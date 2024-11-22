///////////////////////////////////////////////////////////////////////////////
//
// File: RogersRicci1D.cpp
//
// For more information, please see: http://www.nektar.info
//
// The MIT License
//
// Copyright (c) 2006 Division of Applied Mathematics, Brown University (USA),
// Department of Aeronautics, Imperial College London (UK), and Scientific
// Computing and Imaging Institute, University of Utah (USA).
//
// License for the specific language governing rights and limitations under
// Permission is hereby granted, free of charge, to any person obtaining a
// copy of this software and associated documentation files (the "Software"),
// to deal in the Software without restriction, including without limitation
// the rights to use, copy, modify, merge, publish, distribute, sublicense,
// and/or sell copies of the Software, and to permit persons to whom the
// Software is furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included
// in all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
// OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
// THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
// FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
// DEALINGS IN THE SOFTWARE.
//
// Description: Image warping solve routines
//
///////////////////////////////////////////////////////////////////////////////

#include "RogersRicci1D.h"
#include <MultiRegions/ContField.h>
#include <SolverUtils/RiemannSolvers/RiemannSolver.h>

namespace Nektar
{
std::string RogersRicci1D::className =
    GetEquationSystemFactory().RegisterCreatorFunction(
        "RogersRicci1D", RogersRicci1D::create,
        "System for the Rogers-Ricci 1D system of equations.");

RogersRicci1D::RogersRicci1D(
    const LibUtilities::SessionReaderSharedPtr &session,
    const SpatialDomains::MeshGraphSharedPtr &graph)
    : UnsteadySystem(session, graph), AdvectionSystem(session, graph),
      m_advVel(1)
{
}

void RogersRicci1D::v_InitObject(bool DeclareField)
{
    AdvectionSystem::v_InitObject(DeclareField);

    m_npts = m_fields[0]->GetNpoints();

    ASSERTL0(m_fields.size() == 2,
             "Incorrect number of variables detected (expected 2): check your "
             "session file.");

    // Store mesh dimension for easy retrieval later.
    m_ndims = m_graph->GetMeshDimension();
    ASSERTL0(m_ndims == 1, "Eq sys only supports 1D mesh");

    m_intVariables = {n_idx, u_idx};

    m_session->LoadParameter("n_star", this->n_star, 0.06);
    m_session->LoadParameter("T", this->T, 1.0);
    m_session->LoadParameter("tau", this->tau, 1.0);

    m_advVel[0] = Array<OneD, NekDouble>(m_npts);

    switch (m_projectionType)
    {
            // case MultiRegions::eGalerkin:
            // {

            //     break;
            // }

        case MultiRegions::eDiscontinuous:
        {
            // Make u continuous
            m_fields[u_idx] = MemoryManager<MR::ContField>::AllocateSharedPtr(
                m_session, m_graph, m_session->GetVariable(u_idx), true, true);

            // Do not forwards transform initial condition.
            m_homoInitialFwd = false;

            // Define the normal velocity fields.
            if (m_fields[0]->GetTrace())
            {
                m_traceVn = Array<OneD, NekDouble>(GetTraceNpoints());
            }

            // The remainder of this code is fairly generic boilerplate for the
            // DG setup.
            std::string advName, riemName;

            // Load what type of advection we want to use -- in theory we also
            // support flux reconstruction for quad-based meshes, or you can use
            // a standard convective term if you were fully continuous in
            // space. Default is DG.
            m_session->LoadSolverInfo("AdvectionType", advName, "WeakDG");

            // Create an advection object of the type above using the factory
            // pattern.
            m_advObject = SolverUtils::GetAdvectionFactory().CreateInstance(
                advName, advName);

            // The advection object needs to know the flux vector being
            // calculated: this is done with a callback.
            m_advObject->SetFluxVector(&RogersRicci1D::GetFluxVector, this);

            // Repeat the above for the Riemann solver: in this case we use an
            // upwind by default. The solver also needs to know the trace
            // normal, which we again implement using a callback.
            m_session->LoadSolverInfo("UpwindType", riemName, "Upwind");
            m_riemannSolver =
                GetRiemannSolverFactory().CreateInstance(riemName, m_session);
            m_riemannSolver->SetScalar("Vn", &RogersRicci1D::GetNormalVelocity,
                                       this);

            // Tell the advection object about the Riemann solver to use, and
            // then get it set up.
            m_advObject->SetRiemannSolver(m_riemannSolver);
            m_advObject->InitObject(m_session, m_fields);
            break;
        }

        default:
        {
            ASSERTL0(false, "Unsupported projection type.");
            break;
        }
    }

    m_ode.DefineOdeRhs(&RogersRicci1D::ExplicitTimeInt, this);
    m_ode.DefineProjection(&RogersRicci1D::DoOdeProjection, this);

    if (!m_explicitAdvection)
    {
        m_implHelper = std::make_shared<ImplicitHelper>(
            m_session, m_fields, m_ode, m_intVariables.size());
        m_implHelper->InitialiseNonlinSysSolver();
        m_ode.DefineImplicitSolve(&ImplicitHelper::ImplicitTimeInt,
                                  m_implHelper);
    }

    // Set up user-defined boundary conditions, if any were specified
    for (int ifld = 0; ifld < m_fields.size(); ifld++)
    {
        int cnt = 0;
        for (int icnd = 0; icnd < m_fields[ifld]->GetBndConditions().size();
             ++icnd)
        {
            SpatialDomains::BoundaryConditionShPtr cnd =
                m_fields[ifld]->GetBndConditions()[icnd];
            if (cnd->GetBoundaryConditionType() != SpatialDomains::ePeriodic)
            {
                std::string type = cnd->GetUserDefined();
                if (!type.empty())
                {
                    CustomBCsSharedPtr BCs_instance =
                        GetCustomBCsFactory().CreateInstance(
                            type, m_session, m_fields, m_traceNormals, ifld,
                            m_spacedim, icnd, cnt, cnd);
                    this->custom_BCs.push_back(BCs_instance);
                }
                cnt +=
                    m_fields[ifld]->GetBndCondExpansions()[icnd]->GetExpSize();
            }
        }
    }
}

/**
 * @brief Evaluate the right-hand side of the ODE system used to integrate in
 * time.
 */
void RogersRicci1D::ExplicitTimeInt(
    const Array<OneD, const Array<OneD, NekDouble>> &inarray,
    Array<OneD, Array<OneD, NekDouble>> &outarray, const NekDouble time)
{
    /**
     * ndot = - (nu)' + n_star
     * udot = - u u' - tau T n'/n
     *
     * n is DG, u is CG; Advect n only
     */

    Array<OneD, NekDouble> n     = inarray[n_idx];
    Array<OneD, NekDouble> u     = inarray[u_idx];
    Array<OneD, NekDouble> n_out = outarray[n_idx];
    Array<OneD, NekDouble> u_out = outarray[u_idx];

    // Advect n only
    Vmath::Vcopy(m_npts, u, 1, m_advVel[0], 1);
    m_advObject->Advect(1, m_fields, m_advVel, inarray, outarray, time);
    Vmath::Neg(m_npts, outarray[0], 1);

    Array<OneD, NekDouble> grad_n(m_npts), grad_u(m_npts);
    m_fields[n_idx]->PhysDeriv(n, grad_n);
    m_fields[u_idx]->PhysDeriv(u, grad_u);

    for (auto ii = 0; ii < m_npts; ii++)
    {
        n_out[ii] += this->n_star;
        u_out[ii] =
            -u[ii] * grad_u[ii] - this->tau * this->T * grad_n[ii] / n[ii];
    }
}

void RogersRicci1D::DoOdeProjection(
    const Array<OneD, const Array<OneD, NekDouble>> &inarray,
    Array<OneD, Array<OneD, NekDouble>> &outarray, const NekDouble time)
{
    int nvariables = inarray.size(), npoints = GetNpoints();

    switch (m_projectionType)
    {
        case MultiRegions::eDiscontinuous:
        {
            for (int i = 0; i < nvariables; ++i)
            {
                Vmath::Vcopy(npoints, inarray[i], 1, outarray[i], 1);
            }
            break;
        }
        case MultiRegions::eGalerkin:
        {
            int ncoeffs = m_fields[0]->GetNcoeffs();
            Array<OneD, NekDouble> coeffs(ncoeffs, 0.0);
            for (int i = 0; i < nvariables; ++i)
            {
                m_fields[i]->FwdTrans(inarray[i], coeffs);
                m_fields[i]->BwdTrans(coeffs, outarray[i]);
            }
            break;
        }
    }
    SetBoundaryConditions(outarray, time);
}

/**
 * @brief Compute the flux vector for this system.
 */
void RogersRicci1D::GetFluxVector(
    const Array<OneD, Array<OneD, NekDouble>> &physfield,
    Array<OneD, Array<OneD, Array<OneD, NekDouble>>> &flux)
{
    ASSERTL1(flux[0].size() == m_advVel.size(),
             "Dimension of flux array and velocity array do not match");

    int nq = physfield[0].size();

    for (int i = 0; i < flux.size(); ++i)
    {
        for (int j = 0; j < flux[0].size(); ++j)
        {
            Vmath::Vmul(nq, physfield[i], 1, m_advVel[j], 1, flux[i][j], 1);
        }
    }
}

/**
 * @brief Compute the normal advection velocity for this system on the
 * trace/skeleton/edges of the 2D mesh.
 */
Array<OneD, NekDouble> &RogersRicci1D::GetNormalVelocity()
{
    // Number of trace (interface) points
    int nTracePts = GetTraceNpoints();

    // Auxiliary variable to compute the normal velocity
    Array<OneD, NekDouble> tmp(nTracePts);

    // Reset the normal velocity
    Vmath::Zero(nTracePts, m_traceVn, 1);

    // Compute dot product of velocity along trace with trace normals. Store in
    // m_traceVn.
    for (int i = 0; i < m_advVel.size(); ++i)
    {
        m_fields[0]->ExtractTracePhys(m_advVel[i], tmp);

        Vmath::Vvtvp(nTracePts, m_traceNormals[i], 1, tmp, 1, m_traceVn, 1,
                     m_traceVn, 1);
    }

    return m_traceVn;
}

void RogersRicci1D::SetBoundaryConditions(
    Array<OneD, Array<OneD, NekDouble>> &physarray, NekDouble time)
{
    if (!this->custom_BCs.empty())
    {
        for (auto &bc : this->custom_BCs)
        {
            if (time > 1e-3)
            {
                bc->Apply(physarray, time);
            }
        }
    }
}

} // namespace Nektar
