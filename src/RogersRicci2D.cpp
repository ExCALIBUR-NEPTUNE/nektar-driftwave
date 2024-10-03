///////////////////////////////////////////////////////////////////////////////
//
// File: RogersRicci2D.cpp
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

#include "RogersRicci2D.h"
#include <MultiRegions/ContField.h>

namespace Nektar
{

std::string RogersRicci2D::className =
    GetEquationSystemFactory().RegisterCreatorFunction(
        "RogersRicci2D", RogersRicci2D::create,
        "System for the Rogers-Ricci 2D system of equations.");

RogersRicci2D::RogersRicci2D(
    const LibUtilities::SessionReaderSharedPtr &session,
    const SpatialDomains::MeshGraphSharedPtr &graph)
    : UnsteadySystem(session, graph), AdvectionSystem(session, graph),
      m_driftVel(2)
{
    // Set up constants
    m_c["B"] = c("omega_ci") * c("m_i") * c("q_E");
    m_c["c_s0"] = sqrt(c("T_e0") / c("m_i"));
    m_c["rho_s0"] = c("c_s0") / c("omega_ci");
    m_c["S_0n"] = 0.03 * c("n0") * c("c_s0") / c("R");
    m_c["S_0T"] = 0.03 * c("T_e0") * c("c_s0") / c("R");
    m_c["omega"] = 1.5 * c("R") / c("L_z");
}

void RogersRicci2D::v_InitObject(bool DeclareField)
{
    AdvectionSystem::v_InitObject(DeclareField);

    ASSERTL0(m_fields.size() == 4,
             "Incorrect number of variables detected (expected 4): check your "
             "session file.");

    m_fields[3] = MemoryManager<MultiRegions::ContField>::AllocateSharedPtr(
        m_session, m_graph, m_session->GetVariable(3), true, true);
    m_intVariables = {0, 1, 2};

    // Assign storage for drift velocity.
    for (int i = 0; i < 2; ++i)
    {
        m_driftVel[i] = Array<OneD, NekDouble>(m_fields[i]->GetNpoints());
    }

    switch (m_projectionType)
    {
        case MultiRegions::eDiscontinuous:
        {
            m_homoInitialFwd = false;

            if (m_fields[0]->GetTrace())
            {
                m_traceVn = Array<OneD, NekDouble>(GetTraceNpoints());
            }

            std::string advName, riemName;
            m_session->LoadSolverInfo("AdvectionType", advName, "WeakDG");
            m_advObject = SolverUtils::GetAdvectionFactory().CreateInstance(
                advName, advName);
            m_advObject->SetFluxVector(&RogersRicci2D::GetFluxVector, this);

            m_session->LoadSolverInfo("UpwindType", riemName, "Upwind");
            m_riemannSolver =
                GetRiemannSolverFactory().CreateInstance(riemName, m_session);
            m_riemannSolver->SetScalar(
                "Vn", &RogersRicci2D::GetNormalVelocity, this);

            // Tell the advection object about the Riemann solver to use, and
            // then get it set up.
            m_advObject->SetRiemannSolver(m_riemannSolver);
            m_advObject->InitObject(m_session, m_fields);
            break;
        }

        default:
        {
            ASSERTL0(false, "Unsupported projection type: only discontinuous"
                            " projection supported.");
            break;
        }
    }

    m_ode.DefineOdeRhs(&RogersRicci2D::ExplicitTimeInt, this);
    m_ode.DefineProjection(&RogersRicci2D::DoOdeProjection, this);

    if (!m_explicitAdvection)
    {
        m_implHelper = std::make_shared<ImplicitHelper>(
            m_session, m_fields, m_ode, 3);
        m_implHelper->InitialiseNonlinSysSolver();
        m_ode.DefineImplicitSolve(
            &ImplicitHelper::ImplicitTimeInt, m_implHelper);
    }

    const int nPts = m_fields[0]->GetNpoints();
    m_x = Array<OneD, NekDouble>(nPts);
    m_y = Array<OneD, NekDouble>(nPts);
    m_r = Array<OneD, NekDouble>(nPts);

    m_fields[0]->GetCoords(m_x, m_y);

    for (int i = 0; i < nPts; ++i)
    {
        m_r[i] = sqrt(m_x[i] * m_x[i] + m_y[i] * m_y[i]);
    }
}

/**
 * @brief Evaluate the right-hand side of the ODE system used to integrate in
 * time.
 *
 * This routine performs the bulk of the work in this class, and essentially
 * computes the right hand side term of the generalised ODE system
 *
 * \f\[ \frac{\partial \mathbf{u}}{\partial t} = \mathbf{R}(\mathbf{u}) \f\]
 *
 * The order of operations is as follows:
 *
 * - First, compute the electrostatic potential \f$ \phi \f$, given the
 * - Using this, compute the drift velocity \f$ (\partial_y\phi,
 *   -\partial_x\phi).
 * - Then evaluate the \f$ \nabla\cdot\mathbf{F} \f$ operator using the
 *   advection object #m_advObject.
 * - Finally put this on the right hand side and evaluate the source terms for
 *   each field.
 *
 * The assumption here is that fields are ordered inside `m_fields` so that
 * field 0 is vorticity \f$ \zeta \f$, field 1 is number density \f$ n \f$, and
 * field 2 is electrostatic potential. Only \f$ \zeta \f$ and \f$ n \f$ are time
 * integrated.
 *
 * @param inarray    Array containing each field's current state.
 * @param outarray   The result of the right-hand side operator for each field
 *                   being time integrated.
 * @param time       Current value of time.
 */
void RogersRicci2D::ExplicitTimeInt(
    const Array<OneD, const Array<OneD, NekDouble>> &inarray,
    Array<OneD, Array<OneD, NekDouble>> &outarray, const NekDouble time)
{
    // nPts below corresponds to the total number of solution/integration
    // points: i.e. number of elements * quadrature points per element.
    int i, nPts = GetNpoints();

    // Set up factors for electrostatic potential solve. We support a generic
    // Helmholtz solve of the form (\nabla^2 - \lambda) u = f, so this sets
    // \lambda to zero.
    StdRegions::ConstFactorMap factors;
    factors[StdRegions::eFactorLambda] = 0.0;

    // Solve for phi. Output of this routine is in coefficient (spectral) space,
    // so backwards transform to physical space since we'll need that for the
    // advection step & computing drift velocity.
    m_fields[3]->HelmSolve(inarray[0], m_fields[3]->UpdateCoeffs(), factors);
    m_fields[3]->BwdTrans(m_fields[3]->GetCoeffs(), m_fields[3]->UpdatePhys());

    // Calculate drift velocity v_E: PhysDeriv takes input and computes spatial
    // derivatives.
    m_fields[3]->PhysDeriv(m_fields[3]->GetPhys(), m_driftVel[1],
                           m_driftVel[0]);

    // We frequently use vector math (Vmath) routines for one-line operations
    // like negating entries in a vector.
    Vmath::Neg(nPts, m_driftVel[1], 1);

    // Do advection for zeta, n. The hard-coded '2' here indicates that we
    // should only advect the first two components of inarray.
    m_advObject->Advect(3, m_fields, m_driftVel, inarray, outarray, time);

    Array<OneD, NekDouble> n = inarray[0];
    Array<OneD, NekDouble> T_e = inarray[1];
    Array<OneD, NekDouble> w = inarray[2];
    Array<OneD, NekDouble> phi = m_fields[3]->UpdatePhys();

    // Put advection term on the right hand side.
    const NekDouble rho_s0 = 1.2e-2;
    const NekDouble r_s = 1.0; // can't find this in list of constants
    const NekDouble L_s = 18.0; // not sure this is right -- not listed in table?

    for (i = 0; i < nPts; ++i)
    {
        NekDouble et = exp(3 - phi[i] / T_e[i]);
        NekDouble st = 0.03 * (1.0 - tanh(rho_s0 * m_r[i] - r_s) / L_s);
        outarray[0][i] = -40 * outarray[0][i] - 1.0/24.0 * et * n[i] + st; // add source term
        outarray[1][i] = -40 * outarray[1][i] - 1.0/36.0 * (1.71 * et - 0.71) * T_e[i] + st; // add source term
        outarray[2][i] = -40 * outarray[2][i] + 1.0/24.0 * (1 - et); // should be ok
    }
}

/**
 * @brief Perform projection into correct polynomial space.
 *
 * This routine projects the @p inarray input and ensures the @p outarray output
 * lives in the correct space. Since we are hard-coding DG, this corresponds to
 * a simple copy from in to out, since no elemental connectivity is required and
 * the output of the RHS function is polynomial.
 */
void RogersRicci2D::DoOdeProjection(
    const Array<OneD, const Array<OneD, NekDouble>> &inarray,
    Array<OneD, Array<OneD, NekDouble>> &outarray, const NekDouble time)
{
    int nvariables = inarray.size(), npoints = GetNpoints();
    SetBoundaryConditions(time);

    for (int i = 0; i < nvariables; ++i)
    {
        Vmath::Vcopy(npoints, inarray[i], 1, outarray[i], 1);
    }
}

/**
 * @brief Compute the flux vector for this system.
 */
void RogersRicci2D::GetFluxVector(
    const Array<OneD, Array<OneD, NekDouble>> &physfield,
    Array<OneD, Array<OneD, Array<OneD, NekDouble>>> &flux)
{
    ASSERTL1(flux[0].size() == m_driftVel.size(),
             "Dimension of flux array and velocity array do not match");

    int nq = physfield[0].size();

    for (int i = 0; i < flux.size(); ++i)
    {
        for (int j = 0; j < flux[0].size(); ++j)
        {
            Vmath::Vmul(nq, physfield[i], 1, m_driftVel[j], 1, flux[i][j], 1);
        }
    }
}

/**
 * @brief Compute the normal advection velocity for this system on the
 * trace/skeleton/edges of the 2D mesh.
 */
Array<OneD, NekDouble> &RogersRicci2D::GetNormalVelocity()
{
    // Number of trace (interface) points
    int nTracePts = GetTraceNpoints();

    // Auxiliary variable to compute the normal velocity
    Array<OneD, NekDouble> tmp(nTracePts);

    // Reset the normal velocity
    Vmath::Zero(nTracePts, m_traceVn, 1);

    // Compute dot product of velocity along trace with trace normals. Store in
    // m_traceVn.
    for (int i = 0; i < m_driftVel.size(); ++i)
    {
        m_fields[0]->ExtractTracePhys(m_driftVel[i], tmp);

        Vmath::Vvtvp(nTracePts, m_traceNormals[i], 1, tmp, 1, m_traceVn, 1,
                     m_traceVn, 1);
    }

    return m_traceVn;
}

} // namespace Nektar
