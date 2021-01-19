///////////////////////////////////////////////////////////////////////////////
//
// File: DriftWaveSystem.cpp
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

#include "DriftWaveSystem.h"
#include <MultiRegions/ContField2D.h>

namespace Nektar
{

std::string DriftWaveSystem::className =
    GetEquationSystemFactory().RegisterCreatorFunction(
        "DriftWaveSystem", DriftWaveSystem::create,
        "System for the Hasegawa-Wakatani equations.");

DriftWaveSystem::DriftWaveSystem(
    const LibUtilities::SessionReaderSharedPtr& session,
    const SpatialDomains::MeshGraphSharedPtr& graph)
    : UnsteadySystem(session, graph),
      AdvectionSystem(session, graph),
      m_driftVel(2)
{
}

void DriftWaveSystem::v_InitObject()
{
    AdvectionSystem::v_InitObject();

    // Redefine phi field so that it is continuous for Poisson solve.
    m_fields[2] = MemoryManager<MultiRegions::ContField2D>
        ::AllocateSharedPtr(
            m_session, m_graph, m_session->GetVariable(2), true, true);

    // Tell UnsteadySystem to only integrate first two fields in time
    // (i.e. vorticity and density).
    m_intVariables.push_back(0);
    m_intVariables.push_back(1);

    // Assign storage for drift velocity.
    for (int i = 0; i < 2; ++i)
    {
        m_driftVel[i] = Array<OneD, NekDouble>(m_fields[i]->GetNpoints());
    }

    // Load constant alpha.
    ASSERTL0(m_session->DefinesParameter("alpha"),
             "Session file should define parameter alpha.");
    m_session->LoadParameter("alpha", m_alpha, 1.0);

    // Load constant kappa.
    ASSERTL0(m_session->DefinesParameter("kappa"),
             "Session file should define parameter kappa.");
    m_session->LoadParameter("kappa", m_kappa, 1.0);

    // Type of advection class to be used
    switch(m_projectionType)
    {
        // Discontinuous field
        case MultiRegions::eDiscontinuous:
        {
            // Do not forwards transform initial condition
            m_homoInitialFwd = false;

            // Define the normal velocity fields
            if (m_fields[0]->GetTrace())
            {
                m_traceVn = Array<OneD, NekDouble>(GetTraceNpoints());
            }

            std::string advName, riemName;
            m_session->LoadSolverInfo(
                "AdvectionType", advName, "WeakDG");
            m_advObject = SolverUtils::
                GetAdvectionFactory().CreateInstance(advName, advName);
            m_advObject->SetFluxVector(
                &DriftWaveSystem::GetFluxVector, this);
            m_session->LoadSolverInfo(
                "UpwindType", riemName, "Upwind");
            m_riemannSolver =
                GetRiemannSolverFactory().CreateInstance(riemName, m_session);
            m_riemannSolver->SetScalar(
                "Vn", &DriftWaveSystem::GetNormalVelocity, this);

            m_advObject->SetRiemannSolver(m_riemannSolver);
            m_advObject->InitObject(m_session, m_fields);
            break;
        }

        default:
        {
            ASSERTL0(false,
                     "Unsupported projection type: only discontinuous"
                     " projection supported.");
            break;
        }
    }

    ASSERTL0(m_explicitAdvection,
             "This solver only supports explicit-in-time advection.");

    m_ode.DefineOdeRhs    (&DriftWaveSystem::ExplicitTimeInt, this);
    m_ode.DefineProjection(&DriftWaveSystem::DoOdeProjection, this);
}

DriftWaveSystem::~DriftWaveSystem()
{

}

void DriftWaveSystem::ExplicitTimeInt(
    const Array<OneD, const Array<OneD,NekDouble> > &inarray,
    Array<OneD,       Array<OneD,NekDouble> > &outarray,
    const NekDouble time)
{
    // Field 0 vorticity
    // Field 1 density
    // Field 2 electric potential

    int i;
    int nPts = GetNpoints();

    // Solve for electric potential.
    StdRegions::ConstFactorMap factors;
    factors[StdRegions::eFactorLambda] = 0.0;

    // Solve for phi
    m_fields[2]->HelmSolve(inarray[0], m_fields[2]->UpdateCoeffs(),
                           NullFlagList, factors);
    m_fields[2]->BwdTrans (m_fields[2]->GetCoeffs(),
                           m_fields[2]->UpdatePhys());

    // Calculate drift velocity v_E
    m_fields[2]->PhysDeriv(m_fields[2]->GetPhys(),
                           m_driftVel[1], m_driftVel[0]);
    Vmath::Neg(nPts, m_driftVel[1], 1);

    // Do advection for zeta, n.
    m_advObject->Advect(2, m_fields, m_driftVel, inarray, outarray, time);

    // Put advection term on the right hand side.
    for (i = 0; i < 2; ++i)
    {
        Vmath::Neg(nPts, outarray[i], 1);
    }

    // Add source term alpha*(phi - n) to right hand side.
    Array<OneD, NekDouble> sourceTerm(nPts);
    Vmath::Vsub(nPts, m_fields[2]->GetPhys(), 1, inarray[1], 1, sourceTerm, 1);
    Vmath::Smul(nPts, m_alpha, sourceTerm, 1, sourceTerm, 1);
    Vmath::Vadd(nPts, sourceTerm, 1, outarray[0], 1, outarray[0], 1);
    Vmath::Vadd(nPts, sourceTerm, 1, outarray[1], 1, outarray[1], 1);

    // Add source term -kappa * d(phi)/dy to n equation.
    Vmath::Svtvp(nPts, -m_kappa, m_driftVel[0], 1, outarray[1], 1, outarray[1], 1);
}

/**
 *
 */
void DriftWaveSystem::DoOdeProjection(
    const Array<OneD, const Array<OneD, NekDouble> > &inarray,
    Array<OneD,       Array<OneD, NekDouble> > &outarray,
    const NekDouble time)
{
    int nvariables = inarray.num_elements(), npoints = GetNpoints();
    SetBoundaryConditions(time);

    for (int i = 0; i < nvariables; ++i)
    {
        Vmath::Vcopy(npoints, inarray[i], 1, outarray[i], 1);
    }
}

void DriftWaveSystem::GetFluxVector(
    const Array<OneD, Array<OneD, NekDouble> >               &physfield,
          Array<OneD, Array<OneD, Array<OneD, NekDouble> > > &flux)
{
    ASSERTL1(flux[0].num_elements() == m_driftVel.num_elements(),
             "Dimension of flux array and velocity array do not match");

    int i , j;
    int nq = physfield[0].num_elements();

    for (i = 0; i < flux.num_elements(); ++i)
    {
        for (j = 0; j < flux[0].num_elements(); ++j)
        {
            Vmath::Vmul(nq, physfield[i], 1, m_driftVel[j], 1,
                        flux[i][j], 1);
        }
    }
}

/**
 * @brief Get the normal velocity for the linear advection equation.
 */
Array<OneD, NekDouble> &DriftWaveSystem::GetNormalVelocity()
{
    // Number of trace (interface) points
    int i;
    int nTracePts = GetTraceNpoints();

    // Auxiliary variable to compute the normal velocity
    Array<OneD, NekDouble> tmp(nTracePts);

    // Reset the normal velocity
    Vmath::Zero(nTracePts, m_traceVn, 1);

    for (i = 0; i < m_driftVel.num_elements(); ++i)
    {
        m_fields[0]->ExtractTracePhys(m_driftVel[i], tmp);

        Vmath::Vvtvp(nTracePts,
                     m_traceNormals[i], 1,
                     tmp,               1,
                     m_traceVn,         1,
                     m_traceVn,         1);
    }

    return m_traceVn;
}

void DriftWaveSystem::v_GenerateSummary(SolverUtils::SummaryList& s)
{
    UnsteadySystem::v_GenerateSummary(s);
}
}
