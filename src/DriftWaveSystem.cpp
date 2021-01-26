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
#include <MultiRegions/ContField.h>

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

    // Since we are starting from a setup where each field is defined to be a
    // discontinuous field (and thus support DG), the first thing we do is to
    // recreate the phi field so that it is continuous to support the Poisson
    // solve. Note that you can still perform a Poisson solve using a
    // discontinuous field, which is done via the hybridisable discontinuous
    // Galerkin (HDG) approach.
    m_fields[2] = MemoryManager<MultiRegions::ContField>
        ::AllocateSharedPtr(
            m_session, m_graph, m_session->GetVariable(2), true, true);

    // Tell UnsteadySystem to only integrate first two fields in time
    // (i.e. vorticity and density), ignoring electrostatic potential since .
    m_intVariables = { 0, 1 };

    // Assign storage for drift velocity.
    for (int i = 0; i < 2; ++i)
    {
        m_driftVel[i] = Array<OneD, NekDouble>(m_fields[i]->GetNpoints());
    }

    // Load constant alpha from the session file.
    ASSERTL0(m_session->DefinesParameter("alpha"),
             "Session file should define parameter alpha.");
    m_session->LoadParameter("alpha", m_alpha, 1.0);

    // Load constant kappa from the session file.
    ASSERTL0(m_session->DefinesParameter("kappa"),
             "Session file should define parameter kappa.");
    m_session->LoadParameter("kappa", m_kappa, 1.0);

    // Type of advection class to be used. By default, we only support the
    // discontinuous projection, since this is the only approach we're
    // considering for this solver.
    switch(m_projectionType)
    {
        // Discontinuous field
        case MultiRegions::eDiscontinuous:
        {
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
            m_advObject = SolverUtils::
                GetAdvectionFactory().CreateInstance(advName, advName);

            // The advection object needs to know the flux vector being
            // calculated: this is done with a callback.
            m_advObject->SetFluxVector(&DriftWaveSystem::GetFluxVector, this);

            // Repeat the above for the Riemann solver: in this case we use an
            // upwind by default. The solver also needs to know the trace
            // normal, which we again implement using a callback.
            m_session->LoadSolverInfo("UpwindType", riemName, "Upwind");
            m_riemannSolver =
                GetRiemannSolverFactory().CreateInstance(riemName, m_session);
            m_riemannSolver->SetScalar(
                "Vn", &DriftWaveSystem::GetNormalVelocity, this);

            // Tell the advection object about the Riemann solver to use, and
            // then get it set up.
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

    // The m_ode object defines the timestepping to be used, and lives in the
    // SolverUtils::UnsteadySystem class. For explicit solvers, you need to
    // supply a right-hand side function, and a projection function (e.g. for
    // continuous Galerkin this would be an assembly-type operation to ensure
    // C^0 connectivity). These are done again through callbacks.
    m_ode.DefineOdeRhs    (&DriftWaveSystem::ExplicitTimeInt, this);
    m_ode.DefineProjection(&DriftWaveSystem::DoOdeProjection, this);
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
void DriftWaveSystem::ExplicitTimeInt(
    const Array<OneD, const Array<OneD,NekDouble> > &inarray,
          Array<OneD,       Array<OneD,NekDouble> > &outarray,
    const NekDouble time)

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
    m_fields[2]->HelmSolve(inarray[0], m_fields[2]->UpdateCoeffs(),
                           factors);
    m_fields[2]->BwdTrans (m_fields[2]->GetCoeffs(),
                           m_fields[2]->UpdatePhys());

    // Calculate drift velocity v_E: PhysDeriv takes input and computes spatial
    // derivatives.
    m_fields[2]->PhysDeriv(m_fields[2]->GetPhys(), m_driftVel[1], m_driftVel[0]);

    // We frequently use vector math (Vmath) routines for one-line operations
    // like negating entries in a vector.
    Vmath::Neg(nPts, m_driftVel[1], 1);

    // Do advection for zeta, n. The hard-coded '2' here indicates that we
    // should only advect the first two components of inarray.
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
    Vmath::Svtvp(nPts, -m_kappa, m_driftVel[0], 1,
                 outarray[1], 1, outarray[1], 1);
}

/**
 * @brief Perform projection into correct polynomial space.
 *
 * This routine projects the @p inarray input and ensures the @p outarray output
 * lives in the correct space. Since we are hard-coding DG, this corresponds to
 * a simple copy from in to out, since no elemental connectivity is required and
 * the output of the RHS function is polynomial.
 */
void DriftWaveSystem::DoOdeProjection(
    const Array<OneD, const Array<OneD, NekDouble> > &inarray,
          Array<OneD,       Array<OneD, NekDouble> > &outarray,
    const NekDouble time)
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
void DriftWaveSystem::GetFluxVector(
    const Array<OneD, Array<OneD, NekDouble> >               &physfield,
          Array<OneD, Array<OneD, Array<OneD, NekDouble> > > &flux)
{
    ASSERTL1(flux[0].size() == m_driftVel.size(),
             "Dimension of flux array and velocity array do not match");

    int i , j;
    int nq = physfield[0].size();

    for (i = 0; i < flux.size(); ++i)
    {
        for (j = 0; j < flux[0].size(); ++j)
        {
            Vmath::Vmul(nq, physfield[i], 1, m_driftVel[j], 1,
                        flux[i][j], 1);
        }
    }
}

/**
 * @brief Compute the normal advection velocity for this system on the
 * trace/skeleton/edges of the 2D mesh.
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

    for (i = 0; i < m_driftVel.size(); ++i)
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

}
