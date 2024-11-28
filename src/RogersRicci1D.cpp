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
// Description: A greatly simplified version of the Rogers & Ricci system framed
// as a 1D problem.
//
///////////////////////////////////////////////////////////////////////////////

#include "RogersRicci1D.h"
#include <SolverUtils/RiemannSolvers/RiemannSolver.h>

namespace Nektar
{
std::string RogersRicci1D::className =
    GetEquationSystemFactory().RegisterCreatorFunction(
        "RogersRicci1D", RogersRicci1D::create,
        "System for the Rogers-Ricci 1D system of equations.");

class CustomUpwindSolver : public SolverUtils::RiemannSolver
{
public:
    static RiemannSolverSharedPtr create(
        const LibUtilities::SessionReaderSharedPtr &pSession)
    {
        return RiemannSolverSharedPtr(new CustomUpwindSolver(pSession));
    }

    static std::string solver_name;

    CustomUpwindSolver(const LibUtilities::SessionReaderSharedPtr &pSession)
        : SolverUtils::RiemannSolver(pSession)
    {
    }

    void set_pt_types(std::vector<TracePtType> trace_pt_types)
    {
        this->trace_pt_types = trace_pt_types;
    }
protected:
    NekDouble calc_alpha(const NekDouble &u_L, const NekDouble &u_R,
                         const NekDouble &cs)
    {
        NekDouble C_L = std::max(std::abs(u_L - cs), std::abs(u_L + cs));
        NekDouble C_R = std::max(std::abs(u_R - cs), std::abs(u_R + cs));
        return std::max(C_L, C_R);
    }

    NekDouble calc_trace_flux(const NekDouble &q_L, const NekDouble &q_R,
                              const NekDouble &qF_L, const NekDouble &qF_R,
                              const NekDouble &alpha, const NekDouble norm = 1)
    {
        return norm * 0.5 * (qF_L + qF_R) - 0.5 * alpha * (q_R - q_L);
    }

    virtual void v_Solve(
        const int nDim, const Array<OneD, const Array<OneD, NekDouble>> &Fwd,
        const Array<OneD, const Array<OneD, NekDouble>> &Bwd,
        Array<OneD, Array<OneD, NekDouble>> &flux) override final
    {

        ASSERTL0(CheckScalars("flow_trace_norm"),
                 "flow_trace_norm array not defined.");
        const Array<OneD, NekDouble> &flow_norms =
            m_scalars["flow_trace_norm"]();

        const int n_idx       = 0;
        const int u_idx       = 1;
        const int n_trace_pts = Fwd[0].size();
        ASSERTL0(this->trace_pt_types.size() == n_trace_pts,
                 "trace_pt_types vector has the wrong number of points");

        // Isothermal for now
        NekDouble T = 1.0;

        for (int j = 0; j < n_trace_pts; ++j)
        {
            // Extract variables
            NekDouble n_L = Fwd[n_idx][j];
            NekDouble u_L = Fwd[u_idx][j];
            NekDouble n_R = Bwd[n_idx][j];
            NekDouble u_R = Bwd[u_idx][j];

            // Set sound speed
            NekDouble cs = sqrt(T);

            // Special handling for boundary points
            if (this->trace_pt_types[j] == eOutflowBdyLow)
            {
                n_R = n_L;  // (zero Neumann)
                u_R = -1.0; // (Dirichlet)
            }
            else if (this->trace_pt_types[j] == eOutflowBdyHigh)
            {
                n_R = n_L; // (zero Neumann)
                u_R = 1.0; // (Dirichlet)
            }

            // Compute left and right state fluxes, accounting for grad(P) term
            NekDouble flux_n_L = n_L * u_L;
            NekDouble flux_n_R = n_R * u_R;
            NekDouble flux_u_L = 0.5 * u_L * u_L + T * log(n_L);
            NekDouble flux_u_R = 0.5 * u_R * u_R + T * log(n_R);

            NekDouble alpha = calc_alpha(u_L, u_R, cs);
            flux[n_idx][j]  = calc_trace_flux(n_L, n_R, flux_n_L, flux_n_R,
                                             alpha, flow_norms[j]);
            flux[u_idx][j]  = calc_trace_flux(u_L, u_R, flux_u_L, flux_u_R,
                                             alpha, flow_norms[j]);
        }
    }

private:
    std::vector<TracePtType> trace_pt_types;
};

std::string CustomUpwindSolver::solver_name =
    SolverUtils::GetRiemannSolverFactory().RegisterCreatorFunction(
        "CustomUpwind", CustomUpwindSolver::create,
        "Customised upwind solver for 1DRR");

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
        case MultiRegions::eDiscontinuous:
        {

            // Do not forwards transform initial condition.
            m_homoInitialFwd = false;

            // Define the normal velocity fields.
            if (m_fields[0]->GetTrace())
            {
                int nTrace           = GetTraceNpoints();
                m_traceVn            = Array<OneD, NekDouble>(nTrace);
                this->trace_pt_types = std::vector(nTrace, eNormalPt);
            }
            this->pt_types = std::vector(m_npts, eNormalPt);

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
            m_riemannSolver->SetScalar(
                "flow_trace_norm", &RogersRicci1D::GetOutflowTraceNormal, this);

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
    // Also populate trace_pt_types array (passed to the riemann solver) and
    // pt_types array (used in DoOdeProjection)
    const Array<OneD, const int> &traceBndMap = m_fields[0]->GetTraceBndMap();
    Array<OneD, int> ElmtID, EdgeID;
    m_fields[0]->GetBoundaryToElmtMap(ElmtID, EdgeID);
    for (int ifld = 0; ifld < m_fields.size(); ifld++)
    {
        // Arrays of BCs for this field and corresponding expansion lists
        auto BCs     = m_fields[ifld]->GetBndConditions();
        auto BC_exps = m_fields[ifld]->GetBndCondExpansions();
        int cnt      = 0;
        for (int region = 0; region < BCs.size(); ++region)
        {
            auto BC_exp                               = BC_exps[region];
            int nExp                                  = BC_exp->GetExpSize();
            SpatialDomains::BoundaryConditionShPtr BC = BCs[region];

            // For all types other than periodic:
            if (BC->GetBoundaryConditionType() != SpatialDomains::ePeriodic)
            {
                // For user-defined BC types, intantiate and add to the
                // class-level BCs array
                std::string user_type = BC->GetUserDefined();
                if (!user_type.empty())
                {
                    CustomBCsSharedPtr BCs_instance =
                        GetCustomBCsFactory().CreateInstance(
                            user_type, m_session, m_fields, m_traceNormals,
                            ifld, m_spacedim, region, cnt, BC);
                    this->custom_BCs.push_back(BCs_instance);

                    // For outflow BCs type, label the appropriate trace pts and
                    // regular field pts as outflow
                    if (user_type == "Outflow")
                    {
                        // Use Dirichlet "equation" (VALUE string) to identify
                        // which outflow boundary we're on
                        auto dirichlet_bc = std::static_pointer_cast<
                            SpatialDomains::DirichletBoundaryCondition>(BC);
                        ASSERTL0(dirichlet_bc,
                                 "Expected 'Outflow' BCs to have type "
                                 "Dirichlet");
                        LibUtilities::Equation eqn =
                            dirichlet_bc->m_dirichletCondition;
                        TracePtType pt_type = eqn.Evaluate() < 0
                                                  ? eOutflowBdyLow
                                                  : eOutflowBdyHigh;

                        for (int exp_idx = 0; exp_idx < nExp; ++exp_idx)
                        {
                            // Set label in global points array
                            int bnd_cumu_idx_start = cnt + exp_idx;
                            int npts = BC_exp->GetExp(exp_idx)->GetTotPoints();

                            for (int iq = 0; iq < npts; iq++)
                            {
                                int bnd_cumu_idx = bnd_cumu_idx_start + iq;
                                int foo          = traceBndMap[bnd_cumu_idx];
                                int idx_in_trace =
                                    m_fields[ifld]->GetTrace()->GetPhys_Offset(
                                        foo);
                                int elementID = ElmtID[bnd_cumu_idx];
                                int edgeID    = EdgeID[bnd_cumu_idx];
                                auto element =
                                    m_fields[ifld]->GetExp(elementID);
                                int npts_element = element->GetTotPoints();

                                Array<OneD, int> edge_pt_indices;
                                element->GetTracePhysMap(edgeID,
                                                         edge_pt_indices);
                                for (auto &edge_pt_idx : edge_pt_indices)
                                {
                                    int phys_idx =
                                        m_fields[ifld]->GetPhys_Offset(
                                            elementID) +
                                        edge_pt_idx;
                                    pt_types[phys_idx] = pt_type;
                                }
                            }
                            auto trace_bdy_map_offset =
                                m_fields[0]->GetTrace()->GetPhys_Offset(
                                    traceBndMap[cnt + exp_idx]);
                            // Set outflow label for all of the pts
                            for (int pt_idx = 0; pt_idx < npts; ++pt_idx)
                            {
                                int trace_idx = trace_bdy_map_offset + pt_idx;
                                trace_pt_types[trace_idx] = pt_type;
                            }
                        }
                    }
                }
                cnt += nExp;
            }
        }
    }

    std::dynamic_pointer_cast<CustomUpwindSolver>(this->m_riemannSolver)
        ->set_pt_types(this->trace_pt_types);
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
     * All terms other than the density source are handled by setting
     * appropriate fluxes in the advection operation.
     */

    Array<OneD, NekDouble> n     = inarray[n_idx];
    Array<OneD, NekDouble> u     = inarray[u_idx];
    Array<OneD, NekDouble> n_out = outarray[n_idx];
    Array<OneD, NekDouble> u_out = outarray[u_idx];

    int nadvect = 2;
    m_advObject->Advect(nadvect, m_fields, m_advVel, inarray, outarray, time);
    for (auto ii = 0; ii < nadvect; ii++)
    {
        Vmath::Neg(m_npts, outarray[ii], 1);
    }

    // Add density source
    for (auto ii = 0; ii < m_npts; ii++)
    {
        n_out[ii] += this->n_star;
    }
}

void RogersRicci1D::DoOdeProjection(
    const Array<OneD, const Array<OneD, NekDouble>> &inarray,
    Array<OneD, Array<OneD, NekDouble>> &outarray, const NekDouble time)
{
    int nvariables = inarray.size();

    for (int i = 0; i < nvariables; ++i)
    {
        Vmath::Vcopy(m_npts, inarray[i], 1, outarray[i], 1);
    }
    EquationSystem::SetBoundaryConditions(time);

    // Fix u for points on the outflow boundaries
    for (auto jj = 0; jj < m_npts; jj++)
    {
        if (this->pt_types[jj] == eOutflowBdyLow)
        {
            outarray[this->u_idx][jj] = -1.0;
        }
        else if (this->pt_types[jj] == eOutflowBdyHigh)
        {
            outarray[this->u_idx][jj] = 1.0;
        }
    }
}

/**
 * @brief Compute the flux vector for this system.
 */
void RogersRicci1D::GetFluxVector(
    const Array<OneD, Array<OneD, NekDouble>> &physfield,
    Array<OneD, Array<OneD, Array<OneD, NekDouble>>> &flux)
{
    int ndims = flux[0].size();
    ASSERTL1(ndims == m_advVel.size(),
             "Dimension of flux array and velocity array do not match");

    int nq = physfield[0].size();

    for (size_t p = 0; p < this->m_npts; ++p)
    {
        // Flux for the n eqn
        flux[0][0][p] = physfield[0][p] * physfield[1][p];
        // Flux for the u eqn
        flux[1][0][p] = 0.5 * physfield[1][p] * physfield[1][p] +
                        this->tau * this->T * std::log(physfield[0][p]);
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

    for (int i = 0; i < nTracePts; ++i)
    {
        m_traceVn[i] = m_traceNormals[0][i];
    }

    return m_traceVn;
}

Array<OneD, NekDouble> &RogersRicci1D::GetOutflowTraceNormal()
{
    int flow_dir;
    switch (m_ndims)
    {
        case 1:
        {
            flow_dir = 0;
            break;
        }
        case 3:
        {
            flow_dir = 2;
            break;
        }
        default:
        {
            ASSERTL0(false, "GetOutflowTraceNormal: Only set up for 1 or 3 "
                            "dimensional meshes.");
            break;
        }
    }
    return m_traceNormals[flow_dir];
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
