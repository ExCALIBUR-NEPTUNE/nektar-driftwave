///////////////////////////////////////////////////////////////////////////////
//
// File RogersRicci1D.h
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
// Description: Unsteady driftwave solve routines
//
///////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_ROGERSRICCI1D_H
#define NEKTAR_ROGERSRICCI1D_H

#include <SolverUtils/AdvectionSystem.h>
#include <SolverUtils/RiemannSolvers/RiemannSolver.h>

#include "BoundaryConditions/CustomBCs.h"
#include "ImplicitHelper.hpp"

using namespace Nektar::SolverUtils;

namespace Nektar
{

// Enum used to identify trace points on outflow boundaries
enum TracePtType
{
    eOutflowBdyLow,
    eOutflowBdyHigh,
    eNormalPt,
};

/**
 * @brief Equation system for 1D, reduced R&R.
 */
class RogersRicci1D : public AdvectionSystem
{
public:
    // Friend class to allow the memory manager to allocate shared pointers of
    // this class.
    friend class MemoryManager<RogersRicci1D>;

    /// Creates an instance of this class. This static method is registered with
    /// a factory.
    static SolverUtils::EquationSystemSharedPtr create(
        const LibUtilities::SessionReaderSharedPtr &session,
        const SpatialDomains::MeshGraphSharedPtr &graph)
    {
        SolverUtils::EquationSystemSharedPtr p =
            MemoryManager<RogersRicci1D>::AllocateSharedPtr(session, graph);
        p->InitObject();
        return p;
    }

    /// Name of class, used to statically initialise a function pointer for the
    /// create method above.
    static std::string className;

    // Require a fixed variable order; use these indices for clarity
    static constexpr int n_idx = 0;
    static constexpr int u_idx = 1;

    /// Default destructor.
    virtual ~RogersRicci1D() = default;

protected:
    Array<OneD, Array<OneD, NekDouble>> m_advVel;
    /// Storage for the dot product of drift velocity with element edge normals,
    /// required for the DG formulation.
    Array<OneD, NekDouble> m_traceVn;
    /// A SolverUtils::Advection object, which abstracts the calculation of the
    /// \f$ \nabla\cdot\mathbf{F} \f$ operator using different approaches.
    AdvectionSharedPtr m_advObject;
    /// A Riemann solver object to solve numerical fluxes arising from DG: in
    /// this case a simple upwind.
    RiemannSolverSharedPtr m_riemannSolver;
    /// Helper object for fully-implicit solve.
    std::shared_ptr<ImplicitHelper> m_implHelper;

    /// Protected constructor. Since we use a factory pattern, objects should be
    /// constructed via the SolverUtils::EquationSystem factory.
    RogersRicci1D(const LibUtilities::SessionReaderSharedPtr &session,
                  const SpatialDomains::MeshGraphSharedPtr &graph);

    virtual void v_InitObject(bool DeclareField = true) override;

    void InitialiseNonlinSysSolver(void);

    void ExplicitTimeInt(
        const Array<OneD, const Array<OneD, NekDouble>> &inarray,
        Array<OneD, Array<OneD, NekDouble>> &outarray, const NekDouble time);
    void DoOdeProjection(
        const Array<OneD, const Array<OneD, NekDouble>> &inarray,
        Array<OneD, Array<OneD, NekDouble>> &outarray, const NekDouble time);
    void GetFluxVector(const Array<OneD, Array<OneD, NekDouble>> &physfield,
                       Array<OneD, Array<OneD, Array<OneD, NekDouble>>> &flux);

    Array<OneD, NekDouble> &GetNormalVelocity();

    Array<OneD, NekDouble> &GetOutflowTraceNormal();
    void SetBoundaryConditions(Array<OneD, Array<OneD, NekDouble>> &physarray,
                               NekDouble time);
    int m_npts;
    int m_ndims;

private:
    int outflow_dim;
    NekDouble n_star;
    NekDouble tau;
    NekDouble T;

    /// User defined boundary conditions
    std::vector<CustomBCsSharedPtr> custom_BCs;

    // Point types in trace arrays
    std::vector<TracePtType> trace_pt_types;

    // Point types in field arrays
    std::vector<TracePtType> pt_types;
};

} // namespace Nektar

#endif
