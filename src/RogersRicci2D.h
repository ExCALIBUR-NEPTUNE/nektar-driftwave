///////////////////////////////////////////////////////////////////////////////
//
// File RogersRicci2D.h
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

#ifndef NEKTAR_ROGERSRICCI2D_H
#define NEKTAR_ROGERSRICCI2D_H

#include <SolverUtils/AdvectionSystem.h>
#include <SolverUtils/RiemannSolvers/RiemannSolver.h>

#include "ImplicitHelper.hpp"

using namespace Nektar::SolverUtils;

namespace Nektar
{

/**
 * @brief An equation system for the drift-wave solver.
 */
class RogersRicci2D : public AdvectionSystem
{
public:
    // Friend class to allow the memory manager to allocate shared pointers of
    // this class.
    friend class MemoryManager<RogersRicci2D>;

    /// Creates an instance of this class. This static method is registered with
    /// a factory.
    static SolverUtils::EquationSystemSharedPtr create(
        const LibUtilities::SessionReaderSharedPtr &session,
        const SpatialDomains::MeshGraphSharedPtr &graph)
    {
        SolverUtils::EquationSystemSharedPtr p =
            MemoryManager<RogersRicci2D>::AllocateSharedPtr(session, graph);
        p->InitObject();
        return p;
    }

    /// Name of class, used to statically initialise a function pointer for the
    /// create method above.
    static std::string className;

    // Require a fixed variable order; use these indices for clarity
    static constexpr int n_idx   = 0;
    static constexpr int Te_idx  = 1;
    static constexpr int w_idx   = 2;
    static constexpr int phi_idx = 3;

    /// Default destructor.
    virtual ~RogersRicci2D() = default;

protected:
    /// Helper function to define constants.
    NekDouble c(std::string n)
    {
        auto it = m_c.find(n);

        ASSERTL0(it != m_c.end(), "Unknown constant");

        return it->second;
    }

    /// Map of known constants
    std::map<std::string, NekDouble> m_c = {
        {"T_e", 6.0},      {"L_z", 18.0},       {"n_0", 2.0e18},
        {"m_i", 6.67e-27}, {"omega_ci", 9.6e5}, {"lambda", 3.0},
        {"R", 0.5}};

    /// Storage for the drift velocity. The outer index is dimension, and inner
    /// index the solution nodes (in physical space).
    Array<OneD, Array<OneD, NekDouble>> m_driftVel;
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

    Array<OneD, NekDouble> m_r;

    /// Protected constructor. Since we use a factory pattern, objects should be
    /// constructed via the SolverUtils::EquationSystem factory.
    RogersRicci2D(const LibUtilities::SessionReaderSharedPtr &session,
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

    int m_npts;
    int m_ndims;
};

} // namespace Nektar

#endif
