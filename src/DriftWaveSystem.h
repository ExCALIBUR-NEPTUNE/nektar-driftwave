///////////////////////////////////////////////////////////////////////////////
//
// File DriftWaveSystem.h
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

#ifndef NEKTAR_DRIFTWAVESYSTEM_H
#define NEKTAR_DRIFTWAVESYSTEM_H

#include <SolverUtils/AdvectionSystem.h>
#include <SolverUtils/RiemannSolvers/RiemannSolver.h>

#include "ImplicitHelper.hpp"

using namespace Nektar::SolverUtils;

namespace Nektar
{

/**
 * @brief An equation system for the drift-wave solver.
 */
class DriftWaveSystem : public AdvectionSystem
{
public:
    // Friend class to allow the memory manager to allocate shared pointers of
    // this class.
    friend class MemoryManager<DriftWaveSystem>;

    /// Creates an instance of this class. This static method is registered with
    /// a factory.
    static SolverUtils::EquationSystemSharedPtr create(
        const LibUtilities::SessionReaderSharedPtr &session,
        const SpatialDomains::MeshGraphSharedPtr &graph)
    {
        SolverUtils::EquationSystemSharedPtr p =
            MemoryManager<DriftWaveSystem>::AllocateSharedPtr(session, graph);
        p->InitObject();
        return p;
    }

    /// Name of class, used to statically initialise a function pointer for the
    /// create method above.
    static std::string className;

    /// Default destructor.
    virtual ~DriftWaveSystem() = default;

protected:
    /// Storage for the drift velocity. The outer index is dimension, and inner
    /// index the solution nodes (in physical space).
    Array<OneD, Array<OneD, NekDouble>> m_driftVel;
    /// Storage for the dot product of drift velocity with element edge normals,
    /// required for the DG formulation.
    Array<OneD, NekDouble> m_traceVn;
    /// \f$ \alpha \f$ constant for the Hasegawa-Wakatani equations.
    NekDouble m_alpha;
    /// \f$ \kappa \f$ constant for the Hasegawa-Wakatani equations.
    NekDouble m_kappa;
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
    DriftWaveSystem(const LibUtilities::SessionReaderSharedPtr &session,
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
};

} // namespace Nektar

#endif
