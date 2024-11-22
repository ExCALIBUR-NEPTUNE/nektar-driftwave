#ifndef DRIFTWAVE_OUTFLOWBCS_H
#define DRIFTWAVE_OUTFLOWBCS_H

#include "CustomBCs.h"

namespace LR = Nektar::LocalRegions;
namespace LU = Nektar::LibUtilities;
namespace MR = Nektar::MultiRegions;
namespace SD = Nektar::SpatialDomains;
namespace Nektar
{

/**
 * @brief Sheath boundary conditions class
 */
class OutflowBCs : public CustomBCs
{
public:
    friend class MemoryManager<OutflowBCs>;

    /// Creates an instance of this class
    static CustomBCsSharedPtr create(
        const LU::SessionReaderSharedPtr &pSession,
        const Array<OneD, MR::ExpListSharedPtr> &pFields,
        const Array<OneD, Array<OneD, NekDouble>> &pTraceNormals,
        const int field_idx, const int pSpaceDim, const int bcRegion,
        const int cnt, SD::BoundaryConditionShPtr cnd)
    {
        CustomBCsSharedPtr p = MemoryManager<OutflowBCs>::AllocateSharedPtr(
            pSession, pFields, pTraceNormals, field_idx, pSpaceDim, bcRegion,
            cnt, cnd);
        return p;
    }

    /// Name of the class
    static std::string className;

protected:
    virtual void v_Apply(Array<OneD, Array<OneD, NekDouble>> &physarray,
                         const NekDouble &time);

private:
    OutflowBCs(const LU::SessionReaderSharedPtr &session,
               const Array<OneD, MR::ExpListSharedPtr> &fields,
               const Array<OneD, Array<OneD, NekDouble>> &traceNormals,
               const int target_idx, const int spaceDim, const int bcRegion,
               const int cnt, SD::BoundaryConditionShPtr cnd);

    virtual ~OutflowBCs(void){};

    int n_idx;
    int u_idx;
};

} // namespace Nektar

#endif
