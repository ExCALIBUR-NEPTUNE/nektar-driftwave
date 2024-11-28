#include <boost/core/ignore_unused.hpp>

#include "OutflowBCs.h"

namespace Nektar
{

std::string OutflowBCs::className =
    GetCustomBCsFactory().RegisterCreatorFunction("Outflow", OutflowBCs::create,
                                                  "OutflowBCs.");

OutflowBCs::OutflowBCs(const LU::SessionReaderSharedPtr &session,
                       const Array<OneD, MR::ExpListSharedPtr> &fields,
                       const Array<OneD, Array<OneD, NekDouble>> &traceNormals,
                       const int target_idx, const int spaceDim,
                       const int bcRegion, const int cnt,
                       SD::BoundaryConditionShPtr cnd)
    : CustomBCs(session, fields, traceNormals, target_idx, spaceDim, bcRegion,
                cnt, cnd)
{
    this->n_idx = this->field_indices["n"];
    this->u_idx = this->field_indices["u"];
}

void OutflowBCs::v_Apply(Array<OneD, Array<OneD, NekDouble>> &physarray,
                         const NekDouble &time)
{

}
} // namespace Nektar
