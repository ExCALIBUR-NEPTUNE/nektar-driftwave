#include <boost/core/ignore_unused.hpp>

#include "OutflowBCs.h"

namespace Nektar
{

std::string OutflowBCs::className =
    GetCustomBCsFactory().RegisterCreatorFunction(
        "OutflowBCs", OutflowBCs::create, "OutflowBCs.");

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

    // n, u  in boundary region element(s)
    Array<OneD, NekDouble> n_bdy_elt, u_bdy_elt;
    m_fields[this->n_idx]->ExtractPhysToBndElmt(
        m_bcRegion, physarray[this->n_idx], n_bdy_elt);
    m_fields[this->u_idx]->ExtractPhysToBndElmt(
        m_bcRegion, physarray[this->u_idx], u_bdy_elt);

    // flux = n*u on boundary points
    Array<OneD, NekDouble> bdy_n(this->npts_bdy), bdy_u(this->npts_bdy),
        bdy_flux(this->npts_bdy);
    m_fields[this->n_idx]->ExtractElmtToBndPhys(m_bcRegion, n_bdy_elt, bdy_n);
    m_fields[this->u_idx]->ExtractElmtToBndPhys(m_bcRegion, u_bdy_elt, bdy_u);
    Vmath::Vmul(this->npts_bdy, bdy_n, 1, bdy_n, 1, bdy_flux, 1);

    // Update boundary expansion coeffs
    auto bdy_coeffs = this->bdy_explists[this->n_idx]->UpdateCoeffs();
    // segfaults
    // this->bdy_explists[this->n_idx]->IProductWRTBase(bdy_flux, bdy_coeffs);
}
} // namespace Nektar
