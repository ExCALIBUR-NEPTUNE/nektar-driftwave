#include "CustomBCs.h"

namespace Nektar
{
// Static factory singleton
CustomBCsFactory &GetCustomBCsFactory()
{
    static CustomBCsFactory instance;
    return instance;
}

/**
 * Abstract base class for custom boundary conditions.
 */
CustomBCs::CustomBCs(const LU::SessionReaderSharedPtr &session,
                     const Array<OneD, MR::ExpListSharedPtr> &fields,
                     const Array<OneD, Array<OneD, NekDouble>> &trace_normals,
                     const int target_idx, const int space_dim,
                     const int bc_region, const int cnt,
                     SD::BoundaryConditionShPtr cnd)
    : m_session(session), m_fields(fields), m_traceNormals(trace_normals),
      target_idx(target_idx), m_spacedim(space_dim), m_bcRegion(bc_region),
      m_offset(cnt), m_cnd(cnd), field_indices()
{

    // Variable name => index map
    auto vars = session->GetVariables();
    for (auto ii = 0; ii < vars.size(); ii++)
    {
        field_indices[vars[ii]] = ii;
    }

    this->bdy_explists =
        Array<OneD, MultiRegions::ExpListSharedPtr>(m_fields.size());
    for (auto ii = 0; ii < m_fields.size(); ++ii)
    {
        this->bdy_explists[ii] =
            m_fields[0]->GetBndCondExpansions()[m_bcRegion];
    }

    // Number of points on the boundary in this region
    this->npts_bdy =
        m_fields[0]->GetBndCondExpansions()[m_bcRegion]->GetTotPoints();
}

/**
 * @param   bcRegion      id of the boundary region
 * @param   cnt
 * @param   Fwd
 * @param   physarray
 * @param   time
 */
void CustomBCs::Apply(Array<OneD, Array<OneD, NekDouble>> &physarray,
                      const NekDouble &time)
{
    v_Apply(physarray, time);
}

/**
 * @ brief Newly added bc should specify this virtual function
 * if the Bwd/value in m_bndCondExpansions is the target value, like Dirichlet,
 * bc weight should be 1.0.
 * if some average Fwd and Bwd/value in m_bndCondExpansions
 * is the target value like WallViscousBC weight should be 0.5.
 */
void CustomBCs::v_ApplyBwdWeight()
{
    for (int i = 0; i < m_fields.size(); ++i)
    {
        m_fields[i]->SetBndCondBwdWeight(m_bcRegion, m_weight);
    }
}

} // namespace Nektar
