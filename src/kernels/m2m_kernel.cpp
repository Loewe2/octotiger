#include "m2m_kernel.hpp"

#include "defs.hpp"
#include "grid_flattened_indices.hpp"
#include "interaction_types.hpp"

#include <array>
#include <functional>

extern taylor<4, real> factor;

std::vector<interaction_type> ilist_debugging;

namespace octotiger {
namespace fmm {

    void m2m_kernel::operator()(const multiindex& cell_index, const size_t cell_flat_index,
        const multiindex& cell_index_unpadded, const size_t cell_flat_index_unpadded,
        const multiindex& interaction_partner_index, const size_t interaction_partner_flat_index) {
        // if (cell_flat_index_unpadded != 0) {
        //     return;
        // }
        // if (interaction_partner_index.x - cell_index.x < 0 ||
        //     interaction_partner_index.y - cell_index.y < 0 ||
        //     interaction_partner_index.z - cell_index.z < 0) {
        //     return;
        // }

        std::array<real, NDIM>
            X;    // TODO: replace by space_vector for vectorization or get rid of temporary
        for (integer d = 0; d < NDIM; ++d) {
            // TODO: need SoA for this access if vectorized
            X[d] = center_of_masses[cell_flat_index][d];
        }

        // interaction_type current_interaction;
        // current_interaction.first =
        //     cell_flat_index;    // has to be translated to unpadded before comparison
        // current_interaction.second =
        //     interaction_partner_flat_index;    // has to be translated to unpadded before
        //     comparison
        // current_interaction.four = {0};
        // current_interaction.x[0] = (interaction_partner_index.x - cell_index.x);
        // current_interaction.x[1] = (interaction_partner_index.y - cell_index.y);
        // current_interaction.x[2] = (interaction_partner_index.z - cell_index.z);
        // current_interaction.first_index = {
        //     {cell_index_unpadded.x, cell_index_unpadded.y, cell_index_unpadded.z}};
        // current_interaction.second_index = {{interaction_partner_index.x - 8,
        //     interaction_partner_index.y - 8, interaction_partner_index.z - 8}};

        // ilist_debugging.push_back(current_interaction);

        std::array<real, NDIM>
            Y;    // TODO: replace by space_vector for vectorization or get rid of temporary
        for (integer d = 0; d < NDIM; ++d) {
            Y[d] = center_of_masses[interaction_partner_flat_index][d];
        }

        // cell specific taylor series coefficients
        // multipole const& Miii1 = M_ptr[iii1];
        // taylor<4, real> m0;    // TODO: replace by expansion type or get rid of temporary
        //                        // (replace by simd_vector for vectorization)
        // for (integer j = 0; j != taylor_sizes[3]; ++j) {
        //     m0[j] = local_expansions[interaction_partner_flat_index][j];
        // }
        expansion& m_partner = local_expansions[interaction_partner_flat_index];

        // n angular momentum of the cells
        // TODO: replace by expansion type or get rid of temporary
        // (replace by simd_vector for vectorization)
        taylor<4, real> n0;

        // Initalize moments and momentum
        if (type != RHO) {
            for (integer j = taylor_sizes[2]; j != taylor_sizes[3]; ++j) {
                n0[j] = ZERO;
            }
        } else {
            // this branch computes the angular momentum correction, (20) in the
            // paper divide by mass of other cell
            expansion& m_cell = local_expansions[cell_flat_index];

            real const tmp1 = m_partner() / m_cell();
            // calculating the coefficients for formula (M are the octopole moments)
            // the coefficients are calculated in (17) and (18)
            for (integer j = taylor_sizes[2]; j != taylor_sizes[3]; ++j) {
                n0[j] = m_partner[j] - m_cell[j] * tmp1;
            }
        }

        // taylor<4, simd_vector> A0;
        // std::array<simd_vector, NDIM> B0 = {
        //     {simd_vector(ZERO), simd_vector(ZERO), simd_vector(ZERO)}};

        // distance between cells in all dimensions
        // TODO: replace by simd_vector for vectorization or get rid of temporary
        std::array<real, NDIM> dX;
        for (integer d = 0; d < NDIM; ++d) {
            dX[d] = X[d] - Y[d];
        }

        // R_i in paper is the dX in the code
        // D is taylor expansion value for a given X expansion of the gravitational potential
        // (multipole expansion)
        taylor<5, real> D;    // TODO: replace by simd_vector for vectorization

        // std::cout << "dX: ";
        // for (size_t i = 0; i < dX.size(); i++) {
        //     if (i > 0) {
        //         std::cout << ", ";
        //     }
        //     std::cout << dX[i];
        // }
        // std::cout << std::endl;

        // calculates all D-values, calculate all coefficients of 1/r (not the potential),
        // formula (6)-(9) and (19)
        D.set_basis(dX);    // TODO: after vectorization, make sure the vectorized version is called

        // std::cout << "D: ";
        // for (size_t i = 0; i < D.size(); i++) {
        //     if (i > 0) {
        //         std::cout << ", ";
        //     }
        //     std::cout << D[i];
        // }
        // std::cout << std::endl;

        // output variable references
        // TODO: use these again after debugging individual components!
        // expansion& current_potential = potential_expansions[cell_flat_index_unpadded];
        // space_vector& current_angular_correction = angular_corrections[cell_flat_index_unpadded];

        expansion current_potential;
        for (size_t i = 0; i < current_potential.size(); i++) {
            current_potential[i] = 0.0;
        }
        space_vector current_angular_correction;
        for (size_t i = 0; i < current_angular_correction.size(); i++) {
            current_angular_correction[i] = 0.0;
        }

        // the following loops calculate formula (10), potential from B->A
        current_potential[0] += m_partner[0] * D[0];
        if (type != RHO) {
            for (integer i = taylor_sizes[0]; i != taylor_sizes[1]; ++i) {
                current_potential[0] -= m_partner[i] * D[i];
            }
        }
        for (integer i = taylor_sizes[1]; i != taylor_sizes[2]; ++i) {
            const auto tmp = D[i] * (factor[i] * HALF);
            current_potential[0] += m_partner[i] * tmp;
        }
        for (integer i = taylor_sizes[2]; i != taylor_sizes[3]; ++i) {
            const auto tmp = D[i] * (factor[i] * SIXTH);
            current_potential[0] -= m_partner[i] * tmp;
        }

        for (integer a = 0; a < NDIM; ++a) {
            int const* ab_idx_map = to_ab_idx_map3[a];
            int const* abc_idx_map = to_abc_idx_map3[a];

            current_potential(a) += m_partner() * D(a);
            for (integer i = 0; i != 6; ++i) {
                if (type != RHO && i < 3) {
                    current_potential(a) -= m_partner(a) * D[ab_idx_map[i]];
                }
                const integer cb_idx = cb_idx_map[i];
                const auto tmp1 = D[abc_idx_map[i]] * (factor[cb_idx] * HALF);
                current_potential(a) += m_partner[cb_idx] * tmp1;
            }
        }

        if (type == RHO) {
            for (integer a = 0; a != NDIM; ++a) {
                int const* abcd_idx_map = to_abcd_idx_map3[a];
                for (integer i = 0; i != 10; ++i) {
                    const integer bcd_idx = bcd_idx_map[i];
                    const auto tmp = D[abcd_idx_map[i]] * (factor[bcd_idx] * SIXTH);
                    current_angular_correction[a] -= n0[bcd_idx] * tmp;
                }
            }
        }

        for (integer i = 0; i != 6; ++i) {
            int const* abc_idx_map6 = to_abc_idx_map6[i];

            integer const ab_idx = ab_idx_map6[i];
            current_potential[ab_idx] += m_partner() * D[ab_idx];
            for (integer c = 0; c < NDIM; ++c) {
                const auto& tmp = D[abc_idx_map6[c]];
                current_potential[ab_idx] -= m_partner(c) * tmp;
            }
        }

        for (integer i = taylor_sizes[2]; i != taylor_sizes[3]; ++i) {
            const auto& tmp = D[i];
            current_potential[i] += m_partner[0] * tmp;
        }

        // TODO: remove this when switching back to non-copy (reference-based) approach
        expansion& current_potential_result = potential_expansions[cell_flat_index_unpadded];
        space_vector& current_angular_correction_result =
            angular_corrections[cell_flat_index_unpadded];
        for (size_t i = 0; i < current_potential.size(); i++) {
            current_potential_result[i] += current_potential[i];
        }
        for (size_t i = 0; i < current_angular_correction.size(); i++) {
            current_angular_correction_result[i] += current_angular_correction[i];
        }

        // if (cell_index_unpadded.x == 0 && cell_index_unpadded.y == 0 &&
        //     cell_index_unpadded.z == 0) {
        //     std::cout << "-------------------------------------------------------------"
        //               << std::endl;
        //     std::cout << "cell_index: " << cell_index;
        //     std::cout << " interaction_partner_index: " << interaction_partner_index;
        //     std::cout << " diff: " << (interaction_partner_index.x - cell_index.x) << ", "
        //               << (interaction_partner_index.y - cell_index.y) << ", "
        //               << (interaction_partner_index.z - cell_index.z) << std::endl;
        //     // std::cout << " cell_index_unpadded: " << cell_index_unpadded;
        //     std::cout << " cell_flat_index: " << cell_flat_index;
        //     std::cout << " cell_flat_index_unpadded: " << cell_flat_index_unpadded;
        //     std::cout << " interaction_partner_flat_index: " << interaction_partner_flat_index;
        //     std::cout << std::endl;
        //     if (interaction_partner_index.x < 8 || interaction_partner_index.y < 8 ||
        //         interaction_partner_index.z < 8 || interaction_partner_index.x >= 16 ||
        //         interaction_partner_index.y >= 16 || interaction_partner_index.z >= 16) {
        //         multiindex old_partner_index(interaction_partner_index.x % 8,
        //             interaction_partner_index.y % 8, interaction_partner_index.z % 8);
        //         std::cout << "old partner index: " << old_partner_index << std::endl;
        //         std::cout << "old partner flat index: "
        //                   << to_inner_flat_index_not_padded(old_partner_index) << std::endl;
        //     }
        //     std::cout << "X: ";
        //     for (size_t d = 0; d < NDIM; d++) {
        //         if (d > 0) {
        //             std::cout << ", ";
        //         }
        //         std::cout << X[d];
        //     }
        //     std::cout << std::endl;
        //     std::cout << "Y: ";
        //     for (size_t d = 0; d < NDIM; d++) {
        //         if (d > 0) {
        //             std::cout << ", ";
        //         }
        //         std::cout << Y[d];
        //     }
        //     std::cout << std::endl;
        //     std::cout << "m_partner: ";
        //     for (size_t i = 0; i < m_partner.size(); i++) {
        //         if (i > 0) {
        //             std::cout << ", ";
        //         }
        //         std::cout << m_partner[i];
        //     }
        //     std::cout << std::endl;
        //     std::cout << "D: ";
        //     for (size_t i = 0; i < D.size(); i++) {
        //         if (i > 0) {
        //             std::cout << ", ";
        //         }
        //         std::cout << D[i];
        //     }
        //     std::cout << std::endl;
        //     std::cout << "current_potential: ";
        //     for (size_t i = 0; i < current_potential.size(); i++) {
        //         if (i > 0) {
        //             std::cout << ", ";
        //         }
        //         std::cout << current_potential[i];
        //     }
        //     std::cout << std::endl;

        //     std::cout << "current_angular_correction: ";
        //     for (size_t d = 0; d < NDIM; d++) {
        //         if (d > 0) {
        //             std::cout << ", ";
        //         }
        //         std::cout << current_angular_correction[d];
        //     }
        //     std::cout << std::endl;
        // }
    }

    m2m_kernel::m2m_kernel(std::vector<expansion>& local_expansions,
        std::vector<space_vector>& center_of_masses, std::vector<expansion>& potential_expansions,
        std::vector<space_vector>& angular_corrections, gsolve_type type)
      : local_expansions(local_expansions)
      , center_of_masses(center_of_masses)
      , potential_expansions(potential_expansions)
      , angular_corrections(angular_corrections)
      , type(type) {}

    void m2m_kernel::apply_stencils(std::array<std::vector<multiindex>, 8>& stencils) {
        for (int64_t i0 = 0; i0 < 2; i0 += 1) {
            for (int64_t i1 = 0; i1 < 2; i1 += 1) {
                for (int64_t i2 = 0; i2 < 2; i2 += 1) {
                    multiindex offset(i0, i1, i2);
                    size_t stencil_index = i0 * 4 + i1 * 2 + i2;
                    // std::cout << "stencil_index: " << stencil_index << std::endl;
                    std::vector<multiindex>& stencil = stencils[stencil_index];
                    for (multiindex& stencil_element : stencil) {
                        iterate_inner_cells_padded_stridded_stencil(offset, stencil_element, *this);
                    }
                }
            }
        }
    }

    // void m2m_kernel::apply_stencil_element(multiindex& stencil_element) {
    //     // execute kernel for individual interactions for every inner cell
    //     // use this object as a functor for the iteration
    //     iterate_inner_cells_padded_stencil(stencil_element, *this);

    //     // std::cout << "potential_expansions.size(): " << potential_expansions.size() <<
    //     std::endl;
    //     // std::cout << "potential_expansions[0]: " << potential_expansions[0];
    //     // std::cout << std::endl;
    // }

}    // namespace fmm
}    // namespace octotiger
