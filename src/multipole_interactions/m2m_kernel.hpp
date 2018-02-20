#pragma once

#include "../common_kernel/interaction_constants.hpp"
#include "../common_kernel/interactions_iterators.hpp"
#include "../common_kernel/kernel_simd_types.hpp"
#include "../common_kernel/multiindex.hpp"
#include "../common_kernel/struct_of_array_data.hpp"
#include "taylor.hpp"

#include <vector>

namespace octotiger {
namespace fmm {
    namespace multipole_interactions {
    // Only counts in innermost loop!! Should be close enough though
    // 1+1+13+55+4+8+8+8+12+12+12+12+12+12+30+6+12+12+12+10+20
    constexpr size_t non_rho_vc_operations = 272;
    constexpr size_t non_rho_flop = non_rho_vc_operations * m2m_vector::size();
    // 1+4+12+9+6+11+3+11+2+11+2+11+12+8+11+2+11+1+11+7+11+3
    constexpr size_t only_rho_vc_operations = 160;
    constexpr size_t only_rho_flop = only_rho_vc_operations * m2m_vector::size();
    constexpr size_t rho_flop = only_rho_flop + non_rho_flop;

        class m2m_kernel
        {
        private:
            // std::vector<expansion>& local_expansions;
            // struct_of_array_data<expansion, real, 20, ENTRIES, SOA_PADDING>&
            // local_expansions_SoA;

            // com0 = *(com_ptr[0])
            // std::vector<space_vector>& center_of_masses;
            // down from 4 components
            // struct_of_array_data<space_vector, real, 3, ENTRIES, SOA_PADDING>&
            // center_of_masses_SoA;

            // // multipole expansion on this cell (L)
            // // std::vector<expansion>& potential_expansions;
            // struct_of_array_data<expansion, real, 20, ENTRIES, SOA_PADDING>&
            // potential_expansions_SoA;
            // // angular momentum correction on this cell (L_c)
            // // std::vector<space_vector>& angular_corrections;
            // struct_of_array_data<space_vector, real, 3, ENTRIES, SOA_PADDING>&
            // angular_corrections_SoA;

            std::vector<bool>& neighbor_empty;

            // so skip non-existing interaction partners faster, one entry per vector variable
            std::vector<bool> vector_is_empty;

            // std::vector<int32_t> coarse_x;
            // std::vector<int32_t> coarse_y;
            // std::vector<int32_t> coarse_z;

            const m2m_vector theta_rec_squared;
            m2m_int_vector offset_vector;

            // const double theta_rec_squared_scalar;

            // // for superimposed stencil
            // template <typename F>
            // void iterate_inner_cells_padded_stencil(multiindex<>& stencil_element, F& f) {
            //     for (size_t i0 = 0; i0 < INNER_CELLS_PER_DIRECTION; i0++) {
            //         for (size_t i1 = 0; i1 < INNER_CELLS_PER_DIRECTION; i1++) {
            //             // for (size_t i2 = 0; i2 < INNER_CELLS_PER_DIRECTION; i2++) {
            //             for (size_t i2 = 0; i2 < INNER_CELLS_PER_DIRECTION; i2 +=
            //             m2m_vector::Size)
            //             {
            //                 const multiindex<> cell_index(i0 + INNER_CELLS_PADDING_DEPTH,
            //                     i1 + INNER_CELLS_PADDING_DEPTH, i2 + INNER_CELLS_PADDING_DEPTH);
            //                 // BUG: indexing has to be done with uint32_t because of Vc
            //                 limitation
            //                 const int64_t cell_flat_index = to_flat_index_padded(cell_index);
            //                 const multiindex<> cell_index_unpadded(i0, i1, i2);
            //                 const int64_t cell_flat_index_unpadded =
            //                     to_inner_flat_index_not_padded(cell_index_unpadded);
            //                 const multiindex<> interaction_partner_index(
            //                     cell_index.x + stencil_element.x, cell_index.y +
            //                     stencil_element.y,
            //                     cell_index.z + stencil_element.z);

            //                 const int64_t interaction_flat_partner_index =
            //                     to_flat_index_padded(interaction_partner_index);
            //                 // std::cout << "cur: " << cell_index
            //                 //           << " partner: " << interaction_partner_index <<
            //                 std::endl;
            //                 f(cell_index, cell_flat_index, cell_index_unpadded,
            //                     cell_flat_index_unpadded, interaction_partner_index,
            //                     interaction_flat_partner_index);
            //             }
            //         }
            //     }
            // }

            // template <typename F>
            // void iterate_inner_cells_padded_stridded_stencil(
            //     multiindex<>& offset, multiindex<>& stencil_element, F& f) {
            //     for (size_t i0 = offset.x; i0 < INNER_CELLS_PER_DIRECTION; i0 += 2) {
            //         for (size_t i1 = offset.y; i1 < INNER_CELLS_PER_DIRECTION; i1 += 2) {
            //             for (size_t i2 = offset.z; i2 < INNER_CELLS_PER_DIRECTION; i2 += 2) {
            //                 // for (size_t i0 = 0; i0 < INNER_CELLS_PER_DIRECTION; i0 += 1) {
            //                 //     for (size_t i1 = 0; i1 < INNER_CELLS_PER_DIRECTION; i1 += 1) {
            //                 //         for (size_t i2 = 0; i2 < INNER_CELLS_PER_DIRECTION; i2 +=
            //                 1) {

            //                 const multiindex<> cell_index(i0 + INNER_CELLS_PADDING_DEPTH,
            //                     i1 + INNER_CELLS_PADDING_DEPTH, i2 + INNER_CELLS_PADDING_DEPTH);
            //                 const size_t cell_flat_index = to_flat_index_padded(cell_index);
            //                 const multiindex<> cell_index_unpadded(i0, i1, i2);
            //                 const size_t cell_flat_index_unpadded =
            //                     to_inner_flat_index_not_padded(cell_index_unpadded);
            //                 const multiindex<> interaction_partner_index(
            //                     cell_index.x + stencil_element.x, cell_index.y +
            //                     stencil_element.y,
            //                     cell_index.z + stencil_element.z);
            //                 const size_t interaction_partner_flat_index =
            //                     to_flat_index_padded(interaction_partner_index);

            //                 // std::cout << "cur: " << cell_index
            //                 //           << " partner: " << interaction_partner_index <<
            //                 std::endl;
            //                 f(cell_index, cell_flat_index, cell_index_unpadded,
            //                     cell_flat_index_unpadded, interaction_partner_index,
            //                     interaction_partner_flat_index);
            //             }
            //         }
            //     }
            // }

            void blocked_interaction_rho(struct_of_array_data<expansion, real, 20, ENTRIES,
                                             SOA_PADDING>& local_expansions_SoA,
                struct_of_array_data<space_vector, real, 3, ENTRIES, SOA_PADDING>&
                    center_of_masses_SoA,
                struct_of_array_data<expansion, real, 20, INNER_CELLS, SOA_PADDING>&
                    potential_expansions_SoA,
                struct_of_array_data<space_vector, real, 3, INNER_CELLS, SOA_PADDING>&
                    angular_corrections_SoA,
                               std::vector<real> &mons,
                const multiindex<>& cell_index, const size_t cell_flat_index,
                const multiindex<m2m_int_vector>& cell_index_coarse,
                const multiindex<>& cell_index_unpadded, const size_t cell_flat_index_unpadded,
                const two_phase_stencil& stencil, const size_t outer_stencil_index);

            void blocked_interaction_non_rho(struct_of_array_data<expansion, real, 20, ENTRIES,
                                                 SOA_PADDING>& local_expansions_SoA,
                struct_of_array_data<space_vector, real, 3, ENTRIES, SOA_PADDING>&
                    center_of_masses_SoA,
                struct_of_array_data<expansion, real, 20, INNER_CELLS, SOA_PADDING>&
                    potential_expansions_SoA,
                struct_of_array_data<space_vector, real, 3, INNER_CELLS, SOA_PADDING>&
                    angular_corrections_SoA,
                               std::vector<real> &mons,
                const multiindex<>& cell_index, const size_t cell_flat_index,
                const multiindex<m2m_int_vector>& cell_index_coarse,
                const multiindex<>& cell_index_unpadded, const size_t cell_flat_index_unpadded,
                const two_phase_stencil& stencil, const size_t outer_stencil_index);

            void vectors_check_empty();

            // void calculate_coarse_indices();

        public:
            static size_t multipole_callcounter_rho;
            static size_t multipole_callcounter_non_rho;
            m2m_kernel(
                //            struct_of_array_data<expansion, real, 20, ENTRIES, SOA_PADDING>&
                //            local_expansions_SoA,
                // struct_of_array_data<space_vector, real, 3, ENTRIES, SOA_PADDING>&
                // center_of_masses_SoA,
                // struct_of_array_data<expansion, real, 20, ENTRIES, SOA_PADDING>&
                //     potential_expansions_SoA,
                // struct_of_array_data<space_vector, real, 3, ENTRIES, SOA_PADDING>&
                //     angular_corrections_SoA,
                std::vector<bool>& neighbor_empty);

            m2m_kernel(m2m_kernel& other) = delete;

            m2m_kernel(const m2m_kernel& other) = delete;

            m2m_kernel operator=(const m2m_kernel& other) = delete;

            void apply_stencil(struct_of_array_data<expansion, real, 20, ENTRIES, SOA_PADDING>&
                                   local_expansions_SoA,
                struct_of_array_data<space_vector, real, 3, ENTRIES, SOA_PADDING>&
                    center_of_masses_SoA,
                struct_of_array_data<expansion, real, 20, INNER_CELLS, SOA_PADDING>&
                    potential_expansions_SoA,
                struct_of_array_data<space_vector, real, 3, INNER_CELLS, SOA_PADDING>&
                    angular_corrections_SoA,
                               std::vector<real> &mons,
                               const two_phase_stencil& stencil, gsolve_type type);
        };

    }    // namespace multipole_interactions
}    // namespace fmm
}    // namespace octotiger
