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
    namespace monopole_interactions {
        constexpr uint64_t P2M_STENCIL_BLOCKING = 1;
    //  3+13+55+4+12+12+12+12+12+12+12+12+6+8+4
    constexpr size_t non_rho_p2m_vc_operations = 189;
    constexpr size_t non_rho_p2m_flop = non_rho_p2m_vc_operations * m2m_vector::size();
    // 14+9+6+9+3+9+2+9+2+9+10+7+10+11+10+16+6+3
    constexpr size_t only_rho_p2m_vc_operations = 145;
    constexpr size_t only_rho_p2m_flop = only_rho_p2m_vc_operations * m2m_vector::size();
    constexpr size_t rho_p2m_flop = non_rho_p2m_flop + only_rho_p2m_flop;

        class p2m_kernel
        {
        private:
            /// Stores which neighbors are empty so that we can just skip them
            std::vector<bool>& neighbor_empty;

            /// so skip non-existing interaction partners faster, one entry per vector variable
            std::vector<bool> vector_is_empty;

            const m2m_vector theta_rec_squared;
            m2m_int_vector offset_vector;

            /// Calculates the monopole multipole boundary interactions with solve type rho
            void blocked_interaction_rho(struct_of_array_data<expansion, real, 20, ENTRIES,
                                             SOA_PADDING>& local_expansions_SoA,
                struct_of_array_data<space_vector, real, 3, ENTRIES, SOA_PADDING>&
                    center_of_masses_SoA,
                struct_of_array_data<expansion, real, 20, INNER_CELLS, SOA_PADDING>&
                    potential_expansions_SoA,
                struct_of_array_data<space_vector, real, 3, INNER_CELLS, SOA_PADDING>&
                    angular_corrections_SoA,
                const multiindex<>& cell_index, const size_t cell_flat_index,
                const multiindex<m2m_int_vector>& cell_index_coarse,
                const multiindex<>& cell_index_unpadded, const size_t cell_flat_index_unpadded,
                const multiindex<>& interaction_partner_index,
                const size_t interaction_partner_flat_index,
                multiindex<m2m_int_vector>& interaction_partner_index_coarse);

            /// Calculates the monopole multipole boundary interactions without the solve type rho
            void blocked_interaction_non_rho(struct_of_array_data<expansion, real, 20, ENTRIES,
                                                 SOA_PADDING>& local_expansions_SoA,
                struct_of_array_data<space_vector, real, 3, ENTRIES, SOA_PADDING>&
                    center_of_masses_SoA,
                struct_of_array_data<expansion, real, 20, INNER_CELLS, SOA_PADDING>&
                    potential_expansions_SoA,
                struct_of_array_data<space_vector, real, 3, INNER_CELLS, SOA_PADDING>&
                    angular_corrections_SoA,
                const multiindex<>& cell_index, const size_t cell_flat_index,
                const multiindex<m2m_int_vector>& cell_index_coarse,
                const multiindex<>& cell_index_unpadded, const size_t cell_flat_index_unpadded,
                const multiindex<>& interaction_partner_index,
                const size_t interaction_partner_flat_index,
                multiindex<m2m_int_vector>& interaction_partner_index_coarse);

            void vectors_check_empty();

        public:
            static size_t p2m_counter_rho;
            static size_t p2m_counter_non_rho;
            p2m_kernel(std::vector<bool>& neighbor_empty);

            p2m_kernel(p2m_kernel& other) = delete;
            p2m_kernel(const p2m_kernel& other) = delete;
            p2m_kernel operator=(const p2m_kernel& other) = delete;

            void apply_stencil(struct_of_array_data<expansion, real, 20, ENTRIES, SOA_PADDING>&
                                   local_expansions_SoA,
                struct_of_array_data<space_vector, real, 3, ENTRIES, SOA_PADDING>&
                    center_of_masses_SoA,
                struct_of_array_data<expansion, real, 20, INNER_CELLS, SOA_PADDING>&
                    potential_expansions_SoA,
                struct_of_array_data<space_vector, real, 3, INNER_CELLS, SOA_PADDING>&
                    angular_corrections_SoA,
                std::vector<multiindex<>>& stencil, gsolve_type type, bool (&z_skip)[3][3][3],
                bool (&y_skip)[3][3], bool (&x_skip)[3]);
        };

    }    // namespace monopole_interactions
}    // namespace fmm
}    // namespace octotiger
