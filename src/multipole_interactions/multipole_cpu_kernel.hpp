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
        /** Controls the order in which the cpu multipole FMM interactions are calculated
         * (blocking). The actual numeric operations are found in
         * compute_kernel_templates.hpp. This
         * class is mostly responsible for loading data and control the order to increase
         * cache
         * efficieny on the cpu
         */
        class multipole_cpu_kernel
        {
        private:
            const m2m_vector theta_rec_squared;
            m2m_int_vector offset_vector;

            /// Executes a small block of RHO interactions (size is controlled by
            /// STENCIL_BLOCKING)
            void blocked_interaction_rho(const struct_of_array_data<expansion, real, 20, ENTRIES,
                                             SOA_PADDING>& local_expansions_SoA,
                const struct_of_array_data<space_vector, real, 3, ENTRIES, SOA_PADDING>&
                    center_of_masses_SoA,
                struct_of_array_data<expansion, real, 20, INNER_CELLS, SOA_PADDING>&
                    potential_expansions_SoA,
                struct_of_array_data<space_vector, real, 3, INNER_CELLS, SOA_PADDING>&
                    angular_corrections_SoA,
                const std::vector<real>& mons, const multiindex<>& cell_index,
                const size_t cell_flat_index, const multiindex<m2m_int_vector>& cell_index_coarse,
                const multiindex<>& cell_index_unpadded, const size_t cell_flat_index_unpadded,
                const two_phase_stencil& stencil, const size_t outer_stencil_index);

            /// Executes a small block of non-RHO interactions (size is controlled by
            /// STENCIL_BLOCKING)
            void blocked_interaction_non_rho(const struct_of_array_data<expansion, real, 20,
                                                 ENTRIES, SOA_PADDING>& local_expansions_SoA,
                const struct_of_array_data<space_vector, real, 3, ENTRIES, SOA_PADDING>&
                    center_of_masses_SoA,
                struct_of_array_data<expansion, real, 20, INNER_CELLS, SOA_PADDING>&
                    potential_expansions_SoA,
                struct_of_array_data<space_vector, real, 3, INNER_CELLS, SOA_PADDING>&
                    angular_corrections_SoA,
                const std::vector<real>& mons, const multiindex<>& cell_index,
                const size_t cell_flat_index, const multiindex<m2m_int_vector>& cell_index_coarse,
                const multiindex<>& cell_index_unpadded, const size_t cell_flat_index_unpadded,
                const two_phase_stencil& stencil, const size_t outer_stencil_index);

        public:
            multipole_cpu_kernel(void);

            multipole_cpu_kernel(multipole_cpu_kernel& other) = delete;

            multipole_cpu_kernel(const multipole_cpu_kernel& other) = delete;

            multipole_cpu_kernel operator=(const multipole_cpu_kernel& other) = delete;

            /// Calculate all multipole interactions for this kernel (runs the kernel)
            void apply_stencil(const struct_of_array_data<expansion, real, 20, ENTRIES,
                                   SOA_PADDING>& local_expansions_SoA,
                const struct_of_array_data<space_vector, real, 3, ENTRIES, SOA_PADDING>&
                    center_of_masses_SoA,
                struct_of_array_data<expansion, real, 20, INNER_CELLS, SOA_PADDING>&
                    potential_expansions_SoA,
                struct_of_array_data<space_vector, real, 3, INNER_CELLS, SOA_PADDING>&
                    angular_corrections_SoA,
                const std::vector<real>& mons, const two_phase_stencil& stencil, gsolve_type type);
        };

        // Helper classes for simplied access:
        template <size_t startindex, typename mask_t, typename mons_t, typename out_t>
        struct update_mask
        {
            update_mask(mask_t& mask_orig, const mask_t& mask_phase, const mons_t& mons,
                out_t& m_partner, const size_t interaction_partner_flat_index) {
                mask_orig = mask_orig & mask_phase;
            }
        };

        template <typename mask_t, typename mons_t, typename out_t>
        struct update_mask<0, mask_t, mons_t, out_t>
        {
            update_mask(mask_t& mask_orig, const mask_t& mask_phase, const mons_t& mons,
                out_t& m_partner, const size_t interaction_partner_flat_index) {
                Vc::where(mask_orig, m_partner[0]) = m2m_vector(
                    mons.data() + interaction_partner_flat_index, Vc::flags::element_aligned);
                mask_orig = mask_orig & mask_phase;
            }
        };

        template <size_t startindex, size_t endindex, typename in_t, typename out_t,
            typename mask_t>
        struct unrolled_SoA_load : unrolled_SoA_load<startindex, endindex - 1, in_t, out_t, mask_t>
        {
            unrolled_SoA_load(const in_t& local_expansions_SoA, out_t& m_partner,
                const mask_t& mask, const size_t interaction_partner_flat_index)
              : unrolled_SoA_load<startindex, endindex - 1, in_t, out_t, mask_t>(
                    local_expansions_SoA, m_partner, mask, interaction_partner_flat_index) {
                Vc::where(mask, m_partner[endindex]) =
                    local_expansions_SoA.template value<endindex>(interaction_partner_flat_index);
            }
        };

        template <typename in_t, typename out_t, typename mask_t>
        struct unrolled_SoA_load<0, 0, in_t, out_t, mask_t>
        {
            unrolled_SoA_load(const in_t& local_expansions_SoA, out_t& m_partner,
                const mask_t& mask, const size_t interaction_partner_flat_index) {
                Vc::where(mask, m_partner[0]) =
                    local_expansions_SoA.template value<0>(interaction_partner_flat_index);
            }
        };

    }    // namespace multipole_interactions
}    // namespace fmm
}    // namespace octotiger
