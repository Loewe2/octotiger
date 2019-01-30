#include "octotiger/monopole_interactions/p2p_interaction_interface.hpp"
#include "octotiger/monopole_interactions/calculate_stencil.hpp"

#include "octotiger/common_kernel/interactions_iterators.hpp"

#include "octotiger/options.hpp"

#include <algorithm>
#include <array>
#include <octotiger/debug_vector.hpp>

namespace octotiger {
namespace fmm {
    namespace monopole_interactions {
        const thread_local oct::vector<multiindex<>> p2p_interaction_interface::stencil =
            calculate_stencil().first;
        const thread_local oct::vector<bool> p2p_interaction_interface::stencil_masks =
            calculate_stencil_masks(p2p_interaction_interface::stencil).first;
        const thread_local oct::vector<oct::array<real, 4>> p2p_interaction_interface::four =
            calculate_stencil().second;
        const thread_local oct::vector<oct::array<real, 4>> p2p_interaction_interface::stencil_four_constants =
            calculate_stencil_masks(p2p_interaction_interface::stencil).second;
        thread_local oct::vector<real> p2p_interaction_interface::local_monopoles_staging_area(
            ENTRIES);

        p2p_interaction_interface::p2p_interaction_interface(void)
          : neighbor_empty_monopoles(27)
          , kernel_monopoles(neighbor_empty_monopoles) {
            this->p2p_type = opts().p2p_kernel_type;
        }

        void p2p_interaction_interface::compute_p2p_interactions(oct::vector<real>& monopoles,
            oct::vector<neighbor_gravity_type>& neighbors, gsolve_type type, real dx,
            oct::array<bool, geo::direction::count()>& is_direction_empty) {
            update_input(monopoles, neighbors, type, local_monopoles_staging_area);
            compute_interactions(type, is_direction_empty, neighbors, dx);
        }

        void p2p_interaction_interface::compute_interactions(gsolve_type type,
            oct::array<bool, geo::direction::count()>& is_direction_empty,
            oct::vector<neighbor_gravity_type>& all_neighbor_interaction_data, real dx) {
            if (p2p_type == interaction_kernel_type::SOA_CPU) {
                struct_of_array_data<expansion, real, 20, INNER_CELLS, SOA_PADDING>
                    potential_expansions_SoA;
                kernel_monopoles.apply_stencil_non_blocked(
                    local_monopoles_staging_area, potential_expansions_SoA, stencil_masks, stencil_four_constants, dx);
                potential_expansions_SoA.to_non_SoA(grid_ptr->get_L());
            } else {
                grid_ptr->compute_interactions(type);
                // waits for boundary data and then computes boundary interactions
                for (auto const& dir : geo::direction::full_set()) {
                    if (!is_direction_empty[dir]) {
                        neighbor_gravity_type& neighbor_data = all_neighbor_interaction_data[dir];
                        if (neighbor_data.is_monopole) {
                            grid_ptr->compute_boundary_interactions(type, neighbor_data.direction,
                                neighbor_data.is_monopole, neighbor_data.data);
                        }
                    }
                }
            }
        }
    }    // namespace monopole_interactions
}    // namespace fmm
}    // namespace octotiger
