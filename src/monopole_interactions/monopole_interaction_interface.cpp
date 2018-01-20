#include "monopole_interaction_interface.hpp"

#include "../common_kernel/interactions_iterators.hpp"
#include "calculate_stencil.hpp"
#include "p2m_kernel.hpp"
#include "p2p_kernel.hpp"

#include <algorithm>

// Big picture questions:
// - use any kind of tiling?

namespace octotiger {
namespace fmm {
    namespace monopole_interactions {
        size_t total_neighbors = 0;
        size_t missing_neighbors = 0;

        std::vector<multiindex<>> monopole_interaction_interface::stencil;
        std::vector<std::array<real, 4>> monopole_interaction_interface::four;

        monopole_interaction_interface::monopole_interaction_interface(void)
          : neighbor_empty_multipoles(27)
          , neighbor_empty_monopoles(27)
          , kernel(neighbor_empty_multipoles)
          , kernel_monopoles(neighbor_empty_monopoles) {
            // Create our input structure for the compute kernel
            for (auto i = 0; i < 27; ++i) {
                local_expansions_array[i] = std::vector<expansion>(EXPANSION_COUNT_NOT_PADDED);
                center_of_masses_array[i] = std::vector<space_vector>(EXPANSION_COUNT_NOT_PADDED);
            }

            potential_expansions = std::vector<expansion>(EXPANSION_COUNT_NOT_PADDED);
            angular_corrections = std::vector<space_vector>(EXPANSION_COUNT_NOT_PADDED);

            local_monopoles = std::vector<real>(EXPANSION_COUNT_PADDED);
        }

        void monopole_interaction_interface::update_input(std::vector<real>& mons,
            std::vector<multipole>& multipoles,
            std::vector<std::shared_ptr<std::vector<space_vector>>>& com_ptr,
            std::vector<neighbor_gravity_type>& neighbors, gsolve_type t, real dx_arg) {
            std::vector<space_vector> const& com0 = *(com_ptr[0]);
            type = t;
            dx = dx_arg;

            iterate_inner_cells_padded(
                [this, mons, multipoles, com0](const multiindex<>& i, const size_t flat_index,
                    const multiindex<>& i_unpadded, const size_t flat_index_unpadded) {
                    local_expansions_array[13].at(flat_index_unpadded) = 0.0;
                    center_of_masses_array[13].at(flat_index_unpadded) = com0.at(flat_index_unpadded);

                    local_monopoles.at(flat_index) = mons.at(flat_index_unpadded);
                });
            total_neighbors += 27;

            size_t current_missing = 0;
            size_t current_monopole = 0;

            for (size_t i = 0; i < neighbor_empty_multipoles.size(); i++) {
                neighbor_empty_multipoles[i] = false;
                neighbor_empty_monopoles[i] = false;
            }

            multipole_neighbors_exist = false;

            // Now look at neighboring data
            for (const geo::direction& dir : geo::direction::full_set()) {
                // don't use neighbor.direction, is always zero for empty cells!
                neighbor_gravity_type& neighbor = neighbors[dir];
                // Switch x and z dimension since the stencil code and the old octotiger code
                // a different order for some reason
                auto x = dir.operator[](ZDIM) + 1;
                auto y = dir.operator[](YDIM) + 1;
                auto z = dir.operator[](XDIM) + 1;
                // std::cout << dir.flat_index_with_center() << ":" << z << " " << y << " " << x <<
                // std::endl;
                // std::cout << dir.flat_index_with_center() << ":" << z + 1 << " " << y + 1 << " "
                // << x+1 << std::endl;
                // std::cout << dir.flat_index_with_center() << ":" << (z + 1) * 3 * 3 << " " << (y
                // + 1) * 3 << " " << (x+1) << std::endl;
                // std::cout << dir.flat_index_with_center() << ":" << (z + 1) * 3 * 3  + (y + 1) *
                // 3  + (x+1) << std::endl;
                // std::cin.get();
                auto flat = x * 3 * 3 + y * 3 + z;

                // this dir is setup as a multipole - and we only consider multipoles here
                if (!neighbor.is_monopole) {
                    neighbor_empty_monopoles[dir.flat_index_with_center()] = true;
                    // neighbor has no data - input structure just recevices zeros as padding
                    if (!neighbor.data.M) {
                        iterate_inner_cells_padding(
                            dir, [this, flat](const multiindex<>& i, const size_t flat_index,
                                     const multiindex<>&,
                                     const size_t) { local_monopoles.at(flat_index) = 0.0; });
                        neighbor_empty_multipoles[dir.flat_index_with_center()] = true;
                        z_skip[z][y][x] = true;
                    } else {
                        // Get multipole data into our input structure
                        std::vector<multipole>& neighbor_M_ptr = *(neighbor.data.M);
                        std::vector<space_vector>& neighbor_com0 = *(neighbor.data.x);
                        iterate_inner_cells_padding(
                            dir, [this, neighbor_M_ptr, neighbor_com0, flat](const multiindex<>& i,
                                     const size_t flat_index, const multiindex<>& i_unpadded,
                                     const size_t flat_index_unpadded) {

                                local_expansions_array[flat].at(flat_index_unpadded) =
                                    neighbor_M_ptr.at(flat_index_unpadded);
                                center_of_masses_array[flat].at(flat_index_unpadded) =
                                    neighbor_com0.at(flat_index_unpadded);

                                local_monopoles.at(flat_index) = 0.0;
                            });
                        multipole_neighbors_exist = true;
                        z_skip[z][y][x] = false;
                    }
                } else {
                    neighbor_empty_multipoles[dir.flat_index_with_center()] = true;
                    z_skip[z][y][x] = true;
                    if (neighbor.is_monopole) {
                        if (!neighbor.data.m) {
                            // TODO: ask Dominic why !is_monopole and stuff still empty
                            iterate_inner_cells_padding(
                                dir, [this](const multiindex<>& i, const size_t flat_index,
                                         const multiindex<>&, const size_t) {
                                    // initializes whole expansion, relatively expansion
                                    local_monopoles.at(flat_index) = 0.0;
                                });
                            neighbor_empty_monopoles[dir.flat_index_with_center()] = true;
                        } else {
                            std::vector<real>& neighbor_mons = *(neighbor.data.m);
                            iterate_inner_cells_padding(
                                dir, [this, neighbor_mons](const multiindex<>& i,
                                         const size_t flat_index, const multiindex<>&,
                                         const size_t flat_index_unpadded) {
                                    // initializes whole expansion, relatively expansion
                                    local_monopoles.at(flat_index) =
                                        neighbor_mons.at(flat_index_unpadded);
                                });
                        }
                    }
                }
            }

            neighbor_empty_multipoles[13] = true;
            neighbor_empty_monopoles[13] = false;
            center_of_masses_SoA_array[13].update_data(
            center_of_masses_array[13]);

            iterate_inner_cells_not_padded(
                [this](const multiindex<>& i_unpadded, const size_t flat_index_unpadded) {
                    expansion& e = potential_expansions.at(flat_index_unpadded);
                    e = 0.0;
                });
            iterate_inner_cells_not_padded(
                [this](const multiindex<>& i_unpadded, const size_t flat_index_unpadded) {
                    space_vector& s = angular_corrections.at(flat_index_unpadded);
                    s = 0.0;
                });

            z_skip[1][1][1] = true;
            for (auto xi = 0; xi < 3; ++xi) {
                x_skip[xi] = true;
                for (auto yi = 0; yi < 3; ++yi) {
                    y_skip[xi][yi] = true;
                    for (auto zi = 0; zi < 3; ++zi) {
                        if (!z_skip[xi][yi][zi]) {
                            y_skip[xi][yi] = false;
                            auto flat = xi * 3 * 3 + yi * 3 + zi;
                            local_expansions_SoA_array[flat].update_data(
                                local_expansions_array[flat]);
                            center_of_masses_SoA_array[flat].update_data(
                                center_of_masses_array[flat]);
                        }
                    }
                    if (!y_skip[xi][yi])
                        x_skip[xi] = false;
                }
            }

            // if (multipole_neighbors_exist) {
            // for (auto zi = 0; zi < 3; ++zi) {
            //     std::cout << z_skip[zi] << " ";
            // }
            // std::cout << "\n";
            // for (auto zi = 0; zi < 3; ++zi) {
            //     for (auto yi = 0; yi < 3; ++yi) {
            //         std::cout << y_skip[zi][yi] << " ";
            //     }
            // }
            // std::cout << "\n";
            // for (auto zi = 0; zi < 3; ++zi) {
            //     for (auto yi = 0; yi < 3; ++yi) {
            //         for (auto xi = 0; xi < 3; ++xi) {
            //             std::cout << x_skip[zi][yi][xi] << " ";
            //         }
            //     }
            // }
            // std::cout << "\n";
            // std::cin.get();
            // }
        }

        void monopole_interaction_interface::compute_interactions(interaction_kernel_type p2p_type,
            interaction_kernel_type p2m_type,
            std::array<bool, geo::direction::count()>& is_direction_empty,
            std::vector<neighbor_gravity_type>& all_neighbor_interaction_data) {
            // auto start_complete = std::chrono::high_resolution_clock::now();
            // Convert input structure to new datastructure (SoA)
            if (p2p_type == interaction_kernel_type::SOA_CPU &&
                p2m_type == interaction_kernel_type::SOA_CPU) {
                potential_expansions_SoA.update_data(potential_expansions);
                // auto start = std::chrono::high_resolution_clock::now();
                if (multipole_neighbors_exist) {
                    angular_corrections_SoA.update_data(angular_corrections);
                    kernel.apply_stencil(local_expansions_SoA_array, center_of_masses_SoA_array,
                        potential_expansions_SoA, angular_corrections_SoA, stencil, type, z_skip,
                        y_skip, x_skip);
                    if (type == RHO) {
                        angular_corrections_SoA.to_non_SoA(angular_corrections);
                        std::vector<space_vector>& L_c = grid_ptr->get_L_c();
                        for (size_t i = 0; i < L_c.size(); i++) {
                            L_c[i] = angular_corrections[i];
                        }
                    }
                }
                kernel_monopoles.apply_stencil(
                    local_monopoles, potential_expansions_SoA, stencil, four, dx);
                potential_expansions_SoA.to_non_SoA(potential_expansions);

                std::vector<expansion>& L = grid_ptr->get_L();
                for (size_t i = 0; i < L.size(); i++) {
                    L[i] = potential_expansions[i];
                }

            } else if (p2p_type == interaction_kernel_type::SOA_CPU) {
                potential_expansions_SoA.update_data(potential_expansions);
                kernel_monopoles.apply_stencil(
                    local_monopoles, potential_expansions_SoA, stencil, four, dx);
                potential_expansions_SoA.to_non_SoA(potential_expansions);

                std::vector<expansion>& L = grid_ptr->get_L();
                std::vector<space_vector>& L_c = grid_ptr->get_L_c();
                std::fill(std::begin(L), std::end(L), ZERO);
                std::fill(std::begin(L_c), std::end(L_c), ZERO);
                for (auto const& dir : geo::direction::full_set()) {
                    if (!is_direction_empty[dir]) {
                        neighbor_gravity_type& neighbor_data = all_neighbor_interaction_data[dir];
                        if (!neighbor_data.is_monopole) {
                            grid_ptr->compute_boundary_interactions(type, neighbor_data.direction,
                                neighbor_data.is_monopole, neighbor_data.data);
                        }
                    }
                }
                for (size_t i = 0; i < L.size(); i++) {
                    L[i] += potential_expansions[i];
                }
                for (size_t i = 0; i < L_c.size(); i++) {
                    L_c[i] += angular_corrections[i];
                }
            } else if (p2m_type == interaction_kernel_type::SOA_CPU) {
                if (multipole_neighbors_exist) {
                    potential_expansions_SoA.update_data(potential_expansions);
                    angular_corrections_SoA.update_data(angular_corrections);
                    kernel.apply_stencil(local_expansions_SoA_array, center_of_masses_SoA_array,
                        potential_expansions_SoA, angular_corrections_SoA, stencil, type, z_skip,
                        y_skip, x_skip);
                    potential_expansions_SoA.to_non_SoA(potential_expansions);
                    if (type == RHO) {
                        angular_corrections_SoA.to_non_SoA(angular_corrections);
                        std::vector<space_vector>& L_c = grid_ptr->get_L_c();
                        for (size_t i = 0; i < L_c.size(); i++) {
                            L_c[i] = angular_corrections[i];
                        }
                    }
                }

                std::vector<expansion>& L = grid_ptr->get_L();
                std::fill(std::begin(L), std::end(L), ZERO);
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
                for (size_t i = 0; i < L.size(); i++) {
                    L[i] += potential_expansions[i];
                }
            } else {
                // old-style interaction calculation
                // computes inner interactions
                grid_ptr->compute_interactions(type);
                // waits for boundary data and then computes boundary interactions
                for (auto const& dir : geo::direction::full_set()) {
                    if (!is_direction_empty[dir]) {
                        neighbor_gravity_type& neighbor_data = all_neighbor_interaction_data[dir];
                        grid_ptr->compute_boundary_interactions(type, neighbor_data.direction,
                            neighbor_data.is_monopole, neighbor_data.data);
                    }
                }
            }
        }

        std::vector<expansion>& monopole_interaction_interface::get_potential_expansions() {
            return potential_expansions;
        }

        std::vector<space_vector>& monopole_interaction_interface::get_angular_corrections() {
            return angular_corrections;
        }

        void monopole_interaction_interface::print_potential_expansions() {
            print_layered_not_padded(true, [this](const multiindex<>& i, const size_t flat_index) {
                std::cout << " (" << i << ") =[0] " << this->potential_expansions[flat_index][0];
            });
        }

        void monopole_interaction_interface::print_angular_corrections() {
            print_layered_not_padded(true, [this](const multiindex<>& i, const size_t flat_index) {
                std::cout << " (" << i << ") =[0] " << this->angular_corrections[flat_index];
            });
        }

        void monopole_interaction_interface::print_local_expansions() {
            // print_layered_padded(true, [this](const multiindex<>& i, const size_t flat_index) {
            //     std::cout << " " << this->local_expansions[flat_index];
            // });
        }

        void monopole_interaction_interface::print_center_of_masses() {
            // print_layered_padded(true, [this](const multiindex<>& i, const size_t flat_index) {
            //     std::cout << this->center_of_masses[flat_index];
            // });
        }

        void monopole_interaction_interface::add_to_potential_expansions(
            std::vector<expansion>& L) {
            // if (!multipole_neighbors_exist)
            //     return;
            iterate_inner_cells_not_padded([this, &L](multiindex<>& i, size_t flat_index) {
                potential_expansions[flat_index] += L[flat_index];
            });
        }

        void monopole_interaction_interface::add_to_center_of_masses(
             std::vector<space_vector>& L_c) {
            // // if (!multipole_neighbors_exist)
            // //     return;
            // iterate_inner_cells_not_padded([this, &L_c](multiindex<>& i, size_t flat_index) {
            //     center_of_masses[flat_index] += L_c[flat_index];
            // });
        }

    }    // namespace monopole_interactions
}    // namespace fmm
}    // namespace octotiger
