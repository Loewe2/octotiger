#include "multipole_interaction_interface.hpp"

#ifdef OCTOTIGER_WITH_CUDA
#include "m2m_cuda.hpp"
#endif

#include "../common_kernel/interactions_iterators.hpp"
#include "calculate_stencil.hpp"
#include "m2m_kernel.hpp"
#include "m2p_kernel.hpp"

#include <algorithm>

#include "options.hpp"

// required for cuda
extern options opts;
extern taylor<4, real> factor;
extern taylor<4, real> factor_half;
extern taylor<4, real> factor_sixth;

namespace octotiger {
namespace fmm {
    namespace multipole_interactions {

        two_phase_stencil multipole_interaction_interface::stencil;
        multipole_interaction_interface::multipole_interaction_interface(void)
          : neighbor_empty_multipole(27)
          , neighbor_empty_monopole(27)
          , mixed_interactions_kernel(neighbor_empty_monopole) {
            local_monopoles = std::vector<real>(EXPANSION_COUNT_PADDED);
        }

        void multipole_interaction_interface::update_input(std::vector<real>& monopoles,
            std::vector<multipole>& M_ptr,
            std::vector<std::shared_ptr<std::vector<space_vector>>>& com_ptr,
            std::vector<neighbor_gravity_type>& neighbors, gsolve_type t, real dx,
            std::array<real, NDIM> xbase) {
            type = t;
            dX = dx;
            xBase = xbase;
            std::vector<space_vector> const& com0 = *(com_ptr[0]);

            iterate_inner_cells_padded(
                [this, M_ptr, com0](const multiindex<>& i, const size_t flat_index,
                    const multiindex<>& i_unpadded, const size_t flat_index_unpadded) {
                    // local_expansions.at(flat_index) = M_ptr.at(flat_index_unpadded);
                    // center_of_masses.at(flat_index) = com0.at(flat_index_unpadded);
                    local_expansions_SoA.set_AoS_value(
                        std::move(M_ptr.at(flat_index_unpadded)), flat_index);
                    center_of_masses_SoA.set_AoS_value(
                        std::move(com0.at(flat_index_unpadded)), flat_index);
                    local_monopoles.at(flat_index) = 0.0;
                });

            for (size_t i = 0; i < neighbor_empty_multipole.size(); i++) {
                neighbor_empty_multipole[i] = false;
                neighbor_empty_monopole[i] = false;
            }

            monopole_neighbors_exist = false;
            for (const geo::direction& dir : geo::direction::full_set()) {
                // don't use neighbor.direction, is always zero for empty cells!
                neighbor_gravity_type& neighbor = neighbors[dir];
                // Switch x and z dimension since the stencil code and the old octotiger code
                // a different order for some reason
                auto x = dir.operator[](ZDIM) + 1;
                auto y = dir.operator[](YDIM) + 1;
                auto z = dir.operator[](XDIM) + 1;

                // this dir is setup as a multipole
                if (!neighbor.is_monopole) {
                    neighbor_empty_monopole[dir.flat_index_with_center()] = true;
                    x_skip[z][y][x] = true;
                    if (!neighbor.data.M) {
                        // TODO: ask Dominic why !is_monopole and stuff still empty
                        iterate_inner_cells_padding(dir, [this](const multiindex<>& i,
                                                             const size_t flat_index,
                                                             const multiindex<>&, const size_t) {
                            // // initializes whole expansion, relatively expansion
                            // local_expansions.at(flat_index) = 0.0;
                            // // initializes x,y,z vector
                            // center_of_masses.at(flat_index) = 0.0;
                            local_expansions_SoA.set_AoS_value(std::move(expansion()), flat_index);
                            center_of_masses_SoA.set_AoS_value(
                                std::move(space_vector()), flat_index);
                            local_monopoles.at(flat_index) = 0.0;
                        });
                        neighbor_empty_multipole[dir.flat_index_with_center()] = true;
                    } else {
                        std::vector<multipole>& neighbor_M_ptr = *(neighbor.data.M);
                        std::vector<space_vector>& neighbor_com0 = *(neighbor.data.x);
                        iterate_inner_cells_padding(
                            dir, [this, neighbor_M_ptr, neighbor_com0](const multiindex<>& i,
                                     const size_t flat_index, const multiindex<>& i_unpadded,
                                     const size_t flat_index_unpadded) {
                                // local_expansions.at(flat_index) =
                                //     neighbor_M_ptr.at(flat_index_unpadded);
                                // center_of_masses.at(flat_index) =
                                //     neighbor_com0.at(flat_index_unpadded);
                                local_expansions_SoA.set_AoS_value(
                                    std::move(neighbor_M_ptr.at(flat_index_unpadded)), flat_index);
                                center_of_masses_SoA.set_AoS_value(
                                    std::move(neighbor_com0.at(flat_index_unpadded)), flat_index);
                                local_monopoles.at(flat_index) = 0.0;

                            });
                    }
                } else {
                    neighbor_empty_multipole[dir.flat_index_with_center()] = true;
                    // neighbor has no data - input structure just recevices zeros as padding
                    if (!neighbor.data.m) {
                        iterate_inner_cells_padding(dir, [this](const multiindex<>& i,
                                                             const size_t flat_index,
                                                             const multiindex<>&, const size_t) {
                            // initializes whole expansion, relatively expansion
                            // local_expansions.at(flat_index) = 0.0;
                            // // initializes x,y,z vector
                            // center_of_masses.at(flat_index) = 0.0;
                            local_expansions_SoA.set_AoS_value(std::move(expansion()), flat_index);
                            center_of_masses_SoA.set_AoS_value(
                                std::move(space_vector()), flat_index);
                            local_monopoles.at(flat_index) = 0.0;

                        });
                        neighbor_empty_monopole[dir.flat_index_with_center()] = true;
                        x_skip[z][y][x] = true;
                    } else {
                        // Get multipole data into our input structure
                        std::vector<real>& neighbor_mons = *(neighbor.data.m);
                        std::vector<space_vector>& neighbor_com0 = *(neighbor.data.x);
                        iterate_inner_cells_padding(dir, [this, neighbor_mons, xbase, dx](
                                                             const multiindex<>& i,
                                                             const size_t flat_index,
                                                             const multiindex<>& i_unpadded,
                                                             const size_t flat_index_unpadded) {
                            // local_expansions.at(flat_index) = 0.0;
                            // center_of_masses.at(flat_index) = 0.0;

                            space_vector e;
                            e[0] = (i.x) * dx + xbase[0] - INNER_CELLS_PER_DIRECTION * dx;
                            e[1] = (i.y) * dx + xbase[1] - INNER_CELLS_PER_DIRECTION * dx;
                            e[2] = (i.z) * dx + xbase[2] - INNER_CELLS_PER_DIRECTION * dx;
                            center_of_masses_SoA.set_AoS_value(std::move(e), flat_index);
                            // local_monopoles.at(flat_index) =
                            // neighbor_mons.at(flat_index_unpadded);
                            local_expansions_SoA.set_AoS_value(std::move(expansion()), flat_index);
                            // local_expansions_SoA.set_value(
                            //     std::move(neighbor_mons.at(flat_index_unpadded)), flat_index);
                            local_monopoles.at(flat_index) = neighbor_mons.at(flat_index_unpadded);
                        });
                        monopole_neighbors_exist = true;
                        x_skip[z][y][x] = false;
                    }
                }
            }

            neighbor_empty_multipole[13] = false;
            neighbor_empty_monopole[13] = true;

            x_skip[1][1][1] = true;
            for (auto zi = 0; zi < 3; ++zi) {
                z_skip[zi] = true;
                for (auto yi = 0; yi < 3; ++yi) {
                    y_skip[zi][yi] = true;
                    for (auto xi = 0; xi < 3; ++xi) {
                        if (!x_skip[zi][yi][xi]) {
                            y_skip[zi][yi] = false;
                            break;
                        }
                    }
                    if (!y_skip[zi][yi])
                        z_skip[zi] = false;
                }
            }

            // std::fill(std::begin(potential_expansions), std::end(potential_expansions), ZERO);

            // std::cout << "local_expansions:" << std::endl;
            // this->print_local_expansions();
            // std::cout << "center_of_masses:" << std::endl;
            // this->print_center_of_masses();
        }

        void multipole_interaction_interface::compute_interactions(interaction_kernel_type m2m_type,
            interaction_kernel_type m2p_type,
            std::array<bool, geo::direction::count()>& is_direction_empty,
            std::vector<neighbor_gravity_type>& all_neighbor_interaction_data) {
            if (m2m_type == interaction_kernel_type::SOA_CPU &&
                m2p_type == interaction_kernel_type::SOA_CPU) {
#ifdef OCTOTIGER_WITH_CUDA
                struct_of_array_data<real, 20, INNER_CELLS, SOA_PADDING>
                    potential_expansions_SoA;
                for (size_t i = 0; i < potential_expansions_SoA.size(); i++) {
                    potential_expansions_SoA.get_underlying_pointer()[i] = 0.0;
                }
                struct_of_array_data<real, 3, INNER_CELLS, SOA_PADDING>
                    angular_corrections_SoA;
                for (size_t i = 0; i < angular_corrections_SoA.size(); i++) {
                    angular_corrections_SoA.get_underlying_pointer()[i] = 0.0;
                }
                cuda::m2m_cuda cuda_kernel;
                // auto start = std::chrono::high_resolution_clock::now();
                cuda_kernel.compute_interactions(local_expansions_SoA, center_of_masses_SoA,
                    potential_expansions_SoA, angular_corrections_SoA, opts.theta,
                    factor.get_array(), factor_half.get_array(), factor_sixth.get_array(), type);
                // auto end = std::chrono::high_resolution_clock::now();

                if (monopole_neighbors_exist) {
                    mixed_interactions_kernel.apply_stencil(local_monopoles,
                    local_expansions_SoA,
                        center_of_masses_SoA, potential_expansions_SoA, angular_corrections_SoA,
                        stencil.stencil_elements, type, dX, xBase, x_skip, y_skip, z_skip);
                }
#else
                struct_of_array_data<real, 20, INNER_CELLS, SOA_PADDING>
                    potential_expansions_SoA;
                struct_of_array_data<real, 3, INNER_CELLS, SOA_PADDING>
                    angular_corrections_SoA;

                // if (monopole_neighbors_exist) {
                //     mixed_interactions_kernel.apply_stencil(local_monopoles,
                //     local_expansions_SoA,
                //         center_of_masses_SoA, potential_expansions_SoA, angular_corrections_SoA,
                //         stencil_mixed_interactions, type, dX, xBase, x_skip, y_skip, z_skip);
                // }
                m2m_kernel kernel(neighbor_empty_multipole);
                kernel.apply_stencil(local_expansions_SoA, center_of_masses_SoA,
                    potential_expansions_SoA, angular_corrections_SoA, local_monopoles, stencil,
                    type);

#endif
                if (type == RHO) {
                  angular_corrections_SoA.to_non_SoA<space_vector>(grid_ptr->get_L_c());
                }

                potential_expansions_SoA.add_to_non_SoA<expansion>(grid_ptr->get_L());
            } else if (m2m_type == interaction_kernel_type::SOA_CPU) {
                struct_of_array_data<real, 20, INNER_CELLS, SOA_PADDING>
                    potential_expansions_SoA;
                struct_of_array_data<real, 3, INNER_CELLS, SOA_PADDING>
                    angular_corrections_SoA;
#ifdef OCTOTIGER_WITH_CUDA
#else
                m2m_kernel kernel(neighbor_empty_multipole);
                kernel.apply_stencil(local_expansions_SoA, center_of_masses_SoA,
                    potential_expansions_SoA, angular_corrections_SoA, local_monopoles, stencil,
                    type);
#endif

                std::vector<expansion>& L = grid_ptr->get_L();
                std::vector<space_vector>& L_c = grid_ptr->get_L_c();
                std::fill(std::begin(L), std::end(L), ZERO);
                std::fill(std::begin(L_c), std::end(L_c), ZERO);
                for (auto const& dir : geo::direction::full_set()) {
                    if (!is_direction_empty[dir]) {
                        neighbor_gravity_type& neighbor_data = all_neighbor_interaction_data[dir];
                        if (neighbor_data.is_monopole) {
                            grid_ptr->compute_boundary_interactions(type, neighbor_data.direction,
                                neighbor_data.is_monopole, neighbor_data.data);
                        }
                    }
                }
                potential_expansions_SoA.add_to_non_SoA<expansion>(grid_ptr->get_L());
                angular_corrections_SoA.add_to_non_SoA<space_vector>(grid_ptr->get_L_c());

            } else if (m2p_type == interaction_kernel_type::SOA_CPU) {
                struct_of_array_data<real, 20, INNER_CELLS, SOA_PADDING>
                    potential_expansions_SoA;
                struct_of_array_data<real, 3, INNER_CELLS, SOA_PADDING>
                    angular_corrections_SoA;
                // if (monopole_neighbors_exist) {
                //     mixed_interactions_kernel.apply_stencil(local_monopoles,
                //     local_expansions_SoA,
                //         center_of_masses_SoA, potential_expansions_SoA, angular_corrections_SoA,
                //         stencil_mixed_interactions, type, dX, xBase, x_skip, y_skip, z_skip);
                // }
                m2m_kernel kernel(neighbor_empty_multipole);
                kernel.apply_stencil(local_expansions_SoA, center_of_masses_SoA,
                    potential_expansions_SoA, angular_corrections_SoA, local_monopoles, stencil,
                    type);

                std::vector<expansion>& L = grid_ptr->get_L();
                std::vector<space_vector>& L_c = grid_ptr->get_L_c();
                std::fill(std::begin(L), std::end(L), ZERO);
                std::fill(std::begin(L_c), std::end(L_c), ZERO);

                grid_ptr->compute_interactions(type);
                // waits for boundary data and then computes boundary interactions
                for (auto const& dir : geo::direction::full_set()) {
                    if (!is_direction_empty[dir]) {
                        neighbor_gravity_type& neighbor_data = all_neighbor_interaction_data[dir];
                        if (!neighbor_data.is_monopole) {
                            grid_ptr->compute_boundary_interactions(type, neighbor_data.direction,
                                neighbor_data.is_monopole, neighbor_data.data);
                        }
                    }
                }
                potential_expansions_SoA.add_to_non_SoA<expansion>(grid_ptr->get_L());
                angular_corrections_SoA.add_to_non_SoA<space_vector>(grid_ptr->get_L_c());
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
    }    // namespace multipole_interactions
}    // namespace fmm
}    // namespace octotiger
