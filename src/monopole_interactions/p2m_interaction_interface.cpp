//  Copyright (c) 2019 AUTHORS
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#include "octotiger/common_kernel/interactions_iterators.hpp"
#include "octotiger/monopole_interactions/calculate_stencil.hpp"
#include "octotiger/monopole_interactions/p2m_interaction_interface.hpp"
#include "octotiger/options.hpp"
#include "octotiger/real.hpp"

#include <algorithm>
#include <array>
#include <vector>
// #include <string>

#include <boost/serialization/vector.hpp>
#include <boost/serialization/array.hpp>
#include <boost/serialization/extended_type_info.hpp>

// Big picture questions:
// - use any kind of tiling?
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
namespace boost {
namespace serialization {

template<class Archive, typename T>
void serialize(Archive & ar, octotiger::fmm::multiindex<T> & g, const unsigned int version)
{
    ar & g.x;
    ar & g.y;
    ar & g.z;
}

} // namespace seriali
}
template <typename>
class my_trait;

template <typename AoS_type_t, typename component_type_t, size_t num_components_v, size_t entries_v, 
        size_t padding_v, typename backend_t_t>
struct my_trait<octotiger::fmm::struct_of_array_data<AoS_type_t, component_type_t, num_components_v, entries_v, padding_v, backend_t_t>> {
    using aos_type = AoS_type_t;
    using component_type = component_type_t;
    static constexpr size_t num_components = num_components_v;
    static constexpr size_t entries = entries_v;
    static constexpr size_t padding = padding_v;
    using backend_type = backend_t_t;
};


namespace octotiger {
namespace fmm {
    namespace monopole_interactions {
        std::vector<multiindex<>>& p2m_interaction_interface::stencil()
        {
            static thread_local std::vector<multiindex<>> stencil_ =
                calculate_stencil().first;
            return stencil_;
        }
        thread_local std::vector<real> p2m_interaction_interface::local_monopoles_staging_area(
            ENTRIES);
        thread_local struct_of_array_data<expansion, real, 20, ENTRIES, SOA_PADDING>
            p2m_interaction_interface::local_expansions_staging_area;
        thread_local struct_of_array_data<space_vector, real, 3, ENTRIES, SOA_PADDING>
            p2m_interaction_interface::center_of_masses_staging_area;

        p2m_interaction_interface::p2m_interaction_interface()
          : neighbor_empty_multipoles(27)
          , kernel(neighbor_empty_multipoles) {
            this->p2m_type = opts().p2m_kernel_type;
        }

        void p2m_interaction_interface::compute_p2m_interactions(std::vector<real>& monopoles,
            std::vector<multipole>& M_ptr,
            std::vector<std::shared_ptr<std::vector<space_vector>>>& com_ptr,
            std::vector<neighbor_gravity_type>& neighbors, gsolve_type type,
            std::array<bool, geo::direction::count()>& is_direction_empty) {
            update_input(monopoles, M_ptr, com_ptr, neighbors, type, local_monopoles_staging_area,
                local_expansions_staging_area, center_of_masses_staging_area);
            compute_interactions(type, is_direction_empty, neighbors);
        }

        void p2m_interaction_interface::compute_interactions(gsolve_type type,
            std::array<bool, geo::direction::count()>& is_direction_empty,
            std::vector<neighbor_gravity_type>& all_neighbor_interaction_data) {
            if (p2m_type == interaction_kernel_type::SOA_CPU) {
                if (multipole_neighbors_exist) {
                    struct_of_array_data<expansion, real, 20, INNER_CELLS, SOA_PADDING>
                        potential_expansions_SoA;
                    struct_of_array_data<space_vector, real, 3, INNER_CELLS, SOA_PADDING>
                        angular_corrections_SoA;
                        std::ofstream ofs("filename" + std::to_string(hpx::get_worker_thread_num()) );

                            boost::archive::text_oarchive oa(ofs);

                            std::vector<typename my_trait<decltype(local_expansions_staging_area)>::aos_type> local_expansions_staging_area_vector( my_trait<decltype(local_expansions_staging_area)>::entries);
                            std::vector<space_vector_gen<double> > center_of_masses_staging_area_vector(my_trait<decltype( center_of_masses_staging_area)>::entries);
                            std::vector<taylor<4, double> > potential_expansions_SoA_vector(my_trait<decltype(potential_expansions_SoA)>::entries);
                            std::vector<space_vector_gen<double> > angular_corrections_SoA_vector(my_trait<decltype(angular_corrections_SoA)>::entries);
                            local_expansions_staging_area.to_non_SoA(local_expansions_staging_area_vector);
                            center_of_masses_staging_area.to_non_SoA(center_of_masses_staging_area_vector);
                            potential_expansions_SoA.to_non_SoA(potential_expansions_SoA_vector);
                            angular_corrections_SoA.to_non_SoA(angular_corrections_SoA_vector);
                            int size = local_expansions_staging_area_vector.size();
        
                            oa << neighbor_empty_multipoles;
                            oa << local_expansions_staging_area_vector;
                            oa << center_of_masses_staging_area_vector;
    
                            oa << potential_expansions_SoA_vector;
                            oa << angular_corrections_SoA_vector;
     
                            oa << type;
                            oa << x_skip;
                            oa << y_skip;
                            oa << z_skip;
                            std::vector<multiindex<>>& stencil_ = stencil();
                            oa << stencil_;



                    kernel.apply_stencil(local_expansions_staging_area,
                        center_of_masses_staging_area, potential_expansions_SoA,
                        angular_corrections_SoA, stencil_, type, x_skip, y_skip, z_skip);
                    potential_expansions_SoA.to_non_SoA(potential_expansions_SoA_vector);
                    angular_corrections_SoA.to_non_SoA(angular_corrections_SoA_vector);
                    oa << potential_expansions_SoA_vector;
                    oa << angular_corrections_SoA_vector;
                    potential_expansions_SoA.add_to_non_SoA(grid_ptr->get_L());
                    if (type == RHO) {
                        angular_corrections_SoA.to_non_SoA(grid_ptr->get_L_c());
                    }
                }
            } else {
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
            }
        }
    }    // namespace monopole_interactions
}    // namespace fmm
}    // namespace octotiger
