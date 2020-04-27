//  Copyright (c) 2019 AUTHORS
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#pragma once

#include "octotiger/common_kernel/interaction_constants.hpp"
#include "octotiger/common_kernel/interactions_iterators.hpp"
#include "octotiger/common_kernel/kernel_simd_types.hpp"
#include "octotiger/common_kernel/multiindex.hpp"
#include "octotiger/common_kernel/struct_of_array_data.hpp"
#include "octotiger/taylor.hpp"

#include <cstdint>
#include <vector>

#include <fstream>

// include headers that implement a archive in simple text format
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <boost/serialization/vector.hpp>

namespace octotiger {
namespace fmm {
    namespace monopole_interactions {
        constexpr uint64_t P2M_STENCIL_BLOCKING = 1;

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
            void blocked_interaction_rho(const struct_of_array_data<expansion, real, 20, ENTRIES,
                                             SOA_PADDING>& local_expansions_SoA,
                const struct_of_array_data<space_vector, real, 3, ENTRIES, SOA_PADDING>&
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
            void blocked_interaction_non_rho(const struct_of_array_data<expansion, real, 20, ENTRIES,
                                                 SOA_PADDING>& local_expansions_SoA,
                const struct_of_array_data<space_vector, real, 3, ENTRIES, SOA_PADDING>&
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

            friend class boost::serialization::access;
            template<class Archive>
            void serialize(Archive & ar, const unsigned int version)
            {
                ar & neighbor_empty;
                ar & vector_is_empty;
                // ar & theta_rec_squared;
                // ar & offset_vector;
            }
        public:
            p2m_kernel(std::vector<bool>& neighbor_empty);

            p2m_kernel(p2m_kernel& other) = delete;
            p2m_kernel(const p2m_kernel& other) = delete;
            p2m_kernel operator=(const p2m_kernel& other) = delete;

            void apply_stencil(const struct_of_array_data<expansion, real, 20, ENTRIES, SOA_PADDING>&
                                   local_expansions_SoA,
                const struct_of_array_data<space_vector, real, 3, ENTRIES, SOA_PADDING>&
                    center_of_masses_SoA,
                struct_of_array_data<expansion, real, 20, INNER_CELLS, SOA_PADDING>&
                    potential_expansions_SoA,
                struct_of_array_data<space_vector, real, 3, INNER_CELLS, SOA_PADDING>&
                    angular_corrections_SoA,
                const std::vector<multiindex<>>& stencil, gsolve_type type, const bool (&z_skip)[3][3][3],
                const bool (&y_skip)[3][3], const bool (&x_skip)[3]);
        };

    }    // namespace monopole_interactions
}    // namespace fmm
}    // namespace octotiger


// namespace boost {
// namespace serialization {

// template<class Archive>
// void serialize(Archive & ar, octotiger::fmm::monopole_interactions::p2m_kernel & kernel, const unsigned int version)
// {
//     ar & kernel.neighbor_empty;
//     ar & kernel.vector_is_empty;
//     ar & kernel.theta_rec_squared;
//     ar & kernel.offset_vector;
// }

// } // namespace serialization
// } // namespace boost