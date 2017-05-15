#pragma once

#include <memory>
#include <vector>

#include "geometry.hpp"
#include "grid.hpp"
#include "node_server.hpp"
#include "taylor.hpp"

#include "interactions_constants.hpp"
#include "multiindex.hpp"

namespace octotiger {
namespace fmm {

    // for both local and multipole expansion
    // typedef taylor<4, real> expansion;

    class m2m_interactions
    {
    private:
        bool verbose;

        /*
         * logical structure of all arrays:
         * cube of 8x8x8 cells of configurable size,
         * input variables (local_expansions, center_of_masses) have
         * an additional 1-sized layer of cells around it for padding
         */

        // M_ptr
        std::vector<expansion> local_expansions;

        // com0 = *(com_ptr[0])
        std::vector<space_vector> center_of_masses;

        // multipole expansion on this cell (L)
        std::vector<expansion> potential_expansions;
        // angular momentum correction on this cell (L_c)
        std::vector<space_vector> angular_corrections;

        // inline size_t to_outer_flat_index(multiindex& m) {
        //     // skip volumes in the two "slow" dimensions, skip only single row in direction of
        //     // fast dimension
        //     return INNER_CELLS * (m.x * OUTER_CELLS_PER_DIRECTION * OUTER_CELLS_PER_DIRECTION +
        //                              m.y * OUTER_CELLS_PER_DIRECTION) +
        //         m.z * INNER_CELLS_PER_DIRECTION;
        // }

        // stride for multiple outer cells (and/or padding)
        inline size_t to_inner_flat_index_padded(const multiindex& m) {
            return m.x * PADDED_STRIDE * PADDED_STRIDE + m.y * PADDED_STRIDE + m.z;
        }

        // inline size_t to_flat_index_padded(multiindex& outer, multiindex& inner) {
        //     return to_outer_flat_index(outer) + to_inner_flat_index_padded(inner);
        // }

        // strides are only valid for single cell! (no padding)
        inline size_t to_inner_flat_index_not_padded(const multiindex& m) {
            return m.x * INNER_CELLS_PER_DIRECTION * INNER_CELLS_PER_DIRECTION +
                m.y * INNER_CELLS_PER_DIRECTION + m.z;
        }

        // meant to iterate the input data structure
        template <typename F>
        void iterate_inner_cells_padded(F f) {
            for (size_t i0 = 0; i0 < INNER_CELLS_PER_DIRECTION; i0++) {
                for (size_t i1 = 0; i1 < INNER_CELLS_PER_DIRECTION; i1++) {
                    for (size_t i2 = 0; i2 < INNER_CELLS_PER_DIRECTION; i2++) {
                        const multiindex m(i0 + INNER_CELLS_PADDING_DEPTH,
                            i1 + INNER_CELLS_PADDING_DEPTH, i2 + INNER_CELLS_PADDING_DEPTH);
                        const size_t inner_flat_index = to_inner_flat_index_padded(m);
                        const multiindex m_unpadded(i0, i1, i2);
                        const size_t inner_flat_index_unpadded =
                            to_inner_flat_index_not_padded(m_unpadded);
                        f(m, inner_flat_index, m_unpadded, inner_flat_index_unpadded);
                    }
                }
            }
        }

        // meant to iterate the input data structure
        template <typename F>
        void iterate_inner_cells_padding(const geo::direction& dir, F f) {
            // TODO: implementation not finished
            for (size_t i0 = 0; i0 < INNER_CELLS_PER_DIRECTION; i0++) {
                for (size_t i1 = 0; i1 < INNER_CELLS_PER_DIRECTION; i1++) {
                    for (size_t i2 = 0; i2 < INNER_CELLS_PER_DIRECTION; i2++) {
                        const multiindex m(
                            i0 + INNER_CELLS_PADDING_DEPTH + dir[0] * INNER_CELLS_PADDING_DEPTH,
                            i1 + INNER_CELLS_PADDING_DEPTH + dir[1] * INNER_CELLS_PADDING_DEPTH,
                            i2 + INNER_CELLS_PADDING_DEPTH + dir[2] * INNER_CELLS_PADDING_DEPTH);
                        const size_t inner_flat_index = to_inner_flat_index_padded(m);
                        const multiindex m_unpadded(i0, i1, i2);
                        const size_t inner_flat_index_unpadded =
                            to_inner_flat_index_not_padded(m_unpadded);
                        f(m, inner_flat_index, m_unpadded, inner_flat_index_unpadded);
                    }
                }
            }
        }

        // meant to iterate the output data structure
        template <typename F>
        void iterate_inner_cells_not_padded(F f) {
            for (size_t i0 = 0; i0 < INNER_CELLS_PER_DIRECTION; i0++) {
                for (size_t i1 = 0; i1 < INNER_CELLS_PER_DIRECTION; i1++) {
                    for (size_t i2 = 0; i2 < INNER_CELLS_PER_DIRECTION; i2++) {
                        multiindex i(i0, i1, i2);
                        size_t inner_flat_index = to_inner_flat_index_not_padded(i);
                        f(i, inner_flat_index);
                    }
                }
            }
        }

    public:
        // at this point, uses the old datamembers of the grid class as input
        // and converts them to the new data structure
        m2m_interactions(grid& g, std::vector<node_server::neighbor_gravity_type>& neighbors);

        // dummy constructor for debugging
        m2m_interactions();

        ~m2m_interactions() {
            std::cout << "in destructor" << std::endl;
        }

        void get_converted_local_expansions(std::vector<multipole>& M_ptr);

        void get_converted_center_of_masses(
            std::vector<std::shared_ptr<std::vector<space_vector>>> com_ptr);

        void get_converted_potential_expansions(std::vector<expansion>& L);

        void get_converted_angular_corrections(std::vector<space_vector>& L_c);
    };
}    // namespace fmm
}    // namespace octotiger
