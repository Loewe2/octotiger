#pragma once

#include <memory>
#include <vector>

#include "geometry.hpp"
#include "grid.hpp"
#include "node_server.hpp"
#include "taylor.hpp"

#include "interaction_constants.hpp"
#include "multiindex.hpp"

namespace octotiger {
namespace fmm {

    // for both local and multipole expansion
    // typedef taylor<4, real> expansion;

    class m2m_interactions
    {
    private:
        bool verbose;
        std::vector<multiindex> stencil;

        /*
         * logical structure of all arrays:
         * cube of 8x8x8 cells of configurable size,
         * input variables (local_expansions, center_of_masses) have
         * an additional 1-sized layer of cells around it for padding
         */

        // M_ptr
        std::vector<expansion> local_expansions;    // TODO: needs to be converted to SoA

        // com0 = *(com_ptr[0])
        std::vector<space_vector> center_of_masses;    // TODO: needs to be converted to SoA

        // multipole expansion on this cell (L)
        std::vector<expansion> potential_expansions;    // TODO: needs to be converted to SoA
        // angular momentum correction on this cell (L_c)
        std::vector<space_vector> angular_corrections;    // TODO: needs to be converted to SoA

        gsolve_type type;

    public:
        // at this point, uses the old datamembers of the grid class as input
        // and converts them to the new data structure
        m2m_interactions(
            grid& g, std::vector<node_server::neighbor_gravity_type>& neighbors, gsolve_type type);

        // dummy constructor for debugging
        m2m_interactions();

        void compute_interactions();

        void get_converted_local_expansions(std::vector<multipole>& M_ptr);

        void get_converted_center_of_masses(
            std::vector<std::shared_ptr<std::vector<space_vector>>> com_ptr);

        void get_converted_potential_expansions(std::vector<expansion>& L);

        void get_converted_angular_corrections(std::vector<space_vector>& L_c);

        void print_potential_expansions();

        void print_angular_corrections();
    };
}    // namespace fmm
}    // namespace octotiger
