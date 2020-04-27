
#include "octotiger/common_kernel/interactions_iterators.hpp"
#include "octotiger/monopole_interactions/calculate_stencil.hpp"
#include "octotiger/monopole_interactions/p2m_interaction_interface.hpp"
#include "octotiger/options.hpp"
#include "octotiger/real.hpp"
#include "octotiger/monopole_interactions/p2m_kernel.hpp"
// #include "octotiger/taylor.hpp"

// #include <array>
// #include <memory>

#include <boost/serialization/vector.hpp>
#include <boost/serialization/array.hpp>

#include <vector>

#include <gtest/gtest.h>

 std::vector< octotiger::fmm::multiindex<>>& stencil()
        {
            
            static thread_local std::vector< octotiger::fmm::multiindex<>> stencil_ =
                 octotiger::fmm::monopole_interactions::calculate_stencil().first;
            return stencil_;
        }

TEST(ExampleTests, DemonstrateGTestMacros){
    std::vector<bool> neighbor_empty_multipoles(27);
    octotiger::fmm::monopole_interactions::p2m_kernel kernel(neighbor_empty_multipoles);
    std::vector<space_vector_gen<double> > center_of_masses_staging_area_vector;
    std::vector<taylor<4, double> > potential_expansions_SoA_vector;
    std::vector<space_vector_gen<double> > angular_corrections_SoA_vector;

     // create and open an archive for input
    std::ifstream ifs("/home/alexander/Documents/Uni/Masterarbeit/octotiger/build/filename");
    boost::archive::text_iarchive ia(ifs);
    int size;
    ia >> size;
    std::vector<taylor<4, double> > local_expansions_staging_area_vector(size);
    std::cout << "site: " << size<< std::endl;;
    for(int i  = 0; i < size; ++i){
        std::cout << i << std::endl;
        size_t temp;
        ia >> temp;
        // local_expansions_staging_area_vector[i].resize(temp);
        ia >> local_expansions_staging_area_vector[i];

    }
    // ia >> center_of_masses_staging_area_vector;
    // ia >> potential_expansions_SoA_vector;
    // ia >> angular_corrections_SoA_vector;


    gsolve_type type;
    // ia >> type;
    bool x_skip[3][3][3];
    bool y_skip[3][3];
    bool z_skip[3];
    // ia >> x_skip;
    // ia >> y_skip;
    // ia >> z_skip;
    std::cout << z_skip[0] << " " << z_skip[1] << " " << z_skip[2] << std::endl;
   
    octotiger::fmm::struct_of_array_data<expansion, real, 20, octotiger::fmm::ENTRIES, octotiger::fmm::SOA_PADDING> local_expansions_staging_area(local_expansions_staging_area_vector);
    octotiger::fmm::struct_of_array_data<space_vector, real, 3,  octotiger::fmm::ENTRIES,  octotiger::fmm::SOA_PADDING> center_of_masses_staging_area(center_of_masses_staging_area_vector);
    octotiger::fmm::struct_of_array_data<expansion, real, 20, octotiger::fmm::INNER_CELLS, octotiger::fmm::SOA_PADDING> potential_expansions_SoA(potential_expansions_SoA_vector);
    octotiger::fmm::struct_of_array_data<space_vector, real, 3, octotiger::fmm::INNER_CELLS, octotiger::fmm::SOA_PADDING> angular_corrections_SoA(angular_corrections_SoA_vector);

    kernel.apply_stencil(local_expansions_staging_area,
                        center_of_masses_staging_area, 
                        potential_expansions_SoA,
                        angular_corrections_SoA, 
                        stencil(), 
                        type,
                        x_skip, y_skip, z_skip);
    EXPECT_TRUE(true);
}

