
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


TEST(ExampleTests, DemonstrateGTestMacros){
    
    std::vector<space_vector_gen<double> > center_of_masses_staging_area_vector;
    std::vector<space_vector_gen<double> > angular_corrections_SoA_vector;
    std::vector<space_vector_gen<double> > angular_corrections_SoA_result_vector_correct;

     // create and open an archive for input
    std::ifstream ifs("/home/alexander/Documents/Uni/Masterarbeit/octotiger/build/filename1");
    boost::archive::text_iarchive ia(ifs);
    // std::cout << "site: " << size<< std::endl;;
    // for(int i  = 0; i < size; ++i){
    //     std::cout << i << std::endl;
    //     size_t temp;
    //     ia >> temp;
    //     // local_expansions_staging_area_vector[i].resize(temp);

    // }
    std::vector<bool> neighbor_empty_multipoles;
    ia >> neighbor_empty_multipoles;
    octotiger::fmm::monopole_interactions::p2m_kernel kernel(neighbor_empty_multipoles);

    std::vector<taylor<4, double> > local_expansions_staging_area_vector;
    std::vector<taylor<4, double> > potential_expansions_SoA_vector;
    std::vector<taylor<4, double> > potential_expansions_SoA_result_vector_correct;
    ia >> local_expansions_staging_area_vector;
    ia >> center_of_masses_staging_area_vector;
    ia >> potential_expansions_SoA_vector;
    ia >> angular_corrections_SoA_vector;


    gsolve_type type;
    ia >> type;
    bool x_skip[3][3][3];
    bool y_skip[3][3];
    bool z_skip[3];
    ia >> x_skip;
    ia >> y_skip;
    ia >> z_skip;
    std::cout << z_skip[0] << " " << z_skip[1] << " " << z_skip[2] << std::endl;
    std::vector<octotiger::fmm::multiindex<>> stencil_;
    ia >> stencil_;
   
    octotiger::fmm::struct_of_array_data<expansion, real, 20, octotiger::fmm::ENTRIES, octotiger::fmm::SOA_PADDING> local_expansions_staging_area(local_expansions_staging_area_vector);
    octotiger::fmm::struct_of_array_data<space_vector, real, 3,  octotiger::fmm::ENTRIES,  octotiger::fmm::SOA_PADDING> center_of_masses_staging_area(center_of_masses_staging_area_vector);
    octotiger::fmm::struct_of_array_data<expansion, real, 20, octotiger::fmm::INNER_CELLS, octotiger::fmm::SOA_PADDING> potential_expansions_SoA(potential_expansions_SoA_vector);
    octotiger::fmm::struct_of_array_data<space_vector, real, 3, octotiger::fmm::INNER_CELLS, octotiger::fmm::SOA_PADDING> angular_corrections_SoA(angular_corrections_SoA_vector);

    kernel.apply_stencil(local_expansions_staging_area,
                        center_of_masses_staging_area, 
                        potential_expansions_SoA,
                        angular_corrections_SoA, 
                        stencil_, 
                        type,
                        x_skip, y_skip, z_skip);

    std::vector<taylor<4, double> > potential_expansions_SoA_result_vector(my_trait<decltype(potential_expansions_SoA)>::entries);
    std::vector<space_vector_gen<double> > angular_corrections_SoA_result_vector(my_trait<decltype(angular_corrections_SoA)>::entries);
    potential_expansions_SoA.to_non_SoA(potential_expansions_SoA_result_vector);
    angular_corrections_SoA.to_non_SoA(angular_corrections_SoA_result_vector);
    ia >> potential_expansions_SoA_result_vector_correct;
    ia >> angular_corrections_SoA_result_vector_correct;

    ASSERT_EQ(angular_corrections_SoA_result_vector.size(), angular_corrections_SoA_result_vector_correct.size());
    for (std::size_t i = 0; i < angular_corrections_SoA_result_vector.size(); ++i) {
        SCOPED_TRACE(i);
        for(std::size_t j = 0; j < NDIM; ++j){
            SCOPED_TRACE(j);
            EXPECT_NEAR (angular_corrections_SoA_result_vector[i][j], angular_corrections_SoA_result_vector_correct[i][j], 1e-19);
        }
    }

    // // ASSERT_EQ(potential_expansions_SoA_result_vector.size(), potential_expansions_SoA_result_vector_correct.size());
    for (std::size_t i = 0; i < potential_expansions_SoA_result_vector.size(); ++i) {
        SCOPED_TRACE(i);
        for(std::size_t j = 0; j < 4; ++j){
            SCOPED_TRACE(j);
            EXPECT_NEAR(potential_expansions_SoA_result_vector[i][j], potential_expansions_SoA_result_vector_correct[i][j], 1e-14);
        }

    }

    EXPECT_TRUE(true);
}

