#pragma once

#include "interaction_constants.hpp"
#include "struct_of_array_data.hpp"
#include "taylor.hpp"

namespace octotiger {
namespace fmm {
namespace cuda {

class m2m_cuda {
public:
  void compute_interactions(
      octotiger::fmm::struct_of_array_data<expansion, real, 20, ENTRIES,
                                           SOA_PADDING>
          &local_expansions_SoA,
      octotiger::fmm::struct_of_array_data<space_vector, real, 3, ENTRIES,
                                           SOA_PADDING>
          &center_of_masses_SoA,
      octotiger::fmm::struct_of_array_data<expansion, real, 20, ENTRIES,
                                           SOA_PADDING>
          &potential_expansions_SoA,
      octotiger::fmm::struct_of_array_data<space_vector, real, 3, ENTRIES,
                                           SOA_PADDING>
          &angular_corrections_SoA);
};
}
}
}
