#pragma once

#include "octotiger/common_kernel/multiindex.hpp"

namespace octotiger {
namespace fmm {
    namespace multipole_interactions {

        two_phase_stencil calculate_stencil(void);
        std::pair<oct::vector<bool>, oct::vector<bool>> calculate_stencil_masks(two_phase_stencil superimposed_stencil);

    }    // namespace multipole_interactions
}    // namespace fmm
}    // namespace octotiger
