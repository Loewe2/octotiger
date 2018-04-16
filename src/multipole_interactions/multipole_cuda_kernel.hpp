#pragma once
#ifdef OCTOTIGER_CUDA_ENABLED
#include "../common_kernel/interaction_constants.hpp"
#include "../common_kernel/multiindex.hpp"

namespace octotiger {
namespace fmm {
    namespace multipole_interactions {
        __global__ void cuda_multipole_interactions_kernel_rho(
            const double (&local_monopoles)[NUMBER_LOCAL_MONOPOLE_VALUES],
            const double (&center_of_masses)[NUMBER_MASS_VALUES],
            const double (&multipoles)[NUMBER_LOCAL_EXPANSION_VALUES],
            double (&potential_expansions)[3 * NUMBER_POT_EXPANSIONS],
            double (&angular_corrections)[3 * NUMBER_ANG_CORRECTIONS],
            const octotiger::fmm::multiindex<> (&stencil)[STENCIL_SIZE],
            const double (&stencil_phases)[STENCIL_SIZE], const double theta,
            size_t id_x, size_t id_y);
        __global__ void cuda_multipole_interactions_kernel_non_rho(
            const double (&local_monopoles)[NUMBER_LOCAL_MONOPOLE_VALUES],
            const double (&center_of_masses)[NUMBER_MASS_VALUES],
            const double (&multipoles)[NUMBER_LOCAL_EXPANSION_VALUES],
            double (&potential_expansions)[3 * NUMBER_POT_EXPANSIONS],
            const octotiger::fmm::multiindex<> (&stencil)[STENCIL_SIZE],
            const double (&stencil_phases)[STENCIL_SIZE], const double theta,
            size_t id_x, size_t id_y);
        __global__ void cuda_add_multipole_pot_blocks(
            double (&potential_expansions)[3 * NUMBER_POT_EXPANSIONS]);
        __global__ void cuda_add_multipole_ang_blocks(
            double (&angular_corrections)[3 * NUMBER_ANG_CORRECTIONS]);
    }    // namespace multipole_interactions
}    // namespace fmm
}    // namespace octotiger
#endif
