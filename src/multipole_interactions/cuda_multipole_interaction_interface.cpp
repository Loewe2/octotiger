#ifdef OCTOTIGER_HAVE_CUDA
#include "octotiger/multipole_interactions/cuda_multipole_interaction_interface.hpp"
#include "octotiger/multipole_interactions/calculate_stencil.hpp"
#include "octotiger/multipole_interactions/multipole_cuda_kernel.hpp"

#include "octotiger/options.hpp"
#include "octotiger/defs.hpp"

#include <array>
#include <vector>

namespace octotiger {
namespace fmm {
    namespace multipole_interactions {
        thread_local size_t cuda_multipole_interaction_interface::cpu_launch_counter = 0;
        thread_local size_t cuda_multipole_interaction_interface::cuda_launch_counter = 0;
        thread_local size_t cuda_multipole_interaction_interface::cpu_launch_counter_non_rho = 0;
        thread_local size_t cuda_multipole_interaction_interface::cuda_launch_counter_non_rho = 0;


        cuda_multipole_interaction_interface::cuda_multipole_interaction_interface(void)
          : multipole_interaction_interface()
          , theta(opts().theta) {
        }

        void cuda_multipole_interaction_interface::compute_multipole_interactions(
            std::vector<real>& monopoles, std::vector<multipole>& M_ptr,
            std::vector<std::shared_ptr<std::vector<space_vector>>>& com_ptr,
            std::vector<neighbor_gravity_type>& neighbors, gsolve_type type, real dx,
            std::array<bool, geo::direction::count()>& is_direction_empty,
            std::array<real, NDIM> xbase) {
            kernel_scheduler::scheduler.init();
            // Check where we want to run this:
            int slot = kernel_scheduler::scheduler.get_launch_slot();
            if (slot == -1 || m2m_type == interaction_kernel_type::OLD) {    // Run fallback cpu implementation
                if (type == RHO)
                    cpu_launch_counter++;
                else
                    cpu_launch_counter_non_rho++;
                multipole_interaction_interface::compute_multipole_interactions(
                    monopoles, M_ptr, com_ptr, neighbors, type, dx, is_direction_empty, xbase);
            } else {    // run on cuda device
                if (type == RHO)
                    cuda_launch_counter++;
                else
                    cuda_launch_counter_non_rho++;

                auto kernel_launcher = [this, &monopoles, &neighbors, &type, dx, &is_direction_empty]
                       (std::vector<real, cuda_pinned_allocator<real>>& local_monopoles,
                        struct_of_array_data<expansion, real, 20, ENTRIES, SOA_PADDING,
                        std::vector<real, cuda_pinned_allocator<real>>>& local_expansions_SoA,
                        struct_of_array_data<space_vector, real, 3, ENTRIES, SOA_PADDING,
                        std::vector<real, cuda_pinned_allocator<real>>>& center_of_masses_SoA,
                        kernel_device_enviroment &env, util::cuda_helper& gpu_interface) {
                    // Move data into SoA arrays
                    update_input(monopoles, M_ptr, com_ptr, neighbors, type, dx, xbase,
                        local_monopoles, local_expansions_SoA,
                        center_of_masses_SoA);

                    // Queue moving of input data to device
                    gpu_interface.copy_async(env.device_local_monopoles,
                        staging_area.local_monopoles.data(), local_monopoles_size,
                        cudaMemcpyHostToDevice);
                    gpu_interface.copy_async(env.device_local_expansions,
                        staging_area.local_expansions_SoA.get_pod(), local_expansions_size,
                        cudaMemcpyHostToDevice);
                    gpu_interface.copy_async(env.device_center_of_masses,
                        staging_area.center_of_masses_SoA.get_pod(), center_of_masses_size,
                        cudaMemcpyHostToDevice);

                    // Launch kernel and queue copying of results
                    const dim3 grid_spec(INX, 1, 1);
                    const dim3 threads_per_block(1, INX, INX);
                    if (type == RHO) {
                        bool second_phase = false;
                        void* args[] = {&(env.device_local_monopoles), &(env.device_center_of_masses),
                            &(env.device_local_expansions), &(env.device_potential_expansions),
                            &(env.device_angular_corrections), &theta, &second_phase};
                        gpu_interface.execute((const void*)&cuda_multipole_interactions_kernel_rho, grid_spec,
                        threads_per_block, args, 0);
                        gpu_interface.copy_async(angular_corrections_SoA.get_pod(),
                                                env.device_angular_corrections, angular_corrections_size,
                                                cudaMemcpyDeviceToHost);

                    } else {
                        bool second_phase = false;
                        void* args[] = {&(env.device_local_monopoles), &(env.device_center_of_masses),
                            &(env.device_local_expansions), &(env.device_potential_expansions),
                            &theta, &second_phase};
                        gpu_interface.execute((const void*)&cuda_multipole_interactions_kernel_non_rho, grid_spec,
                            threads_per_block, args, 0);
                    }
                    gpu_interface.copy_async(potential_expansions_SoA.get_pod(),
                        env.device_potential_expansions, potential_expansions_size,
                        cudaMemcpyDeviceToHost);
                }
                auto fut = kernel_scheduler::scheduler.launch(slot, kernel_launcher);
                fut.get();

                // Copy results back into non-SoA array
                potential_expansions_SoA.add_to_non_SoA(grid_ptr->get_L());
                if (type == RHO)
                    angular_corrections_SoA.to_non_SoA(grid_ptr->get_L_c());
            }
        }

    }    // namespace multipole_interactions
}    // namespace fmm
}    // namespace octotiger
#endif
