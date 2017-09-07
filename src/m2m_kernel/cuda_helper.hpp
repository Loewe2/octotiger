#pragma once

#include <iostream>
// #include "struct_of_array_data.hpp"

namespace octotiger {
namespace fmm {
    namespace cuda {

        inline void check_cuda_error(cudaError_t error, const char* file, int line) {
            if (error != cudaSuccess) {
                std::cerr << "CUDA error: " << cudaGetErrorString(error) << " file: " << file
                          << " line: " << line << std::endl;
            }
        }
    }
}
}

#define CUDA_CHECK_ERROR(error) \
    { octotiger::fmm::cuda::check_cuda_error((error), __FILE__, __LINE__); }
namespace octotiger {
namespace fmm {
    namespace cuda {

        template <typename cuda_component_type>
        class cuda_buffer
        {
        private:
            cuda_component_type* d_ptr;

        public:
            cuda_buffer(cuda_component_type* d_ptr_)
              : d_ptr(d_ptr_) {}
            ~cuda_buffer() {
                if (d_ptr != nullptr) {
                    CUDA_CHECK_ERROR(cudaFree(d_ptr));
                }
            }
            cuda_buffer(cuda_buffer<cuda_component_type>&& other) {
                d_ptr = other.d_ptr;
                other.d_ptr = nullptr;
            };
            cuda_buffer(const cuda_buffer<cuda_component_type>&) = delete;
            void operator=(const cuda_buffer<cuda_component_type>& other) = delete;

            cuda_component_type* get_device_ptr() {
                return d_ptr;
            }

            template <typename component_type, size_t num_components, size_t entries,
                size_t padding>
            void move_to_host(octotiger::fmm::struct_of_array_data<component_type, num_components,
                entries, padding>& host_buffer) {
                size_t SoA_bytes = host_buffer.get_size_bytes();
                component_type* h_SoA = host_buffer.get_underlying_pointer();
                CUDA_CHECK_ERROR(cudaMemcpy(h_SoA, d_ptr, SoA_bytes, cudaMemcpyDeviceToHost));
            }
        };

        template <typename component_type, size_t num_components, size_t entries, size_t padding>
        cuda_buffer<component_type> move_to_device(
            octotiger::fmm::struct_of_array_data<component_type, num_components, entries, padding>&
                data) {
            component_type* d_SoA;
            size_t SoA_bytes = data.get_size_bytes();
            std::cout << "bytes: " << SoA_bytes << std::endl;
            CUDA_CHECK_ERROR(cudaMalloc(&d_SoA, SoA_bytes));
            component_type* h_SoA = data.get_underlying_pointer();
            CUDA_CHECK_ERROR(cudaMemcpy(d_SoA, h_SoA, SoA_bytes, cudaMemcpyHostToDevice));
            return cuda_buffer<component_type>(d_SoA);
        }

        template <typename component_type>
        cuda_buffer<component_type> move_to_device(std::vector<component_type>& data) {
            component_type* d_SoA;
            size_t SoA_bytes = data.size() * sizeof(component_type);
            CUDA_CHECK_ERROR(cudaMalloc(&d_SoA, SoA_bytes));
            component_type* h_SoA = data.data();
            CUDA_CHECK_ERROR(cudaMemcpy(d_SoA, h_SoA, SoA_bytes, cudaMemcpyHostToDevice));
            return cuda_buffer<component_type>(d_SoA);
        }

        template <typename component_type, size_t N>
        cuda_buffer<component_type> move_to_device(std::array<component_type, N>& data) {
            component_type* d_SoA;
            size_t SoA_bytes = N * sizeof(component_type);
            CUDA_CHECK_ERROR(cudaMalloc(&d_SoA, SoA_bytes));
            component_type* h_SoA = data.data();
            CUDA_CHECK_ERROR(cudaMemcpy(d_SoA, h_SoA, SoA_bytes, cudaMemcpyHostToDevice));
            return cuda_buffer<component_type>(d_SoA);
        }

        // component type has to be indexable, and has to have a size() operator
        template <typename component_type, size_t num_components, size_t entries, size_t padding>
        class struct_of_array_data
        {
        private:
            // data in SoA form
            static constexpr size_t padded_entries_per_component = entries + padding;

            component_type* const data;

            struct_of_array_data()
              : data(nullptr) {}

        public:
            template <size_t component_access>
            __device__ inline component_type* pointer(const size_t flat_index) const {
                constexpr size_t component_array_offset =
                    component_access * padded_entries_per_component;
                // should result in single move instruction, indirect addressing: reg + reg +
                // constant
                return data + flat_index + component_array_offset;
            }

            // careful, this returns a copy!
            template <size_t component_access>
            __device__ inline double value(const size_t flat_index) const {
                return *this->pointer<component_access>(flat_index);
            }

            template <typename AoS_type>
            __device__ static struct_of_array_data<component_type, num_components, entries, padding>
            from_vector(const std::vector<AoS_type>& org) {
                component_type* data =
                    new component_type[num_components * padded_entries_per_component];
                for (size_t component = 0; component < num_components; component++) {
                    for (size_t entry = 0; entry < org.size(); entry++) {
                        data[component * padded_entries_per_component + entry] =
                            org[entry][component];
                    }
                }
                return struct_of_array_data(data);
            }

            // constructor that works on preallocated and initialized data
            __device__ struct_of_array_data(component_type* preallocated_data)
              : data(preallocated_data) {}

            __device__ struct_of_array_data(const struct_of_array_data& other) = delete;

            __device__ struct_of_array_data(const struct_of_array_data&& other) = delete;

            __device__ struct_of_array_data& operator=(const struct_of_array_data& other) = delete;

            // write back into non-SoA style array
            template <typename AoS_type>
            __device__ void to_non_SoA(std::vector<AoS_type>& org) {
                // constexpr size_t padded_entries_per_component = entries + padding;
                for (size_t component = 0; component < num_components; component++) {
                    for (size_t entry = 0; entry < org.size(); entry++) {
                        org[entry][component] =
                            data[component * padded_entries_per_component + entry];
                    }
                }
            }

            __device__ double* get_underlying_pointer() {
                return data;
            }

            __device__ size_t get_size_bytes() {
                return num_components * padded_entries_per_component * sizeof(component_type);
            }

            __device__ inline size_t size() {
                return num_components * padded_entries_per_component;
            }

            __device__ component_type& at(size_t component_access, size_t flat_index) {
                size_t component_array_offset = component_access * padded_entries_per_component;
                return *(data + flat_index + component_array_offset);
            }
        };
    }
}
}
