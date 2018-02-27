#include "m2m_cuda.hpp"

#include <array>
#include <vector>
#include "common_kernel/multiindex.hpp"
#include "common_kernel/struct_of_array_data.hpp"
#include "cuda_helper.hpp"

namespace octotiger {
namespace fmm {

    template <>
    __device__ void multiindex<int32_t>::transform_coarse() {
        // const T patch_size = static_cast<typename T::value_type>(INX);
        // const T subtract = static_cast<typename T::value_type>(INX / 2);
        const int32_t patch_size = static_cast<int32_t>(INX);
        const int32_t subtract = static_cast<int32_t>(INX / 2);

        x = ((x + patch_size) >> 1) - subtract;
        y = ((y + patch_size) >> 1) - subtract;
        z = ((z + patch_size) >> 1) - subtract;
    }

    extern std::vector<octotiger::fmm::multiindex<>> stencil;

    namespace cuda {

        template <typename T>
        __device__ constexpr inline T sqr(T const& val) {
            return val * val;
        }

        class D_split
        {
        public:
            double (&X)[NDIM];

            double X_00;
            double X_11;
            double X_22;

            double r2;
            double r2inv;

            double d0;
            double d1;
            double d2;
            double d3;

        public:
            __device__ D_split(double (&X)[NDIM])
              : X(X) {
                X_00 = X[0] * X[0];
                X_11 = X[1] * X[1];
                X_22 = X[2] * X[2];

                r2 = X_00 + X_11 + X_22;
                r2inv = 1.0 / max(r2, 1.0e-20);

                // parts of formula (6)
                d0 = -sqrt(r2inv);
                // parts of formula (7)
                d1 = -d0 * r2inv;
                // parts of formula (8)
                d2 = -3.0 * d1 * r2inv;
                // parts of  formula (9)
                d3 = -5.0 * d2 * r2inv;
            }

            // overload for kernel-specific simd type
            __device__ inline void calculate_D_lower(double (&A)[20]) {
                // formula (6)
                A[0] = d0;

                A[1] = X[0] * d1;
                A[2] = X[1] * d1;
                A[3] = X[2] * d1;

                double X_12 = X[1] * X[2];
                double X_01 = X[0] * X[1];
                double X_02 = X[0] * X[2];

                A[4] = d2 * X_00;
                A[4] += d1;
                A[5] = d2 * X_01;
                A[6] = d2 * X_02;

                A[7] = d2 * X_11;
                A[7] += d1;
                A[8] = d2 * X_12;

                A[9] = d2 * X_22;
                A[9] += d1;

                A[10] = d3 * X_00 * X[0];
                double d2_X0 = d2 * X[0];
                A[10] += 3.0 * d2_X0;
                A[11] = d3 * X_00 * X[1];
                A[11] += d2 * X[1];
                A[12] = d3 * X_00 * X[2];
                A[12] += d2 * X[2];

                A[13] = d3 * X[0] * X_11;
                A[13] += d2 * X[0];
                A[14] = d3 * X[0] * X_12;

                A[15] = d3 * X[0] * X_22;
                A[15] += d2_X0;

                A[16] = d3 * X_11 * X[1];
                double d2_X1 = d2 * X[1];
                A[16] += 3.0 * d2_X1;

                A[17] = d3 * X_11 * X[2];
                A[17] += d2 * X[2];

                A[18] = d3 * X[1] * X_22;
                A[18] += d2 * X[1];

                A[19] = d3 * X_22 * X[2];
                double d2_X2 = d2 * X[2];
                A[19] += 3.0 * d2_X2;
            }
        };

        namespace detail {
            __device__ inline int32_t distance_squared_reciprocal(
                const multiindex<>& i, const multiindex<>& j) {
                return (sqr(i.x - j.x) + sqr(i.y - j.y) + sqr(i.z - j.z));
            }
        }

        __global__ void blocked_interaction_rho(double* local_expansions_SoA_ptr,
            double* center_of_masses_SoA_ptr, double* potential_expansions_SoA_ptr,
            double* angular_corrections_SoA_ptr, octotiger::fmm::multiindex<>* stencil,
            size_t stencil_elements, double theta, double* factor, double* factor_half,
            double* factor_sixth) {
            // TODO: make sure you pick up the right expansion type
            octotiger::fmm::cuda::struct_of_array_data<double, 20, octotiger::fmm::ENTRIES,
                octotiger::fmm::SOA_PADDING>
                local_expansions_SoA(local_expansions_SoA_ptr);
            octotiger::fmm::cuda::struct_of_array_data<double, 3, octotiger::fmm::ENTRIES,
                octotiger::fmm::SOA_PADDING>
                center_of_masses_SoA(center_of_masses_SoA_ptr);
            octotiger::fmm::cuda::struct_of_array_data<double, 20, octotiger::fmm::INNER_CELLS,
                octotiger::fmm::SOA_PADDING>
                potential_expansions_SoA(potential_expansions_SoA_ptr);
            octotiger::fmm::cuda::struct_of_array_data<double, 3, octotiger::fmm::INNER_CELLS,
                octotiger::fmm::SOA_PADDING>
                angular_corrections_SoA(angular_corrections_SoA_ptr);

            const double theta_rec_squared = sqr(1.0 / theta);

            octotiger::fmm::multiindex<> cell_index(threadIdx.x + INNER_CELLS_PADDING_DEPTH,
                threadIdx.y + INNER_CELLS_PADDING_DEPTH, threadIdx.z + INNER_CELLS_PADDING_DEPTH);

            octotiger::fmm::multiindex<> cell_index_coarse(cell_index);
            cell_index_coarse.transform_coarse();

            size_t cell_flat_index = octotiger::fmm::to_flat_index_padded(cell_index);
            octotiger::fmm::multiindex<> cell_index_unpadded(threadIdx.x, threadIdx.y, threadIdx.z);
            size_t cell_flat_index_unpadded =
                octotiger::fmm::to_inner_flat_index_not_padded(cell_index_unpadded);

            double X[NDIM];
            X[0] = center_of_masses_SoA.value<0>(cell_flat_index);
            X[1] = center_of_masses_SoA.value<1>(cell_flat_index);
            X[2] = center_of_masses_SoA.value<2>(cell_flat_index);
            double Y[NDIM];
            double dX[NDIM];
            double m_partner[20];
            double tmp[20];
            tmp[0] = *potential_expansions_SoA.pointer<0>(cell_flat_index_unpadded);
            tmp[1] = *potential_expansions_SoA.pointer<1>(cell_flat_index_unpadded);
            tmp[2] = *potential_expansions_SoA.pointer<2>(cell_flat_index_unpadded);
            tmp[3] = *potential_expansions_SoA.pointer<3>(cell_flat_index_unpadded);
            tmp[4] = *potential_expansions_SoA.pointer<4>(cell_flat_index_unpadded);
            tmp[5] = *potential_expansions_SoA.pointer<5>(cell_flat_index_unpadded);
            tmp[6] = *potential_expansions_SoA.pointer<6>(cell_flat_index_unpadded);
            tmp[7] = *potential_expansions_SoA.pointer<7>(cell_flat_index_unpadded);
            tmp[8] = *potential_expansions_SoA.pointer<8>(cell_flat_index_unpadded);
            tmp[9] = *potential_expansions_SoA.pointer<9>(cell_flat_index_unpadded);
            tmp[10] = *potential_expansions_SoA.pointer<10>(cell_flat_index_unpadded);
            tmp[11] = *potential_expansions_SoA.pointer<11>(cell_flat_index_unpadded);
            tmp[12] = *potential_expansions_SoA.pointer<12>(cell_flat_index_unpadded);
            tmp[13] = *potential_expansions_SoA.pointer<13>(cell_flat_index_unpadded);
            tmp[14] = *potential_expansions_SoA.pointer<14>(cell_flat_index_unpadded);
            tmp[15] = *potential_expansions_SoA.pointer<15>(cell_flat_index_unpadded);
            tmp[16] = *potential_expansions_SoA.pointer<16>(cell_flat_index_unpadded);
            tmp[17] = *potential_expansions_SoA.pointer<17>(cell_flat_index_unpadded);
            tmp[18] = *potential_expansions_SoA.pointer<18>(cell_flat_index_unpadded);
            tmp[19] = *potential_expansions_SoA.pointer<19>(cell_flat_index_unpadded);

            double tmp_cor[NDIM];
            tmp_cor[0] = *angular_corrections_SoA.value<0>(cell_flat_index_unpadded);
            tmp_cor[1] = *angular_corrections_SoA.value<1>(cell_flat_index_unpadded);
            tmp_cor[2] = *angular_corrections_SoA.value<2>(cell_flat_index_unpadded);

            for (size_t i = 0; i < stencil_elements; i++) {
                octotiger::fmm::multiindex<>& cur_stencil = stencil[i];
                const octotiger::fmm::multiindex<> interaction_partner_index(
                    cell_index.x + cur_stencil.x, cell_index.y + cur_stencil.y,
                    cell_index.z + cur_stencil.z);

                const size_t interaction_partner_flat_index =
                    octotiger::fmm::to_flat_index_padded(interaction_partner_index);

                // implicitly broadcasts to vector
                octotiger::fmm::multiindex<> interaction_partner_index_coarse(
                    interaction_partner_index);
                // note that this is the same for groups of 2x2x2 elements
                // -> maps to the same for some SIMD lanes
                interaction_partner_index_coarse.transform_coarse();

                const double theta_c_rec_squared =
                    static_cast<double>(detail::distance_squared_reciprocal(
                        cell_index_coarse, interaction_partner_index_coarse));

                const bool mask_b = theta_rec_squared > theta_c_rec_squared;
                const double mask = mask_b ? 1.0 : 0.0;

                Y[0] = center_of_masses_SoA.value<0>(interaction_partner_flat_index);
                Y[1] = center_of_masses_SoA.value<1>(interaction_partner_flat_index);
                Y[2] = center_of_masses_SoA.value<2>(interaction_partner_flat_index);

                dX[0] = X[0] - Y[0];
                dX[1] = X[1] - Y[1];
                dX[2] = X[2] - Y[2];

                m_partner[0] = local_expansions_SoA.value<0>(interaction_partner_flat_index);
                m_partner[1] = local_expansions_SoA.value<1>(interaction_partner_flat_index);
                m_partner[2] = local_expansions_SoA.value<2>(interaction_partner_flat_index);
                m_partner[3] = local_expansions_SoA.value<3>(interaction_partner_flat_index);
                m_partner[4] = local_expansions_SoA.value<4>(interaction_partner_flat_index);
                m_partner[5] = local_expansions_SoA.value<5>(interaction_partner_flat_index);
                m_partner[6] = local_expansions_SoA.value<6>(interaction_partner_flat_index);
                m_partner[7] = local_expansions_SoA.value<7>(interaction_partner_flat_index);
                m_partner[8] = local_expansions_SoA.value<8>(interaction_partner_flat_index);
                m_partner[9] = local_expansions_SoA.value<9>(interaction_partner_flat_index);
                m_partner[10] = local_expansions_SoA.value<10>(interaction_partner_flat_index);
                m_partner[11] = local_expansions_SoA.value<11>(interaction_partner_flat_index);
                m_partner[12] = local_expansions_SoA.value<12>(interaction_partner_flat_index);
                m_partner[13] = local_expansions_SoA.value<13>(interaction_partner_flat_index);
                m_partner[14] = local_expansions_SoA.value<14>(interaction_partner_flat_index);
                m_partner[15] = local_expansions_SoA.value<15>(interaction_partner_flat_index);
                m_partner[16] = local_expansions_SoA.value<16>(interaction_partner_flat_index);
                m_partner[17] = local_expansions_SoA.value<17>(interaction_partner_flat_index);
                m_partner[18] = local_expansions_SoA.value<18>(interaction_partner_flat_index);
                m_partner[19] = local_expansions_SoA.value<19>(interaction_partner_flat_index);

                // R_i in paper is the dX in the code
                // D is taylor expansion value for a given X expansion of the gravitational
                // potential
                // (multipole expansion)

                // calculates all D-values, calculate all coefficients of 1/r (not the
                // potential),
                // formula (6)-(9) and (19)
                D_split D_calculator(dX);
                // TODO: this translates to an initialization loop, bug!
                double D_lower[20];
                D_calculator.calculate_D_lower(D_lower);

                // expansion_v current_potential;
                // TODO: this translates to an initialization loop, bug!
                double cur_pot[10];

                // 10-19 are not cached!

                // the following loops calculate formula (10), potential from B->A

                cur_pot[0] = m_partner[0] * D_lower[0];
                cur_pot[1] = m_partner[0] * D_lower[1];
                cur_pot[2] = m_partner[0] * D_lower[2];
                cur_pot[3] = m_partner[0] * D_lower[3];

                cur_pot[0] += m_partner[4] * (D_lower[4] * factor_half[4]);
                cur_pot[1] += m_partner[4] * (D_lower[10] * factor_half[4]);
                cur_pot[2] += m_partner[4] * (D_lower[11] * factor_half[4]);
                cur_pot[3] += m_partner[4] * (D_lower[12] * factor_half[4]);

                cur_pot[0] += m_partner[5] * (D_lower[5] * factor_half[5]);
                cur_pot[1] += m_partner[5] * (D_lower[11] * factor_half[5]);
                cur_pot[2] += m_partner[5] * (D_lower[13] * factor_half[5]);
                cur_pot[3] += m_partner[5] * (D_lower[14] * factor_half[5]);

                cur_pot[0] += m_partner[6] * (D_lower[6] * factor_half[6]);
                cur_pot[1] += m_partner[6] * (D_lower[12] * factor_half[6]);
                cur_pot[2] += m_partner[6] * (D_lower[14] * factor_half[6]);
                cur_pot[3] += m_partner[6] * (D_lower[15] * factor_half[6]);

                cur_pot[0] += m_partner[7] * (D_lower[7] * factor_half[7]);
                cur_pot[1] += m_partner[7] * (D_lower[13] * factor_half[7]);
                cur_pot[2] += m_partner[7] * (D_lower[16] * factor_half[7]);
                cur_pot[3] += m_partner[7] * (D_lower[17] * factor_half[7]);

                cur_pot[0] += m_partner[8] * (D_lower[8] * factor_half[8]);
                cur_pot[1] += m_partner[8] * (D_lower[14] * factor_half[8]);
                cur_pot[2] += m_partner[8] * (D_lower[17] * factor_half[8]);
                cur_pot[3] += m_partner[8] * (D_lower[18] * factor_half[8]);

                cur_pot[0] += m_partner[9] * (D_lower[9] * factor_half[9]);
                cur_pot[1] += m_partner[9] * (D_lower[15] * factor_half[9]);
                cur_pot[2] += m_partner[9] * (D_lower[18] * factor_half[9]);
                cur_pot[3] += m_partner[9] * (D_lower[19] * factor_half[9]);

                cur_pot[0] -= m_partner[10] * (D_lower[10] * factor_sixth[10]);
                cur_pot[0] -= m_partner[11] * (D_lower[11] * factor_sixth[11]);
                cur_pot[0] -= m_partner[12] * (D_lower[12] * factor_sixth[12]);
                cur_pot[0] -= m_partner[13] * (D_lower[13] * factor_sixth[13]);
                cur_pot[0] -= m_partner[14] * (D_lower[14] * factor_sixth[14]);
                cur_pot[0] -= m_partner[15] * (D_lower[15] * factor_sixth[15]);
                cur_pot[0] -= m_partner[16] * (D_lower[16] * factor_sixth[16]);
                cur_pot[0] -= m_partner[17] * (D_lower[17] * factor_sixth[17]);
                cur_pot[0] -= m_partner[18] * (D_lower[18] * factor_sixth[18]);
                cur_pot[0] -= m_partner[19] * (D_lower[19] * factor_sixth[19]);

                cur_pot[4] = m_partner[0] * D_lower[4];
                cur_pot[5] = m_partner[0] * D_lower[5];
                cur_pot[6] = m_partner[0] * D_lower[6];
                cur_pot[7] = m_partner[0] * D_lower[7];
                cur_pot[8] = m_partner[0] * D_lower[8];
                cur_pot[9] = m_partner[0] * D_lower[9];

                cur_pot[4] -= m_partner[1] * D_lower[10];
                cur_pot[5] -= m_partner[1] * D_lower[11];
                cur_pot[6] -= m_partner[1] * D_lower[12];
                cur_pot[7] -= m_partner[1] * D_lower[13];
                cur_pot[8] -= m_partner[1] * D_lower[14];
                cur_pot[9] -= m_partner[1] * D_lower[15];

                cur_pot[4] -= m_partner[2] * D_lower[11];
                cur_pot[5] -= m_partner[2] * D_lower[13];
                cur_pot[6] -= m_partner[2] * D_lower[14];
                cur_pot[7] -= m_partner[2] * D_lower[16];
                cur_pot[8] -= m_partner[2] * D_lower[17];
                cur_pot[9] -= m_partner[2] * D_lower[18];

                cur_pot[4] -= m_partner[3] * D_lower[12];
                cur_pot[5] -= m_partner[3] * D_lower[14];
                cur_pot[6] -= m_partner[3] * D_lower[15];
                cur_pot[7] -= m_partner[3] * D_lower[17];
                cur_pot[8] -= m_partner[3] * D_lower[18];
                cur_pot[9] -= m_partner[3] * D_lower[19];

                tmp[0] +=
                    (potential_expansions_SoA.value<0>(cell_flat_index_unpadded) + cur_pot[0]) *
                    mask;
                tmp[1] +=
                    (potential_expansions_SoA.value<1>(cell_flat_index_unpadded) + cur_pot[1]) *
                    mask;
                tmp[2] +=
                    (potential_expansions_SoA.value<2>(cell_flat_index_unpadded) + cur_pot[2]) *
                    mask;
                tmp[3] +=
                    (potential_expansions_SoA.value<3>(cell_flat_index_unpadded) + cur_pot[3]) *
                    mask;
                tmp[4] +=
                    (potential_expansions_SoA.value<4>(cell_flat_index_unpadded) + cur_pot[4]) *
                    mask;
                tmp[5] +=
                    (potential_expansions_SoA.value<5>(cell_flat_index_unpadded) + cur_pot[5]) *
                    mask;
                tmp[6] +=
                    (potential_expansions_SoA.value<6>(cell_flat_index_unpadded) + cur_pot[6]) *
                    mask;
                tmp[7] +=
                    (potential_expansions_SoA.value<7>(cell_flat_index_unpadded) + cur_pot[7]) *
                    mask;
                tmp[8] +=
                    (potential_expansions_SoA.value<8>(cell_flat_index_unpadded) + cur_pot[8]) *
                    mask;
                tmp[9] +=
                    (potential_expansions_SoA.value<9>(cell_flat_index_unpadded) + cur_pot[9]) *
                    mask;
                tmp[10] += (potential_expansions_SoA.value<10>(cell_flat_index_unpadded) +
                               m_partner[0] * D_lower[10]) *
                    mask;
                tmp[11] += (potential_expansions_SoA.value<11>(cell_flat_index_unpadded) +
                               m_partner[0] * D_lower[11]) *
                    mask;
                tmp[12] += (potential_expansions_SoA.value<12>(cell_flat_index_unpadded) +
                               m_partner[0] * D_lower[12]) *
                    mask;
                tmp[13] += (potential_expansions_SoA.value<13>(cell_flat_index_unpadded) +
                               m_partner[0] * D_lower[13]) *
                    mask;
                tmp[14] += (potential_expansions_SoA.value<14>(cell_flat_index_unpadded) +
                               m_partner[0] * D_lower[14]) *
                    mask;
                tmp[15] += (potential_expansions_SoA.value<15>(cell_flat_index_unpadded) +
                               m_partner[0] * D_lower[15]) *
                    mask;
                tmp[16] += (potential_expansions_SoA.value<16>(cell_flat_index_unpadded) +
                               m_partner[0] * D_lower[16]) *
                    mask;
                tmp[17] += (potential_expansions_SoA.value<17>(cell_flat_index_unpadded) +
                               m_partner[0] * D_lower[17]) *
                    mask;
                tmp[18] += (potential_expansions_SoA.value<18>(cell_flat_index_unpadded) +
                               m_partner[0] * D_lower[18]) *
                    mask;
                tmp[19] += (potential_expansions_SoA.value<19>(cell_flat_index_unpadded) +
                               m_partner[0] * D_lower[19]) *
                    mask;

                ////////////// angular momentum correction, if enabled /////////////////////////

                // Initalize moments and momentum
                // this branch computes the angular momentum correction, (20) in the
                // paper divide by mass of other cell
                // struct_of_array_iterator<expansion, double, 20> m_cell_iterator(
                //     local_expansions_SoA, cell_flat_index);

                double const n0_constant =
                    m_partner[0] / local_expansions_SoA.value<0>(cell_flat_index);

                // m_cell_iterator.increment(10);

                // calculating the coefficients for formula (M are the octopole moments)
                // the coefficients are calculated in (17) and (18)
                // struct_of_array_iterator<space_vector, double, 3>
                // current_angular_correction_result(
                //     angular_corrections_SoA, cell_flat_index_unpadded);

                double D_upper[15];
                // D_calculator.calculate_D_upper(D_upper);

                double current_angular_correction[NDIM];
                current_angular_correction[0] = 0.0;
                current_angular_correction[1] = 0.0;
                current_angular_correction[2] = 0.0;

                D_upper[0] =
                    D_calculator.X[0] * D_calculator.X[0] * D_calculator.d3 + 2.0 * D_calculator.d2;
                double d3_X00 = D_calculator.d3 * D_calculator.X_00;
                D_upper[0] += D_calculator.d2;
                D_upper[0] += 5.0 * d3_X00;
                double d3_X01 = D_calculator.d3 * D_calculator.X[0] * D_calculator.X[1];
                D_upper[1] = 3.0 * d3_X01;
                double X_02 = D_calculator.X[0] * D_calculator.X[2];
                double d3_X02 = D_calculator.d3 * X_02;
                D_upper[2] = 3.0 * d3_X02;
                double n0_tmp =
                    m_partner[10] - local_expansions_SoA.value<10>(cell_flat_index) * n0_constant;
                // m_cell_iterator++; // 11

                current_angular_correction[0] -= n0_tmp * (D_upper[0] * factor_sixth[10]);
                current_angular_correction[1] -= n0_tmp * (D_upper[1] * factor_sixth[10]);
                current_angular_correction[2] -= n0_tmp * (D_upper[2] * factor_sixth[10]);

                D_upper[3] = D_calculator.d2;
                double d3_X11 = D_calculator.d3 * D_calculator.X_11;
                D_upper[3] += d3_X11;
                D_upper[3] += D_calculator.d3 * D_calculator.X_00;
                double d3_X12 = D_calculator.d3 * D_calculator.X[1] * D_calculator.X[2];
                D_upper[4] = d3_X12;

                n0_tmp =
                    m_partner[11] - local_expansions_SoA.value<11>(cell_flat_index) * n0_constant;

                current_angular_correction[0] -= n0_tmp * (D_upper[1] * factor_sixth[11]);
                current_angular_correction[1] -= n0_tmp * (D_upper[3] * factor_sixth[11]);
                current_angular_correction[2] -= n0_tmp * (D_upper[4] * factor_sixth[11]);

                D_upper[5] = D_calculator.d2;
                double d3_X22 = D_calculator.d3 * D_calculator.X_22;
                D_upper[5] += d3_X22;
                D_upper[5] += d3_X00;

                n0_tmp =
                    m_partner[12] - local_expansions_SoA.value<12>(cell_flat_index) * n0_constant;

                current_angular_correction[0] -= n0_tmp * (D_upper[2] * factor_sixth[12]);
                current_angular_correction[1] -= n0_tmp * (D_upper[4] * factor_sixth[12]);
                current_angular_correction[2] -= n0_tmp * (D_upper[5] * factor_sixth[12]);

                D_upper[6] = 3.0 * d3_X01;
                D_upper[7] = D_calculator.d3 * X_02;

                n0_tmp =
                    m_partner[13] - local_expansions_SoA.value<13>(cell_flat_index) * n0_constant;

                current_angular_correction[0] -= n0_tmp * (D_upper[3] * factor_sixth[13]);
                current_angular_correction[1] -= n0_tmp * (D_upper[6] * factor_sixth[13]);
                current_angular_correction[2] -= n0_tmp * (D_upper[7] * factor_sixth[13]);

                D_upper[8] = D_calculator.d3 * D_calculator.X[0] * D_calculator.X[1];

                n0_tmp =
                    m_partner[14] - local_expansions_SoA.value<14>(cell_flat_index) * n0_constant;

                current_angular_correction[0] -= n0_tmp * (D_upper[4] * factor_sixth[14]);
                current_angular_correction[1] -= n0_tmp * (D_upper[7] * factor_sixth[14]);
                current_angular_correction[2] -= n0_tmp * (D_upper[8] * factor_sixth[14]);

                D_upper[9] = 3.0 * d3_X02;

                n0_tmp =
                    m_partner[15] - local_expansions_SoA.value<15>(cell_flat_index) * n0_constant;

                current_angular_correction[0] -= n0_tmp * (D_upper[5] * factor_sixth[15]);
                current_angular_correction[1] -= n0_tmp * (D_upper[8] * factor_sixth[15]);
                current_angular_correction[2] -= n0_tmp * (D_upper[9] * factor_sixth[15]);

                D_upper[10] =
                    D_calculator.X[1] * D_calculator.X[1] * D_calculator.d3 + 2.0 * D_calculator.d2;
                D_upper[10] += D_calculator.d2;
                D_upper[10] += 5.0 * d3_X11;

                D_upper[11] = 3.0 * d3_X12;

                n0_tmp =
                    m_partner[16] - local_expansions_SoA.value<16>(cell_flat_index) * n0_constant;

                current_angular_correction[0] -= n0_tmp * (D_upper[6] * factor_sixth[16]);
                current_angular_correction[1] -= n0_tmp * (D_upper[10] * factor_sixth[16]);
                current_angular_correction[2] -= n0_tmp * (D_upper[11] * factor_sixth[16]);

                D_upper[12] = D_calculator.d2;
                D_upper[12] += d3_X22;
                D_upper[12] += d3_X11;

                n0_tmp =
                    m_partner[17] - local_expansions_SoA.value<17>(cell_flat_index) * n0_constant;

                current_angular_correction[0] -= n0_tmp * (D_upper[7] * factor_sixth[17]);
                current_angular_correction[1] -= n0_tmp * (D_upper[11] * factor_sixth[17]);
                current_angular_correction[2] -= n0_tmp * (D_upper[12] * factor_sixth[17]);

                D_upper[13] = 3.0 * d3_X12;

                n0_tmp =
                    m_partner[18] - local_expansions_SoA.value<18>(cell_flat_index) * n0_constant;

                current_angular_correction[0] -= n0_tmp * (D_upper[8] * factor_sixth[18]);
                current_angular_correction[1] -= n0_tmp * (D_upper[12] * factor_sixth[18]);
                current_angular_correction[2] -= n0_tmp * (D_upper[13] * factor_sixth[18]);

                D_upper[14] =
                    D_calculator.X[2] * D_calculator.X[2] * D_calculator.d3 + 2.0 * D_calculator.d2;
                D_upper[14] += D_calculator.d2;
                D_upper[14] += 5.0 * d3_X22;

                n0_tmp =
                    m_partner[19] - local_expansions_SoA.value<19>(cell_flat_index) * n0_constant;

                current_angular_correction[0] -= n0_tmp * (D_upper[9] * factor_sixth[19]);
                current_angular_correction[1] -= n0_tmp * (D_upper[13] * factor_sixth[19]);
                current_angular_correction[2] -= n0_tmp * (D_upper[14] * factor_sixth[19]);

                tmp_cor[0] += (angular_corrections_SoA.value<0>(cell_flat_index_unpadded) +
                                 current_angular_correction[0]) *
                    mask;

                tmp_cor[1] += (angular_corrections_SoA.value<1>(cell_flat_index_unpadded) +
                                 current_angular_correction[1]) *
                    mask;

                tmp_cor[2] += (angular_corrections_SoA.value<2>(cell_flat_index_unpadded) +
                                 current_angular_correction[2]) *
                    mask;
            }
            *angular_corrections_SoA.pointer<0>(cell_flat_index_unpadded) = tmp_cor[0];
            *angular_corrections_SoA.pointer<1>(cell_flat_index_unpadded) = tmp_cor[1];
            *angular_corrections_SoA.pointer<2>(cell_flat_index_unpadded) = tmp_cor[2];

            *potential_expansions_SoA.pointer<0>(cell_flat_index_unpadded) = tmp[0];
            *potential_expansions_SoA.pointer<1>(cell_flat_index_unpadded) = tmp[1];
            *potential_expansions_SoA.pointer<2>(cell_flat_index_unpadded) = tmp[2];
            *potential_expansions_SoA.pointer<3>(cell_flat_index_unpadded) = tmp[3];
            *potential_expansions_SoA.pointer<4>(cell_flat_index_unpadded) = tmp[4];
            *potential_expansions_SoA.pointer<5>(cell_flat_index_unpadded) = tmp[5];
            *potential_expansions_SoA.pointer<6>(cell_flat_index_unpadded) = tmp[6];
            *potential_expansions_SoA.pointer<7>(cell_flat_index_unpadded) = tmp[7];
            *potential_expansions_SoA.pointer<8>(cell_flat_index_unpadded) = tmp[8];
            *potential_expansions_SoA.pointer<9>(cell_flat_index_unpadded) = tmp[9];
            *potential_expansions_SoA.pointer<10>(cell_flat_index_unpadded) = tmp[10];
            *potential_expansions_SoA.pointer<11>(cell_flat_index_unpadded) = tmp[11];
            *potential_expansions_SoA.pointer<12>(cell_flat_index_unpadded) = tmp[12];
            *potential_expansions_SoA.pointer<13>(cell_flat_index_unpadded) = tmp[13];
            *potential_expansions_SoA.pointer<14>(cell_flat_index_unpadded) = tmp[14];
            *potential_expansions_SoA.pointer<15>(cell_flat_index_unpadded) = tmp[15];
            *potential_expansions_SoA.pointer<16>(cell_flat_index_unpadded) = tmp[16];
            *potential_expansions_SoA.pointer<17>(cell_flat_index_unpadded) = tmp[17];
            *potential_expansions_SoA.pointer<18>(cell_flat_index_unpadded) = tmp[18];
            *potential_expansions_SoA.pointer<19>(cell_flat_index_unpadded) = tmp[19];
        }

        __global__ void blocked_interaction_non_rho(double* local_expansions_SoA_ptr,
            double* center_of_masses_SoA_ptr, double* potential_expansions_SoA_ptr,
            double* angular_corrections_SoA_ptr, octotiger::fmm::multiindex<>* stencil,
            size_t stencil_elements, double theta, double* factor, double* factor_half,
            double* factor_sixth) {
            const double theta_rec_squared = sqr(1.0 / theta);

            // TODO: make sure you pick up the right expansion type
            octotiger::fmm::cuda::struct_of_array_data<double, 20, octotiger::fmm::ENTRIES,
                octotiger::fmm::SOA_PADDING>
                local_expansions_SoA(local_expansions_SoA_ptr);
            octotiger::fmm::cuda::struct_of_array_data<double, 3, octotiger::fmm::ENTRIES,
                octotiger::fmm::SOA_PADDING>
                center_of_masses_SoA(center_of_masses_SoA_ptr);
            octotiger::fmm::cuda::struct_of_array_data<double, 20, octotiger::fmm::INNER_CELLS,
                octotiger::fmm::SOA_PADDING>
                potential_expansions_SoA(potential_expansions_SoA_ptr);
            octotiger::fmm::cuda::struct_of_array_data<double, 3, octotiger::fmm::INNER_CELLS,
                octotiger::fmm::SOA_PADDING>
                angular_corrections_SoA(angular_corrections_SoA_ptr);

            octotiger::fmm::multiindex<> cell_index(threadIdx.x + INNER_CELLS_PADDING_DEPTH,
                threadIdx.y + INNER_CELLS_PADDING_DEPTH, threadIdx.z + INNER_CELLS_PADDING_DEPTH);

            octotiger::fmm::multiindex<> cell_index_coarse(cell_index);
            cell_index_coarse.transform_coarse();

            size_t cell_flat_index = octotiger::fmm::to_flat_index_padded(cell_index);
            octotiger::fmm::multiindex<> cell_index_unpadded(threadIdx.x, threadIdx.y, threadIdx.z);
            size_t cell_flat_index_unpadded =
                octotiger::fmm::to_inner_flat_index_not_padded(cell_index_unpadded);

            double X[NDIM];
            X[0] = center_of_masses_SoA.value<0>(cell_flat_index);
            X[1] = center_of_masses_SoA.value<1>(cell_flat_index);
            X[2] = center_of_masses_SoA.value<2>(cell_flat_index);
            double Y[NDIM];
            double dX[NDIM];
            double m_partner[20];
            double tmp[20];
            tmp[0] = *potential_expansions_SoA.pointer<0>(cell_flat_index_unpadded);
            tmp[1] = *potential_expansions_SoA.pointer<1>(cell_flat_index_unpadded);
            tmp[2] = *potential_expansions_SoA.pointer<2>(cell_flat_index_unpadded);
            tmp[3] = *potential_expansions_SoA.pointer<3>(cell_flat_index_unpadded);
            tmp[4] = *potential_expansions_SoA.pointer<4>(cell_flat_index_unpadded);
            tmp[5] = *potential_expansions_SoA.pointer<5>(cell_flat_index_unpadded);
            tmp[6] = *potential_expansions_SoA.pointer<6>(cell_flat_index_unpadded);
            tmp[7] = *potential_expansions_SoA.pointer<7>(cell_flat_index_unpadded);
            tmp[8] = *potential_expansions_SoA.pointer<8>(cell_flat_index_unpadded);
            tmp[9] = *potential_expansions_SoA.pointer<9>(cell_flat_index_unpadded);
            tmp[10] = *potential_expansions_SoA.pointer<10>(cell_flat_index_unpadded);
            tmp[11] = *potential_expansions_SoA.pointer<11>(cell_flat_index_unpadded);
            tmp[12] = *potential_expansions_SoA.pointer<12>(cell_flat_index_unpadded);
            tmp[13] = *potential_expansions_SoA.pointer<13>(cell_flat_index_unpadded);
            tmp[14] = *potential_expansions_SoA.pointer<14>(cell_flat_index_unpadded);
            tmp[15] = *potential_expansions_SoA.pointer<15>(cell_flat_index_unpadded);
            tmp[16] = *potential_expansions_SoA.pointer<16>(cell_flat_index_unpadded);
            tmp[17] = *potential_expansions_SoA.pointer<17>(cell_flat_index_unpadded);
            tmp[18] = *potential_expansions_SoA.pointer<18>(cell_flat_index_unpadded);
            tmp[19] = *potential_expansions_SoA.pointer<19>(cell_flat_index_unpadded);

            for (size_t i = 0; i < stencil_elements; i++) {
                octotiger::fmm::multiindex<>& cur_stencil = stencil[i];
                const octotiger::fmm::multiindex<> interaction_partner_index(
                    cell_index.x + cur_stencil.x, cell_index.y + cur_stencil.y,
                    cell_index.z + cur_stencil.z);

                const size_t interaction_partner_flat_index =
                    octotiger::fmm::to_flat_index_padded(interaction_partner_index);

                // implicitly broadcasts to vector
                octotiger::fmm::multiindex<> interaction_partner_index_coarse(
                    interaction_partner_index);
                // note that this is the same for groups of 2x2x2 elements
                // -> maps to the same for some SIMD lanes
                interaction_partner_index_coarse.transform_coarse();

                double theta_c_rec_squared =
                    static_cast<double>(detail::distance_squared_reciprocal(
                        cell_index_coarse, interaction_partner_index_coarse));

                const bool mask_b = theta_rec_squared > theta_c_rec_squared;
                const double mask = mask_b ? 1.0 : 0.0;
                // mask = !mask;

                Y[0] = center_of_masses_SoA.value<0>(interaction_partner_flat_index);
                Y[1] = center_of_masses_SoA.value<1>(interaction_partner_flat_index);
                Y[2] = center_of_masses_SoA.value<2>(interaction_partner_flat_index);

                dX[0] = X[0] - Y[0];
                dX[1] = X[1] - Y[1];
                dX[2] = X[2] - Y[2];

                m_partner[0] = local_expansions_SoA.value<0>(interaction_partner_flat_index);
                m_partner[1] = local_expansions_SoA.value<1>(interaction_partner_flat_index);
                m_partner[2] = local_expansions_SoA.value<2>(interaction_partner_flat_index);
                m_partner[3] = local_expansions_SoA.value<3>(interaction_partner_flat_index);
                m_partner[4] = local_expansions_SoA.value<4>(interaction_partner_flat_index);
                m_partner[5] = local_expansions_SoA.value<5>(interaction_partner_flat_index);
                m_partner[6] = local_expansions_SoA.value<6>(interaction_partner_flat_index);
                m_partner[7] = local_expansions_SoA.value<7>(interaction_partner_flat_index);
                m_partner[8] = local_expansions_SoA.value<8>(interaction_partner_flat_index);
                m_partner[9] = local_expansions_SoA.value<9>(interaction_partner_flat_index);
                m_partner[10] = local_expansions_SoA.value<10>(interaction_partner_flat_index);
                m_partner[11] = local_expansions_SoA.value<11>(interaction_partner_flat_index);
                m_partner[12] = local_expansions_SoA.value<12>(interaction_partner_flat_index);
                m_partner[13] = local_expansions_SoA.value<13>(interaction_partner_flat_index);
                m_partner[14] = local_expansions_SoA.value<14>(interaction_partner_flat_index);
                m_partner[15] = local_expansions_SoA.value<15>(interaction_partner_flat_index);
                m_partner[16] = local_expansions_SoA.value<16>(interaction_partner_flat_index);
                m_partner[17] = local_expansions_SoA.value<17>(interaction_partner_flat_index);
                m_partner[18] = local_expansions_SoA.value<18>(interaction_partner_flat_index);
                m_partner[19] = local_expansions_SoA.value<19>(interaction_partner_flat_index);

                // R_i in paper is the dX in the code
                // D is taylor expansion value for a given X expansion of the gravitational
                // potential
                // (multipole expansion)

                // calculates all D-values, calculate all coefficients of 1/r (not the
                // potential),
                // formula (6)-(9) and (19)
                D_split D_calculator(dX);
                // TODO: this translates to an initialization loop, bug!
                double D_lower[20];
                D_calculator.calculate_D_lower(D_lower);

                // expansion_v current_potential;
                // TODO: this translates to an initialization loop, bug!
                double cur_pot[10];

                // 10-19 are not cached!

                // the following loops calculate formula (10), potential from B->A

                cur_pot[0] = m_partner[0] * D_lower[0];
                cur_pot[1] = m_partner[0] * D_lower[1];
                cur_pot[2] = m_partner[0] * D_lower[2];
                cur_pot[3] = m_partner[0] * D_lower[3];

                cur_pot[0] -= m_partner[1] * D_lower[1];

                cur_pot[1] -= m_partner[1] * D_lower[4];
                cur_pot[1] -= m_partner[1] * D_lower[5];
                cur_pot[1] -= m_partner[1] * D_lower[6];

                cur_pot[0] -= m_partner[2] * D_lower[2];

                cur_pot[2] -= m_partner[2] * D_lower[5];
                cur_pot[2] -= m_partner[2] * D_lower[7];
                cur_pot[2] -= m_partner[2] * D_lower[8];

                cur_pot[0] -= m_partner[3] * D_lower[3];

                cur_pot[3] -= m_partner[3] * D_lower[6];
                cur_pot[3] -= m_partner[3] * D_lower[8];
                cur_pot[3] -= m_partner[3] * D_lower[9];

                cur_pot[0] += m_partner[4] * (D_lower[4] * factor_half[4]);
                cur_pot[1] += m_partner[4] * (D_lower[10] * factor_half[4]);
                cur_pot[2] += m_partner[4] * (D_lower[11] * factor_half[4]);
                cur_pot[3] += m_partner[4] * (D_lower[12] * factor_half[4]);

                cur_pot[0] += m_partner[5] * (D_lower[5] * factor_half[5]);
                cur_pot[1] += m_partner[5] * (D_lower[11] * factor_half[5]);
                cur_pot[2] += m_partner[5] * (D_lower[13] * factor_half[5]);
                cur_pot[3] += m_partner[5] * (D_lower[14] * factor_half[5]);

                cur_pot[0] += m_partner[6] * (D_lower[6] * factor_half[6]);
                cur_pot[1] += m_partner[6] * (D_lower[12] * factor_half[6]);
                cur_pot[2] += m_partner[6] * (D_lower[14] * factor_half[6]);
                cur_pot[3] += m_partner[6] * (D_lower[15] * factor_half[6]);

                cur_pot[0] += m_partner[7] * (D_lower[7] * factor_half[7]);
                cur_pot[1] += m_partner[7] * (D_lower[13] * factor_half[7]);
                cur_pot[2] += m_partner[7] * (D_lower[16] * factor_half[7]);
                cur_pot[3] += m_partner[7] * (D_lower[17] * factor_half[7]);

                cur_pot[0] += m_partner[8] * (D_lower[8] * factor_half[8]);
                cur_pot[1] += m_partner[8] * (D_lower[14] * factor_half[8]);
                cur_pot[2] += m_partner[8] * (D_lower[17] * factor_half[8]);
                cur_pot[3] += m_partner[8] * (D_lower[18] * factor_half[8]);

                cur_pot[0] += m_partner[9] * (D_lower[9] * factor_half[9]);
                cur_pot[1] += m_partner[9] * (D_lower[15] * factor_half[9]);
                cur_pot[2] += m_partner[9] * (D_lower[18] * factor_half[9]);
                cur_pot[3] += m_partner[9] * (D_lower[19] * factor_half[9]);

                cur_pot[0] -= m_partner[10] * (D_lower[10] * factor_sixth[10]);
                cur_pot[0] -= m_partner[11] * (D_lower[11] * factor_sixth[11]);
                cur_pot[0] -= m_partner[12] * (D_lower[12] * factor_sixth[12]);
                cur_pot[0] -= m_partner[13] * (D_lower[13] * factor_sixth[13]);
                cur_pot[0] -= m_partner[14] * (D_lower[14] * factor_sixth[14]);
                cur_pot[0] -= m_partner[15] * (D_lower[15] * factor_sixth[15]);
                cur_pot[0] -= m_partner[16] * (D_lower[16] * factor_sixth[16]);
                cur_pot[0] -= m_partner[17] * (D_lower[17] * factor_sixth[17]);
                cur_pot[0] -= m_partner[18] * (D_lower[18] * factor_sixth[18]);
                cur_pot[0] -= m_partner[19] * (D_lower[19] * factor_sixth[19]);

                cur_pot[4] = m_partner[0] * D_lower[4];
                cur_pot[5] = m_partner[0] * D_lower[5];
                cur_pot[6] = m_partner[0] * D_lower[6];
                cur_pot[7] = m_partner[0] * D_lower[7];
                cur_pot[8] = m_partner[0] * D_lower[8];
                cur_pot[9] = m_partner[0] * D_lower[9];

                cur_pot[4] -= m_partner[1] * D_lower[10];
                cur_pot[5] -= m_partner[1] * D_lower[11];
                cur_pot[6] -= m_partner[1] * D_lower[12];
                cur_pot[7] -= m_partner[1] * D_lower[13];
                cur_pot[8] -= m_partner[1] * D_lower[14];
                cur_pot[9] -= m_partner[1] * D_lower[15];

                cur_pot[4] -= m_partner[2] * D_lower[11];
                cur_pot[5] -= m_partner[2] * D_lower[13];
                cur_pot[6] -= m_partner[2] * D_lower[14];
                cur_pot[7] -= m_partner[2] * D_lower[16];
                cur_pot[8] -= m_partner[2] * D_lower[17];
                cur_pot[9] -= m_partner[2] * D_lower[18];

                cur_pot[4] -= m_partner[3] * D_lower[12];
                cur_pot[5] -= m_partner[3] * D_lower[14];
                cur_pot[6] -= m_partner[3] * D_lower[15];
                cur_pot[7] -= m_partner[3] * D_lower[17];
                cur_pot[8] -= m_partner[3] * D_lower[18];
                cur_pot[9] -= m_partner[3] * D_lower[19];

                tmp[0] +=
                    (potential_expansions_SoA.value<0>(cell_flat_index_unpadded) + cur_pot[0]) *
                    mask;
                tmp[1] +=
                    (potential_expansions_SoA.value<1>(cell_flat_index_unpadded) + cur_pot[1]) *
                    mask;
                tmp[2] +=
                    (potential_expansions_SoA.value<2>(cell_flat_index_unpadded) + cur_pot[2]) *
                    mask;
                tmp[3] +=
                    (potential_expansions_SoA.value<3>(cell_flat_index_unpadded) + cur_pot[3]) *
                    mask;
                tmp[4] +=
                    (potential_expansions_SoA.value<4>(cell_flat_index_unpadded) + cur_pot[4]) *
                    mask;
                tmp[5] +=
                    (potential_expansions_SoA.value<5>(cell_flat_index_unpadded) + cur_pot[5]) *
                    mask;
                tmp[6] +=
                    (potential_expansions_SoA.value<6>(cell_flat_index_unpadded) + cur_pot[6]) *
                    mask;
                tmp[7] +=
                    (potential_expansions_SoA.value<7>(cell_flat_index_unpadded) + cur_pot[7]) *
                    mask;
                tmp[8] +=
                    (potential_expansions_SoA.value<8>(cell_flat_index_unpadded) + cur_pot[8]) *
                    mask;
                tmp[9] +=
                    (potential_expansions_SoA.value<9>(cell_flat_index_unpadded) + cur_pot[9]) *
                    mask;
                tmp[10] += (potential_expansions_SoA.value<10>(cell_flat_index_unpadded) +
                               m_partner[0] * D_lower[10]) *
                    mask;
                tmp[11] += (potential_expansions_SoA.value<11>(cell_flat_index_unpadded) +
                               m_partner[0] * D_lower[11]) *
                    mask;
                tmp[12] += (potential_expansions_SoA.value<12>(cell_flat_index_unpadded) +
                               m_partner[0] * D_lower[12]) *
                    mask;
                tmp[13] += (potential_expansions_SoA.value<13>(cell_flat_index_unpadded) +
                               m_partner[0] * D_lower[13]) *
                    mask;
                tmp[14] += (potential_expansions_SoA.value<14>(cell_flat_index_unpadded) +
                               m_partner[0] * D_lower[14]) *
                    mask;
                tmp[15] += (potential_expansions_SoA.value<15>(cell_flat_index_unpadded) +
                               m_partner[0] * D_lower[15]) *
                    mask;
                tmp[16] += (potential_expansions_SoA.value<16>(cell_flat_index_unpadded) +
                               m_partner[0] * D_lower[16]) *
                    mask;
                tmp[17] += (potential_expansions_SoA.value<17>(cell_flat_index_unpadded) +
                               m_partner[0] * D_lower[17]) *
                    mask;
                tmp[18] += (potential_expansions_SoA.value<18>(cell_flat_index_unpadded) +
                               m_partner[0] * D_lower[18]) *
                    mask;
                tmp[19] += (potential_expansions_SoA.value<19>(cell_flat_index_unpadded) +
                               m_partner[0] * D_lower[19]) *
                    mask;
            }
            *potential_expansions_SoA.pointer<0>(cell_flat_index_unpadded) = tmp[0];
            *potential_expansions_SoA.pointer<1>(cell_flat_index_unpadded) = tmp[1];
            *potential_expansions_SoA.pointer<2>(cell_flat_index_unpadded) = tmp[2];
            *potential_expansions_SoA.pointer<3>(cell_flat_index_unpadded) = tmp[3];
            *potential_expansions_SoA.pointer<4>(cell_flat_index_unpadded) = tmp[4];
            *potential_expansions_SoA.pointer<5>(cell_flat_index_unpadded) = tmp[5];
            *potential_expansions_SoA.pointer<6>(cell_flat_index_unpadded) = tmp[6];
            *potential_expansions_SoA.pointer<7>(cell_flat_index_unpadded) = tmp[7];
            *potential_expansions_SoA.pointer<8>(cell_flat_index_unpadded) = tmp[8];
            *potential_expansions_SoA.pointer<9>(cell_flat_index_unpadded) = tmp[9];
            *potential_expansions_SoA.pointer<10>(cell_flat_index_unpadded) = tmp[10];
            *potential_expansions_SoA.pointer<11>(cell_flat_index_unpadded) = tmp[11];
            *potential_expansions_SoA.pointer<12>(cell_flat_index_unpadded) = tmp[12];
            *potential_expansions_SoA.pointer<13>(cell_flat_index_unpadded) = tmp[13];
            *potential_expansions_SoA.pointer<14>(cell_flat_index_unpadded) = tmp[14];
            *potential_expansions_SoA.pointer<15>(cell_flat_index_unpadded) = tmp[15];
            *potential_expansions_SoA.pointer<16>(cell_flat_index_unpadded) = tmp[16];
            *potential_expansions_SoA.pointer<17>(cell_flat_index_unpadded) = tmp[17];
            *potential_expansions_SoA.pointer<18>(cell_flat_index_unpadded) = tmp[18];
            *potential_expansions_SoA.pointer<19>(cell_flat_index_unpadded) = tmp[19];
        }

        void m2m_cuda::compute_interactions(octotiger::fmm::struct_of_array_data<double, 20,
                                                ENTRIES, SOA_PADDING>& local_expansions_SoA,

            octotiger::fmm::struct_of_array_data<double, 3, ENTRIES, SOA_PADDING>&
                center_of_masses_SoA,
            octotiger::fmm::struct_of_array_data<double, 20, INNER_CELLS, SOA_PADDING>&
                potential_expansions_SoA,
            octotiger::fmm::struct_of_array_data<double, 3, INNER_CELLS, SOA_PADDING>&
                angular_corrections_SoA,
            double theta, std::array<double, 20>& factor, std::array<double, 20>& factor_half,
            std::array<double, 20>& factor_sixth, gsolve_type type) {
            int nDevices;

            // cudaGetDeviceCount(&nDevices);
            // for (int i = 0; i < nDevices; i++) {
            //   cudaDeviceProp prop;
            //   cudaGetDeviceProperties(&prop, i);
            //   printf("Device Number: %d\n", i);
            //   printf("  Device name: %s\n", prop.name);
            //   printf("  Memory Clock Rate (KHz): %d\n",
            //          prop.memoryClockRate);
            //   printf("  Memory Bus Width (bits): %d\n",
            //          prop.memoryBusWidth);
            //   printf("  Peak Memory Bandwidth (GB/s): %f\n\n",
            //          2.0*prop.memoryClockRate*(prop.memoryBusWidth/8)/1.0e6);
            // }
            // std::cin.get();
            CUDA_CHECK_ERROR(cudaSetDevice(0));
            std::cout << "moving local_expansions_SoA" << std::endl;
            cuda_buffer<double> d_local_expansions_SoA =
                octotiger::fmm::cuda::move_to_device(local_expansions_SoA);
            std::cout << "moving center_of_masses_SoA" << std::endl;
            cuda_buffer<double> d_center_of_masses_SoA =
                octotiger::fmm::cuda::move_to_device(center_of_masses_SoA);
            std::cout << "moving potential_expansions_SoA" << std::endl;
            cuda_buffer<double> d_potential_expansions_SoA =
                octotiger::fmm::cuda::move_to_device(potential_expansions_SoA);
            std::cout << "moving angular_corrections_SoA" << std::endl;
            cuda_buffer<double> d_angular_corrections_SoA =
                octotiger::fmm::cuda::move_to_device(angular_corrections_SoA);
            // std::cin.get();

            octotiger::fmm::struct_of_array_data<double, 20, INNER_CELLS, SOA_PADDING>
                potential_expansions_SoA_copy;
            d_potential_expansions_SoA.move_to_host(potential_expansions_SoA_copy);

            octotiger::fmm::struct_of_array_data<double, 3, INNER_CELLS, SOA_PADDING>
                angular_corrections_SoA_copy;
            d_angular_corrections_SoA.move_to_host(angular_corrections_SoA_copy);

            std::cout << "moving factor" << std::endl;
            cuda_buffer<double> d_factor = octotiger::fmm::cuda::move_to_device(factor);
            std::cout << "moving factor_half" << std::endl;
            cuda_buffer<double> d_factor_half = octotiger::fmm::cuda::move_to_device(factor_half);
            std::cout << "moving factor_sixth" << std::endl;
            cuda_buffer<double> d_factor_sixth = octotiger::fmm::cuda::move_to_device(factor_sixth);
            std::cout << "moving stencil" << std::endl;
            cuda_buffer<multiindex<>> d_stencil = octotiger::fmm::cuda::move_to_device(stencil);
            // for (multiindex<> s_e : stencil) {
            //     std::cout << "x: " << s_e.x << " y: " << s_e.y << " z: " << s_e.z << std::endl;
            // }

            // int num_blocks = 1;
            dim3 grid_spec(1, 1, 1);
            dim3 threads_per_block(8, 8, 8);
            if (type == RHO) {
                blocked_interaction_rho<<<grid_spec, threads_per_block>>>(
                    d_local_expansions_SoA.get_device_ptr(),
                    d_center_of_masses_SoA.get_device_ptr(),
                    d_potential_expansions_SoA.get_device_ptr(),
                    d_angular_corrections_SoA.get_device_ptr(), d_stencil.get_device_ptr(),
                    stencil.size(), theta, d_factor.get_device_ptr(),
                    d_factor_half.get_device_ptr(), d_factor_sixth.get_device_ptr());
                CUDA_CHECK_ERROR(cudaThreadSynchronize());
            } else {
                blocked_interaction_non_rho<<<grid_spec, threads_per_block>>>(
                    d_local_expansions_SoA.get_device_ptr(),
                    d_center_of_masses_SoA.get_device_ptr(),
                    d_potential_expansions_SoA.get_device_ptr(),
                    d_angular_corrections_SoA.get_device_ptr(), d_stencil.get_device_ptr(),
                    stencil.size(), theta, d_factor.get_device_ptr(),
                    d_factor_half.get_device_ptr(), d_factor_sixth.get_device_ptr());
                CUDA_CHECK_ERROR(cudaThreadSynchronize());
            }

            d_potential_expansions_SoA.move_to_host(potential_expansions_SoA);
            d_angular_corrections_SoA.move_to_host(angular_corrections_SoA);
        }
    }
}
}
