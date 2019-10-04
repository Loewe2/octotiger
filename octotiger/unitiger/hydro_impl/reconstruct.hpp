//  Copyright (c) 2019 AUTHORS
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#pragma once
#include "octotiger/unitiger/physics.hpp"
#include "octotiger/unitiger/physics_impl.hpp"

#include <octotiger/cuda_util/cuda_helper.hpp>
#include <octotiger/cuda_util/cuda_scheduler.hpp>
#include <octotiger/unitiger/hydro_impl/hydro_cuda/hydro_cuda_kernels.hpp>
#include <octotiger/common_kernel/struct_of_array_data.hpp>

//#ifdef OCTOTIGER_WITH_CUDA
template<int NDIM, int INX>
const hydro::recon_type<NDIM>& hydro_computer<NDIM, INX>::reconstruct_cuda(hydro::state_type &U_, const hydro::x_type &X, safe_real omega) {

//	static thread_local octotiger::fmm::struct_of_array_data<std::array<safe_real, geo::NDIR>, safe_real, geo::NDIR, geo::H_N3, 19>
//		D1_SoA;
	static thread_local auto D1 = std::vector<std::array<safe_real, geo::NDIR>>(geo::H_N3);
	static thread_local auto Q = std::vector < std::vector<std::array<safe_real, geo::NDIR>> > (nf_, std::vector<std::array<safe_real, geo::NDIR>>(geo::H_N3));

#ifdef OCTOTIGER_HAVE_CUDA
	if constexpr (geo::NDIR == 27) {
		static thread_local octotiger::fmm::struct_of_array_data<std::array<safe_real, geo::NDIR>, safe_real, geo::NDIR, geo::H_N3, 0, octotiger::fmm::pinned_vector<safe_real>>
		D1_SoA;
		static thread_local std::vector<octotiger::fmm::struct_of_array_data<std::array<safe_real, geo::NDIR>, safe_real, geo::NDIR, geo::H_N3, 0, octotiger::fmm::pinned_vector<safe_real>>> Q_SoA(nf_);

		octotiger::fmm::struct_of_array_data<std::array<safe_real, geo::NDIR>, safe_real, geo::NDIR, geo::H_N3, 0, octotiger::fmm::pinned_vector<safe_real>>
			U_SoA;
			U_SoA.concatenate_vectors(U_);
		octotiger::fmm::struct_of_array_data<std::array<safe_real, NDIM>, safe_real, NDIM, geo::H_N3, 0, octotiger::fmm::pinned_vector<safe_real>>
			X_SoA;
			X_SoA.concatenate_vectors(X);

		//reconstruct_kernel_interface(D1_SoA, Q_SoA, U_SoA, X_SoA);
	} else {
		std::cerr << "CUDA is currently only supported for 3D problems" << std::endl;

	}
#endif
	return Q;
}

const std::vector<bool> alltrue = {true, true, true, true, true, true, true, true, true, true, true, true, true, true, true};
const std::vector<bool> allfalse = {false, false, false, false, false, false, false, false, false, false, false, false, false, false, false};

template<int NDIM, int INX>
void hydro_computer<NDIM, INX>::reconstruct_ppm(std::vector<std::vector<std::vector<safe_real>>> &Q_SoA,
 hydro::state_type &U_SoA, const hydro::x_type &X, safe_real omega, int face_offset, int faces, const std::vector<bool> &smooth) {
	if constexpr (geo::NDIR != 27) {
		reconstruct_ppm_cpu(Q_SoA, U_SoA, omega, face_offset, faces, smooth);
	} else {
		octotiger::fmm::kernel_scheduler::scheduler().init();
		// Get Slot
		int slot = octotiger::fmm::kernel_scheduler::scheduler().get_launch_slot();
		if (slot == -1) {
			reconstruct_ppm_cpu(Q_SoA, U_SoA, omega, face_offset, faces, smooth);
		} else {
			// Get interface
			octotiger::util::cuda_helper& gpu_interface =
				octotiger::fmm::kernel_scheduler::scheduler().get_launch_interface(slot);

			// Get staging area
			auto staging_area =
				octotiger::fmm::kernel_scheduler::scheduler().get_hydro_staging_area(slot);
			auto env =
				octotiger::fmm::kernel_scheduler::scheduler().get_hydro_device_enviroment(slot);

			std::vector<octotiger::fmm::struct_of_array_data<std::array<safe_real, 27>, safe_real, 27, 2744, 0, octotiger::fmm::pinned_vector<safe_real>>> &D1_SoA =
				staging_area.D1_SoA;

			octotiger::fmm::struct_of_array_data<std::array<safe_real, 15>, safe_real, 15, geo::H_N3, 0, octotiger::fmm::pinned_vector<safe_real>>
				&U_SoA2 = staging_area.U_SoA;
			U_SoA2.concatenate_vectors(U_SoA);

			std::vector<octotiger::fmm::struct_of_array_data<std::vector<safe_real>, safe_real, geo::NDIR, geo::H_N3, 0, octotiger::fmm::pinned_vector<safe_real>>>
				&Q_SoA2 = staging_area.Q1_SoA;
			// for (auto f = face_offset; f < faces; f++) {
			// Q_SoA2[f].concatenate_vectors(Q_SoA[f]);
			// }

			reconstruct_ppm_interface(D1_SoA, Q_SoA2, U_SoA2, slot, face_offset, faces);

			for (auto f = face_offset; f < faces; f++) {
				for (auto component = 0; component < geo::NDIR; component++) {
					Q_SoA2[f].copy_component(component, Q_SoA[f][component]);
				}
			}

			reconstruct_ppm_partial(Q_SoA, U_SoA, omega, face_offset, faces, smooth);
		}


	}

}
template<int NDIM, int INX>
void hydro_computer<NDIM, INX>::reconstruct_ppm_partial(std::vector<std::vector<std::vector<safe_real>>> &Q_SoA,
 hydro::state_type &U_SoA, safe_real omega, int face_offset, int faces, const std::vector<bool> &smooth) {

	static constexpr auto dir = geo::direction();
	// static thread_local auto D1_SoA = std::vector < std::vector < safe_real >> (geo::NDIR, std::vector < safe_real > (geo::H_N3));

 	for (int f = face_offset; f < faces; f++) {
		/*for (int j = 0; j < geo::H_NX_XM2; j++) {
			for (int k = 0; k < geo::H_NX_YM2; k++) {
				for (int l = 0; l < geo::H_NX_ZM2; l++) {
					const int i = geo::to_index(j + 1, k + 1, l + 1);
					for (int d = 0; d < geo::NDIR; d++) {
						const auto di = dir[d];
						D1_SoA[d][i] = minmod_theta(U_SoA[f][i + di] - U_SoA[f][i], U_SoA[f][i] - U_SoA[f][i - di], 2.0);
					}
				}
			}
		}

		for (int d = 0; d < geo::NDIR; d++) {
			const auto di = dir[d];
			for (int j = 0; j < geo::H_NX_XM2; j++) {
			for (int k = 0; k < geo::H_NX_YM2; k++) {
				for (int l = 0; l < geo::H_NX_ZM2; l++) {
						const int i = geo::to_index(j + 1, k + 1, l + 1);
						Q_SoA[f][d][i] = 0.5 * (U_SoA[f][i] + U_SoA[f][i + di]);
						Q_SoA[f][d][i] += (1.0 / 6.0) * (D1_SoA[d][i] - D1_SoA[d][i + di]);
					}
				}
			}
		} */
#ifndef DISABLE_VERTEX_AVG
		for (int j = 0; j < geo::H_NX_XM2; j++) {
			for (int k = 0; k < geo::H_NX_YM2; k++) {
				for (int l = 0; l < geo::H_NX_ZM2; l++) {
					const int i = geo::to_index(j + 1, k + 1, l + 1);
					for (int gi = 0; gi < geo::group_count(); gi++) {
						safe_real sum = 0.0;
						for (int n = 0; n < geo::group_size(gi); n++) {
							const auto pair = geo::group_pair(gi, n);
							sum += Q_SoA[f][pair.second][i + pair.first];
						}
						sum /= safe_real(geo::group_size(gi));
						for (int n = 0; n < geo::group_size(gi); n++) {
							const auto pair = geo::group_pair(gi, n);
							Q_SoA[f][pair.second][i + pair.first] = sum;
						}
					}
				}
			}
		}
		for (int d = 0; d < geo::NDIR; d++) {
			if (d != geo::NDIR / 2) {
				const auto di = dir[d];
				for (int j = 0; j < geo::H_NX_XM2; j++) {
					for (int k = 0; k < geo::H_NX_YM2; k++) {
						for (int l = 0; l < geo::H_NX_ZM2; l++) {
							const int i = geo::to_index(j + 1, k + 1, l + 1);
							const auto M = std::max(U_SoA[f][i], U_SoA[f][i + di]);
							const auto m = std::min(U_SoA[f][i], U_SoA[f][i + di]);
							Q_SoA[f][d][i] = std::max(Q_SoA[f][d][i], m);
							Q_SoA[f][d][i] = std::min(Q_SoA[f][d][i], M);
						}
					}
				}
			}
		}
#endif

		if (!smooth[f]) {
			for (int d = 0; d < geo::NDIR / 2; d++) {
				for (int j = 0; j < geo::H_NX_XM4; j++) {
					for (int k = 0; k < geo::H_NX_YM4; k++) {
						for (int l = 0; l < geo::H_NX_ZM4; l++) {
							const int i = geo::to_index(j + 2, k + 2, l + 2);
							auto &qp = Q_SoA[f][geo::flip(d)][i];
							auto &qm = Q_SoA[f][d][i];
							limit_slope(qm, U_SoA[f][i], qp);
						}
					}
				}
			}
		}
	}
}


template<int NDIM, int INX>
void hydro_computer<NDIM, INX>::reconstruct_ppm_cpu(std::vector<std::vector<std::vector<safe_real>>> &Q_SoA,
 hydro::state_type &U_SoA, safe_real omega, int face_offset, int faces, const std::vector<bool> &smooth) {

	static thread_local auto D1_SoA = std::vector < std::vector < safe_real >> (geo::NDIR, std::vector < safe_real > (geo::H_N3));
	static constexpr auto dir = geo::direction();

	for (int f = face_offset; f < faces; f++) {
		for (int j = 0; j < geo::H_NX_XM2; j++) {
			for (int k = 0; k < geo::H_NX_YM2; k++) {
				for (int l = 0; l < geo::H_NX_ZM2; l++) {
					const int i = geo::to_index(j + 1, k + 1, l + 1);
					for (int d = 0; d < geo::NDIR; d++) {
						const auto di = dir[d];
						D1_SoA[d][i] = minmod_theta(U_SoA[f][i + di] - U_SoA[f][i], U_SoA[f][i] - U_SoA[f][i - di], 2.0);
					}
				}
			}
		}

		for (int d = 0; d < geo::NDIR; d++) {
			const auto di = dir[d];
			for (int j = 0; j < geo::H_NX_XM2; j++) {
			for (int k = 0; k < geo::H_NX_YM2; k++) {
				for (int l = 0; l < geo::H_NX_ZM2; l++) {
						const int i = geo::to_index(j + 1, k + 1, l + 1);
						Q_SoA[f][d][i] = 0.5 * (U_SoA[f][i] + U_SoA[f][i + di]);
						Q_SoA[f][d][i] += (1.0 / 6.0) * (D1_SoA[d][i] - D1_SoA[d][i + di]);
					}
				}
			}
		}

#ifndef DISABLE_VERTEX_AVG
		for (int j = 0; j < geo::H_NX_XM2; j++) {
			for (int k = 0; k < geo::H_NX_YM2; k++) {
				for (int l = 0; l < geo::H_NX_ZM2; l++) {
					const int i = geo::to_index(j + 1, k + 1, l + 1);
					for (int gi = 0; gi < geo::group_count(); gi++) {
						safe_real sum = 0.0;
						for (int n = 0; n < geo::group_size(gi); n++) {
							const auto pair = geo::group_pair(gi, n);
							sum += Q_SoA[f][pair.second][i + pair.first];
						}
						sum /= safe_real(geo::group_size(gi));
						for (int n = 0; n < geo::group_size(gi); n++) {
							const auto pair = geo::group_pair(gi, n);
							Q_SoA[f][pair.second][i + pair.first] = sum;
						}
					}
				}
			}
		}
		for (int d = 0; d < geo::NDIR; d++) {
			if (d != geo::NDIR / 2) {
				const auto di = dir[d];
				for (int j = 0; j < geo::H_NX_XM2; j++) {
					for (int k = 0; k < geo::H_NX_YM2; k++) {
						for (int l = 0; l < geo::H_NX_ZM2; l++) {
							const int i = geo::to_index(j + 1, k + 1, l + 1);
							const auto M = std::max(U_SoA[f][i], U_SoA[f][i + di]);
							const auto m = std::min(U_SoA[f][i], U_SoA[f][i + di]);
							Q_SoA[f][d][i] = std::max(Q_SoA[f][d][i], m);
							Q_SoA[f][d][i] = std::min(Q_SoA[f][d][i], M);
						}
					}
				}
			}
		}
#endif

		if (!smooth[f]) {
			for (int d = 0; d < geo::NDIR / 2; d++) {
				for (int j = 0; j < geo::H_NX_XM4; j++) {
					for (int k = 0; k < geo::H_NX_YM4; k++) {
						for (int l = 0; l < geo::H_NX_ZM4; l++) {
							const int i = geo::to_index(j + 2, k + 2, l + 2);
							auto &qp = Q_SoA[f][geo::flip(d)][i];
							auto &qm = Q_SoA[f][d][i];
							limit_slope(qm, U_SoA[f][i], qp);
						}
					}
				}
			}
		}
	}
}
//#endif

template<int NDIM, int INX>
const hydro::recon_type<NDIM>& hydro_computer<NDIM, INX>::reconstruct(hydro::state_type &U_, const hydro::x_type &X, safe_real omega) {

	static thread_local auto Q = std::vector < std::vector<std::array<safe_real, geo::NDIR>> > (nf_, std::vector<std::array<safe_real, geo::NDIR>>(geo::H_N3));

	static thread_local auto Q_SoA = std::vector < std::vector<std::vector<safe_real>>
			> (nf_, std::vector < std::vector < safe_real >> (geo::NDIR, std::vector < safe_real > (geo::H_N3)));
	static thread_local auto D1_SoA = std::vector < std::vector < safe_real >> (geo::NDIR, std::vector < safe_real > (geo::H_N3));

	static const auto SoA2AoS = [](int f1, int f2) {
		for (int f = f1; f < f2; f++) {
			for (int i = 0; i < geo::H_N3; i++) {
				for (int d = 0; d < geo::NDIR; d++) {
					Q[f][i][d] = Q_SoA[f][d][i];
				}
			}
		}
	};

	static const auto AoS2SoA = [](int f1, int f2) {
		for (int f = f1; f < f2; f++) {
			for (int i = 0; i < geo::H_N3; i++) {
				for (int d = 0; d < geo::NDIR; d++) {
					Q_SoA[f][d][i] = Q[f][i][d];
				}
			}
		}
	};
	reconstruct_cuda(U_, X, omega);
	static constexpr auto xloc = geo::xloc();
	static constexpr auto kdelta = geo::kronecker_delta();
	static constexpr auto vw = geo::volume_weight();
	static constexpr auto dir = geo::direction();

	const auto dx = X[0][geo::H_DNX] - X[0][0];
	auto U = physics < NDIM > ::template pre_recon<INX>(U_, X, omega, angmom_count_ > 0);

	const auto measure_angmom = [dx](const std::array<std::array<safe_real, geo::NDIR>, NDIM> &C) {
		std::array < safe_real, geo::NANGMOM > L;
		for (int n = 0; n < geo::NANGMOM; n++) {
			L[n] = 0.0;
			for (int m = 0; m < NDIM; m++) {
				for (int l = 0; l < NDIM; l++) {
					for (int d = 0; d < geo::NDIR; d++) {
						if (d != geo::NDIR / 2) {
							L[n] += vw[d] * kdelta[n][m][l] * 0.5 * xloc[d][m] * C[l][d] * dx;
						}
					}
				}
			}
		}
		return L;
	};

	const auto add_angmom = [dx](std::array<std::array<safe_real, geo::NDIR>, NDIM> &C, std::array<safe_real, geo::NANGMOM> &Z) {
		for (int d = 0; d < geo::NDIR; d++) {
			if (d != geo::NDIR / 2) {
				for (int n = 0; n < geo::NANGMOM; n++) {
					for (int m = 0; m < NDIM; m++) {
						for (int l = 0; l < NDIM; l++) {
							const auto tmp = 6.0 * Z[n] / dx;
							C[l][d] += kdelta[n][m][l] * 0.5 * xloc[d][m] * tmp;
						}
					}
				}
			}
		}
	};

	if (angmom_count_ == 0 || NDIM == 1) {
		reconstruct_ppm(Q_SoA, U, X, omega, 0, nf_, smooth_field_);

	} else {
		reconstruct_ppm(Q_SoA, U, X, omega, 0, angmom_index_, smooth_field_);

		int sx_i = angmom_index_;
		int zx_i = sx_i + NDIM;

		for (int angmom_pair = 0; angmom_pair < angmom_count_; angmom_pair++) {
			reconstruct_ppm(Q_SoA, U, X, omega, sx_i, sx_i + NDIM, alltrue);

			SoA2AoS(rho_i, rho_i + 1);
			SoA2AoS(sx_i, sx_i + NDIM);

			for (int j = 0; j < geo::H_NX_XM4; j++) {
				for (int k = 0; k < geo::H_NX_YM4; k++) {
					for (int l = 0; l < geo::H_NX_ZM4; l++) {
						const int i = geo::to_index(j + 2, k + 2, l + 2);

						std::array < safe_real, geo::NANGMOM > Z;
						std::array<std::array<safe_real, geo::NDIR>, NDIM> S;
						for (int dim = 0; dim < geo::NANGMOM; dim++) {
							Z[dim] = U[zx_i + dim][i];
						}
						for (int dim = 0; dim < NDIM; dim++) {
							for (int d = 0; d < geo::NDIR; d++) {
								S[dim][d] = Q[sx_i + dim][i][d];
							}
						}

						physics < NDIM > ::template pre_angmom<INX>(U, Q, Z, S, i, dx);
						auto am1 = measure_angmom(S);
						decltype(Z) am2;
						for (int dim = 0; dim < geo::NANGMOM; dim++) {
							am2[dim] = Z[dim] - am1[dim];
						}
						add_angmom(S, am2);
						physics < NDIM > ::template post_angmom<INX>(U, Q, Z, S, i, dx);

						for (int dim = 0; dim < NDIM; dim++) {
							for (int d = 0; d < geo::NDIR; d++) {
								if (d != geo::NDIR / 2) {
									auto &s = S[dim][d];
									const auto &q = U[sx_i + dim][i + dir[d]];
									const auto &u0 = U[sx_i + dim][i];
									const auto M = std::max(u0, q);
									const auto m = std::min(u0, q);
									s = std::min(s, M);
									s = std::max(s, m);
								}
							}
						}
						for (int f = sx_i; f < sx_i + NDIM; f++) {
							const auto dim = f - sx_i;
							for (int d = 0; d < geo::NDIR / 2; d++) {
								limit_slope(S[dim][d], U[f][i], S[dim][geo::flip(d)]);
							}
						}

						for (int dim = 0; dim < NDIM; dim++) {
							for (int d = 0; d < geo::NDIR; d++) {
								Q[sx_i + dim][i][d] = S[dim][d];
							}
						}
					}
				}
			}

			AoS2SoA(sx_i, zx_i + geo::NANGMOM);

			reconstruct_ppm(Q_SoA, U, X, omega, zx_i, zx_i + geo::NANGMOM, allfalse);
			sx_i += geo::NANGMOM + NDIM;
			zx_i += geo::NANGMOM + NDIM;
		}
		reconstruct_ppm(Q_SoA, U, X, omega, angmom_index_ + angmom_count_ * (geo::NANGMOM + NDIM), nf_, smooth_field_);
	}

	Q_SoA = physics < NDIM > ::template post_recon<INX>(Q_SoA, X, omega, angmom_count_ > 0);

	SoA2AoS(0, nf_);
	return Q;
}

