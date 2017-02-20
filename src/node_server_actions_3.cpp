#include "node_server.hpp"
#include "node_client.hpp"
#include "future.hpp"
#include "options.hpp"
#include "util.hpp"

#include <hpx/include/run_as.hpp>
#include <hpx/include/lcos.hpp>
#include <hpx/include/util.hpp>

extern options opts;

typedef node_server::send_gravity_boundary_action send_gravity_boundary_action_type;
HPX_REGISTER_ACTION(send_gravity_boundary_action_type);

void node_client::send_gravity_boundary(gravity_boundary_type&& data,
    const geo::direction& dir, bool monopole) const {
    hpx::apply<typename node_server::send_gravity_boundary_action>(get_unmanaged_gid(),
        std::move(data), dir, monopole);
}

void node_server::recv_gravity_boundary(gravity_boundary_type&& bdata,
    const geo::direction& dir, bool monopole) {
    neighbor_gravity_type tmp;
    tmp.data = std::move(bdata);
    tmp.is_monopole = monopole;
    tmp.direction = dir;
    neighbor_gravity_channels[dir].set_value(std::move(tmp));
}

typedef node_server::send_gravity_expansions_action send_gravity_expansions_action_type;
HPX_REGISTER_ACTION(send_gravity_expansions_action_type);

void node_server::recv_gravity_expansions(expansion_pass_type&& v) {
    parent_gravity_channel.set_value(std::move(v));
}

void node_client::send_gravity_expansions(expansion_pass_type&& data) const {
    hpx::apply<typename node_server::send_gravity_expansions_action>(get_unmanaged_gid(),
        std::move(data));
}

typedef node_server::send_gravity_multipoles_action send_gravity_multipoles_action_type;
HPX_REGISTER_ACTION(send_gravity_multipoles_action_type);

void node_client::send_gravity_multipoles(multipole_pass_type&& data,
    const geo::octant& ci) const {
    hpx::apply<typename node_server::send_gravity_multipoles_action>(get_unmanaged_gid(),
        std::move(data), ci);
}

void node_server::recv_gravity_multipoles(multipole_pass_type&& v,
    const geo::octant& ci) {
    child_gravity_channels[ci].set_value(std::move(v));
}

typedef node_server::send_hydro_boundary_action send_hydro_boundary_action_type;
HPX_REGISTER_ACTION(send_hydro_boundary_action_type);

void node_client::send_hydro_boundary(std::vector<real>&& data,
    const geo::direction& dir) const {
    hpx::apply<typename node_server::send_hydro_boundary_action>(get_unmanaged_gid(),
        std::move(data), dir);
}

void node_server::recv_hydro_boundary(std::vector<real>&& bdata,
    const geo::direction& dir) {
    sibling_hydro_type tmp;
    tmp.data = std::move(bdata);
    tmp.direction = dir;
    sibling_hydro_channels[dir].set_value(std::move(tmp));
}

typedef node_server::send_hydro_children_action send_hydro_children_action_type;
HPX_REGISTER_ACTION(send_hydro_children_action_type);

void node_server::recv_hydro_children(std::vector<real>&& data, const geo::octant& ci) {
    child_hydro_channels[ci].set_value(std::move(data));
}

void node_client::send_hydro_children(std::vector<real>&& data,
    const geo::octant& ci) const {
    hpx::apply<typename node_server::send_hydro_children_action>(get_unmanaged_gid(),
        std::move(data), ci);
}

typedef node_server::send_hydro_flux_correct_action send_hydro_flux_correct_action_type;
HPX_REGISTER_ACTION(send_hydro_flux_correct_action_type);

void node_client::send_hydro_flux_correct(std::vector<real>&& data,
    const geo::face& face,
    const geo::octant& ci) const {
    hpx::apply<typename node_server::send_hydro_flux_correct_action>(get_unmanaged_gid(),
        std::move(data), face, ci);
}

void node_server::recv_hydro_flux_correct(std::vector<real>&& data, const geo::face& face,
    const geo::octant& ci) {
    const geo::quadrant index(ci, face.get_dimension());
    niece_hydro_channels[face][index].set_value(std::move(data));
}

typedef node_server::line_of_centers_action line_of_centers_action_type;
HPX_REGISTER_ACTION(line_of_centers_action_type);

hpx::future<line_of_centers_t> node_client::line_of_centers(
    const std::pair<space_vector, space_vector>& line) const {
    return hpx::async<typename node_server::line_of_centers_action>(get_unmanaged_gid(), line);
}

void output_line_of_centers(FILE* fp, const line_of_centers_t& loc) {
    for (integer i = 0; i != loc.size(); ++i) {
        fprintf(fp, "%e ", loc[i].first);
        for (integer j = 0; j != NF + NGF; ++j) {
            fprintf(fp, "%e ", loc[i].second[j]);
        }
        fprintf(fp, "\n");
    }
}

line_of_centers_t node_server::line_of_centers(
    const std::pair<space_vector, space_vector>& line) const {
    line_of_centers_t return_line;
    if (is_refined) {
        std::array<hpx::future<line_of_centers_t>, NCHILD> futs;
        for (integer ci = 0; ci != NCHILD; ++ci) {
            futs[ci] = children[ci].line_of_centers(line);
        }
        std::map<real, std::vector<real>> map;
        for (auto&& fut : futs) {
            auto tmp = fut.get();
            for (integer ii = 0; ii != tmp.size(); ++ii) {
                if (map.find(tmp[ii].first) == map.end()) {
                    map.emplace(std::move(tmp[ii]));
                }
            }
        }
        return_line.resize(map.size());
        std::move(map.begin(), map.end(), return_line.begin());
    } else {
        return_line = grid_ptr->line_of_centers(line);
    }

    return return_line;
}

void line_of_centers_analyze(const line_of_centers_t& loc, real omega,
    std::pair<real, real>& rho1_max, std::pair<real, real>& rho2_max,
    std::pair<real, real>& l1_phi, std::pair<real, real>& l2_phi,
    std::pair<real, real>& l3_phi, real& rho1_phi, real& rho2_phi) {

    for (auto& l : loc) {
        ASSERT_NONAN(l.first);
        for (integer f = 0; f != NF + NGF; ++f) {
            ASSERT_NONAN(l.second[f]);
        }
    }

    rho1_max.second = rho2_max.second = 0.0;
    integer rho1_maxi, rho2_maxi;
    ///	printf( "LOCSIZE %i\n", loc.size());
    for (integer i = 0; i != loc.size(); ++i) {
        const real x = loc[i].first;
        const real rho = loc[i].second[rho_i];
        const real pot = loc[i].second[pot_i];
        //	printf( "%e %e\n", x, rho);
        if (rho1_max.second < rho) {
            //	printf( "!\n");
            rho1_max.second = rho;
            rho1_max.first = x;
            rho1_maxi = i;
            real phi_eff = pot / ASSERT_POSITIVE(rho) - 0.5 * x * x * omega * omega;
            rho1_phi = phi_eff;
        }
    }
    for (integer i = 0; i != loc.size(); ++i) {
        const real x = loc[i].first;
        if (x * rho1_max.first < 0.0) {
            const real rho = loc[i].second[rho_i];
            const real pot = loc[i].second[pot_i];
            if (rho2_max.second < rho) {
                rho2_max.second = rho;
                rho2_max.first = x;
                rho2_maxi = i;
                real phi_eff = pot / ASSERT_POSITIVE(rho) - 0.5 * x * x * omega * omega;
                rho2_phi = phi_eff;
            }
        }
    }
    l1_phi.second = -std::numeric_limits < real > ::max();
    l2_phi.second = -std::numeric_limits < real > ::max();
    l3_phi.second = -std::numeric_limits < real > ::max();
    for (integer i = 0; i != loc.size(); ++i) {
        const real x = loc[i].first;
        const real rho = loc[i].second[rho_i];
        const real pot = loc[i].second[pot_i];
        real phi_eff = pot / ASSERT_POSITIVE(rho) - 0.5 * x * x * omega * omega;
        if (x > std::min(rho1_max.first, rho2_max.first)
            && x < std::max(rho1_max.first, rho2_max.first)) {
            if (phi_eff > l1_phi.second) {
                l1_phi.second = phi_eff;
                l1_phi.first = x;
            }
        } else if (std::abs(x) > std::abs(rho2_max.first) && x * rho2_max.first > 0.0) {
            if (phi_eff > l2_phi.second) {
                l2_phi.second = phi_eff;
                l2_phi.first = x;
            }
        } else if (std::abs(x) > std::abs(rho1_max.first)) {
            if (phi_eff > l3_phi.second) {
                l3_phi.second = phi_eff;
                l3_phi.first = x;
            }
        }
    }
}

typedef node_server::start_run_action start_run_action_type;
HPX_REGISTER_ACTION(start_run_action_type);

hpx::future<void> node_client::start_run(bool b) const {
    return hpx::async<typename node_server::start_run_action>(get_unmanaged_gid(), b);
}

void node_server::start_run(bool scf)
{
    timings_.times_[timings::time_regrid] = 0.0;
    timings::scope ts(timings_, timings::time_total);
    integer output_cnt;

    if (!hydro_on) {
        save_to_file("X.chk");
        diagnostics();
        return;
    }
    printf("%e %e\n", grid::get_A(), grid::get_B());
    if (scf) {
        run_scf();
        set_pivot();
        printf("Adjusting velocities:\n");
        auto diag = diagnostics();
        space_vector dv;
        dv[XDIM] = -diag.grid_sum[sx_i] / diag.grid_sum[rho_i];
        dv[YDIM] = -diag.grid_sum[sy_i] / diag.grid_sum[rho_i];
        dv[ZDIM] = -diag.grid_sum[sz_i] / diag.grid_sum[rho_i];
        this->velocity_inc(dv);
        save_to_file("scf.chk");
    }

    printf("Starting...\n");
    solve_gravity(false);
    int ngrids = regrid(me.get_gid(), false);

    real output_dt = opts.output_dt;

    printf("OMEGA = %e, output_dt = %e\n", grid::get_omega(), output_dt);
    real& t = current_time;
    integer step_num = 0;

    auto fut_ptr = me.get_ptr();
    node_server* root_ptr = fut_ptr.get();

    output_cnt = root_ptr->get_rotation_count() / output_dt;

    profiler_output(stdout);

    real bench_start, bench_stop;
    while (current_time < opts.stop_time) {
        if (step_num > opts.stop_step)
            break;

        auto time_start = std::chrono::high_resolution_clock::now();
        if (!opts.disable_output && root_ptr->get_rotation_count() / output_dt >= output_cnt) {
            //	if (step_num != 0) {

            char fname[33];    // 21 bytes for int (max) + some leeway
            sprintf(fname, "X.%i.chk", int(output_cnt));
            save_to_file(fname);

            sprintf(fname, "X.%i.silo", int(output_cnt));
            output(fname, output_cnt, false);

            //	SYSTEM(std::string("cp *.dat ./dat_back/\n"));
            //	}
            ++output_cnt;

        }
        if (step_num == 0) {
            bench_start = hpx::util::high_resolution_clock::now() / 1e9;
        }

        //	break;
        real dt = 0;

        integer next_step = (std::min)(step_num + refinement_freq(), opts.stop_step + 1);
        dt = step(next_step - step_num).get();

        real omega_dot = 0.0, omega = 0.0, theta = 0.0, theta_dot = 0.0;
        omega = grid::get_omega();
        if ((opts.problem == DWD) && (step_num % refinement_freq() == 0)) {
            auto diags = diagnostics();

            const real dx = diags.secondary_com[XDIM] - diags.primary_com[XDIM];
            const real dy = diags.secondary_com[YDIM] - diags.primary_com[YDIM];
            const real dx_dot = diags.secondary_com_dot[XDIM]
                - diags.primary_com_dot[XDIM];
            const real dy_dot = diags.secondary_com_dot[YDIM]
                - diags.primary_com_dot[YDIM];
            theta = atan2(dy, dx);
            omega = grid::get_omega();
            theta_dot = (dy_dot * dx - dx_dot * dy) / (dx * dx + dy * dy) - omega;
            const real w0 = grid::get_omega() * 100.0;
            const real theta_dot_dot = (2.0 * w0 * theta_dot + w0 * w0 * theta);
            omega_dot = theta_dot_dot;
            omega += omega_dot * dt;
//            omega_dot += theta_dot_dot*dt;
            grid::set_omega(omega);
        }
        double time_elapsed = std::chrono::duration_cast<std::chrono::duration<double>>(
            std::chrono::high_resolution_clock::now() - time_start).count();

        // run output on separate thread
        if (!opts.disable_output)
        {
            hpx::threads::run_as_os_thread([=]()
            {
                FILE* fp = fopen( "step.dat", "at");
                fprintf(fp, "%i %e %e %e %e %e %e %e %e %i\n",
                    int(next_step - 1), double(t), double(dt), time_elapsed, rotational_time,
                    theta, theta_dot, omega, omega_dot, int(ngrids));
                fclose(fp);
            });     // do not wait for it to finish
        }

        hpx::threads::run_as_os_thread([=]()
        {
            printf("%i %e %e %e %e %e %e %e %e\n", int(next_step - 1), double(t), double(dt),
                time_elapsed, rotational_time, theta, theta_dot, omega, omega_dot);
        });     // do not wait for output to finish

//		t += dt;
        step_num += next_step;

        if (step_num % refinement_freq() == 0) {
            ngrids = regrid(me.get_gid(), false);

            // run output on separate thread
            auto need_break = hpx::threads::run_as_os_thread([&]()
            {
                FILE* fp = fopen("profile.txt", "wt");
                profiler_output(fp);
                fclose(fp);

                //		set_omega_and_pivot();
                bench_stop = hpx::util::high_resolution_clock::now() / 1e9;
                if (scf || opts.bench) {
                    printf("Total time = %e s\n", double(bench_stop - bench_start));
                    FILE* fp = fopen("bench.dat", "at");
                    fprintf(fp, "%i %e\n", int(hpx::find_all_localities().size()),
                        double(bench_stop - bench_start));
                    fclose(fp);
                    return true;
                }
                return false;
            });
            if (need_break.get())
                break;
        }
        //		set_omega_and_pivot();
        if (scf) {
            bench_stop = hpx::util::high_resolution_clock::now() / 1e9;
             printf("Total time = %e s\n", double(bench_stop - bench_start));
            //	FILE* fp = fopen( "bench.dat", "at" );
            //	fprintf( fp, "%i %e\n", int(hpx::find_all_localities().size()), double(bench_stop - bench_start));
            //	fclose(fp);
            break;
        }
    }
    bench_stop = hpx::util::high_resolution_clock::now() / 1e9;
    {
        timings::scope ts(timings_, timings::time_compare_analytic);
        compare_analytic();
        if (!opts.disable_output)
            output("final.silo", output_cnt, true);
    }

    if( opts.bench ) {
        FILE* fp = fopen( "scaling.dat", "at");
        const auto nproc = hpx::find_all_localities().size();
        fprintf( fp, "%i %e\n", int(nproc), float(bench_stop - bench_start));
        fclose( fp );
    }
}

typedef node_server::step_action step_action_type;
HPX_REGISTER_ACTION(step_action_type);

hpx::future<void> node_client::step(integer steps) const {
    return hpx::async<typename node_server::step_action>(get_unmanaged_gid(), steps);
}

void node_server::refined_step() {
#if HPX_HAVE_ITTNOTIFY != 0 && !defined(HPX_HAVE_APEX)
    static hpx::util::itt::string_handle sh("node_server::refined_step");
    hpx::util::itt::task t(hpx::get_thread_itt_domain(), sh);
#endif

    timings::scope ts(timings_, timings::time_computation);
    const real dx = TWO * grid::get_scaling_factor() / real(INX << my_location.level());
    real cfl0 = cfl;

    // FIXME: is this correct? ('a' was never re-initialized for refined == true)
    real a = std::numeric_limits<real>::min();

    dt_ = cfl0 * dx / a;

    all_hydro_bounds();
    local_timestep_channels[NCHILD].set_value(dt_);
    auto dt_fut = global_timestep_channel.get_future();

#ifdef RADIATION
    dt_ = dt_fut.get();
    compute_radiation(dt_/2.0);
    all_hydro_bounds();
#endif

    for (integer rk = 0; rk < NRK; ++rk) {

        compute_fmm(DRHODT, false);

        compute_fmm(RHO, true);
        all_hydro_bounds();

    }
#ifdef RADIATION
    compute_radiation(dt_/2.0);
    all_hydro_bounds();
#else
    dt_ = dt_fut.get();
#endif

    update();
}

hpx::future<void> node_server::nonrefined_step() {
#if HPX_HAVE_ITTNOTIFY != 0 && !defined(HPX_HAVE_APEX)
    static hpx::util::itt::string_handle sh("node_server::nonrefined_step");
    hpx::util::itt::task t(hpx::get_thread_itt_domain(), sh);
#endif

    timings::scope ts(timings_, timings::time_computation);

    real cfl0 = cfl;
    dt_ = ZERO;

    all_hydro_bounds();

    grid_ptr->store();
    hpx::future<void> fut = hpx::make_ready_future();

    hpx::shared_future<real> dt_fut = global_timestep_channel.get_future();

    for (integer rk = 0; rk < NRK; ++rk) {

        fut = fut.then(
            hpx::util::annotated_function(
                [rk, cfl0, this, dt_fut](hpx::future<void> f)
                {
                    f.get();        // propagate exceptions

                    grid_ptr->reconstruct();
                    real a = grid_ptr->compute_fluxes();
#ifdef RADIATION
                    if( rk == 0 ) {
                        const real dx = TWO * grid::get_scaling_factor() /
                            real(INX << my_location.level());
                        dt_ = cfl0 * dx / a;
                        local_timestep_channels[NCHILD].set_value(dt_);
                        dt_ = dt_fut.get();
                    	compute_radiation(dt_/2.0);
                    	all_hydro_bounds();
                        grid_ptr->reconstruct();
                        grid_ptr->compute_fluxes();
                    }

                    hpx::future<void> fut_flux = exchange_flux_corrections();
#else
                    hpx::future<void> fut_flux = exchange_flux_corrections();

                    if (rk == 0) {
                        const real dx = TWO * grid::get_scaling_factor() /
                            real(INX << my_location.level());
                        dt_ = cfl0 * dx / a;
                        local_timestep_channels[NCHILD].set_value(dt_);
                    }
#endif

                    return fut_flux.then(
                        hpx::launch::async(hpx::threads::thread_priority_boost),
                        hpx::util::annotated_function(
                            [rk, this, dt_fut](hpx::future<void> f)
                            {
                                f.get();        // propagate exceptions

                                grid_ptr->compute_sources(current_time);
                                grid_ptr->compute_dudt();

                                compute_fmm(DRHODT, false);

                                if (rk == 0) {
                                    dt_ = dt_fut.get();
                                }
                                grid_ptr->next_u(rk, current_time, dt_);

                                compute_fmm(RHO, true);
                                all_hydro_bounds();

#ifdef RADIATION
                                if(rk == NRK - 1) {
                                	all_hydro_bounds();
                                }
#endif
                            }, "node_server::nonrefined_step::compute_fmm"
                        ));
                }, "node_server::nonrefined_step::compute_fluxes"
            )
        );
    }

    return fut.then(hpx::launch::sync,
        [this](hpx::future<void>&& f)
        {
            f.get(); // propagate exceptions...
            update();
        }
    );
}

void node_server::update()
{
    grid_ptr->dual_energy_update();
    current_time += dt_;
    if (grid::get_omega() != 0.0)
    {
        rotational_time += grid::get_omega() * dt_;
    }
    else
    {
        rotational_time = current_time;
    }
}

hpx::future<real> node_server::step(integer steps) {
	grid_ptr->set_coordinates();

    std::array<hpx::future<void>, NCHILD> child_futs;
    if (is_refined)
    {
        for (integer ci = 0; ci != NCHILD; ++ci) {
            child_futs[ci] = children[ci].step(steps);
        }
    }

    hpx::future<real> dt_fut = hpx::make_ready_future(0.0);
    for (integer i = 0; i != steps; ++i)
    {
        dt_fut = dt_fut.then(
            [this, i, steps](hpx::future<real> dt_fut) -> hpx::future<real>
            {
                auto time_start = std::chrono::high_resolution_clock::now();
                auto next_dt = timestep_driver_descend();

                if (is_refined)
                {
                    refined_step();
                }
                else
                {
                    nonrefined_step().get();
                }

                real dt = dt_fut.get();
                if (my_location.level() == 0)
                {
                    double time_elapsed = std::chrono::duration_cast<std::chrono::duration<double>>(
                        std::chrono::high_resolution_clock::now() - time_start).count();

                    if (i + 1 != steps)
                    {
                        hpx::threads::run_as_os_thread([=]()
                        {
                            printf("%i %e %e %e %e\n", int(step_num), double(current_time), double(dt),
                                time_elapsed, rotational_time);
                        });     // do not wait for output to finish
                    }
                }
                ++step_num;

                return next_dt;
            });
    }

    if (is_refined)
    {
        return hpx::dataflow(hpx::launch::sync,
            [this](hpx::future<real> dt_fut, hpx::future<void>&& f)
            {
                f.get(); // propagate exceptions
                return dt_fut;
            },
            std::move(dt_fut),
            hpx::when_all(std::move(child_futs))
            );
    }
    else
    {
        return dt_fut;
    }
}

typedef node_server::timestep_driver_ascend_action timestep_driver_ascend_action_type;
HPX_REGISTER_ACTION(timestep_driver_ascend_action_type);

void node_client::timestep_driver_ascend(real dt) const {
    hpx::apply<typename node_server::timestep_driver_ascend_action>(get_unmanaged_gid(), dt);
}

void node_server::timestep_driver_ascend(real dt) {
    global_timestep_channel.set_value(dt);
    if (is_refined) {
        for(auto& child: children) {
            child.timestep_driver_ascend(dt);
        }
    }
}

typedef node_server::set_local_timestep_action set_local_timestep_action_type;
HPX_REGISTER_ACTION(set_local_timestep_action_type);

void node_client::set_local_timestep(integer idx, real dt) const {
    hpx::apply<typename node_server::set_local_timestep_action>(get_unmanaged_gid(), idx, dt);
}

void node_server::set_local_timestep(integer idx, real dt)
{
    local_timestep_channels[idx].set_value(dt);
}

hpx::future<real> node_server::timestep_driver_descend() {
    if (is_refined) {
        std::array<hpx::future<real>, NCHILD+1> futs;
        integer index = 0;
        for(auto& local_timestep: local_timestep_channels)
        {
            futs[index++] = local_timestep.get_future();
        }

        return hpx::dataflow(hpx::launch::sync,
            hpx::util::annotated_function(
                [this](std::array<hpx::future<real>, NCHILD+1> dts_fut) -> double
                {
                    auto dts = hpx::util::unwrapped(dts_fut);
                    real dt = *std::min_element(dts.begin(), dts.end());

                    if (my_location.level() == 0)
                    {
                        timestep_driver_ascend(dt);
                    }
                    else
                    {
                        parent.set_local_timestep(my_location.get_child_index(), dt);
                    }

                    return dt;
                },
                "node_server::timestep_driver_descend"),
            futs);
    } else {
        return local_timestep_channels[NCHILD].get_future().then(hpx::launch::sync,
            [this](hpx::future<real>&& f)
            {
                real dt = f.get();
                parent.set_local_timestep(my_location.get_child_index(), dt);
                return dt;
            });
    }
}

typedef node_server::velocity_inc_action velocity_inc_action_type;
HPX_REGISTER_ACTION(velocity_inc_action_type);

hpx::future<void> node_client::velocity_inc(const space_vector& dv) const {
    return hpx::async<typename node_server::velocity_inc_action>(get_gid(), dv);
}

void node_server::velocity_inc(const space_vector& dv) {
    if (is_refined) {
        std::array<hpx::future<void>, NCHILD> futs;
        integer index = 0;
        for (auto& child : children) {
            futs[index++] = child.velocity_inc(dv);
        }
        wait_all_and_propagate_exceptions(futs);
    } else {
        grid_ptr->velocity_inc(dv);
    }
}

