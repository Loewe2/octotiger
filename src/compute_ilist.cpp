#include "defs.hpp"
#include "grid.hpp"
#include "options.hpp"
#include "simd.hpp"

// non-root and non-leaf node (David)
extern std::vector<interaction_type> ilist_n;
// monopole-monopole interactions? (David)
extern std::vector<interaction_type> ilist_d;
// interactions of root node (empty list? boundary stuff?) (David)
extern std::vector<interaction_type> ilist_r;
extern std::vector<std::vector<boundary_interaction_type>> ilist_d_bnd;
extern std::vector<std::vector<boundary_interaction_type>> ilist_n_bnd;
extern std::vector<std::vector<boundary_interaction_type>> ilist_n_bnd_new;

extern options opts;

// calculates 1/distance between i and j
real reciprocal_distance(const integer i0, const integer i1, const integer i2, const integer j0,
    const integer j1, const integer j2) {
    real tmp = (sqr(i0 - j0) + sqr(i1 - j1) + sqr(i2 - j2));
    // protect against sqrt(0)
    if (tmp > 0.0) {
        return 1.0 / (std::sqrt(tmp));
    } else {
        return 1.0e+10;    // dummy value
    }
}

// checks whether the index tuple is inside the current node
bool is_interior(const integer i0, const integer i1, const integer i2) {
    bool check = true;
    if (i0 < 0 || i0 >= G_NX) {
        check = false;
    } else if (i1 < 0 || i1 >= G_NX) {
        check = false;
    } else if (i2 < 0 || i2 >= G_NX) {
        check = false;
    }
    return check;
}

geo::direction get_neighbor_dir(const integer i, const integer j, const integer k) {
    integer i0 = 0, j0 = 0, k0 = 0;
    if (i < 0) {
        i0 = -1;
    } else if (i >= G_NX) {
        i0 = +1;
    }
    if (j < 0) {
        j0 = -1;
    } else if (j >= G_NX) {
        j0 = +1;
    }
    if (k < 0) {
        k0 = -1;
    } else if (k >= G_NX) {
        k0 = +1;
    }
    geo::direction d;
    d.set(i0, j0, k0);
    return d;
}

void compute_ilist() {
    // temporaries for insertion, why are there two of them?
    interaction_type np;    // M2M
    interaction_type dp;    // non-M2M
    // temporaries that store the boundary interactions in their non-converted form (see below)
    std::array<std::vector<interaction_type>, geo::direction::count()> ilist_n0_bnd;
    std::array<std::vector<interaction_type>, geo::direction::count()> ilist_d0_bnd;

    const real theta0 = opts.theta;    // TODO: verify that this is the radius of the sphere
    integer width = INX;
    // iterate node inner points
    for (integer i0 = 0; i0 != INX; ++i0) {
        for (integer i1 = 0; i1 != INX; ++i1) {
            for (integer i2 = 0; i2 != INX; ++i2) {
                // lower and upper boundaries for iteration
                // gives inner cube of interaction candidates plus bounding box where the boundaries
                // of the node are left
                // (iterate points in cube, iterate cube around every inner point and bounding box
                // of out-of-bounds points)
                integer i0_lower_bound = std::min(i0 - width, integer(0));
                integer i1_lower_bound = std::min(i1 - width, integer(0));
                integer i2_lower_bound = std::min(i2 - width, integer(0));
                integer i0_upper_bound = std::max(i0 + width + 1, G_NX);
                integer i1_upper_bound = std::max(i1 + width + 1, G_NX);
                integer i2_upper_bound = std::max(i2 + width + 1, G_NX);
                for (integer j0 = i0_lower_bound; j0 < i0_upper_bound; ++j0) {
                    for (integer j1 = i1_lower_bound; j1 < i1_upper_bound; ++j1) {
                        for (integer j2 = i2_lower_bound; j2 < i2_upper_bound; ++j2) {
                            if (i0 == j0 && i1 == j1 && i2 == j2) {
                                continue;
                            }
                            // in here: now i and j point to possible interaction partners
                            // calculate distance between inner nodes
                            const real x = i0 - j0;
                            const real y = i1 - j1;
                            const real z = i2 - j2;

                            // sqrt(0) cannot happen, as i != j
                            const real r = std::sqrt(sqr(x) + sqr(y) + sqr(z));
                            const real r3 = r * r * r;
                            // TODO: what is the purpose of four (David), speculation:
                            // - precomputation for interaction calculation
                            // - uses a little less data compared to storing i and j
                            v4sd four;

                            four[0] = -1.0 / r;
                            four[1] = x / r3;
                            four[2] = y / r3;
                            four[3] = z / r3;

                            // speculation: could be the inner sphere cutoff?
                            // arithmetic for round towards negative infinitive (relevant for
                            // negative numbers)
                            // i / 2 -> round-to-zero
                            // indices on a coarser level?
                            const integer i0_c = (i0 + INX) / 2 - INX / 2;
                            const integer i1_c = (i1 + INX) / 2 - INX / 2;
                            const integer i2_c = (i2 + INX) / 2 - INX / 2;

                            const integer j0_c = (j0 + INX) / 2 - INX / 2;
                            const integer j1_c = (j1 + INX) / 2 - INX / 2;
                            const integer j2_c = (j2 + INX) / 2 - INX / 2;

                            const real theta_f = reciprocal_distance(i0, i1, i2, j0, j1, j2);
                            const real theta_c =
                                reciprocal_distance(i0_c, i1_c, i2_c, j0_c, j1_c, j2_c);
                            // linear index of cell TODO move to outer loop?
                            const integer iii0 = gindex(i0, i1, i2);
                            // calculates linear index in same or neighbor cell
                            // (+ INX to map to positive range for mod)
                            const integer iii1n =
                                gindex((j0 + INX) % INX, (j1 + INX) % INX, (j2 + INX) % INX);
                            // linear index of interaction partner
                            // (only used for saving work using the reflective property)
                            const integer iii1 = gindex(j0, j1, j2);
                            // not in inner sphere (theta_c > theta0), but in outer sphere
                            if (theta_c > theta0 && theta_f <= theta0) {
                                np.first = iii0;
                                np.second = iii1n;
                                np.four = four;
                                np.x[XDIM] = j0;
                                np.x[YDIM] = j1;
                                np.x[ZDIM] = j2;
                                // TODO: i should always be in interior?
                                if (is_interior(j0, j1, j2)) {
                                    ilist_n.push_back(np);
                                } else {
                                    integer neighbor_index =
                                        get_neighbor_dir(j0, j1, j2).flat_index();
                                    ilist_n0_bnd[neighbor_index].push_back(np);
                                }
                            }
                            // not in inner sphere, and maybe not in outer sphere?
                            if (theta_c > theta0) {
                                dp.first = iii0;
                                dp.second = iii1n;
                                dp.x[XDIM] = j0;
                                dp.x[YDIM] = j1;
                                dp.x[ZDIM] = j2;
                                dp.four = four;
                                if (is_interior(j0, j1, j2)) {
                                    if (iii1 > iii0) {
                                        ilist_d.push_back(dp);
                                    }
                                } else {
                                    integer neighbor_index =
                                        get_neighbor_dir(j0, j1, j2).flat_index();
                                    ilist_d0_bnd[neighbor_index].push_back(dp);
                                }
                            }
                            // in outer sphere, and maybe not in inner sphere
                            if (theta_f <= theta0) {
                                np.first = iii0;
                                np.second = iii1n;
                                np.x[XDIM] = j0;
                                np.x[YDIM] = j1;
                                np.x[ZDIM] = j2;
                                np.four = four;
                                if (is_interior(j0, j1, j2)) {
                                    ilist_r.push_back(np);
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    // non_M2M d
    // M2M n

    // Probable functionality:
    // Uses the interaction_type-lists created above and converts them into a list of list
    // that contains the boundary interactions for each list entry
    // (0 are temporaries without the innerlists aggregated)

    for (auto& dir : geo::direction::full_set()) {    // iterate 26 neighbors
        std::vector<boundary_interaction_type>& non_M2M = ilist_d_bnd[dir];
        std::vector<interaction_type>& non_M2M_0 = ilist_d0_bnd[dir];
        std::vector<boundary_interaction_type>& M2M = ilist_n_bnd[dir];
        std::vector<boundary_interaction_type>& M2M_new = ilist_n_bnd_new[dir];
        std::vector<interaction_type>& M2M_0 = ilist_n0_bnd[dir];

        // add entry for every non-M2M boundary interaction
        for (interaction_type& interaction_0 : non_M2M_0) {
            boundary_interaction_type boundary_interaction;
            boundary_interaction.second.push_back(interaction_0.second);
            boundary_interaction.x = interaction_0.x;
            non_M2M.push_back(boundary_interaction);
        }
        // fill up list of actual non-M2M boundary interactions
        for (interaction_type& interaction_0 : non_M2M_0) {
            for (boundary_interaction_type& boundary_interaction : non_M2M) {
                if (boundary_interaction.second[0] == interaction_0.second) {
                    boundary_interaction.first.push_back(interaction_0.first);
                    boundary_interaction.four.push_back(interaction_0.four);
                    break;
                }
            }
        }

        // add entry for every M2M boundary interaction
        for (interaction_type& interaction_0 : M2M_0) {
            {
                boundary_interaction_type boundary_interaction;
                boundary_interaction.second.push_back(interaction_0.second);
                boundary_interaction.x = interaction_0.x;
                M2M.push_back(boundary_interaction);
            }

            {
                // new version
                boundary_interaction_type boundary_interaction;
                boundary_interaction.first.push_back(interaction_0.first);
                boundary_interaction.x = interaction_0.x;
                M2M_new.push_back(boundary_interaction);
            }
        }
        // fill up list of actual M2M boundary interactions
        for (interaction_type& interaction_0 : M2M_0) {
            for (boundary_interaction_type& boundary_interaction : M2M) {
                if (boundary_interaction.second[0] == interaction_0.second) {
                    boundary_interaction.first.push_back(interaction_0.first);
                    boundary_interaction.four.push_back(interaction_0.four);
                    break;
                }
            }

            // new version
            for (boundary_interaction_type& boundary_interaction : M2M_new) {
                if (boundary_interaction.first[0] == interaction_0.first) {
                    boundary_interaction.second.push_back(interaction_0.second);
                    boundary_interaction.four.push_back(interaction_0.four);
                    break;
                }
            }
            // boundary_interaction_type boundary_interaction;
            // boundary_interaction.first.push_back(interaction_0.first);
            // boundary_interaction.four.push_back(interaction_0.four);
            // boundary_interaction.second.push_back(interaction_0.second);
            // boundary_interaction.x = interaction_0.x;
            // M2M_new.push_back(boundary_interaction);
        }
    }

    {
        std::int32_t cur_index = ilist_r[0].first;
        std::int32_t cur_index_start = 0;
        for (std::int32_t li = 0; li != ilist_r.size(); ++li) {
            // std::cout << ilist_r[li].first << " ?? " << cur_index << std::endl;
            if (ilist_r[li].first != cur_index) {
                // std::cout << "range size: " << (li - cur_index_start) << std::endl;
                ilist_r[cur_index_start].inner_loop_stop = li;
                cur_index = ilist_r[li].first;
                cur_index_start = li;
            }
        }
        // make sure the last element is handled correctly as well
        ilist_r[cur_index_start].inner_loop_stop = ilist_r.size();
    }

    {
        std::int32_t cur_index = ilist_n[0].first;
        std::int32_t cur_index_start = 0;
        for (std::int32_t li = 0; li != ilist_n.size(); ++li) {
            // std::cout << ilist_n[li].first << " ?? " << cur_index << std::endl;
            if (ilist_n[li].first != cur_index) {
                // std::cout << "range size: " << (li - cur_index_start) << std::endl;
                ilist_n[cur_index_start].inner_loop_stop = li;
                cur_index = ilist_n[li].first;
                cur_index_start = li;
            }
        }
        // make sure the last element is handled correctly as well
        ilist_n[cur_index_start].inner_loop_stop = ilist_n.size();
    }

    {
        std::int32_t cur_index = ilist_d[0].first;
        std::int32_t cur_index_start = 0;
        for (std::int32_t li = 0; li != ilist_d.size(); ++li) {
            // std::cout << ilist_d[li].first << " ?? " << cur_index << std::endl;
            if (ilist_d[li].first != cur_index) {
                // std::cout << "range size: " << (li - cur_index_start) << std::endl;
                ilist_d[cur_index_start].inner_loop_stop = li;
                cur_index = ilist_d[li].first;
                cur_index_start = li;
            }
        }
        // make sure the last element is handled correctly as well
        ilist_d[cur_index_start].inner_loop_stop = ilist_d.size();
    }
}
