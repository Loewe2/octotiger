#include "geometry.hpp"
#include "interaction_types.hpp"
#include "interactions_iterators.hpp"

extern std::vector<interaction_type> ilist_debugging;

extern std::vector<interaction_type> ilist_n;
extern std::array<std::vector<interaction_type>, geo::direction::count()> ilist_n0_bnd;

void compare_interaction_lists() {
    if (ilist_debugging.size() != ilist_n.size() + ilist_n0_bnd.size()) {
        std::cout << "error: sizes of interaction lists do not match" << std::endl;
    }

    // std::cout << "check which interactions of new kernel cannot be found in either old list"
    //           << std::endl;
    // for (interaction_type& mine : ilist_debugging) {
    //     bool found = false;
    //     for (interaction_type& ref : ilist_n) {
    //         if (mine.first == ref.first && mine.second == ref.second) {
    //             found = true;
    //             std::cout << "found inner!" << std::endl;
    //             break;
    //         }
    //     }
    //     if (found) {
    //         continue;
    //     }
    //     for (geo::direction& dir : geo::direction::full_set()) {
    //         for (interaction_type& ref : ilist_n0_bnd[dir]) {
    //             if (mine.first == ref.first && mine.second == ref.second) {
    //                 found = true;
    //                 std::cout << "found boundary!" << std::endl;
    //                 break;
    //             }
    //         }
    //     }
    //     if (found) {
    //         continue;
    //     }

    //     octotiger::fmm::multiindex i =
    //     octotiger::fmm::inner_flat_index_to_multiindex(mine.first);
    //     octotiger::fmm::multiindex partner =
    //         octotiger::fmm::inner_flat_index_to_multiindex(mine.second);
    //     octotiger::fmm::multiindex diff(partner.x - i.x, partner.y - i.y, partner.z - i.z);
    //     std::cout << "interaction not found i_f: " << mine.first << " partner_f: " << mine.second
    //               << " i: " << i << " partner: " << partner << " diff: " << diff << std::endl;
    // }

    std::cout << "check whether old inner list -> in new list" << std::endl;
    for (interaction_type& mine : ilist_n) {
        octotiger::fmm::multiindex i =
            octotiger::fmm::flat_index_to_multiindex_not_padded(mine.first);
        octotiger::fmm::multiindex partner =
            octotiger::fmm::flat_index_to_multiindex_not_padded(mine.second);
        octotiger::fmm::multiindex diff(partner.x - i.x, partner.y - i.y, partner.z - i.z);

        bool found = false;
        for (interaction_type& ref : ilist_debugging) {
            octotiger::fmm::multiindex ref_first_padded =
                octotiger::fmm::flat_index_to_multiindex_padded(ref.first);
            octotiger::fmm::multiindex ref_second_padded =
                octotiger::fmm::flat_index_to_multiindex_padded(ref.second);
            octotiger::fmm::multiindex ref_first(
                ref_first_padded.x - 8, ref_first_padded.y - 8, ref_first_padded.z - 8);
            octotiger::fmm::multiindex ref_second(
                ref_second_padded.x - 8, ref_second_padded.y - 8, ref_second_padded.z - 8);
            size_t ref_first_flat = octotiger::fmm::to_inner_flat_index_not_padded(ref_first);
            size_t ref_second_flat = octotiger::fmm::to_inner_flat_index_not_padded(ref_second);

            if (mine.first == ref_first_flat && mine.second == ref_second_flat) {
                found = true;
                if (i.x == 0 && i.y == 0 && i.z == 1) {
                    octotiger::fmm::multiindex partner_debug(
                        i.x + ref.x[0], i.y + ref.x[1], i.z + ref.x[2]);
                    octotiger::fmm::multiindex diff_debug(ref.x[0], ref.x[1], ref.x[2]);
                    const integer second_check = gindex((ref_second_padded.x + INX) % INX,
                        (ref_second_padded.y + INX) % INX, (ref_second_padded.z + INX) % INX);
                    std::cout << "found!"
                              << " i: " << i << " partner: " << partner << " diff: " << diff
                              << " diff_len: " << diff.length() << std::endl;
                    std::cout << "diff_debug: " << diff_debug << std::endl;
                    if (mine.second != second_check) {
                        throw;
                    }
                }
                break;
            }
        }
        if (found) {
            continue;
        }

        if (i.x == 0 && i.y == 0 && i.z == 1) {
            std::cout << "old inner interaction not found o_f: " << mine.first
                      << " o_partner_f: " << mine.second << " i: " << i << " partner: " << partner
                      << " diff: " << diff << " diff_len: " << diff.length() << std::endl;
        }
    }

    // std::cout << "check whether old boundary list -> in new list" << std::endl;
    // for (geo::direction& dir : geo::direction::full_set()) {
    //     for (interaction_type& mine : ilist_n0_bnd[dir]) {
    //         bool found = false;
    //         for (interaction_type& ref : ilist_debugging) {
    //             if (mine.first == ref.first && mine.second == ref.second) {
    //                 found = true;
    //                 std::cout << "found!" << std::endl;
    //                 break;
    //             }
    //         }
    //         if (found) {
    //             continue;
    //         }

    //         std::cout << "old boundary interaction not found o_f: " << mine.first
    //                   << " o_partner_f: " << mine.second << std::endl;
    //     }
    // }
}
