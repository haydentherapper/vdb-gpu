/**
 * \file coord.hpp
 * \brief Provides a coordinate struct
 *
 * \author Hayden Blauzvern
 */


#ifndef COORD_HPP_INCLUDED
#define COORD_HPP_INCLUDED

#include <array>

struct Coord {

    int x_, y_, z_;

    Coord() = delete;

    // Constructor
    Coord(const int x, const int y, const int z) :
        x_(x), y_(y), z_(z) {
        // Nothing
    }

    std::array<int, 3> get_rootkey(uint64_t level2_slog2) const {
        return {{x_ & ~((1 << level2_slog2) - 1),
                 y_ & ~((1 << level2_slog2) - 1),
                 z_ & ~((1 << level2_slog2) - 1)}};
    }

    static uint64_t hash(std::array<int, 3> rootkey,
                         uint64_t hashmap_log_size) {
        return ((1 << hashmap_log_size) - 1) &
                    (rootkey[0] * 73856093 ^
                     rootkey[1] * 19349663 ^
                     rootkey[2] * 83492791);
    }
};

#endif
