/**
 * \file coord.hpp
 * \brief Provides a coordinate struct
 *
 * \author Hayden Blauzvern
 */


#ifndef COORD_HPP_INCLUDED
#define COORD_HPP_INCLUDED

#include <array>
#include <iostream>

struct Coord {

    int x_, y_, z_;

    Coord() = delete;

    // Constructor
    Coord(const int x, const int y, const int z) :
        x_(x), y_(y), z_(z) {
        // Nothing
    }

    std::array<uint64_t, 3> get_rootkey(uint64_t level2_slog2) const {
        return {{x_ & ~((1LLU << level2_slog2) - 1),
                 y_ & ~((1LLU << level2_slog2) - 1),
                 z_ & ~((1LLU << level2_slog2) - 1)}};
    }

    uint64_t get_compressed_coord(uint64_t compressed) const {
        // Divide number of bits by number of axes to get compressed size
        uint64_t size_of_axis = 20 - compressed; // Likely 8
        uint64_t amount_to_pad = 64 - (size_of_axis * 3); // Likely 40
        uint64_t padding = (1LLU << amount_to_pad) - 1;
        // Each coordinate is at most 2^20, so right-shift by compressed
        return ((x_ >> compressed) << (amount_to_pad + size_of_axis * 2)) +
               ((y_ >> compressed) << (amount_to_pad + size_of_axis)) +
               ((z_ >> compressed) << amount_to_pad) +
               padding;
    }

    static uint64_t hash(std::array<uint64_t, 3> rootkey,
                         uint64_t hashmap_log_size) {
        return ((1LLU << hashmap_log_size) - 1) &
                    (rootkey[0] * 73856093LLU ^
                     rootkey[1] * 19349663LLU ^
                     rootkey[2] * 83492791LLU);
    }
};

#endif
