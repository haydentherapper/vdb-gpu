/**
 * \file utils.hpp
 * \brief
 *
 * \author Hayden Blauzvern
 */

#ifndef UTILS_HPP_INCLUDED
#define UTILS_HPP_INCLUDED

#include <stdint.h>

bool extract_bit(uint64_t index, uint64_t offset) {
    uint64_t constant = 1LLU << offset;
    // return if the nth bit is set
    return index & constant;
}

#endif /* end of include guard: UTILS_HPP_INCLUDED */
