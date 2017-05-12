/**
    NCG++
    basic_maths.hpp

    Purpose:
    Provided some simple maths utility functions
    and possibly wrappers


    @author Peter Labus
    @version 0.1
    12.05.2017
*/

#pragma once

#include <vector>



// Calculates n!
uint64_t factorial(uint64_t const n);

// Calculates a (choose) b
uint64_t binomial(uint64_t const a, uint64_t const b);

/**
 *  \brief Generate [num_comb]th combination of [num_elems]
 *         in the range of integers [ 0, 1, ..., upper - 1 ]
 *
 *  [num+comb] runs from 0 ... (upper choose num_elems) - 1
 */
std::vector<int> combination(int const upper, int const num_elems, int const num_comb);
