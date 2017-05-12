/**
    NCG++
    basic_maths.cpp

    Purpose:
    Provided some simple maths utility functions
    and possibly wrappers


    @author Peter Labus
    @version 0.1
    12.05.2017
*/

#include <vector>
#include <cstdint>
#include "basic_maths.hpp"



uint64_t factorial(uint64_t const n)
{
  if( n == 0 || n == 1 ) return 1;
  else return ( n * factorial(n-1) );
}


uint64_t binomial(uint64_t const a, uint64_t const b)
{
  return factorial(a) / ( factorial(b) * factorial(a-b) );
}


// FIXME: This is hard to understand legacy code.
// Can we make it more understandable?
std::vector<int> combination(int const upper, int const num_elems, int const num_comb)
{
  /**
   *  \brief Generate [num_comb]th combination of [num_elems]
   *         in the range of integers [ 0, 1, ..., upper - 1 ]
   *
   *  [num+comb] runs from 0 ... (upper choose num_elems) - 1
   */

  if( num_elems == 0 ) return std::vector<int>(0);

  std::vector<int> combinations(num_elems, 0);
  int k = 0;
  int r = 0;

  for(auto i = 0; i < num_elems - 1; ++i)
  {
    combinations[i] = (i != 0) ? combinations[i-1] : 0;

    do
    {
      ++combinations[i];
      r = binomial( upper - combinations[i], num_elems-(i+1) );
      k = k + r;
    } while(k <= num_comb);

    k = k - r;
  }

  combinations[num_elems-1] = combinations[num_elems-2] + num_comb + 1 - k;

  // FIXME: Brute-Force solution to decrease all values by one
  for(auto& c : combinations) --c;

  return combinations;
}
