/**
    NCG++
    rng.hpp

    Purpose:
    Provide an interface for a
    Random Number Generator with
    all the functionality needed.

    @author Peter Labus
    @version 0.1
    08.05.2017
*/

#pragma once

template<typename FT>
class RNG
{
  public:

    virtual FT uniform(FT lower, FT uppper);
    virtual FT signed_uniform(FT symetric_boundary);
    virtual int uniform_int(int lower, int upper);
};
