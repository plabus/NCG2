
/**
    NCG++
    model_parameters.hpp

    Purpose:
    Structure all parameters needed for
    the matrix models.

    @author Peter Labus
    @version 0.1
    08.05.2017
*/

#pragma once

#include <cmath>

class ModelParameters
{
  public:

    explicit ModelParameters(
        int P,
        int Q,
        int N
    ) :
      p_(P),
      q_(Q),
      n_(N),
      d_( p_ + q_ ),
      s_( (q_-p_+8) % 8 ),
      k_( d_%2 ? (int)pow(2,(d_-1)/2) : (int)pow(2,d_/2) )
    {}

    int p(void) const { return p_; }
    int q(void) const { return q_; }
    int n(void) const { return n_; }
    int d(void) const { return d_; }
    int s(void) const { return s_; }
    int k(void) const { return k_; }


  private:

    const int p_;
    const int q_;
    const int n_;

    const int d_;
    const int s_;
    const int k_;
};
