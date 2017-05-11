
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

#include <complex>
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
      k_( d_%2 ? (int)pow(2,(d_-1)/2) : (int)pow(2,d_/2) ),
      exponent_ ( (s_ * (s_+1)/2) % 4 )
    {
      const std::complex<int> I({0,1});
      if(exponent_ == 0)      gamma5_prefactor_ =  1; // I^(4n+0)
      else if(exponent_ == 1) gamma5_prefactor_ =  I; // I^(4n+1)
      else if(exponent_ == 2) gamma5_prefactor_ = -1; // I^(4n+2)
      else if(exponent_ == 3) gamma5_prefactor_ = -I; // I^(4n+3)
    }

    int p(void) const { return p_; }
    int q(void) const { return q_; }
    int n(void) const { return n_; }
    int d(void) const { return d_; }
    int s(void) const { return s_; }
    int k(void) const { return k_; }
    std::complex<int> gamma5_prefactor(void) const { return gamma5_prefactor_; }


  private:

    const int p_;
    const int q_;
    const int n_;

    const int d_;
    const int s_;
    const int k_;
    const int exponent_;
    std::complex<int> gamma5_prefactor_;
};
