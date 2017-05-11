/**
    NCG++
    clifford_algebra.hpp

    Purpose:
    Given the model parameters (p, q) initialise
    the entire ODD gamma matrices

    \Gamma_\mu = { 1, \gamma_mu, \gamma_{\mu \nu}, ..., \gamma_5 }

    and provide a print / read function.


    @author Peter Labus
    @version 0.1
    09.05.2017
*/

#pragma once

#include <iostream>
#include <vector>
#include "model_parameters.hpp"
#include "gamma_matrix.hpp"


// Forword declaration
std::vector<GammaMatrix> generate_gammas(const int d);


class CliffordAlgebra
{
  public:


    explicit CliffordAlgebra(ModelParameters const pqn)
      :
    pqn_(pqn),
    Gammas_(generate_gammas_())
    {
    }

    friend std::ostream& operator<<(std::ostream&, CliffordAlgebra const& A);

  private:

    const ModelParameters pqn_;
    std::vector<GammaMatrix> Gammas_;

    std::vector<GammaMatrix> generate_gammas_(void);
};


// Member Functions:
// =================

std::vector<GammaMatrix> CliffordAlgebra::generate_gammas_(void)
{
  const int d = pqn_.d();
  const int p = pqn_.p();
  const std::complex<int> I({0,1});

  // Generate all small gamma matrices with signiture (d,0)
  // and then multiply the last q matrices with I
  auto gammas = generate_gammas(d);
  for(auto i = p; i != gammas.size(); ++i) gammas[i] = I * gammas[i];
  return gammas;
}


// Non-Member Functions:
// =====================

std::ostream& operator<<(std::ostream& os, CliffordAlgebra const& A)
{
  for(auto i= 0; i < A.Gammas_.size(); ++i)
  {
    os << " Gamma " << i+1 << ":" << std::endl;
    os << A.Gammas_[i] << std::endl;
  }
  return os;
}

std::vector<GammaMatrix> generate_gammas(const int d)
{
  std::vector<GammaMatrix> gammas;
  PauliMatrices Pauli;

  if( d == 1 )
  {
    gammas.push_back( Unity(1) );
  }
  else if( d == 2 )
  {
    gammas.push_back( Pauli.sigma1 );
    gammas.push_back( Pauli.sigma2 );
  }
  else if( d > 2 && d % 2 == 0 )
  {
    auto small_gammas = generate_gammas(d-2);
    auto one = Unity(d-2);

    // The following are outer (tensor) products
    for(auto gamma : small_gammas)
    {
      gammas.push_back( gamma % Pauli.sigma1 );
    }
    gammas.push_back( one % Pauli.sigma2 );
    gammas.push_back( one % Pauli.sigma3 );
  }
  else
  {
    // FIXME: How to better calculate s?
    const int s = (8*8*8-d) % 8;
    const int exponent = ( s * (s+1)/2 ) % 4;
    const std::complex<int> I({0,1});

    std::complex<int> prefactor;
    if(exponent==0)      prefactor =  1; // I^(4n+0)
    else if(exponent==1) prefactor =  I; // I^(4n+1)
    else if(exponent==2) prefactor = -1; // I^(4n+2)
    else if(exponent==3) prefactor = -I; // I^(4n+3)

    gammas = generate_gammas(d-1);
    auto gamma_5 = prefactor * gammas.front();
    for(auto i = 1; i != gammas.size(); ++i) gamma_5 = gamma_5 * gammas[i];
    gammas.push_back( gamma_5 );
  }

  return gammas;
}
