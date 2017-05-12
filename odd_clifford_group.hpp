/**
    NCG++
    odd_clifford_group.hpp

    Purpose:
    Given the model parameters (p, q) initialise
    the entire odd big gamma matrices

    \Gamma_\mu = { \gamma_mu, \gamma_{\mu \nu \rho}, ... }

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



class OddCliffordGroup
{
  public:

    OddCliffordGroup(void) = delete;
    explicit OddCliffordGroup(ModelParameters const pqn);

    friend std::ostream& operator<<(std::ostream&, OddCliffordGroup const& A);


  private:

    ModelParameters const pqn_;
    std::vector<GammaMatrix> const Gammas_;

    // Generates the Odd Clifford Group
    //   \Gamma_\mu = { \gamma_mu, \gamma_{\mu \nu \rho}, ... }
    // with signature (p,q)
    std::vector<GammaMatrix> generate_odd_clifford_group_(void);
};
