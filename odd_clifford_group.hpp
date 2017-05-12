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



// Non-Member Utility Functions:
// =============================

uint64_t factorial( uint64_t const n);

uint64_t binomial(uint64_t const a, uint64_t const b);

std::vector<int> combination(
    int const upper,
    int const num_elems,
    int const num_comb
);

void count_Hs_and_Ls(
    std::vector<int> const& sequence,
    int const p,
    int const q,
    int& num_H,
    int& num_L
);

GammaMatrix antisymmetrise(
    std::vector<GammaMatrix> const& gammas,
    std::vector<int> const& sequence
);
