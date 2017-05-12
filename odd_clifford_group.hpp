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

    const ModelParameters pqn_;
    std::vector<GammaMatrix> Gammas_;
};



// Non-Member Functions:
// =====================

// Generates the Odd Clifford Group
//   \Gamma_\mu = { \gamma_mu, \gamma_{\mu \nu \rho}, ... }
// with signature (p,q)
std::vector<GammaMatrix> generate_odd_clifford_group(
    const ModelParameters pqn
);

// Generates the Clifford Algebra { \gamma_\mu } with signature (p,q)
std::vector<GammaMatrix> generate_small_gammas(
    const ModelParameters pqn
);

// Generates the Euclidean Clifford Algebra { \gamma_\mu } with signature (d=p+q,0)
std::vector<GammaMatrix> generate_euclidean_gammas(
    const int d,                                        // dimensionality of Clifford algebra = # of gamma's
    const ModelParameters pqn = ModelParameters(0,0,0)  // provides gamma5_prefactor for odd dimensions,
                                                        // set to zero by default
);


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
