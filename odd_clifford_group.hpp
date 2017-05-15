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
#include <utility>
#include "model_parameters.hpp"
#include "gamma_matrix.hpp"



class OddCliffordGroup
{
  public:

    // Public Member Functions:
    // ========================

    // Constructors:
    OddCliffordGroup(void) = delete;
    explicit OddCliffordGroup(ModelParameters const pqn);

    // Getters:
    GammaMatrix const& Gamma(int index) const { return Gammas_[index]; }
    int num_H(void) const { return num_matrices_.first; }
    int num_L(void) const { return num_matrices_.second; }
    int size(void) const { return Gammas_.size(); }

    // Setters:
    GammaMatrix& Gamma(int index) { return Gammas_[index]; }


  private:

    // Private Member Variables:
    // =========================
    ModelParameters const pqn_;
    std::vector<GammaMatrix> Gammas_;
    std::pair<int,int> const num_matrices_;

    // Private Member Functions:
    // =========================

    // Generates the Odd Clifford Group
    //   \Gamma_\mu = { \gamma_mu, \gamma_{\mu \nu \rho}, ... }
    // with signature (p,q)
    std::vector<GammaMatrix> generate_odd_clifford_group_(void);
};



// Non-Member Functions:
// =====================

std::ostream& operator<<(std::ostream& os, OddCliffordGroup const& A);

// Reshuffles a vector of GammaMatrices such that
// all hermitian matrices come first and all anti-hermitian
// matrices come second
void reshuffle_gammas(std::vector<GammaMatrix>& gammas);

// Count the number of H's and L's in a vector of GammaMatrices
// assuming all are either hermitian or antihermitian and
// the hermitian matrices come first
std::pair<int,int> count_Hs_and_Ls(std::vector<GammaMatrix> const& gammas);
