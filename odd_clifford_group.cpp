/**
    NCG++
    odd_clifford_group.cpp

    Purpose:
    Given the model parameters (p, q) initialise
    the entire odd big gamma matrices

    \Gamma_\mu = { \gamma_mu, \gamma_{\mu \nu \rho}, ... }

    and provide a print / read function.


    @author Peter Labus
    @version 0.1
    11.05.2017
*/

#include <iostream>
#include <iterator>
#include <vector>
#include <utility>
#include "model_parameters.hpp"
#include "basic_maths.hpp"
#include "gamma_matrix.hpp"
#include "clifford_algebra.hpp"
#include "odd_clifford_group.hpp"



// Member Functions:
// =================

OddCliffordGroup::OddCliffordGroup(ModelParameters const pqn)
  :
    pqn_(pqn),
    Gammas_(generate_odd_clifford_group_()),
    num_matrices_(count_Hs_and_Ls(Gammas_))
{
}


std::vector<GammaMatrix> OddCliffordGroup::generate_odd_clifford_group_()
{
  std::vector<GammaMatrix> big_gammas;
  const CliffordAlgebra small_gammas(pqn_);
  const auto d = small_gammas.size();

  for(auto num_indices = 1; num_indices <= d; num_indices += 2) // odd numbers of indices
  {
    // Calculate number of Gamma matrices with fixed number of indices:
    //   # Gamma matrices = (d choose num_indices)
    const auto num_matrices = binomial(d, num_indices);

    // Iterations over big_gammas with fixed number of indices
    for(auto num_comb = 0; num_comb < num_matrices; ++num_comb)
    {
      // 1. Generate the [num_comb]th combination with num_indices elements
      //    out of the range [0, 1, ..., d-1]. This will represent the
      //    indices of the antisymmetric product of small gamma matrices.
      // 2. Calculate the antisymmetric product and add it to the big
      //    gamma matrices.
      const auto index_sequence = combination(d, num_indices, num_comb);
      const auto matrix = small_gammas.antisymmetric_product(index_sequence);
      big_gammas.push_back(matrix);
    }
  }

  reshuffle_gammas( big_gammas );

  return big_gammas;
}



// Non-Member Functions:
// =====================

std::ostream& operator<<(std::ostream& os, OddCliffordGroup const& A)
{
  for(auto i= 0; i < A.size(); ++i)
  {
    os << " \u0393_" << i+1 << ":" << std::endl;
    os << A.Gamma(i) << std::endl;
  }
  return os;
}


void reshuffle_gammas(std::vector<GammaMatrix>& gammas)
{
  /**
   *  \brief Reshuffles a vector of GammaMatrices such that
   *         all hermitian matrices come first and all anti-hermitian
   *         matrices come second
   */

  std::vector<GammaMatrix> herm;
  std::vector<GammaMatrix> anti;

  for( auto const& gamma : gammas )
  {
    if( is_hermitian(gamma) )
    {
      herm.push_back(gamma);
    }
    else if( is_antihermitian(gamma) )
    {
      anti.push_back(gamma);
    }
    else
    {
      std::cout << " reshuffle_gammas :: ERROR: non-valid matrix " << std::endl;
    }
  }

  gammas.clear();
  gammas.insert(std::end(gammas), std::begin(herm), std::end(herm));
  gammas.insert(std::end(gammas), std::begin(anti), std::end(anti));
}


std::pair<int,int> count_Hs_and_Ls(std::vector<GammaMatrix> const& gammas)
{
  /**
   *  \brief Count the number of H's and L's in a vector of GammaMatrices
   *         assuming all are either hermitian or antihermitian and
   *         the hermitian matrices come first
   */

  int num_H = 0;
  for( auto const& gamma : gammas )
  {
    if( is_hermitian(gamma) ) ++num_H;
    else                      break;
  }

  int num_L = gammas.size() - num_H;
  return std::make_pair( num_H, num_L );
}
