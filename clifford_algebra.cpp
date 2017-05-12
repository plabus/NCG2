/**
    NCG++
    clifford_algebra.cpp

    Purpose:
    Given the model parameters (p, q) initialise
    the small gamma matrices

      \gamma_\mu

    and provide functionality to calculate anti-symmetrised
    products given a vector of indices.


    @author Peter Labus
    @version 0.1
    12.05.2017
*/

#include <iostream>
#include <vector>
#include <cmath>
#include <cassert>
#include "model_parameters.hpp"
#include "basic_maths.hpp"
#include "gamma_matrix.hpp"
#include "clifford_algebra.hpp"



// Member Functions:
// =================

CliffordAlgebra::CliffordAlgebra(ModelParameters const pqn)
  :
    pqn_(pqn),
    gammas_(generate_small_gammas_())
{
}


GammaMatrix CliffordAlgebra::antisymmetric_product(std::vector<int> const& indices) const
{
  /**
   *  \brief Return antisymmetrised product of gamma matrices
   *         with the indices given.
   *
   *  1. Generate the sum of (signed) products of all
   *     possible permutations
   *  2. Normalise the sum by dividing through the number
   *     of permutations
   */

  auto matrix = antisymmetrise(gammas_, indices);
  auto const n = factorial( indices.size() );
  matrix = matrix / n;
  return matrix;
}


std::vector<GammaMatrix> CliffordAlgebra::generate_small_gammas_(void)
{
  /**
   *  \brief Generate all small gamma matrices.
   *
   *  Generate all small gamma matrices with signiture (d,0)
   *  and then multiply the last q matrices with I to
   *  get small gamma matrices of signature (p,q)
   */

  auto gammas = generate_euclidean_gammas(pqn_.d(), pqn_);

  const std::complex<int> I({0,1});
  for(auto i = pqn_.p(); i != gammas.size(); ++i)
  {
    gammas[i] *= I;
  }

  return gammas;
}



// Non-Member Functions:
// =====================

std::ostream& operator<<(std::ostream& os, CliffordAlgebra const& A)
{
  for(auto i= 0; i < A.size(); ++i)
  {
    os << " gamma " << i+1 << ":" << std::endl;
    os << A.gammas_[i] << std::endl;
  }
  return os;
}


std::vector<GammaMatrix> generate_euclidean_gammas(
    int const d,                // dimensionality of Clifford algebra = # of gamma's
    ModelParameters const pqn   // provides gamma5_prefactor for odd dimensions,
                                // set to zero by default
)
{
  /**
   *  \brief Generate all small Euclidean gamma matrices in d dimensions.
   *
   *  FIXME: Add description of algorithm
   */

  const PauliMatrices Pauli;
  std::vector<GammaMatrix> gammas;

  if( d == 1 ) // Start of recursion
  {
    gammas.push_back( Unity(1) );
  }
  else if( d == 2 ) // Start of recursion
  {
    gammas.push_back( Pauli.sigma1 );
    gammas.push_back( Pauli.sigma2 );
  }
  else if( d > 2 && d % 2 == 0 ) // Even-dim case
  {
    /**
     *  Gamma matrices in d even dimensions
     *
     *  { gamma_mu (x) sigma1, 1I (x) sigma2, 1I (x) sigma3 },
     *
     *  where gamma_mu and 1I are the matrices in (d-2) dimensions,
     *  i.e. k x k square complex matrices with
     *
     *    k = 2^{(d-2)/2}
     *
     */
    auto d_recursive = d-2;
    auto k_recursive = static_cast<int>( pow(2, d_recursive/2) );
    auto small_gammas = generate_euclidean_gammas(d_recursive);
    auto one = Unity(k_recursive);

    // The following are outer (tensor) products
    for(auto gamma : small_gammas)
    {
      gammas.push_back( gamma % Pauli.sigma1 );
    }
    gammas.push_back( one % Pauli.sigma2 );
    gammas.push_back( one % Pauli.sigma3 );
  }
  else // Odd-dim case
  {
    /**
     *  Gamma matrices in d odd dimensions
     *
     *    { gamma_mu, gamma_5 },
     *
     *  where gamma_mu are the matrices d-1.
     */
    gammas = generate_euclidean_gammas(d-1);
    auto prefactor = pqn.gamma5_prefactor();
    auto gamma_5 = prefactor * gammas.front();

    for(auto i = 1; i != gammas.size(); ++i)
    {
      gamma_5 *= gammas[i];
    }
    gammas.push_back( gamma_5 );
  }

  return gammas;
}


GammaMatrix antisymmetrise(
    std::vector<GammaMatrix> const& gammas,
    std::vector<int> const& sequence
)
{
  // Base Variables and Assertions:
  // ==============================

  const auto num_indices = sequence.size();
  const auto d = gammas.size();
  const auto k = gammas.front().size();

  assert(num_indices <= d && "antisymmetrise: ERROR: Number of indices bigger than dimension!");
  for(auto i = 0; i < num_indices; ++i)
  {
    assert(sequence[i] <= d && "antisymmetrise: ERROR: Number of indices bigger than dimension!");
  }


  // Base Cases:
  // ===========

  if( num_indices == 0 )
  {
    return Unity(k);
  }
  else if( num_indices == 1 )
  {
    auto index = sequence[0];
    return gammas[index];
  }
  else if( num_indices == 2 )
  {
    auto index1 = sequence[0];
    auto index2 = sequence[1];
    return commutator(gammas[index1], gammas[index2]);
  }


  // Recursive Step:
  // ===============

  GammaMatrix matrix(k);

  // Iterate over all elements in sequence:
  //
  // 1. delete one element out of sequence and construct a new one
  //    sequence_new = [a1, a2, ..., (an), ..., a_num_indices]
  // 2. call antisymmetrise recursively and save in buffer1
  // 3. take the product between buffer1 and the remaining
  //    gamma matrix and save in buffer2
  // 4. matrix += (-1)^(pos-1) * buffer2

  for(auto i = 0; i < num_indices; ++i)
  {
    // 1. Copy new sequence without the ith element
    auto new_sequence = sequence;
    new_sequence.erase( new_sequence.begin() + i );

    // 2. Recursive step
    auto buffer1 = antisymmetrise(gammas, new_sequence);

    // 3. Take the product
    auto index = sequence[i];
    auto buffer2 = gammas[index] * buffer1;

    // 4. Add to the result matrix
    if(i % 2 == 0) matrix = matrix + buffer2;
    else           matrix = matrix - buffer2;
  }

  return matrix;
}
