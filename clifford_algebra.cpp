/**
    NCG++
    clifford_algebra.cpp

    Purpose:
    Given the model parameters (p, q) initialise
    the entire ODD gamma matrices

    \Gamma_\mu = { \gamma_mu, \gamma_{\mu \nu \rho}, ... }

    and provide a print / read function.


    @author Peter Labus
    @version 0.1
    11.05.2017
*/

#include <iostream>
#include <vector>
#include <cmath>
#include <cassert>
#include "model_parameters.hpp"
#include "gamma_matrix.hpp"
#include "clifford_algebra.hpp"



// Member Functions:
// =================

CliffordAlgebra::CliffordAlgebra(ModelParameters const pqn)
  :
    pqn_(pqn),
    // Gammas_(generate_small_gammas(pqn_))
    Gammas_(generate_odd_clifford_group(pqn_))
{
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


// TODO: Move to a better place
uint64_t factorial(uint64_t const n)
{
  if( n == 0 || n == 1 ) return 1;
  else return ( n * factorial(n-1) );
}


// TODO: Move to a better place
uint64_t binomial(uint64_t const a, uint64_t const b)
{
  return factorial(a) / ( factorial(b) * factorial(a-b) );
}


// FIXME: This is hard to understand legacy code.
// Can we make it more understandable?
std::vector<int> combination(int const upper, int const num_elems, int const num_comb)
{
  /**
   *  \brief Generate [num_comb]th combination of [num_elems]
   *         in the range of integers [ 1, 2, ..., upper ]
   *
   *  [num+comb] runs from 0 ... (upper choose num_elems) - 1
   */

  std::vector<int> combinations(num_elems, 0);
  int k = 0;
  int r = 0;

  for(auto i = 0; i < num_elems - 1; ++i)
  {
    combinations[i] = (i != 0) ? combinations[i-1] : 0;

    do
    {
      combinations[i]++;
      r = binomial( upper - combinations[i], num_elems-(i+1) );
      k = k + r;
    } while(k <= num_comb);

    k = k - r;
  }

  combinations[num_elems-1] = combinations[num_elems-2] + num_comb + 1 - k;
  return combinations;
}


void count_Hs_and_Ls(std::vector<int> const& sequence, int const p, int const q, int& num_H, int& num_L)
{
  // Assume n >= 1
  auto n = sequence.size();

  if(n == 1)
  {
    if(sequence[0] <= p) num_H += 1;
    else                 num_L += 1;
  }
  else
  {
    auto exp = ( n - 1 ) * n / 2;
    for(int i = 0; i < n; ++i) if(sequence[i] > p) exp++;

    if(exp%2==0) num_H += 1;
    else         num_L += 1;
  }
}


GammaMatrix antisymmetrise(
    std::vector<GammaMatrix> const& gammas,
    std::vector<int> const& sequence,
    int const d,
    int const k
)
{
  // Number of indices we have antisymmetrise over
  auto num_indices = sequence.size();

  assert(num_indices <= d && "antisymmetrise: ERROR: Number of indices bigger than dimension!");
  for(auto i = 0; i < num_indices; ++i)
  {
    assert(sequence[i] <= d && "antisymmetrise: ERROR: Number of indices bigger than dimension!");
  }

  // If there is one index return relevant gamma matrix
  if( num_indices == 1 )
  {
    /* Note n1 in [1,2,...,dim] */
    auto index = sequence[0] - 1;
    return gammas[index];
  }

  // If there are two indices calculate and return commutator
  if( num_indices == 2 )
  {
    auto index1 = sequence[0] - 1;
    auto index2 = sequence[1] - 1;
    return commutator(gammas[index1], gammas[index2]);
  }

  // If there are more than two indices anti-symmetrise recursively
  if( num_indices > 2 )
  {
    std::cout << " CAN'T DO THAT YET!" << std::endl;
  }
  // {
  //   /* Iterate over all elements in sequence:               *
  //    * 0. choose element n                                  *
  //    * 1. seq_new = [a1, a2, ..., (an), ..., a_num_indices] *
  //    * 2. call antisymmetrise recursively and save in buff1 *
  //    * 3. take the product between buff1 and the remaining  *
  //    *    gamma matrix and save in buff2                    *
  //    * 4. matrix += (-1)^(pos-1) * buff2                    *
  //    * 5. after iteration matrix = matrix / num_indices     */
  //    int new_seq[num_indices-1];
  //    float complex *buff1 = (float complex*) malloc(k*k*sizeof(float complex));
  //    float complex *buff2 = (float complex*) malloc(k*k*sizeof(float complex));

  //    for(int i=0;i<num_indices;++i) { /* Iterator over sequence */
  //      /* 0. Initialise buff2 to zero */
  //      for(int ii=0;ii<k*k;++ii) buff2[ii] = 0. + 0.*I;
  //      /* 1. Copy new sequence w/o ith element */
  //      index = 0;
  //      for(int j=0;j<num_indices;++j) {
  //        if(i!=j) {
  //          new_seq[index] = sequence[j];
  //          index++;
  //        }
  //      }
  //      /* 2. Recursive step */
  //      antisymmetrise(gammas, dim, k, num_indices-1, new_seq, buff1);
  //      /* 3. Take the product */
  //      for(int ii=0;ii<k;++ii) {
  //        for(int jj=0;jj<k;++jj) {
  //          for(int ll=0;ll<k;++ll) {
  //            buff2[ii*k+jj] += buff1[ii*k+ll]*gammas[ll*k+jj + (sequence[i]-1)*k*k];
  //          }
  //        }
  //      }
  //      /* 4. Add to the result matrix */
  //      for(int ii=0;ii<k*k;++ii)
  //        if(i%2==0) matrix[ii] += buff2[ii];
  //        else       matrix[ii] -= buff2[ii];
  //    }
  //    /* 5. Divide by #indices */
  //    for(int i=0;i<k*k;++i) matrix[i] /= num_indices;

  //    free(buff1);
  //    free(buff2);
  // }

  return GammaMatrix(k);
}


std::vector<GammaMatrix> generate_odd_clifford_group(
    const ModelParameters pqn
)
{
  const auto p = pqn.p();
  const auto q = pqn.q();
  const auto d = pqn.d();
  const auto k = pqn.k();
  int num_H = 0;
  int num_L = 0;

  std::vector<GammaMatrix> Gammas;
  const auto gammas = generate_small_gammas(pqn);

  for(auto num_indices = 1; num_indices <= d; num_indices += 2)
  {
    // Calculate number of Gamma matrices with fixed number of indices:
    //   # matrices = (d choose num_indices)
    auto num_matrices = binomial(d, num_indices);

    // Iterations over Gammas with fixed number of indices */
    for(auto num_comb = 0; num_comb < num_matrices; ++num_comb)
    {
      // 1. Generate the [num_comb]th combination with num_indices elements
      //    out of the range [1..d]
      auto const index_sequence = combination(d, num_indices, num_comb);
      count_Hs_and_Ls(index_sequence, p, q, num_H, num_L);
      auto matrix = antisymmetrise(gammas, index_sequence, d, k);
      Gammas.push_back(matrix);
    }
  }

  // TODO:
  // add reshuffling!
  // copy num_H & num_L into CliffordAlgebra class object
  std::cout << " H's = " << num_H << ", L's = " << num_L << std::endl;

  return Gammas;
}


std::vector<GammaMatrix> generate_small_gammas(ModelParameters const pqn)
{
  /**
   *  \brief Generate all small gamma matrices.
   *
   *  Generate all small gamma matrices with signiture (d,0)
   *  and then multiply the last q matrices with I to
   *  get small gamma matrices of signature (p,q)
   */

  auto gammas = generate_euclidean_gammas(pqn.d(), pqn);

  const std::complex<int> I({0,1});
  for(auto i = pqn.p(); i != gammas.size(); ++i)
  {
    gammas[i] *= I;
  }

  return gammas;
}


std::vector<GammaMatrix> generate_euclidean_gammas(
    const int d,                // dimensionality of Clifford algebra = # of gamma's
    const ModelParameters pqn   // provides gamma5_prefactor for odd dimensions,
                                // set to zero by default
)
{
  /**
   *  \brief Generate all small Euclidean gamma matrices in d dimensions.
   *
   *  Generate all small gamma matrices with
   *  signiture (d,0) recursively.
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
