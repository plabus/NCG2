/**
    NCG++
    clifford_algebra.hpp

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

#pragma once

#include <iostream>
#include <vector>
#include "model_parameters.hpp"
#include "gamma_matrix.hpp"



class CliffordAlgebra
{
  public:

    CliffordAlgebra(void) = delete;
    explicit CliffordAlgebra(ModelParameters const pqn);

    friend std::ostream& operator<<(std::ostream&, CliffordAlgebra const& A);


  private:

    const ModelParameters pqn_;
    std::vector<GammaMatrix> gammas_;

    // Generates the Clifford Algebra { \gamma_\mu } with signature (p,q)
    std::vector<GammaMatrix> generate_small_gammas_(void);
};



// Non-Member Functions:
// =====================

// Generates the Euclidean Clifford Algebra { \gamma_\mu } with signature (d=p+q,0)
std::vector<GammaMatrix> generate_euclidean_gammas(
    const int d,                                        // dimensionality of Clifford algebra = # of gamma's
    const ModelParameters pqn = ModelParameters(0,0,0)  // provides gamma5_prefactor for odd dimensions,
                                                        // set to zero by default
);


// Non-Member Utility Functions:
// =============================

GammaMatrix antisymmetrise(
    std::vector<GammaMatrix> const& gammas,
    std::vector<int> const& sequence
);
