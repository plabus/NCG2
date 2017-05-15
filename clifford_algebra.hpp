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

    // Public Member Functions:
    // ========================

    // Constructors:
    CliffordAlgebra(void) = delete;
    explicit CliffordAlgebra(ModelParameters const pqn);

    // Getters:
    GammaMatrix const& gamma(int index) const { return gammas_[index]; }
    int size(void) const { return gammas_.size(); }

    // Setters:
    GammaMatrix& gamma(int index) { return gammas_[index]; }

    // Functionality:
    GammaMatrix antisymmetric_product(std::vector<int> const& indices) const;


  private:

    // Private Member Variables:
    // =========================
    ModelParameters const pqn_;
    std::vector<GammaMatrix> gammas_;

    // Private Member Functions:
    // =========================

    // Generates the Clifford Algebra { \gamma_\mu } with signature (p,q)
    std::vector<GammaMatrix> generate_small_gammas_(void);
};



// Non-Member Functions:
// =====================

std::ostream& operator<<(std::ostream&, CliffordAlgebra const& A);

// Generates the Euclidean Clifford Algebra { \gamma_\mu } with signature (d=p+q,0)
std::vector<GammaMatrix> generate_euclidean_gammas(
    int const d,                                        // dimensionality of Clifford algebra = # of gamma's
    ModelParameters const pqn = ModelParameters(0,0,0)  // provides gamma5_prefactor for odd dimensions,
                                                        // set to zero by default
);


// Non-Member Utility Functions:
// =============================

GammaMatrix antisymmetrise(
    std::vector<GammaMatrix> const& gammas,
    std::vector<int> const& indices
);
