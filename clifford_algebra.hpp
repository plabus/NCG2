/**
    NCG++
    clifford_algebra.hpp

    Purpose:
    Given the model parameters (p, q) initialise
    the entire ODD gamma matrices

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



class CliffordAlgebra
{
  public:

    CliffordAlgebra(void) = delete;
    explicit CliffordAlgebra(ModelParameters const pqn);

    friend std::ostream& operator<<(std::ostream&, CliffordAlgebra const& A);


  private:

    const ModelParameters pqn_;
    std::vector<GammaMatrix> Gammas_;

    std::vector<GammaMatrix> generate_gammas_(void);
};
