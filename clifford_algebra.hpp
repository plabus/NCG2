/**
    NCG++
    clifford_algebra.hpp

    Purpose:
    Given the model parameters (p, q) initialise
    the entire ODD gamma matrices

    \Gamma_\mu = { 1, \gamma_mu, \gamma_{\mu \nu}, ..., \gamma_5 }

    and provide a print / read function.


    @author Peter Labus
    @version 0.1
    09.05.2017
*/

#pragma once

#include <vector>


class CliffordAlgebra
{
  public:



  private:

    const ModelParameters pqn_;
    std::vector< GammaMatrix > gammas_;
};
