/**
    NCG++
    dirac_matrix.hpp

    Purpose:
    Provide a container/wrapper for the Dirac Matrix
    in terms of hermitian and anti-hermitian
    matrices.

    @author Peter Labus
    @version 0.1
    08.05.2017
*/

#pragma once

#include <vector>
#include "math_libraries.hpp"

template<typename FT>
class DiracMatrix
{
  public:

    FT operator()(int row, int column) const;
    FT operator()(int row, int column);


  private:

    const MathLibrary matrix_impl_;
    const ModelParameters pqn_;
    const CliffordAlgebra cliff_;
    const int num_H_;
    const int num_L_;
    std::vector< HermitianMatrix > H_;
    std::vector< AntiHermitianMatrix > L_;
};
