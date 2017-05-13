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
#include "model_parameters.hpp"
#include "clifford_algebra.hpp"

template<typename FT>
class DiracMatrix
{
  public:

    DiracMatrix(void) = delete;
    DiracMatrix(
        ModelParameters const pqn,
        MathLibrary const matrix_impl = MathLibrary::Custom
    )
      :
    pqn_(pqn)
    , matrix_impl_(matrix_impl)
    , cliff_(CliffordAlgebra(pqn_))
    , H_( std::vector< HermitianMatrix<FT> >(0) )
    , L_( std::vector< AntiHermitianMatrix<FT> >(0) )
    {
    }

    FT operator()(int row, int column) const;
    FT operator()(int row, int column);

    int num_H(void) const { return H_.size(); }
    int num_L(void) const { return L_.size(); }
    int size(void) const { return H_.size() + L_.size(); }

  private:

    const ModelParameters pqn_;
    const MathLibrary matrix_impl_;
    const CliffordAlgebra cliff_;
    std::vector< HermitianMatrix<FT> > H_;
    std::vector< AntiHermitianMatrix<FT> > L_;
};
