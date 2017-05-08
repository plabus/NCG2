/**
    NCG++
    dirac_matrix.hpp

    Purpose:
    Provide a container/wrapper for the Dirac
    Matrix in terms of hermitian and anti-hermitian
    matrices.

    @author Peter Labus
    @version 0.1
    08.05.2017
*/

#pragma once

template<typename FT>
class DiracMatrix
{
  public:

    virtual FT operator()(int row, int column) const;
    virtual FT operator()(int row, int column);
};
