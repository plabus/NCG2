/**
    NCG++
    gamma_matrix.hpp

    Purpose:
    Matrix class especially for
    gamma matrices

    @author Peter Labus
    @version 0.1
    10.05.2017
*/

#pragma once

#include <iostream>
#include <vector>
#include <complex>



//===================//
//                   //
//    GammaMatrix    //
//                   //
//===================//

class GammaMatrix
{
  public:

    // Public Member Functions:
    // ========================

    // Constructors:
    GammaMatrix(void) = delete;
    explicit GammaMatrix(const int size);
    explicit GammaMatrix(std::initializer_list< std::complex<int> > const& list);
    GammaMatrix& operator=(GammaMatrix const& other);

    // Assigment Operators:
    GammaMatrix& operator*=(GammaMatrix const& other);
    GammaMatrix& operator*=(std::complex<int> c);

    // Logical Operators:
    bool operator==(GammaMatrix const& other) const;

    // Miscellaneous:
    GammaMatrix operator-(void) const; //!< Return the negative of this as a new GammaMatrix

    // Getters:
    std::complex<int> operator()(const int row, const int col) const { return M_[row*size_+col]; }
    int size(void) const { return size_; }

    // Setters:
    std::complex<int>& operator()(const int row, const int col) { return M_[row*size_+col]; }


  private:

    // Private Member Variables:
    // =========================

    const int size_;
    std::vector< std::complex<int>  > M_;
};


// Non-Member Functions:
// =====================

std::ostream& operator<<(std::ostream& os, GammaMatrix const& A);

// These are addition, substraction, multiplication and outer (tensor) multiplication
GammaMatrix operator+(GammaMatrix const& A, GammaMatrix const& B);
GammaMatrix operator-(GammaMatrix const& A, GammaMatrix const& B);
GammaMatrix operator*(GammaMatrix const& A, GammaMatrix const& B);
GammaMatrix operator%(GammaMatrix const& A, GammaMatrix const& B);

// These allow to multiply a GammaMatrix with a complex integer
// and divide by an integer (integer divison)
GammaMatrix operator*(std::complex<int> c, GammaMatrix const& A);
GammaMatrix operator/(GammaMatrix const& A, int c);

// These are commutators and anticommutators of two GammaMatrices
GammaMatrix commutator(GammaMatrix const& A, GammaMatrix const& B);
GammaMatrix anticommutator(GammaMatrix const& A, GammaMatrix const& B);

// These return true is M is hermitain / antihermitian respectively
bool is_hermitian(GammaMatrix const& M);
bool is_antihermitian(GammaMatrix const& M);



//=====================//
//                     //
//    PauliMatrices    //
//                     //
//=====================//

struct PauliMatrices
{
  PauliMatrices(void);
  const GammaMatrix sigma1;
  const GammaMatrix sigma2;
  const GammaMatrix sigma3;
};

// Generates a d x d unity (gamma) matrix
GammaMatrix Unity(const int d);
