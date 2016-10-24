/***************************************************************************************
 *                                                                                      *
 *                                     Deformetrica                                     *
 *                                                                                      *
 *    Copyright Inria and the University of Utah.  All rights reserved. This file is    *
 *    distributed under the terms of the Inria Non-Commercial License Agreement.        *
 *                                                                                      *
 *                                                                                      *
 ****************************************************************************************/


#ifndef _ArmadilloMatrixWrapper_friend_h
#define _ArmadilloMatrixWrapper_friend_h


#ifndef _ArmadilloMatrixWrapper_h
#error Do not include ArmadilloMatrixWrapper_friend.h : include ArmadilloMatrixWrapper.h instead
#endif

#include "ArmadilloVectorWrapper.h"
#include <iostream>




//
// Left Matrix / Right Matrix :
//
template<class ScalarType>
inline ArmadilloMatrixWrapper<ScalarType> operator-(const ArmadilloMatrixWrapper<ScalarType> & left, const ArmadilloMatrixWrapper<ScalarType> & right) {
	return ArmadilloMatrixWrapper<ScalarType>(left.m_Matrix - right.m_Matrix);
}

template<class ScalarType>
inline ArmadilloMatrixWrapper<ScalarType> operator+(const ArmadilloMatrixWrapper<ScalarType> & left, const ArmadilloMatrixWrapper<ScalarType> & right) {
	return ArmadilloMatrixWrapper<ScalarType>(left.m_Matrix + right.m_Matrix);
}

template<class ScalarType>
inline ArmadilloMatrixWrapper<ScalarType> operator*(const ArmadilloMatrixWrapper<ScalarType> & left, const ArmadilloMatrixWrapper<ScalarType> & right) {
	return ArmadilloMatrixWrapper<ScalarType>(left.m_Matrix * right.m_Matrix);
}



//
// (Left-Right) Matrix / Scalar :
//
template<class ScalarType>
inline ArmadilloMatrixWrapper<ScalarType> operator*(const ArmadilloMatrixWrapper<ScalarType> & leftMatrix, ScalarType const& rightScalar) {
	return ArmadilloMatrixWrapper<ScalarType>(leftMatrix.m_Matrix * rightScalar);
}

template<class ScalarType>
inline ArmadilloMatrixWrapper<ScalarType> operator*(ScalarType const& leftScalar, const ArmadilloMatrixWrapper<ScalarType> & rightMatrix) {
	return ArmadilloMatrixWrapper<ScalarType>(leftScalar * rightMatrix.m_Matrix);
}


template<class ScalarType>
inline ArmadilloMatrixWrapper<ScalarType> operator/(const ArmadilloMatrixWrapper<ScalarType> & leftMatrix, ScalarType const& rightScalar) {
	return ArmadilloMatrixWrapper<ScalarType>(leftMatrix.m_Matrix / rightScalar);
}

template<class ScalarType>
inline ArmadilloMatrixWrapper<ScalarType> operator/(ScalarType const& leftScalar, const ArmadilloMatrixWrapper<ScalarType> & rightMatrix) {
	return ArmadilloMatrixWrapper<ScalarType>(leftScalar / rightMatrix.m_Matrix);
}


//
// Left Matrix / Right Vector :
//
template<class ScalarType>
ArmadilloVectorWrapper<ScalarType> operator*(ArmadilloMatrixWrapper<ScalarType> const &leftMatrix, ArmadilloVectorWrapper<ScalarType> const &rightVector) {
	return ArmadilloVectorWrapper<ScalarType>(leftMatrix.m_Matrix * rightVector.m_Vector);
}



//
// Other operations :
//

template <class ScalarType>
ArmadilloMatrixWrapper<ScalarType> diagonal_matrix(unsigned N, ScalarType const & value) {
	typename ArmadilloMatrixWrapper<ScalarType>::ArmadilloMatrixType result;
	return ArmadilloMatrixWrapper<ScalarType>( result.eye(N, N)*value );
}

template <class ScalarType>
ScalarType trace(ArmadilloMatrixWrapper<ScalarType> const& M) {
	return arma::trace(M.m_Matrix);
}

template <class ScalarType>
ArmadilloMatrixWrapper<ScalarType> chol(ArmadilloMatrixWrapper<ScalarType> const& M) {
	return ArmadilloMatrixWrapper<ScalarType>(arma::chol(M.m_Matrix));
}

template <class ScalarType>
ArmadilloMatrixWrapper<ScalarType> inverse(ArmadilloMatrixWrapper<ScalarType> const& M) {
	return ArmadilloMatrixWrapper<ScalarType>(arma::inv(M.m_Matrix));
}

template <class ScalarType>
ArmadilloMatrixWrapper<ScalarType> inverse_sympd(ArmadilloMatrixWrapper<ScalarType> const& M) {
	return ArmadilloMatrixWrapper<ScalarType>(arma::inv_sympd(M.m_Matrix));
}

template <class ScalarType>
ArmadilloVectorWrapper<ScalarType> eigenvalues_sym(ArmadilloMatrixWrapper<ScalarType> const& M) {
	//	arma::Col<ScalarType> eigen_values;
	typename ArmadilloVectorWrapper<ScalarType>::ArmadilloVectorType eigen_values;
	arma::eig_sym(eigen_values, M.m_Matrix);
	return ArmadilloVectorWrapper<ScalarType>(eigen_values);
}

template <class ScalarType>
std::ostream & operator<<(std::ostream& os, ArmadilloMatrixWrapper<ScalarType> const& rhs) {
	return (os << rhs.m_Matrix);
}




#endif /* _ArmadilloMatrixWrapper_friend_h */