/***************************************************************************************
 *                                                                                      *
 *                                     Deformetrica                                     *
 *                                                                                      *
 *    Copyright Inria and the University of Utah.  All rights reserved. This file is    *
 *    distributed under the terms of the Inria Non-Commercial License Agreement.        *
 *                                                                                      *
 *                                                                                      *
 ****************************************************************************************/


#ifndef _ArmadilloVectorWrapper_friend_h
#define _ArmadilloVectorWrapper_friend_h

#ifndef _ArmadilloVectorWrapper_h
#error Do not include ArmadilloVectorWrapper_friend.h : include ArmadilloVectorWrapper.h instead
#endif

#include "ArmadilloVectorWrapper.h"
#include <iostream>


//
// Left Vector / Right Vector :
//
template<class ScalarType>
inline ArmadilloVectorWrapper<ScalarType> operator-(const ArmadilloVectorWrapper<ScalarType> & left, const ArmadilloVectorWrapper<ScalarType> & right) {
	return ArmadilloVectorWrapper<ScalarType>(left.m_Vector - right.m_Vector);
}

template<class ScalarType>
inline ArmadilloVectorWrapper<ScalarType> operator+(const ArmadilloVectorWrapper<ScalarType> & left, const ArmadilloVectorWrapper<ScalarType> & right) {
	return ArmadilloVectorWrapper<ScalarType>(left.m_Vector + right.m_Vector);
}



//
// Left Vector / Right Scalar :
//
template<class ScalarType>
inline ArmadilloVectorWrapper<ScalarType> operator+(const ArmadilloVectorWrapper<ScalarType> & leftVector, const ScalarType & rightScalar) {
	return ArmadilloVectorWrapper<ScalarType>(leftVector.m_Vector + rightScalar);
}

template<class ScalarType>
inline ArmadilloVectorWrapper<ScalarType> operator-(const ArmadilloVectorWrapper<ScalarType> & leftVector, const ScalarType & rightScalar) {
	return ArmadilloVectorWrapper<ScalarType>(leftVector.m_Vector - rightScalar);
}

template<class ScalarType>
inline ArmadilloVectorWrapper<ScalarType> operator/(const ArmadilloVectorWrapper<ScalarType> & leftVector, const ScalarType & rightScalar) {
	return ArmadilloVectorWrapper<ScalarType>(leftVector.m_Vector / rightScalar);
}



//
// Left Scalar / Right Vector :
//
template<class ScalarType>
inline ArmadilloVectorWrapper<ScalarType> operator*(const ScalarType & leftScalar, const ArmadilloVectorWrapper<ScalarType> & rightVector) {
	return ArmadilloVectorWrapper<ScalarType>(leftScalar * rightVector.m_Vector);
}

//
// Other operations :
//
template <class ScalarType>
std::ostream & operator<<(std::ostream& os, ArmadilloVectorWrapper<ScalarType> const& rhs) {
	return (os << rhs.m_Vector.t());
}


/// ADDED BY IGOR:
template <class ScalarType>
bool operator==(const ArmadilloVectorWrapper<ScalarType> & left, const ArmadilloVectorWrapper<ScalarType> & right) {
    arma::umat a =  (left.m_Vector == right.m_Vector);
    for(auto it = a.begin(); it != a.end(); ++it)
    {
        if(*it == 0) { return false; }
    }
    return true;
}



#endif /* _ArmadilloVectorWrapper_friend_h */