/***************************************************************************************
*                                                                                      *
*                                     Deformetrica                                     *
*                                                                                      *
*    Copyright Inria and the University of Utah.  All rights reserved. This file is    *
*    distributed under the terms of the Inria Non-Commercial License Agreement.        *
*                                                                                      *
*                                                                                      *
****************************************************************************************/

#ifndef _ArmadilloLinearAlgebra_h
#define _ArmadilloLinearAlgebra_h

#include "ArmadilloVectorWrapper.h"
#include "ArmadilloMatrixWrapper.h"


/**
 *  \brief      Armadillo linear algebra class.
 *
 *  \copyright  Inria and the University of Utah
 *  \version    Deformetrica X.X
 *
 *  \details    The LinearAlgebra class only contains typedefs for armadillo vector and matrix types.
 */

template<class ScalarType>
class ArmadilloLinearAlgebra {

public :

	/// Vector type.
	typedef ArmadilloVectorWrapper<ScalarType> VectorType;
	/// Matrix type.
	typedef ArmadilloMatrixWrapper<ScalarType> MatrixType;

};


#endif /* _ArmadilloLinearAlgebra_h */