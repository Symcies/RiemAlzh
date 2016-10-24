/***************************************************************************************
*                                                                                      *
*                                     Deformetrica                                     *
*                                                                                      *
*    Copyright Inria and the University of Utah.  All rights reserved. This file is    *
*    distributed under the terms of the Inria Non-Commercial License Agreement.        *
*                                                                                      *
*                                                                                      *
****************************************************************************************/

#ifndef _LinearAlgebra_h
#define _LinearAlgebra_h

#include "Armadillo/ArmadilloLinearAlgebra.h"


/**
 *  \brief      Linear algebra class.
 *
 *  \copyright  Inria and the University of Utah
 *  \version    Deformetrica X.X
 *
 *  \details    The LinearAlgebra class only contains typedefs for vector, matrix, matrix lists, linear variable and
 *  			linear variable map types.
 */

template<class ScalarType>
class LinearAlgebra {

public :

	/// Vector type.
	typedef typename ArmadilloLinearAlgebra<ScalarType>::VectorType VectorType;
	/// Matrix type.
	typedef typename ArmadilloLinearAlgebra<ScalarType>::MatrixType MatrixType;
    

};


#endif /* _LinearAlgebra_h */