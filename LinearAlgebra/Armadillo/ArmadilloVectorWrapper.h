/***************************************************************************************
 *                                                                                      *
 *                                     Deformetrica                                     *
 *                                                                                      *
 *    Copyright Inria and the University of Utah. All rights reserved. This file is     *
 *    distributed under the terms of the Inria Non-Commercial License Agreement.        *
 *                                                                                      *
 *                                                                                      *
 ****************************************************************************************/

#ifndef _ArmadilloVectorWrapper_h
#define _ArmadilloVectorWrapper_h

/// Libraries files.
#include <vector>
#include <iostream>
#include <armadillo>
#include <cstdlib>


template<class ScalarType>
class ArmadilloMatrixWrapper;

template<class ScalarType>
class ArmadilloVectorWrapper;

template<class ScalarType>
ArmadilloVectorWrapper<ScalarType> operator*(ArmadilloMatrixWrapper<ScalarType> const &leftMatrix, ArmadilloVectorWrapper<ScalarType> const &rightVector);

template<class ScalarType>
ArmadilloVectorWrapper<ScalarType> operator-(const ArmadilloVectorWrapper<ScalarType> & left, const ArmadilloVectorWrapper<ScalarType> & right);

template<class ScalarType>
ArmadilloVectorWrapper<ScalarType> operator+(const ArmadilloVectorWrapper<ScalarType> & left, const ArmadilloVectorWrapper<ScalarType> & right);

template<class ScalarType>
ArmadilloVectorWrapper<ScalarType> operator+(const ArmadilloVectorWrapper<ScalarType> & leftVector, const ScalarType & rightScalar);

template<class ScalarType>
ArmadilloVectorWrapper<ScalarType> operator-(const ArmadilloVectorWrapper<ScalarType> & leftVector, const ScalarType & rightScalar);

template<class ScalarType>
ArmadilloVectorWrapper<ScalarType> operator/(const ArmadilloVectorWrapper<ScalarType> & leftVector, const ScalarType & rightScalar);

template<class ScalarType>
ArmadilloVectorWrapper<ScalarType> operator*(const ScalarType & leftScalar, const ArmadilloVectorWrapper<ScalarType> & rightVector);


template<class ScalarType>
std::ostream & operator<<(std::ostream& os, ArmadilloVectorWrapper<ScalarType> const& rhs);


// ADDED BY IGOR
template <class ScalarType>    
bool operator==(const ArmadilloVectorWrapper<ScalarType> & left, const ArmadilloVectorWrapper<ScalarType> & right);


/**
 *  \brief      Mathematical vector class (Armadillo).
 *
 *  \copyright  Inria and the University of Utah
 *  \version    Deformetrica X.X
 *
 *  \details    The ArmadilloVectorWrapper class is a wrapping for mathematical vector class, templated by type of element
 *              and using the Armadillo library.
 */
template<class ScalarType>
class ArmadilloVectorWrapper {

	////////////////////////////////////////////////////////////////////////////////////////////////////
	// Friend methods :
	////////////////////////////////////////////////////////////////////////////////////////////////////

	friend ArmadilloVectorWrapper<ScalarType> operator*<>(ArmadilloMatrixWrapper<ScalarType> const &leftMatrix, ArmadilloVectorWrapper<ScalarType> const &rightVector);

	friend ArmadilloVectorWrapper<ScalarType> operator-<>(const ArmadilloVectorWrapper<ScalarType> & left, const ArmadilloVectorWrapper<ScalarType> & right);
	friend ArmadilloVectorWrapper<ScalarType> operator+<>(const ArmadilloVectorWrapper<ScalarType> & left, const ArmadilloVectorWrapper<ScalarType> & right);

	friend ArmadilloVectorWrapper<ScalarType> operator+<>(const ArmadilloVectorWrapper<ScalarType> & leftVector, const ScalarType & rightScalar);
	friend ArmadilloVectorWrapper<ScalarType> operator-<>(const ArmadilloVectorWrapper<ScalarType> & leftVector, const ScalarType & rightScalar);
	friend ArmadilloVectorWrapper<ScalarType> operator/<>(const ArmadilloVectorWrapper<ScalarType> & leftVector, const ScalarType & rightScalar);

	friend ArmadilloVectorWrapper<ScalarType> operator*<>(const ScalarType & leftScalar, const ArmadilloVectorWrapper<ScalarType> & rightVector);
    // ADDED BY IGOR
    friend bool operator==<>(const ArmadilloVectorWrapper<ScalarType> & left, const ArmadilloVectorWrapper<ScalarType> & right);

	/// Overload of the stream operator.
	friend std::ostream & operator<<<>(std::ostream& os, ArmadilloVectorWrapper<ScalarType> const& rhs);


public :

	////////////////////////////////////////////////////////////////////////////////////////////////////
	// typedef :
	////////////////////////////////////////////////////////////////////////////////////////////////////

	/// Vector type
	typedef typename arma::Col<ScalarType> ArmadilloVectorType;
	/// Matrix type.
	typedef typename arma::Mat<ScalarType> ArmadilloMatrixType;

    /// Iterator on a raw Armadillo vector type.
    typedef typename ArmadilloVectorType::iterator iterator;
    typedef typename ArmadilloVectorType::const_iterator const_iterator;


	////////////////////////////////////////////////////////////////////////////////////////////////////
	// Constructors / Destructor :
	////////////////////////////////////////////////////////////////////////////////////////////////////

	/// Default constructor.
	inline ArmadilloVectorWrapper() : m_Vector() {}

 	/// Creates vector of \e len elements.
 	explicit inline ArmadilloVectorWrapper(unsigned len) : m_Vector(len) { }

 	/// Creates vector of \e len elements, all set to \e v0.
 	explicit inline ArmadilloVectorWrapper(unsigned len, ScalarType const& v0) : m_Vector(len) { m_Vector.fill(v0); }

	/// Constructor from a raw Armadillo vector (do not use it in Deformetrica!).
	inline ArmadilloVectorWrapper(const ArmadilloVectorType& v) { m_Vector = v; }

    /// Copy constructor.
    inline ArmadilloVectorWrapper(const ArmadilloVectorWrapper<ScalarType>& other) { m_Vector = other.m_Vector; }

	/// Access directly to the Armadillo vector (do not use it in Deformetrica!).
	inline ArmadilloVectorType const & toArmadillo() const { return m_Vector; }



	////////////////////////////////////////////////////////////////////////////////////////////////////
	// Methods :
	////////////////////////////////////////////////////////////////////////////////////////////////////

	/// Returns the length, number of elements, dimension of this vector.
	inline unsigned size() const { return m_Vector.n_rows; }
	/// Resizes to \e n elements.
	inline bool set_size(unsigned n) { m_Vector.set_size(n); return true; }

	/// Returns the number of elements.
	inline unsigned n_elem() const { return m_Vector.n_rows; }

	/// Gets value at element \e i.
	inline ScalarType get(unsigned int i) const { return m_Vector.at(i); }
	/// Gets contiguous values between \e first_index and \e last_index.
	inline ArmadilloVectorWrapper<ScalarType> subvec(const unsigned int first_index, const unsigned int last_index) const {
		return ArmadilloVectorWrapper<ScalarType>(m_Vector.subvec(first_index, last_index)); }

	/// Sets all elements of matrix to specified value.
	inline ArmadilloVectorWrapper<ScalarType>& fill(ScalarType const& scalar) { m_Vector.fill(scalar); return *this; }

	/// Converts a (column) vector to a matrix of size \e rows x \e cols, in a column-wise fashion.
	inline ArmadilloMatrixWrapper<ScalarType> convert_to_matrix_col_wise(unsigned rows, unsigned cols) const {
		ArmadilloMatrixType result(m_Vector);
		result.reshape(rows, cols);
		return ArmadilloMatrixWrapper<ScalarType>(result); }
	/// Converts a (column) vector to a matrix of size \e rows x \e cols, in row-wise fashion.
	inline ArmadilloMatrixWrapper<ScalarType> convert_to_matrix_row_wise(unsigned rows, unsigned cols) const {
		ArmadilloMatrixType result(m_Vector);
		result.reshape(cols, rows);
		return ArmadilloMatrixWrapper<ScalarType>(arma::trans(result)); }

	/// Concatenates vertically.
	inline ArmadilloVectorWrapper<ScalarType> join_cols(const ArmadilloVectorWrapper<ScalarType>& v) const {
		return ArmadilloVectorWrapper<ScalarType>(arma::join_cols(m_Vector, v.toArmadillo())); }

    /// Returns an iterator on the first element of the raw Armadillo vector.
    inline       iterator begin()       { return m_Vector.begin(); }
    inline const_iterator begin() const { return m_Vector.begin(); }

    /// Returns an iterator on the end of the raw Armadillo vector.
    inline       iterator end()       { return m_Vector.end(); }
    inline const_iterator end() const { return m_Vector.end(); }
    
    /// Probably only for debugging : should be deleted at the end if possible
   // TODO : Address the previous concern and delete if possible
   inline std::vector<double> convert_to_stl() { 
       return arma::conv_to<std::vector<double>>::from(m_Vector); }


	////////////////////////////////////////////////////////////////////////////////////////////////////
	// Arithmetic operations :
	////////////////////////////////////////////////////////////////////////////////////////////////////

	/// Returns magnitude (length) of vector.
	inline ScalarType magnitude() const { return arma::norm(m_Vector, 2); }

	/// Returns sum of squares of elements.
	inline ScalarType squared_magnitude() const { ScalarType result = arma::norm(m_Vector, 2); return result*result; }

	/// Returns sum of elements.
	inline ScalarType sum() const { return arma::sum(m_Vector); }
	/// Returns sum of squares of elements.
	inline ScalarType sum_of_squares() const { return arma::sum(arma::square(m_Vector)); }

	/// Smallest value.
	inline ScalarType min_value() const { return m_Vector.min(); }
	/// Largest value.
	inline ScalarType max_value() const { return m_Vector.max(); }

	/// Smallest value index.
//	inline unsigned int index_min() const { return m_Vector.index_min(); } // REQUIRES ARMADILLO 7.400
	inline unsigned int index_min() const {
		unsigned int index = 0;
		for (unsigned int k = 1 ; k < m_Vector.size() ; ++k) {
			if (m_Vector(k) < m_Vector(index)) { index = k; } }
		return index; }
	/// Largest value index.
//	inline unsigned int index_max() const { return m_Vector.index_max(); } // REQUIRES ARMADILLO 7.400
	inline unsigned int index_max() const {
		unsigned int index = 0;
		for (unsigned int k = 1 ; k < m_Vector.size() ; ++k) {
			if (m_Vector(k) > m_Vector(index)) { index = k; } }
		return index; }

	/// Access the contiguous block storing the elements in the vector.
	inline ScalarType const* memptr() const { return m_Vector.memptr(); }

	/// Access the contiguous block storing the elements in the vector.
	inline ScalarType      * memptr() { return m_Vector.memptr(); }



	////////////////////////////////////////////////////////////////////////////////////////////////////
	// Operator overloading :
	////////////////////////////////////////////////////////////////////////////////////////////////////

    /// Returns reference to the element at index \e i for reading.
    inline ScalarType const & operator()(unsigned int i) const { return m_Vector(i); }
	/// Returns reference to the element at index \e i for reading or writing.
	inline ScalarType       & operator()(unsigned int i)       { return m_Vector(i); }

	/// Returns reference to the element at index \e i.
	/// \warning Please use operator() in Deformetrica instead of operator[].
	inline ScalarType       & operator[](unsigned int i)       { return m_Vector(i); }
	/// Returns reference to the element at index \e i.
	/// \warning Please use operator() in Deformetrica instead of operator[].
	inline ScalarType const & operator[](unsigned int i) const { return m_Vector(i); }

	/// Scalar multiplication of lhs vector in situ by \e rhsScalar.
	inline ArmadilloVectorWrapper<ScalarType>& operator*=(ScalarType rhsScalar) { m_Vector *= rhsScalar; return *this; }
	/// Scalar division of lhs vector in situ by \e rhsScalar.
	inline ArmadilloVectorWrapper<ScalarType>& operator/=(ScalarType rhsScalar) { m_Vector /= rhsScalar; return *this; }

	/// Adds \e rhs to lhs vector in situ.
	inline ArmadilloVectorWrapper<ScalarType>& operator+=(ArmadilloVectorWrapper<ScalarType> const& rhs) { m_Vector += rhs.m_Vector; return *this; }

	/// Unary minus operator.
	inline ArmadilloVectorWrapper<ScalarType> operator-() const { return ArmadilloVectorWrapper<ScalarType>(- m_Vector); }
	/// Unary plus operator.
	inline ArmadilloVectorWrapper<ScalarType> operator+() const { return ArmadilloVectorWrapper<ScalarType>(+ m_Vector); }

	/// Scalar multiplication of lhs vector by \e scalar.
	inline ArmadilloVectorWrapper<ScalarType> operator*(ScalarType scalar) const { return ArmadilloVectorWrapper<ScalarType>(m_Vector*scalar); }
	/// Scalar division of lhs vector by \e scalar.
	inline ArmadilloVectorWrapper<ScalarType> operator/(ScalarType scalar) const { return ArmadilloVectorWrapper<ScalarType>(m_Vector/scalar); }
    

private :

	ArmadilloVectorType m_Vector;

};


template<class ScalarType>
inline ScalarType dot_product (ArmadilloVectorWrapper<ScalarType> const & v1, ArmadilloVectorWrapper<ScalarType> const & v2) {
	return arma::dot(v1.toArmadillo(), v2.toArmadillo()); }

template<class ScalarType>
inline ArmadilloVectorWrapper<ScalarType> cross_3d (ArmadilloVectorWrapper<ScalarType> const & v1, ArmadilloVectorWrapper<ScalarType> const & v2) {
	return ArmadilloVectorWrapper<ScalarType>(arma::cross(v1.toArmadillo(), v2.toArmadillo())); }


#include "ArmadilloVectorWrapper_friend.h"



#endif /* _ArmadilloVectorWrapper_h */