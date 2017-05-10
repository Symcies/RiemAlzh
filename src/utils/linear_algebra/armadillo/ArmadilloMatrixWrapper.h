/***************************************************************************************
 *                                                                                      *
 *                                     Deformetrica                                     *
 *                                                                                      *
 *    Copyright Inria and the University of Utah.  All rights reserved. This file is    *
 *    distributed under the terms of the Inria Non-Commercial License Agreement.        *
 *                                                                                      *
 *                                                                                      *
 ****************************************************************************************/

#ifndef _ArmadilloMatrixWrapper_h
#define _ArmadilloMatrixWrapper_h

#include <armadillo>


template<class TScalar>
class ArmadilloMatrixWrapper;

template<class TScalar>
class ArmadilloVectorWrapper;



template<class TScalar>
ArmadilloMatrixWrapper<TScalar> operator-(const ArmadilloMatrixWrapper<TScalar> & left, const ArmadilloMatrixWrapper<TScalar> & right);
template<class TScalar>
ArmadilloMatrixWrapper<TScalar> operator+(const ArmadilloMatrixWrapper<TScalar> & left, const ArmadilloMatrixWrapper<TScalar> & right);
template<class TScalar>
ArmadilloMatrixWrapper<TScalar> operator*(const ArmadilloMatrixWrapper<TScalar> & left, const ArmadilloMatrixWrapper<TScalar> & right);

template<class TScalar>
ArmadilloMatrixWrapper<TScalar> operator*(const ArmadilloMatrixWrapper<TScalar> & leftMatrix, TScalar const& rightScalar);
template<class TScalar>
ArmadilloMatrixWrapper<TScalar> operator*(TScalar const& leftScalar, const ArmadilloMatrixWrapper<TScalar> & rightMatrix);

template<class TScalar>
ArmadilloMatrixWrapper<TScalar> operator/(const ArmadilloMatrixWrapper<TScalar> & leftMatrix, unsigned int const& rightScalar);
template<class TScalar>
ArmadilloMatrixWrapper<TScalar> operator/(TScalar const& leftScalar, const ArmadilloMatrixWrapper<TScalar> & rightMatrix);

template<class TScalar>
ArmadilloVectorWrapper<TScalar> operator*(ArmadilloMatrixWrapper<TScalar> const &leftMatrix, ArmadilloVectorWrapper<TScalar> const &rightVector);



template<class TScalar>
ArmadilloMatrixWrapper<TScalar> diagonal_matrix(unsigned N, TScalar const & value);

template<class TScalar>
ArmadilloMatrixWrapper<TScalar> identity_matrix(unsigned N);

template<class TScalar>
TScalar trace(ArmadilloMatrixWrapper<TScalar> const& M);
template<class TScalar>
ArmadilloMatrixWrapper<TScalar> inverse(ArmadilloMatrixWrapper<TScalar> const& M);
template<class TScalar>
ArmadilloMatrixWrapper<TScalar> inverse_sympd(ArmadilloMatrixWrapper<TScalar> const& M);
template<class TScalar>
ArmadilloVectorWrapper<TScalar> eigenvalues_sym(ArmadilloMatrixWrapper<TScalar> const& M);

template <class TScalar>
std::ostream & operator<<(std::ostream& os, ArmadilloMatrixWrapper<TScalar> const& rhs);



/**
 *  \brief      Mathematical matrix class (Armadillo).
 *
 *  \copyright  Inria and the University of Utah
 *  \version    Deformetrica X.X
 *
 *  \details    The ArmadilloMatrixWrapper class is a wrapping for mathematical matrix class, templated by type of element
 *              and using Armadillo library.
 *  \warning    Elements are stored with column-major ordering (ie. column by column).
 *
 */
template<class TScalar>
class ArmadilloMatrixWrapper {

	////////////////////////////////////////////////////////////////////////////////////////////////////
	// Friend methods :
	////////////////////////////////////////////////////////////////////////////////////////////////////

	friend ArmadilloMatrixWrapper<TScalar> operator-<>(const ArmadilloMatrixWrapper<TScalar> & left, const ArmadilloMatrixWrapper<TScalar> & right);
	friend ArmadilloMatrixWrapper<TScalar> operator+<>(const ArmadilloMatrixWrapper<TScalar> & left, const ArmadilloMatrixWrapper<TScalar> & right);
	friend ArmadilloMatrixWrapper<TScalar> operator*<>(const ArmadilloMatrixWrapper<TScalar> & left, const ArmadilloMatrixWrapper<TScalar> & right);

	friend ArmadilloMatrixWrapper<TScalar> operator*<>(const ArmadilloMatrixWrapper<TScalar> & leftMatrix, TScalar const& rightScalar);
	friend ArmadilloMatrixWrapper<TScalar> operator*<>(TScalar const& leftScalar, const ArmadilloMatrixWrapper<TScalar> & rightMatrix);

	friend ArmadilloMatrixWrapper<TScalar> operator/<>(const ArmadilloMatrixWrapper<TScalar> & leftMatrix, unsigned int const& rightScalar);
	friend ArmadilloMatrixWrapper<TScalar> operator/<>(TScalar const& leftScalar, const ArmadilloMatrixWrapper<TScalar> & rightMatrix);

	friend ArmadilloVectorWrapper<TScalar> operator*<>(ArmadilloMatrixWrapper<TScalar> const &leftMatrix, ArmadilloVectorWrapper<TScalar> const &rightVector);

	/// Creates a diagonal matrix size \e N x \e N with with all diagonal elements set to \e value.
	friend ArmadilloMatrixWrapper<TScalar> diagonal_matrix<>(unsigned N, TScalar const & value);
    
    /// Creates a identity matrix size \e N x \e N with with all diagonal elements set to 1.
	friend ArmadilloMatrixWrapper<TScalar> identity_matrix<>(unsigned N);
    
	/// Sum of the diagonal elements of square matrix \e M.
	friend TScalar trace<>(ArmadilloMatrixWrapper<TScalar> const& M);

	/// Inverse of square matrix \e M.
	/// \warning A std::logic_error exception is thrown if \e M does not evaluate to a square matrix.
	friend ArmadilloMatrixWrapper<TScalar> inverse<>(ArmadilloMatrixWrapper<TScalar> const& M);

	/// Inverse of symmetric positive definite matrix \e M.
	/// \warning There is currently (i.e. 04-2015) no explicit check whether \e M is symmetric or positive definite.
	friend ArmadilloMatrixWrapper<TScalar> inverse_sympd<>(ArmadilloMatrixWrapper<TScalar> const& M);

	/// Inverse of symmetric positive definite matrix \e M.
	/// \warning There is currently (i.e. 04-2015) no explicit check whether \e M is symmetric.
	friend ArmadilloVectorWrapper<TScalar> eigenvalues_sym<>(ArmadilloMatrixWrapper<TScalar> const& M);

	/// Overload of the stream operator.
	friend std::ostream & operator<<<>(std::ostream& os, ArmadilloMatrixWrapper<TScalar> const& rhs);



public :

	////////////////////////////////////////////////////////////////////////////////////////////////////
	// typedef :
	////////////////////////////////////////////////////////////////////////////////////////////////////

	/// Vector type
	typedef typename arma::Col<TScalar> ArmadilloVectorType;
	/// Matrix type.
	typedef typename arma::Mat<TScalar> ArmadilloMatrixType;



	////////////////////////////////////////////////////////////////////////////////////////////////////
	// Constructors / Destructor :
	////////////////////////////////////////////////////////////////////////////////////////////////////

	/// Default constructor.
	inline ArmadilloMatrixWrapper() : m_Matrix() {}

	/// Matrix constructor of size \e r rows by \e c columns.
	explicit inline ArmadilloMatrixWrapper(unsigned r, unsigned c) : m_Matrix(r, c) {}

	/// Construct a matrix of size \e r rows by \e c columns, and all elements equal to \e v0.
	explicit inline ArmadilloMatrixWrapper(unsigned r, unsigned c, TScalar const& v0) : m_Matrix(r, c) { m_Matrix.fill(v0); }

	/// Construct a matrix of size \e r rows by \e c columns, initialized by a (column-wise) memory block.
	/// \warning The auxiliary memory \e data_block is not copied !
	explicit inline ArmadilloMatrixWrapper(TScalar *data_block, unsigned r, unsigned c) : m_Matrix(data_block, r, c, false, true) {}

	/// Constructor which converts a (column) vector to a matrix.
	inline ArmadilloMatrixWrapper(const ArmadilloVectorWrapper<TScalar> & v) : m_Matrix(v.toArmadillo()) {}

	/// Special constructor (do not use it in Deformetrica!).
	inline ArmadilloMatrixWrapper(const ArmadilloMatrixType & M) : m_Matrix(M) {}

//	/// Access directly to the Armadillo matrix (do not use it in Deformetrica!).
//	inline ArmadilloMatrixType const & toArmadillo() const { return m_Matrix; }



	////////////////////////////////////////////////////////////////////////////////////////////////////
	// Methods :
	////////////////////////////////////////////////////////////////////////////////////////////////////

	/// Returns the number of rows.
	inline unsigned rows()    const { return m_Matrix.n_rows; }

	/// Returns the number of columns.
	inline unsigned columns()  const { return m_Matrix.n_cols; }

	/// Returns the number of columns.
	inline unsigned cols()    const { return m_Matrix.n_cols; }

	/// Returns the number of elements (This equals to rows() * cols()).
	inline unsigned size()    const { return m_Matrix.n_elem; }

	/// Access the contiguous block storing the elements in the matrix row-wise.
	inline TScalar      * memory_pointer_row_wise() { ArmadilloMatrixType result = arma::vectorise(m_Matrix, 1); return result.memptr(); }

	/// Access the contiguous block storing the elements in the matrix row-wise.
	inline TScalar const* memory_pointer_row_wise() const { ArmadilloMatrixType result = arma::vectorise(m_Matrix, 1); return result.memptr(); }

	/// Resize to r rows by c columns. Old data lost.
	inline bool set_size(unsigned r, unsigned c) { m_Matrix.set_size(r, c); return true; }

	/// Set all elements of matrix to specified value.
	inline ArmadilloMatrixWrapper<TScalar>& fill(TScalar const& scalar) { m_Matrix.fill(scalar); return *this; }

	/// Gets \e n rows beginning at \e rowstart.
	inline ArmadilloMatrixWrapper<TScalar> get_n_rows(unsigned rowstart, unsigned n) const {
		return ArmadilloMatrixWrapper<TScalar>(m_Matrix.rows(rowstart, rowstart+n-1));
	}

	/// Gets \e n columns beginning at \e colstart.
	inline ArmadilloMatrixWrapper<TScalar> get_n_columns(unsigned colstart, unsigned n) const {
		return ArmadilloMatrixWrapper<TScalar>(m_Matrix.cols(colstart, colstart+n-1));
	}

	/// Sets columns to those in \e M, starting at \e starting_column, then return *this.
	inline ArmadilloMatrixWrapper<TScalar> & set_columns(unsigned starting_column, ArmadilloMatrixWrapper<TScalar> const& M) {
		m_Matrix(arma::span(0, 0+M.rows()-1), arma::span(starting_column, starting_column+M.cols()-1)) = M.m_Matrix;
		return *this;
	}

	/// Sets values of this matrix to those of \e M, starting at (\e top, \e left).
	inline ArmadilloMatrixWrapper<TScalar>& update(ArmadilloMatrixWrapper<TScalar> const& M, unsigned top=0, unsigned left=0) {
		m_Matrix(arma::span(top, top+M.rows()-1), arma::span(left, left+M.cols()-1)) = M.m_Matrix;
		return *this;
	}

	/// Gets a vector equal to the given row.
	inline ArmadilloVectorWrapper<TScalar> get_row(unsigned r) const { return ArmadilloVectorWrapper<TScalar>(m_Matrix.row(r).t()); }

	/// Gets a vector equal to the given column.
	inline ArmadilloVectorWrapper<TScalar> get_column(unsigned c) const { return ArmadilloVectorWrapper<TScalar>(m_Matrix.col(c)); }

	/// Sets the \e i-th row to \e v.
	inline void set_row(unsigned i, ArmadilloVectorWrapper<TScalar> const& v) { m_Matrix.row(i) = v.toArmadillo().t(); }

	/// Sets the elements of the i'th row to \e value, then return *this.
	inline ArmadilloMatrixWrapper<TScalar> & set_row(unsigned i, TScalar value ) { m_Matrix.row(i).fill(value); return *this; }

	/// Sets \e j-th column to \e v.
	inline void set_column(unsigned j, ArmadilloVectorWrapper<TScalar> const& v) { m_Matrix.col(j) = v.toArmadillo(); }

	/// Sets this matrix to an identity matrix.
	inline void set_identity() { m_Matrix.eye(); }

	/// Converts a matrix to vector with a row-wise order.
	inline ArmadilloVectorWrapper<TScalar> vectorise_row_wise() const { return ArmadilloVectorWrapper<TScalar>(arma::vectorise(m_Matrix, 1)); }



	////////////////////////////////////////////////////////////////////////////////////////////////////
	// Arithmetic operations :
	////////////////////////////////////////////////////////////////////////////////////////////////////

	/// Returns Frobenius norm of matrix (sqrt of sum of squares of its elements).
	inline TScalar frobenius_norm() const { return arma::norm(m_Matrix, "fro"); }

	/// Returns transpose.
	inline ArmadilloMatrixWrapper<TScalar> transpose() const { return ArmadilloMatrixWrapper<TScalar>(m_Matrix.t()); }



	////////////////////////////////////////////////////////////////////////////////////////////////////
	// Operator overloading :
	////////////////////////////////////////////////////////////////////////////////////////////////////

	/// Access an element for reading.
	inline TScalar const & operator()(unsigned r, unsigned c) const { return m_Matrix(r, c); }
	/// Access an element for reading or writing.
	inline TScalar       & operator()(unsigned r, unsigned c)       { return m_Matrix(r, c); }

	/// Scalar division of lhs matrix in situ by \e rhsScalar.
	inline ArmadilloMatrixWrapper<TScalar>& operator/=(TScalar rhsScalar) { m_Matrix /= rhsScalar; return *this; }
	/// Scalar multiplication in situ of lhs matrix by \e rhsScalar.
	inline ArmadilloMatrixWrapper<TScalar>& operator*=(TScalar rhsScalar) { m_Matrix *= rhsScalar; return *this; }
	/// Add rhs to lhs matrix in situ.
	inline ArmadilloMatrixWrapper<TScalar>& operator+=(TScalar rhsScalar) { m_Matrix += rhsScalar; return *this; }
	/// Subtract rhs from lhs matrix in situ.
	inline ArmadilloMatrixWrapper<TScalar>& operator-=(TScalar rhsScalar) { m_Matrix -= rhsScalar; return *this; }

	/// Multiply lhs matrix in situ by rhs.
	inline ArmadilloMatrixWrapper<TScalar>& operator*=(ArmadilloMatrixWrapper<TScalar> const& rhsMatrix)
							{ m_Matrix *= rhsMatrix.m_Matrix; return *this; }
	/// Add rhs to lhs matrix in situ.
	inline ArmadilloMatrixWrapper<TScalar>& operator+=(ArmadilloMatrixWrapper<TScalar> const& rhsMatrix)
							{ m_Matrix += rhsMatrix.m_Matrix; return *this; }
	/// Subtract rhs from lhs matrix in situ.
	inline ArmadilloMatrixWrapper<TScalar>& operator-=(ArmadilloMatrixWrapper<TScalar> const& rhsMatrix)
							{ m_Matrix -= rhsMatrix.m_Matrix; return *this; }

	/// Unary minus operator.
	inline ArmadilloMatrixWrapper<TScalar> operator-() const { return ArmadilloMatrixWrapper<TScalar>(- m_Matrix); }
	/// Unary plus operator.
	inline ArmadilloMatrixWrapper<TScalar> operator+() const { return ArmadilloMatrixWrapper<TScalar>(+ m_Matrix); }



private :

	ArmadilloMatrixType m_Matrix;

};


#include "ArmadilloMatrixWrapper_friend.h"


#endif /* _ArmadilloMatrixWrapper_h */