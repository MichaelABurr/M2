// Copyright 2014  Michael E. Stillman

#ifndef M2_BASIC_MUTMAT_DMAT_LU_INPLACE_HPP_
#define M2_BASIC_MUTMAT_DMAT_LU_INPLACE_HPP_

#include "basic-mutable-matrices/dmat.hpp"
#include "basic-mutable-matrices/mat-arith.hpp"
#include "basic-mutable-matrices/mat-elem-ops.hpp"
#include "basic-mutable-matrices/mat-util.hpp"

#include <iostream>
#include <memory>
#include <vector>

template <typename RT>
class LUUtil
{
 public:
  typedef RT RingType;
  typedef DMat<RingType> Mat;

  static void setUpperLower(const Mat& LU, Mat& lower, Mat& upper);
  static void computePivotColumns(const Mat& LU,
                                  std::vector<size_t>& result_columns);
  // not written yet:
  // static bool computeSign(const std::vector<size_t>& perm);  // true: 1,
  // false: -1

  void debug_out(const Mat& M)
  {
    buffer o;
    displayMat(o, M);
    std::cout << o.str() << std::endl;
  }
};

template <typename RT>
class DMatLUinPlace;

template <typename RT>
class DMatLUinPlace
{
 public:
  typedef RT RingType;
  typedef DMat<RingType> Mat;

 public:
  DMatLUinPlace(const Mat& A);

  const RingType& ring() const { return mLU->ring(); }
  long numRows() const { return mLU->numRows(); }
  long numColumns() const { return mLU->numColumns(); }
  const Mat& LUinPlace()
  {
    computeLU();
    return *mLU;
  }  // raises an exception if there is an error
  // Can be called repeatedly: the result is remembered once done.
  // Returns a constant ref to the internal "in place" LU.

  bool signOfPermutation() { return mSign; }
  const std::vector<size_t>& permutation() { return mPerm; }
  const std::vector<size_t>& pivotColumns() { return mPivotColumns; }
 private:
  void computeLU();
  size_t findPivot(size_t row, size_t col);

 private:
  std::unique_ptr<Mat> mLU;
  std::vector<size_t> mPerm;
  bool mSign;
  bool mIsDone;
  std::vector<size_t> mPivotColumns;
};

template <class RingType>
DMatLUinPlace<RingType>::DMatLUinPlace(const Mat& A)
    : mLU(std::make_unique<Mat>(A)),  // copies A
      mSign(true),  // sign = 1
      mIsDone(false)
{
  for (size_t i = 0; i < A.numRows(); i++) mPerm.push_back(i);
}

template <class RingType>
size_t DMatLUinPlace<RingType>::findPivot(size_t row, size_t col)
{
  // Look at elements A[row,col], A[row+1,col], ..., A[nrows-1, col]
  // Return the index r s.y. abs(A[r,col]) is maximum over all of these

  for (size_t i = row; i < mLU->numRows(); i++)
    {
      if (!ring().is_zero(mLU->entry(i, col))) return i;
    }
  return static_cast<size_t>(-1);
}

template <class RingType>
void DMatLUinPlace<RingType>::computeLU()
{
  if (mIsDone) return;

  //  std::cout << "computing LU decomposition generic version" << std::endl;
  typename RingType::Element tmp(mLU->ring());

  size_t col = 0;  // current column we are working on
  size_t row = 0;  // current row we are working on
  size_t nrows = mLU->numRows();
  size_t ncols = mLU->numColumns();

  while (col < ncols && row < nrows)
    {
      // printf("*** in naive row,col = (%ld, %ld) ***\n", row, col);
      // debug_out();

      // Step 1: Set the 'upper' values: (row,col)..(nrows-1,col)
      for (size_t r = row; r < nrows; r++)
        {
          for (size_t i = 0; i < row; i++)
            {
              mLU->ring().mult(tmp, mLU->entry(r, i), mLU->entry(i, col));
              mLU->ring().subtract(
                  mLU->entry(r, col), mLU->entry(r, col), tmp);
            }
        }

      // printf("after step 1\n");
      // debug_out();

      // Step 2: Find a pivot among the elements in step 1.
      //  If one: swap rows if needed
      //  If none, increment 'col', and continue at start of loop
      size_t k = findPivot(row, col);
      if (k == static_cast<size_t>(-1))
        {
          col = col + 1;
          continue;
        }
      // printf("pivot is in row %ld\n", k);
      std::swap(mPerm[row], mPerm[k]);
      if (k != row)
        {
          MatElementaryOps<Mat>::interchange_rows(*mLU, k, row);
          mSign = !mSign;
        }
      mPivotColumns.push_back(col);
      const typename RingType::ElementType& pivot = mLU->entry(row, col);

      // printf("after step 2:\n");
      // debug_out();

      // Step 3A: Set the 'upper' elements in (row,col+1), ..., (row,ncols-1).
      for (size_t c = col + 1; c < ncols; c++)
        {
          for (size_t i = 0; i < row; i++)
            {
              mLU->ring().mult(tmp, mLU->entry(row, i), mLU->entry(i, c));
              mLU->ring().subtract(
                  mLU->entry(row, c), mLU->entry(row, c), tmp);
            }
        }
      // printf("after step 3A:\n");
      // debug_out();

      // Step 3B: Set the 'lower' elements in (row+1,row), ..., (nrows-1,row)
      //  from (row+1,col), ..., (nrows-1,col)
      // This just means dividing then by the pivot
      // except, if we have skipped columns for pivots, we must set these
      // elements
      // in column 'row', not 'col'...
      // Step 3C: if row != col, then set these elements to 0:
      //  (row+1,col), ..., (nrows-1,col)
      for (size_t r = row + 1; r < nrows; r++)
        {
          mLU->ring().divide(mLU->entry(r, row), mLU->entry(r, col), pivot);
          if (row < col) ring().set_zero(mLU->entry(r, col));
        }

      // printf("after step 3B:\n");
      // debug_out();

      row++;
      col++;
    }

  mIsDone = true;
}

template <class RingType>
void LUUtil<RingType>::setUpperLower(const Mat& LU, Mat& lower, Mat& upper)
{
  size_t min = std::min(LU.numRows(), LU.numColumns());
  lower.resize(LU.numRows(), min);
  upper.resize(min, LU.numColumns());

  // At this point, lower and upper should be zero matrices.
  assert(MatrixOps::isZero(lower));
  assert(MatrixOps::isZero(upper));

  for (size_t c = 0; c < LU.numColumns(); c++)
    {
      if (c < min) LU.ring().set(lower.entry(c, c), 1);
      for (size_t r = 0; r < LU.numRows(); r++)
        {
          if (r <= c)
            LU.ring().copy(upper.entry(r, c), LU.entry(r, c));
          else if (c < lower.numRows())
            {
              LU.ring().copy(lower.entry(r, c), LU.entry(r, c));
            }
        }
    }
}

template <class RingType>
void LUUtil<RingType>::computePivotColumns(const Mat& LU,
                                           std::vector<size_t>& result_columns)
{
  result_columns.clear();
  size_t thiscol = 0;
  size_t thisrow = 0;
  while (thisrow < LU.numRows() and thiscol < LU.numColumns())
    {
      if (not LU.ring().is_zero(LU.entry(thisrow, thiscol)))
        {
          result_columns.push_back(thiscol);
          thisrow++;
        }
      thiscol++;
    }
}

#endif

// Local Variables:
// compile-command: "make -C $M2BUILDDIR/Macaulay2/e "
// indent-tabs-mode: nil
// End:
