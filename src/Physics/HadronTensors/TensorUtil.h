/*
 * =====================================================================================
 *
 *       Filename:  TensorUtil.h
 *
 *    Description:  Tiny fixed-size linear-algebra library with an Eigen-shaped
 *                  API.  Header-only, no external dependencies beyond the C++
 *                  standard library.  Lets the hadron-tensor code in this
 *                  directory be written in the usual Eigen idiom (Matrix4cd,
 *                  Vector4d, .block<>(), .head<>(), .transposeInPlace(),
 *                  comma-initialization, etc.) without pulling Eigen in as a
 *                  new GENIE dependency.
 *
 *                  Provides (in namespace genie::TensorUtil):
 *
 *                     Matrix2cd     - 2x2  complex<double>
 *                     Matrix4cd     - 4x4  complex<double>
 *                     Vector2cd     - len-2 complex<double>
 *                     Vector4cd     - len-4 complex<double>
 *                     RowVector2cd  - len-2 complex<double> row vector
 *                     Vector4d      - len-4 double
 *
 *                  Supported subset of the Eigen-shaped API:
 *                     m.setZero(), v.setZero()
 *                     M::Zero()
 *                     m(i,j), v(i)
 *                     m1+m2, m1-m2, m1*m2 (matrix*matrix)
 *                     M*v       (matrix * column)
 *                     row*M     (row * matrix)
 *                     scalar*m, m*scalar, m/scalar, +=, -=, *=, /=
 *                     m.block<R,C>(i,j) = src                   (lvalue)
 *                     m.block<R,C>(i,j)                         (rvalue)
 *                     v.head<N>(), v.segment<N>(i), v.tail<N>() (lvalue + rvalue)
 *                     v.transpose()                              (column -> row)
 *                     M.transposeInPlace()
 *                     m << a, b, c, d ... (comma-initializer, row-major)
 *
 *                  Storage is std::array<T, R*C> (or std::array<T, N> for
 *                  vectors), inline and stack-allocated.  No heap, no virtual
 *                  dispatch.  Requires C++11.
 *
 *        Version:  1.0
 *        Created:  05/07/2026 02:40:52 PM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Liang Liu (L. Liu), liangliu@fnal.gov
 *		    Fermi National Accelerator Laboratory
 *  Collaboration:  GENIE
 *
 * =====================================================================================
 */


#ifndef _TENSOR_UTIL_H_
#define _TENSOR_UTIL_H_

#include <array>
#include <complex>
#include <cstddef>
#include <type_traits>
#include <utility>

namespace genie {
  namespace TensorUtil {

    // =============================================================================
    // Forward declarations
    // =============================================================================
    template <typename T, int R, int C> class SmallMatrix;
    template <typename T, int N>        class SmallVector;
    template <typename T, int N>        class SmallRowVector;


    // =============================================================================
    // Comma-initializer (m << a, b, c, ...)
    //
    // The host class is responsible for translating linear index k into an
    // element write through the call host.set_linear(k, value).
    // =============================================================================
    template <typename Host>
      class CommaInit
      {
        public:
          using Scalar = typename Host::Scalar;

          CommaInit(Host & host, const Scalar & first)
            : m_host(host), m_index(0)
          {
            m_host.set_linear(m_index++, first);
          }

          CommaInit & operator,(const Scalar & v)
          {
            m_host.set_linear(m_index++, v);
            return *this;
          }

        private:
          Host &      m_host;
          std::size_t m_index;
      };


    // =============================================================================
    // Block proxy: lvalue / rvalue view over a sub-block of a matrix
    // =============================================================================
    template <typename Mat, int Rows, int Cols>
      class Block
      {
        public:
          using Scalar = typename Mat::Scalar;

          Block(Mat & host, int r0, int c0)
            : m_host(host), m_r0(r0), m_c0(c0) {}

          // Assign from any source providing operator()(i,j).
          template <typename Src>
            Block & operator=(const Src & src)
            {
              for(int i = 0; i < Rows; ++i)
                for(int j = 0; j < Cols; ++j)
                  m_host(m_r0 + i, m_c0 + j) = src(i, j);
              return *this;
            }

          Scalar operator()(int i, int j) const
          {
            return m_host(m_r0 + i, m_c0 + j);
          }

          // Materialize as a standalone SmallMatrix when used in expressions.
          operator SmallMatrix<Scalar, Rows, Cols>() const
          {
            SmallMatrix<Scalar, Rows, Cols> out;
            for(int i = 0; i < Rows; ++i)
              for(int j = 0; j < Cols; ++j)
                out(i, j) = m_host(m_r0 + i, m_c0 + j);
            return out;
          }

        private:
          Mat & m_host;
          int   m_r0;
          int   m_c0;
      };


    // =============================================================================
    // Vector segment proxy (head<N>, segment<N>, tail<N>)
    // =============================================================================
    template <typename Vec, int N>
      class Segment
      {
        public:
          using Scalar = typename Vec::Scalar;

          Segment(Vec & host, int i0) : m_host(host), m_i0(i0) {}

          template <typename Src>
            Segment & operator=(const Src & src)
            {
              for(int i = 0; i < N; ++i)
                m_host(m_i0 + i) = src(i);
              return *this;
            }

          Scalar operator()(int i) const { return m_host(m_i0 + i); }

          operator SmallVector<Scalar, N>() const
          {
            SmallVector<Scalar, N> out;
            for(int i = 0; i < N; ++i) out(i) = m_host(m_i0 + i);
            return out;
          }

        private:
          Vec & m_host;
          int   m_i0;
      };


    // =============================================================================
    // SmallMatrix<T, Rows, Cols>
    //
    // Storage: std::array<T, Rows*Cols>, row-major, inline (no heap).
    // =============================================================================
    template <typename T, int Rows, int Cols>
      class SmallMatrix
      {
        public:
          using Scalar = T;
          static constexpr int RowsAtCompileTime = Rows;
          static constexpr int ColsAtCompileTime = Cols;

          // -- ctors ---------------------------------------------------------------
          SmallMatrix() { m_data.fill(T(0)); }

          // -- factories -----------------------------------------------------------
          static SmallMatrix Zero() { return SmallMatrix(); }

          // -- mutation ------------------------------------------------------------
          void setZero() { m_data.fill(T(0)); }

          // -- element access ------------------------------------------------------
          T &       operator()(int i, int j)       { return m_data[i * Cols + j]; }
          const T & operator()(int i, int j) const { return m_data[i * Cols + j]; }

          // For the comma-initializer (row-major linear index).
          void set_linear(std::size_t k, const T & v) { m_data[k] = v; }

          // -- arithmetic ----------------------------------------------------------
          SmallMatrix & operator+=(const SmallMatrix & rhs)
          {
            for(std::size_t k = 0; k < m_data.size(); ++k) m_data[k] += rhs.m_data[k];
            return *this;
          }
          SmallMatrix & operator-=(const SmallMatrix & rhs)
          {
            for(std::size_t k = 0; k < m_data.size(); ++k) m_data[k] -= rhs.m_data[k];
            return *this;
          }
          SmallMatrix & operator*=(const T & s)
          {
            for(auto & x : m_data) x *= s;
            return *this;
          }
          SmallMatrix & operator/=(const T & s)
          {
            for(auto & x : m_data) x /= s;
            return *this;
          }

          // -- in-place transpose (square only) -----------------------------------
          void transposeInPlace()
          {
            static_assert(Rows == Cols, "transposeInPlace requires a square matrix");
            for(int i = 0; i < Rows; ++i)
              for(int j = i + 1; j < Cols; ++j)
                std::swap((*this)(i, j), (*this)(j, i));
          }

          // -- block view ---------------------------------------------------------
          template <int R, int C>
            Block<SmallMatrix, R, C> block(int r0, int c0)
            {
              return Block<SmallMatrix, R, C>(*this, r0, c0);
            }
          template <int R, int C>
            Block<const SmallMatrix, R, C> block(int r0, int c0) const
            {
              return Block<const SmallMatrix, R, C>(*this, r0, c0);
            }

          // -- comma-initializer ---------------------------------------------------
          CommaInit<SmallMatrix> operator<<(const T & first)
          {
            return CommaInit<SmallMatrix>(*this, first);
          }

        private:
          std::array<T, Rows * Cols> m_data;
      };


    // =============================================================================
    // SmallVector<T, N> (column vector)
    // =============================================================================
    template <typename T, int N>
      class SmallVector
      {
        public:
          using Scalar = T;
          static constexpr int SizeAtCompileTime = N;

          SmallVector() { m_data.fill(T(0)); }

          static SmallVector Zero() { return SmallVector(); }

          void setZero() { m_data.fill(T(0)); }

          // element access
          T &       operator()(int i)       { return m_data[i]; }
          const T & operator()(int i) const { return m_data[i]; }
          T &       operator[](int i)       { return m_data[i]; }
          const T & operator[](int i) const { return m_data[i]; }

          void set_linear(std::size_t k, const T & v) { m_data[k] = v; }

          // arithmetic
          SmallVector & operator+=(const SmallVector & rhs)
          {
            for(int i = 0; i < N; ++i) m_data[i] += rhs.m_data[i];
            return *this;
          }
          SmallVector & operator-=(const SmallVector & rhs)
          {
            for(int i = 0; i < N; ++i) m_data[i] -= rhs.m_data[i];
            return *this;
          }
          SmallVector & operator*=(const T & s)
          {
            for(auto & x : m_data) x *= s;
            return *this;
          }
          SmallVector & operator/=(const T & s)
          {
            for(auto & x : m_data) x /= s;
            return *this;
          }

          // segment views
          template <int K>
            Segment<SmallVector, K> head()             { return Segment<SmallVector, K>(*this, 0); }
          template <int K>
            Segment<SmallVector, K> segment(int i0)    { return Segment<SmallVector, K>(*this, i0); }
          template <int K>
            Segment<SmallVector, K> tail()             { return Segment<SmallVector, K>(*this, N - K); }

          template <int K>
            Segment<const SmallVector, K> head() const          { return Segment<const SmallVector, K>(*this, 0); }
          template <int K>
            Segment<const SmallVector, K> segment(int i0) const { return Segment<const SmallVector, K>(*this, i0); }
          template <int K>
            Segment<const SmallVector, K> tail() const          { return Segment<const SmallVector, K>(*this, N - K); }

          // transpose -> row vector
          SmallRowVector<T, N> transpose() const
          {
            SmallRowVector<T, N> r;
            for(int i = 0; i < N; ++i) r(i) = m_data[i];
            return r;
          }

          // comma-init
          CommaInit<SmallVector> operator<<(const T & first)
          {
            return CommaInit<SmallVector>(*this, first);
          }

        private:
          std::array<T, N> m_data;
      };


    // =============================================================================
    // SmallRowVector<T, N>
    // =============================================================================
    template <typename T, int N>
      class SmallRowVector
      {
        public:
          using Scalar = T;

          SmallRowVector() { m_data.fill(T(0)); }

          T &       operator()(int i)       { return m_data[i]; }
          const T & operator()(int i) const { return m_data[i]; }

          void set_linear(std::size_t k, const T & v) { m_data[k] = v; }

          CommaInit<SmallRowVector> operator<<(const T & first)
          {
            return CommaInit<SmallRowVector>(*this, first);
          }

        private:
          std::array<T, N> m_data;
      };


    // =============================================================================
    // Free-function arithmetic
    // =============================================================================

    // matrix +/- matrix
    template <typename T, int R, int C>
      SmallMatrix<T, R, C> operator+(SmallMatrix<T, R, C> a, const SmallMatrix<T, R, C> & b)
      {
        a += b; return a;
      }
    template <typename T, int R, int C>
      SmallMatrix<T, R, C> operator-(SmallMatrix<T, R, C> a, const SmallMatrix<T, R, C> & b)
      {
        a -= b; return a;
      }
    template <typename T, int R, int C>
      SmallMatrix<T, R, C> operator-(const SmallMatrix<T, R, C> & a)
      {
        SmallMatrix<T, R, C> out;
        for(int i = 0; i < R; ++i)
          for(int j = 0; j < C; ++j)
            out(i, j) = -a(i, j);
        return out;
      }

    // matrix * matrix  (R x K) * (K x C) -> (R x C)
    template <typename T, int R, int K, int C>
      SmallMatrix<T, R, C> operator*(const SmallMatrix<T, R, K> & a,
          const SmallMatrix<T, K, C> & b)
      {
        SmallMatrix<T, R, C> out;
        for(int i = 0; i < R; ++i)
        {
          for(int j = 0; j < C; ++j)
          {
            T s = T(0);
            for(int k = 0; k < K; ++k) s += a(i, k) * b(k, j);
            out(i, j) = s;
          }
        }
        return out;
      }

    // matrix * column  (R x C) * (C) -> (R)
    template <typename T, int R, int C>
      SmallVector<T, R> operator*(const SmallMatrix<T, R, C> & a,
          const SmallVector<T, C> & v)
      {
        SmallVector<T, R> out;
        for(int i = 0; i < R; ++i)
        {
          T s = T(0);
          for(int j = 0; j < C; ++j) s += a(i, j) * v(j);
          out(i) = s;
        }
        return out;
      }

    // row * matrix  (1 x R) * (R x C) -> (1 x C)
    template <typename T, int R, int C>
      SmallRowVector<T, C> operator*(const SmallRowVector<T, R> & v,
          const SmallMatrix<T, R, C> & a)
      {
        SmallRowVector<T, C> out;
        for(int j = 0; j < C; ++j)
        {
          T s = T(0);
          for(int i = 0; i < R; ++i) s += v(i) * a(i, j);
          out(j) = s;
        }
        return out;
      }

    // scalar * matrix and matrix * scalar (S convertible to T)
    template <typename T, int R, int C, typename S>
      auto operator*(const S & s, SmallMatrix<T, R, C> a)
      -> typename std::enable_if<std::is_convertible<S, T>::value, SmallMatrix<T, R, C>>::type
      {
        a *= T(s); return a;
      }
    template <typename T, int R, int C, typename S>
      auto operator*(SmallMatrix<T, R, C> a, const S & s)
      -> typename std::enable_if<std::is_convertible<S, T>::value, SmallMatrix<T, R, C>>::type
      {
        a *= T(s); return a;
      }
    template <typename T, int R, int C, typename S>
      auto operator/(SmallMatrix<T, R, C> a, const S & s)
      -> typename std::enable_if<std::is_convertible<S, T>::value, SmallMatrix<T, R, C>>::type
      {
        a /= T(s); return a;
      }

    // scalar * vector, vector * scalar, vector / scalar
    template <typename T, int N, typename S>
      auto operator*(const S & s, SmallVector<T, N> v)
      -> typename std::enable_if<std::is_convertible<S, T>::value, SmallVector<T, N>>::type
      {
        v *= T(s); return v;
      }
    template <typename T, int N, typename S>
      auto operator*(SmallVector<T, N> v, const S & s)
      -> typename std::enable_if<std::is_convertible<S, T>::value, SmallVector<T, N>>::type
      {
        v *= T(s); return v;
      }
    template <typename T, int N, typename S>
      auto operator/(SmallVector<T, N> v, const S & s)
      -> typename std::enable_if<std::is_convertible<S, T>::value, SmallVector<T, N>>::type
      {
        v /= T(s); return v;
      }

    // Segment + Segment / Segment - Segment (returns concrete SmallVector).
    // Used e.g. for `pp.tail<3>() = p.tail<3>() + q.tail<3>()`.
    template <typename VecA, typename VecB, int N>
      SmallVector<typename VecA::Scalar, N>
      operator+(const Segment<VecA, N> & a, const Segment<VecB, N> & b)
      {
        SmallVector<typename VecA::Scalar, N> out;
        for(int i = 0; i < N; ++i) out(i) = a(i) + b(i);
        return out;
      }
    template <typename VecA, typename VecB, int N>
      SmallVector<typename VecA::Scalar, N>
      operator-(const Segment<VecA, N> & a, const Segment<VecB, N> & b)
      {
        SmallVector<typename VecA::Scalar, N> out;
        for(int i = 0; i < N; ++i) out(i) = a(i) - b(i);
        return out;
      }


    // =============================================================================
    // Eigen-compatible typedefs (the names client code uses)
    // =============================================================================
    using Matrix2cd    = SmallMatrix<std::complex<double>, 2, 2>;
    using Matrix4cd    = SmallMatrix<std::complex<double>, 4, 4>;
    using Vector2cd    = SmallVector<std::complex<double>, 2>;
    using Vector4cd    = SmallVector<std::complex<double>, 4>;
    using RowVector2cd = SmallRowVector<std::complex<double>, 2>;
    using Vector4d     = SmallVector<double, 4>;

  }  // namespace TensorUtil
}  // namespace genie


#endif  // _TENSOR_UTIL_H_
