/*
  This file is part of MADNESS.
  
  Copyright (C) 2007-10 Oak Ridge National Laboratory
  
  This program is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2 of the License, or
  (at your option) any later version.
  
  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
  GNU General Public License for more details.
  
  You should have received a copy of the GNU General Public License
  along with this program; if not, write to the Free Software
  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
  
  For more information please contact:

  Robert J. Harrison
  Oak Ridge National Laboratory
  One Bethel Valley Road
  P.O. Box 2008, MS-6367

  email: harrisonrj@ornl.gov 
  tel:   865-241-3937
  fax:   865-572-0680

  
  $Id$
*/

#ifndef MADNESS_LINALG_GMRES_H__INCLUDED
#define MADNESS_LINALG_GMRES_H__INCLUDED

#include <tensor/tensor.h>
#include <world/print.h>
#include <iostream>
#include <linalg/tensor_lapack.h>

/// \file gmres.h
/// \brief Defines a general operator interface and a templated GMRES solver
/// for solving linear equations.
///
/// \ingroup solvers
///
/// An AbstractVectorSpace class is also defined to guarantee a uniform way
/// for linear algebra routines to access common operations like norm, inner
/// product etc.  Implementations for MADNESS Vectors and Functions are
/// provided in the VectorSpace<floating_point, NDIM> and
/// FunctionSpace<floating_point, NDIM> classes.

namespace madness {

    /// A generic operator: takes in one T and produces another T.
    /// Override the protected action() function to implement an Operator.
    template <typename T>
    class Operator {
        protected:
            /// The action of the operator
            virtual void action(const T &in, T &out) const = 0;

        public:
            /// public access to the operator's action, returns out for
            /// convenience
            T & applyOp(const T &in, T &out) const {
                action(in, out);
                return out;
            }

	    virtual ~Operator(){}
    };

    /// A generic vector space which provides common operations needed by
    /// linear algebra routines (norm, inner product, etc.)
    ///
    /// Most of these routines will presumedly be defined, but the names may
    /// not be the same (i.e. madness:Function uses norm2, madness::Tensor
    /// uses normf)
    ///
    /// T is the vector type
    /// real_type the real type used for norms, etc.
    /// scalar_type is the type used for inner products and gaxpy coefficients,
    /// etc.
    ///
    /// When implementing a child class, real_type and scalar_type will
    /// probably be obtained from TensorTypeData<T> (see the VectorSpace
    /// and FunctionSpace classes below)
    template <typename T, typename real_type, typename scalar_type>
    class AbstractVectorSpace {
        public:
            virtual ~AbstractVectorSpace() {}

            /// the norm of a vector
            virtual real_type norm(const T &) const = 0;

            /// scale the vector by a constant
            /// should return the scaled vector
            virtual T & scale(T &, const scalar_type &) const = 0;

            /// standard bilinear gaxpy
            ///   x <- a*x + b*y
            /// should return the new vector, x
            virtual T & gaxpy(T &x, const scalar_type &a, const T &y,
                const scalar_type &b) const = 0;

            /// the inner product between two vectors
            virtual scalar_type inner(const T &, const T &) const = 0;

            /// any special instructions to be executed when a vector is no
            /// longer needed
            ///
            /// unless otherwise specified, do nothing
            virtual void destroy(T &) const {}
    };

    /// Real part of a real number
    /// This is necessary to make templating nightmares disappear in
    /// VectorSpace::norm
    static inline double real(double x) { return x; }

    /// Imaginary part of a real number
    static inline double imag(double x) { return 0.0; }

    /// A vector space using MADNESS Vectors
    template <typename T, int NDIM>
    class VectorSpace : public AbstractVectorSpace<Vector<T, NDIM>,
        typename TensorTypeData<T>::float_scalar_type, T> {

        protected:
            static const bool iscplx = TensorTypeData<T>::iscomplex;

        public:
            typedef typename TensorTypeData<T>::float_scalar_type real_type;
            typedef T scalar_type;

            virtual ~VectorSpace() {}

            virtual real_type norm(const Vector<scalar_type, NDIM> &vec)
                const {

                real_type ret, mag;
                int i;

                mag = real(vec[0]);
                ret = mag*mag;
                if(iscplx) {
                    mag = imag(vec[0]);
                    ret += mag*mag;
                }

                for(i = 1; i < NDIM; ++i) {
                    mag = real(vec[i]);
                    ret += mag*mag;
                    if(iscplx) {
                        mag = imag(vec[i]);
                        ret += mag*mag;
                    }
                }

                return sqrt(ret);
            }

            virtual Vector<scalar_type, NDIM> & scale(
                Vector<scalar_type, NDIM> &vec, const scalar_type &c) const {

                vec *= c;
                return vec;
            }

            virtual Vector<scalar_type, NDIM> & gaxpy(
                Vector<scalar_type, NDIM> &x, const scalar_type &a,
                const Vector<scalar_type, NDIM> &y, const scalar_type &b)
                const {

                for(int i = 0; i < NDIM; ++i)
                    x[i] = a * x[i] + b * y[i];
                return x;
            }

            virtual scalar_type inner(const Vector<scalar_type, NDIM> &l,
                const Vector<scalar_type, NDIM> &r) const {

                scalar_type ret;

                if(iscplx) {
                    ret = conj(l[0]) * r[0];
                    for(int i = 1; i < NDIM; ++i)
                        ret += conj(l[i]) * r[i];
                }
                else {
                    ret  = l[0] * r[0];
                    for(int i = 1; i < NDIM; ++i)
                        ret += l[i] * r[i];
                }

                return ret;
            }
    };

    /// A vector space using MADNESS Functions
    template <typename T, int NDIM>
    class FunctionSpace : public AbstractVectorSpace<Function<T, NDIM>,
        typename TensorTypeData<T>::float_scalar_type, T> {

        public:
            typedef typename TensorTypeData<T>::float_scalar_type real_type;
            typedef T scalar_type;

            virtual ~FunctionSpace() {}

            virtual real_type norm(const Function<scalar_type, NDIM> &vec)
                const {

                return vec.norm2();
            }

            virtual Function<scalar_type, NDIM> & scale(
                Function<scalar_type, NDIM> &vec, const scalar_type &c) const {

                vec.scale(c);
                return vec;
            }

            virtual Function<scalar_type, NDIM> & gaxpy(
                Function<scalar_type, NDIM> &x, const scalar_type &a,
                const Function<scalar_type, NDIM> &y, const scalar_type &b)
                const {

                x.gaxpy(a, y, b);
                return x;
            }

            virtual scalar_type inner(const Function<scalar_type, NDIM> &l,
                const Function<scalar_type, NDIM> &r) const {

                return l.inner(r);
            }

            virtual void destroy(Function<scalar_type, NDIM> &f) const {
                f.clear();
            }
    };

    /// A vector space using MADNESS Vectors of MADNESS Functions
    template <typename T, int VDIM, int FDIM>
    class VectorOfFunctionsSpace : public AbstractVectorSpace<
        Vector<Function<T, FDIM>, VDIM>,
        typename TensorTypeData<T>::float_scalar_type, T> {

        public:
            typedef typename TensorTypeData<T>::float_scalar_type real_type;
            typedef T scalar_type;

            virtual ~VectorOfFunctionsSpace() {}

            virtual real_type norm(
                    const Vector<Function<scalar_type, FDIM>, VDIM> &vec)
                const {

                real_type ret = vec[0].norm2();
                for(int i = 1; i < VDIM; ++i)
                    ret += vec[i].norm2();
                return ret;
            }

            virtual Vector<Function<scalar_type, FDIM>, VDIM> & scale(
                    Vector<Function<scalar_type, FDIM>, VDIM> &vec,
                    const scalar_type &c) const {

                for(int i = 0; i < VDIM; ++i)
                    vec[i].scale(c);
                return vec;
            }

            virtual Vector<Function<scalar_type, FDIM>, VDIM> & gaxpy(
                    Vector<Function<scalar_type, FDIM>, VDIM> &x,
                    const scalar_type &a,
                    const Vector<Function<scalar_type, FDIM>, VDIM> &y,
                    const scalar_type &b) const {

                for(int i = 0; i < VDIM; ++i)
                    x[i].gaxpy(a, y[i], b);
                return x;
            }

            virtual scalar_type inner(
                    const Vector<Function<scalar_type, FDIM>, VDIM> &l,
                    const Vector<Function<scalar_type, FDIM>, VDIM> &r) const {

                scalar_type ret = l[0].inner(r[0]);
                for(int i = 0; i < VDIM; ++i)
                    ret += l[i].inner(r[i]);
                return ret;
            }

            virtual void destroy(Vector<Function<scalar_type, FDIM>, VDIM> &f)
                    const {
                for(int i = 0; i < VDIM; ++i)
                    f[i].clear();
            }
    };

    /// A GMRES solver routine for linear systems, Ax == b.
    /// Requires the vector space, the operator, A, the inhomogeneity b, an
    /// initial guess x, the maximum number of iterations, the convergence
    /// threshold, and a flag for producing output.
    ///
    /// The vector space object provides a standard way to compute the needed
    /// linear algebra operations (norm, inner product, etc.).
    ///
    /// The printed output, if desired, shows iteration number, the residual,
    /// and the effective rank of the GMRES matrix at that iteration.
    ///
    /// On output, x has the computed solution, thresh has its residual,
    /// and maxiters has the number of iterations needed
    template <typename T, typename real_type, typename scalar_type>
    void GMRES(const AbstractVectorSpace<T, real_type, scalar_type> &space,
        const Operator<T> &op, const T &b, T &x, int &maxiters,
        real_type &thresh, const bool outp = true) {

        int iter, i;
        long rank;
        std::vector<T> V;
        T r;
        Tensor<scalar_type> H(maxiters+1, maxiters);
        Tensor<scalar_type> betae(maxiters+1);
        Tensor<scalar_type> y;
        Tensor<real_type> s;
        Tensor<real_type> sumsq;
        real_type resid, norm;

        // initialize
        H = 0.0;
        betae = 0.0;

        // construct the first subspace basis vector
        iter = 0;
        op.applyOp(x, r);
        space.gaxpy(r, -1.0, b, 1.0);
        betae[0] = resid = space.norm(r);
        if(outp)
            printf("itr rnk resid\n%.3d N/A %.6e\n", iter, resid);
        if(resid < thresh) {
            maxiters = 1;
            thresh = resid;
            return;
        }
        space.scale(r, 1.0 / resid);
        V.push_back(r);
        ++iter;

        do {
            // compute the new vector
            op.applyOp(V[iter - 1], r);

            // orthogonalize the new vector
            for(i = 0; i < iter; ++i) {
                H(i, iter-1) = space.inner(V[i], r);
                space.gaxpy(r, 1.0, V[i], -H(i, iter-1));
            }
            H(iter, iter-1) = norm = space.norm(r);

            // solve Hy == betae for y
            gelss(H(Slice(0, iter), Slice(0, iter-1)), betae(Slice(0, iter)),
                1.0e-12, y, s, rank, sumsq);

            resid = sumsq[0];

            if(norm < thresh) {
                // we just got the zero vector -> no more progress

                // if this is the first iteration, compute the residual
                if(iter == 1) {
                    space.destroy(r);
                    space.gaxpy(x, 1.0, V[0], betae[0]);
                    op.applyOp(x, r);
                    space.gaxpy(r, -1.0, b, 1.0);
                    thresh = space.norm(r);
                    space.destroy(r);
                    space.destroy(V[0]);

                    maxiters = 1;
                    if(outp)
                        printf("%.3d N/A %.6e ** Zero Vector Encountered **\n",
                            iter, thresh);
                    return;
                }

                if(outp)
                    printf("%.3d %.3ld %.6e ** Zero Vector Encountered **\n",
                        iter, rank, resid);
                break;
            }

            if(outp) {
                printf("%.3d %.3ld %.6e", iter, rank, resid);
                if(iter != rank)
                    printf(" ** Questionable Progress **");
                printf("\n");
            }

            // normalize the vector and save it
            space.scale(r, 1.0 / norm);
            V.push_back(r);
            ++iter;
       } while(iter < maxiters && resid > thresh);

       // build the solution and destroy the basis vectors
       for(i = 0; i < iter; ++i) {
           space.gaxpy(x, 1.0, V[i], y[i]);
           space.destroy(V[i]);
       }

       thresh = resid;
       maxiters = iter;
    }

} // end of madness namespace

#endif // MADNESS_LINALG_GMRES_H__INCLUDED
