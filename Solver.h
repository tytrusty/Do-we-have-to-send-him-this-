#ifndef SOLVER_H
#define SOLVER_H

#include <omp.h>
#include <Eigen/Core>

using namespace Eigen;
using Eigen::internal::BandMatrix;

namespace Solver
{
    const int preIters=3, postIters=2, coarseIters=10;
    int iters_ = 400;
    float tolerance = 1e-7;

    void updateIters(int iters) { iters_ = iters; }
    void updateTolerance(int tol) { tolerance = tol; }
    /**
     * Fast(ish) error computation
     * This is orphan omp function so we're assuming it's in parallel already
     */
    void fast_error(const VectorXd& a, const VectorXd& b, double& norm) {
        // Assert a.size == b.size TODO
        norm = 0.0;
        #pragma omp for reduction(+:norm)
        for (int i = 0; i < a.size(); i++) {
            double diff = a[i]-b[i];
            norm += diff*diff;
        }
        norm = std::sqrt(norm);
    }

    void gauss_seidel(const SparseMatrix<double>& A, const VectorXd& b, VectorXd& x,
            int max_iters = iters_)
    {
        std::cout << "gauss_seidel size(): " << x.size() << std::endl;
        int iters = 0;
        x.setZero();
        double dxi = 0.0;
        int size = x.size();
        double error = 1.0;
        VectorXd xold;
        while (iters < max_iters && error >= tolerance) {
            xold = x;
            // #pragma omp parallel for private(dxi)  
            for (int i = 0; i < A.outerSize(); ++i) {
                dxi = b.coeff(i);
                for(SparseMatrix<double>::InnerIterator it(A,i); it; ++it) {
                    if (it.index() != i)
                        dxi -= it.value()*x[it.index()];
                }
                x[i] = dxi/A.coeff(i,i);
            }
            ++iters;
            fast_error(x, xold, error);
            std::cout << "error: " << error << std::endl;
            std::cout << "tol: " << tolerance << std::endl;
        }
        std::cout << "iters : " << iters << std::endl;
    }

    void restrict() {}
    void prolong() {}
    VectorXd restrict_residual(const SparseMatrix<double>& A, const VectorXd& b, const VectorXd& x) 
    {
        VectorXd r_old = b - A*x;
        SparseMatrix<double, RowMajor> restrict(r_old.size()/2, r_old.size());

        std::vector<Triplet<double>> triplets;
        triplets.reserve(r_old.size()/2 * 3);
        for (int i = 0; i < restrict.rows()-1; ++i) {
            triplets.push_back(Triplet<double>(i, (i*2),   0.25));
            triplets.push_back(Triplet<double>(i, (i*2)+1, 0.50));
            triplets.push_back(Triplet<double>(i, (i*2)+2, 0.25));
        }
        // handling edge case for even # of cols
        int last = restrict.rows()-1;
        triplets.push_back(Triplet<double>(last, (last*2),   0.25));
        triplets.push_back(Triplet<double>(last, (last*2)+1, 0.50));
        if (last*2+2 < restrict.cols())
            triplets.push_back(Triplet<double>(last, (last*2)+2, 0.25));
        restrict.setFromTriplets(triplets.begin(), triplets.end());
        return restrict * r_old;
    }

    VectorXd prolong_residual(const SparseMatrix<double>& A, const VectorXd& b, const VectorXd& x) 
    {
        VectorXd r_old = b - A*x;
        SparseMatrix<double, RowMajor> prolong(r_old.size()*2+1, r_old.size());

        std::vector<Triplet<double>> triplets;
        triplets.reserve(r_old.size()/2 * 3);
        for (int i = 0; i < prolong.cols()-1; ++i) {
            triplets.push_back(Triplet<double>((i*2),   i, 0.5));
            triplets.push_back(Triplet<double>((i*2)+1, i, 1.0));
            triplets.push_back(Triplet<double>((i*2)+2, i, 0.5));
        }
        // handling edge case for even # of cols
        int last = prolong.cols()-1;
        triplets.push_back(Triplet<double>(last*2, last,   0.5));
        triplets.push_back(Triplet<double>(last*2+1, last, 1.0));
        if (last*2+2 < prolong.cols())
            triplets.push_back(Triplet<double>(last*2+2, last, 0.5));

        prolong.setFromTriplets(triplets.begin(), triplets.end());
        return prolong * r_old;
    }
    void multigrid(const SparseMatrix<double>& A, const VectorXd& b, VectorXd& x, int depth)
    {
        if (depth == 0) {
            gauss_seidel(A, b, x, coarseIters); 
        } else {
            gauss_seidel(A, b, x, preIters); 
            VectorXd r = restrict_residual(A, b, x);
            gauss_seidel(A, b, x, postIters); 
        }
    }
}

#endif // SOLVER_H
