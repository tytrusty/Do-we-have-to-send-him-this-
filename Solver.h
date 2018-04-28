#ifndef SOLVER_H
#define SOLVER_H

#include <omp.h>

using namespace Eigen;

namespace Solver
{
    int max_iters = 40;
    float tolerance = 1e-7;

    /**
     * Fast(ish) error computation
     * This is orphan omp function so we're assuming it's in parallel already
     */
    void fast_error(const VectorXd& a, const VectorXd& b, double& norm) {
        // Assert a.size == b.size TODO
        // double norm = 0.0;
        #pragma omp for reduction(+:norm)
        for (int i = 0; i < a.size(); i++) {
            double diff = a[i]-b[i];
            norm += diff*diff;
        }
        norm = std::sqrt(norm);
    }

    void gauss_seidel_serial(const MatrixXd& A, const VectorXd& b, VectorXd& x)
    {
        int iters = 0;
        x.setZero();
        VectorXd dx(x.size());
        double dxi = 0.0;
        int size = x.size();
        double maxdiff = 1.0;
        while (iters < max_iters && maxdiff >= tolerance) {
            
            dx = b;
            // gauss-seidel red/black iteration
            for (int i = 0; i < size; i++) {
                dxi = 0.0;
                for (int j = 0; j < size; j++) 
                    dxi += A.coeffRef(i,j)*x[j];
                dx[i] = (b[i]-dxi) / A.coeffRef(i,i);
                x[i] += dx[i];
            }
            ++iters;
        }
        std::cout << "maxdiff: " <<  maxdiff << std::endl;
        std::cout << "iters : " << iters << std::endl;
    }

    void gauss_seidel_parallel(const SparseMatrix<double>& A, const VectorXd& b, VectorXd& x)
    {
        int iters = 0;
        x.setZero();
        double dxi = 0.0;
        int size = x.size();
        double error = 1.0;
        VectorXd xold;
        while (iters < max_iters && error >= tolerance) {
            xold = x;
            #pragma omp parallel for private(dxi)  
            for (int i = 0; i < A.outerSize(); ++i) {
                dxi = b[i];
                for(SparseMatrix<double>::InnerIterator it(A,i); it; ++it) {
                    if (it.index() != i)
                        dxi -= it.value()*x[it.index()];
                }
                x[i] = dxi/A.coeff(i,i);
            }
            ++iters;
            fast_error(x, xold, error);
        }
        std::cout << "iters : " << iters << std::endl;
    }
}

#endif // SOLVER_H
