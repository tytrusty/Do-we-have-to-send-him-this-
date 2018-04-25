#ifndef SOLVER_H
#define SOLVER_H

using namespace Eigen;

namespace Solver
{
    int max_iters = 40;
    float tolerance = 1e-7;

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

    void gauss_seidel_parallel(const MatrixXd& A, const VectorXd& b, VectorXd& x)
    {
        int iters = 0;
        x.setZero();
        double dxi = 0.0;
        int size = x.size();
        double error = 1.0;
        VectorXd xold;
        while (iters < max_iters && error >= tolerance) {
            xold = x;
            // gauss-seidel red/black iteration
            for (int i = 0; i < size; i++) {
                dxi = b[i];
                for (int j = 0; j < size; j++) {
                    if (i != j)
                        dxi -= A.coeffRef(i,j)*x[j];
                }
                x[i] = dxi/A.coeffRef(i,i);
            }
            ++iters;
            error = (x - xold).norm();
        }
        std::cout << "iters : " << iters << std::endl;
    }

}

#endif // SOLVER_H