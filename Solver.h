#ifndef SOLVER_H
#define SOLVER_H

#include <omp.h>
#include <algorithm> 
#include <Eigen/Core>
#include <igl/slice.h>
#include <igl/edge_flaps.h>
#include <igl/collapse_edge.h>
#include <igl/adjacency_list.h>
#include <igl/remove_unreferenced.h>
#include <igl/shortest_edge_and_midpoint.h>
#include <igl/max_faces_stopping_condition.h>
#include <igl/decimate.h>

using namespace Eigen;
using std::vector;

class Solver
{
public:
    const int preIters=3, postIters=2, coarseIters=10;
    int iters_ = 400;
    float tolerance = 1e-7;

    void updateIters(int iters) { iters_ = iters; } void updateTolerance(int tol) { tolerance = tol; }
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
            int max_iters)
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

    std::vector<SparseMatrix<double>,Eigen::aligned_allocator<SparseMatrix<double>> > operators;

    void multigrid_init(const MatrixXd& V, const MatrixXi& F)
    {
        MatrixXd fineV = V, coarseV;
        MatrixXi fineF = F, coarseF;
        VectorXi J,I;
        // Coarsen mesh
        while (igl::decimate(fineV,fineF,fineF.rows()/2,coarseV,coarseF, J, I))
        {
            const int fsize = fineV.rows();
            const int csize = coarseV.rows();
            SparseMatrix<double, RowMajor> p(fsize, csize);
            std::vector<Triplet<double>> triplets;
            // Form Two blocks of prolongation operator -- an Identity block coarse x coarse
            // and a mapping from fine -> coarse vertices
            for(int i = 0; i < csize; ++i)
                triplets.push_back(Triplet<double>(i,i,1));

            vector<vector<int>> A;
            igl::adjacency_list(fineF,A);
            // loop through invalid vertices and add rows to prolongation matrix that
            // map a fine vertex to its neighbors in the coarse mesh.
            int cnt = 0;
            for (int i = 0; i < I.size(); i++) {
                if (I[i] == -1) {
                    int ringCnt = 0;
                    // First sweep to count number of "valid" neighbors
                    for (int j = 0; j < A[i].size(); j++) {
                        if (I[A[i][j]] != -1) ringCnt++;
                    }
                    // Second sweep to add prolongation row
                    for (int j = 0; j < A[i].size(); j++) {
                        int fineToCoarse = I[A[i][j]];
                        if (fineToCoarse != -1)
                            triplets.push_back(Triplet<double>(cnt,fineToCoarse,1.0/ringCnt));
                    }
                    cnt++;
                }
            }
            p.setFromTriplets(triplets.begin(), triplets.end());
            operators.push_back(p);
            std::cout << " coarseV.size() : " << coarseV.rows() << std::endl;
            fineV = coarseV;
            fineF = coarseF;
        }
    }


    void multigrid(const SparseMatrix<double>& A, const VectorXd& b, VectorXd& x, int depth)
    {
        if (depth == 0) {
            gauss_seidel(A, b, x, coarseIters); 
        } else {
            gauss_seidel(A, b, x, preIters); 
            // VectorXd r = restrict_residual(A, b, x);
            gauss_seidel(A, b, x, postIters); 
        }
    }
};

#endif // SOLVER_H
