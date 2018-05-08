#ifndef SOLVER_H
#define SOLVER_H

#include <omp.h>
#include <algorithm> 
#include <Eigen/Core>
#include <igl/decimate.h>

using namespace Eigen;
using std::vector;

class Solver
{
public:
    const int preIters=3, postIters=2, coarseIters=10;
    // const int preIters=30, postIters=20, coarseIters=100;
    int iters_ = 400;
    float tolerance = 1e-7;

    // https://github.com/libigl/libigl/blob/master/include/igl/decimate.cpp 
    // Edited so that the I map is not spliced. (I need the -1 values so that I know which
    // vertices were destroyed). 
    bool decimate(
        const Eigen::MatrixXd & V,
        const Eigen::MatrixXi & F,
        const size_t max_m,
        Eigen::MatrixXd & U,
        Eigen::MatrixXi & G,
        Eigen::VectorXi & J,
        Eigen::VectorXi & I)
    {
        // Original number of faces
        const int orig_m = F.rows();
        // Tracking number of faces
        int m = F.rows();
        typedef Eigen::MatrixXd DerivedV;
        typedef Eigen::MatrixXi DerivedF;
        DerivedV VO;
        DerivedF FO;
        igl::connect_boundary_to_infinity(V,F,VO,FO);
        // decimate will not work correctly on non-edge-manifold meshes. By extension
        // this includes meshes with non-manifold vertices on the boundary since these
        // will create a non-manifold edge when connected to infinity.
        if(!igl::is_edge_manifold(FO))
        {
          return false;
        }
        bool ret = igl::decimate(
          VO,
          FO,
          igl::shortest_edge_and_midpoint,
          igl::max_faces_stopping_condition(m,orig_m,max_m),
          U,
          G,
          J,
          I);
        const Eigen::Array<bool,Eigen::Dynamic,1> keep = (J.array()<orig_m);
        igl::slice_mask(Eigen::MatrixXi(G),keep,1,G);
        igl::slice_mask(Eigen::VectorXi(J),keep,1,J);
        VectorXi _1, I2;
        igl::remove_unreferenced(Eigen::MatrixXd(U),Eigen::MatrixXi(G),U,G,_1,I2);
        //igl::slice(Eigen::VectorXi(I),I2,1,I);
        return ret;

    }

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

    double gauss_seidel(const SparseMatrix<double>& A, const VectorXd& b, VectorXd& x,
            int max_iters)
    {
        // std::cout << "gauss_seidel size(): " << x.size() << std::endl;
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
            //std::cout << "error: " << error << std::endl;
            //std::cout << "tol: " << tolerance << std::endl;
        }
        // std::cout << "iters : " << iters << std::endl;
        return error;
    }

    std::vector<SparseMatrix<double>,Eigen::aligned_allocator<SparseMatrix<double>> > P;
    int max_depth;

    void multigrid_init(const MatrixXd& V, const MatrixXi& F)
    {
        MatrixXd fineV = V, coarseV;
        MatrixXi fineF = F, coarseF;
        VectorXi J,I;
        // Coarsen mesh
        while (decimate(fineV,fineF,fineF.rows()/2,coarseV,coarseF, J, I))
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
            int cnt = csize;
            for (int i = 0; i < A.size(); i++) {
                if (I[i] == -1) {
                    int ringCnt = 0;
                    // First sweep to count number of "valid" neighbors
                    for (int j = 0; j < A[i].size(); j++) {
                         if (I[A[i][j]] != -1) ringCnt++;
                    }
                    //// Second sweep to add prolongation row
                    for (int j = 0; j < A[i].size(); j++) {
                        int fineToCoarse = I[A[i][j]];
                        if (fineToCoarse != -1) // wrong index
                            triplets.push_back(Triplet<double>(cnt,fineToCoarse,1.0/ringCnt));
                    }
                    cnt++;
                }
            }
            std::cout << std::endl;
            p.setFromTriplets(triplets.begin(), triplets.end());
            P.push_back(p);
            std::cout << " coarseV.size() : " << coarseV.rows() << std::endl;
            fineV = coarseV;
            fineF = coarseF;
        }
        max_depth = P.size();
    }


    void multigrid(const SparseMatrix<double>& A, const VectorXd& b, VectorXd& x, int depth = 0)
    {
        if (depth == max_depth) {
            double error = gauss_seidel(A, b, x, coarseIters); 
            std::cout << "size: " << x.size() << " error: " <<  error << std::endl;
        } else {
            // Presmoothing
            double preerror = gauss_seidel(A, b, x, preIters); 
            std::cout << "size: " << x.size() << " preerror: " << preerror << std::endl;

            // Restriction
            SparseMatrix<double> R = P[depth].transpose();
            VectorXd r = R * (b - A*x);
            SparseMatrix<double> restrictA = R * A * P[depth];
            if (depth == max_depth-1) {
                // std::cout << "P Mat: \n" << P[depth] << std::endl;
                std::cout << "r: " << r << std::endl;
            }

            VectorXd cx(r.size());
            multigrid(restrictA, r, cx, depth+1); 
            if (depth == max_depth-1) {
                // std::cout << "P Mat: \n" << P[depth] << std::endl;
                std::cout << "cx: " << cx << std::endl;
            }

            // Prolong error and correct fine solution
            x += P[depth]*cx;

            // Post smoothing
            double posterror = gauss_seidel(A, b, x, postIters); 
            std::cout << "size: " << x.size() << " posterror: " << posterror << std::endl;
        }

    }
};

#endif // SOLVER_H
