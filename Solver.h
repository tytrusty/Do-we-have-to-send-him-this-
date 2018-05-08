#ifndef SOLVER_H
#define SOLVER_H

#include <omp.h>
#include <Eigen/Core>
#include <igl/edge_flaps.h>
#include <igl/collapse_edge.h>
#include <igl/max_faces_stopping_condition.h>
#include <igl/shortest_edge_and_midpoint.h>
#include <igl/slice.h>
#include <igl/remove_unreferenced.h>
#include <algorithm>

using namespace Eigen;

class Solver
{
public:
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

    VectorXi EMAP;
    MatrixXi E,EF,EI;
    typedef std::set<std::pair<double,int> > PriorityQueue;
    PriorityQueue Q;
    std::vector<PriorityQueue::iterator > Qit;
    std::vector<SparseMatrix<double>,Eigen::aligned_allocator<SparseMatrix<double>> > operators;
    MatrixXd C;

    void setup(const MatrixXd& V, const MatrixXi& F) {
        // If an edge were collapsed, we'd collapse it to these points:
        igl::edge_flaps(F, E, EMAP, EF, EI);
        Qit.resize(E.rows());

        // Cost matrix
        C.resize(E.rows(),V.cols());
        VectorXd costs(E.rows());

        // Build cost matrix
        for(int e = 0;e<E.rows();e++)
        {
          double cost = e;
          RowVectorXd p(1,3);
          igl::shortest_edge_and_midpoint(e,V,F,E,EMAP,EF,EI,cost,p);
          C.row(e) = p;
          Qit[e] = Q.insert(std::pair<double,int>(cost,e)).first;
        }

    }
    
    void construct_p(const MatrixXd& V, const MatrixXi& F, const size_t max_m, MatrixXd& VO, MatrixXi& FO) {
        // TODO Check if is_edge_manifold
        MatrixXd newV = V;
        MatrixXi newF = F;

        const int max_iter = std::ceil(0.5*Q.size());
        int count = 0;
        if (Q.empty()) std::cout << "Q EMPTY" << std::endl;
          for(int j = 0;j<max_iter;j++)
          {
            if (Q.empty())
                break;

            if (Q.begin()->first == std::numeric_limits<double>::infinity())
                break; // min cost edge is infinite cost


            if (!igl::collapse_edge(igl::shortest_edge_and_midpoint,newV,newF,E,EMAP,EF,EI,Q,Qit,C))
            {
                  break; // failed to remove
            } 
            ++count;
        }

        // https://github.com/libigl/libigl/blob/master/include/igl/decimate.cpp
        std::cout << "Num collapsed: " << count << std::endl;
        int m = 0;
        for(int f = 0;f<F.rows();f++)
        {
            if(newF(f,0) != IGL_COLLAPSE_EDGE_NULL || 
               newF(f,1) != IGL_COLLAPSE_EDGE_NULL || 
               newF(f,2) != IGL_COLLAPSE_EDGE_NULL)
            {
                if (m != f)
                    newF.row(m) = newF.row(f);
                m++;
            }
        }
        newF.conservativeResize(m,3);

        VectorXi I; // I maps from Old to New Vertices
        MatrixXd NV;
        MatrixXi NF;
        igl::remove_unreferenced(newV, newF, NV, NF, I);
        std::cout << " I size: " << I.size() << std::endl;
        int c = 0;
        for (int i = 0; i < I.size(); i++) {
            if (I[i] == -1) c++;
        }
        std::cout << "invalid cnt: " << c << std::endl;

        VO = NV;
        FO = NF;

        // course and fine matrix sizes
        const int fsize = V.rows();
        const int csize = NV.rows();
        SparseMatrix<double, RowMajor> p(fsize, csize);
        std::vector<Triplet<double>> triplets;
        // Form Two blocks of prolongation operator -- an Identity block coarse x coarse
        // and a mapping from fine -> coarse vertices
        for(int i = 0; i < csize; ++i)
            triplets.push_back(Triplet(i,i,1);

        // Make map 
        // 1. Adjacency list for fine vertices
        // 2. loop through invalid vertices
        // 3. for-each vertex, find neighbors
        // 4. Count the number of neighbors for prolongation value 1/R
        // 5. map each fine vertex (that is removed) to a coarse, using I(fine) -> coarse(Index)
        // 6. Add triplet
        

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
