#ifndef SOLVER_H
#define SOLVER_H #include <omp.h>
#include <algorithm> 
#include <Eigen/Core>
#include <igl/decimate.h>

using namespace Eigen;
using std::vector;

class Solver
{
public:
    Solver()
    {

    }
    // const int preIters=3, postIters=2, coarseIters=10;
    const int preIters=30, postIters=20, coarseIters=100;
    int iters_ = 400;
    float tolerance = 1e-7;
    std::map<int,vector<int> > fineMap;

	bool decimate(
	  const Eigen::MatrixXd & OV,
	  const Eigen::MatrixXi & OF,
	  const std::function<void(
	    const int,
	    const Eigen::MatrixXd &,
	    const Eigen::MatrixXi &,
	    const Eigen::MatrixXi &,
	    const Eigen::VectorXi &,
	    const Eigen::MatrixXi &,
	    const Eigen::MatrixXi &,
	    double &,
	    Eigen::RowVectorXd &)> & cost_and_placement,
	  const std::function<bool(
	      const Eigen::MatrixXd &,
	      const Eigen::MatrixXi &,
	      const Eigen::MatrixXi &,
	      const Eigen::VectorXi &,
	      const Eigen::MatrixXi &,
	      const Eigen::MatrixXi &,
	      const std::set<std::pair<double,int> > &,
	      const std::vector<std::set<std::pair<double,int> >::iterator > &,
	      const Eigen::MatrixXd &,
	      const int,
	      const int,
	      const int,
	      const int,
	      const int)> & stopping_condition,
	    const std::function<bool(
	      const Eigen::MatrixXd &                                         ,/*V*/
	      const Eigen::MatrixXi &                                         ,/*F*/
	      const Eigen::MatrixXi &                                         ,/*E*/
	      const Eigen::VectorXi &                                         ,/*EMAP*/
	      const Eigen::MatrixXi &                                         ,/*EF*/
	      const Eigen::MatrixXi &                                         ,/*EI*/
	      const std::set<std::pair<double,int> > &                        ,/*Q*/
	      const std::vector<std::set<std::pair<double,int> >::iterator > &,/*Qit*/
	      const Eigen::MatrixXd &                                         ,/*C*/
	      const int                                                        /*e*/
	      )> & pre_collapse,
	    const std::function<void(
	      const Eigen::MatrixXd &                                         ,   /*V*/
	      const Eigen::MatrixXi &                                         ,   /*F*/
	      const Eigen::MatrixXi &                                         ,   /*E*/
	      const Eigen::VectorXi &                                         ,/*EMAP*/
	      const Eigen::MatrixXi &                                         ,  /*EF*/
	      const Eigen::MatrixXi &                                         ,  /*EI*/
	      const std::set<std::pair<double,int> > &                        ,   /*Q*/
	      const std::vector<std::set<std::pair<double,int> >::iterator > &, /*Qit*/
	      const Eigen::MatrixXd &                                         ,   /*C*/
	      const int                                                       ,   /*e*/
	      const int                                                       ,  /*e1*/
	      const int                                                       ,  /*e2*/
	      const int                                                       ,  /*f1*/
	      const int                                                       ,  /*f2*/
	      const bool                                                  /*collapsed*/
	      )> & post_collapse,
	  Eigen::MatrixXd & U,
	  Eigen::MatrixXi & G,
	  Eigen::VectorXi & J,
	  Eigen::VectorXi & I
	  )
	{
	  // Working copies
	  Eigen::MatrixXd V = OV;
	  Eigen::MatrixXi F = OF;
      VectorXi EMAP;
      MatrixXi E,EF,EI;
      igl::edge_flaps(OF,E,EMAP,EF,EI);
      MatrixXi Ecopy = E;
	  typedef std::set<std::pair<double,int> > PriorityQueue;
	  PriorityQueue Q;
	  std::vector<PriorityQueue::iterator > Qit;
	  Qit.resize(E.rows());
	  // If an edge were collapsed, we'd collapse it to these points:
	  MatrixXd C(E.rows(),V.cols());
	  for(int e = 0;e<E.rows();e++)
	  {
	    double cost = e;
	    RowVectorXd p(1,3);
	    cost_and_placement(e,V,F,E,EMAP,EF,EI,cost,p);
	    C.row(e) = p;
	    Qit[e] = Q.insert(std::pair<double,int>(cost,e)).first;
	  }
	  int prev_e = -1;
	  bool clean_finish = false;
	
	  while(true)
	  {
	    if(Q.empty())
	    {
	      break;
	    }
	    if(Q.begin()->first == std::numeric_limits<double>::infinity())
	    {
	      // min cost edge is infinite cost
	      break;
	    }
	    int e,e1,e2,f1,f2;
	    if(igl::collapse_edge(
	       cost_and_placement, pre_collapse, post_collapse,
	       V,F,E,EMAP,EF,EI,Q,Qit,C,e,e1,e2,f1,f2))
	    {
            const int eflip = Ecopy(e,0)>Ecopy(e,1);
            // source and destination
            const int s = eflip?Ecopy(e,1):Ecopy(e,0); 
            const int d = eflip?Ecopy(e,0):Ecopy(e,1); 
            std::cout << "fineMap push back " << s << " d : " << d << std::endl;
/* DOESNT HANDLE ALL EDGE CASES
fineMap push back 11 d : 14
fineMap push back 5 d : 15
fineMap push back 0 d : 12
fineMap push back 0 d : 10
fineMap push back 2 d : 13
fineMap push back 1 d : 3
fineMap push back 1 d : 4
fineMap push back 0 d : 2
fineMap push back 8 d : 15
fineMap push back 10 d : 16
*/
            if (fineMap.find(d) != fineMap.end())
                fineMap[s].insert(fineMap[s].end(), fineMap[d].begin(), fineMap[d].end());
            fineMap[s].push_back(d);
            if(stopping_condition(V,F,E,EMAP,EF,EI,Q,Qit,C,e,e1,e2,f1,f2))
	        {
	            clean_finish = true;
	            break;
	        }
	    }else
	    {
	      if(prev_e == e)
	      {
	        assert(false && "Edge collapse no progress... bad stopping condition?");
	        break;
	      }
	      // Edge was not collapsed... must have been invalid. collapse_edge should
	      // have updated its cost to inf... continue
	    }
	    prev_e = e;
	  }
	  // remove all IGL_COLLAPSE_EDGE_NULL faces
	  MatrixXi F2(F.rows(),3);
	  J.resize(F.rows());
	  int m = 0;
	  for(int f = 0;f<F.rows();f++)
	  {
	    if(
	      F(f,0) != IGL_COLLAPSE_EDGE_NULL || 
	      F(f,1) != IGL_COLLAPSE_EDGE_NULL || 
	      F(f,2) != IGL_COLLAPSE_EDGE_NULL)
	    {
	      F2.row(m) = F.row(f);
	      J(m) = f;
	      m++;
	    }
	  }
	  F2.conservativeResize(m,F2.cols());
	  J.conservativeResize(m);
	  VectorXi _1;
      igl::remove_unreferenced(V,F2,U,G,I,_1);
	  return clean_finish;
	}

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
        const auto always_try = [](
          const Eigen::MatrixXd &, const Eigen::MatrixXi &, const Eigen::MatrixXi &,
          const Eigen::VectorXi &, const Eigen::MatrixXi &, const Eigen::MatrixXi &,
          const std::set<std::pair<double,int> > &                        ,
          const std::vector<std::set<std::pair<double,int> >::iterator > &,
          const Eigen::MatrixXd &, const int) -> bool { return true;};
        const auto never_care = [](
          const Eigen::MatrixXd &, const Eigen::MatrixXi &, const Eigen::MatrixXi &,
          const Eigen::VectorXi &, const Eigen::MatrixXi &, const Eigen::MatrixXi &,
          const std::set<std::pair<double,int> > &,
          const std::vector<std::set<std::pair<double,int> >::iterator > &,
          const Eigen::MatrixXd &, const int, const int, const int, const int, 
          const int, const bool )-> void { };
        bool ret = decimate(
          VO,
          FO,
          igl::shortest_edge_and_midpoint,
          igl::max_faces_stopping_condition(m,orig_m,max_m),
          always_try,
          never_care,
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
        //std::cout << " I: " << I << std::endl;
        //std::cout << "G: \n" << G << std::endl;
        //std::cout << "F: \n" << F << std::endl;
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
        }
        return error;
    }

    std::vector<SparseMatrix<double>,Eigen::aligned_allocator<SparseMatrix<double>> > P;
    std::vector<SparseMatrix<double>,Eigen::aligned_allocator<SparseMatrix<double>> > R;
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
            const int csize = coarseV.rows(); // don't account for infinite vertex
            SparseMatrix<double, RowMajor> r(csize, fsize);
            SparseMatrix<double, RowMajor> p(fsize, csize);
            std::vector<Triplet<double>> rtriplets;
            std::vector<Triplet<double>> ptriplets;

            for (int v = 0; v < I.size() - 1; ++v) {
                int coarsev = I[v];
                // std::cout << "v: " << v << " coarsev: " << coarsev << std::endl;
                if (coarsev != -1) {
                    if (fineMap.find(v) != fineMap.end()) {
                        double avg = 1.0 / (1 + fineMap[v].size()); 
                        for (int i = 0; i < fineMap[v].size(); i++) {
                            rtriplets.push_back(Triplet<double>(coarsev,fineMap[v][i],avg));
                            ptriplets.push_back(Triplet<double>(fineMap[v][i], coarsev, 1)); // 
                        }
                        rtriplets.push_back(Triplet<double>(coarsev,v,avg));
                    } else {
                        rtriplets.push_back(Triplet<double>(coarsev,v,1));
                    }
                    ptriplets.push_back(Triplet<double>(v,coarsev,1));
                }

            }
            p.setFromTriplets(ptriplets.begin(), ptriplets.end());
            P.push_back(p);
            r.setFromTriplets(rtriplets.begin(), rtriplets.end());
            R.push_back(r);
            if (csize == 13) {
                std::cout << "I:  \n" << I << std::endl;
                std::cout << " p : \n" << MatrixXd(p) << std::endl;
                std::cout << " r : \n" << MatrixXd(r) << std::endl;

            }
            std::cout << " coarseV.size() : " << coarseV.rows() << std::endl;
            fineV = coarseV;
            fineF = coarseF;
            fineMap.clear();
        }
        max_depth = P.size();
    }

    void multigrid(const SparseMatrix<double>& A, const VectorXd& b, VectorXd& x, int depth = 0)
    {
        if (depth == max_depth) {
            std::cout << "End condition " << std::endl;
            double error = gauss_seidel(A, b, x, coarseIters); 
            std::cout << "size: " << x.size() << " error: " <<  error << std::endl;
            std::cout << "A dim : " << A.rows() << " x " << A.cols() << std::endl;
            std::cout << "A : \n" << A << std::endl;
            std::cout << "r: \n " << b << std::endl;
            std::cout << "e: \n " << x << std::endl;
        } else {
            // Presmoothing
            double preerror = gauss_seidel(A, b, x, preIters); 
            std::cout << "size: " << x.size() << " preerror: " << preerror << std::endl;

            // Restriction
            VectorXd r = R[depth] * (b - A*x);
            std::cout << "R dim : " << R[depth].rows() << " x " << R[depth].cols() << std::endl;
            std::cout << "P dim : " << P[depth].rows() << " x " << P[depth].cols() << std::endl;
            std::cout << "A dim : " << A.rows() << " x " << A.cols() << std::endl;
            //std::cout << "A : \n" << A << std::endl;
            //std::cout << "r: \n " <<  r << std::endl;
            SparseMatrix<double> restrictA = R[depth] * A * P[depth];

            VectorXd cx(r.size());
            multigrid(restrictA, r, cx, depth+1); 
            if (depth == max_depth-1) {
                std::cout << "P Mat: \n" << P[depth] << std::endl;
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
