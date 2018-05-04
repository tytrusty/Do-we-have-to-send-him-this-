#include "HeatHook.h"
#include <igl/unproject_onto_mesh.h>
#include <igl/cotmatrix.h>
#include "igl/opengl/glfw/imgui/ImGuiHelpers.h"
#include <map>
#include <omp.h>
#include "Solver.h"
#include <string>
#include <igl/readPLY.h>
#include "Cluster.h"
using namespace Eigen;

HeatHook::HeatHook() : PhysicsHook() 
{
    clickedVertex = -1;
    dt = 1e0;
    heat_dt = 1e-5;
    mcf_dt = 5e-4;
    mass_fixed = false;
    cotan_fixed = true;
    explicit_mcf = false;
    meshFile_ = "cube.obj";
    solverIters = 40;
    solverTol = 1e-7;
}

double HeatHook::computeVolume() {
    double volume_ = 0;
    for (int i = 0; i < F.rows(); i++)
    {
        Vector3d pts[3];
        Vector3d centroid(0, 0, 0);
        for (int j = 0; j < 3; j++)
        {
            pts[j] = V.row(F(i, j));
            centroid += pts[j];
        }
        Vector3d normal = (pts[1] - pts[0]).cross(pts[2] - pts[0]);
        double area = 0.5 * normal.norm();
        normal /= normal.norm();

        centroid /= 3.0;
        volume_ += centroid.dot(normal) * area / 3.0;
    }
	return volume_;
}

Vector3d HeatHook::computeCenterOfMass(double volume_) {
    Vector3d cm(0, 0, 0);
    for (int i = 0; i < F.rows(); i++)
    {
        Vector3d pts[3];
        for (int j = 0; j < 3; j++)
        {
            pts[j] = V.row(F(i, j));
        }
        Vector3d normal = (pts[1] - pts[0]).cross(pts[2] - pts[0]);
        double area = 0.5 * normal.norm();
        normal /= normal.norm();

        Vector3d term(0, 0, 0);
        for (int j = 0; j < 3; j++)
        {
            for (int k = 0; k < 3; k++)
            {
                for (int l = k; l < 3; l++)
                    term[j] += pts[k][j] * pts[l][j];
            }
            term[j] *= area*normal[j] / 12.0;
        }

        cm += term;
    }

    return cm / volume_;
}

void HeatHook::drawGUI(igl::opengl::glfw::imgui::ImGuiMenu &menu)
{
	if (ImGui::CollapsingHeader("Mesh", ImGuiTreeNodeFlags_DefaultOpen))
    {
        ImGui::InputText("Filename", meshFile_);
	}
	if (ImGui::CollapsingHeader("Simulation Options", ImGuiTreeNodeFlags_DefaultOpen))
    {
        ImGui::InputFloat("MCF timestep", &mcf_dt, 0, 0, 5);
        ImGui::Checkbox("Mass fixed", &mass_fixed);
        ImGui::Checkbox("Cotan fixed", &cotan_fixed);
        ImGui::Checkbox("MCF is explicit", &explicit_mcf);
        ImGui::InputFloat("Heat flow timestep", &heat_dt, 0, 0, 5);
        ImGui::InputFloat("Timestep", &dt, 0, 0, 3);
        ImGui::InputInt("Solver Iters", &solverIters);
        ImGui::InputFloat("Solver Tolerance", &solverTol, 0, 0, 12);
	}
}

bool HeatHook::mouseClicked(igl::opengl::glfw::Viewer &viewer, int button)
{
    if(button != 0)
        return false;
    render_mutex.lock();
    double x = viewer.current_mouse_x;
    double y = viewer.core.viewport(3) - viewer.current_mouse_y;
    int fid;
    Eigen::Vector3f bc;
    MouseEvent me;
    bool ret = true;
    if (igl::unproject_onto_mesh(Eigen::Vector2f(x, y), viewer.core.view * viewer.core.model,
        viewer.core.proj, viewer.core.viewport, V, F, fid, bc))
    {
        int bestvert = -1;
        double bestcoord = 2.0;
        for (int j = 0; j < 3; j++)
        {
            if (bc[j] < bestcoord)
            {
                bestcoord = bc[j];
                bestvert = j;
            }
        }
        me.type = MouseEvent::ME_CLICKED;
        me.vertex = F(fid, bestvert);        
        
        Eigen::Vector3f proj;
        Eigen::Vector3f pt;
        for (int i = 0; i < 3; i++)
            pt[i] = V(me.vertex, i);
        Eigen::Matrix4f modelview = viewer.core.view * viewer.core.model;
        proj = igl::project(pt, modelview,
            viewer.core.proj, viewer.core.viewport);

        clickedz = proj[2];
        Eigen::Vector3f pos;
        pos[0] = float(x);
        pos[1] = float(y);
        pos[2] = float(clickedz);
        Eigen::Vector3f unproj = igl::unproject(pos, modelview,
            viewer.core.proj, viewer.core.viewport);
        for (int i = 0; i < 3; i++)
            me.pos[i] = unproj[i];
        ret = true;
    }
    else
    {
        me.type = MouseEvent::ME_RELEASED;
        ret = false;
    }       
    render_mutex.unlock();

    mouseMutex.lock();
    mouseEvents.push_back(me);
    mouseMutex.unlock();
    return ret;
}

bool HeatHook::mouseReleased(igl::opengl::glfw::Viewer &viewer, int button)
{
    MouseEvent me;
    me.type = MouseEvent::ME_RELEASED;
    mouseMutex.lock();
    mouseEvents.push_back(me);
    mouseMutex.unlock();
    return false;
}

bool HeatHook::mouseMoved(igl::opengl::glfw::Viewer &viewer, int button)
{
    MouseEvent me;
    me.type = MouseEvent::ME_DRAGGED;    
    double x = viewer.current_mouse_x;
    double y = viewer.core.viewport(3) - viewer.current_mouse_y;
    int fid;
    Eigen::Vector3f bc;
    if (igl::unproject_onto_mesh(Eigen::Vector2f(x, y), viewer.core.view * viewer.core.model,
        viewer.core.proj, viewer.core.viewport, V, F, fid, bc))
    {
        int bestvert = -1;
        double bestcoord = 2.0;
        for (int j = 0; j < 3; j++)
        {
            if (bc[j] < bestcoord)
            {
                bestcoord = bc[j];
                bestvert = j;
            }
        }
        me.type = MouseEvent::ME_CLICKED;
        me.vertex = F(fid, bestvert);        
        
        Eigen::Vector3f proj;
        Eigen::Vector3f pt;
        for (int i = 0; i < 3; i++)
            pt[i] = V(me.vertex, i);
        Eigen::Matrix4f modelview = viewer.core.view * viewer.core.model;
        proj = igl::project(pt, modelview,
            viewer.core.proj, viewer.core.viewport);

        clickedz = proj[2];
        Eigen::Vector3f pos;
        pos[0] = float(x);
        pos[1] = float(y);
        pos[2] = float(clickedz);
        Eigen::Vector3f unproj = igl::unproject(pos, modelview,
            viewer.core.proj, viewer.core.viewport);
        for (int i = 0; i < 3; i++)
            me.pos[i] = unproj[i];
    }

    mouseMutex.lock();
    mouseEvents.push_back(me);
    mouseMutex.unlock();
    return false;
}

void HeatHook::tick()
{
    mouseMutex.lock();
    for (MouseEvent me : mouseEvents)
    {
        if (me.type == MouseEvent::ME_CLICKED)
        {            
            curPos = me.pos;
            clickedVertex = me.vertex;         
        }
        if (me.type == MouseEvent::ME_RELEASED)
        {
            clickedVertex = -1;
        }
        if (me.type == MouseEvent::ME_DRAGGED)
        {
            curPos = me.pos;
        }
    }
    mouseEvents.clear();
    mouseMutex.unlock();
}

void HeatHook::buildCotanLaplacian(Eigen::SparseMatrix<double>& L_)
{
    std::cout << "L size " << L_.size() << std::endl;
//    igl::cotmatrix(V, F, L_);
//    std::cout << "igl\n" << L_ << std::endl << std::endl;
    L_.setZero();
    std::vector<Triplet<double>> trips;
    int nfaces = F.rows();
    for (int i = 0; i < nfaces; i++) {
        Vector3i face = F.row(i);
        for (int j = 0; j < 3; j++) {
            Vector3d e1 = V.row(face[(j+1)%3]) - V.row(face[j]);
            Vector3d e2 = V.row(face[(j+2)%3]) - V.row(face[j]);
            Vector3d e3 = e2 - e1;

            double cot1 = e2.dot(e3)/(-e2).cross(-e3).norm();
            double cot2 = (-e1).dot(e3)/(e1).cross(-e3).norm();

            trips.push_back(Triplet<double>(face[j], face[(j+1)%3], .5*(true ? cot1 : (cot1 + cot2))));
            trips.push_back(Triplet<double>(face[(j+1)%3], face[j], .5*(true ? cot1 : (cot1 + cot2))));
            trips.push_back(Triplet<double>(face[j], face[j], -.5*(true? cot1 : (cot1 + cot2))));
            trips.push_back(Triplet<double>(face[(j+1)%3], face[(j+1)%3], -.5*(true? cot1 : (cot1 + cot2))));
        }
    }
    L_.setFromTriplets(trips.begin(), trips.end());
//    std::cout << L_ << std::endl;
}
void HeatHook::integrateHeat(MatrixXd& ugrad)
{
    double start;
    VectorXd u = VectorXd(V.rows());

    start = omp_get_wtime();
    // solve heat flow
    SparseMatrix<double> A = M - dt*L;
    ConjugateGradient<SparseMatrix<double>, Lower|Upper> cg;
    cg.compute(A);
    u = cg.solve(source);
    // Solver::gauss_seidel(A, source, u);
    std::cout << "solve heat time (s): " << omp_get_wtime() - start << std::endl;

    start = omp_get_wtime();
    // evaluate flow field
    int nfaces = F.rows();
    ugrad.setZero();
    // get u gradient
    #pragma omp parallel for 
    for (int i = 0; i < nfaces; i++) { 
        Vector3i face = F.row(i);
        Vector3d e1 = V.row(face[1]) - V.row(face[0]);
        Vector3d e2 = V.row(face[2]) - V.row(face[0]);
        Vector3d normal = e1.cross(e2);
        for (int j = 0; j < 3; j++) {
            Vector3d oppEdge = V.row(face[(j+2)%3]) - V.row(face[(j+1)%3]);
            ugrad.row(i) += u[face[j]]*(normal.cross(oppEdge));
        }
        ugrad.row(i) /= (e1.cross(e2)).norm(); // divide by area of face
        ugrad.row(i) /= -ugrad.row(i).norm();  // normalize and flip direction
    }
    std::cout << "compute ugrad time (s): " << omp_get_wtime() - start << std::endl;
}

void HeatHook::solveDistance(const MatrixXd& ugrad) 
{
    // divergence
    double start = omp_get_wtime();
    VectorXd div = VectorXd::Zero(V.rows());
    phi = VectorXd::Zero(V.rows());
    int nfaces = F.rows();
    for (int i = 0; i < nfaces; i++) {
        Vector3i face = F.row(i);
        for (int j = 0; j < 3; j++) {
            Vector3d e1 = V.row(face[(j+1)%3]) - V.row(face[j]);
            Vector3d e2 = V.row(face[(j+2)%3]) - V.row(face[j]);
            Vector3d e3 = e2 - e1;
            double cot1 = e2.dot(e3)/(-e2).cross(-e3).norm();
            double cot2 = (-e1).dot(e3)/(-e1).cross(e3).norm();
            div[face[j]] += 0.5*(cot1*(e1.dot(ugrad.row(i))) + cot2*(e2.dot(ugrad.row(i))));
        }
    }
    std::cout << "div computation time (s): " << omp_get_wtime() - start << std::endl;
    // Solve for distance
    start = omp_get_wtime();
    Solver::gauss_seidel(L, div, phi);
    // Finds a different solution ???
    //ConjugateGradient<SparseMatrix<double>, Lower|Upper> cg;
    //cg.compute(L);
    //phi = cg.solve(div);
    std::cout << "solve dist time (s): " << omp_get_wtime() - start << std::endl;
}


void HeatHook::initSimulation()
{
    Eigen::initParallel();
    std::cout << "Num threads: " << Eigen::nbThreads() << std::endl;
    #pragma omp parallel
    {
        std::cout << "Thread num: " << omp_get_thread_num() << std::endl;
    }

	std::string meshfname = std::string("../meshes/") + meshFile_;
    std::string ext = meshFile_.substr(meshFile_.find_last_of(".")+1);
    if (ext == "obj") {
        std::cout << "Reading OBJ file" << std::endl;
        if(!igl::readOBJ(meshfname, V, F))
        {
            std::cerr << "Couldn't read mesh file" << std::endl;
            exit(-1);
        }
    } else if (ext == "ply") {
        std::cout << "Reading PLY file" << std::endl;
        if(!igl::readPLY(meshfname, V, F))
        {
            std::cerr << "Couldn't read mesh file" << std::endl;
            exit(-1);
        }
        std::cout << "Done reading PLY" << std::endl;

    }

    V *= 10.0; 
    // V /= 10.0;
    prevClicked = -1;

    L = SparseMatrix<double>(V.rows(), V.rows());
    double start;
    start = omp_get_wtime();
    //igl::cotmatrix(V, F, L);
    buildCotanLaplacian(L);

    igl::massmatrix(V, F, igl::MASSMATRIX_TYPE_DEFAULT, Morig);
    V /= std::sqrt(Morig.diagonal().sum());
    igl::massmatrix(V, F, igl::MASSMATRIX_TYPE_DEFAULT, Morig);
    igl::massmatrix(V, F, igl::MASSMATRIX_TYPE_DEFAULT, M);
    std::cout << "cot & mass matrix time (s): " << omp_get_wtime() - start << std::endl;
    source = VectorXd::Zero(V.rows());
    phi    = VectorXd::Zero(V.rows());

    double volume = computeVolume();
    std::cout << "vol: " << volume << std::endl;
    Vector3d cm = computeCenterOfMass(volume);
    std::cout << "CM: " << cm << std::endl;
    for (int i = 0; i < V.rows(); i++)
        V.row(i) -= cm;

//    Cluster cluster(F, 69, V.rows());
//    cluster.BFS();
//    int n = 0;
//    for(unordered_set<int> set : cluster.cluster_map_true) {
//        for (int i : set) {
//            phi[i] = n/4.0;
//        }
//        ++n;
//    }
}

bool HeatHook::simulateOneStep()
{

    if (true) // curvature flow 
    {
        SparseMatrix<double>& useM = mass_fixed ? M : Morig;
        if (explicit_mcf) {
            SparseMatrix<double> Minv(V.rows(), V.rows());
            std::vector<Triplet<double>> trips;
            VectorXd diag = M.diagonal().cwiseQuotient(M.diagonal().cwiseProduct(M.diagonal()));
            int nverts = V.rows();
            for (int i = 0; i < nverts; i++) {
                trips.push_back(Triplet<double>(i, i, diag[i]));
            }
            Minv.setFromTriplets(trips.begin(), trips.end());
            V +=  mcf_dt * Minv * L  * V;
        } else {
            SimplicialLDLT<SparseMatrix<double>> solver;
            solver.compute(useM - mcf_dt*L);
            V = solver.solve(useM * V);
        }
        // normalize area

        if (!mass_fixed) {
            igl::massmatrix(V, F, igl::MASSMATRIX_TYPE_DEFAULT, Morig);
            V /= std::sqrt(Morig.diagonal().sum());
        }
        else {
            igl::massmatrix(V, F, igl::MASSMATRIX_TYPE_DEFAULT, M);
            V /= std::sqrt(M.diagonal().sum());
        }

        if (!cotan_fixed)
            buildCotanLaplacian(L);
            //igl::cotmatrix(V, F, L);

        if (false) {} // if non-conformal, update the cotan matrix

        // account for new c.o.m
        double volume = computeVolume();
        Vector3d cm = computeCenterOfMass(volume);
        for (int i = 0; i < V.rows(); i++)
            V.row(i) -= cm;
    }

    if (false) // heat flow
    {
        VectorXd u = VectorXd(V.rows());
        SparseMatrix<double> A = M - heat_dt*L;
        SimplicialLDLT<SparseMatrix<double>> solver;
        solver.compute(A);
        u = solver.solve(M*source);
        phi = u;
        source = u;
    }

    if (clickedVertex != -1 && clickedVertex != prevClicked)
    {
        if (false) // calc geodesics is checked 
        {
            // Update params
            Solver::updateIters(solverIters);
            Solver::updateTolerance(solverTol);
            source[clickedVertex] = 1.0;
            //std::cout << "clicked: " << clickedVertex << std::endl;
            int nfaces = F.rows();
            MatrixXd ugrad(nfaces, 3);
            integrateHeat(ugrad);
            solveDistance(ugrad);
            prevClicked = clickedVertex;
        } else if (true) { // show heat flow
            //std::cout << "clicked: " << clickedVertex << std::endl;
            source[clickedVertex] = 1.0;
            prevClicked = clickedVertex;
        }
    }
    return false;
}
