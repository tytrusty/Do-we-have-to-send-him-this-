#include "HeatHook.h"
#include <igl/unproject_onto_mesh.h>
#include <igl/cotmatrix.h>
#include "igl/opengl/glfw/imgui/ImGuiHelpers.h"
#include <map>
#include <omp.h>
#include "Solver.h"

using namespace Eigen;

HeatHook::HeatHook() : PhysicsHook() 
{
    clickedVertex = -1;
    dt = 1e0;
    mcf_dt = 5e-4;
    meshFile_ = "rect-coarse.obj";
    solverIters = 40;
    solverTol = 1e-7;
}

void HeatHook::drawGUI(igl::opengl::glfw::imgui::ImGuiMenu &menu)
{
	if (ImGui::CollapsingHeader("Mesh", ImGuiTreeNodeFlags_DefaultOpen))
    {
        ImGui::InputText("Filename", meshFile_);
	}
	if (ImGui::CollapsingHeader("Simulation Options", ImGuiTreeNodeFlags_DefaultOpen))
    {
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

//bool HeatHook::mouseMoved(igl::opengl::glfw::Viewer &viewer, int button)
//{
//    MouseEvent me;
//    me.type = MouseEvent::ME_DRAGGED;    
//    double x = viewer.current_mouse_x;
//    double y = viewer.core.viewport(3) - viewer.current_mouse_y;
//    Eigen::Vector3d pos(x, y, clickedz);
//    igl::unproject(pos, viewer.core.view * viewer.core.model,
//        viewer.core.proj, viewer.core.viewport, me.pos);
//    mouseMutex.lock();
//    mouseEvents.push_back(me);
//    mouseMutex.unlock();
//    return false;
//}

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

	std::string meshfname = std::string("meshes/") + meshFile_;
    if(!igl::readOBJ(meshfname, V, F))
        meshfname = std::string("../meshes/") + meshFile_;
        if (!igl::readOBJ(meshfname, V, F))
        {
            std::cerr << "Couldn't read mesh file" << std::endl;
            exit(-1);
        }
    // V *= 40.0; 
    // V /= 10.0;
    prevClicked = -1;

    double start;
    start = omp_get_wtime();
    igl::cotmatrix(V, F, L);
    igl::massmatrix(V, F, igl::MASSMATRIX_TYPE_DEFAULT, Morig);
    V /= std::sqrt(Morig.diagonal().sum());
    igl::massmatrix(V, F, igl::MASSMATRIX_TYPE_DEFAULT, Morig);
    M = Morig;
    Vdot = MatrixXd::Zero(V.rows(), V.cols());
    std::cout << "cot & mass matrix time (s): " << omp_get_wtime() - start << std::endl;
    source = VectorXd::Zero(V.rows());
    phi    = VectorXd::Zero(V.rows());
}

bool HeatHook::simulateOneStep()
{
    //V += mcf_dt * Vdot;
    //std::cout << std::sqrt(Morig.diagonal().prod() * 1./M.diagonal().prod()) << std::endl;
    //Vdot = -std::sqrt(Morig.diagonal().prod() * 1./M.diagonal().prod()) * L * V;

    SimplicialLDLT<SparseMatrix<double>> solver;
    double prefactor = std::sqrt(Morig.diagonal().cwiseQuotient(M.diagonal()).prod());
    //std::cout << Morig.diagonal() << std::endl;
    std::cout << prefactor << std::endl;
    solver.compute(M - prefactor*mcf_dt*L);
    V = solver.solve(M * V);

    V /= std::sqrt(M.diagonal().sum()); // rescale
    igl::massmatrix(V, F, igl::MASSMATRIX_TYPE_DEFAULT, M);
    // igl::cotmatrix(V, F, L);

    if (clickedVertex != -1 && clickedVertex != prevClicked)
    {
        // Update params
        Solver::updateIters(solverIters);
        Solver::updateTolerance(solverTol);

        source[clickedVertex] = 1.0;
        std::cout << "clicked: " << clickedVertex << std::endl;
        int nfaces = F.rows();
        MatrixXd ugrad(nfaces, 3);
        integrateHeat(ugrad);
        solveDistance(ugrad);
        prevClicked = clickedVertex;
    }
    return false;
}
