#include "HeatHook.h"
#include <limits>
#include <map>
#include <igl/unproject_onto_mesh.h>
#include <igl/cotmatrix.h>
#include <limits>
#include "igl/opengl/glfw/imgui/ImGuiHelpers.h"

using namespace Eigen;

HeatHook::HeatHook() : PhysicsHook() 
{
    clickedVertex = -1;
    dt = 1e0;
    meshFile_ = "rect-coarse.obj";
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
        viewer.core.proj, viewer.core.viewport, renderV, renderF, fid, bc))
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
        me.vertex = renderF(fid, bestvert);        
        
        Eigen::Vector3f proj;
        Eigen::Vector3f pt;
        for (int i = 0; i < 3; i++)
            pt[i] = renderV(me.vertex, i);
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
    Eigen::Vector3d pos(x, y, clickedz);
    igl::unproject(pos, viewer.core.view * viewer.core.model,
        viewer.core.proj, viewer.core.viewport, me.pos);
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

void HeatHook::integrateHeat(MatrixXd& ugrad)
{
    VectorXd b = VectorXd::Zero(V.rows());
    b[0] = 1.0;
    // b[1] = 1.0;
    VectorXd u = VectorXd(V.rows());

    // solve heat flow
    SparseMatrix<double> A = M - dt*L;
    ConjugateGradient<SparseMatrix<double>, Lower|Upper> cg;
    cg.compute(A);
    u = cg.solve(b);

    // evaluate flow field
    int nfaces = F.rows();
    ugrad.setZero();
    // get u gradient
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
}


void HeatHook::initSimulation()
{
    Eigen::initParallel();
    std::cout << "Num threads: " << Eigen::nbThreads();
	std::string meshfname = std::string("meshes/") + meshFile_;
    if(!igl::readOBJ(meshfname, V, F))
        meshfname = std::string("../meshes/") + meshFile_;
        if (!igl::readOBJ(meshfname, V, F))
        {
            std::cerr << "Couldn't read mesh file" << std::endl;
            exit(-1);
        }
    // V *= 40.0; 
    igl::cotmatrix(V, F, L);
    igl::massmatrix(V, F, igl::MASSMATRIX_TYPE_DEFAULT, M);
    int nfaces = F.rows();
    MatrixXd ugrad(nfaces, 3);
    integrateHeat(ugrad);

    // divergence
    VectorXd div = VectorXd::Zero(V.rows());
    phi = VectorXd::Zero(V.rows());
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

    ConjugateGradient<SparseMatrix<double>, Lower|Upper> cg;
    cg.compute(L);
    phi = cg.solve(div);

}

bool HeatHook::simulateOneStep()
{
}
