#include "HeatHook.h"
#include <limits>
#include <map>
#include <igl/unproject_onto_mesh.h>
#include <limits>

using namespace Eigen;

HeatHook::HeatHook() : PhysicsHook() 
{
    clickedVertex = -1;

    dt = 1e-3;
    constraintIters = 5;

    gravityEnabled = true;
    gravityG = -9.8;

    pinEnabled = true;
    pinWeight = 1.0;

    stretchEnabled = true;
    stretchWeight = 0.5;

    bendingEnabled = true;
    bendingWeight = 0.5;

    pullingEnabled = true;
    pullingWeight = 0.5;
}

void HeatHook::drawGUI(igl::opengl::glfw::imgui::ImGuiMenu &menu)
{
    if (ImGui::CollapsingHeader("Simulation Options", ImGuiTreeNodeFlags_DefaultOpen))
    {
        ImGui::InputFloat("Timestep", &dt, 0, 0, 3);        
        ImGui::InputInt("Constraint Iters", &constraintIters);
    }
    if (ImGui::CollapsingHeader("Forces", ImGuiTreeNodeFlags_DefaultOpen))
    {
        ImGui::Checkbox("Gravity Enabled", &gravityEnabled);
        ImGui::InputFloat("Gravity G", &gravityG, 0, 0, 3);        
        ImGui::Checkbox("Pins Enabled", &pinEnabled);
        ImGui::InputFloat("Pin Weight", &pinWeight, 0, 0, 3);        
        ImGui::Checkbox("Stretching Enabled", &stretchEnabled);
        ImGui::InputFloat("Stretching Weight", &stretchWeight, 0, 0, 3);    
        ImGui::Checkbox("Bending Enabled", &bendingEnabled);
        ImGui::InputFloat("Bending Weight", &bendingWeight, 0, 0, 3);   
        ImGui::Checkbox("Pulling Enabled", &pullingEnabled);
        ImGui::InputFloat("Pulling Weight", &pullingWeight, 0, 0, 3);   
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
        viewer.core.proj, viewer.core.viewport, renderQ, renderF, fid, bc))
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
            pt[i] = renderQ(me.vertex, i);
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


void HeatHook::initSimulation()
{
    if(!igl::readOBJ("meshes/cube.obj", origQ, F))
        if (!igl::readOBJ("../meshes/cube.obj", origQ, F))
        {
            std::cerr << "Couldn't read mesh file" << std::endl;
            exit(-1);
        }
    //mesh is tiny for some reason
    origQ /= 10;
    Q     = origQ;  
    prevQ = origQ;
    Qdot.resize(Q.rows(), 3);
    Qdot.setZero();

    // precompute diamonds
    for (int i = 0; i < F.rows(); i++)
    {
        VectorXi face = F.row(i);
        for (int j = 0; j < 2; j++)
            for (int k = j + 1; k < 3; k++) {
                int first = std::min(face[j], face[k]);
                int second = std::max(face[j], face[k]);
                Edge e = std::make_tuple(first, second);
                std::map<Edge, std::vector<int>>::iterator got = diamonds.find(e);
                if (got == diamonds.end()) {
                    std::vector<int> vertices;
                    for (int l = 0; l < 3; l++)
                        vertices.push_back(face[l]);
                    std::pair<Edge, std::vector<int>> new_pair = std::make_pair(e, vertices);
                    diamonds.insert(new_pair);
                } else {
                    for (int l = 0; l < 3; l++)
                    {
                        if (std::get<0>(got->first) != face[l] and
                            std::get<1>(got->first) != face[l])
                            got->second.push_back(face[l]);
                    }
                }
            }
    }


    // get top corners
    double xmin = Q.col(0).minCoeff();
    double xmax = Q.col(0).maxCoeff(); 
    double ymax = Q.col(1).maxCoeff();
    for (int i = 0; i < Q.rows(); ++i) {
        if (Q.row(i).y()==ymax && (Q.row(i).x()==xmin || Q.row(i).x()==xmax)) {
            pinned_ids.push_back(i);
            pinned.push_back(Q.row(i));
        }
    }
}

bool HeatHook::simulateOneStep()
{
    prevQ = Q;
    Q += Qdot*dt;
    for (int i = 0; i < constraintIters; ++i) {
        if (pinEnabled) {
            Q.row(pinned_ids[0]) = pinWeight*pinned[0].transpose() + (1-pinWeight)*Q.row(pinned_ids[0]);
            Q.row(pinned_ids[1]) = pinWeight*pinned[1].transpose() + (1-pinWeight)*Q.row(pinned_ids[1]);
        }
        if (stretchEnabled) {
            for (int j = 0; j < F.rows(); j++)
            {
                Vector3i face = F.row(j);
                Vector3d c0 = 1./3 * (origQ.row(face[0]) + origQ.row(face[1]) + origQ.row(face[2]));
                Vector3d c = 1./3 * (Q.row(face[0]) + Q.row(face[1]) + Q.row(face[2]));
                Matrix3d A, B;
                for (int k = 0; k < 3; k++) {
                    A.col(k) = Q.row(face[k]).transpose() - c;
                    B.col(k) = origQ.row(face[k]).transpose() - c0;
                }
                Matrix3d ABT = A*B.transpose();
                JacobiSVD<Matrix3d> svd(ABT, ComputeFullU | ComputeFullV);
                MatrixXd R = svd.matrixU() * svd.matrixV().transpose();
                for (int k = 0; k < 3; k++) {
                    Vector3d qtransform = R*(origQ.row(face[k]).transpose() - c0) + c; 
                    Q.row(face[k]) = stretchWeight*qtransform.transpose() + (1-stretchWeight)*Q.row(face[k]);
                }
            }
        }

        if (bendingEnabled) {
            //std::cout << " wtf xd lol " << std::endl;
            for (auto const& pair : diamonds)
            {
                std::vector<int> vertices = pair.second;
                if (vertices.size() < 4)
                    continue;
                Vector3d c, c0; c.setZero(); c0.setZero();
                for (int i : vertices) {
                    c  += Q.row(i)/4.0;
                    c0 += origQ.row(i)/4.0;
                }

                MatrixXd A(3,4), B(3,4);
                for (int i = 0; i < vertices.size(); ++i) {
                    A.col(i) = Q.row(vertices[i]).transpose() - c;
                    B.col(i) = origQ.row(vertices[i]).transpose() - c0;
                }
                MatrixXd ABT = A*B.transpose();
                JacobiSVD<MatrixXd> svd(ABT, ComputeFullU | ComputeFullV);
                MatrixXd R = svd.matrixU()*svd.matrixV().transpose();
                for (int i: vertices) {
                    Vector3d qtransform = R*(origQ.row(i).transpose() - c0) + c; 
                    Q.row(i) = stretchWeight*qtransform.transpose() + (1-stretchWeight)*Q.row(i);
                }
            }
        }
        if (pullingEnabled and clickedVertex != -1) {
            Q.row(clickedVertex) = pullingWeight*curPos.transpose() + (1-pinWeight)*Q.row(clickedVertex);
        }
    }
    Qdot = ((Q - prevQ) / dt);
    if (gravityEnabled) {
        for (int i = 0; i < Qdot.rows(); ++i)
            Qdot.row(i).y() += dt*gravityG;
    }
    return false;
}
