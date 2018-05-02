#ifndef HEATHOOK_H
#define HEATHOOK_H
#include "PhysicsHook.h"
#include <igl/readOBJ.h>
#include <iostream>
#include <igl/isolines.h>
#include <igl/edges.h>
#include <igl/components.h>
#include <Eigen/Core>

typedef std::tuple<int, int> Edge;

struct MouseEvent
{
    enum METype {
        ME_CLICKED,
        ME_RELEASED,
        ME_DRAGGED
    };

    METype type;
    int vertex;
    Eigen::Vector3d pos;
};

class HeatHook : public PhysicsHook
{
public:
    HeatHook();

    virtual void drawGUI(igl::opengl::glfw::imgui::ImGuiMenu &menu);
    
    virtual void tick();

    virtual void initSimulation();

    virtual void updateRenderGeometry()
    {
        renderV = V;
        renderF = F;        
    }

    virtual bool simulateOneStep();
    
    virtual void renderRenderGeometry(igl::opengl::glfw::Viewer &viewer)
    {
        Eigen::MatrixXd isoV;
        Eigen::MatrixXi isoE;
        igl::isolines(V, F, phi, 30, isoV, isoE);

        viewer.data().clear();
        viewer.data().set_mesh(renderV, renderF);
        igl::jet(phi, true, C);
        viewer.data().set_colors(C);
        viewer.data().set_edges(isoV, isoE, Eigen::RowVector3d(0.,0.,0.));

    }

    virtual bool mouseClicked(igl::opengl::glfw::Viewer &viewer, int button);
    virtual bool mouseReleased(igl::opengl::glfw::Viewer &viewer,  int button);
    // virtual bool mouseMoved(igl::opengl::glfw::Viewer &viewer,  int button);
    
private:
    void integrateHeat(Eigen::MatrixXd&);
    void solveDistance(const Eigen::MatrixXd& ugrad);

    Eigen::MatrixXd V;
    Eigen::MatrixXi F;
    Eigen::VectorXd phi;
    std::string meshFile_;


    std::mutex mouseMutex;
    std::vector<MouseEvent> mouseEvents;
    int clickedVertex; // the currently selected vertex (-1 if no vertex)
    int prevClicked;
    double clickedz;
    Eigen::Vector3d curPos; // the current position of the mouse cursor in 3D

    float dt;
    double mcf_dt;
    int solverIters;
    float solverTol;

    Eigen::SparseMatrix<double> L;
    Eigen::SparseMatrix<double> M;
    Eigen::SparseMatrix<double> Morig;
    Eigen::VectorXd u;
    Eigen::VectorXd source;
    Eigen::MatrixXd C;
    Eigen::MatrixXd Vdot;

    Eigen::MatrixXd renderV;
    Eigen::MatrixXi renderF;
};

#endif // HEATHOOK_H
