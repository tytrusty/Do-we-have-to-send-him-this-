#ifndef HEATHOOK_H
#define HEATHOOK_H
#include "PhysicsHook.h"
#include "Solver.h"
#include <igl/readOBJ.h>
#include <iostream>
#include <igl/isolines.h>
#include <igl/edges.h>
#include <igl/components.h>
#include <igl/colormap.h>
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

        viewer.data().clear();
        viewer.data().set_mesh(renderV, renderF);
        //viewer.data().set_mesh(V, F);
        
        igl::colormap(igl::ColorMapType(render_color), phi, normalize_color, C);
        viewer.data().set_colors(C);
        if (enable_iso) 
        {
            Eigen::MatrixXd isoV;
            Eigen::MatrixXi isoE;
            igl::isolines(renderV, renderF, phi, 30, isoV, isoE);
            viewer.data().set_edges(isoV, isoE, Eigen::RowVector3d(0.,0.,0.));
        }

    }

    virtual bool mouseClicked(igl::opengl::glfw::Viewer &viewer, int button);
    virtual bool mouseReleased(igl::opengl::glfw::Viewer &viewer,  int button);
    virtual bool mouseMoved(igl::opengl::glfw::Viewer &viewer,  int button);
    
private:
    void buildCotanLaplacian(Eigen::SparseMatrix<double>&);
    void integrateHeat(Eigen::MatrixXd&);
    void solveDistance(const Eigen::MatrixXd& ugrad);
    double computeVolume();
    Eigen::Vector3d computeCenterOfMass(double volume_);
    Solver solver;

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
    bool mouseDown;

    bool mass_fixed;
    bool cotan_fixed;
    bool explicit_mcf;
    float geodesic_dt;
    float mcf_dt;
    float heat_dt;
    int solverIters;
    float solverTol;
    bool enable_geodesics;
    bool enable_heat;
    bool enable_mcf;

    bool enable_iso;
    bool normalize_color;
    int render_color;


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
