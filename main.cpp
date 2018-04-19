#include <igl/opengl/glfw/Viewer.h>
#include <thread>
#include "PhysicsHook.h"
#include "HeatHook.h"
#include <igl/opengl/glfw/imgui/ImGuiMenu.h>
#include <igl/opengl/glfw/imgui/ImGuiHelpers.h>

static PhysicsHook *hook = NULL;

void toggleSimulation()
{
    if (!hook)
        return;

    if (hook->isPaused())
        hook->run();
    else
        hook->pause();
}

void resetSimulation()
{
    if (!hook)
        return;

    hook->reset();
}

bool drawCallback(igl::opengl::glfw::Viewer &viewer)
{
    if (!hook)
        return false;

    hook->render(viewer);
    return false;
}

bool keyCallback(igl::opengl::glfw::Viewer &viewer, unsigned int key, int modifiers)
{
    if (key == ' ')
    {
        toggleSimulation();
        return true;
    }
    return false;
}

bool mouseDownCallback(igl::opengl::glfw::Viewer &viewer, int button, int modifier)
{
    return hook->mouseClicked(viewer, button);    
}

bool mouseUpCallback(igl::opengl::glfw::Viewer &viewer, int button, int modifier)
{
    return hook->mouseReleased(viewer, button);    
}

bool mouseMoveCallback(igl::opengl::glfw::Viewer &viewer, int button, int modifier)
{
    return hook->mouseMoved(viewer, button);    
}


bool drawGUI(igl::opengl::glfw::imgui::ImGuiMenu &menu)
{
    if (ImGui::CollapsingHeader("Simulation Control", ImGuiTreeNodeFlags_DefaultOpen))
    {
        if (ImGui::Button("Run/Pause Sim", ImVec2(-1, 0)))
        {
            toggleSimulation();
        }
        if (ImGui::Button("Reset Sim", ImVec2(-1, 0)))
        {
            resetSimulation();
        }
    }
    hook->drawGUI(menu);
    return false;
}

int main(int argc, char *argv[])
{
    igl::opengl::glfw::Viewer viewer;

    hook = new HeatHook();
    hook->reset();

    viewer.data().set_face_based(true);
    viewer.core.is_animating = true;
    viewer.callback_mouse_down = mouseDownCallback;
    viewer.callback_mouse_up = mouseUpCallback;
    viewer.callback_mouse_move = mouseMoveCallback;
    viewer.callback_key_pressed = keyCallback;
    viewer.callback_pre_draw = drawCallback;

    igl::opengl::glfw::imgui::ImGuiMenu menu;
    viewer.plugins.push_back(&menu);

    menu.callback_draw_viewer_menu = [&]() {drawGUI(menu); };

    viewer.launch();
}
