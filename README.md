# Heat & Laplacian

A simple simulation framework using libigl and cmake. Based on Alec Jacobson's libigl example project. This project contains some boilerplate that sets up a physical simulation to run in its own thread, with rendering provided by libigl.

## Compile

Compile this project using the standard cmake routine:

    mkdir build
    cd build
    cmake ..
    make

This should find and build the dependencies and create a `example_bin` binary.

## Run

From within the `build` directory just issue:

    ./heat_bin

A glfw app should launch displaying a GUI. The code will try to load mesh files from ../meshes folder, so you need to run from a build subdirectory.

## Simulation Parameters

- Filename: .obj or .ply mesh files accepted

### Simulation Options
For heat flow and heat geodesics, only one of the two should be enabled.
- Enable heat flow: enable heat flow on the surface. Click and drag to add heat to a vertex
- Enable MCF: enables mean curvature flow to smooth the surface 
- Enable heat geodesics: enables computing geodesics with the heat method. Click on the
mesh to set a source.

### MCF Options
- MCF timestep:
- Mass fixed:
- Cotan fixed:
- MCF is explicit:

### Heat Options:
- geodesic timestep: timestep to be used when computing geodesics with heat method. 
- flow timestep: timestep to be used when heat flow is enabled

* Render Options:
- Render color: color presets for rendering. 
- isocontour: enables isolines on the surface
## Dependencies

The only dependencies are stl, eigen, [libigl](libigl.github.io/libigl/) and
the dependencies of the `igl::viewer::Viewer` (glfw and opengl).

We recommend you to install libigl using git via:

    git clone --recursive https://github.com/libigl/libigl.git
