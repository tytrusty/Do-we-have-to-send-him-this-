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

- Solver Iters and Tolerance: Ignore this right now. Code currently uses eigen's conjugate gradient solver. Later on I will add a multigrid solver with openmp&mpi support.

### MCF Options
- MCF timestep: timestep to be used for MCF
- Mass fixed: fix the mass matrix. If not fixed, mass matrix will be modified after each timestep due to mesh deformation. 
- Cotan fixed: fix the cotangent laplacian (For conformalized MCF, fix cotan and don't fix the mass)
- MCF is explicit: explicit time integration for MCF (super unstable, dont use lol)

### Heat Options:
- geodesic timestep: timestep to be used when computing geodesics with heat method. For larger values, the result is a smooth approximation of geodesics, and for large t's there are numerical problems. The default timestep used is the surface area/number of faces multiplied by a scalar 5 (this is the suggested value in the referenced Geodesics in Heat paper) 
- flow timestep: timestep to be used when heat flow is enabled

## Render Options:
- Render color: color presets for rendering. 
- isocontour: enables isolines on the surface. Good for demoing the geodesics

## Dependencies

The only dependencies are stl, eigen, [libigl](libigl.github.io/libigl/) and
the dependencies of the `igl::viewer::Viewer` (glfw and opengl).

We recommend you to install libigl using git via:

    git clone --recursive https://github.com/libigl/libigl.git

## References
Crane, K., Weischedel, C., and Wardetzky, M. 2013. Geodesics in heat: A
new approach to computing distance based on heat flow. ACM Trans. Graph.
32, 5. Article 152 (September 2013), 11 pages.
DOI: [http://dx.doi.org/10.1145/2516971.2516977](http://dx.doi.org/10.1145/2516971.2516977).

Kazhdan et al. 2012. Can Mean-Curvature Flow Be Made Non-Singular? [http://www.cs.jhu.edu/~misha/MyPapers/SGP12.pdf](http://www.cs.jhu.edu/~misha/MyPapers/SGP12.pdf)


