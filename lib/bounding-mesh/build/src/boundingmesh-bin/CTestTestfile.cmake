# CMake generated Testfile for 
# Source directory: /home/ty/Physim/bounding-mesh-0.2/src/boundingmesh-bin
# Build directory: /home/ty/Physim/bounding-mesh-0.2/build/src/boundingmesh-bin
# 
# This file includes the relevant testing commands required for 
# testing this directory and lists subdirectories to be tested as well.
add_test(testBoundingmeshTorus "/home/ty/Physim/bounding-mesh-0.2/build/boundingmesh" "-v" "100" "/home/ty/Physim/bounding-mesh-0.2/examples/torus/torus.off" "/home/ty/Physim/bounding-mesh-0.2/build/torus_decimated.off")
add_test(testBoundingmeshTeapot "/home/ty/Physim/bounding-mesh-0.2/build/boundingmesh" "-v" "1000" "/home/ty/Physim/bounding-mesh-0.2/examples/teapot/teapot.off" "/home/ty/Physim/bounding-mesh-0.2/build/teapot_decimated.off")
add_test(testBoundingConvexDecompositionBunny "/home/ty/Physim/bounding-mesh-0.2/build/bounding-convex-decomposition" "-x" "20000" "/home/ty/Physim/bounding-mesh-0.2/examples/bunny/bunny.off" "/home/ty/Physim/bounding-mesh-0.2/build/bunny_decomposition.wrl")
add_test(testBoundingConvexDecompositionTeapot "/home/ty/Physim/bounding-mesh-0.2/build/bounding-convex-decomposition" "/home/ty/Physim/bounding-mesh-0.2/examples/teapot/teapot.off" "/home/ty/Physim/bounding-mesh-0.2/build/teapot_decomposition.wrl")
