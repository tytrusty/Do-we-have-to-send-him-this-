set(BOUNDINGMESH_VERSION "0.2")

include("${CMAKE_CURRENT_LIST_DIR}/boundingmesh-export.cmake")
set(BOUNDINGMESH_INCLUDE_DIRS "/usr/local/include;/usr/include/eigen3;/usr/include")
set(BOUNDINGMESH_CXX_FLAGS " -std=c++0x")
set(BOUNDINGMESH_LIBRARY_DIRS "/usr/local/lib")
set(BOUNDINGMESH_LIBRARIES "boundingmesh;/usr/lib/x86_64-linux-gnu/libCGAL.so;/usr/lib/x86_64-linux-gnu/libboost_thread.so;/usr/lib/x86_64-linux-gnu/libboost_system.so;/usr/lib/x86_64-linux-gnu/libboost_chrono.so;/usr/lib/x86_64-linux-gnu/libboost_date_time.so;/usr/lib/x86_64-linux-gnu/libboost_atomic.so;/usr/lib/x86_64-linux-gnu/libpthread.so;/usr/lib/x86_64-linux-gnu/libgmp.so")

include_directories(${BOUNDINGMESH_INCLUDE_DIRS})
