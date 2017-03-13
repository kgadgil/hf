cmake_minimum_required (VERSION 3.1)
project (HartreeFock)
# The version number.
set (HartreeFock_VERSION_MAJOR 1)
set (HartreeFock_VERSION_MINOR 0)
 
set(LIBINT2_DIR "/Users/kg/Software/kg_libint/install")
set(EIGEN3_DIR "/usr/local/Cellar/eigen/3.3.1/include/eigen3")

#this makes the compiler use the c++11 version
set(CMAKE_CXX_STANDARD 11)

# this finds and includes the header directories for my specific project
#the first one searches this folder, the second go into libints libraries and so on
include_directories("${PROJECT_BINARY_DIR}" ${LIBINT2_DIR}/include ${LIBINT2_DIR}/include/libint2 ${EIGEN3_DIR})

#this portion finds and links the directory that contains libints information
set(LIBINT2_LIB "${LIBINT2_DIR}/lib/")
link_directories(${LIBINT2_LIB})
       
# The executables are made using the first command from the source file 
#the second links the directory we found earlier to the executaable to run the code
add_executable(run mp2.cpp)
target_link_libraries(run int2)


