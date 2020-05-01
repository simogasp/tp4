cmake_minimum_required(VERSION 3.1)

project( OpenGLlabs LANGUAGES CXX VERSION 0.1)

option(BUILD_TESTS "Build tests" OFF)

#########################################################
#
# EXTERNAL LIBRARIES
#
#########################################################

#########################################################
# FIND OPENGL
#########################################################
set(OpenGL_GL_PREFERENCE "GLVND")
find_package(OpenGL REQUIRED)
message(STATUS "OPENGL_gl_LIBRARY: ${OPENGL_gl_LIBRARY}")

#########################################################
# FIND GLUT
#########################################################
if(MSVC)
    set(GLUT_ROOT_PATH "${CMAKE_SOURCE_DIR}/freeglut")
    message(STATUS "GLUT_ROOT_PATH: ${GLUT_ROOT_PATH}")
endif()
find_package(GLUT REQUIRED)
message(STATUS "GLUT_FOUND: ${GLUT_FOUND}")
message(STATUS "GLUT_INCLUDE_DIR: ${GLUT_INCLUDE_DIR}")
message(STATUS "GLUT_LIBRARIES: ${GLUT_LIBRARIES}")

if(CMAKE_SYSTEM_NAME STREQUAL Linux)
    #########################################################
    # FIND Threads, not used but necessary for linkning on linux
    # funny story, even if it is not used it is needed by other
    # dependencies which do not propagate their dependencies (?)
    #########################################################
    find_package(Threads REQUIRED)
    # this force all libs to be include even if not directly used
    set(CMAKE_EXE_LINKER_FLAGS  "-Wl,--no-as-needed")
endif()


#########################################################
# SET COMPILATION FLAGS FOR C++11 
#########################################################
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)


# PUT ALL THE BINARY IN A BIN DIRECTORY
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR})

add_subdirectory(objVisualizer)

if(BUILD_TESTS)
    enable_testing()
    add_test( testEdgeList bin/testEdgeList 1000 )
endif()
