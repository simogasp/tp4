cmake_minimum_required(VERSION 2.8)
include(ExternalProject)

project( OpenGLlabs )

#########################################################
#
# EXTERNAL LIBRARIES
#
#########################################################

#########################################################
# FIND OPENGL
#########################################################
find_package(OpenGL REQUIRED)
include_directories(${OpenGL_INCLUDE_DIRS})
link_directories(${OpenGL_LIBRARY_DIRS})
add_definitions(${OpenGL_DEFINITIONS})
# MESSAGE( "${OPENGL_LIBRARIES}" )
# MESSAGE( "${OPENGL_LIBRARIES_DIRS}" )
# MESSAGE( "${OPENGL_DEFINITIONS}" )
if(NOT OPENGL_FOUND)
    message(ERROR " OPENGL not found!")
endif(NOT OPENGL_FOUND)

#########################################################
# FIND GLUT
#########################################################
find_package(GLUT REQUIRED)
include_directories(${GLUT_INCLUDE_DIR})
link_directories(${GLUT_LIBRARY_DIRS})
add_definitions(${GLUT_DEFINITIONS})
# message( "glut ${GLUT_INCLUDE_DIR}" )
if(NOT GLUT_FOUND)
    message(ERROR " GLUT not found!")
else()
    message( "-- Found GLUT" )
endif(NOT GLUT_FOUND)


set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/bin)

# #Detect c++11 support on Unix system
# if (UNIX)
#   if ("${CMAKE_CXX_COMPILER_ID}" MATCHES "GNU")
#       execute_process(
#           COMMAND ${CMAKE_CXX_COMPILER} -dumpversion OUTPUT_VARIABLE GCC_VERSION)
#       if (NOT (GCC_VERSION VERSION_GREATER 4.7 OR GCC_VERSION VERSION_EQUAL 4.7))
#         if (GCC_VERSION VERSION_GREATER 4.7 OR GCC_VERSION VERSION_EQUAL 4.7)
#           message(STATUS "C++11 activated.")
#           add_definitions("-std=gnu++11")
#         elseif (GCC_VERSION VERSION_GREATER 4.3 OR GCC_VERSION VERSION_EQUAL 4.3)
#           message(WARNING "C++0x activated. If you get any errors update to a compiler which fully supports C++11")
#           add_definitions("-std=gnu++0x")
#         endif ()
#       endif ()
#   elseif ("${CMAKE_CXX_COMPILER_ID}" MATCHES "Clang")
#       set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -stdlib=libc++")
#   else ()
#       message(FATAL_ERROR "Your C++ compiler does not support C++11.")
#   endif ()
# endif (UNIX)

# include(CheckCXXCompilerFlag)
# CHECK_CXX_COMPILER_FLAG("-std=c++11" COMPILER_SUPPORTS_CXX11)
# CHECK_CXX_COMPILER_FLAG("-std=c++0x" COMPILER_SUPPORTS_CXX0X)
# if(COMPILER_SUPPORTS_CXX11)
#         set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
# elseif(COMPILER_SUPPORTS_CXX0X)
#         set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++0x")
#         message("only support c++0x")
# else()
#         message(STATUS "The compiler ${CMAKE_CXX_COMPILER} has no C++11 support. Please use a different C++ compiler.")
# endif()

add_subdirectory(objVisualizer)


