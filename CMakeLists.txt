cmake_minimum_required(VERSION 3.0)
include(ExternalProject)

project( OpenGLlabs LANGUAGES CXX VERSION 0.1)

option(WITH_UNORDERED_MAP "Use unordered map for efficiency" ON)
option(BUILD_TESTS "Build tests" OFF)

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
# message( "${OPENGL_LIBRARIES}" )
# message( "${OPENGL_LIBRARIES_DIRS}" )
# message( "${OPENGL_DEFINITIONS}" )
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

if (CMAKE_VERSION VERSION_LESS "3.1")
    set(CMAKE_CXX_FLAGS "-std=gnu++11 ${CMAKE_CXX_FLAGS}")
else ()
    set(CMAKE_CXX_STANDARD 11)
    set(CMAKE_CXX_STANDARD_REQUIRED ON)
endif ()

set(VISUALIZER_DEFINITIONS "")
if (WITH_UNORDERED_MAP)
#########################################################
# FIND UNORDERED MAP
#########################################################

    INCLUDE(CheckIncludeFileCXX)
    check_include_file_cxx(unordered_map HAVE_STD_UNORDERED_MAP_HEADER)
    if (HAVE_STD_UNORDERED_MAP_HEADER)
        # Even so we've found unordered_map header file it doesn't
        # mean unordered_map and unordered_set will be declared in
        # std namespace.
        #
        # Namely, MSVC 2008 have unordered_map header which declares
        # unordered_map class in std::tr1 namespace. In order to support
        # this, we do extra check to see which exactly namespace is
        # to be used.

        INCLUDE(CheckCXXSourceCompiles)
        check_cxx_source_compiles("#include <unordered_map>
            int main() {
            std::unordered_map<int, int> map;
            return 0;
            }"
            HAVE_UNORDERED_MAP_IN_STD_NAMESPACE)
        if (HAVE_UNORDERED_MAP_IN_STD_NAMESPACE)
            list(APPEND VISUALIZER_DEFINITIONS -DHAVE_STD_UNORDERED_MAP)
            message("-- Found unordered_map/set in std namespace.")
        else (HAVE_UNORDERED_MAP_IN_STD_NAMESPACE)
            check_cxx_source_compiles("#include <unordered_map>
                int main() {
                std::tr1::unordered_map<int, int> map;
                return 0;
                }"
                HAVE_UNORDERED_MAP_IN_TR1_NAMESPACE)
            if (HAVE_UNORDERED_MAP_IN_TR1_NAMESPACE)
                list(APPEND VISUALIZER_DEFINITIONS -DHAVE_STD_UNORDERED_MAP_IN_TR1_NAMESPACE)
                message("-- Found unordered_map/set in std::tr1 namespace.")
            else (HAVE_UNORDERED_MAP_IN_TR1_NAMESPACE)
                message("-- Found <unordered_map> but cannot find either std::unordered_map "
                    "or std::tr1::unordered_map.")
                message("-- Replacing unordered_map/set with map/set (warning: slower!)")
                list(APPEND VISUALIZER_DEFINITIONS -DHAVE_NO_UNORDERED_MAP)
            endif (HAVE_UNORDERED_MAP_IN_TR1_NAMESPACE)
        endif (HAVE_UNORDERED_MAP_IN_STD_NAMESPACE)
    else (HAVE_STD_UNORDERED_MAP_HEADER)
        check_include_file_cxx("tr1/unordered_map" UNORDERED_MAP_IN_TR1_NAMESPACE)
        if (UNORDERED_MAP_IN_TR1_NAMESPACE)
            list(APPEND VISUALIZER_DEFINITIONS -DHAVE_TR1_UNORDERED_MAP)
            message("-- Found tr1/unordered_map/set in std::tr1 namespace.")
        else (UNORDERED_MAP_IN_TR1_NAMESPACE)
            message("-- Unable to find <unordered_map> or <tr1/unordered_map>. ")
            message("-- Replacing unordered_map/set with map/set (warning: slower!)")
            list(APPEND VISUALIZER_DEFINITIONS -DHAVE_NO_UNORDERED_MAP)
        endif (UNORDERED_MAP_IN_TR1_NAMESPACE)
    endif (HAVE_STD_UNORDERED_MAP_HEADER)

    #########################################################
    # FIND HASH
    #########################################################
    check_include_file_cxx(functional HAVE_STD_FUNCTIONAL_HEADER)
    if (HAVE_STD_FUNCTIONAL_HEADER)
        # Even so we've found functional header file it doesn't
        # mean functional and functional will be declared in
        # std namespace.
        #
        # Namely, MSVC 2008 have functional header which declares
        # functional class in std::tr1 namespace. In order to support
        # this, we do extra check to see which exactly namespace is
        # to be used.

        include(CheckCXXSourceCompiles)
        check_cxx_source_compiles("#include <functional>
            int main() {
            std::hash<int> fun;
            return 0;
            }"
            HAVE_FUNCTIONAL_IN_STD_NAMESPACE)
        if (HAVE_FUNCTIONAL_IN_STD_NAMESPACE)
            list(APPEND VISUALIZER_DEFINITIONS -DHAVE_STD_FUNCTIONAL)
            message("-- Found functional/hash in std namespace.")
        else (HAVE_FUNCTIONAL_IN_STD_NAMESPACE)
            check_cxx_source_compiles("#include <functional>
                int main() {
                std::tr1::hash<int> fun;
                return 0;
                }"
                HAVE_FUNCTIONAL_IN_TR1_NAMESPACE)
            if (HAVE_FUNCTIONAL_IN_TR1_NAMESPACE)
                list(APPEND VISUALIZER_DEFINITIONS -DHAVE_STD_HASH_IN_TR1_NAMESPACE)
                message("-- Found functional/hash in std::tr1 namespace.")
            else (HAVE_FUNCTIONAL_IN_TR1_NAMESPACE)
                message( FATAL_ERROR "-- Found <functional> but cannot find either std::hash "
                    "or std::tr1::hash.")
            endif (HAVE_FUNCTIONAL_IN_TR1_NAMESPACE)
        endif (HAVE_FUNCTIONAL_IN_STD_NAMESPACE)
    else (HAVE_STD_FUNCTIONAL_HEADER)
        check_include_file_cxx("tr1/functional" FUNCTIONAL_IN_TR1_NAMESPACE)
        if (FUNCTIONAL_IN_TR1_NAMESPACE)
            list(APPEND VISUALIZER_DEFINITIONS -DHAVE_TR1_FUNCTIONAL)
            message("-- Found tr1/functional/set in std::tr1 namespace.")
        else (FUNCTIONAL_IN_TR1_NAMESPACE)
            message(FATAL_ERROR "-- Unable to find <functional> or <tr1/functional>. ")
        endif (FUNCTIONAL_IN_TR1_NAMESPACE)
    endif (HAVE_STD_FUNCTIONAL_HEADER)

else(WITH_UNORDERED_MAP)
    list(APPEND VISUALIZER_DEFINITIONS -DHAVE_NO_UNORDERED_MAP)
endif(WITH_UNORDERED_MAP)

message(STATUS "VISUALIZER_DEFINITIONS: ${VISUALIZER_DEFINITIONS}")

# PUT ALL THE BINARY IN A BIN DIRECTORY
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR})



if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
    # suppress the hideous deprecated warning
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-deprecated")
endif()


add_subdirectory(objVisualizer)

if(BUILD_TESTS)
    enable_testing()
    add_test( testEdgeList bin/testEdgeList 1000 )
endif()