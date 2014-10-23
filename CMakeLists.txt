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


#########################################################
# SET COMPILATION FLAGS FOR C++11 
#########################################################

if (CMAKE_COMPILER_IS_GNUCC)
    execute_process(COMMAND ${CMAKE_C_COMPILER} -dumpversion
                    OUTPUT_VARIABLE GCC_VERSION)
    if (GCC_VERSION VERSION_LESS 4.7 )
        message(STATUS "Version < 4.7")
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++0x")
    else(GCC_VERSION VERSION_LESS 4.7 )
        message(STATUS "Version >= 4.7")
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
    endif(GCC_VERSION VERSION_LESS 4.7 )
endif(CMAKE_COMPILER_IS_GNUCC)


#########################################################
# FIND UNORDERED MAP
#########################################################

INCLUDE(CheckIncludeFileCXX)
CHECK_INCLUDE_FILE_CXX(unordered_map HAVE_STD_UNORDERED_MAP_HEADER)
IF (HAVE_STD_UNORDERED_MAP_HEADER)
	# Even so we've found unordered_map header file it doesn't
	# mean unordered_map and unordered_set will be declared in
	# std namespace.
	#
	# Namely, MSVC 2008 have unordered_map header which declares
	# unordered_map class in std::tr1 namespace. In order to support
	# this, we do extra check to see which exactly namespace is
	# to be used.

	INCLUDE(CheckCXXSourceCompiles)
	CHECK_CXX_SOURCE_COMPILES("#include <unordered_map>
		int main() {
		std::unordered_map<int, int> map;
		return 0;
		}"
		HAVE_UNORDERED_MAP_IN_STD_NAMESPACE)
	IF (HAVE_UNORDERED_MAP_IN_STD_NAMESPACE)
		ADD_DEFINITIONS(-DHAVE_STD_UNORDERED_MAP)
		MESSAGE("-- Found unordered_map/set in std namespace.")
	ELSE (HAVE_UNORDERED_MAP_IN_STD_NAMESPACE)
		CHECK_CXX_SOURCE_COMPILES("#include <unordered_map>
			int main() {
			std::tr1::unordered_map<int, int> map;
			return 0;
			}"
			HAVE_UNORDERED_MAP_IN_TR1_NAMESPACE)
		IF (HAVE_UNORDERED_MAP_IN_TR1_NAMESPACE)
			ADD_DEFINITIONS(-DHAVE_STD_UNORDERED_MAP_IN_TR1_NAMESPACE)
			MESSAGE("-- Found unordered_map/set in std::tr1 namespace.")
		ELSE (HAVE_UNORDERED_MAP_IN_TR1_NAMESPACE)
			MESSAGE("-- Found <unordered_map> but cannot find either std::unordered_map "
				"or std::tr1::unordered_map.")
			MESSAGE("-- Replacing unordered_map/set with map/set (warning: slower!)")
			ADD_DEFINITIONS(-DHAVE_NO_UNORDERED_MAP)
		ENDIF (HAVE_UNORDERED_MAP_IN_TR1_NAMESPACE)
	ENDIF (HAVE_UNORDERED_MAP_IN_STD_NAMESPACE)
ELSE (HAVE_STD_UNORDERED_MAP_HEADER)
	CHECK_INCLUDE_FILE_CXX("tr1/unordered_map" UNORDERED_MAP_IN_TR1_NAMESPACE)
	IF (UNORDERED_MAP_IN_TR1_NAMESPACE)
		ADD_DEFINITIONS(-DHAVE_TR1_UNORDERED_MAP)
		MESSAGE("-- Found tr1/unordered_map/set in std::tr1 namespace.")
	ELSE (UNORDERED_MAP_IN_TR1_NAMESPACE)
		MESSAGE("-- Unable to find <unordered_map> or <tr1/unordered_map>. ")
		MESSAGE("-- Replacing unordered_map/set with map/set (warning: slower!)")
		ADD_DEFINITIONS(-DHAVE_NO_UNORDERED_MAP)
	ENDIF (UNORDERED_MAP_IN_TR1_NAMESPACE)
ENDIF (HAVE_STD_UNORDERED_MAP_HEADER)

#########################################################
# FIND HASH
#########################################################
CHECK_INCLUDE_FILE_CXX(functional HAVE_STD_FUNCTIONAL_HEADER)
IF (HAVE_STD_FUNCTIONAL_HEADER)
	# Even so we've found functional header file it doesn't
	# mean functional and functional will be declared in
	# std namespace.
	#
	# Namely, MSVC 2008 have functional header which declares
	# functional class in std::tr1 namespace. In order to support
	# this, we do extra check to see which exactly namespace is
	# to be used.

	INCLUDE(CheckCXXSourceCompiles)
	CHECK_CXX_SOURCE_COMPILES("#include <functional>
		int main() {
		std::hash<int> fun;
		return 0;
		}"
		HAVE_FUNCTIONAL_IN_STD_NAMESPACE)
	IF (HAVE_FUNCTIONAL_IN_STD_NAMESPACE)
		ADD_DEFINITIONS(-DHAVE_STD_FUNCTIONAL)
		MESSAGE("-- Found functional/hash in std namespace.")
	ELSE (HAVE_FUNCTIONAL_IN_STD_NAMESPACE)
		CHECK_CXX_SOURCE_COMPILES("#include <functional>
			int main() {
			std::tr1::hash<int> fun;
			return 0;
			}"
			HAVE_FUNCTIONAL_IN_TR1_NAMESPACE)
		IF (HAVE_FUNCTIONAL_IN_TR1_NAMESPACE)
			ADD_DEFINITIONS(-DHAVE_STD_HASH_IN_TR1_NAMESPACE)
			MESSAGE("-- Found functional/hash in std::tr1 namespace.")
		ELSE (HAVE_FUNCTIONAL_IN_TR1_NAMESPACE) 
			MESSAGE( FATAL_ERROR "-- Found <functional> but cannot find either std::hash "
				"or std::tr1::hash.")
		ENDIF (HAVE_FUNCTIONAL_IN_TR1_NAMESPACE)
	ENDIF (HAVE_FUNCTIONAL_IN_STD_NAMESPACE)
ELSE (HAVE_STD_FUNCTIONAL_HEADER)
	CHECK_INCLUDE_FILE_CXX("tr1/functional" FUNCTIONAL_IN_TR1_NAMESPACE)
	IF (FUNCTIONAL_IN_TR1_NAMESPACE)
		ADD_DEFINITIONS(-DHAVE_TR1_FUNCTIONAL)
		MESSAGE("-- Found tr1/functional/set in std::tr1 namespace.")
	ELSE (FUNCTIONAL_IN_TR1_NAMESPACE)
		MESSAGE(FATAL_ERROR "-- Unable to find <functional> or <tr1/functional>. ")
	ENDIF (FUNCTIONAL_IN_TR1_NAMESPACE)
ENDIF (HAVE_STD_FUNCTIONAL_HEADER)

# PUT ALL THE BINARY IN A BIN DIRECTORY
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/bin)



if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
	# suppress the hideous deprecated warning
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-deprecated")
endif()


add_subdirectory(objVisualizer)


