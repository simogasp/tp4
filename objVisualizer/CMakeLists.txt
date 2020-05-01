add_executable( visualizer main.cpp ObjModelsolution.cpp core.cpp)
target_include_directories(visualizer PUBLIC ${OPENGL_INCLUDE_DIR} ${GLUT_INCLUDE_DIR})
target_link_libraries( visualizer ${OPENGL_opengl_LIBRARY} ${OPENGL_glu_LIBRARY} ${GLUT_LIBRARIES} )
target_compile_options(visualizer PRIVATE -Wall -Wextra -pedantic -Wno-comment)

if("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
    target_compile_options(visualizer PRIVATE -Wno-deprecated-declarations)
endif()

if(CMAKE_SYSTEM_NAME STREQUAL Linux)
    target_link_libraries( visualizer ${CMAKE_THREAD_LIBS_INIT} )
endif()

if(BUILD_TESTS)
    add_executable( testEdge testEdge.cpp core.cpp)
    target_include_directories(testEdge PUBLIC ${OPENGL_INCLUDE_DIR} ${GLUT_INCLUDE_DIR})
    target_link_libraries( testEdge ${OPENGL_opengl_LIBRARY} ${OPENGL_glu_LIBRARY} ${GLUT_LIBRARIES} )

    add_executable( testEdgeList testEdgeList.cpp core.cpp)
    target_include_directories(testEdgeList PUBLIC ${OPENGL_INCLUDE_DIR} ${GLUT_INCLUDE_DIR})
    target_link_libraries( testEdgeList ${OPENGL_opengl_LIBRARY} ${OPENGL_glu_LIBRARY} ${GLUT_LIBRARIES} )
    if(CMAKE_SYSTEM_NAME STREQUAL Linux)
        target_link_libraries( testEdge ${CMAKE_THREAD_LIBS_INIT} )
        target_link_libraries( testEdgeList ${CMAKE_THREAD_LIBS_INIT} )
    endif()
endif()

