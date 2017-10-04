add_executable( visualizer main.cpp ObjModelsolution.cpp core.cpp)
target_link_libraries( visualizer ${OPENGL_LIBRARIES} ${GLUT_LIBRARIES} )

add_executable( testEdge testEdge.cpp core.cpp)
target_link_libraries( testEdge ${OPENGL_LIBRARIES} ${GLUT_LIBRARIES} )

add_executable( testEdgeList testEdgeList.cpp core.cpp)
target_link_libraries( testEdgeList ${OPENGL_LIBRARIES} ${GLUT_LIBRARIES} )

if(CMAKE_SYSTEM_NAME STREQUAL Linux)
    target_link_libraries( visualizer ${CMAKE_THREAD_LIBS_INIT} )
    target_link_libraries( testEdge ${CMAKE_THREAD_LIBS_INIT} )
    target_link_libraries( testEdgeList ${CMAKE_THREAD_LIBS_INIT} )
endif()
