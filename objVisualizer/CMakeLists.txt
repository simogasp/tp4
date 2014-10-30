add_executable( visualizer main.cpp ObjModelsolution.cpp core.cpp)
target_link_libraries( visualizer ${OPENGL_LIBRARIES} ${GLUT_LIBRARIES} )

add_executable( testEdge testEdge.cpp core.cpp)
target_link_libraries( testEdge ${OPENGL_LIBRARIES} ${GLUT_LIBRARIES} )

add_executable( testEdgeList testEdgeList.cpp core.cpp)
target_link_libraries( testEdgeList ${OPENGL_LIBRARIES} ${GLUT_LIBRARIES} )
