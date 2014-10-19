add_executable( visualizer glutobjsolution.cpp ObjModelsolution.cpp core.cpp)
target_link_libraries( visualizer ${OPENGL_LIBRARIES} ${GLUT_LIBRARIES} )

add_executable( testEdge testEdge.cpp core.cpp)
target_link_libraries( testEdge ${OPENGL_LIBRARIES} ${GLUT_LIBRARIES} )
