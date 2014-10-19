add_executable( visualizer glutobjsolution.cpp ObjModelsolution.cpp core.cpp)
target_link_libraries( visualizer ${OPENGL_LIBRARIES} ${GLUT_LIBRARIES} )
