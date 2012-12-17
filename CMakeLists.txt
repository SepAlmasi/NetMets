#Specify the version being used aswell as the language
cmake_minimum_required(VERSION 2.8)
#Name your project here
project(NetMets)

#set the module directory
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}")

#find OpenGL
find_package(OpenGL REQUIRED)

#find GLUT
set(GLUT_ROOT_PATH $ENV{GLUT_ROOT_PATH})
find_package(GLUT REQUIRED)

#find GLEW
find_package(GLEW REQUIRED)

#find BOOST
find_package(Boost REQUIRED)

#find the Approximate Nearest Neighbor Library
find_package(ANN REQUIRED)

#ask the user for the RTS location
set(RTS_ROOT_PATH $ENV{RTS_ROOT_PATH})
find_package(RTS REQUIRED)

#set the include directories
include_directories(
	${CMAKE_CURRENT_BINARY_DIR}
	${OPENGL_INCLUDE_DIR}
	${GLEW_INCLUDE_PATH}
	${GLUT_INCLUDE_DIR}
	${RTS_INCLUDE_DIR}
	${ANN_INCLUDE_DIR}
	${Boost_INCLUDE_DIR}
)

#Assign source files to the appropriate variables
file(GLOB SRC_CPP "*.cpp")
file(GLOB SRC_H "*.h")

#set up copying data files
configure_file(00_GT.obj ${CMAKE_CURRENT_BINARY_DIR}/00_GT.obj @ONLY)
configure_file(00_T.obj ${CMAKE_CURRENT_BINARY_DIR}/00_T.obj @ONLY)
configure_file(SmoothShader_Vertex.glsl ${CMAKE_CURRENT_BINARY_DIR}/SmoothShader_Vertex.glsl @ONLY)
configure_file(SmoothShader_Fragment.glsl ${CMAKE_CURRENT_BINARY_DIR}/SmoothShader_Fragment.glsl @ONLY)
configure_file(ErrorMap_Fragment.glsl ${CMAKE_CURRENT_BINARY_DIR}/ErrorMap_Fragment.glsl @ONLY)

#create an executable
add_executable(NetMets ${SRC_CPP} ${SRC_H})

#set the link libraries
target_link_libraries(NetMets ${ANN_LIBRARY} ${OPENGL_gl_LIBRARY} ${OPENGL_glu_LIBRARY} ${GLEW_LIBRARY} ${GLUT_glut_LIBRARY})



