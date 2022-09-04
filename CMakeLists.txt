cmake_minimum_required(VERSION 3.16)
project(mean_curvature)

list(PREPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

# Libigl
option(LIBIGL_WITH_OPENGL            "Use OpenGL"         ON)
option(LIBIGL_WITH_OPENGL_GLFW       "Use GLFW"           ON)
include(libigl)

# Add your project files
file(GLOB SRC_FILES *.cpp)
add_executable(${PROJECT_NAME} src/main_mean_curvature.cpp)
add_executable(laplacian_smooth src/laplacian.cpp)
target_link_libraries(${PROJECT_NAME} PUBLIC igl::core igl::opengl_glfw)
target_link_libraries(laplacian_smooth PUBLIC igl::core igl::opengl_glfw)