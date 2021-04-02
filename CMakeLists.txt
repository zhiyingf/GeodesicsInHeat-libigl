# CMakeList.txt : CMake project for GeodesicInHeat, include source and define
# project specific logic here.3.8
#
cmake_minimum_required (VERSION 3.8)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

# libigl
option(LIBIGL_WITH_OPENGL            "Use OpenGL"         ON)
option(LIBIGL_WITH_OPENGL_GLFW       "Use GLFW"           ON)
option(LIBIGL_WITH_OPENGL_GLFW_IMGUI "Use ImGui"          ON)
option(LIBIGL_WITH_PNG               "Use PNG"            ON)

project (GeodesicInHeat)

find_package(LIBIGL REQUIRED QUIET)

# Add source to this project's executable.
file(GLOB SRCFILES src/*.cpp)
file(GLOB HEADERSFILES include/*.h)
add_executable(${PROJECT_NAME} ${HEADERSFILES} ${SRCFILES})


# TODO: Add tests and install targets if needed.
target_link_libraries(${PROJECT_NAME} igl::core igl::opengl_glfw igl::opengl igl::opengl_glfw_imgui igl::png)
include_directories(include/ src/)

