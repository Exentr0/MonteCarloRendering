cmake_minimum_required(VERSION 3.30)
project(OpenGLProject)

set(CMAKE_CXX_STANDARD 20)

find_package(OpenGL REQUIRED)
find_package(glfw3 REQUIRED)
find_package(GLEW REQUIRED)
find_package(glm REQUIRED)

set(SOURCES
        main.cpp
)

add_executable(${PROJECT_NAME} ${SOURCES})

target_link_libraries(${PROJECT_NAME}
        glfw
        GLEW::GLEW
        OpenGL::GL
        glm::glm
)
