cmake_minimum_required(VERSION 2.8)

project(tBall)

find_package(VTK REQUIRED)
find_package(OpenCV REQUIRED)

INCLUDE(${VTK_USE_FILE})

add_executable(tBall main.cpp myInteractorStyle.h myInteractorStyle.cpp 
	Pool.h Pool.cpp drawing.h drawing.cpp Cover.h Cover.cpp connectivity_tumor.h
	connectivity_tumor.cpp Miniball.hpp)

target_link_libraries(tBall ${VTK_LIBRARIES} ${OpenCV_LIBS})