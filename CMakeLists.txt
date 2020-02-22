cmake_minimum_required(VERSION 3.00)

# set the project name
project(MkSPH)

add_compile_options(-Wno-format-security -Wno-format)
add_definitions(-D__GL_H__ -DMKDEBUG)

SET(CMAKE_CXX_FLAGS "-std=c++11") 
SET(CMAKE_CXX_STANDARD 11)

# SRC_FILES   Test.cpp  MkSparseMatrix.cpp MkKrig.cpp    MkSpline.cpp    MkBaseGenetic.cpp  MkDimUnit.cpp MkBndCon.cpp MkDOF.cpp MkEntity.cpp MkExcavStep.cpp MkGrid.cpp MkMesh.cpp  MkSmoothFunc.cpp rngs.c rvgs.c MkSPH.cpp MkStress.cpp MkFailure.cpp MkStrain.cpp  

# add library
#add_library(libMkSPH.a ${SRC_FILES}) 

SET ( CMAKE_BUILD_TYPE Debug)
SET ( SRC_FILES  MkLiuTest.cpp MkArray.cpp MkContainer.cpp MkMatrix.cpp MkPoint.cpp MkPointData.cpp  MkOpenGL.cpp  MkMisc.cpp  MkColor.cpp MkShape.cpp MkLine.cpp MkPolygon.cpp MkPlane.cpp MkSphere.cpp MkCube.cpp MkCircle.cpp MkArc.cpp MkCylinder.cpp MkRect.cpp MkTriangle.cpp MkLiuSPH_mk2.cpp MkLiuParam.cpp MkLiuParticle.cpp MkLiuKernel.cpp MkLiuPair.cpp  MkLiuGrid.cpp MkLiuBound.cpp )

link_directories(.)
link_libraries( -lGL -lGLU -lglut)

# add the executable
add_executable(Test  ${SRC_FILES})
