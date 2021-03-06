PROJECT(VACORESS)
CMAKE_MINIMUM_REQUIRED(VERSION 2.8)

OPTION(VACORESS_MAKE_TESTS "Compile applications that use VESimulator?" ON)

# SET(CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake/Modules" CACHE STRING
#     "Path to cmake modules")
INCLUDE(FindMatlab)

SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pg")
SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -pg")

IF(NOT CMAKE_BUILD_TYPE)
    SET(CMAKE_BUILD_TYPE "Release" CACHE STRING "Build type (Release or Debug)"
        FORCE)
ENDIF(NOT CMAKE_BUILD_TYPE)

IF(NOT LIBRARY_OUTPUT_PATH)
    SET(LIBRARY_OUTPUT_PATH ${PROJECT_BINARY_DIR}/lib)
ENDIF(NOT LIBRARY_OUTPUT_PATH)

IF(NOT EXECUTABLE_OUTPUT_PATH)
    SET(EXECUTABLE_OUTPUT_PATH ${PROJECT_BINARY_DIR}/bin)
ENDIF(NOT EXECUTABLE_OUTPUT_PATH)

IF(NOT CMAKE_SIZEOF_VOID_P EQUAL 4)
    IF($ENV{MATLAB_ARCH} MATCHES "glnx86")
        SET(FORCE_32_BIT ON CACHE STRING "Force 32 bit compilation")
    ENDIF($ENV{MATLAB_ARCH} MATCHES "glnx86")
ENDIF(NOT CMAKE_SIZEOF_VOID_P EQUAL 4)

FIND_PACKAGE(OpenMP REQUIRED)
FIND_PACKAGE(Matlab REQUIRED)

INCLUDE_DIRECTORIES(${Matlab_INCLUDE_DIRS})

ADD_SUBDIRECTORY(VESimulator)
ADD_SUBDIRECTORY(mexInterface)

