cmake_minimum_required(VERSION 3.0.0)
project(llbp VERSION 0.1.0 LANGUAGES C CXX)


set(CMAKE_CXX_STANDARD 20)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

set(CMAKE_CXX_FLAGS "-Wall -Wextra -Wno-unused-parameter")
set(CMAKE_CXX_FLAGS_DEBUG "-g")
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -g")



SET(BASEPATH ${CMAKE_SOURCE_DIR})
INCLUDE_DIRECTORIES(${BASEPATH})



add_subdirectory(utils)
add_subdirectory(bpmodels)


## Create executable
set(TARGET_NAME predictor)
set(SOURCES
    main.cc
    btb/btb.cc
    )

add_executable(${TARGET_NAME} ${SOURCES})

target_link_libraries(${TARGET_NAME}
                       utils
                       bpmodels)
target_compile_options(${TARGET_NAME} PRIVATE -Wall)


## Add boost library
find_package(Boost COMPONENTS program_options REQUIRED)
include_directories( ${Boost_INCLUDE_DIRS} )
target_link_libraries(${TARGET_NAME} ${Boost_LIBRARIES})




