cmake_minimum_required(VERSION 3.16)

project(
  racon
  VERSION 2.0.0
  LANGUAGES C CXX
  DESCRIPTION "Racon2 is a consensus module for de novo genome assembly.")

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -pedantic")
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/bin)

include(FetchContent)
include(GNUInstallDirs)

if(CMAKE_SOURCE_DIR STREQUAL PROJECT_SOURCE_DIR)
  set(racon_main_project ON)
endif()
option(racon_build_wrapper "Build wrapper" OFF)

find_package(bioparser 3.0.15 QUIET)
if(NOT bioparser_FOUND)
  FetchContent_Declare(
    bioparser
    GIT_REPOSITORY https://github.com/rvaser/bioparser
    GIT_TAG 3.0.15)

  FetchContent_GetProperties(bioparser)
  if(NOT bioparser_POPULATED)
    FetchContent_Populate(bioparser)
    add_subdirectory(${bioparser_SOURCE_DIR} ${bioparser_BINARY_DIR}
                     EXCLUDE_FROM_ALL)
  endif()
endif()

find_package(biosoup 0.10.0 QUIET)
if(NOT biosoup_FOUND)
  FetchContent_Declare(
    biosoup
    GIT_REPOSITORY https://github.com/tbrekalo/biosoup
    GIT_TAG master)

  FetchContent_GetProperties(biosoup)
  if(NOT biosoup_POPULATED)
    FetchContent_Populate(biosoup)
    add_subdirectory(${biosoup_SOURCE_DIR} ${biosoup_BINARY_DIR}
                     EXCLUDE_FROM_ALL)
  endif()
endif()

find_package(cxxopts 3.1.1 QUIET)
if(NOT cxxopts_FOUND)
  FetchContent_Declare(
    cxxopts
    GIT_REPOSITORY https://github.com/jarro2783/cxxopts
    GIT_TAG v3.1.1)

  FetchContent_GetProperties(cxxopts)
  if(NOT cxxopts_POPULATED)
    FetchContent_Populate(cxxopts)
    add_subdirectory(${cxxopts_SOURCE_DIR} ${cxxopts_BINARY_DIR}
                     EXCLUDE_FROM_ALL)
  endif()
endif()

find_package(edlib 1.2.7 QUIET)
if(NOT edlib_FOUND)
  FetchContent_Declare(
    edlib
    GIT_REPOSITORY https://github.com/martinsos/edlib
    GIT_TAG v1.2.7)

  FetchContent_GetProperties(edlib)
  if(NOT edlib_POPULATED)
    FetchContent_Populate(edlib)
    add_subdirectory(${edlib_SOURCE_DIR} ${edlib_BINARY_DIR} EXCLUDE_FROM_ALL)
  endif()
endif()

find_package(fmt 9.1.0 QUIET)
if(NOT fmt_FOUND)
  FetchContent_Declare(
    fmt
    GIT_REPOSITORY https://github.com/fmtlib/fmt
    GIT_TAG 9.1.0)

  FetchContent_GetProperties(fmt)
  if(NOT fmt_POPULATED)
    FetchContent_Populate(fmt)
    add_subdirectory(${fmt_SOURCE_DIR} ${fmt_BINARY_DIR} EXCLUDE_FROM_ALL)
  endif()
endif()

find_package(spdlog 4.0.8 QUIET)
if(NOT spdlog_FOUND)
  FetchContent_Declare(
    spdlog
    GIT_REPOSITORY https://github.com/gabime/spdlog
    GIT_TAG v1.11.0)

  FetchContent_GetProperties(spdlog)
  if(NOT spdlog_POPULATED)
    FetchContent_Populate(spdlog)
    add_subdirectory(${spdlog_SOURCE_DIR} ${spdlog_BINARY_DIR} EXCLUDE_FROM_ALL)
  endif()
endif()

find_package(spoa 4.0.8 QUIET)
if(NOT spoa_FOUND)
  FetchContent_Declare(
    spoa
    GIT_REPOSITORY https://github.com/rvaser/spoa
    GIT_TAG 4.0.8)

  FetchContent_GetProperties(spoa)
  if(NOT spoa_POPULATED)
    FetchContent_Populate(spoa)
    add_subdirectory(${spoa_SOURCE_DIR} ${spoa_BINARY_DIR} EXCLUDE_FROM_ALL)
  endif()
endif()

FetchContent_Declare(
  wfa2
  GIT_REPOSITORY https://github.com/smarco/WFA2-lib
  GIT_TAG development)

FetchContent_GetProperties(wfa2)
if(NOT wfa2_POPULATED)
  FetchContent_Populate(wfa2)
  add_subdirectory(${wfa2_SOURCE_DIR} ${wfa2_BINARY_DIR} EXCLUDE_FROM_ALL)
endif()

find_package(TBB REQUIRED)
set(racon_sources
    src/logger.cpp
    src/polisher.cpp
    src/overlap.cpp
    src/sequence.cpp
    src/window.cpp
    src/data.cpp
    src/io.cpp)

add_library(racon ${racon_sources})
target_link_libraries(
  racon
  PUBLIC bioparser::bioparser biosoup::biosoup edlib::edlib spoa::spoa
  PRIVATE fmt::fmt tbb wfa2cpp_static)

target_include_directories(
  racon PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/src>)

add_executable(racon_exe src/main.cpp)
target_link_libraries(racon_exe cxxopts::cxxopts fmt::fmt racon)

target_compile_options(
  racon_exe PUBLIC $<$<CONFIG:Debug,RelWithDebInfo>:-fsanitize=address>
                   $<$<CONFIG:Debug,RelWithDebInfo>:-fno-omit-frame-pointer>)
target_link_options(racon_exe PUBLIC
                    $<$<CONFIG:Debug,RelWithDebInfo>:-fsanitize=address>)

set_target_properties(racon_exe PROPERTIES VERSION ${racon_VERSION}
                                           SOVERSION ${racon_VERSION_MAJOR})
configure_file(${PROJECT_SOURCE_DIR}/src/version.h.in
               ${CMAKE_CURRENT_BINARY_DIR}/include/version.h)
target_include_directories(racon_exe
                           PRIVATE ${CMAKE_CURRENT_BINARY_DIR}/include)

set_property(TARGET racon_exe PROPERTY OUTPUT_NAME racon2)
install(TARGETS racon_exe DESTINATION ${CMAKE_INSTALL_BINDIR})

# Add Debian packaging
set(CPACK_GENERATOR "DEB")
set(CPACK_DEBIAN_PACKAGE_MAINTAINER "Robert Vaser")
set(CPACK_PACKAGE_VERSION "${PROJECT_VERSION}")
include(CPack)
