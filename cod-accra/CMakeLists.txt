cmake_minimum_required(VERSION "3.14")

project(CODACCRA LANGUAGES Fortran)
set(CMAKE_Fortran_MODULE_DIRECTORY "${CMAKE_BINARY_DIR}/modules")

include(cmake/CPM.cmake)

CPMAddPackage(
  NAME clima
  VERSION 0.4.3
  OPTIONS
    "BUILD_EXECUTABLES OFF"
  GITHUB_REPOSITORY "nicholaswogan/clima"
  GIT_TAG "v0.4.3"
  EXCLUDE_FROM_ALL ON
)

set(COD_TESTS test1 test2)

foreach(COD_TEST ${COD_TESTS})
  add_executable(${COD_TEST} "${CMAKE_SOURCE_DIR}/${COD_TEST}/${COD_TEST}.f90")
  target_link_libraries(${COD_TEST} clima)
  target_include_directories(${COD_TEST} PUBLIC ${CMAKE_Fortran_MODULE_DIRECTORY})
  if ("${CMAKE_Fortran_COMPILER_ID}" MATCHES "GNU")
    target_compile_options(${COD_TEST} PRIVATE -Wunused -Wimplicit-interface -fimplicit-none)
    if (CMAKE_BUILD_TYPE STREQUAL "Debug")
      target_compile_options(${COD_TEST} PRIVATE -fcheck=all,no-array-temps)
    endif()
  endif()
endforeach()
