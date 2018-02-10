cmake_minimum_required(VERSION 2.8.12)
project (KMapSolver)
set(CMAKE_CXX_FLAGS "-std=c++11 -stdlib=libc++ ${CMAKE_CXX_FLAGS}")

if (MSVC)
  # /wd4005  macro-redefinition
  # /wd4068  unknown pragma
  # /wd4244  conversion from 'type1' to 'type2'
  # /wd4267  conversion from 'size_t' to 'type2'
  # /wd4800  force value to bool 'true' or 'false' (performance warning)
  add_compile_options(/wd4005 /wd4068 /wd4244 /wd4267 /wd4800)
  add_definitions(/DNOMINMAX /DWIN32_LEAN_AND_MEAN=1 /D_CRT_SECURE_NO_WARNINGS)
endif()

add_subdirectory(cctz)
add_subdirectory(abseil-cpp)
include_directories(${abseil_SOURCE_DIR}/include ${abseil_SOURCE_DIR})
# Prevent overriding the parent project's compiler/linker
# settings on Windows
set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
add_subdirectory(googletest)

add_library(KMapSolverLib kmap_solver.cc kmap_solver.h)
target_link_libraries(KMapSolverLib absl::strings)
add_executable(KMapSolver main.cc)
target_link_libraries(KMapSolver KMapSolverLib absl::strings)

add_executable(KmapSolverTest kmap_solver_test.cc)
include_directories(${gtest_SOURCE_DIR}/include ${gtest_SOURCE_DIR})
add_dependencies(KmapSolverTest gtest)
target_link_libraries(KmapSolverTest KMapSolverLib gtest gmock gtest_main absl::strings)