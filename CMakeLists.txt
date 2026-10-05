cmake_minimum_required(VERSION 3.16)
project(BEVA VERSION 1.0.0 LANGUAGES C)

option(BEVA_BUILD_TESTS "Build tests" ON)
option(BEVA_BUILD_EXAMPLES "Build examples" ON)
option(BEVA_BUILD_SHARED "Build shared library" ON)
option(BEVA_ENABLE_OPENMP "Enable OpenMP when available" ON)
option(BEVA_ENABLE_WERROR "Treat warnings as errors" OFF)

set(CMAKE_C_STANDARD 11)
set(CMAKE_C_STANDARD_REQUIRED ON)
set(CMAKE_C_EXTENSIONS ON)

set(BEVA_SOURCES
    src/status.c src/vector.c src/matrix.c src/lattice.c
    src/qr.c src/number_theory.c src/decomposition.c src/physics.c src/quantum.c
)

add_library(beva_static STATIC ${BEVA_SOURCES})
add_library(BEVA::beva_static ALIAS beva_static)
set_target_properties(beva_static PROPERTIES OUTPUT_NAME beva)
target_include_directories(beva_static PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_compile_definitions(beva_static PRIVATE _CRT_SECURE_NO_WARNINGS)
if(MSVC)
  target_compile_options(beva_static PRIVATE /W4)
else()
  target_compile_options(beva_static PRIVATE -Wall -Wextra -Wpedantic)
  if(BEVA_ENABLE_WERROR)
    target_compile_options(beva_static PRIVATE -Werror)
  endif()
endif()
if(BEVA_ENABLE_OPENMP)
  find_package(OpenMP QUIET)
  if(OpenMP_C_FOUND)
    target_link_libraries(beva_static PUBLIC OpenMP::OpenMP_C)
    target_compile_definitions(beva_static PUBLIC BEVA_HAS_OPENMP=1)
  endif()
endif()
if(NOT MSVC)
  target_link_libraries(beva_static PUBLIC m)
endif()

if(BEVA_BUILD_SHARED)
  add_library(beva_shared SHARED ${BEVA_SOURCES})
  add_library(BEVA::beva_shared ALIAS beva_shared)
  set_target_properties(beva_shared PROPERTIES OUTPUT_NAME beva)
  target_include_directories(beva_shared PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)
  if(NOT MSVC)
    target_link_libraries(beva_shared PUBLIC m)
  endif()
  if(BEVA_ENABLE_OPENMP AND OpenMP_C_FOUND)
    target_link_libraries(beva_shared PUBLIC OpenMP::OpenMP_C)
    target_compile_definitions(beva_shared PUBLIC BEVA_HAS_OPENMP=1)
  endif()
endif()

if(BEVA_BUILD_TESTS)
  enable_testing()
  add_executable(beva_tests tests/test_beva.c)
  target_link_libraries(beva_tests PRIVATE beva_static)
  add_test(NAME beva_tests COMMAND beva_tests)
endif()

if(BEVA_BUILD_EXAMPLES)
  add_executable(beva_example examples/example.c)
  target_link_libraries(beva_example PRIVATE beva_static)
endif()

install(DIRECTORY include/ DESTINATION include)
install(TARGETS beva_static
  ARCHIVE DESTINATION lib
  LIBRARY DESTINATION lib
  RUNTIME DESTINATION bin)
if(BEVA_BUILD_SHARED)
  install(TARGETS beva_shared
    ARCHIVE DESTINATION lib
    LIBRARY DESTINATION lib
    RUNTIME DESTINATION bin)
endif()
install(DIRECTORY docs/ DESTINATION share/doc/beva)
