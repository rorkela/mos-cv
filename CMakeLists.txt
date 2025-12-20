cmake_minimum_required(VERSION 3.16)

project(moscv C)

set(CMAKE_C_STANDARD 11)
set(CMAKE_C_STANDARD_REQUIRED ON)

# OpenMP
find_package(OpenMP REQUIRED)

add_executable(moscv
    src/main.c

    src/common/common.c

    src/thomas/thomas.c

    src/carrier_continuity/carrier.c
    src/parameter_fetch/parameter_fetch.c
    src/poisson/poisson.c
    src/solve_c/solve_c.c
    src/fileio/fileio.c
)

target_include_directories(moscv PRIVATE
    src
)

target_compile_options(moscv PRIVATE
    -Wall
    -Wextra
    -O2
)

target_link_libraries(moscv PRIVATE
    m
    OpenMP::OpenMP_C
)
