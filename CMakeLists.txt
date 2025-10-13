cmake_minimum_required(VERSION 3.10)

project(mos_cv C)

# Compiler flags
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wall -O2 -fopenmp")

# Explicitly list source files
set(SRC_FILES
    main.c
    parameter_fetch/parameter_fetch.c
    carrier_continuity/carrier.c
    poisson/poisson.c
    solve_c/solve_c.c
    fileio/fileio.c
)

# Create executable
add_executable(mos_cv ${SRC_FILES})

# Link with math library
target_link_libraries(mos_cv m)

# Link OpenMP if found
if(OpenMP_C_FOUND)
    target_link_libraries(mos_cv PUBLIC OpenMP::OpenMP_C)
endif()
