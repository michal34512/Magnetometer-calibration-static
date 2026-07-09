cmake_minimum_required(VERSION 3.16)

project(mag_callib_static C)

set(CMAKE_C_STANDARD 11)
set(CMAKE_C_STANDARD_REQUIRED ON)

add_executable(mag_callib_static
    main.c
    vector.c
    calibration.c
    matrix.c
    eig.c
)

if(MSVC)
    target_compile_options(mag_callib_static PRIVATE /W4)
else()
    target_compile_options(mag_callib_static PRIVATE -Wall -Wextra -Wpedantic)
    target_link_libraries(mag_callib_static PRIVATE m)
endif()
