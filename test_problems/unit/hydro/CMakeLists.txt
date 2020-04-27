cmake_minimum_required(VERSION 3.8)

set(Sources
    Test.cpp
)

add_executable(Hydro_test  ${Sources})

target_link_libraries(Hydro_test PUBLIC
    gtest_main
    hydrolib
)

add_test(
    NAME Hydro_test
    COMMAND Hydro_test
)