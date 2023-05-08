set(SATSUMA_TEST_WITH_GCOV OFF CACHE BOOL "Build satsuma tests with gcov support for coverage reports")

if(SATSUMA_TEST_WITH_GCOV)
    set(CMAKE_BUILD_TYPE "Debug")
    add_compile_options("-O0;-g;-fprofile-arcs;-ftest-coverage")
    add_link_options("-lgcov;--coverage")
endif()

include(FetchContent)

FetchContent_Declare(
  googletest
  GIT_REPOSITORY https://github.com/google/googletest.git
  GIT_TAG        release-1.12.1
)
if(WIN32)
    # avoid linking errors, cf https://stackoverflow.com/questions/12540970/how-to-make-gtest-build-mdd-instead-of-mtd-by-default-using-cmake
    set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
endif()
FetchContent_MakeAvailable(googletest)


if (NOT TARGET satsuma::satsuma)
    find_package(satsuma REQUIRED)
endif()


enable_testing()
include(GoogleTest)

add_executable(unittests 
    basic.cc
    )
target_link_libraries(unittests PRIVATE
    satsuma::satsuma
    gtest gtest_main
    ${EXTRA_LIBRARIES}
)
gtest_add_tests(TARGET unittests
    EXTRA_ARGS "--gtest_color=yes"
)
