cmake_minimum_required(VERSION 3.18 FATAL_ERROR)
cmake_policy(VERSION 3.18.0)

project (satsuma
    VERSION 0.0.1
    LANGUAGES C CXX
    )

option(SATSUMA_ENABLE_BLOSSOM5 "Enable Blossom-V (non-free license)" OFF)

if(${CMAKE_SOURCE_DIR} STREQUAL ${CMAKE_BINARY_DIR})
    message(FATAL_ERROR "In-source builds not allowed. Please make a new directory (called a build directory) and run CMake from there.")
endif()

if(${CMAKE_PROJECT_NAME} STREQUAL ${PROJECT_NAME})
    set(STANDALONE_BUILD TRUE)
else()
    set(STANDALONE_BUILD FALSE)
endif()

set(SATSUMA_ENABLE_TESTS ${STANDALONE_BUILD} CACHE BOOL "Build satsuma tests")
set(SATSUMA_ENABLE_TOOLS ${STANDALONE_BUILD} CACHE BOOL "Build satsuma examples/tools")

set(CMAKE_CXX_STANDARD 20)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
set(CMAKE_FIND_PACKAGE_NO_PACKAGE_REGISTRY ON) # deprecated since CMake 3.16
set(CMAKE_FIND_USE_PACKAGE_REGISTRY OFF) # since CMake 3.16

include(FetchContent)
set(FETCHCONTENT_QUIET OFF)
set(EXTERNAL_DIR "${CMAKE_CURRENT_SOURCE_DIR}/external")
set(FETCHCONTENT_UPDATES_DISCONNECTED TRUE)

if(STANDALONE_BUILD)
    if (NOT TARGET lemon::lemon)
        FetchContent_Declare(
        lemon
        GIT_REPOSITORY https://github.com/cgg-bern/lemon
        GIT_TAG        cgg
        SOURCE_DIR "${EXTERNAL_DIR}/lemon"
        )
        FetchContent_MakeAvailable(lemon)
    endif()


    if(SATSUMA_ENABLE_BLOSSOM5)
        if (NOT TARGET Blossom5::Blossom5)
            FetchContent_Declare(blossom5-cmake
                GIT_REPOSITORY https://github.com/cgg-bern/blossom5-cmake
                GIT_TAG main
                SOURCE_DIR "${EXTERNAL_DIR}/blossom5-cmake"
                )
            FetchContent_MakeAvailable(blossom5-cmake)
        endif()
    endif()

    if (NOT TARGET Timekeeper::libTimekeeper)
        FetchContent_Declare(timekeeper
            GIT_REPOSITORY https://github.com/cgg-bern/libTimekeeper
            GIT_TAG main
            SOURCE_DIR "${EXTERNAL_DIR}/timekeeper"
            )
        FetchContent_Populate(timekeeper)
        add_subdirectory("${timekeeper_SOURCE_DIR}/src")
    endif()
endif()

add_subdirectory(src)

if(${SATSUMA_ENABLE_TESTS})
    enable_testing()
    add_subdirectory(tests)
endif()
if(${SATSUMA_ENABLE_TOOLS})
    add_subdirectory(tools)
endif()
