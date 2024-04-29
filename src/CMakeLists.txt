cmake_minimum_required(VERSION 3.18 FATAL_ERROR)
cmake_policy(VERSION 3.18.0)

project (libsatsuma
    VERSION 0.0.2
    LANGUAGES C CXX
    )

configure_file (
  "${CMAKE_CURRENT_SOURCE_DIR}/libsatsuma/Config/Version.hh.in"
  "${CMAKE_CURRENT_BINARY_DIR}/libsatsuma/Config/Version.hh"
)
if (TARGET Gurobi::GurobiCXX) 
    set (HAVE_GUROBI 1)
else()
    set (HAVE_GUROBI 0)
endif()

configure_file (
  "${CMAKE_CURRENT_SOURCE_DIR}/libsatsuma/Config/Gurobi.hh.in"
  "${CMAKE_CURRENT_BINARY_DIR}/libsatsuma/Config/Gurobi.hh"
)

if (SATSUMA_ENABLE_BLOSSOM5)
    if (TARGET Blossom5::Blossom5) 
        message("libSatsuma: building WITH blossom-v")
        set (HAVE_BLOSSOM5 1)
    else()
        message(WARNING "libSatsuma: SATSUMA_ENABLE_BLOSSOM5 specified, but target not found.")
        set (HAVE_BLOSSOM5 0)
    endif()
else()
    message("libSatsuma: building WITHOUT blossom-v")
    set (HAVE_BLOSSOM5 0)
endif()
configure_file (
  "${CMAKE_CURRENT_SOURCE_DIR}/libsatsuma/Config/Blossom5.hh.in"
  "${CMAKE_CURRENT_BINARY_DIR}/libsatsuma/Config/Blossom5.hh"
)

add_library(satsuma
    ./libsatsuma/Extra/Highlevel.cc
    ./libsatsuma/Problems/BMatching.cc
    ./libsatsuma/Problems/MCF.cc
    ./libsatsuma/Problems/CostFunction.cc
    ./libsatsuma/Problems/BiFlowGraph.cc
    ./libsatsuma/Problems/Matching.cc
    ./libsatsuma/Problems/BiMCF.cc
    ./libsatsuma/Problems/BiMDF.cc
    ./libsatsuma/Reductions/BiMCF_to_MCF.cc
    ./libsatsuma/Reductions/BiMCF_to_BMatching.cc
    ./libsatsuma/Reductions/BMatching_to_Matching.cc
    ./libsatsuma/Reductions/BiMDF_to_BiMCF.cc
    ./libsatsuma/Reductions/BiMDF_Simplification.cc
    ./libsatsuma/Reductions/OrientableBiMCF_to_MCF.cc
    ./libsatsuma/Solvers/BiMDFDoubleCover.cc
    ./libsatsuma/Solvers/BiMDFGuess.cc
    ./libsatsuma/Solvers/BiMDFRefinement.cc
    ./libsatsuma/Solvers/EvenBiMDF.cc
    ./libsatsuma/Solvers/Matching.cc
    ./libsatsuma/Solvers/MCF.cc
    ./libsatsuma/Solvers/TJoinMST.cc
    ./libsatsuma/Solvers/OrientBinet.cc
    ./libsatsuma/IO/read_bimdf.cc
    ./libsatsuma/IO/write_bimdf.cc
    )
add_library(satsuma::satsuma ALIAS satsuma)

target_link_libraries(satsuma PUBLIC
    lemon::lemon
    Timekeeper::libTimekeeper
    )
if (HAVE_BLOSSOM5)
    target_link_libraries(satsuma PUBLIC Blossom5::Blossom5)
endif()

if (TARGET Gurobi::GurobiCXX) 
    target_link_libraries(satsuma PUBLIC Gurobi::GurobiCXX)
    target_sources(satsuma PRIVATE
        ./libsatsuma/Solvers/BiMCFGurobi.cc
        ./libsatsuma/Solvers/BiMDFGurobi.cc
        )

endif()

target_compile_features(satsuma INTERFACE cxx_std_20)
set_property(TARGET satsuma PROPERTY POSITION_INDEPENDENT_CODE ON)

target_include_directories(satsuma
  PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/>
        $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/>
        $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
)

include(CMakePackageConfigHelpers)
include(GNUInstallDirs)

#[[=
# for this to work, we would need to enhance lemon's cmake
#  export called with target "satsuma" which requires target "lemon" that is not in any export set.

set(INSTALL_CONFIGDIR ${CMAKE_INSTALL_LIBDIR}/cmake/satsuma)

configure_package_config_file(
    "${CMAKE_CURRENT_SOURCE_DIR}/.cmake/satsumaConfig.cmake.in"
    "${CMAKE_CURRENT_BINARY_DIR}/satsumaConfig.cmake"
    INSTALL_DESTINATION ${INSTALL_CONFIGDIR}
    )

write_basic_package_version_file(
    satsumaConfigVersion.cmake
    VERSION ${PROJECT_VERSION}
    COMPATIBILITY SameMinorVersion
    )

install(TARGETS satsuma
    EXPORT satsumaTargets
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
    INCLUDES DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
)

install(EXPORT satsumaTargets
    FILE satsumaTargets.cmake
    NAMESPACE satsuma::
    DESTINATION ${INSTALL_CONFIGDIR}
)

export(EXPORT satsumaTargets
    NAMESPACE satsuma::)

# Install Header Files
install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/satsuma
    DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
    FILES_MATCHING PATTERN "*.hh"
    )

install(FILES
    "${CMAKE_CURRENT_BINARY_DIR}/satsumaConfigVersion.cmake"
    "${CMAKE_CURRENT_BINARY_DIR}/satsumaConfig.cmake"
    DESTINATION ${INSTALL_CONFIGDIR})


install(DIRECTORY
    "${CMAKE_CURRENT_BINARY_DIR}/satsuma/Config"
    DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/satsuma)

#=]]
