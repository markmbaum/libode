
cmake_minimum_required(VERSION 3.18)

project("ode"
		VERSION 0.5.0
		DESCRIPTION "Easy to compile ODE integrators as C++ classes"
		HOMEPAGE_URL "https://github.com/markmbaum/libode"
		LANGUAGES CXX)

include(GNUInstallDirs)
include(CMakePackageConfigHelpers)

if(CMAKE_PROJECT_NAME STREQUAL PROJECT_NAME)
	include(CTest)
endif()

# source files
file (GLOB_RECURSE ODE_SRC CONFIGURE_DEPENDS
       ${PROJECT_SOURCE_DIR}/src/*.cc
       )

# header files
file (GLOB_RECURSE ODE_HDR CONFIGURE_DEPENDS
       ${PROJECT_SOURCE_DIR}/include/ode/*.h
       )

# message("ODE_SRC: ${ODE_HDR}")
# message("ODE_SRC: ${ODE_SRC}")

add_library(${PROJECT_NAME} "${ODE_SRC}")
add_library(${PROJECT_NAME}::${PROJECT_NAME} ALIAS ${PROJECT_NAME})

target_compile_features(${PROJECT_NAME} INTERFACE cxx_std_11)

set_target_properties(${PROJECT_NAME} PROPERTIES PUBLIC_HEADER "${ODE_HDR}")

target_include_directories(
  ${PROJECT_NAME}
  PUBLIC 
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include/${PROJECT_NAME}>
  $<INSTALL_INTERFACE:include/${PROJECT_NAME}>
)

install(TARGETS ${PROJECT_NAME}
        EXPORT ${PROJECT_NAME}_Targets
        ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
        RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
        PUBLIC_HEADER DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/${PROJECT_NAME}
        )

write_basic_package_version_file("${PROJECT_NAME}ConfigVersion.cmake"
                                 VERSION ${PROJECT_VERSION}
                                 COMPATIBILITY SameMajorVersion)
configure_package_config_file(
  "${PROJECT_SOURCE_DIR}/cmake/${PROJECT_NAME}Config.cmake.in"
  "${PROJECT_BINARY_DIR}/${PROJECT_NAME}Config.cmake"
  INSTALL_DESTINATION
  ${CMAKE_INSTALL_DATAROOTDIR}/${PROJECT_NAME}/cmake)

install(EXPORT ${PROJECT_NAME}_Targets
        FILE ${PROJECT_NAME}Targets.cmake
        NAMESPACE ${PROJECT_NAME}::
        DESTINATION ${CMAKE_INSTALL_DATAROOTDIR}/${PROJECT_NAME}/cmake)

install(FILES "${PROJECT_BINARY_DIR}/${PROJECT_NAME}Config.cmake"
              "${PROJECT_BINARY_DIR}/${PROJECT_NAME}ConfigVersion.cmake"
        DESTINATION ${CMAKE_INSTALL_DATAROOTDIR}/${PROJECT_NAME}/cmake)

#install(FILES ${ODE_HDR} 
#        DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/${PROJECT_NAME}
#        )

# =========================================================================
SET (PROJECT_VERSION_MAJOR "0")
SET (PROJECT_VERSION_MINOR "1")
SET (PROJECT_VERSION_PATCH "0")
# ----------------------------

SET (CPACK_RESOURCE_FILE_README "${CMAKE_SOURCE_DIR}/README.md")
set(CPACK_RESOURCE_FILE_LICENSE "${PROJECT_SOURCE_DIR}/LICENSE")

# NEW
SET (CPACK_DEBIAN_PACKAGE_SHLIBDEPS ON)
set (CPACK_DEBIAN_PACKAGE_DEPENDS " ")
 
set(CPACK_PACKAGE_NAME "libode")
SET (CPACK_SOURCE_IGNORE_FILES build/*,./.kdev4,.git/*,.gitignore )

# build a CPack driven installer package
SET (CPACK_RESOURCE_FILE_LICENSE "${PROJECT_SOURCE_DIR}/LICENSE")
SET (CPACK_PACKAGE_DESCRIPTION_SUMMARY "Description of ODE library")
SET (CPACK_PACKAGE_VENDOR "opherden")
SET (CPACK_PACKAGE_DESCRIPTION_FILE "${CPACK_RESOURCE_FILE_README}")
SET (CPACK_PACKAGE_VENDOR "Team")
SET (CPACK_PACKAGE_SECTION "library")
SET (CPACK_PACKAGE_CONTACT "Knut Opherden <knut.opherden@hotmail.com>")
SET (CPACK_DEBIAN_PACKAGE_MAINTAINER "Knut Opherden")
 
SET (CPACK_PACKAGE_VERSION_MAJOR "${PROJECT_VERSION_MAJOR}")
SET (CPACK_PACKAGE_VERSION_MINOR "${PROJECT_VERSION_MINOR}")
SET (CPACK_PACKAGE_VERSION_PATCH "${PROJECT_VERSION_PATCH}")

IF (UNIX)
SET (CPACK_SOURCE_GENERATOR "TGZ;ZIP")
SET (CPACK_GENERATOR "TGZ;TBZ2")
ELSE (UNIX)
SET (CPACK_SOURCE_GENERATOR "ZIP")
SET (CPACK_GENERATOR "NSIS")
ENDIF (UNIX)

# include(Installing)
include(CPack)

# =========================================================================

