
file(GLOB SOURCES "*.cc")

add_executable(iprins ${SOURCES})

include_directories(caller)
add_subdirectory(caller)
add_subdirectory(fasta)

target_link_libraries(iprins
    caller
  )

if (NOT ${INSTALL_BIN_PREFIX} STREQUAL "")
    message("INSTALL BINARY PATH : " ${INSTALL_BIN_PREFIX})
    set(INSTALL_BIN_DIR ${INSTALL_BIN_PREFIX} CACHE PATH "Installation directory for executables")
    install(TARGETS iprins DESTINATION ${INSTALL_BIN_DIR})
else()
    message("INSTALL BINARY PATH : " "/usr/local/bin")
    set(INSTALL_BIN_DIR "/usr/local/bin" CACHE PATH "Installation directory for executables")
    install(TARGETS iprins DESTINATION ${INSTALL_BIN_DIR})
endif()