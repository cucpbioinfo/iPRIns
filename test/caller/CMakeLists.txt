include_directories(${PROJECT_SOURCE_DIR}/src)

add_executable(
    test_caller
    tests.cc
    )

target_link_libraries(
    test_caller
    gtest_main
    caller
    )

add_test(
  NAME
    test_caller_1
  COMMAND
    ${CMAKE_BINARY_DIR}/${CMAKE_INSTALL_BINDIR}/test_caller
  )