project(gabp-tests)

add_executable(gabp-tests main.cc matrix.cc)
# target_link_libraries(gabp-tests PRIVATE gabp-lib)
target_include_directories(gabp-tests PUBLIC ../include)

add_test(gabp-tests gabp-tests)
