add_library(LinAlg INTERFACE)
target_include_directories(LinAlg INTERFACE ${CMAKE_CURRENT_SOURCE_DIR})

add_executable(test main.cpp)

target_link_libraries(test PRIVATE LinAlg)
