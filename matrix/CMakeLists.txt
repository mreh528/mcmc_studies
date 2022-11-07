add_library(matrix Covariance.cxx Metrics.h)
target_link_libraries(matrix ${ROOT_LIBRARIES})

target_include_directories(matrix PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

