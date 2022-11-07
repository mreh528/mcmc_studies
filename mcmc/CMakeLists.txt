add_library(mcmc mcmc.cxx)
target_link_libraries(mcmc ${ROOT_LIBRARIES})

target_include_directories(mcmc PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

