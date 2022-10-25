add_library(utils mcmc.cxx Covariance.cxx CommandHandler.cxx)
target_link_libraries(utils ${ROOT_LIBRARIES})

target_include_directories(utils PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

