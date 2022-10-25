# Make sure ROOT loaded before proceeding
if (ROOT_FOUND)
    # add executables and link libraries
    add_executable(mcmcPerformance MCMCPerformance.cxx)
    target_link_libraries(mcmcPerformance utils ${ROOT_LIBRARIES})
    
    add_executable(generateCovmat GenerateCovmat.cxx)
    target_link_libraries(generateCovmat utils ${ROOT_LIBRARIES})
else()
    message(STATUS "ROOT not found! Not building executables!")
endif()

