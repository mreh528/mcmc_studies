# Make sure ROOT loaded before proceeding
if (ROOT_FOUND)
    # add executables and link libraries
    add_executable(mcmcPerformance MCMCPerformance.cxx)
    target_link_libraries(mcmcPerformance utils ${ROOT_LIBRARIES})

    add_executable(generateTargetCovmat GenerateTargetCovmat.cxx)
    target_link_libraries(generateTargetCovmat utils ${ROOT_LIBRARIES})

    add_executable(calcPosteriorCovmat CalcPosteriorCovmat.cxx)
    target_link_libraries(calcPosteriorCovmat utils ${ROOT_LIBRARIES})

    add_executable(plotTraces PlotTraces.cxx)
    target_link_libraries(plotTraces utils ${ROOT_LIBRARIES})
else()
    message(STATUS "ROOT not found! Not building executables!")
endif()

