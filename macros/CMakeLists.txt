# Make sure ROOT loaded before proceeding
if (ROOT_FOUND)
    # add executables and link libraries
    add_executable(runMCMC RunMCMC.cxx)
    target_link_libraries(runMCMC utils matrix mcmc ${ROOT_LIBRARIES})

    add_executable(generateTargetCovmat GenerateTargetCovmat.cxx)
    target_link_libraries(generateTargetCovmat utils matrix ${ROOT_LIBRARIES})

    add_executable(calcPosteriorCovmat CalcPosteriorCovmat.cxx)
    target_link_libraries(calcPosteriorCovmat utils matrix ${ROOT_LIBRARIES})

    add_executable(plotTraces PlotTraces.cxx)
    target_link_libraries(plotTraces utils matrix mcmc ${ROOT_LIBRARIES})

    add_executable(calcSuboptimality CalcSuboptimality.cxx)
    target_link_libraries(calcSuboptimality utils matrix mcmc ${ROOT_LIBRARIES})
else()
    message(STATUS "ROOT not found! Not building executables!")
endif()

