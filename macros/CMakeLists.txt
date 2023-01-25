# Make sure ROOT loaded before proceeding
if (ROOT_FOUND)
    # add executables and link libraries
    add_executable(analyzePosteriorCov AnalyzePosteriorCov.cxx)
    target_link_libraries(analyzePosteriorCov utils matrix mcmc ${ROOT_LIBRARIES})

    add_executable(analyzePosteriorChains AnalyzePosteriorChains.cxx)
    target_link_libraries(analyzePosteriorChains utils mcmc ${ROOT_LIBRARIES})

    add_executable(calcPosteriorCovmat CalcPosteriorCovmat.cxx)
    target_link_libraries(calcPosteriorCovmat utils matrix ${ROOT_LIBRARIES})

    add_executable(findMCMCStart FindMCMCStart.cxx)
    target_link_libraries(findMCMCStart utils mcmc ${ROOT_LIBRARIES})

    add_executable(generateTargetCovmat GenerateTargetCovmat.cxx)
    target_link_libraries(generateTargetCovmat utils matrix ${ROOT_LIBRARIES})

    add_executable(runMCMC RunMCMC.cxx)
    target_link_libraries(runMCMC utils mcmc ${ROOT_LIBRARIES})

else()
    message(STATUS "ROOT not found! Not building executables!")
endif()

