cmake_minimum_required(VERSION 3.10)
project(WaterPaths)

set(CMAKE_CXX_STANDARD 14)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp -O2 -march=native")
#Debug: não pode ter otimização (usar O0)
#Otimização: O2 -> utilizar esse
#Otimização mais agressiva: O3

include_directories(src/ContinuityModels)
include_directories(src/ContinuityModels/Base)
include_directories(src/Controls)
include_directories(src/Controls/Base)
include_directories(src/Controls/Custom)
include_directories(src/DataCollector)
include_directories(src/DataCollector/Base)
include_directories(src/DroughtMitigationInstruments)
include_directories(src/DroughtMitigationInstruments/Base)
include_directories(src/Problem)
include_directories(src/Problem/Base)
include_directories(src/Simulation)
include_directories(src/SystemComponents)
include_directories(src/SystemComponents/Bonds)
include_directories(src/SystemComponents/Bonds/Base)
include_directories(src/SystemComponents/Utility)
include_directories(src/SystemComponents/WaterSources)
include_directories(src/SystemComponents/WaterSources/Base)
include_directories(src/Utils)
include_directories(src/Utils/Graph)
include_directories(src/Utils/QPSolver)
#include_directories(Tests)
include_directories(Utils)

add_executable(WaterPaths
        src/ContinuityModels/Base/ContinuityModel.cpp
        src/ContinuityModels/Base/ContinuityModel.h
        src/ContinuityModels/ContinuityModelRealization.cpp
        src/ContinuityModels/ContinuityModelRealization.h
        src/ContinuityModels/ContinuityModelROF.cpp
        src/ContinuityModels/ContinuityModelROF.h
        src/Controls/Base/ControlRules.h
        src/Controls/Base/MinEnvFlowControl.cpp
        src/Controls/Base/MinEnvFlowControl.h
        src/Controls/Base/MinEnvironFlowControl.cpp
        src/Controls/Base/MinEnvironFlowControl.h
        src/Controls/Custom/FallsLakeMinEnvFlowControl.cpp
        src/Controls/Custom/FallsLakeMinEnvFlowControl.h
        src/Controls/Custom/JordanLakeMinEnvFlowControl.cpp
        src/Controls/Custom/JordanLakeMinEnvFlowControl.h
        src/Controls/EvaporationSeries.cpp
        src/Controls/EvaporationSeries.h
        src/Controls/FixedMinEnvFlowControl.cpp
        src/Controls/FixedMinEnvFlowControl.h
        src/Controls/InflowMinEnvFlowControl.cpp
        src/Controls/InflowMinEnvFlowControl.h
        src/Controls/SeasonalMinEnvFlowControl.cpp
        src/Controls/SeasonalMinEnvFlowControl.h
        src/Controls/StorageMinEnvFlowControl.cpp
        src/Controls/StorageMinEnvFlowControl.h
        src/Controls/WwtpDischargeRule.cpp
        src/Controls/WwtpDischargeRule.h
        src/DataCollector/Base/DataCollector.cpp
        src/DataCollector/Base/DataCollector.h
        src/DataCollector/AllocatedReservoirDataCollector.cpp
        src/DataCollector/AllocatedReservoirDataCollector.h
        src/DataCollector/EmptyDataCollector.cpp
        src/DataCollector/EmptyDataCollector.h
        src/DataCollector/IntakeDataCollector.cpp
        src/DataCollector/IntakeDataCollector.h
        src/DataCollector/MasterDataCollector.cpp
        src/DataCollector/MasterDataCollector.h
        src/DataCollector/QuaryDataCollector.cpp
        src/DataCollector/QuaryDataCollector.h
        src/DataCollector/ReservoirDataCollector.cpp
        src/DataCollector/ReservoirDataCollector.h
        src/DataCollector/RestrictionsDataCollector.cpp
        src/DataCollector/RestrictionsDataCollector.h
        src/DataCollector/TransfersDataCollector.cpp
        src/DataCollector/TransfersDataCollector.h
        src/DataCollector/UtilitiesDataCollector.cpp
        src/DataCollector/UtilitiesDataCollector.h
        src/DataCollector/WaterReuseDataCollector.cpp
        src/DataCollector/WaterReuseDataCollector.h
        src/DroughtMitigationInstruments/Base/DroughtMitigationPolicy.cpp
        src/DroughtMitigationInstruments/Base/DroughtMitigationPolicy.h
        src/DroughtMitigationInstruments/InsuranceStorageToROF.cpp
        src/DroughtMitigationInstruments/InsuranceStorageToROF.h
        src/DroughtMitigationInstruments/Restrictions.cpp
        src/DroughtMitigationInstruments/Restrictions.h
        src/DroughtMitigationInstruments/Transfers.cpp
        src/DroughtMitigationInstruments/Transfers.h
        src/Problem/Base/Problem.cpp
        src/Problem/Base/Problem.h
        src/Problem/Caesb.cpp
        src/Problem/Caesb.h
        src/Simulation/Simulation.cpp
        src/Simulation/Simulation.h
        src/SystemComponents/Bonds/Base/Bond.cpp
        src/SystemComponents/Bonds/Base/Bond.h
        src/SystemComponents/Bonds/BalloonPaymentBond.cpp
        src/SystemComponents/Bonds/BalloonPaymentBond.h
        src/SystemComponents/Bonds/FloatingInterestBalloonPaymentBond.cpp
        src/SystemComponents/Bonds/FloatingInterestBalloonPaymentBond.h
        src/SystemComponents/Bonds/LevelDebtServiceBond.cpp
        src/SystemComponents/Bonds/LevelDebtServiceBond.h
        src/SystemComponents/Utility/InfrastructureManager.cpp
        src/SystemComponents/Utility/InfrastructureManager.h
        src/SystemComponents/Utility/Utility.cpp
        src/SystemComponents/Utility/Utility.h
        src/SystemComponents/WaterSources/Base/WaterSource.cpp
        src/SystemComponents/WaterSources/Base/WaterSource.h
        src/SystemComponents/WaterSources/AllocatedReservoir.cpp
        src/SystemComponents/WaterSources/AllocatedReservoir.h
        src/SystemComponents/WaterSources/Intake.cpp
        src/SystemComponents/WaterSources/Intake.h
        src/SystemComponents/WaterSources/Quarry.cpp
        src/SystemComponents/WaterSources/Quarry.h
        src/SystemComponents/WaterSources/Relocation.cpp
        src/SystemComponents/WaterSources/Relocation.h
        src/SystemComponents/WaterSources/Reservoir.cpp
        src/SystemComponents/WaterSources/Reservoir.h
        src/SystemComponents/WaterSources/ReservoirExpansion.cpp
        src/SystemComponents/WaterSources/ReservoirExpansion.h
        src/SystemComponents/WaterSources/SequentialJointTreatmentExpansion.cpp
        src/SystemComponents/WaterSources/SequentialJointTreatmentExpansion.h
        src/SystemComponents/WaterSources/WaterReuse.cpp
        src/SystemComponents/WaterSources/WaterReuse.h
        src/SystemComponents/Catchment.cpp
        src/SystemComponents/Catchment.h
        src/Utils/Graph/Graph.cpp
        src/Utils/Graph/Graph.h
        src/Utils/QPSolver/Array.cpp
        src/Utils/QPSolver/Array.h
        src/Utils/QPSolver/QuadProg++.cpp
        src/Utils/QPSolver/QuadProg++.h
        src/Utils/Constants.h
        src/Utils/DataSeries.cpp
        src/Utils/DataSeries.h
        src/Utils/Matrices.h
        src/Utils/ObjectivesCalculator.cpp
        src/Utils/ObjectivesCalculator.h
        src/Utils/Solutions.h
        src/Utils/Utils.cpp
        src/Utils/Utils.h
        src/main.cpp
#        Tests/catch.hpp
#        Tests/Tests.cpp
        src/DroughtMitigationInstruments/TransfersBilateral.cpp src/DroughtMitigationInstruments/TransfersBilateral.h src/DataCollector/TransfersBilateralDataCollector.cpp src/DataCollector/TransfersBilateralDataCollector.h)
