set(MODULE_NAME "ITSSimulation")

O2_SETUP(NAME ${MODULE_NAME})

set(SRCS
    src/V11Geometry.cxx
    src/V1Layer.cxx
    src/Detector.cxx
    src/GeometryHandler.cxx
    src/Digitizer.cxx
    src/DigitizerTask.cxx
    src/DigitWriteoutBuffer.cxx
    src/DigitLayer.cxx
    src/DigitStave.cxx
    )
set(HEADERS
    include/${MODULE_NAME}/Detector.h
    include/${MODULE_NAME}/DigitLayer.h
    include/${MODULE_NAME}/DigitStave.h
    include/${MODULE_NAME}/DigitWriteoutBuffer.h
    include/${MODULE_NAME}/Digitizer.h
    include/${MODULE_NAME}/DigitizerTask.h
    include/${MODULE_NAME}/GeometryHandler.h
    include/${MODULE_NAME}/V1Layer.h
    include/${MODULE_NAME}/V11Geometry.h
    )

Set(LINKDEF src/ITSSimulationLinkDef.h)
Set(LIBRARY_NAME ${MODULE_NAME})
set(BUCKET_NAME its_simulation_bucket)

O2_GENERATE_LIBRARY()
