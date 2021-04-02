set(MODULE_NAME "O2QcMuonChambers")


# ---- Files ----

set(SRCS
  src/GlobalHistogram.cxx
  src/PedestalsTask.cxx
  src/PhysicsTaskDigits.cxx
  src/PhysicsTaskPreClusters.cxx
  src/PedestalsCheck.cxx
  src/PhysicsCheck.cxx
  src/MCHPostProcessing.cxx
  src/TH1MCHReductor.cxx
  src/TH2MCHReductor.cxx
  src/GlobalHistogramReductor.cxx
)

set(HEADERS
  include/MCH/GlobalHistogram.h
  include/MCH/PedestalsTask.h
  include/MCH/PhysicsTaskDigits.h
  include/MCH/PhysicsTaskPreClusters.h
  include/MCH/PedestalsCheck.h
  include/MCH/PhysicsCheck.h
  include/MCH/MCHPostProcessing.h
  include/MCH/TH1MCHReductor.h
  include/MCH/TH2MCHReductor.h
  include/MCH/GlobalHistogramReductor.h
  include/MCH/CustomMergeableTH2Quotient.h
)

# ---- Library ----

add_library(${MODULE_NAME} SHARED ${SRCS})

target_include_directories(
  ${MODULE_NAME}
  PUBLIC $<INSTALL_INTERFACE:include> $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
  PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/src
)

target_link_libraries(${MODULE_NAME} PUBLIC O2QualityControl O2::CommonDataFormat O2::GPUCommon
        $<TARGET_NAME_IF_EXISTS:O2::MCHMappingFactory> O2::MCHMappingImpl4 O2::MCHMappingSegContour O2::MCHBase O2::DataFormatsMCH O2::MCHRawDecoder O2::MCHCalibration O2::MCHPreClustering O2::Mergers)

target_compile_definitions(${MODULE_NAME} PRIVATE $<$<TARGET_EXISTS:O2::MCHMappingFactory>:MCH_HAS_MAPPING_FACTORY>)

# Digit.h is moving from MCHBase to DataFormatsMCH : let's handle both
# gracefully for the moment...
get_target_property(O2_INCLUDE_DIRS O2::CommonDataFormat INTERFACE_INCLUDE_DIRECTORIES)
get_target_property(ROOT_INCLUDE_DIRS ROOT::Core INTERFACE_INCLUDE_DIRECTORIES)

set(CMAKE_REQUIRED_INCLUDES ${O2_INCLUDE_DIRS} ${ROOT_INCLUDE_DIRS})
check_include_file_cxx("DataFormatsMCH/Digit.h" HAVE_DIGIT_IN_DATAFORMATS)

if(HAVE_DIGIT_IN_DATAFORMATS)
    target_compile_definitions(${MODULE_NAME} PRIVATE HAVE_DIGIT_IN_DATAFORMATS)
endif()

install(
        TARGETS ${MODULE_NAME}
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
        ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
        RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
)

# ---- ROOT dictionary ----

add_root_dictionary(${MODULE_NAME}
                    HEADERS include/MCH/GlobalHistogram.h
                            include/MCH/PedestalsTask.h
                            include/MCH/PhysicsTaskDigits.h
                            include/MCH/PhysicsTaskPreClusters.h
                            include/MCH/PedestalsCheck.h
                            include/MCH/PhysicsCheck.h
                            include/MCH/MCHPostProcessing.h
                            include/MCH/TH1MCHReductor.h
                            include/MCH/TH2MCHReductor.h
                            include/MCH/GlobalHistogramReductor.h
                            include/MCH/CustomMergeableTH2Quotient.h
                    LINKDEF include/MCH/LinkDef.h)
                    
# ---- Executables ----

set(EXE_SRCS
     src/TestMerging.cxx
     src/MergeMCH.cxx)
set(EXE_NAMES
     o2-qc-run-test-mergin
     o2-qc-run-mch-merger)

list(LENGTH EXE_SRCS count)
math(EXPR count "${count}-1")
foreach(i RANGE ${count})
  list(GET EXE_SRCS ${i} src)
  list(GET EXE_NAMES ${i} name)
  add_executable(${name} ${src})
  target_link_libraries(${name} PRIVATE QualityControl QcMuonChambers)
endforeach()

install(
  TARGETS o2-qc-run-test-mergin
  o2-qc-run-mch-merger
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
)

# ---- Tests ----

set(
  TEST_SRCS
)

foreach(test ${TEST_SRCS})
  get_filename_component(test_name ${test} NAME)
  string(REGEX REPLACE ".cxx" "" test_name ${test_name})

  add_executable(${test_name} ${test})
  target_link_libraries(${test_name} PRIVATE ${MODULE_NAME} Boost::unit_test_framework)
  add_test(NAME ${test_name} COMMAND ${test_name})
  set_tests_properties(${test_name} PROPERTIES TIMEOUT 60)
endforeach()
