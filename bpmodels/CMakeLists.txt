
## Branch Predictor Model library

set(SOURCES
    base_predictor.cc

    tage/tage.cc
    tage/tage_scl.cc
    llbp/llbp.cc
)

add_library(bpmodels STATIC ${SOURCES})

target_include_directories(bpmodels PUBLIC "${BASEPATH}")
target_link_libraries(bpmodels PUBLIC utils)
# INCLUDE_DIRECTORIES("${BASEPATH}")