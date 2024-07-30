
## Utils

set(SOURCES
    common.cc
    fileutils.cc
)

add_library(utils ${SOURCES})
target_link_libraries(utils -lz)