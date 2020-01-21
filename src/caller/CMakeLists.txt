file(GLOB SOURCES "*.cc")

add_library(caller
        STATIC
        ${SOURCES}
)

include_directories(${PROJECT_SOURCE_DIR})
include_directories(${PROJECT_SOURCE_DIR}/src)

if(APPLE)
    target_link_libraries(
        caller
        fasta
        ${HTSLIB_LIBRARIES}
        ${ZLIB_LIBRARIES}
        ${TBB_IMPORTED_TARGETS}
    )
elseif(UNIX AND NOT APPLE)
    target_link_libraries(
        caller
        fasta
        ${HTSLIB_LIBRARIES}
        ${ZLIB_LIBRARIES}
        ${TBB_IMPORTED_TARGETS}
    )
endif()
