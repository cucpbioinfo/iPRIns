file(GLOB SOURCES "*.cc")


add_library(fasta
     STATIC
    ${SOURCES}
)