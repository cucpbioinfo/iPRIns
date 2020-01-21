ExternalProject_Add(htslib
	PREFIX ${CMAKE_BINARY_DIR}/external/htslib
    GIT_REPOSITORY "https://github.com/samtools/htslib.git"
	GIT_TAG 1.9
	UPDATE_COMMAND ""
	BUILD_IN_SOURCE 1
    CONFIGURE_COMMAND autoheader && autoconf && ./configure --prefix=${CMAKE_BINARY_DIR}/external/htslib
	BUILD_COMMAND make lib-static
	INSTALL_COMMAND make install prefix=${CMAKE_BINARY_DIR}/external/htslib
	LOG_DOWNLOAD 0
    LOG_UPDATE 0
    LOG_CONFIGURE 0
    LOG_BUILD 0
    LOG_TEST 0
    LOG_INSTALL 0
)

set(HTSLIB_INCLUDE_DIRS ${CMAKE_BINARY_DIR}/external/htslib/include)
if(APPLE)
    set(HTSLIB_LIBRARIES ${CMAKE_BINARY_DIR}/external/htslib/lib/libhts.dylib)
elseif(UNIX AND NOT APPLE)
    set(HTSLIB_LIBRARIES ${CMAKE_BINARY_DIR}/external/htslib/lib/libhts.so)
    #set(HTSLIB_LIBRARIES ${CMAKE_BINARY_DIR}/external/htslib/lib/libhts.a)
endif()