# Build external library
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -I/usr/include")


find_package( ZLIB REQUIRED )
if (ZLIB_FOUND )
	message( ${ZLIB_INCLUDE_DIRS} )
	message( ${ZLIB_LIBRARIES} )
else()
ExternalProject_Add(zlib
	PREFIX ${CMAKE_BINARY_DIR}/external/zlib
	GIT_REPOSITORY "https://github.com/madler/zlib.git"
	GIT_TAG "v1.2.11"
	UPDATE_COMMAND ""
	BUILD_IN_SOURCE 1
	CONFIGURE_COMMAND ${CMAKE_BINARY_DIR}/external/zlib/src/zlib/configure --prefix=${CMAKE_BINARY_DIR}/external/zlib
	INSTALL_DIR ${CMAKE_BINARY_DIR}/external/zlib
	LOG_DOWNLOAD 0
    LOG_UPDATE 0
    LOG_CONFIGURE 0
    LOG_BUILD 0
    LOG_TEST 0
    LOG_INSTALL 0
	)
	
	set(ZLIB_INCLUDE_DIRS ${CMAKE_BINARY_DIR}/external/zlib/src/zlib/include)
	set(ZLIB_LIBRARIES ${CMAKE_BINARY_DIR}/external/zlib/lib/libz.so.1)
endif( ZLIB_FOUND )
