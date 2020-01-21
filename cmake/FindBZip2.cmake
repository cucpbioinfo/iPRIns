find_package( BZip2 REQUIRED )
if ( BZIP2_FOUND )
	message( ${BZIP2_INCLUDE_DIR} )
	message( ${BZIP2_LIBRARIES} )
else()
	ExternalProject_Add(bzip2
	PREFIX ${CMAKE_BINARY_DIR}/external/bzip2
  	URL https://jaist.dl.sourceforge.net/project/bzip2/bzip2-1.0.6.tar.gz
  	UPDATE_COMMAND ""
	BUILD_IN_SOURCE TRUE
	UPDATE_COMMAND ""
	BUILD_IN_SOURCE 1
    CONFIGURE_COMMAND ""
	CONFIGURE_COMMAND ""
	BUILD_COMMAND make lib-static
	INSTALL_COMMAND make install prefix=${CMAKE_BINARY_DIR}/external/bzip2
	LOG_DOWNLOAD 0
    LOG_UPDATE 0
    LOG_CONFIGURE 0
    LOG_BUILD 0
    LOG_TEST 0
    LOG_INSTALL 0
	)

endif( BZIP2_FOUND )