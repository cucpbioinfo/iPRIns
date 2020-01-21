find_package( LibLZMA REQUIRED )
if ( LIBLZMA_FOUND )
	message( ${LIBLZMA_INCLUDE_DIRS} )
	message( ${LIBLZMA_LIBRARIES} )
else()

endif( LIBLZMA_FOUND )
