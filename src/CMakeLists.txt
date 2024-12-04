set(REMAT_dependencies
    PARENT_SCOPE
    )

set(REMAT_sources
    )

set(REMAT_headers
    ${CMAKE_CURRENT_LIST_DIR}/types.h
    ${CMAKE_CURRENT_LIST_DIR}/Fixed.h
    ${CMAKE_CURRENT_LIST_DIR}/Rational.h
    ${CMAKE_CURRENT_LIST_DIR}/Dual.h
    PARENT_SCOPE
    )
  
#blt_add_library(
#    NAME        REMAT
#    SOURCES     ${REMAT_sources}
#    HEADERS     ${REMAT_headers}
#    DEPENDS_ON  ${REMAT_dependencies}
#    SHARED      FALSE
#    )
