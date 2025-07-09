set(REMAT_dependencies
    PARENT_SCOPE
    )

set(REMAT_sources
    ${CMAKE_CURRENT_LIST_DIR}/REMAT_API.cpp
    )

set(REMAT_python
    ${CMAKE_CURRENT_LIST_DIR}/__init__.py
    ${CMAKE_CURRENT_LIST_DIR}/REMAT.py
    PARENT_SCOPE
    )

set(REMAT_headers
    ${CMAKE_CURRENT_LIST_DIR}/Dual.h
    ${CMAKE_CURRENT_LIST_DIR}/Element.h
    ${CMAKE_CURRENT_LIST_DIR}/Fixed.h
    ${CMAKE_CURRENT_LIST_DIR}/Integrator.h
    ${CMAKE_CURRENT_LIST_DIR}/Material.h
    ${CMAKE_CURRENT_LIST_DIR}/Matrix.h
    ${CMAKE_CURRENT_LIST_DIR}/Parameters.h
    ${CMAKE_CURRENT_LIST_DIR}/Rational.h
    ${CMAKE_CURRENT_LIST_DIR}/System.h
    ${CMAKE_CURRENT_LIST_DIR}/types.h
    PARENT_SCOPE
    )
  
blt_add_library(
    NAME        REMAT
    SOURCES     ${REMAT_sources}
    HEADERS     ${REMAT_headers}
    DEPENDS_ON  ${REMAT_dependencies}
    SHARED      TRUE
    )
set_target_properties(REMAT PROPERTIES SUFFIX ".so")
