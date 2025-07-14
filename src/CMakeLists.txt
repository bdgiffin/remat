set(REMAT_dependencies
#    Boost
    PARENT_SCOPE
    )

set(REMAT_include_dirs
#    ${Boost_INCLUDE_DIR}
    PARENT_SCOPE
    )

set(REMAT_sources
    ${CMAKE_CURRENT_LIST_DIR}/REMAT_API.cpp
    )

set(REMAT_python
    ${CMAKE_CURRENT_LIST_DIR}/__init__.py
    ${CMAKE_CURRENT_LIST_DIR}/REMAT.py
    ${CMAKE_CURRENT_LIST_DIR}/MeshUtils/FaceSet.py
    ${CMAKE_CURRENT_LIST_DIR}/MeshUtils/Geometry.py
    ${CMAKE_CURRENT_LIST_DIR}/MeshUtils/GeometryFactory.py
    ${CMAKE_CURRENT_LIST_DIR}/MeshUtils/Material.py
    ${CMAKE_CURRENT_LIST_DIR}/MeshUtils/Model.py
    ${CMAKE_CURRENT_LIST_DIR}/MeshUtils/MeshUtils.py
    ${CMAKE_CURRENT_LIST_DIR}/MeshUtils/NodeSet.py
    ${CMAKE_CURRENT_LIST_DIR}/MeshUtils/Part.py
    ${CMAKE_CURRENT_LIST_DIR}/MeshUtils/Selector.py
    PARENT_SCOPE
    )

set(REMAT_headers
    ${CMAKE_CURRENT_LIST_DIR}/Dual.h
    ${CMAKE_CURRENT_LIST_DIR}/ContactInteraction.h
    ${CMAKE_CURRENT_LIST_DIR}/Element.h
    ${CMAKE_CURRENT_LIST_DIR}/Fixed.h
    ${CMAKE_CURRENT_LIST_DIR}/Integrator.h
    ${CMAKE_CURRENT_LIST_DIR}/Material.h
    ${CMAKE_CURRENT_LIST_DIR}/Matrix.h
    ${CMAKE_CURRENT_LIST_DIR}/Parameters.h
    ${CMAKE_CURRENT_LIST_DIR}/Rational.h
    ${CMAKE_CURRENT_LIST_DIR}/System.h
    ${CMAKE_CURRENT_LIST_DIR}/Truss.h
    ${CMAKE_CURRENT_LIST_DIR}/UniaxialMaterial.h
    ${CMAKE_CURRENT_LIST_DIR}/types.h
    PARENT_SCOPE
    )
  
blt_add_library(
    NAME        REMAT
    SOURCES     ${REMAT_sources}
    HEADERS     ${REMAT_headers}
    INCLUDES    ${REMAT_include_dirs}
    DEPENDS_ON  ${REMAT_dependencies}
    SHARED      TRUE
    )
set_target_properties(REMAT PROPERTIES SUFFIX ".so")
#target_link_libraries(REMAT LINK_PUBLIC ${Boost_LIBRARIES})
