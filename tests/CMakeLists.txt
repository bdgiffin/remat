#set(test_dependencies REMAT ${REMAT_dependencies} gtest)
set(test_dependencies ${REMAT_dependencies} gtest)
set(test_include_dirs ${REMAT_include_dirs} ../src)

blt_add_executable(NAME       test_Fixed
                   SOURCES    ${CMAKE_CURRENT_LIST_DIR}/test_Fixed.cpp
                   INCLUDES   ${test_include_dirs}
                   DEPENDS_ON ${test_dependencies}
		   )
		 
blt_add_test(NAME    test_Fixed
	     COMMAND test_Fixed)

blt_add_executable(NAME       test_Rational
                   SOURCES    ${CMAKE_CURRENT_LIST_DIR}/test_Rational.cpp
                   INCLUDES   ${test_include_dirs}
                   DEPENDS_ON ${test_dependencies}
		   )
		 
blt_add_test(NAME    test_Rational
	     COMMAND test_Rational)

blt_add_executable(NAME       test_Dual
                   SOURCES    ${CMAKE_CURRENT_LIST_DIR}/test_Dual.cpp
                   INCLUDES   ${test_include_dirs}
                   DEPENDS_ON ${test_dependencies}
		   )
		 
blt_add_test(NAME    test_Dual
	     COMMAND test_Dual)

blt_add_executable(NAME       test_Matrix
                   SOURCES    ${CMAKE_CURRENT_LIST_DIR}/test_Matrix.cpp
                   INCLUDES   ${test_include_dirs}
                   DEPENDS_ON ${test_dependencies}
		   )
		 
blt_add_test(NAME    test_Matrix
	     COMMAND test_Matrix)

blt_add_executable(NAME       test_Integrator
                   SOURCES    ${CMAKE_CURRENT_LIST_DIR}/test_Integrator.cpp
                   INCLUDES   ${test_include_dirs}
                   DEPENDS_ON ${test_dependencies}
		   )
		 
blt_add_test(NAME    test_Integrator
	     COMMAND test_Integrator)

blt_add_executable(NAME       test_System
                   SOURCES    ${CMAKE_CURRENT_LIST_DIR}/test_System.cpp
                   INCLUDES   ${test_include_dirs}
                   DEPENDS_ON ${test_dependencies}
		   )
		 
blt_add_test(NAME    test_System
	     COMMAND test_System)

