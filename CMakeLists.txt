cmake_minimum_required(VERSION 2.8.12.2)

project(cfd_3d)

set(CMAKE_CXX_COMPILER mpic++)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

#set(MPI_INCLUDE "/usr/include/mpich")
#set(MPI_LIB "/usr/lib/x86_64-linux-gnu")


set(CONTRIB_DIR "${CMAKE_SOURCE_DIR}/contrib")


set(LOG_INCLUDE "${CONTRIB_DIR}/log4cpp/include")
set(LOG_LIB "${CONTRIB_DIR}/log4cpp/lib")

set(HYPRE_INCLUDE "${CONTRIB_DIR}/hypre/include")
set(HYPRE_LIB "${CONTRIB_DIR}/hypre/lib")


include_directories(
	SYSTEM ${MPI_INCLUDE_PATH}
	#${MPI_INCLUDE}
	${LOG_INCLUDE}
	${HYPRE_INCLUDE}
	src
	src/mesh
	src/tinyxml
	src/iterators
	src/mesh_properties
	src/mesh_readers
	src/solvers
	src/solvers/MatrixSolver
)

link_directories(
	#${MPI_LIB}
	${LOG_LIB}
	${HYPRE_LIB}
)


set(METHODS_SRC
	src/solvers/Fvm_tvd_implicit.cpp
	src/solvers/Solver.cpp
)

set(SOLVERS_SRC
	src/solvers/MatrixSolver/CSR.cpp
	src/solvers/MatrixSolver/MatrixSolver.cpp
	src/solvers/MatrixSolver/SolverZeidel.cpp
	src/solvers/MatrixSolver/SolverJacobi.cpp
	src/solvers/MatrixSolver/SolverHypre.cpp
	src/solvers/MatrixSolver/SolverHypreGmres.cpp
	src/solvers/MatrixSolver/SolverHypreFlexGmres.cpp
		src/solvers/MatrixSolver/SolverHypreBoomerAmg.cpp )

set(MESH_SRC
	src/mesh/Grid.cpp
	src/mesh/Mesh.cpp
)

set(PROPS_SRC
	src/mesh_properties/BndFaceTemperature.cpp
	src/mesh_properties/CellTemperature.cpp
	src/mesh_properties/CellFluidDynamicsProps.cpp
)

set(READERS_SRC
	src/mesh_readers/MeshReader.cpp
	src/mesh_readers/MeshReaderUnv.cpp
)

set(TINYXML_SRC
	src/tinyxml/tinystr.cpp
	src/tinyxml/tinyxml.cpp
	src/tinyxml/tinyxmlerror.cpp
	src/tinyxml/tinyxmlparser.cpp
)


add_executable(
	cfd-3d
	src/main.cpp
	src/global.cpp
	${METHODS_SRC}
	${SOLVERS_SRC}
	${MESH_SRC}
	${PROPS_SRC}
	${READERS_SRC}
	${TINYXML_SRC}
)

target_link_libraries(cfd-3d ${MPI_LIBRARIES} HYPRE log4cpp pthread)
install(TARGETS cfd-3d RUNTIME DESTINATION .)
