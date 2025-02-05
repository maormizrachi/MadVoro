SUBDIRS = include
AM_CPPFLAGS = -I$(top_srcdir)/include

lib_LTLIBRARIES = libmadvoro.la
libmadvoro_la_LDFLAGS = -version-info 1 -no-undefined

libmadvoro_la_SOURCES = delaunay/Delaunay3D.cpp \
                        elementary/Face.cpp \
                        elementary/Face3D.cpp \
                        elementary/Tetrahedron.cpp \
                        elementary/Point3D.cpp \
                        elementary/Vector3D.cpp \
                        exception/InvalidArgumentException.cpp \
                        exception/SizeException.cpp \
                        exception/MadVoroException.cpp \
                        geometry/Intersections.cpp \
                        hilbert/HilbertOrder3D.cpp \
                        hilbert/HilbertOrder3D_Utils.cpp \
                        hilbert/rectangular/HilbertRectangularConvertor3D.cpp \
                        io/io3D.cpp \
                        io/simple_io.cpp \
                        io/hdf5/hdf5_utils.cpp \
                        io/hdf5/WritePointsHDF5.cpp \
                        io/hdf5/WriteVoronoiHDF5.cpp \
                        io/vtk/write_vtu_3d.cpp \
                        io/vtk/WriteVoronoiVTK.cpp \
                        utils/Predicates3D.cpp \
                        utils/utils.cpp \
                        voronoi/Voronoi3D.cpp \
                        voronoi/pointsManager/HilbertPointsManager.cpp \
                        voronoi/pointsManager/PointsManager.cpp

dist_noinst_HEADERS = delaunay/Delaunay3D.hpp \
                        ds/BinaryTree.hpp \
                        ds/DistributedOctTree/DistributedOctTree.hpp \
                        ds/DistributedOctTree/RankedValue.hpp \
                        ds/GroupRangeTree/GroupRangeTree.hpp \
                        ds/GroupTree/GroupTree.hpp \
                        ds/KDTree/KDTree.hpp \
                        ds/OctTree/OctTree.hpp \
                        ds/RangeTree/RangeTree.hpp \
                        ds/utils/BoundingBox.hpp \
                        ds/utils/geometry.hpp \
                        ds/utils/raw_type.h \
                        ds/utils/Sphere.hpp \
                        elementary/Face3D.hpp \
                        elementary/Face.hpp \
                        elementary/Mat33.hpp \
                        elementary/Point3D.hpp \
                        elementary/Tetrahedron.hpp \
                        elementary/Vector3D.hpp \
                        environment/CurveEnvAgent.hpp \
                        environment/EnvironmentAgent.h \
                        environment/hilbert/DistributedOctEnvAgent.hpp \
                        environment/hilbert/HilbertCurveEnvAgent.hpp \
                        environment/hilbert/HilbertEnvAgent.hpp \
                        environment/hilbert/HilbertTreeEnvAgent.hpp \
                        exception/InvalidArgumentException.hpp \
                        exception/MadVoroException.hpp \
                        exception/SizeException.hpp \
                        geometry/Intersections.hpp \
                        hilbert/HilbertConvertor3D.hpp \
                        hilbert/HilbertOrder3D.hpp \
                        hilbert/HilbertOrder3D_Utils.hpp \
                        hilbert/hilbertTypes.h \
                        hilbert/ordinary/HilbertOrdinaryConvertor3D.hpp \
                        hilbert/rectangular/HilbertRectangularConvertor3D.hpp \
                        hilbert/rectangular/HilbertRectangularTree3D.hpp \
                        io/hdf5/hdf5_utils.hpp \
                        io/hdf5/WritePointsHDF5.hpp \
                        io/hdf5/WriteVoronoiHDF5.hpp \
                        io/io3D.hpp \
                        io/simple_io.hpp \
                        io/vtk/WriteVoronoiVTK.hpp \
                        io/vtk/write_vtu_3d.hpp \
                        mpi/balance/balance.hpp \
                        mpi/exchange/exchange.hpp \
                        mpi/serialize/mpi_commands.hpp \
                        mpi/serialize/Serializable.hpp \
                        mpi/serialize/Serializer.hpp \
                        queryAgent/AnswerAgent.hpp \
                        queryAgent/BusyWaitQueryAgent.hpp \
                        queryAgent/QueryAgent.hpp \
                        queryAgent/QueryData.hpp \
                        queryAgent/TalkAgent.hpp \
                        queryAgent/WaitUntilAnsweredQueryAgent.hpp \
                        range/BigRangeAgent.hpp \
                        range/finders/BruteForce.hpp \
                        range/finders/GroupRangeTree.hpp \
                        range/finders/KDTree.hpp \
                        range/finders/OctTree.hpp \
                        range/finders/RangeFinder.hpp \
                        range/finders/RangeTree.hpp \
                        range/finders/utils/IndexedVector.hpp \
                        range/RangeQueryData.h \
                        range/SentPointsContainer.hpp \
                        range/SmallRangeAgent.hpp \
                        utils/compiler.h \
                        utils/point/3DPoint.hpp \
                        utils/Predicates3D.hpp \
                        utils/print/all.hpp \
                        utils/print/printTuples.hpp \
                        utils/print/printVectors.hpp \
                        utils/utils.hpp \
                        voronoi/pointsManager/HilbertPointsManager.hpp \
                        voronoi/pointsManager/PointsManager.hpp \
                        voronoi/Voronoi3D.hpp

get_sources:
	echo $(libmadvoro_la_SOURCES)
