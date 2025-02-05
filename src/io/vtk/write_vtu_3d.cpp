#ifdef MADVORO_WITH_VTK

#include "write_vtu_3d.hpp"

void MadVoro::IO::write_vtu3d::write_vtu_3d(std::filesystem::path const& file_name,
				std::vector<std::string> const& cell_variable_names,
				std::vector<std::vector<double>> const& cell_variables,
				std::vector<std::string> const& cell_vectors_names,
				std::vector<std::vector<Vector3D>> const& cell_vectors,
				double const time,
				std::size_t cycle,
				Voronoi3D const& tess)
{
	std::vector<Vector3D> const& vertices = tess.GetFacePoints();
	std::size_t const num_vertices = vertices.size();
	std::size_t const num_cells = tess.GetPointNo();

	#ifdef MADVORO_WITH_MPI
		int mpi_rank;
		int mpi_size;
		MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
		MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);
	#endif
	
	vtkNew<vtkUnstructuredGrid> ugrid;

	// ----------- time and cycle
	// see https://www.visitusers.org/index.php?title=Time_and_Cycle_in_VTK_files

	// time
	if(time != std::numeric_limits<double>::max())
	{
		vtkNew<vtkDoubleArray> td;
		td->SetName("TIME");
		td->SetNumberOfTuples(1);
		td->SetTuple1(0, time);
		ugrid->GetFieldData()->AddArray(td);
	}
	
	// cycle
	if(cycle != std::numeric_limits<size_t>::max())
	{
		vtkNew<vtkIntArray> cd;
		cd->SetName("CYCLE");
		cd->SetNumberOfTuples(1);
		cd->SetTuple1(0, cycle);
		ugrid->GetFieldData()->AddArray(cd);
	}

	// ----------- grid vertices coordinates
	// see https://kitware.github.io/vtk-examples/site/Cxx/GeometricObjects/QuadraticHexahedron/
	vtkNew<vtkPoints> points;
	points->SetNumberOfPoints(num_vertices);

	// ----------- cells connectivity
	std::set<size_t> real_vertices;
	std::vector<face_vec > const& cell_faces = tess.GetAllCellFaces();
	std::vector<point_vec > const& points_in_face = tess.GetAllPointsInFace();
	ugrid->Allocate(num_cells);
	std::vector<vtkIdType> point_array_in_cell;
	for(std::size_t cell=0; cell<num_cells; ++cell){
		vtkNew<vtkIdList> faces;
		point_array_in_cell.clear();
		size_t const Nfaces = cell_faces[cell].size();
		//vtkIdType ptIds[currect_cell_num_vertices];
		for(size_t i = 0; i < Nfaces; ++i)
		{
			size_t const face_index = cell_faces[cell][i];
			size_t const Nvertices_in_face = points_in_face[face_index].size();
			faces->InsertNextId(Nvertices_in_face);
			for(size_t j = 0; j < Nvertices_in_face; ++j)
			{
				real_vertices.insert(points_in_face[face_index][j]);
				faces->InsertNextId(points_in_face[face_index][j]);
				point_array_in_cell.push_back(points_in_face[face_index][j]);
			}
		}
		std::sort(point_array_in_cell.begin(), point_array_in_cell.end());
		auto it = std::unique(point_array_in_cell.begin(), point_array_in_cell.end());
		point_array_in_cell = std::vector(point_array_in_cell.begin(), it);
		vtkIdType* ptIds = &point_array_in_cell[0];
		ugrid->InsertNextCell(VTK_POLYHEDRON, point_array_in_cell.size(), ptIds, Nfaces, faces->GetPointer(0));
	}

	Vector3D mid_vertice = 0.5 * (tess.GetBoxCoordinates().first + tess.GetBoxCoordinates().second);
	for(std::size_t p=0; p<num_vertices; ++p){
		if(real_vertices.count(p) > 0)
			points->SetPoint(p, vertices[p].x, vertices[p].y, vertices[p].z);
		else
			points->SetPoint(p, mid_vertice.x, mid_vertice.y, mid_vertice.z);
	}
	ugrid->SetPoints(points);

	// ------------cell centered data
	// see here:
	// https://paraview.paraview.narkive.com/GUqqKIK1/assign-scalars-vectors-to-mesh-points-in-vtkrectilineargrid-c

		// cell scalars data
		for(std::size_t var_index=0; var_index<cell_variable_names.size(); ++var_index){
			vtkNew<vtkDoubleArray> var_data;
			var_data->SetName(cell_variable_names[var_index].c_str());
			var_data->SetNumberOfComponents(1);
			var_data->SetNumberOfValues(num_cells);
			auto const& var = cell_variables[var_index];
			assert(var.size() == num_cells);
			for(std::size_t cell=0; cell<num_cells; ++cell){
				var_data->SetValue(cell, var[cell]);
			}
			ugrid->GetCellData()->AddArray(var_data);
		}

		// write mpi rank as a cell centered data
		#ifdef MADVORO_WITH_MPI
			vtkNew<vtkIntArray> var_mpi_rank;
			var_mpi_rank->SetName(std::string("mpi_rank").c_str());
			var_mpi_rank->SetNumberOfComponents(1);
			var_mpi_rank->SetNumberOfValues(num_cells);
			for(std::size_t cell=0; cell<num_cells; ++cell){
				var_mpi_rank->SetValue(cell, mpi_rank);
			}
			ugrid->GetCellData()->AddArray(var_mpi_rank);
		#endif

		// cell vector data
		for(std::size_t var_index=0; var_index<cell_vectors_names.size(); ++var_index){
			vtkNew<vtkDoubleArray> var_data;
			var_data->SetName(cell_vectors_names[var_index].c_str());
			var_data->SetNumberOfComponents(3);
			var_data->SetNumberOfTuples(num_cells);
			auto const& var = cell_vectors[var_index];
			for(std::size_t cell=0; cell<num_cells; ++cell){
				var_data->SetTuple3(cell, var[cell].x, var[cell].y, var[cell].z);
			}
			ugrid->GetCellData()->AddArray(var_data);
		}


	// ------------ write the 'pieced' .pvtu file
	/*
	Write a 'pvtu' file which holds metadata about all per-cpu 'vtu' files.
	this 'pvtu' file should be opened in paraview/VisIt etc.
	see for example:
	https://gerstrong.github.io/blog/2016/08/20/hacking-vtk-for-parallelisation
	https://www.steinzone.de/wordpress/hacking-vtk-for-parallelisation-mpi-and-c/
	http://rotorbit.blogspot.com/2017/02/how-to-write-vtk-files-in-parallel.html
	https://github.com/rotorbit/RotorbitTutorials/blob/master/rotorbit_WriteVTK.tar.gz
	
	HOWEVER! there is a bug in all of those links: not all .vtu files paths are listed
	in the resulting .pvtu file. Instead, only the root rank's .vtu file is listed.
	this can be solved by adding manually the .vtu file paths of all processes to the .pvtu file.
	
	Instead, a solution to this problem is obtained by passing the MPI communicator 
	to the VTK library, as offered here:
	https://discourse.vtk.org/t/distributed-i-o-with-vtkxmlpunstructuredgridwriter/7418/3
	which gives a link to these lines libmesh implementation, which I use here:
	https://github.com/libMesh/libmesh/blob/484c8652977b6504e1613633b4c7d87297f94957/src/mesh/vtk_io.C#L50-L54
	https://github.com/libMesh/libmesh/blob/484c8652977b6504e1613633b4c7d87297f94957/src/mesh/vtk_io.C#L275-L286
	but it only works if VTK library was compiled with MPI.
	*/
	#ifdef MADVORO_WITH_MPI
		vtkNew<vtkXMLPUnstructuredGridWriter> pwriter;

		// Set VTK library to use the same MPI communicator as we do
		vtkNew<vtkMPICommunicator> vtk_comm;
		MPI_Comm mpi_comm(MPI_COMM_WORLD);
		vtkMPICommunicatorOpaqueComm vtk_opaque_comm(&mpi_comm);
		vtk_comm ->InitializeExternal(&vtk_opaque_comm);
		vtkNew<vtkMPIController> vtk_mpi_ctrl;
		vtk_mpi_ctrl->SetCommunicator(vtk_comm );
		pwriter->SetController(vtk_mpi_ctrl);

		// Tell the writer how many partitions ('pieces') exist and on which processor we are currently
		pwriter->SetNumberOfPieces(mpi_size);
		pwriter->SetStartPiece(mpi_rank);
		pwriter->SetEndPiece(mpi_rank);
		pwriter->SetUseSubdirectory(true);

		// ----- write file to disk
		std::filesystem::path pname(file_name);
		pname.replace_extension("pvtu");
		//printf("%d/%d writing vtu file '%s'\n", mpi_rank+1, mpi_size, file_pvtu.c_str());
		pwriter->SetFileName(pname.c_str());
		pwriter->SetInputData(ugrid);
		pwriter->Write();
	#else
		// ------------ write the vtu file to disk
		vtkNew<vtkXMLUnstructuredGridWriter> writer;
		writer->SetCompressionLevel(9);
		writer->SetFileName(file_name.c_str());
		writer->SetInputData(ugrid);
		writer->Write();
	#endif
}

void write_vtu_3d_points(std::filesystem::path const& file_name,
			std::vector<std::string> const& point_variable_names,
			std::vector<std::vector<double>> const& point_variables,
			std::vector<std::string> const& point_vectors_names,
			std::vector<std::vector<Vector3D>> const& point_vectors,
			double const time,
			std::size_t cycle,
			Voronoi3D const& tess){
	std::vector<Vector3D> mesh_points = tess.getMeshPoints();
	std::size_t const num_points = tess.GetPointNo();
	mesh_points.resize(num_points);
	std::vector<Vector3D> cm_points = tess.GetAllCM();
	cm_points.resize(num_points);

	#ifdef MADVORO_WITH_MPI
		int mpi_rank;
		int mpi_size;
		MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
		MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);
	#endif
	
	vtkNew<vtkUnstructuredGrid> ugrid;

	// ----------- time and cycle
	// see https://www.visitusers.org/index.php?title=Time_and_Cycle_in_VTK_files

	// time
	vtkNew<vtkDoubleArray> td;
	td->SetName("TIME");
	td->SetNumberOfTuples(1);
	td->SetTuple1(0, time);
	ugrid->GetFieldData()->AddArray(td);
	
	// cycle
	vtkNew<vtkIntArray> cd;
	cd->SetName("CYCLE");
	cd->SetNumberOfTuples(1);
	cd->SetTuple1(0, cycle);
	ugrid->GetFieldData()->AddArray(cd);

	// ----------- grid vertices coordinates
	// see https://kitware.github.io/vtk-examples/site/Cxx/GeometricObjects/QuadraticHexahedron/
	vtkNew<vtkPoints> points;
	points->SetNumberOfPoints(num_points);
	for(std::size_t p=0; p<num_points; ++p){
		points->SetPoint(p, mesh_points[p].x, mesh_points[p].y, mesh_points[p].z);
	}
	ugrid->SetPoints(points);

	
	// ------------cell centered data
	// see here:
	// https://paraview.paraview.narkive.com/GUqqKIK1/assign-scalars-vectors-to-mesh-points-in-vtkrectilineargrid-c

	// cell scalars data
	for(std::size_t var_index=0; var_index<point_variable_names.size(); ++var_index){
		vtkNew<vtkDoubleArray> var_data;
		var_data->SetName(point_variable_names[var_index].c_str());
		var_data->SetNumberOfComponents(1);
		var_data->SetNumberOfValues(num_points);
		auto const& var = point_variables[var_index];
		// assert(var.size() == num_cells);
		for(std::size_t cell=0; cell<num_points; ++cell){
			var_data->SetValue(cell, var[cell]);
		}
		ugrid->GetPointData()->AddArray(var_data);
	}

	// write mpi rank as a cell centered data
	#ifdef MADVORO_WITH_MPI
		vtkNew<vtkIntArray> var_mpi_rank;
		var_mpi_rank->SetName(std::string("mpi_rank").c_str());
		var_mpi_rank->SetNumberOfComponents(1);
		var_mpi_rank->SetNumberOfValues(num_points);
		for(std::size_t cell=0; cell<num_points; ++cell){
			var_mpi_rank->SetValue(cell, mpi_rank);
		}
		ugrid->GetPointData()->AddArray(var_mpi_rank);
	#endif

	// point vector data
	for(std::size_t var_index=0; var_index<point_vectors_names.size(); ++var_index){
		vtkNew<vtkDoubleArray> var_data;
		var_data->SetName(point_vectors_names[var_index].c_str());
		var_data->SetNumberOfComponents(3);
		var_data->SetNumberOfTuples(num_points);
		auto const& var = point_vectors[var_index];
		for(std::size_t cell=0; cell<num_points; ++cell){
			var_data->SetTuple3(cell, var[cell].x, var[cell].y, var[cell].z);
		}
		ugrid->GetPointData()->AddArray(var_data);
	}
	// Write CM
	{
		vtkNew<vtkDoubleArray> var_data;
		var_data->SetName(std::string("CM").c_str());
		var_data->SetNumberOfComponents(3);
		var_data->SetNumberOfTuples(num_points);
		auto const& var = cm_points;
		for(std::size_t cell=0; cell<num_points; ++cell){
			var_data->SetTuple3(cell, var[cell].x, var[cell].y, var[cell].z);
		}
		ugrid->GetPointData()->AddArray(var_data);
	}


	// ------------ write the 'pieced' .pvtu file
	/*
	Write a 'pvtu' file which holds metadata about all per-cpu 'vtu' files.
	this 'pvtu' file should be opened in paraview/VisIt etc.
	see for example:
	https://gerstrong.github.io/blog/2016/08/20/hacking-vtk-for-parallelisation
	https://www.steinzone.de/wordpress/hacking-vtk-for-parallelisation-mpi-and-c/
	http://rotorbit.blogspot.com/2017/02/how-to-write-vtk-files-in-parallel.html
	https://github.com/rotorbit/RotorbitTutorials/blob/master/rotorbit_WriteVTK.tar.gz
	
	HOWEVER! there is a bug in all of those links: not all .vtu files paths are listed
	in the resulting .pvtu file. Instead, only the root rank's .vtu file is listed.
	this can be solved by adding manually the .vtu file paths of all processes to the .pvtu file.
	
	Instead, a solution to this problem is obtained by passing the MPI communicator 
	to the VTK library, as offered here:
	https://discourse.vtk.org/t/distributed-i-o-with-vtkxmlpunstructuredgridwriter/7418/3
	which gives a link to these lines libmesh implementation, which I use here:
	https://github.com/libMesh/libmesh/blob/484c8652977b6504e1613633b4c7d87297f94957/src/mesh/vtk_io.C#L50-L54
	https://github.com/libMesh/libmesh/blob/484c8652977b6504e1613633b4c7d87297f94957/src/mesh/vtk_io.C#L275-L286
	but it only works if VTK library was compiled with MPI.
	*/
	#ifdef MADVORO_WITH_MPI
		vtkNew<vtkXMLPUnstructuredGridWriter> pwriter;

		// Set VTK library to use the same MPI communicator as we do
		vtkNew<vtkMPICommunicator> vtk_comm;
		MPI_Comm mpi_comm(MPI_COMM_WORLD);
		vtkMPICommunicatorOpaqueComm vtk_opaque_comm(&mpi_comm);
		vtk_comm ->InitializeExternal(&vtk_opaque_comm);
		vtkNew<vtkMPIController> vtk_mpi_ctrl;
		vtk_mpi_ctrl->SetCommunicator(vtk_comm );
		pwriter->SetController(vtk_mpi_ctrl);

		// Tell the writer how many partitions ('pieces') exist and on which processor we are currently
		pwriter->SetNumberOfPieces(mpi_size);
		pwriter->SetStartPiece(mpi_rank);
		pwriter->SetEndPiece(mpi_rank);
		pwriter->SetUseSubdirectory(true);

		// ----- write file to disk
		std::filesystem::path pname(file_name);
		pname.replace_extension("pvtu");
		//printf("%d/%d writing vtu file '%s'\n", mpi_rank+1, mpi_size, file_pvtu.c_str());
		pwriter->SetFileName(pname.c_str());
		pwriter->SetInputData(ugrid);
		pwriter->Write();
	#else
		// ------------ write the vtu file to disk
		vtkNew<vtkXMLUnstructuredGridWriter> writer;
		writer->SetCompressionLevel(9);
		writer->SetFileName(file_name.c_str());
		writer->SetInputData(ugrid);
		writer->Write();
	#endif
}

#endif // MADVORO_WITH_VTK