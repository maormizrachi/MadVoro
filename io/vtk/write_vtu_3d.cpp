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

	if(time != std::numeric_limits<double>::max())
	{
		vtkNew<vtkDoubleArray> td;
		td->SetName("TIME");
		td->SetNumberOfTuples(1);
		td->SetTuple1(0, time);
		ugrid->GetFieldData()->AddArray(td);
	}
	
	if(cycle != std::numeric_limits<size_t>::max())
	{
		vtkNew<vtkIntArray> cd;
		cd->SetName("CYCLE");
		cd->SetNumberOfTuples(1);
		cd->SetTuple1(0, cycle);
		ugrid->GetFieldData()->AddArray(cd);
	}

	vtkNew<vtkPoints> points;
	points->SetNumberOfPoints(num_vertices);

	std::set<size_t> real_vertices;
	std::vector<face_vec > const& cell_faces = tess.GetAllCellFaces();
	std::vector<point_vec > const& points_in_face = tess.GetAllPointsInFace();
	ugrid->Allocate(num_cells);
	std::vector<vtkIdType> point_array_in_cell;
	for(std::size_t cell=0; cell<num_cells; ++cell){
		vtkNew<vtkIdList> faces;
		point_array_in_cell.clear();
		size_t const Nfaces = cell_faces[cell].size();
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

	#ifdef MADVORO_WITH_MPI
		vtkNew<vtkXMLPUnstructuredGridWriter> pwriter;
		vtkNew<vtkMPICommunicator> vtk_comm;
		MPI_Comm mpi_comm(MPI_COMM_WORLD);
		vtkMPICommunicatorOpaqueComm vtk_opaque_comm(&mpi_comm);
		vtk_comm->InitializeExternal(&vtk_opaque_comm);
		vtkNew<vtkMPIController> vtk_mpi_ctrl;
		vtk_mpi_ctrl->SetCommunicator(vtk_comm);
		pwriter->SetController(vtk_mpi_ctrl);
		pwriter->SetNumberOfPieces(mpi_size);
		pwriter->SetStartPiece(mpi_rank);
		pwriter->SetEndPiece(mpi_rank);
		pwriter->SetUseSubdirectory(true);
		std::filesystem::path pname(file_name);
		pname.replace_extension("pvtu");
		pwriter->SetFileName(pname.c_str());
		pwriter->SetInputData(ugrid);
		pwriter->Write();
	#else
		vtkNew<vtkXMLUnstructuredGridWriter> writer;
		writer->SetCompressionLevel(9);
		writer->SetFileName(file_name.c_str());
		writer->SetInputData(ugrid);
		writer->Write();
	#endif
}

#endif // MADVORO_WITH_VTK
