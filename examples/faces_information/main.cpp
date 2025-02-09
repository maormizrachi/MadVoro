#include <madvoro/Voronoi3D.hpp>

int main(int argc, char *argv[])
{
    MadVoro::Vector3D ll(-10, -10, -10); // lower left corner
    MadVoro::Vector3D ur(10, 10, 10); // upper right corner

    MadVoro::Voronoi3D diagram(ll, ur);

    std::vector<MadVoro::Vector3D> points {{ -0.115, 4.484, 0.688 },
                                            { -6.121, 6.335, -9.093 },
                                            { -5.529, -7.464, 2.459 },
                                            { -3.448, -7.794, -6.939 },
                                            { -3.231, 8.213, -7.829 },
                                            { 5.204, 9.943, 6.175 },
                                            { -7.018, 7.155, -0.226 },
                                            { 0.567, 2.263, -0.962 },
                                            { 9.998, -3.030, 7.772 },
                                            { 3.811, 6.672, -8.231 },
                                            { 4.252, -8.460, -2.693 },
                                            { 5.184, 0.449, -1.665 },
                                            { -3.333, 3.909, 5.544 },
                                            { 3.495, -0.703, 5.494 },
                                            { -4.844, 4.808, -4.166 },
                                            { 3.016, -0.006, -1.534 },
                                            { -3.990, -1.579, 0.691 },
                                            { 6.928, 1.059, 6.088 },
                                            { 9.439, -6.873, 5.741 },
                                            { 6.238, -6.793, -4.822 }
                                        };

    // the box coordinate should be equal to (ll, ur)
    std::pair<MadVoro::Vector3D, MadVoro::Vector3D> boundingBox = diagram.GetBoxCoordinates();

    // Build the diagram serially
    diagram.Build(points);

    std::cout << "ll is " << boundingBox.first << ", ur is " << boundingBox.second << std::endl;
    std::cout << "Number of points (sites/cells) is: " << diagram.GetPointNo() << std::endl;

    size_t examplePointIdx = 3; // index 3

    // the vector `faces` contains the indices of the faces
    const MadVoro::face_vec &faces = diagram.GetCellFaces(examplePointIdx);
    std::cout << "Number of faces of point number " << examplePointIdx << " (" << diagram.GetMeshPoint(examplePointIdx) << ") is " << faces.size() << std::endl;
    std::cout << "The number of neighbors of point number " << examplePointIdx << " is " << faces.size() << std::endl;


    std::vector<MadVoro::Vector3D> facesVertices = diagram.GetFacePoints();

    std::cout << "And these are the faces: " << std::endl;
    for(const size_t &faceIdx : faces)
    {
        MadVoro::Vector3D centerOfMass = diagram.FaceCM(faceIdx);
        const std::pair<size_t, size_t> neighbors = diagram.GetFaceNeighbors(faceIdx); // get neighbors of the face
        
        // one of the values of the pair is our point (`examplePointIdx`), the other is the neighbor
        size_t neighborIdx = (neighbors.first == examplePointIdx)? neighbors.second : neighbors.first;
        std::cout << "Face number " << faceIdx << ", the other side neighbor is point number " << neighborIdx << " (which is " << diagram.GetMeshPoint(neighborIdx) << ")" << std::endl;;
        std::cout << "Vertices of the face are: ";
        for(const size_t &vertexIdx : diagram.GetPointsInFace(faceIdx))
        {
            std::cout << facesVertices[vertexIdx] << ", ";
        }
        std::cout << "\b\b " << std::endl;
    }

    return 0;
}