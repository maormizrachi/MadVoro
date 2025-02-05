#include "hilbert/HilbertOrder3D.hpp"

using namespace MadVoro;

// Constructor - uses the reference Hilbert curve shape:
MadVoro::HilbertCurve3D_shape::HilbertCurve3D_shape() : m_vShapePoints(boost::array<Point3D, 7>())
{
	m_vShapePoints[0] = Point3D(0, 0, -1);
	m_vShapePoints[1] = Point3D(0, 1, 0);
	m_vShapePoints[2] = Point3D(0, 0, 1);
	m_vShapePoints[3] = Point3D(-1, 0, 0);
	m_vShapePoints[4] = Point3D(0, 0, -1);
	m_vShapePoints[5] = Point3D(0, -1, 0);
	m_vShapePoints[6] = Point3D(0, 0, 1);
}

// Compare to a given Hilbert curve shape by comparing pairs of shape points:
bool MadVoro::HilbertCurve3D_shape::operator==(const HilbertCurve3D_shape & shape) const
{
	bool b = true;
	for (std::size_t ii = 0; ii < m_vShapePoints.size(); ++ii)
	{
		b = b && (m_vShapePoints[ii] == shape.m_vShapePoints[ii]);
	}

	return b;
}

// Constructor - performs all required initiallizations and preprocessing:
MadVoro::HilbertCurve3D::HilbertCurve3D() :m_vRotatedShapes(boost::array<HilbertCurve3D_shape, NUMBER_OF_SHAPES> ()), m_vRotations(boost::array < vector<int>, NUMBER_OF_SHAPES >()),
m_vShapeRecursion(boost::array< boost::array<int, 8> , NUMBER_OF_SHAPES>())
{
	int rot[MAX_ROTATION_LENGTH];
	//	int iRotLength;
	for (size_t iRotIndex = 1; iRotIndex < NUMBER_OF_SHAPES; ++iRotIndex)
	{
	  const int iRotLength = GetRotation(rot, static_cast<int>(iRotIndex));
		m_vRotations[iRotIndex].assign(rot, rot + iRotLength);
	}

	for (size_t ii = 1; ii < NUMBER_OF_SHAPES; ++ii)
	{
	  RotateShape(static_cast<int>(ii), m_vRotations[ii]);
	}

	BuildRecursionRule();
	BuildShapeOrder();
}

// FindShapeIndex - returns the index of a shape:
int MadVoro::HilbertCurve3D::FindShapeIndex(const HilbertCurve3D_shape & roShape)
{
  for (size_t ii = 0; ii < NUMBER_OF_SHAPES; ++ii)
	{
		if (roShape == m_vRotatedShapes[ii])
		{
		  return static_cast<int>(ii);
		}
	}
	// TODO - manage this kind of return value (error)
	return -1;
}

void MadVoro::HilbertCurve3D::BuildRecursionRule()
{
	// Reference recursion rule:
	m_vShapeRecursion[0][0] = 12;
	m_vShapeRecursion[0][1] = 16;
	m_vShapeRecursion[0][2] = 16;
	m_vShapeRecursion[0][3] = 2;
	m_vShapeRecursion[0][4] = 2;
	m_vShapeRecursion[0][5] = 14;
	m_vShapeRecursion[0][6] = 14;
	m_vShapeRecursion[0][7] = 10;

	HilbertCurve3D_shape oTempShape;
	// What about ii=0? not necessary 
	for (size_t ii = 0; ii < NUMBER_OF_SHAPES; ++ii)
	{
		for (size_t jj = 0; jj < 8; ++jj)
		{
			// Rotate the appropriate block of the reference recursion rule, according to the ii rotation scheme:
		  RotateShape(m_vRotatedShapes[static_cast<size_t>(m_vShapeRecursion[0][jj])], oTempShape, static_cast<int>(ii));
			// Find the shape index of the rotated shape:
			m_vShapeRecursion[ii][jj] = FindShapeIndex(oTempShape);
		}
	}

	return;
}

// Return the rotation scheme of the ii rotation (manually precalculated):
int MadVoro::HilbertCurve3D::GetRotation(int * piRotation, int iRotationIndex)
{
	switch (iRotationIndex)
	{
	case 1:
		piRotation[0] = 1;
		return 1;
	case 2:
		piRotation[0] = 1;
		piRotation[1] = 1;
		return 2;
	case 3:
		piRotation[0] = 1;
		piRotation[1] = 1;
		piRotation[2] = 1;
		return 3;
	
	case 4:
		piRotation[0] = 2;
		return 1;
	case 5:
		piRotation[0] = 2;
		piRotation[1] = 2;
		return 2;
	case 6:
		piRotation[0] = 2;
		piRotation[1] = 2;
		piRotation[2] = 2;
		return 3;

	case 7:
		piRotation[0] = 3;
		return 1;
	case 8:
		piRotation[0] = 3;
		piRotation[1] = 3;
		return 2;
	case 9:
		piRotation[0] = 3;
		piRotation[1] = 3;
		piRotation[2] = 3;
		return 3;

	case 10:
		piRotation[0] = 1;
		piRotation[1] = 2;
		return 2;
	case 11:
		piRotation[0] = 2;
		piRotation[1] = 1;
		return 2;
	case 12:
		piRotation[0] = 3;
		piRotation[1] = 1;
		return 2;
	case 13:
		piRotation[0] = 2;
		piRotation[1] = 3;
		return 2;

	case 14:
		piRotation[0] = 1;
		piRotation[1] = 1;
		piRotation[2] = 1;
		piRotation[3] = 3;
		return 4;
	case 15:
		piRotation[0] = 2;
		piRotation[1] = 2;
		piRotation[2] = 2;
		piRotation[3] = 1;
		return 4;
	case 16:
		piRotation[0] = 3;
		piRotation[1] = 3;
		piRotation[2] = 3;
		piRotation[3] = 2;
		return 4;
	case 17:
		piRotation[0] = 3;
		piRotation[1] = 2;
		piRotation[2] = 3;
		piRotation[3] = 2;
		return 4;
	case 18:
		piRotation[0] = 3;
		piRotation[1] = 2;
		piRotation[2] = 2;
		return 3;

	case 19:
		piRotation[0] = 1;
		piRotation[1] = 1;
		piRotation[2] = 1;
		piRotation[3] = 2;
		piRotation[4] = 2;
		return 5;
	case 20:
		piRotation[0] = 3;
		piRotation[1] = 3;
		piRotation[2] = 2;
		return 3;
	case 21:
		piRotation[0] = 3;
		piRotation[1] = 3;
		piRotation[2] = 1;
		return 3;
	case 22:
		piRotation[0] = 2;
		piRotation[1] = 2;
		piRotation[2] = 3;
		return 3;
	case 23:
		piRotation[0] = 1;
		piRotation[1] = 1;
		piRotation[2] = 2;
		return 3;

	default:
		return 0;
		//		break;
	}
}

// Rotate a shape:
void MadVoro::HilbertCurve3D::RotateShape(int iShapeIndex, vector<int> vAxes)
{
  int iSign;

	for (std::size_t ii = 0; ii < 7; ++ii)
	{
		for (std::size_t iAx = 0; iAx < vAxes.size(); ++iAx)
		{
			// A trick to find the sign of vAxes[iAx]:
			iSign = (vAxes[iAx] > 0) - (vAxes[iAx] < 0);

			switch(std::abs(vAxes[iAx]))
			{
			case 1:
			  m_vRotatedShapes[static_cast<size_t>(iShapeIndex)].m_vShapePoints[ii].RotateX( iSign * PI / 2 );
				break;
			case 2:
			  m_vRotatedShapes[static_cast<size_t>(iShapeIndex)].m_vShapePoints[ii].RotateY( iSign * PI / 2 );
				break;
			case 3:
			  m_vRotatedShapes[static_cast<size_t>(iShapeIndex)].m_vShapePoints[ii].RotateZ( iSign * PI / 2 );
				break;
			default:
				break;
			}
		}
		// Round off the results:
		m_vRotatedShapes[static_cast<size_t>(iShapeIndex)].m_vShapePoints[ii].Round();
	}
}

void MadVoro::HilbertCurve3D::RotateShape(HilbertCurve3D_shape const & roShape, HilbertCurve3D_shape & roShapeOut , int iRotationIndex)
{
  int iSign;

	std::vector<int> vAxes = m_vRotations[static_cast<size_t>(iRotationIndex)];
	roShapeOut = roShape;

	for (int ii = 0; ii < 7; ++ii)
	{
		for (std::size_t iAx = 0; iAx < vAxes.size(); ++iAx)
		{
			// A trick to find the sign of vAxes[iAx]:
			iSign = (vAxes[iAx] > 0) - (vAxes[iAx] < 0);

			switch(std::abs(vAxes[iAx]))
			{
			case 1:
			  roShapeOut.m_vShapePoints[static_cast<size_t>(ii)].RotateX(iSign * PI / 2);
				break;
			case 2:
			  roShapeOut.m_vShapePoints[static_cast<size_t>(ii)].RotateY(iSign * PI / 2);
				break;
			case 3:
			  roShapeOut.m_vShapePoints[static_cast<size_t>(ii)].RotateZ(iSign * PI / 2);
				break;
			default:
				break;
			}
		}
		// Round off the results:
		roShapeOut.m_vShapePoints[static_cast<size_t>(ii)].Round();
	}

	return;
}

void MadVoro::HilbertCurve3D::BuildShapeOrder()
{
	boost::array<int, 8> vShapeVerticesX;
	boost::array<int, 8> vShapeVerticesY;
	boost::array<int, 8> vShapeVerticesZ;

	vShapeVerticesX[0] = 0;
	vShapeVerticesY[0] = 0;
	vShapeVerticesZ[0] = 0;

	for (std::size_t iShapeInd = 0; iShapeInd < NUMBER_OF_SHAPES; ++iShapeInd)
	{
		for (std::size_t ii = 0; ii < m_vRotatedShapes[iShapeInd].m_vShapePoints.size(); ++ii)
		{
		  vShapeVerticesX[ii + 1] = static_cast<int>(vShapeVerticesX[ii] + m_vRotatedShapes[iShapeInd].m_vShapePoints[ii].x);
			vShapeVerticesY[ii + 1] = static_cast<int>(vShapeVerticesY[ii] + m_vRotatedShapes[iShapeInd].m_vShapePoints[ii].y);
			vShapeVerticesZ[ii + 1] = static_cast<int>(vShapeVerticesZ[ii] + m_vRotatedShapes[iShapeInd].m_vShapePoints[ii].z);
		}

		int iMinX = *std::min_element(vShapeVerticesX.begin(), vShapeVerticesX.end());
		int iMinY = *std::min_element(vShapeVerticesY.begin(), vShapeVerticesY.end());
		int iMinZ = *std::min_element(vShapeVerticesZ.begin(), vShapeVerticesZ.end());

		for (std::size_t jj = 0; jj < vShapeVerticesX.size(); ++jj)
		{
			vShapeVerticesX[jj] -= iMinX;
			vShapeVerticesY[jj] -= iMinY;
			vShapeVerticesZ[jj] -= iMinZ;
		}

		for (std::size_t kk = 0; kk < vShapeVerticesX.size(); ++kk)
		{
			m_mShapeOrder[iShapeInd][vShapeVerticesX[kk]][vShapeVerticesY[kk]][vShapeVerticesZ[kk]] = static_cast<int>(kk);
		}
	}

	return;
}

unsigned long long int MadVoro::HilbertCurve3D::Hilbert3D_xyz2d(Point3D const & rvPoint, int numOfIterations) const
{
	// Extract the coordinates:
	double x = rvPoint.x;
	double y = rvPoint.y;
	double z = rvPoint.z;

	// The output distance along the 3D-Hilbert Curve:
	unsigned long long int d = 0;

	// The current shape index:
	int iCurrentShape = 0;
	// The octant number:
	//	int iOctantNum = 0;
	// A temp variable - storing the current (negative) power of 2
	//	double dbPow2;
	// Variables indicating the current octant:
	//bool bX, bY, bZ;
	for (int iN = 1; iN <= numOfIterations; ++iN)
	{
		// Calculate the current power of 0.5:
	  const double dbPow2 = 1.0 / (1 << iN);
	  const bool bX = x > dbPow2;
	  const bool bY = y > dbPow2;
	  const bool bZ = z > dbPow2;

		x -= dbPow2*bX;
		y -= dbPow2*bY;
		z -= dbPow2*bZ;

		// Multiply the distance by 8 (for every recursion iteration):
		d = d << 3;
		const int iOctantNum = m_mShapeOrder[iCurrentShape][bX][bY][bZ];
		d = d + static_cast<size_t>(iOctantNum);
		iCurrentShape = m_vShapeRecursion[static_cast<size_t>(iCurrentShape)][static_cast<size_t>(iOctantNum)];
	}

	return d;
}

std::vector<std::size_t> MadVoro::GetGlobalHibertIndeces(std::vector<Point3D> const& cor, Point3D const& ll, Point3D const& ur,size_t &Hmax)
{
	std::vector<std::size_t> res;
	std::size_t Niter = 20;
	Hmax = static_cast<size_t>(std::pow(static_cast<size_t>(2), static_cast<size_t>(3*Niter)));
	HilbertCurve3D oHilbert;
	Point3D dx = ur - ll,vtemp;
	std::size_t Ncor = cor.size();
	res.resize(Ncor);
	for (size_t i = 0; i < Ncor; ++i)
	{
		vtemp.x = (cor[i].x - ll.x) / dx.x;
		vtemp.y = (cor[i].y - ll.y) / dx.y;
		vtemp.z = (cor[i].z - ll.z) / dx.z;
		res[i]=static_cast<size_t>(oHilbert.Hilbert3D_xyz2d(vtemp,static_cast<int>(Niter)));
	}
	return res;
}

std::vector<std::size_t> MadVoro::HilbertOrder3D(std::vector<Point3D> const& cor)
{
	static int recursiveCalls = 0;
	recursiveCalls++;

	if(recursiveCalls >= MAX_HILBERT_RECURSIVE_CALLS)
	{
		MadVoro::Exception::MadVoroException eo("HilbertOrder3D: too many recursive calls");
		for(size_t i = 0; i < cor.size(); ++i)
		{
			for(size_t j = 0; j < cor.size(); ++j)
			{
				if(i != j and cor[i] == cor[j])
				{
					eo.Append2ErrorMessage(" - Duplicated point found");
					eo.addEntry("Point1", cor[i]);
					eo.addEntry("Point2", cor[j]);
					throw eo;
				}
			}
		}
		eo.Append2ErrorMessage(" - Though no duplicated points found");
		throw eo;
	}
	// If only 1 or 2 points are provided - do not reorder them
	if ( 2 >= cor.size() )
	{
		std::vector<std::size_t> vIndSort( cor.size() );
		for (std::size_t ii = 0; ii < cor.size(); ++ii)
		{
			vIndSort[ii] = ii;
		}
		recursiveCalls--;
		return vIndSort;
	}
	// Create a 3D-Hilbert Curve Object:
	HilbertCurve3D oHilbert;

	// Allocate an output vector:
	size_t N = cor.size();
	std::vector<unsigned long long int> vOut;
	vOut.resize(N);
	
	// Estimate the number of required iterations:
	int numOfIterations = EstimateHilbertIterationNum(cor);

	std::vector<Point3D> vAdjustedPoints;

	// Adjust the points coordinates to the unit cube:
	AdjustPoints(cor, vAdjustedPoints);

	// Run throught the points, and calculate the Hilbert distance of each:
	for (size_t ii = 0; ii < N; ++ii)
	{
		vOut[ii]=oHilbert.Hilbert3D_xyz2d(vAdjustedPoints[ii], numOfIterations+8);
		//vOut.push_back(oHilbert.Hilbert3D_xyz2d(vAdjustedPoints[ii], 2));
	}
	// Get the sorting indices:
	std::vector<std::size_t> vIndSort;
	ordered(vOut, vIndSort);
	// Reorder the Hilbert distances vector (according to the sorting indices):
	reorder( vOut, vIndSort );

	// Find indices with repeated Hilbert distance:
	std::vector<std::vector<std::size_t>> vEqualIndices;
	FindEqualIndices(vOut, vEqualIndices);
	
	// If all points have different Hilbert distances, return the sorting indices:
	if (vEqualIndices.empty())
	{
		recursiveCalls--;
		return vIndSort;
	}
	else
	{
		for (const std::vector<size_t> &equalList : vEqualIndices)
		{
			std::vector<Point3D> vPointsInner(equalList.size() );
			std::vector<std::size_t> vIndInner(equalList.size());
			std::vector<std::size_t> vIndSortInner(equalList.size());
			std::vector<std::size_t> vIndSortInner_cpy(equalList.size());

			// Store the points with the equal indices
			for (std::size_t jj = 0; jj < equalList.size() ; ++jj)
			{
				vIndInner[jj] = vIndSort[equalList[jj]];
				vPointsInner[jj] = cor[vIndInner[jj]];
				vIndSortInner_cpy[jj] = vIndSort[vIndInner[jj]];
			}
			
			// Sort the repeated points:
			vIndSortInner = HilbertOrder3D(vPointsInner);
	//		vector<std::size_t> vIndSortTemp = vIndSort;
			for (std::size_t kk = 0; kk < vIndSortInner.size(); ++kk)
			{
				vIndSort[vIndInner[kk]] = vIndSortInner_cpy[vIndSortInner[kk]];
			}
		}

		// Return the sorting indices:
		recursiveCalls--;
		return vIndSort;
	}
}
