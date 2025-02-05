#include "hilbert/HilbertOrder3D_Utils.hpp"

int MadVoro::EstimateHilbertIterationNum(vector<Point3D> const& cor)
{
  return static_cast<int>(ceil(log(pow(static_cast<double>(cor.size()), (1.0 / 3.0))) / log(2.0)));
}

// Adjust the vector of points so all coordinate values will be in the range [0,1]
void MadVoro::AdjustPoints(vector<Point3D> const & vPointsIn, vector<Point3D> & vPointsOut)
{
	// The output vector:
	vPointsOut.resize(vPointsIn.size());
	// Vectors holding the X,Y,Z coordinates:
	vector<double> vPointsX, vPointsY, vPointsZ;
	
	Split(vPointsIn, vPointsX, vPointsY, vPointsZ);

	// TOOD - what if all points have the same x (or y, or z) values?

	double dbMinX = *min_element(vPointsX.begin(), vPointsX.end());
	double dbMinY = *min_element(vPointsY.begin(), vPointsY.end());
	double dbMinZ = *min_element(vPointsZ.begin(), vPointsZ.end());

	double dbMaxX = *max_element(vPointsX.begin(), vPointsX.end());
	double dbMaxY = *max_element(vPointsY.begin(), vPointsY.end());
	double dbMaxZ = *max_element(vPointsZ.begin(), vPointsZ.end());
	// The scale factor:
	double dbScaleX = dbMaxX - dbMinX;
	double dbScaleY = dbMaxY - dbMinY;
	double dbScaleZ = dbMaxZ - dbMinZ;
	// To prevent division by zero (very unlikely - double precision!)
	bool bFlagX = MadVoro::Utils::close2zero(dbScaleX);
	bool bFlagY = MadVoro::Utils::close2zero(dbScaleY);
	bool bFlagZ = MadVoro::Utils::close2zero(dbScaleZ);

	// X coordinate:
	if (!bFlagX)
	{
		for (std::size_t ii = 0; ii < vPointsIn.size(); ++ii)
		{
			// Scale the X coordinate:
			vPointsX[ii] -= dbMinX;
			vPointsX[ii] /= dbScaleX;
		}
	}
	else
	{
		fill(vPointsX.begin(), vPointsX.end(), 0);
	}
	// Y coordinate:
	if (!bFlagY)
	{
		for (std::size_t ii = 0; ii < vPointsIn.size(); ++ii)
		{
			// Scale the X coordinate:
			vPointsY[ii] -= dbMinY;
			vPointsY[ii] /= dbScaleY;
		}
	}
	else
	{
		fill(vPointsY.begin(), vPointsY.end(), 0);
	}
	// Z coordinate:
	if (!bFlagZ)
	{
		for (std::size_t ii = 0; ii < vPointsIn.size(); ++ii)
		{
			vPointsZ[ii] -= dbMinZ;
			vPointsZ[ii] /= dbScaleZ;
		}
	}
	else
	{
		fill(vPointsZ.begin(), vPointsZ.end(), 0);
	}

	// Store back the coordiantes in a single output vector:
	for (std::size_t ii = 0; ii < vPointsIn.size(); ++ii)
	{
		vPointsOut[ii].x = vPointsX[ii];
		vPointsOut[ii].y = vPointsY[ii];
		vPointsOut[ii].z = vPointsZ[ii];
	}

	return;
}

void MadVoro::FindEqualIndices(vector<unsigned long long int> const & vD_sorted, vector<vector<std::size_t> > & vOut)
{
	vector<unsigned long long int> vD_sorted_cpy = vD_sorted;
	vector<unsigned long long int> vD_sorted_unq = vD_sorted;
	//vector<unsigned long long int>::iterator it = adjacent_find(vD_sorted.begin(), vD_sorted.end());

	vector<unsigned long long int>::iterator it1, itPrev, itCur;
	it1 = unique(vD_sorted_unq.begin(), vD_sorted_unq.end());

	vD_sorted_unq.resize(static_cast<size_t>(distance(vD_sorted_unq.begin(), it1)));
	
	if (vD_sorted.size() == vD_sorted_unq.size())
	{
		return;
	}

	vOut.reserve(vD_sorted.size() - vD_sorted_unq.size());

	int iCurPrevDist = 0;

	itPrev = vD_sorted_cpy.begin();
	for (it1 = vD_sorted_unq.begin()+1; it1 != vD_sorted_unq.end(); ++it1)
	{
		if (distance( it1 , vD_sorted_unq.end() ) == 0)
		{
			itCur = vD_sorted_cpy.end();
		}
		else
		{
			itCur = find(itPrev, vD_sorted_cpy.end(), *it1);
		}

		iCurPrevDist = static_cast<int>(distance(itPrev, itCur));
		if (1 < iCurPrevDist)
		{
			int iBase = static_cast<int>(distance(vD_sorted_cpy.begin(), itPrev));
			vector<std::size_t> vInd( static_cast<size_t>(iCurPrevDist) );
			// C++11
			// iota(vInd.begin(), vInd.end(), iBase);
			for (int ii = 0; ii < iCurPrevDist; ++ii)
			{
			  vInd[static_cast<size_t>(ii)] = static_cast<size_t>(iBase + ii);
			}
			vOut.push_back(vInd);
		}

		itPrev = itCur;
	}

	iCurPrevDist = static_cast<int>(distance(itPrev, vD_sorted_cpy.end()));
	if (1 < iCurPrevDist )
	{
		int iBase = static_cast<int>(distance(vD_sorted_cpy.begin(), itPrev));

		vector<std::size_t> vInd(static_cast<size_t>(iCurPrevDist));
		for (int ii = 0; ii < iCurPrevDist; ++ii)
		{
		  vInd[static_cast<size_t>(ii)] = static_cast<size_t>(iBase + ii);
		}
		vOut.push_back(vInd);
	}

	return;
}
