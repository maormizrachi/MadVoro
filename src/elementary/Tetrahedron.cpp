#include "Tetrahedron.hpp"

using namespace MadVoro;

MadVoro::Tetrahedron::Tetrahedron(): points(), neighbors(), newTetra(true)
{}

MadVoro::Tetrahedron::Tetrahedron(Tetrahedron const & other)
{
	this->newTetra = other.newTetra;
#ifdef __INTEL_COMPILER
#pragma omp simd
#endif
	for (int i = 0; i < 4; i++)
	{
		points[i] = other.points[i];
		neighbors[i] = other.neighbors[i];
	}
}

MadVoro::Tetrahedron::~Tetrahedron()
{}

Tetrahedron &MadVoro::Tetrahedron::operator=(Tetrahedron const & other)
{
	this->newTetra = other.newTetra;
	if (&other == this)
		return *this;
#ifdef __INTEL_COMPILER
#pragma omp simd
#endif
	for (int i = 0; i < 4; ++i)
	{
		points[i] = other.points[i];
		neighbors[i] = other.neighbors[i];
	}
	return *this;
}
