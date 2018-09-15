#include "mesh.h"

int Classify(const xmesh& m, double4 p) {
	// Accurate: shoot a ray and count crossings
	// if ray hits any edge or vertex, shoot a new ray and repeat
	THROW(not_implemented);
}

int Classify(const xmesh& m, const segment3& s) {
	// ??
	THROW(not_implemented);
}

// This will be a ground truth function. Slow and accurate.
int Classify(const xmesh3& ma, const xmesh3& mb) {
	if (Classify(ma, AnyVertex(mb)) == -1 || Classify(mb, AnyVertex(ma)) == -1)
		return -1;

	// TODO check if any edge of one is penetrating the other

	THROW(not_implemented);
}
