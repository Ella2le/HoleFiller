#include "Mesh.h"
#include "HoleFiller.h"

int main(int, char ** argv)
{
	const char* filename = "../in_out/bunny_with_holes.off";
	const char* result = "../in_out/bunny_without_holes.off";

	Mesh * m = new Mesh();
	m->loadOff(filename);
	HoleFiller hFiller(*m);
	hFiller.fillHoles();
	m->exportOff(result);
	return 0;
}