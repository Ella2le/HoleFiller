#include "Mesh.h"
#include "HoleFiller.h"

int main(int, char ** argv)
{
	char* filename = "input/bunny_with_holes.off";
	char* result = "output/bunny_without_holes.off";

	Mesh * m = new Mesh();
	m->loadOff(filename);
	HoleFiller hFiller(*m);
	hFiller.fillHoles();
	m->exportOff(result);
	return 0;
}