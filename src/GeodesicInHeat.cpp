#include "GeodesicInHeat.h"
#include <igl/is_border_vertex.h>

GeodesicInHeat::GeodesicInHeat(TriangleMesh& triMesh, boundaryCondition boundaryCon):triMesh(triMesh),boundaryCon(boundaryCon)
{
	/*computeTimeStep(smooth);
	computeLeftSource();
	solveHeat.compute(leftSource);
	solveDistance.compute(triMesh.getCotLaplacian());*/
}

GeodesicInHeat::~GeodesicInHeat()
{
}

void GeodesicInHeat::findBoundaryPoint(){
	boundaryPoint = igl::is_border_vertex(triMesh.getFaces());
}

void GeodesicInHeat::updateSource(size_t sourceVi) {
	size_t VN = triMesh.getVN();
	source.setConstant(VN, 0);
	source.coeffRef(sourceVi) = 1;
}


void GeodesicInHeat::computeTimeStep(double smooth) {
	double avgL = triMesh.getAvgEdgeLength();
	timeStep = avgL * avgL * smooth;
}

void GeodesicInHeat::computeHeatNeuman() {
	solveHeat.compute(leftSource);
	heatNeuman = solveHeat.solve(source);
}
void GeodesicInHeat::computeHeatDirichlet() {//the heat of boundary point = 0
	solveHeat.compute(leftSource);
	heatDirichlet = solveHeat.solve(source);
	size_t cou = 0;
	for (auto it:boundaryPoint) {
		if (it) {
			heatDirichlet(cou) = 0;
		}
		cou++;
	}
}

void GeodesicInHeat::computeHeat() {
	switch (boundaryCon)
	{
	case NEUMAN:
		computeHeatNeuman();
		heat = heatNeuman;
		break;
	case DIRICHLET:
		computeHeatDirichlet();
		heat = heatDirichlet;
		break;
	case NEUMAN_DIRICHLET:
		computeHeatNeuman();
		computeHeatDirichlet();
		heat = 0.5 * (heatNeuman + heatDirichlet);
		break;
	default:
		assert(false);
		break;
	}
}
void GeodesicInHeat::computeDivergence() {
	triMesh.computeDivergence(heat, divergenceHeat);
}
void GeodesicInHeat::computeDistance(size_t& sourceVi, Eigen::VectorXd& geodesicDistance, double smooth) {
	assert(smooth > 0);
	computeTimeStep(smooth);
	computeLeftSource();
	updateSource(sourceVi);
	computeHeat();

	computeDivergence();
	solveDistance.compute(triMesh.getCotLaplacian());
	distance = solveDistance.solve(divergenceHeat);
	geodesicDistance = distance;
}


void GeodesicInHeat::computeLeftSource() {
	leftSource = triMesh.getMass() - timeStep * triMesh.getCotLaplacian();
}
