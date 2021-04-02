#ifndef GEODESICINHEAT_H
#define GEODESICINHEAT_H


#pragma once
#include "TriangleMesh.h"
#include <Eigen/SparseCholesky>

class GeodesicInHeat
{
public:
	enum boundaryCondition { NEUMAN, DIRICHLET, NEUMAN_DIRICHLET };

	GeodesicInHeat(TriangleMesh& triMesh, boundaryCondition boundaryCon = NEUMAN_DIRICHLET);
	~GeodesicInHeat();
	void computeTimeStep(double smooth);
	void computeDistance(size_t& sourceVi, Eigen::VectorXd& geodesicDistance, double smooth = 1.0);
	void computeLeftSource();
	double getTimeStep()
	{
		return timeStep;
	}
	

private:
	TriangleMesh triMesh;
	boundaryCondition boundaryCon;
	double timeStep;
	Eigen::SparseMatrix<double> leftSource;//leftSource = M - t*Lc
	std::vector<bool> boundaryPoint;
	Eigen::VectorXd source;
	Eigen::VectorXd heatNeuman;
	Eigen::VectorXd heatDirichlet;
	Eigen::VectorXd heat;
	Eigen::VectorXd divergenceHeat;
	Eigen::VectorXd distance;
	Eigen::SimplicialCholesky<Eigen::SparseMatrix<double>> solveHeat;
	Eigen::SimplicialCholesky<Eigen::SparseMatrix<double>> solveDistance;

	void findBoundaryPoint();
	void updateSource(size_t sourceVi);
	void computeHeatNeuman();
	void computeHeatDirichlet();
	void computeHeat();
	void computeDivergence();
};


#endif // !GEODESICINHEAT_H
