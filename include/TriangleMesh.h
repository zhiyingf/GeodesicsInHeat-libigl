//#pragma once
#ifndef TRIANGLEMESH_H
#define TRIANGLEMESH_H

#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <igl/readOBJ.h>
#include <igl/readOFF.h>
#include <igl/readPLY.h>
#include <igl/readSTL.h>
#include <igl/cotmatrix.h>
#include <igl/doublearea.h>
#include <igl/grad.h>
#include <igl/massmatrix.h>
#include <igl/per_vertex_normals.h>
#include <igl/avg_edge_length.h>
#include <igl/invert_diag.h>

class TriangleMesh
{
public:
	TriangleMesh(Eigen::MatrixXd* Vertices, Eigen::MatrixXi* Faces);
	~TriangleMesh();
	Eigen::MatrixXd getVertices() {
		return *Vertices;
	}
	Eigen::MatrixXi getFaces() {
		return *Faces;
	}
	Eigen::SparseMatrix<double> getCotLaplacian() {
		return cotLaplacian;
	}
	Eigen::SparseMatrix<double> getCotLaplacianOperator(){
		return cotLaplacianOperator;
	}
	Eigen::SparseMatrix<double> getMass() {
		return mass;
	}
	Eigen::SparseMatrix<double> getGrad() {
		return grad;
	}
	double getAvgEdgeLength() {
		return avgEdgeLength;                                     
	}
	size_t getVN() {
		return VN;
	}
	size_t getFN() {
		return FN;
	}


private:
	Eigen::MatrixXd* Vertices = nullptr;
	Eigen::MatrixXi* Faces = nullptr;
	Eigen::SparseMatrix<double> cotLaplacian;
	Eigen::SparseMatrix<double> mass;
	Eigen::SparseMatrix<double> cotLaplacianOperator;
	Eigen::SparseMatrix<double> grad;
	Eigen::VectorXd dblA;
	Eigen::MatrixXd normals;
	//Eigen::MatrixXd tangents;
	double avgEdgeLength;
	size_t VN;
	size_t FN;

public:
	void computeDivergence(Eigen::VectorXd& U, Eigen::VectorXd& divergenceU);
	void computeGradient(Eigen::VectorXd& U, Eigen::MatrixXd& gradU_norm);
	
};

#endif // !TRIANGLEMESH_H
