#include "TriangleMesh.h"
TriangleMesh::TriangleMesh(Eigen::MatrixXd* Vertices, Eigen::MatrixXi* Faces)
{
	Eigen::MatrixXd& V = *Vertices; 
	Eigen::MatrixXi& F = *Faces;

	//size
	VN = V.rows();
	FN = F.rows();

	igl::cotmatrix(V, F, cotLaplacian);
	igl::massmatrix(V, F, igl::MASSMATRIX_TYPE_BARYCENTRIC, mass);
	igl::grad(V, F, grad);
	igl::doublearea(V, F, dblA);
	igl::per_vertex_normals(V, F, normals);
	avgEdgeLength = igl::avg_edge_length(V, F);

	//compute cotLaplacianOperator
	Eigen::SparseMatrix<double> massInv;
	igl::invert_diag(mass, massInv);
	cotLaplacianOperator = 0.5 * massInv * cotLaplacian;

}

TriangleMesh::~TriangleMesh()
{
}

void TriangleMesh::computeDivergence(Eigen::VectorXd& U, Eigen::VectorXd& divergenceU) {
	Eigen::MatrixXd gradU;
	computeGradient(U, gradU);

	//X = -gradU
	gradU = -1.0 * gradU;

	Eigen::MatrixXd gradUrow = gradU.cwiseProduct(dblA.replicate(1, 3));
	gradUrow.resize(3 * FN, 1);
	Eigen::SparseMatrix<double> gradT = grad.transpose();
	divergenceU = 0.5 * gradT * gradUrow;
}

void TriangleMesh::computeGradient(Eigen::VectorXd& U, Eigen::MatrixXd& gradU_norm)
{
	// Compute gradient of U
	Eigen::MatrixXd gradU = Eigen::Map<const Eigen::MatrixXd>((grad * U).eval().data(), FN, 3);
	//Compute gradient magnitude
	//Eigen::VectorXd& gradU_mag
	//const VectorXd gradU_mag = gradU.rowwise().norm();

	//compute unit gradient 
	gradU_norm = gradU.rowwise().normalized();
}

