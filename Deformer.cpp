#include "Deformer.h"
#include <iostream>
#include <limits>

// Pseudo inverse of a matrix in Eigen
// Used in approximated local rotation Laplacian editing
// http://eigen.tuxfamily.org/bz/show_bug.cgi?id=257
// http://eigendobetter.com/
template < typename _Matrix_Type_ >
inline _Matrix_Type_ eigenPinv(const _Matrix_Type_& a,
                               double epsilon = std::numeric_limits< double >::epsilon()) {
    if (a.rows() < a.cols()) {
        Eigen::JacobiSVD< _Matrix_Type_ > svd(a.transpose(),
                                              Eigen::ComputeThinU | Eigen::ComputeThinV);

        double tolerance = epsilon * std::max((double)a.cols(), (double)a.rows()) *
                                     svd.singularValues().array().abs().maxCoeff();

        return (svd.matrixV() *
            (svd.singularValues().array().abs() > tolerance)
            .select(svd.singularValues().array().inverse(), 0).matrix().asDiagonal() *
            svd.matrixU().adjoint()
        ).transpose();
    }

    Eigen::JacobiSVD< _Matrix_Type_ > svd(a,
                                          Eigen::ComputeThinU | Eigen::ComputeThinV);

    double tolerance = epsilon * std::max((double)a.cols(), (double)a.rows()) *
                                 svd.singularValues().array().abs().maxCoeff();

    return svd.matrixV() *
            (svd.singularValues().array().abs() > tolerance)
            .select(svd.singularValues().array().inverse(), 0).matrix().asDiagonal() *
            svd.matrixU().adjoint();
}

// rot: true for estimating local rotations, false for naive Laplacian editing
Deformer::Deformer(Mesh * mesh, bool rot) : solver(NULL), handleWeight(1e3), localRotation(rot) {
    this->mesh = mesh;

    roi_list.clear();
    // keep the location of roi vertices
    for (int i = 0; i < mesh->vList.size(); ++i) {
        if (mesh->vList[i]->Flag())
            roi_list.push_back(mesh->vList[i]);
    }

    // build the diff coord operator with soft constraint
    BuildSystemMatrix();
}

Deformer::~Deformer() {
    if (solver != NULL) {
        delete solver;
		delete Lap;
    }
}

// this is the place where the editing techniques take place
void Deformer::Deform() {
	
	/*please use the SparseLinearSystemSolver to solve  Mx=A.transpose()b*/
 
	VertexList vList = mesh->Vertices();
	Eigen::VectorXd deta_Vertex_X = Eigen::VectorXd::Zero(vList.size() + roi_list.size());
	Eigen::VectorXd deta_Vertex_Y = Eigen::VectorXd::Zero(vList.size() + roi_list.size());
	Eigen::VectorXd deta_Vertex_Z = Eigen::VectorXd::Zero(vList.size() + roi_list.size());

	//add the extra set of entrices on the right hand side
	for (int i = 0; i < vList.size(); i++) {
		deta_Vertex_X[i] = vx[i];
		deta_Vertex_Y[i] = vy[i];
		deta_Vertex_Z[i] = vz[i];
	}
	for (int i = 0; i < roi_list.size(); i++) {
		/*************************/
		/* insert your code here */
		/*************************/

		deta_Vertex_X[i + vList.size()] = roi_list[i]->Position()(0);
		deta_Vertex_Y[i + vList.size()] = roi_list[i]->Position()(1);
		deta_Vertex_Z[i + vList.size()] = roi_list[i]->Position()(2);
	}
	//get the right hand side by left multiplying matrix A. And A is not changed over translation
	Eigen::SparseMatrix<double> A = Lap->ToSparseMatrix(vList.size() + roi_list.size(), vList.size());
	Eigen::VectorXd bx = A.transpose()*deta_Vertex_X;
	Eigen::VectorXd by = A.transpose()*deta_Vertex_Y;
	Eigen::VectorXd bz = A.transpose()*deta_Vertex_Z;




	//get vectors of x, y, z coord after the deformation by solver
	Eigen::VectorXd vx2 = Eigen::VectorXd::Zero(vList.size());
	Eigen::VectorXd vy2 = Eigen::VectorXd::Zero(vList.size());
	Eigen::VectorXd vz2 = Eigen::VectorXd::Zero(vList.size());

	// hint: call the solver here
	vx2 = solver->Solve(bx);
	vy2 = solver->Solve(by);
	vz2 = solver->Solve(bz);
	//set the position of vertices
	for (int i = 0; i < vList.size(); i++){
		Eigen::Vector3d newPosition = Eigen::Vector3d(vx2[i], vy2[i], vz2[i]);
		vList[i]->SetPosition(newPosition);
	}
	
}

// build the differential operator matrix and do factorization
void Deformer::BuildSystemMatrix() {

	VertexList vList = mesh->Vertices();

	//compute the L matrix using cotangent weight
	for (int i = 0; i < vList.size(); i++){
		
		//compute the summation of weight around a vertex;
		int k = vList[i]->Valence();
		HEdge* nextHedge = vList[i]->HalfEdge()->Twin();
		Eigen::Vector3d preVertex, currVertex, nexVertex;
		double cotAphla, cotBeta, sumWeight = 0.0;

		for (int r = 0; r < k; r++){
			currVertex = nextHedge->Start()->Position();
			preVertex = nextHedge->Twin()->Prev()->Start()->Position();
			nexVertex = nextHedge->Next()->Twin()->Start()->Position();
			cotAphla = Mesh::Cot(vList[i]->Position(), nexVertex, currVertex);
			cotBeta = Mesh::Cot(vList[i]->Position(), preVertex, currVertex);
			sumWeight += cotAphla + cotBeta;
			nextHedge = nextHedge->Next()->Twin();
		}

		// compute weight for each one-ring vertex		
		for (int r = 0; r < k; r++){
			cotAphla = Mesh::Cot(vList[i]->Position(), nextHedge->Next()->Twin()->Start()->Position(), nextHedge->Start()->Position());
			cotBeta = Mesh::Cot(vList[i]->Position(), nextHedge->Twin()->Prev()->Start()->Position(), nextHedge->Start()->Position());
			Lap->AddEntry(i, nextHedge->Start()->Index(),  (cotAphla + cotBeta) / sumWeight);	
			nextHedge = nextHedge->Next()->Twin();
		}
		Lap->AddEntry(i, i, -1.0);
	}
	
	//get delta vectors of x, y, z coord for each vertex before transformation
	vx = Eigen::VectorXd::Zero(vList.size());
	vy = Eigen::VectorXd::Zero(vList.size());
	vz = Eigen::VectorXd::Zero(vList.size());

	for (int i = 0; i < vList.size(); i++){
		vx[i] = vList[i]->Position()(0);
		vy[i] = vList[i]->Position()(1);
		vz[i] = vList[i]->Position()(2);
	}

	Eigen::SparseMatrix<double> Laplacian(vList.size(), vList.size());
	Laplacian = Lap->ToSparseMatrix(vList.size(), vList.size());
	
	vx = Laplacian*vx;
	vy = Laplacian*vy;
	vz = Laplacian*vz;
	
	//get the rest part of right hand side matrix and initial the solver
	for (int i = 0; i < roi_list.size(); i++){
		Lap->AddEntry(vList.size() + i, roi_list[i]->Index(), 1.0);
		
	}
	Eigen::SparseMatrix<double> A(vList.size() + roi_list.size(), vList.size());
	Eigen::SparseMatrix<double> M(vList.size(), vList.size());	
	
	A = Lap->ToSparseMatrix(vList.size() + roi_list.size(), vList.size());	
	M = A.transpose()*A;	
	solver = new SparseLinearSystemSolver(M);//δ֪

    // hint: initialize the solver here
    // you may employ the flag variable "localRotation" and the weighting variable "handleWeight".
}
