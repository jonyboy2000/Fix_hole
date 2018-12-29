#ifndef __DEFORMER_H__
#define __DEFORMER_H__
#include "Mesh.h"
#include "LinearSystemSolver.h"

////////// Laplacian Deformer //////////
class Deformer
{
public:
    // constructor
    Deformer(Mesh * mesh, bool rot = false);
    ~Deformer();

    // editing
    void Deform();

private:
    // Build left hand side matrix and pre-factorize it
    void BuildSystemMatrix();
    
    SparseLinearSystemSolver* solver;
	Mesh * mesh;
    VertexList roi_list;

	//store the original delta_X (detail)
	Eigen::VectorXd vx, vy, vz;
	SparseMatrixBuilder* Lap = new SparseMatrixBuilder();

    double handleWeight; // The handle weight in the linear system to be solved
    bool localRotation; // true for estimating local rotations, false for naive Laplacian editing
};


#endif // __DEFORMER_H__