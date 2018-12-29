#ifndef __LINEAR_SYSTEM_SOLVER_H__
#define __LINEAR_SYSTEM_SOLVER_H__
#include <vector>
#include "Eigen/Dense"
#include "Eigen/Sparse"

// Building Eigen sparse matrix by adding entries one by one
// Eigen Documentation: http://eigen.tuxfamily.org/dox/group__TutorialSparse.html
class SparseMatrixBuilder {
public:
    void AddEntry(int row, int col, double value);
    void AddEntries(const std::vector< Eigen::Triplet< double > >& e);
    std::vector< Eigen::Triplet< double > > OffsetEntries(int row_offset, int col_offset, double coefficient = 1);

    Eigen::SparseMatrix< double > ToSparseMatrix(int num_of_rows, int num_of_cols) const;
    const std::vector< Eigen::Triplet< double > >& GetEntries() const;

private:
    std::vector< Eigen::Triplet< double > > spmatEntries;  //tripletList;
};

// Solving sparse linear systems
class SparseLinearSystemSolver {
public:
    SparseLinearSystemSolver(const Eigen::SparseMatrix< double >& A);
    ~SparseLinearSystemSolver();

    /* Solve by Cholesky factorizations */
    Eigen::VectorXd Solve(const Eigen::VectorXd& b) const;

    /* Solve the sparse linear system Ax = b by using Conjugete Gradient Method. */
    // Note that this is is a static function.
    static void SolveByConjugateGradient(const Eigen::SparseMatrix< double >& A,
                                         const Eigen::VectorXd& b,
                                         int maxIter, double tolerance,
                                         Eigen::VectorXd& x);

private:
    // Eigen Documentation: http://eigen.tuxfamily.org/dox/group__TopicSparseSystems.html
    Eigen::SimplicialLDLT< Eigen::SparseMatrix< double > >* choleskySolver;
};

#endif __LINEAR_SYSTEM_SOLVER_H__
