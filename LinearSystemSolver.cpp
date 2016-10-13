#include <iostream>
#include "LinearSystemSolver.h"

void SparseMatrixBuilder::AddEntry(int row, int col, double value) {
    spmatEntries.push_back(Eigen::Triplet< double >(row, col, value));
}

void SparseMatrixBuilder::AddEntries(const std::vector< Eigen::Triplet< double > >& e) {
    spmatEntries.insert(spmatEntries.end(), e.begin(), e.end());
}

std::vector< Eigen::Triplet< double > > SparseMatrixBuilder::OffsetEntries(int row_offset, int col_offset, double coefficient) {
    if (!row_offset && !col_offset && coefficient == 1) {
        return spmatEntries;
    }

    std::vector< Eigen::Triplet< double > > newSpmatEntries;
    for (int i = 0; i < spmatEntries.size(); ++i) {
        Eigen::Triplet< double >& e = spmatEntries[i];
        newSpmatEntries.push_back(Eigen::Triplet< double >(e.row() + row_offset,
                                                           e.col() + col_offset,
                                                           coefficient * e.value()));
    }

    return newSpmatEntries;
}

Eigen::SparseMatrix< double > SparseMatrixBuilder::ToSparseMatrix(int num_of_rows, int num_of_cols) const {
    Eigen::SparseMatrix< double > spmat(num_of_rows, num_of_cols);
    spmat.setFromTriplets(spmatEntries.begin(), spmatEntries.end());
    return spmat;
}

const std::vector< Eigen::Triplet< double > >& SparseMatrixBuilder::GetEntries() const {
    return spmatEntries;
}

SparseLinearSystemSolver::SparseLinearSystemSolver(const Eigen::SparseMatrix< double >& A) :
    choleskySolver(NULL) {

    choleskySolver = new Eigen::SimplicialLDLT< Eigen::SparseMatrix< double > >();
    // Cholesky factorizations
    choleskySolver->compute(A);
    if (choleskySolver->info() != Eigen::Success) {
        // decomposition failed
        std::cout << "sparse decomposition failed" << std::endl;
    } else {
        std::cout << "sparse decomposition succeeded" << std::endl;
    }
}

SparseLinearSystemSolver::~SparseLinearSystemSolver() {
    if (choleskySolver != NULL) {
        delete choleskySolver;
    }
}

Eigen::VectorXd SparseLinearSystemSolver::Solve(const Eigen::VectorXd& b) const {
    Eigen::VectorXd x = choleskySolver->solve(b);

    if (choleskySolver->info() != Eigen::Success) {
        // decomposition failed
        std::cout << "solving the linear system failed" << std::endl;
    } else {
        std::cout << "solving the linear system succeeded" << std::endl;
    }

    return x;
}

void SparseLinearSystemSolver::SolveByConjugateGradient(const Eigen::SparseMatrix< double >& A,
                                                        const Eigen::VectorXd& b,
                                                        int maxIter, double tolerance,
                                                        Eigen::VectorXd& x) {
    /*************************/
    /* insert your code here */
    /*************************/

    /*====== Programming Assignment 1 ======*/
    Eigen::VectorXd r = b - A*x;
    Eigen::VectorXd d = r;
    Eigen::VectorXd r_1;
    double alpha;
    double beta;
    for (int i = 0; i < maxIter; i++)
    {
        if (r.norm() > tolerance)
        {
            alpha = r.transpose().dot(r)/(d.transpose()*A*d);
            x = x + alpha * d;
            r_1 = r - alpha * A * d;
            beta = r_1.reverse().dot(r_1)/r.reverse().dot(r);
            r = r_1;
            d = r + beta*d;
        }
    }
}
