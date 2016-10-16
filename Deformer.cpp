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
    }
}




// this is the place where the editing techniques take place
void Deformer::Deform() {
    /*************************/
    /* insert your code here */
    /*************************/
    
    /*====== Programming Assignment 2 ======*/
    
    VertexList vList = mesh->vList;
    int N = (int)vList.size();
    int K = (int)roi_list.size();
    
    // update constraints' positions
    for (int i = 0; i < K; ++i)
    {
        Eigen::Vector3d p = roi_list[i]->Position();
        b_.coeffRef(i + 3 * N, 0) = p(0)*handleWeight;
        b_.coeffRef(i + 3 * N + K, 0) = p(1)*handleWeight;
        b_.coeffRef(i + 3 * N + 2*K, 0) = p(2)*handleWeight;
    }
    
    Eigen::SparseMatrix<double> AT_b = A_Trans*b_;
    Eigen::VectorXd A_b(3*N);
    for (int i = 0; i < 3*N; ++i)
    {
        A_b(i) = AT_b.coeff(i, 0);
    }
    Eigen::VectorXd V_ = solver->Solve(A_b);
    
    for (int i = 0; i < N; ++i)
    {
        Eigen::Vector3d p;
        p(0) = V_(i);
        p(1) = V_(i + N);
        p(2) = V_(i + 2*N);
        vList[i]->SetPosition(p);
    }
    
}


void Deformer::BuildSystemMatrix() {
    /*************************/
    /* insert your code here */
    /*************************/
    
    /*====== Programming Assignment 2 ======*/
    
    
    // hint: initialize the solver here
    // you may employ the flag variable "localRotation" and the weighting variable "handleWeight".
    
    // First implement the naive method
    // Solve the x,y,z at one time using (number of vertices + number of selected vertices) * (3 * number of vertices) matrix.
    
    
    VertexList vList = mesh->vList;
    int numberVertice = vList.size();
    int numberConstrian = roi_list.size();
    
    SparseMatrixBuilder A;   //This is to build the Matrix A
    SparseMatrixBuilder b;   //This is to build the vector b
    Eigen::SparseMatrix<double> A_Original;
    Eigen::SparseMatrix<double> A_Modified;
    
    
    SparseMatrixBuilder L; //This is the matrix caculate the meancurvature
    
    //Here we constructe L for calculating mean curvature
    //And we consturcte A via adding some 1 under L...
    for (int i = 0; i < numberVertice; ++i)
    {
        std::vector<Eigen::Triplet<double> > rowWeights(0);
        double weightSum = 0;
        Vertex *pV = vList[i];
        if (pV->Valence() < 3)
            continue;
        OneRingHEdge ring(pV);
        HEdge *curr = NULL;
        while ( (curr = ring.NextHEdge()) )
        {
            int index = curr->End()->Index();
            Eigen::Vector3d v1 = curr->Next()->End()->Position();
            Eigen::Vector3d v2 = curr->Twin()->Next()->End()->Position();
            double cot1 = mesh->Cot(pV->Position(), v1, curr->End()->Position());
            double cot2 = mesh->Cot(pV->Position(), v2, curr->End()->Position());
            rowWeights.push_back(Eigen::Triplet<double>(i, index, cot1 + cot2));
            weightSum += (cot1 + cot2);
        }
        for (int j = 0; j < rowWeights.size(); ++j)
        {
            L.AddEntry(i, rowWeights[j].col(), -rowWeights[j].value() / weightSum);
            A.AddEntry(i, rowWeights[j].col(), -rowWeights[j].value() / weightSum);
            A.AddEntry(i + numberVertice, rowWeights[j].col() + numberVertice, -rowWeights[j].value() / weightSum);
            A.AddEntry(i + 2*numberVertice, rowWeights[j].col() + 2*numberVertice, -rowWeights[j].value() / weightSum);
        }
        L.AddEntry(i, i, 1);
        A.AddEntry(i, i, 1);
        A.AddEntry(i + numberVertice, i + numberVertice, 1);
        A.AddEntry(i + 2*numberVertice, i + 2*numberVertice, 1);
    }
    Eigen::SparseMatrix<double> Laplacin = L.ToSparseMatrix(numberVertice, numberVertice);
    
    SparseMatrixBuilder V;
    Eigen::SparseMatrix<double> V_; //A vector store meancurvature and handles
    Eigen::SparseMatrix<double> MeanCurve; //MeanCurvature.
    for (int i = 0; i < numberVertice; ++i)
    {
        Eigen::Vector3d p = vList[i]->Position();
        V.AddEntry(i, 0, p(0));
        V.AddEntry(i, 1, p(1));
        V.AddEntry(i, 2, p(2));
    }
    V_ = V.ToSparseMatrix(numberVertice, 3);
    MeanCurve = Laplacin*V_;
    
    
    for(int i=0; i<numberVertice; ++i)
    {
        double x = MeanCurve.coeff(i,0);
        double y = MeanCurve.coeff(i,0);
        double z = MeanCurve.coeff(i,0);
        /*
        //Here we consturct Ai for each vertex
        //And we calculate s, vector h, vector t.
        if(true)
        {
            Vertex *pV = vList[i];
            int K = pV->Valence();
            Eigen::MatrixXd Ati(3*(K+1), 7);
            Eigen::VectorXd bti(3*(K+1));
            Eigen::VectorXd sht(7);
            double vx = pV->Position()[0];
            double vy = pV->Position()[1];
            double vz = pV->Position()[2];
            Ati(0, 0) = vx;
            Ati(0, 1) = 0;
            Ati(0, 2) = vz;
            Ati(0, 3) = -vy;
            Ati(0, 4) = 1;
            Ati(0, 5) = 0;
            Ati(0, 6) = 0;
            bti(0) = vx;
            
            Ati(1, 0) = vy;
            Ati(1, 1) = -vz;
            Ati(1, 2) = 0;
            Ati(1, 3) = vx;
            Ati(1, 4) = 0;
            Ati(1, 5) = 1;
            Ati(1, 6) = 0;
            bti(1) = vy;
            
            Ati(2, 0) = vz;
            Ati(2, 1) = vy;
            Ati(2, 2) = -vx;
            Ati(2, 3) = 0;
            Ati(2, 4) = 0;
            Ati(2, 5) = 0;
            Ati(2, 6) = 1;
            bti(2) = vz;
            
            OneRingVertex ring(pV);
            Vertex *pCurr = NULL;
            int counter = 0;
            while( (pCurr=ring.NextVertex()) )
            {
                counter++;
                vx = pCurr->Position()[0];
                vy = pCurr->Position()[1];
                vz = pCurr->Position()[2];
                
                Ati(counter*3, 0) = vx;
                Ati(counter*3, 1) = 0;
                Ati(counter*3, 2) = vz;
                Ati(counter*3, 3) = -vy;
                Ati(counter*3, 4) = 1;
                Ati(counter*3, 5) = 0;
                Ati(counter*3, 6) = 0;
                bti(counter*3) = vx;
                
                Ati(counter*3+1, 0) = vy;
                Ati(counter*3+1, 1) = -vz;
                Ati(counter*3+1, 2) = 0;
                Ati(counter*3+1, 3) = vx;
                Ati(counter*3+1, 4) = 0;
                Ati(counter*3+1, 5) = 1;
                Ati(counter*3+1, 6) = 0;
                bti(counter*3+1) = vy;
                
                Ati(counter*3+2, 0) = vz;
                Ati(counter*3+2, 1) = vy;
                Ati(counter*3+2, 2) = -vx;
                Ati(counter*3+2, 3) = 0;
                Ati(counter*3+2, 4) = 0;
                Ati(counter*3+2, 5) = 0;
                Ati(counter*3+2, 0) = 1;
                bti(counter*3+2) = vz;
            }
            Eigen::MatrixXd AT = Ati.transpose();
            Eigen::MatrixXd AA = AT*Ati;
            sht = AA.reverse()*AT*bti;
            bool isNan = false;
            for(int j = 0; j < 7; j++)
            {
                if ((sht(j) != sht(j)))
                {
                    isNan = true;
                }
            }
            if (!isNan)
            {
                //std::cout << sht << std::endl;
                x = sht(0)*x - sht(3)*y + sht(2)*z +sht(4);
                y = sht(3)*x + sht(0)*y - sht(1)*z +sht(5);
                z = -sht(2)*x +sht(1)*y + sht(0)*z +sht(6);
            }
        }*/
        b.AddEntry(i, 0, x);
        b.AddEntry(i + numberVertice, 0, y);
        b.AddEntry(i + 2*numberVertice, 0, z);
    }
    
    //Add handles to the vector
    for(int i = 0; i < numberConstrian; ++i)
    {
        Vertex *pV = roi_list[i];
        A.AddEntry(i + 3*numberVertice, pV->Index(), handleWeight);
        A.AddEntry(i + 3*numberVertice + numberConstrian, pV->Index() + numberVertice, handleWeight);
        A.AddEntry(i + 3*numberVertice + 2*numberConstrian, pV->Index() + 2*numberVertice, handleWeight);
    }
    
    //Solve the linear system
    A_Original = A.ToSparseMatrix(3*numberVertice + 3*numberConstrian, 3*numberVertice);
    A_Trans = A_Original.transpose();
    A_Modified = A_Trans * A_Original;
    
    solver = new SparseLinearSystemSolver(A_Modified);
    b_ = b.ToSparseMatrix(3 * numberVertice + 3 * numberConstrian, 1);
    
}
