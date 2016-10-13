#include <list>
#include <cstring>
#include <iostream>
#include <fstream>
#include <cmath>
#include "Mesh.h"
#include "LinearSystemSolver.h"
#include <queue>

inline Eigen::Vector3d MinVector3d(Eigen::Vector3d v1, Eigen::Vector3d v2) {
    return Eigen::Vector3d(std::min(v1(0), v2(0)),
                           std::min(v1(1), v2(1)),
                           std::min(v1(2), v2(2)));
}

inline Eigen::Vector3d MaxVector3d(Eigen::Vector3d v1, Eigen::Vector3d v2) {
    return Eigen::Vector3d(std::max(v1(0), v2(0)),
                           std::max(v1(1), v2(1)),
                           std::max(v1(2), v2(2)));
}

OneRingHEdge::OneRingHEdge(const Vertex* v) {
    if (v == NULL) start = next = NULL;
    else start = next = v->HalfEdge();
}

HEdge* OneRingHEdge::NextHEdge() {
    HEdge* ret = next;
    if (next && next->Prev()->Twin() != start)
        next = next->Prev()->Twin();
    else
        next = NULL;
    return ret;
}

Mesh::~Mesh() {
    Clear();
}

const HEdgeList& Mesh::Edges() const {
    return heList;
}

const HEdgeList& Mesh::BoundaryEdges() const {
    return bheList;
}

const VertexList& Mesh::Vertices() const {
    return vList;
}

const FaceList& Mesh::Faces() const {
    return fList;
}

// load a .obj mesh definition file
bool Mesh::LoadObjFile(const char* filename) {
    if (filename == NULL || strlen(filename) == 0) return false;
    std::ifstream ifs(filename);
    if (ifs.fail())
    {
        return false;
    }
    Clear();

    std::string line;
    while (std::getline(ifs, line)) {
        if (line.empty()) continue;

        std::istringstream iss(line);
        std::string type;
        iss >> type;
        // vertex
        if (type.compare("v") == 0) {
            double x, y, z;
            iss >> x >> y >> z;
            AddVertex(new Vertex(x, y, z));
        }
        // face
        else if (type.compare("f") == 0) {
            int index[3];
            iss >> index[0] >> index[1] >> index[2];
            AddFace(index[0] - 1, index[1] - 1, index[2] - 1);
        }
    }
    ifs.close();

    size_t i;
    Eigen::Vector3d box = this->MaxCoord() - this->MinCoord();
    for (i = 0; i < vList.size(); i++) vList[i]->SetPosition(vList[i]->Position() / box(0));

    Eigen::Vector3d tot = Eigen::Vector3d::Zero();
    for (i = 0; i < vList.size(); i++) tot += vList[i]->Position();
    Eigen::Vector3d avg = tot / vList.size();
    for (i = 0; i < vList.size(); i++) vList[i]->SetPosition(vList[i]->Position() - avg);

    HEdgeList list;
    for (i = 0; i < bheList.size(); i++)
        if (bheList[i]->Start()) list.push_back(bheList[i]);
    bheList = list;

    for (i = 0; i < vList.size(); i++) {
        vList[i]->adjHEdges.clear();
        vList[i]->SetIndex((int)i);
        vList[i]->SetFlag(0);
    }

    return true;
}

void Mesh::AddVertex(Vertex* v) {
    vList.push_back(v);
}

void Mesh::AddFace(int v1, int v2, int v3) {
    int i;
    HEdge *he[3], *bhe[3];
    Vertex* v[3];
    Face* f;

    // obtain objects
    for (i = 0; i < 3; i++) he[i] = new HEdge();
    for (i = 0; i < 3; i++) bhe[i] = new HEdge(true);
    v[0] = vList[v1];
    v[1] = vList[v2];
    v[2] = vList[v3];
    f = new Face();

    // connect prev-next pointers
    SetPrevNext(he[0], he[1]);
    SetPrevNext(he[1], he[2]);
    SetPrevNext(he[2], he[0]);
    SetPrevNext(bhe[0], bhe[1]);
    SetPrevNext(bhe[1], bhe[2]);
    SetPrevNext(bhe[2], bhe[0]);

    // connect twin pointers
    SetTwin(he[0], bhe[0]);
    SetTwin(he[1], bhe[2]);
    SetTwin(he[2], bhe[1]);

    // connect start pointers for bhe
    bhe[0]->SetStart(v[1]);
    bhe[1]->SetStart(v[0]);
    bhe[2]->SetStart(v[2]);
    for (i = 0; i < 3; i++) he[i]->SetStart(v[i]);

    // connect start pointers
    // connect face-hedge pointers
    for (i = 0; i < 3; i++) {
        v[i]->SetHalfEdge(he[i]);
        v[i]->adjHEdges.push_back(he[i]);
        SetFace(f, he[i]);
    }
    v[0]->adjHEdges.push_back(bhe[1]);
    v[1]->adjHEdges.push_back(bhe[0]);
    v[2]->adjHEdges.push_back(bhe[2]);

    // mearge boundary if in need
    for (i = 0; i < 3; i++) {
        Vertex* start = bhe[i]->Start();
        Vertex* end = bhe[i]->End();
        for (size_t j = 0; j < end->adjHEdges.size(); j++) {
            HEdge* curr = end->adjHEdges[j];
            if (curr->IsBoundary() && curr->End() == start) {
                SetPrevNext(bhe[i]->Prev(), curr->Next());
                SetPrevNext(curr->Prev(), bhe[i]->Next());
                SetTwin(bhe[i]->Twin(), curr->Twin());
                bhe[i]->SetStart(NULL); // mark as unused
                curr->SetStart(NULL); // mark as unused
                break;
            }
        }
    }

    // finally add hedges and faces to list
    for (i = 0; i < 3; i++) heList.push_back(he[i]);
    for (i = 0; i < 3; i++) bheList.push_back(bhe[i]);
    fList.push_back(f);
}

void Mesh::Clear() {
    size_t i;
    for (i = 0; i < heList.size(); i++) delete heList[i];
    for (i = 0; i < bheList.size(); i++) delete bheList[i];
    for (i = 0; i < vList.size(); i++) delete vList[i];
    for (i = 0; i < fList.size(); i++) delete fList[i];
    heList.clear();
    bheList.clear();
    vList.clear();
    fList.clear();
}

Eigen::Vector3d Mesh::MinCoord() const {
    Eigen::Vector3d minCoord = Eigen::Vector3d::Zero();
    for (size_t i = 0; i < vList.size(); i++)
        minCoord = MinVector3d((vList[i])->Position(), minCoord);
    return minCoord;
}

Eigen::Vector3d Mesh::MaxCoord() const {
    Eigen::Vector3d maxCoord = Eigen::Vector3d::Zero();
    for (size_t i = 0; i < vList.size(); i++)
        maxCoord = MaxVector3d((vList[i])->Position(), maxCoord);
    return maxCoord;
}

void Mesh::DisplayMeshInfo() {
    /*************************/
    /* insert your code here */
    /*************************/

    /*====== Programming Assignment 0 ======*/
    int NumberBoundary, NumberComponents;
    NumberBoundary = CountBoundaryLoops();
    NumberComponents = CountConnectedComponents();
    std::cout << "*====== Programming Assignment 0 ======*" << std::endl;
    std::cout << "#1. Number of Vertices: " << vList.size() << std::endl;
    std::cout << "#2. Number of Half-Edges: " << heList.size() + bheList.size() << std::endl;
    std::cout << "#3. Number of Faces: " << fList.size() << std::endl;
    std::cout << "#4. Number of Loop Boundary: " << NumberBoundary << std::endl;
    std::cout << "#5. Number of Connect Components: " << NumberComponents << std::endl;
    std::cout << "#6. Number of Genus: " <<NumberComponents - (vList.size() - (heList.size() + bheList.size())/2 +fList.size() + NumberBoundary)/2 <<  std::endl;
}

void Mesh::ComputeVertexNormals() {
    /*************************/
    /* insert your code here */
    /*************************/
    
    /*====== Programming Assignment 0 ======*/
    for (int i = 0; i < vList.size(); i++)
    {
        Vertex* p = vList[i];
        int valence = p->Valence();
        Eigen::Vector3d t1 = Eigen::Vector3d::Zero();
        Eigen::Vector3d t2 = Eigen::Vector3d::Zero();
        std::vector<Vertex*> pSet;
        OneRingVertex ring = OneRingVertex(p);
        
        for (int j = 0; j < valence; j++)
        {
            pSet.push_back(ring.NextVertex());
        }
        if (!p->IsBoundary())
        {
            for (int j = 0; j < valence; j++)
            {
                //curr = ring.NextVertex();
                t1 += cos((3.14*2*j)/valence) * pSet[j]->Position();
                t2 += sin((3.14*2*j)/valence) * pSet[j]->Position();
            }
        }
        else
        {
            t1 = pSet[0]->Position() - pSet[valence-1]->Position();
            if (valence == 2)
            {
                t2 = pSet[0]->Position() + pSet[1]->Position() - p->Position();
            }
            if (valence == 3)
            {
                t2 = pSet[1]->Position() - p->Position();
            }
            if (valence > 3)
            {
                double theta = 3.14/(valence-1);
                Eigen::Vector3d temp = Eigen::Vector3d::Zero();
                for (int j = 0; j < valence - 2; j++)
                {
                    temp += sin(j*theta) * pSet[j]->Position();
                }
                t2 = sin(theta) * (pSet[0]->Position() + pSet[valence-1]->Position()) + 2*(cos(theta) - 1) * temp;
            }
        }
        Eigen::Vector3d normal = t1.cross(t2);
        normal.normalize();
        p->SetNormal(normal);
    }
}

// compute the vertex curvature of the graph
double Mesh::CalculateSpace(Vertex* p, Vertex* p1, Vertex* p2)
{
    double space = 0;
    Eigen::Vector3d t1,t2,t3;
    t1 = p1->Position() - p->Position();
    t2 = p2->Position() - p->Position();
    t3 = p1->Position() - p2->Position();
    double d1 = t1.norm();
    double d2 = t2.norm();
    double d3 = t3.norm();
    double halfLength = (d1 + d2 + d3)/2;
    space = sqrt(halfLength * (halfLength - d1) * (halfLength - d2) * (halfLength - d3));
    return space;
}

double Mesh::CalculateCot(Vertex* p, Vertex* pm, Vertex* pj)
{
    Eigen::Vector3d t1 = p->Position() - pm->Position();
    Eigen::Vector3d t2 = pj->Position() - pm->Position();
    double d1 = t1.norm();
    double d2 = t2.norm();
    double cos = t1.dot(t2) / (d1*d2);
    double sin = sqrt(1-cos*cos);
    return cos/sin;
}

void Mesh::MeanCurvatureNormalise(std::vector<double>& meanList)
{
    double min = 100000000;
    double max = 0;
    double length;
    for (int i = 0; i < meanList.size(); i++)
    {
        meanList[i] = log2(meanList[i]+1) + 0.1;
    }
    for (int i = 0; i < meanList.size(); i++)
    {
        if (meanList[i] > max)
        {
            max = meanList[i];
        }
        if (meanList[i] < min)
        {
            min = meanList[i];
        }
    }
    length = max - min;
    for (int i = 0; i < meanList.size(); i++)
    {
        meanList[i] = (((meanList[i] - min)/length)-0.5)*2;
    }
}
// used for produce colormap;
double interpolate( double val, double y0, double x0, double y1, double x1 ) {
    return (val-x0)*(y1-y0)/(x1-x0) + y0;
}

double base( double val ) {
    if ( val <= -0.75 ) return 0;
    else if ( val <= -0.25 ) return interpolate( val, 0.0, -0.75, 1.0, -0.25 );
    else if ( val <= 0.25 ) return 1.0;
    else if ( val <= 0.75 ) return interpolate( val, 1.0, 0.25, 0.0, 0.75 );
    else return 0.0;
}

double red( double gray ) {
    return base( gray - 0.5 );
}
double green( double gray ) {
    return base( gray );
}
double blue( double gray ) {
    return base( gray + 0.5 );
}

void Mesh::ComputeVertexCurvatures() {
    /*************************/
    /* insert your code here */
    /*************************/

    /*====== Programming Assignment 1 ======*/
    //calculate space
    std::vector<double> meanList;
    for (int i = 0; i < vList.size(); i++)
    {
        Vertex* p = vList[i];
        int valence = vList[i]->Valence();
        std::vector<Vertex*> ringVertex;
        OneRingVertex ring = OneRingVertex(p);
        for (int j = 0; j < valence; j++)
        {
            ringVertex.push_back(ring.NextVertex());
        }
        
        //Calculate Space
        double space = 0;
        for (int j = 0; j < valence-1; j++)
        {
            space += CalculateSpace(p, ringVertex[j], ringVertex[j+1]);
        }
        if (!p->IsBoundary())
        {
            space += CalculateSpace(p, ringVertex[valence-1], ringVertex[0]);
        }
        
        //Calculate Cot
        Eigen::Vector3d sumCot(0,0,0);
        for (int j = 1; j < valence - 1; j++)
        {
            double cot = CalculateCot(p, ringVertex[j-1], ringVertex[j]) + CalculateCot(p, ringVertex[j+1], ringVertex[j]);
            sumCot += cot * (ringVertex[j]->Position() - p->Position());
        }
        if (!p->IsBoundary())
        {
            double cot1 = CalculateCot(p, ringVertex[1], ringVertex[0]) + CalculateCot(p, ringVertex[valence-1], ringVertex[0]);
            double cot2 = CalculateCot(p, ringVertex[valence-2], ringVertex[valence-1]) + CalculateCot(p, ringVertex[0], ringVertex[valence-1]);
            sumCot += (cot1 * (ringVertex[0]->Position() - p->Position()) + cot2 * (ringVertex[valence-1]->Position() - p->Position()));
        }
        sumCot = -(1/(4*space))*sumCot;
        double meanCurve = sumCot.norm();
        meanList.push_back(meanCurve);
    }
    MeanCurvatureNormalise(meanList);
    
    //Five color heat map;
    for (int i = 0; i < vList.size(); i++)
    {
        Eigen::Vector3d Color = Eigen::Vector3d(red(meanList[i]), green(meanList[i]), blue(meanList[i]));
        vList[i]->SetColor(Color);
    }
}

// umbrella smoothing
// uniformWeights: true for uniform-weight Laplacian, false for cotangent-weight Laplacian
Eigen::SparseMatrix<double> Mesh::BuilMatrix(bool uniformWeights)
{
    SparseMatrixBuilder matrixBuilder;
    if (uniformWeights)
    {
        for (int i = 0; i < vList.size(); i++)
        {
            Vertex* p = vList[i];
            if (!p->IsBoundary())
            {
                int valence = p->Valence();
                double weight = valence;
                weight = 1/weight;
                OneRingVertex ring = OneRingVertex(p);
                Vertex* curr = NULL;
                while(curr = ring.NextVertex())
                {
                    matrixBuilder.AddEntry(i, curr->Index(), weight);
                }
                matrixBuilder.AddEntry(i, i, -1);
            }
        }
    }
    else
    {
        for (int i = 0; i < vList.size(); i++)
        {
            Vertex* p = vList[i];
            if (!p->IsBoundary())
            {
                int valence = p->Valence();
                std::vector<Vertex*> pSet;
                std::vector<double> weightSet;
                double weightSum = 0;
                OneRingVertex ring = OneRingVertex(p);
                Vertex* curr = NULL;
                while(curr = ring.NextVertex())
                {
                    pSet.push_back(curr);
                }
                //test
                weightSet.push_back(CalculateCot(p, pSet[1], pSet[0]) + CalculateCot(p, pSet[valence-1], pSet[0]));
                for (int j = 1; j < valence - 1; j++)
                {
                    weightSet.push_back(CalculateCot(p, pSet[j+1], pSet[j]) + CalculateCot(p, pSet[j-1], pSet[j]));
                    weightSum += weightSet[j];
                }
                weightSet.push_back(CalculateCot(p, pSet[0], pSet[valence-1]) + CalculateCot(p, pSet[valence-2], pSet[valence-1]));
                weightSum += (weightSet[0] + weightSet[valence-1]);
                for (int j = 0; j < valence; j++)
                {
                    matrixBuilder.AddEntry(i, pSet[j]->Index(), weightSet[j]/weightSum);
                }
                matrixBuilder.AddEntry(i, i, -1);
            }
        }
    }
    return matrixBuilder.ToSparseMatrix(vList.size(), vList.size());
}


void Mesh::UmbrellaSmooth(bool uniformWeights) {
    /*************************/
    /* insert your code here */
    /*************************/

    /*====== Programming Assignment 1 ======*/
    uniformWeights = 0;
    double Lambda = 1;
    Eigen::SparseMatrix<double> Laplace(vList.size(),vList.size());
    Laplace = BuilMatrix(uniformWeights);
    //std::cout << Laplace;
    Eigen::Matrix<double, Eigen::Dynamic, 3> m1;
    m1.resize(vList.size(), 3);
    for (int i =0; i < vList.size(); i++)
    {
        m1(i,0) = vList[i]->Position().x();
        m1(i,1) = vList[i]->Position().y();
        m1(i,2) = vList[i]->Position().z();
    }
    
    m1 = m1 + Lambda*(Laplace*m1);
    for (int i = 0; i < vList.size(); i++)
    {
        Eigen::Vector3d position(m1(i,0),m1(i,1),m1(i,2));
        vList[i]->SetPosition(position);
    }
    ComputeVertexNormals();
    ComputeVertexCurvatures();
}

// implicit umbrella smoothing
// uniformWeights: true for uniform-weight Laplacian, false for cotangent-weight Laplacian

Eigen::SparseMatrix<double> Mesh::BuilMatrix(bool uniformWeights, double Lambda)
{
    SparseMatrixBuilder matrixBuilder;
    if (uniformWeights)
    {
        for (int i = 0; i < vList.size(); i++)
        {
            Vertex* p = vList[i];
            if (!p->IsBoundary())
            {
                int valence = p->Valence();
                double weight = valence;
                weight = 1/weight;
                OneRingVertex ring = OneRingVertex(p);
                Vertex* curr = NULL;
                while(curr = ring.NextVertex())
                {
                    matrixBuilder.AddEntry(i, curr->Index(), -Lambda*weight);
                }
                matrixBuilder.AddEntry(i, i, 1+Lambda*1);
            }
            else
            {
                matrixBuilder.AddEntry(i, i, 1);
            }
        }
    }
    else
    {
        for (int i = 0; i < vList.size(); i++)
        {
            Vertex* p = vList[i];
            if (!p->IsBoundary())
            {
                int valence = p->Valence();
                std::vector<Vertex*> pSet;
                std::vector<double> weightSet;
                double weightSum = 0;
                OneRingVertex ring = OneRingVertex(p);
                Vertex* curr = NULL;
                while(curr = ring.NextVertex())
                {
                    pSet.push_back(curr);
                }
                //test
                weightSet.push_back(CalculateCot(p, pSet[1], pSet[0]) + CalculateCot(p, pSet[valence-1], pSet[0]));
                for (int j = 1; j < valence - 1; j++)
                {
                    weightSet.push_back(CalculateCot(p, pSet[j+1], pSet[j]) + CalculateCot(p, pSet[j-1], pSet[j]));
                    weightSum += weightSet[j];
                }
                weightSet.push_back(CalculateCot(p, pSet[0], pSet[valence-1]) + CalculateCot(p, pSet[valence-2], pSet[valence-1]));
                weightSum += (weightSet[0] + weightSet[valence-1]);
                for (int j = 0; j < valence; j++)
                {
                    matrixBuilder.AddEntry(i, pSet[j]->Index(), -Lambda*weightSet[j]/weightSum);
                }
                matrixBuilder.AddEntry(i, i, 1+Lambda);
            }
            else
            {
                matrixBuilder.AddEntry(i, i, 1);
            }
        }
    }
    return matrixBuilder.ToSparseMatrix(vList.size(), vList.size());
}

void Mesh::ImplicitUmbrellaSmooth(bool uniformWeights) {
    double Lambda = 1;
    uniformWeights = 0;
    Eigen::SparseMatrix<double> Laplace(vList.size(),vList.size());
    Laplace = BuilMatrix(uniformWeights, Lambda);
    //std::cout << Laplace;
    ComputeVertexNormals();
    ComputeVertexCurvatures();
    Eigen::VectorXd x[3],b[3];
    
    for(int i = 0; i < 3; i++)
    {
        x[i].resize(vList.size(), 1);
        b[i].resize(vList.size(), 1);
    }
    for(int i = 0; i < vList.size(); i++)
    {
        for(int j = 0; j < 3; j++)
        {
            b[j](i,0) = vList[i]->Position()(j);
            x[j](i,0) = 0;
        }
    }
    SparseLinearSystemSolver solver(Laplace);
    for(int i = 0; i < 3; i++)
    {
        solver.SolveByConjugateGradient(Laplace, b[i], 200, 0.0001, x[i]);
        std::cout << (Laplace*x[i] - b[i]).norm() << std::endl;
    }
    for(int i = 0; i < vList.size(); i++)
    {
        Eigen::Vector3d position(x[0](i,0),x[1](i,0),x[2](i,0));
        vList[i]->SetPosition(position);
    }
    ComputeVertexNormals();
    ComputeVertexCurvatures();
}

void Mesh::LoopBhe(int index)
{
    bheList[index]->SetFlag(1);
    HEdge* curr = bheList[index]->Next();
    while (curr != bheList[index])
    {
        curr->SetFlag(1);
        curr = curr->Next();
    }
}

void Mesh::SetBheBack()
{
    for (int i = 0; i < bheList.size(); i++)
    {
        bheList[i]->SetFlag(0);
    }
}

int Mesh::CountBoundaryLoops() {
    /*************************/
    /* insert your code here */
    /*************************/

    /*====== Programming Assignment 0 ======*/
    int NumberLoops = 0;
    for (int i = 0; i < bheList.size(); i++)
    {
        if (!bheList[i]->Flag())
        {
            NumberLoops++;
            LoopBhe(i);
        }
    }
    SetBheBack();
    return NumberLoops;
}

void Mesh::BFSVertex(int index)
{
    std::queue<Vertex*> VertexQueue;
    VertexQueue.push(vList[index]);
    while(!VertexQueue.empty())
    {
        VertexQueue.front()->SetFlag_(1);
        OneRingVertex ring = OneRingVertex(VertexQueue.front());
        Vertex* curr = NULL;
        while (curr = ring.NextVertex())
        {
            if (!curr->Flag_())
            {
                VertexQueue.push(curr);
            }
        }
        VertexQueue.pop();
    }
}

void Mesh::SetVertexBack()
{
    for (int i = 0; i < vList.size(); i++)
    {
        vList[i]->SetFlag_(0);
    }
}

int Mesh::CountConnectedComponents() {
    /*************************/
    /* insert your code here */
    /*************************/

    /*====== Programming Assignment 0 ======*/
    int NumberComponents = 0;
    for (int i = 0; i < vList.size(); i++)
    {
        if (!vList[i]->Flag_())
        {
            NumberComponents ++;
            BFSVertex(i);
        }
    }
    SetVertexBack();
    return NumberComponents;
}

void Mesh::GroupingVertexFlags() {
    // set vertex flag to be 255 initially
    for (size_t i = 0; i < vList.size(); i++)
        if (vList[i]->Flag() != 0)
            vList[i]->SetFlag(255);

    int id = 0;
    VertexList tmpList;
    for (int i = 0; i < vList.size(); i++)
        if (vList[i]->Flag() == 255) {
            id++;
            vList[i]->SetFlag(id);
            tmpList.push_back(vList[i]);
            while (! tmpList.empty()) {
                Vertex* v = tmpList.back();
                tmpList.pop_back();
                OneRingVertex ring = OneRingVertex(v);
                while (Vertex* v2 = ring.NextVertex()) {
                    if (v2->Flag() == 255) {
                        v2->SetFlag(id);
                        tmpList.push_back(v2);
                    }
                }
            }
        }
}


void Mesh::SetPrevNext(HEdge* e1, HEdge* e2) {
    e1->SetNext(e2);
    e2->SetPrev(e1);
}

void Mesh::SetTwin(HEdge* e1, HEdge* e2) {
    e1->SetTwin(e2);
    e2->SetTwin(e1);
}

void Mesh::SetFace(Face* f, HEdge* e) {
    f->SetHalfEdge(e);
    e->SetFace(f);
}

double Mesh::Cot(const Eigen::Vector3d& p1, const Eigen::Vector3d& p2, const Eigen::Vector3d& p3) {
    Eigen::Vector3d v1 = p1 - p2;
    Eigen::Vector3d v2 = p3 - p2;

    double _dot_res = v1.normalized().dot(v2.normalized());
    if (_dot_res < -1.0) {
        _dot_res = -1.0;
    }
    else if (_dot_res >  1.0) {
        _dot_res = 1.0;
    }
    return 1.0 / std::tan(std::acos(_dot_res));
}

double Mesh::TriArea(const Eigen::Vector3d& p1, const Eigen::Vector3d& p2, const Eigen::Vector3d& p3) {
    Eigen::Vector3d v1 = p2 - p1;
    Eigen::Vector3d v2 = p3 - p1;
    return v1.cross(v2).norm() / 2.0;
}
