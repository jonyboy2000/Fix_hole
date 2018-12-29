#include <list>
#include <cstring>
#include <iostream>
#include <fstream>
#include <cmath>
#include "Mesh.h"
#include "LinearSystemSolver.h"

#define pi 3.141592654


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



inline double distance(const Eigen::Vector3d& v1,const Eigen::Vector3d& v2)//两个点的距离
{
	double a;
	a = sqrt((v1(0) - v2(0))*(v1(0) - v2(0)) + (v1(1) - v2(1))*(v1(1) - v2(1)) + (v1(2) - v2(2))*(v1(2) - v2(2)));
	return a;
}

inline double Angle(const Eigen::Vector3d& p1, const Eigen::Vector3d& p2, const Eigen::Vector3d& p3)
{
	Eigen::Vector3d v1 = p1 - p2;
	Eigen::Vector3d v2 = p3 - p2;

	double _dot_res = v1.normalized().dot(v2.normalized());
	
	//Eigen::Vector3d cross = v1.normalized().cross(v2.normalized());
	//if (p2.dot(cross) < 0)
	//{
	//	return pi*2- std::acos(_dot_res);
	//}
	//else
	//{
	//	return std::acos(_dot_res);
	//}
	return std::acos(_dot_res);
}


inline double basefun(const Eigen::Vector3d& p1, const Eigen::Vector3d& p2)//basefun函数为基函数，这边取得是r^2*ln(r)，其中r是欧拉距离
{
	double l = (p1 - p2).norm();
	if (abs(l) < 0.000001)
		return 0;
	//return pow(l, 3);//r^3
	return pow(l, 2)*log(l);
}

Eigen::Vector3d Mesh::sandu(const Eigen::Vector3d& p11)//计算函数散度
{
	double x = 0;
	double y = 0;
	double z = 0;

	//for (int i = 0; i < all_list.size(); i++)//r^3
	//{
	//	x = x + 3 * (lamada_value[i] * (p11(0) - all_list[i]->Position()(0))*distance(p11, all_list[i]->Position()));
	//	y = y + 3 * (lamada_value[i] * (p11(1) - all_list[i]->Position()(1))*distance(p11, all_list[i]->Position()));
	//	z = z + 3 * (lamada_value[i] * (p11(2) - all_list[i]->Position()(2))*distance(p11, all_list[i]->Position()));
	//}
	for (int i = 0; i < all_list.size(); i++)
	{

		double m = 2 * log(distance(p11, all_list[i]->Position())) + 1.0;
		x = x + (lamada_value[i] * (p11(0) - all_list[i]->Position()(0))*m);
		y = y + (lamada_value[i] * (p11(1) - all_list[i]->Position()(1))*m);
		z = z + (lamada_value[i] * (p11(2) - all_list[i]->Position()(2))*m);
	}

	x = x + p1;
	y = y + p2;
	z = z + p3;

	
	Eigen::Vector3d newpoint = Eigen::Vector3d(x, y, z);
	return newpoint;

}
double Mesh::RBF_network(const Eigen::Vector3d& p11)//函数fx
{
	double x = 0;
	for (int i = 0; i < all_list.size(); i++)
	{
		x = x + lamada_value[i] * basefun(p11, all_list[i]->Position());

	}
	x = x + p0 + p1*p11(0) + p2*p11(1) + p3*p11(2);

	return x;
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


const FaceList& Mesh::holefzcelist() const
{
	return myFaceList;
}

const VertexList& Mesh::holevertexlist() const
{
	return myvertexList;
}

const VertexList& Mesh::bianjiedianlist() const
{
	return bianjiedian;
}

const VertexList& Mesh::normai() const {
	return normal_list;
}


// load a .obj mesh definition file
bool Mesh::LoadObjFile(const char* filename) {
    if (filename == NULL || strlen(filename) == 0) return false;
    std::ifstream ifs(filename);
    if (ifs.fail()) return false;

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
       // vList[i]->adjHEdges.clear();
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

	for (i = 0; i < myFaceList.size(); i++) delete myFaceList[i];
	for (i = 0; i < myvertexList.size(); i++) delete myvertexList[i];
	for (i = 0; i < bianjiedian.size(); i++) delete bianjiedian[i];
	for (i = 0; i < normal_list.size(); i++) delete normal_list[i];
	for (i = 0; i < all_list.size(); i++) delete all_list[i];


    heList.clear();
    bheList.clear();
    vList.clear();
    fList.clear();
	myFaceList.clear();
	myvertexList.clear();
	bianjiedian.clear();
	normal_list.clear();
	all_list.clear();
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
	int Num_vertices = vList.size();
	int Num_hedges = heList.size();
	int Num_bhedges = bheList.size();
	int Num_edges;
	int Num_fcaes = fList.size();

	int boundaries;
	int component;
	int genus;
	int Euler;

	Num_edges = (Num_hedges + Num_bhedges) / 2; 
	boundaries = CountBoundaryLoops();
	component = CountConnectedComponents();
	Euler = Num_vertices - Num_edges + Num_fcaes; 
	genus = (2 - Euler - boundaries) / 2; 

	std::cout << "Number of Vertices: " << Num_vertices << std::endl;
	std::cout << "Number of Edges: " << Num_edges << std::endl;
	std::cout << "Number of Faces: " << Num_fcaes << std::endl;
	std::cout << "Number of Boundaries: " << boundaries << std::endl;
	std::cout << "Number of Genus: " << genus << std::endl;
	std::cout << "Number of Components: " << component << std::endl;
	std::cout << "Euler characteristic: " << Euler << std::endl;
	std::cout << "bhelist" << Num_bhedges << std::endl;


}

// compute the normal of each vertex
void Mesh::ComputeVertexNormals() {
	
	const double PI = 3.14159265;
	Eigen::Vector3d t1, t2, nor, Ver1;

	int i, k, j;
	for (i = 0; i < vList.size(); i++) {

		//OneRingVertex ring(vList[i]);  // I didn't use the one ring function. it's a little bit confused especially for new comer of C++;

		k = vList[i]->Valence(); //get the degree of this vertex;
		HEdge* nextHedge;
		Vertex* nextVertex;

		t1 = Eigen::Vector3d(0, 0, 0); //define two vectors;
		t2 = Eigen::Vector3d(0, 0, 0);

		if (vList[i]->IsBoundary()) {

			nextHedge = vList[i]->HalfEdge()->Twin();  //get the first half edge and its twin, one ring vertex;
			nextVertex = nextHedge->Start();

			//make sure the starting half edge is boundary edge, so the next vertex is also on the boundary;
			//while (!(nextVertex->IsBoundary())){
			//nextHedge = nextHedge->Next()->Twin();
			//nextVertex = nextHedge->Start();
			//}

			Ver1 = nextVertex->Position(); //now, the starting vertex is denoted as Ver1;

			if (k == 2) {
				nextVertex = nextHedge->Next()->Twin()->Start();
				t1 = Ver1 - nextVertex->Position();
				t2 = Ver1 + nextVertex->Position() - vList[i]->Position();
			}
			else if (k == 3) {
				nextHedge = nextHedge->Next()->Twin();
				nextVertex = nextHedge->Start();
				t2 = nextVertex->Position() - vList[i]->Position();
				nextHedge = nextHedge->Next()->Twin();
				nextVertex = nextHedge->Start();
				t1 = Ver1 - nextVertex->Position();
			}
			else {
				for (j = 1; j < k - 1; j++) {
					nextHedge = nextHedge->Next()->Twin();
					nextVertex = nextHedge->Start();
					t2[0] += (2 * cos(PI / (k - 1)) - 2)* sin(j*PI / (k - 1))* nextVertex->Position()[0];
					t2[1] += (2 * cos(PI / (k - 1)) - 2)* sin(j*PI / (k - 1))* nextVertex->Position()[1];
					t2[2] += (2 * cos(PI / (k - 1)) - 2)* sin(j*PI / (k - 1))* nextVertex->Position()[2];

				}
				nextHedge = nextHedge->Next()->Twin();
				nextVertex = nextHedge->Start();
				t1 = Ver1 - nextVertex->Position();
				t2 = t2 + sin(PI / (k - 1))* t1;
			}

		}
		else {
			nextHedge = vList[i]->HalfEdge()->Twin();  //get the first half edge that points to this vertex;
			nextVertex = nextHedge->Start();

			for (j = 0; j < k; j++) {

				t1 += cos(2 * PI*j / k)*(nextVertex->Position() - vList[i]->Position());
				t2 += sin(2 * PI*j / k)*(nextVertex->Position() - vList[i]->Position());

				nextHedge = nextHedge->Next()->Twin(); //get next half edge that points to this vertex;
				nextVertex = nextHedge->Start();

			}

		}

		//set the normalied product of two vertors as vertex normal
		nor = -t1.cross(t2); // why it nagates itself??
		nor.normalize();
		vList[i]->SetNormal(nor);
	}


}

// compute the vertex curvature of the graph
void Mesh::ComputeVertexCurvatures() {
    /*************************/
    /* insert your code here */
    /*************************/

    /*====== Programming Assignment 1 ======*/
	
}

// umbrella smoothing
// uniformWeights: true for uniform-weight Laplacian, false for cotangent-weight Laplacian
void Mesh::UmbrellaSmooth(bool uniformWeights) {
    /*************************/
    /* insert your code here */
    /*************************/
    
	/*====== Programming Assignment 1 ======*/
    
}

// implicit umbrella smoothing
// uniformWeights: true for uniform-weight Laplacian, false for cotangent-weight Laplacian
void Mesh::ImplicitUmbrellaSmooth(bool uniformWeights) {
    /*************************/
    /* insert your code here */
    /*************************/		
	
    /*====== Programming Assignment 1 ======*/

}

int Mesh::CountBoundaryLoops() {
    /*************************/
    /* insert your code here */
    /*************************/

    /*====== Programming Assignment 0 ======*/
	return 0;
}

void Mesh::TraverseComponents(HEdge* he)
{
	return;
}


int Mesh::CountConnectedComponents() 
{
    /*************************/
    /* insert your code here */
    /*************************/

    /*====== Programming Assignment 0 ======*/
	return 0;
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
}//

double Mesh::TriArea(const Eigen::Vector3d& p1, const Eigen::Vector3d& p2, const Eigen::Vector3d& p3) {
    Eigen::Vector3d v1 = p2 - p1;
    Eigen::Vector3d v2 = p3 - p1;
    return v1.cross(v2).norm() / 2.0;
}//面积



void Mesh::Mesh_fix()
{
	int vbhe_nub = bheList.size();
	double angle[1000];
	HEdgeList boundary_list;
	boundary_list.push_back(bheList[0]);

	int nNum = 1;
	
	while (true)//boundary_list 是有顺序的bhelist 的排列，首尾相连
	{
		for (int i = 0; i < vbhe_nub; i++)
		{
			if (boundary_list[nNum - 1]->End()->Position() == bheList[i]->Start()->Position())
			{
				boundary_list.push_back(bheList[i]);
			}
		}
		nNum++;
		if (nNum == vbhe_nub)
		{
			break;
		}
	}

	int a = boundary_list.size();
	std::cout << "bhelist:  "<<a<<std::endl;

	double length_boundary_all;
	for (int i = 0; i < boundary_list.size(); i++)
	{
		length_boundary_all = length_boundary_all + distance(boundary_list[i]->Start()->Position(), boundary_list[i]->End()->Position());
	}
	double avg_length = length_boundary_all / boundary_list.size();//得到边界边的平均长度



	for (int i = 0; i < boundary_list.size()-1; i++)
	{
		angle[i] = Angle(boundary_list[i]->Start()->Position(), boundary_list[i]->End()->Position(), boundary_list[i + 1]->End()->Position());
	}
	angle[boundary_list.size() -1] = Angle(boundary_list[boundary_list.size() -1]->Start()->Position(), boundary_list[boundary_list.size() -1]->End()->Position(), boundary_list[0]->End()->Position());

	int min_angle = 0;
	for (int i = 1; i < boundary_list.size(); i++)
	{
		if (angle[min_angle] > angle[i])
		{
			min_angle = i;//得到边界点中最小角度的序号（在boundary_list中）
		}
	}

	Eigen::Vector3d p1;
	Eigen::Vector3d p2;
	Eigen::Vector3d p3;

	if (min_angle == vbhe_nub - 1)
	{
		 p1 = boundary_list[min_angle]->Start()->Position();
		 p2 = boundary_list[min_angle]->End()->Position();
		 p3 = boundary_list[0]->End()->Position();
	}
	else
	{
		 p1 = boundary_list[min_angle]->Start()->Position();
		 p2 = boundary_list[min_angle]->End()->Position();
		 p3 = boundary_list[min_angle+1]->End()->Position();
	}

	double tri_dis = distance(p1, p3);
	if (tri_dis > 2 * avg_length)//如果最小角度点的两边点连接线长度大于平均长度，加上一个点，此点是连接线中点
	{
		Eigen::Vector3d p13 = (p1+p3)/2;
		AddVertex(new Vertex(p13));
		vList[vList.size() - 1]->SetIndex(vList.size() - 1);
		AddFace(boundary_list[min_angle]->Start()->Index(), boundary_list[min_angle]->End()->Index(), vList[vList.size() - 1]->Index());
		if (min_angle == vbhe_nub - 1)//判断是不是boundary_list序号中最后一个点
		{
			AddFace(boundary_list[min_angle]->End()->Index(), boundary_list[0]->End()->Index(), vList[vList.size() - 1]->Index());
		}
		else
		{
			AddFace(boundary_list[min_angle]->End()->Index(), boundary_list[min_angle + 1]->End()->Index(), vList[vList.size() - 1]->Index());
		}
		HEdgeList list;
		for (int i = 0; i < bheList.size(); i++)
			if (bheList[i]->Start()) list.push_back(bheList[i]);
		bheList = list;//新的边界边
	}
	else//如果最小角度点的两边点连接线长度小于平均长度，直接连接即可
	{
		if (min_angle == vbhe_nub - 1)
		{
			AddFace(boundary_list[min_angle]->Start()->Index(), boundary_list[min_angle]->End()->Index(),boundary_list[0]->End()->Index());
		}
		else
		{
			AddFace(boundary_list[min_angle]->Start()->Index(), boundary_list[min_angle]->End()->Index(), boundary_list[min_angle + 1]->End()->Index());
		}
		HEdgeList list;
		for (int i = 0; i < bheList.size(); i++)
			if (bheList[i]->Start()) list.push_back(bheList[i]);
		bheList = list;
	}

}

void Mesh::Least_squares_1()
{

	int m = fList.size();//未修补时的facelist的数量
	int d = vList.size();//未修补时的vlistlist的数量

	for (int i = 0; i < bheList.size(); i++)
	{
			bianjiedian.push_back(bheList[i]->End());
	}

	int bj_ver = bianjiedian.size();


	while (bheList.size() > 0)
	{
		Mesh_fix();
	}


	int n = fList.size();// 修补后的facelist的数量
	for (int i = m; i < n; i++)
	{
		myFaceList.push_back(fList[i]);
	}

	int e = vList.size();//修补后的vlistlist的数量
	for (int i = d; i < e; i++)
	{
		myvertexList.push_back(vList[i]);
	}


	Eigen::VectorXd deta_Vertex_X = Eigen::VectorXd::Zero(e + bj_ver);
	Eigen::VectorXd deta_Vertex_Y = Eigen::VectorXd::Zero(e + bj_ver);
	Eigen::VectorXd deta_Vertex_Z = Eigen::VectorXd::Zero(e + bj_ver);


	for (int i = 0; i < e; i++)
	{
		deta_Vertex_X[i] = 0;
		deta_Vertex_Y[i] = 0;
		deta_Vertex_Z[i] = 0;
	}
	for (int i = e; i < e + bj_ver; i++)
	{
		deta_Vertex_X[i] = bianjiedian[i - e]->Position()(0);
		deta_Vertex_Y[i] = bianjiedian[i - e]->Position()(1);
		deta_Vertex_Z[i] = bianjiedian[i - e]->Position()(2);
	}


	Eigen::SparseMatrix<double> L = Eigen::SparseMatrix<double>(e + bj_ver,e);
	std::vector<Eigen::Triplet<double>> triple;

	for (int i = 0; i < e; i++)
	{
		int k = vList[i]->Valence();
		double avg = 1 / double(k);
		HEdge* nextHedge = vList[i]->HalfEdge()->Twin();

		for (int r = 0; r < k; r++) {
			nextHedge = nextHedge->Next()->Twin();
		}

		// compute weight for each one-ring vertex		
		for (int r = 0; r < k; r++) {

			triple.push_back(Eigen::Triplet<double>(i, nextHedge->Start()->Index(), -avg));
			nextHedge = nextHedge->Next()->Twin();
		}
		triple.push_back(Eigen::Triplet<double>(i, i, 1.0));
	}
	for (int i = e; i < e + bj_ver; i++)
	{

		triple.push_back(Eigen::Triplet<double>(i, bianjiedian[i-e]->Index(), 1.0));
	}

	L.setFromTriplets(triple.begin(), triple.end());
	Eigen::SparseLU < Eigen::SparseMatrix<double>> solver;
	solver.compute(L.transpose()*L);

	Eigen::VectorXd vx2 = Eigen::VectorXd::Zero(e);
	Eigen::VectorXd vy2 = Eigen::VectorXd::Zero(e);
	Eigen::VectorXd vz2 = Eigen::VectorXd::Zero(e);

	vx2 = solver.solve(L.transpose()*deta_Vertex_X);
	vy2 = solver.solve(L.transpose()*deta_Vertex_Y);
	vz2 = solver.solve(L.transpose()*deta_Vertex_Z);
	//set the position of vertices
	for (int i = d; i < e; i++) 
	{
		Eigen::Vector3d newPosition = Eigen::Vector3d(vx2[i], vy2[i], vz2[i]);
		vList[i]->SetPosition(newPosition);
	}

}


void Mesh::Mesh_RBF_fix()
{
	
	Least_squares_1();//得到最小二乘网格点
	int bj_ver = bianjiedian.size();

	ComputeVertexNormals();//求得法线，这边的法线是有了最小二乘网格之后，求得边界点的法线
	for (int i = 0; i < bj_ver; i++)
	{
		int n = bianjiedian[i]->Index();

		Eigen::Vector3d ne = 0.01*vList[n]->Normal() + vList[n]->Position();
		normal_list.push_back(new Vertex(ne));
	}


//all_list是所有点的集合，这些点包括(有顺序)边界点外法线方向的外部点，边界点，边界点内法线方向的内部点，长度为0.01
	for (int i = 0; i < bj_ver; i++)
	{
		int n = bianjiedian[i]->Index();
		Eigen::Vector3d wai = 0.01*vList[n]->Normal() + vList[n]->Position();

		all_list.push_back(new Vertex(wai));//先增加所有的外部点，
	}
	for (int i = 0; i < bj_ver; i++)
	{
		int n = bianjiedian[i]->Index();
		Eigen::Vector3d origin = vList[n]->Position();

		all_list.push_back(new Vertex(origin));//再曾加边界点
	}
	for (int i = 0; i < bj_ver; i++)
	{
		int n = bianjiedian[i]->Index();
		Eigen::Vector3d nei = -0.01*vList[n]->Normal() + vList[n]->Position();

		all_list.push_back(new Vertex(nei));//最后增加内部点
	}



	int k = bj_ver * 3;
	Eigen::MatrixXd L_new = Eigen::MatrixXd::Zero(k + 4, k + 4);


	for (int i = 0; i < k; i++)
	{
		for (int j = 0; j < k; j++)
		{
			L_new(i, j) = basefun(all_list[i]->Position(), all_list[j]->Position());
		}
	}
	for (int i = 0; i < k - 1; i++)
	{
		L_new(i, k) = 1;
		L_new(i, k + 1) = all_list[i]->Position()(0);
		L_new(i, k + 2) = all_list[i]->Position()(1);
		L_new(i, k + 3) = all_list[i]->Position()(2);

		L_new(k, i) = 1;
		L_new(k + 1, i) = all_list[i]->Position()(0);
		L_new(k + 2, i) = all_list[i]->Position()(1);
		L_new(k + 3, i) = all_list[i]->Position()(2);
	}

	Eigen::VectorXd fx = Eigen::VectorXd::Zero(3 * bj_ver + 4);
	for (int i = 0; i < bj_ver; i++)
	{
		double a = 1.0;
		fx[i] = -1;                  //外部点取fx=-1
		fx[i + bj_ver] = 0.0;        //原点取fx=0
		fx[i + 2 * bj_ver] = 1;      //内部点取fx=1
	}


	//lamada_p是[λ1,λ2,λ3……λk，p0,p1,p2,p3].tranpose()。λn是每个基函数的系数
	Eigen::VectorXd lamada_p = Eigen::VectorXd::Zero(k + 4);


	//解方程L_new*lamada_p=fx
	lamada_p = L_new.lu().solve(fx);


	lamada_value = lamada_p;//得到λ的值的vector，只要前k个值

	std::cout << "lamada" << lamada_p << std::endl;


	p0 = lamada_p[k];
	p1 = lamada_p[k + 1];
	p2 = lamada_p[k + 2];
	p3 = lamada_p[k + 3];


	int count = 0;
	double error = 0;//误差值
	double jingdu = 0.000000001;//定义精度

	do{

		for (int i = vList.size() - myvertexList.size(); i < vList.size(); i++)
		{
			Eigen::Vector3d newPosition = Eigen::Vector3d::Zero(3);
			Eigen::Vector3d oldPosition = vList[i]->Position();
			Eigen::Vector3d sandu_vec = sandu(oldPosition);


			newPosition = oldPosition - RBF_network(oldPosition)*sandu_vec / pow(sandu_vec.norm(), 2);
			vList[i]->SetPosition(newPosition);
		}

		for (int i = vList.size() - myvertexList.size(); i < vList.size(); i++)
		{
			Eigen::Vector3d oldPosition = vList[i]->Position();
			double x = RBF_network(oldPosition);
			error = error + x;
		}
		error = error / myvertexList.size();
		count++;
	} while (abs(error) > jingdu);

	std::cout << "count" << count << std::endl;

}