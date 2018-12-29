/**
* Mesh 
*
* using half-edge data structure
*/

#ifndef __MESH_H__
#define __MESH_H__
#include <vector>

// Use the library Eigen for linear algebra operations
// Documentation: http://eigen.tuxfamily.org/dox/
#include "Eigen/Dense"

// classes
class HEdge;
class Vertex;
class Face;
class Mesh;
class OneRingHEdge;
class OneRingVertex;

// types
typedef std::vector< HEdge* > HEdgeList;
typedef std::vector< Vertex* > VertexList;
typedef std::vector< Face* > FaceList;

////////// class HEdge //////////
class HEdge {
private:
    HEdge *twin, *prev, *next; // twin/previous/next half edges
    Vertex* start; // start vertex
    Face* face; // left face
    bool boundary; // flag for boundary edge

    bool flag; // it's for count boundary loop use, 
    // you can use it freely for marking in boundary loop counting 
    // and connected component counting

    bool valid;

public:
    /////////////////////////////////////
    // constructor
    HEdge(bool b = false) {
        boundary = b;
        twin = prev = next = NULL;
        start = NULL;
        face = NULL;
        flag = false;
        valid = true;
    }

    HEdge* Twin() const {
        return twin;
    }

    HEdge* Prev() const {
        return prev;
    }

    HEdge* Next() const {
        return next;
    }

    Vertex* Start() const {
        return start;
    }

    Vertex* End() const {
        return next->start;
    } // for convenience

    Face* LeftFace() const {
        return face;
    }

    bool Flag() const {
        return flag;
    }

    HEdge* SetTwin(HEdge* e) {
        return twin = e;
    }

    HEdge* SetPrev(HEdge* e) {
        return prev = e;
    }

    HEdge* SetNext(HEdge* e) {
        return next = e;
    }

    Vertex* SetStart(Vertex* v) {
        return start = v;
    }

    Face* SetFace(Face* f) {
        return face = f;
    }

    bool SetFlag(bool b) {
        return flag = b;
    }

    bool SetValid(bool b) {
        return valid = b;
    }

    bool IsBoundary() const {
        return boundary;
    }

    bool IsValid() const {
        return valid;
    }
};


////////// class OneRingHEdge //////////
// this class is use for access the neighbor HALF EDGES
// of a given vertex, please see Vertex::IsBoundary() for its usage
class OneRingHEdge {
private:
    HEdge *start, *next;
public:
    OneRingHEdge(const Vertex* v); // constructor
    HEdge* NextHEdge(); // iterator
};


////////// class OneRingVertex //////////
// this class is use for access the neighbor VERTICES 
// of a given vertex, please see Vertex::Valence() for its usage
class OneRingVertex {
private:
    OneRingHEdge ring;
public:
    OneRingVertex(const Vertex* v) : ring(v) { // constructor
    }

    Vertex* NextVertex() { // iterator
        HEdge* he = ring.NextHEdge();
        return (he) ? he->End() : NULL;
    }
};

////////// class Vertex //////////
class Vertex {
private:
    Eigen::Vector3d position; // position (x,y,z) in space
    Eigen::Vector3d normal; // normal vector for smooth shading rendering
    Eigen::Vector3d color; // color value for curvature displaying
    HEdge* he; // one of half edge starts with this vertex
    int index; // index in the Mesh::vList, DO NOT UPDATE IT
    int flag; // 0 for moved, 1 for constrained
    bool valid;

public:
    std::vector< HEdge* > adjHEdges; // for reading object only, do not use it in other place


    Vertex() : he(NULL), flag(0) { // constructors
        color = Eigen::Vector3d::Zero();
        normal = Eigen::Vector3d::Zero();
    }

    Vertex(const Eigen::Vector3d& v) : he(NULL), position(v), flag(0) {
        color = Eigen::Vector3d::Zero();
        normal = Eigen::Vector3d::Zero();
    }

    Vertex(double x, double y, double z) : he(NULL), position(x, y, z), flag(0) {
        color = Eigen::Vector3d::Zero();
        normal = Eigen::Vector3d::Zero();
    }

    // access functions
    const Eigen::Vector3d& Position() const {
        return position;
    }

    const Eigen::Vector3d& Normal() const {
        return normal;
    }

    const Eigen::Vector3d& Color() const {
        return color;
    }

    HEdge* HalfEdge() const {
        return he;
    }

    int Index() const {
        return index;
    }

    int Flag() const {
        return flag;
    }

    const Eigen::Vector3d& SetPosition(const Eigen::Vector3d& p) {
        return position = p;
    }

    const Eigen::Vector3d& SetNormal(const Eigen::Vector3d& n) {
        return normal = n;
    }

    const Eigen::Vector3d& SetColor(const Eigen::Vector3d& c) {
        return color = c;
    }

    HEdge* SetHalfEdge(HEdge* he) {
        return Vertex::he = he;
    }

    int SetIndex(int index) {
        return Vertex::index = index;
    }

    int SetFlag(int value) {
        return flag = value;
    }

    bool IsValid() const {
        return valid;
    }

    bool SetValid(bool b) {
        return valid = b;
    }

    bool IsBoundary() const {
        OneRingHEdge ring(this);
        HEdge* curr = NULL;
        while (curr = ring.NextHEdge())
            if (curr->IsBoundary())
                return true;
        return false;
    }

    int Valence() const {
        int count = 0;
        OneRingVertex ring(this);
        Vertex* curr = NULL;
        while (curr = ring.NextVertex()) count++;
        return count;
    }
};

////////// class Face //////////
class Face {
private:
    HEdge* he;
    bool valid;
public:
    // constructor
    Face() : he(NULL), valid(true) {
    }

    // access function
    HEdge* HalfEdge() const {
        return he;
    }

    HEdge* SetHalfEdge(HEdge* he) {
        return Face::he = he;
    }

    // check for boundary face
    bool IsBoundary() const {
        HEdge* curr = he;
        do {
            if (curr->Twin()->IsBoundary()) return true;
            curr = curr->Next();
        } while (curr != he);
        return false;
    }

    bool SetValid(bool b) {
        return valid = b;
    }

    bool IsValid() const {
        return valid;
    }
};


class Mesh {
public:
    HEdgeList heList;       // list of NON-boundary half edges
    HEdgeList bheList;      // list of boundary half egdes
    VertexList vList;       // list of vertices
    FaceList fList;         // list of faces



	FaceList myFaceList; //新加的facelist
	VertexList myvertexList;//新加的vertexlist
	VertexList bianjiedian;//初始的边界点
	VertexList normal_list;//初始边界点的外法线方向的点，
	VertexList all_list;


	double p0;
	double p1;
	double p2;
	double p3;
	Eigen::VectorXd lamada_value;

    // constructor & destructors
    Mesh() {
    }
    ~Mesh();

    // access functions
    const HEdgeList& Edges() const;
    const HEdgeList& BoundaryEdges() const;
    const VertexList& Vertices() const;
    const FaceList& Faces() const;


	//added
	const FaceList& holefzcelist() const;
	const VertexList& holevertexlist() const;
	const VertexList& bianjiedianlist() const;
	const VertexList& normai() const;



    // functions for loading obj files,
    // you DO NOT need to understand and use them
    bool LoadObjFile(const char* filename);
    void AddVertex(Vertex* v);
    void AddFace(int v1, int v2, int v3);
    void Clear();

    Eigen::Vector3d MinCoord() const;
    Eigen::Vector3d MaxCoord() const;

    /************************************************************************/
    /* please implement the following functions */
    void DisplayMeshInfo();
    void ComputeVertexNormals();
    void ComputeVertexCurvatures();
    void UmbrellaSmooth(bool uniformWeights = true);
    void ImplicitUmbrellaSmooth(bool uniformWeights = true);
    /************************************************************************/

    // additional helper functions
    int CountBoundaryLoops();
    int CountConnectedComponents();
	void TraverseComponents(HEdge*);

    void GroupingVertexFlags();

    static void SetPrevNext(HEdge* e1, HEdge* e2);
    static void SetTwin(HEdge* e1, HEdge* e2);
    static void SetFace(Face* f, HEdge* e);

    static double Cot(const Eigen::Vector3d& p1, const Eigen::Vector3d& p2, const Eigen::Vector3d& p3);
    static double TriArea(const Eigen::Vector3d& p1, const Eigen::Vector3d& p2, const Eigen::Vector3d& p3);

	 Eigen::Vector3d sandu(const Eigen::Vector3d& p11);//计算我们函数的散度vector3d
	 double RBF_network(const Eigen::Vector3d& p11);//得到原始函数

	void Mesh_fix();//步骤0，得到简单初始补洞网格
	void Least_squares_1();//步骤1得到最小二乘网格
	void Mesh_RBF_fix();//步骤2，使用RBF基函数得到fx，最后使用牛顿迭代法得到迭代点


};

#endif // __MESH_H__
