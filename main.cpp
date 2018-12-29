#include <iostream>
#include "Mesh.h"
#include "GLProjector.h"
#include "Deformer.h"

// the way the mesh is rendered
enum EnumDisplayMode {
    WIREFRAME,
    HIDDENLINE,
    FLATSHADED,
    SMOOTHSHADED,
    COLORSMOOTHSHADED,
	FIXHOLE
};

// variables
int displayMode = SMOOTHSHADED;                     // current display mode
int mainMenu, displayMenu;                          // glut menu handlers
int winWidth, winHeight;                            // window width and height
double winAspect;                                   // winWidth / winHeight;
int lastX, lastY;                                   // last mouse motion position
bool leftDown, middleDown, shiftDown;               // mouse down and shift down flags
double sphi = 90.0, stheta = 45.0, sdepth = 10;     // for simple trackball
double xpan = 0.0, ypan = 0.0;                      // for simple trackball
double zNear = 1.0, zFar = 100.0;                   // clipping
double g_fov = 45.0;
Eigen::Vector3d g_center;
double g_sdepth;
Mesh mesh; // our mesh

// editing mode
enum Mode {
    Viewing,
    Selection,
    Moving
};

Mode currentMode = Viewing;
int downX, downY; // mouse down position
int selectedHandleIndex = -1; // the index of the handle
Deformer* deformer = NULL; // Deformer, make a natural deformation

// functions
void SetBoundaryBox(const Eigen::Vector3d& bmin, const Eigen::Vector3d& bmax);
void InitGL();
void InitMenu();
void InitGeometry();

// window related 
void MenuCallback(int value);
void ReshapeFunc(int width, int height);

// rendering functions
void DisplayFunc();
void DrawWireframe();
void DrawHiddenLine();
void DrawFlatShaded();
void DrawSmoothShaded();
void DrawColorSmoothShaded();
void DrawSelectionRect();//draw white rectangle
void DrawSelectedVertices();
void Fix_Hole();

// input related glut functions
void KeyboardFunc(unsigned char ch, int x, int y);
void MouseFunc(int button, int state, int x, int y);
void MotionFunc(int x, int y);

void SelectVertexByRect(); // select ROI for editing
int StartMove(); // a single editing operation


void SetBoundaryBox(const Eigen::Vector3d& bmin, const Eigen::Vector3d& bmax) {
    double PI = 3.14159265358979323846;
    double radius = (bmax - bmin).norm();
    g_center = 0.5 * (bmin + bmax);
    zNear = 0.2 * radius / sin(0.5 * g_fov * PI / 180.0);
    zFar = zNear + 2.0 * radius;
    g_sdepth = zNear + radius;
    zNear *= 0.1;
    zFar *= 10;
    sdepth = g_sdepth;
}

// init openGL environment
void InitGL() {
    GLfloat light0Position[] = {0, 1, 0, 1.0};

    // initialize GLUT stuffs
    glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH);
    glutInitWindowSize(500, 500);
    glutCreateWindow("Comp5411 Mesh Viewer");

    glClearColor(0.0, 0.0, 0.0, 0.0);
    glPolygonOffset(1.0, 1.0);
    glDepthFunc(GL_LEQUAL);
    glEnable(GL_DEPTH_TEST);
    glEnable(GL_CULL_FACE);
    glEnable(GL_COLOR_MATERIAL);
    glColorMaterial(GL_FRONT, GL_DIFFUSE);
    glLightfv(GL_LIGHT0, GL_POSITION, light0Position);
    glEnable(GL_LIGHT0);

    glutReshapeFunc(ReshapeFunc);
    glutDisplayFunc(DisplayFunc);
    glutKeyboardFunc(KeyboardFunc);
    glutMouseFunc(MouseFunc);
    glutMotionFunc(MotionFunc);
}

// init right-click menu
void InitMenu() {
    displayMenu = glutCreateMenu(MenuCallback);
    glutAddMenuEntry("Wireframe", WIREFRAME);
    glutAddMenuEntry("Hidden Line", HIDDENLINE);
    glutAddMenuEntry("Flat Shaded", FLATSHADED);
    glutAddMenuEntry("Smooth Shaded", SMOOTHSHADED);
    glutAddMenuEntry("Color Smooth Shaded", COLORSMOOTHSHADED);
	glutAddMenuEntry("Fix the hole", 111);

    mainMenu = glutCreateMenu(MenuCallback);
    glutAddSubMenu("Display", displayMenu);
    glutAddMenuEntry("Build Deformer", 101);
    glutAddMenuEntry("Exit", 99);
    glutAttachMenu(GLUT_RIGHT_BUTTON);
}

// init geometry (if no input argument is provided)
void InitGeometry() {
    const int VSIZE = 4;
    const int HESIZE = 12;
    const int FSIZE = 4;
    int i;
    Vertex* v[VSIZE];
    HEdge* he[HESIZE];
    Face* f[FSIZE];

    for (i = 0; i < VSIZE; i++) {
        v[i] = new Vertex();
        mesh.vList.push_back(v[i]);
    }
    v[0]->SetPosition(Eigen::Vector3d(0.0, 0.0, 0.0));
    v[1]->SetPosition(Eigen::Vector3d(10.0, 0.0, 0.0));
    v[2]->SetPosition(Eigen::Vector3d(0.0, 10.0, 0.0));
    v[3]->SetPosition(Eigen::Vector3d(0.0, 0.0, 10.0));

    v[0]->SetNormal(Eigen::Vector3d(-0.577, -0.577, -0.577));
    v[1]->SetNormal(Eigen::Vector3d(0.0, -0.7, -0.7));
    v[2]->SetNormal(Eigen::Vector3d(-0.7, 0.0, -0.7));
    v[3]->SetNormal(Eigen::Vector3d(-0.7, -0.7, 0.0));

    for (i = 0; i < FSIZE; i++) {
        f[i] = new Face();
        mesh.fList.push_back(f[i]);
    }

    for (i = 0; i < HESIZE; i++) {
        he[i] = new HEdge();
        mesh.heList.push_back(he[i]);
    }
    for (i = 0; i < FSIZE; i++) {
        int base = i * 3;
        Mesh::SetPrevNext(he[base], he[base + 1]);
        Mesh::SetPrevNext(he[base + 1], he[base + 2]);
        Mesh::SetPrevNext(he[base + 2], he[base]);
        Mesh::SetFace(f[i], he[base]);
    }
    Mesh::SetTwin(he[0], he[4]);
    Mesh::SetTwin(he[1], he[7]);
    Mesh::SetTwin(he[2], he[10]);
    Mesh::SetTwin(he[3], he[8]);
    Mesh::SetTwin(he[5], he[9]);
    Mesh::SetTwin(he[6], he[11]);
    he[0]->SetStart(v[1]);
    he[1]->SetStart(v[2]);
    he[2]->SetStart(v[3]);
    he[3]->SetStart(v[0]);
    he[4]->SetStart(v[2]);
    he[5]->SetStart(v[1]);
    he[6]->SetStart(v[0]);
    he[7]->SetStart(v[3]);
    he[8]->SetStart(v[2]);
    he[9]->SetStart(v[0]);
    he[10]->SetStart(v[1]);
    he[11]->SetStart(v[3]);
    v[0]->SetHalfEdge(he[3]);
    v[1]->SetHalfEdge(he[0]);
    v[2]->SetHalfEdge(he[1]);
    v[3]->SetHalfEdge(he[2]);
}


// GLUT menu callback function
void MenuCallback(int value) {
    switch (value) {
    case 99:
        exit(0);
        break;
    case 101:
        deformer = new Deformer(&mesh);
        std::cout << "Deformer Building Finished\n";
        break;
	case 111:
		Fix_Hole();
		glutPostRedisplay();
		break;
    default:
        displayMode = value;
        glutPostRedisplay();
        break;
    }
}


// GLUT reshape callback function
void ReshapeFunc(int width, int height) {
    winWidth = width;
    winHeight = height;
    winAspect = (double)width / (double)height;
    glViewport(0, 0, width, height);
    glutPostRedisplay();
}


// GLUT display callback function
void DisplayFunc() {
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    gluPerspective(g_fov, winAspect, zNear, zFar);
    //glOrtho(-2.0, 2.0, -2.0, 2.0, zNear, zFar);

    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
    glTranslatef(xpan, ypan, -sdepth);
    glRotatef(-stheta, 1.0, 0.0, 0.0);
    glRotatef(sphi, 0.0, 1.0, 0.0);
    glTranslatef(-g_center[0], -g_center[1], -g_center[2]);

    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    switch (displayMode) {
    case WIREFRAME:
        DrawWireframe();
        break;
    case HIDDENLINE:
        DrawHiddenLine();
        break;
    case FLATSHADED:
        DrawFlatShaded();
        break;
    case SMOOTHSHADED:
        DrawSmoothShaded();
        break;
    case COLORSMOOTHSHADED:
        DrawColorSmoothShaded();
        break;

    }

    DrawSelectedVertices();

    if (currentMode == Selection && downX != lastX && downX != lastY) {
        DrawSelectionRect();// judege if we didn't move it by mouse.
    }

    glutSwapBuffers();
}

// Wireframe render function
void DrawWireframe() {
    HEdgeList heList = mesh.Edges();
    HEdgeList bheList = mesh.BoundaryEdges();

	FaceList myholelist = mesh.holefzcelist();
	VertexList my_Vertex_list = mesh.holevertexlist();
	VertexList bianver = mesh.bianjiedianlist();
	VertexList my_normal_ver = mesh.normai();

    glColor3f(1.0f, 1.0f, 1.0f);
    glBegin(GL_LINES);
    size_t i;
    for (i = 0; i < heList.size(); i++) {
        glVertex3dv(heList[i]->Start()->Position().data());
        glVertex3dv(heList[i]->End()->Position().data());
    }
    glColor3f(1.0f, 0.0f, 0.0f);
    for (i = 0; i < bheList.size(); i++) {
        glVertex3dv(bheList[i]->Start()->Position().data());
        glVertex3dv(bheList[i]->End()->Position().data());
    }
    glEnd();



	glColor3f(0, 0, 1);
	glBegin(GL_TRIANGLES);
	for (int i = 0; i < myholelist.size(); i++)
	{
		Face* f = myholelist[i];
		const Eigen::Vector3d& pos1 = f->HalfEdge()->Start()->Position();
		const Eigen::Vector3d& pos2 = f->HalfEdge()->End()->Position();
		const Eigen::Vector3d& pos3 = f->HalfEdge()->Next()->End()->Position();
		glVertex3dv(pos1.data());
		glVertex3dv(pos2.data());
		glVertex3dv(pos3.data());
	}
	glEnd();

	//glColor3f(1,0,0);
	//glBegin(GL_LINES);
	//for (int m = 0; m < my_normal_ver.size(); m++)
	//{
	//	glVertex3dv(bianver[m]->Position().data());
	//	glVertex3dv(my_normal_ver[m]->Position().data());

	//}
	//glEnd();

}

// Hidden Line render function
void DrawHiddenLine() {
    FaceList fList = mesh.Faces();
    glShadeModel(GL_FLAT);
    glEnable(GL_POLYGON_OFFSET_FILL);
    glColor3f(0, 0, 0);
    glBegin(GL_TRIANGLES);
    for (size_t i = 0; i < fList.size(); i++) {
        Face* f = fList[i];
        const Eigen::Vector3d& pos1 = f->HalfEdge()->Start()->Position();
        const Eigen::Vector3d& pos2 = f->HalfEdge()->End()->Position();
        const Eigen::Vector3d& pos3 = f->HalfEdge()->Next()->End()->Position();
        glVertex3dv(pos1.data());
        glVertex3dv(pos2.data());
        glVertex3dv(pos3.data());
    }
    glEnd();

	

    glDisable(GL_POLYGON_OFFSET_FILL);

    DrawWireframe();
}

// Flat Shaded render function
void DrawFlatShaded() {
    FaceList fList = mesh.Faces();
    glShadeModel(GL_FLAT);
    glEnable(GL_LIGHTING);
    glColor3f(0.4f, 0.4f, 1.0f);
    glBegin(GL_TRIANGLES);
    for (size_t i = 0; i < fList.size(); i++) {
        Face* f = fList[i];
        const Eigen::Vector3d& pos1 = f->HalfEdge()->Start()->Position();
        const Eigen::Vector3d& pos2 = f->HalfEdge()->End()->Position();
        const Eigen::Vector3d& pos3 = f->HalfEdge()->Next()->End()->Position();
        Eigen::Vector3d normal = (pos2 - pos1).cross(pos3 - pos1);
        normal.normalize();
        glNormal3dv(normal.data());
        glVertex3dv(pos1.data());
        glVertex3dv(pos2.data());
        glVertex3dv(pos3.data());
    }
    glEnd();
    glDisable(GL_LIGHTING);
}

// Smooth Shaded render function
void DrawSmoothShaded() {
    FaceList fList = mesh.Faces();
    glShadeModel(GL_SMOOTH);
    glEnable(GL_LIGHTING);
    glColor3f(0.4f, 0.4f, 1.0f);
    glBegin(GL_TRIANGLES);
    for (size_t i = 0; i < fList.size(); i++) {
        Face* f = fList[i];
        Vertex* v1 = f->HalfEdge()->Start();
        Vertex* v2 = f->HalfEdge()->End();
        Vertex* v3 = f->HalfEdge()->Next()->End();

        glNormal3dv(v1->Normal().data());
        glVertex3dv(v1->Position().data());
        glNormal3dv(v2->Normal().data());
        glVertex3dv(v2->Position().data());
        glNormal3dv(v3->Normal().data());
        glVertex3dv(v3->Position().data());
    }
    glEnd();
    glDisable(GL_LIGHTING);
}

void DrawColorSmoothShaded() {
    FaceList fList = mesh.Faces();
    glShadeModel(GL_SMOOTH);
    glEnable(GL_LIGHTING);
    glColor3f(0.4f, 0.4f, 1.0f);
    glBegin(GL_TRIANGLES);
    for (size_t i = 0; i < fList.size(); i++) {
        Face* f = fList[i];
        Vertex* v1 = f->HalfEdge()->Start();
        Vertex* v2 = f->HalfEdge()->End();
        Vertex* v3 = f->HalfEdge()->Next()->End();
        glNormal3dv(v1->Normal().data());
        glColor3dv(v1->Color().data());
        glVertex3dv(v1->Position().data());
        glNormal3dv(v2->Normal().data());
        glColor3dv(v2->Color().data());
        glVertex3dv(v2->Position().data());
        glNormal3dv(v3->Normal().data());
        glColor3dv(v3->Color().data());
        glVertex3dv(v3->Position().data());
    }
    glEnd();
    glDisable(GL_LIGHTING);
}

void DrawSelectionRect() {
    glMatrixMode(GL_PROJECTION);
    glPushMatrix();
    glLoadIdentity();
    gluOrtho2D(0, winWidth, 0, winHeight);
    glMatrixMode(GL_MODELVIEW);
    glPushMatrix();
    glLoadIdentity();

    glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
    glDisable(GL_CULL_FACE);
    glDisable(GL_DEPTH_TEST);
    glColor3f(1.0f, 1.0f, 1.0f);
    glRectd(downX, winHeight - downY, lastX, winHeight - lastY);
    glEnable(GL_CULL_FACE);
    glEnable(GL_DEPTH_TEST);
    glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);

    glMatrixMode(GL_PROJECTION);
    glPopMatrix();
    glMatrixMode(GL_MODELVIEW);
    glPopMatrix();
}


// draw the selected ROI vertices on the mesh
void DrawSelectedVertices() {
    VertexList vList = mesh.Vertices();
    glColor3f(1.0, 0.0, 0.0);
    glPointSize(2.0);
    glBegin(GL_POINTS);
    size_t i;
    for (i = 0; i < vList.size(); i++) {
        if (vList[i]->Flag()) {
            switch (vList[i]->Flag() % 3) {
            case 0: // handle vertices
                glColor3f(1.0, 0.3, 0.3);//red
                break;
            case 1:
                glColor3f(0.3, 1.0, 0.3);//green
                break;
            case 2:
                glColor3f(0.3, 0.3, 1.0);//blue
                break;
            }
            glVertex3dv(vList[i]->Position().data());
        }
    }
    glEnd();
}


// GLUT keyboard callback function
void KeyboardFunc(unsigned char ch, int x, int y) {
	
	int IsUniform; 
	int IterNum;
	switch (ch) {
    case 'u':
    case 'U':
        /************************************************************************/
        /* activate the following code if you finish the corresponding functions*/
        
        /*====== Programming Assignment 1 ======*/
		std::cout << "Choose mode for explicit Umbrellasmooth (1 for Uniform, 2 for cotangent): ";
		std::cin >> IsUniform;
		std::cout << "input the number of iteration: ";
		std::cin >> IterNum;
		if (IsUniform == 1){
			for (int i = 0; i < IterNum; i++)
				mesh.UmbrellaSmooth(true);
		}
		else if (IsUniform == 2){
			for (int i = 0; i < IterNum; i++)
				mesh.UmbrellaSmooth(false);
		}

        mesh.ComputeVertexNormals();
        mesh.ComputeVertexCurvatures();
        /************************************************************************/
        break;

    case 's':
    case 'S':
        /*====== Programming Assignment 1 ======*/
		std::cout << "Choose mode for implicit Umbrellasmooth (1 for Uniform, 2 for contangent):";
		std::cin >> IsUniform;
		std::cout << "input the number of iteration: ";
		std::cin >> IterNum;
		if (IsUniform == 1){
			for (int i = 0; i < IterNum; i++)
				mesh.ImplicitUmbrellaSmooth(true);
		}
		else if (IsUniform == 2){
			for (int i = 0; i < IterNum; i++)
				mesh.ImplicitUmbrellaSmooth(false);
		}
        mesh.ComputeVertexNormals();
        mesh.ComputeVertexCurvatures();
        break;

    case '1': // key '1'
        currentMode = Viewing;
        std::cout << "Viewing mode" << std::endl;
        break;
    case '2': // key '2'
        currentMode = Selection;
        std::cout << "Selection mode" << std::endl;
        break;
    case '3': // key '3'
        currentMode = Moving;
        std::cout << "Moving mode" << std::endl;
        break;
    case 27: // Esc char
        exit(0);
        break;
    }
    glutPostRedisplay();
}


// GLUT mouse callback function
void MouseFunc(int button, int state, int x, int y) {
    if (state == GLUT_DOWN) {
        downX = x;
        downY = y;
    }
    lastX = x;
    lastY = y;
    leftDown = (button == GLUT_LEFT_BUTTON) && (state == GLUT_DOWN);
    middleDown = (button == GLUT_MIDDLE_BUTTON) && (state == GLUT_DOWN);
    shiftDown = (glutGetModifiers() & GLUT_ACTIVE_SHIFT);

    if (currentMode == Selection && state == GLUT_UP) {
        SelectVertexByRect();
        downX = downY = lastX = lastY = 0;
        glutPostRedisplay();
        mesh.GroupingVertexFlags();
    }

    if (currentMode == Moving && state == GLUT_DOWN) {
        selectedHandleIndex = StartMove();
        std::cout << "handle " << selectedHandleIndex << " selected\n";
    }

    if (currentMode == Moving && state == GLUT_UP) {
        if (deformer) {
            // use the deformer to get a natural deformation
            deformer->Deform();
            mesh.ComputeVertexNormals();
            mesh.ComputeVertexCurvatures();
            glutPostRedisplay();
        }
    }

}


// GLUT mouse motion callback function
void MotionFunc(int x, int y) {
    switch (currentMode) {
    case Viewing:
        if (leftDown) {
            if (!shiftDown) { // rotate
                sphi += (double)(x - lastX) / 4.0;
                stheta += (double)(lastY - y) / 4.0;
            } else { // pan with shift key
                xpan += (double)(x - lastX) * sdepth / zNear / winWidth;
                ypan += (double)(lastY - y) * sdepth / zNear / winHeight;
            }
        }
        // scale
        if (middleDown)
            sdepth += (double)(lastY - y) / 10.0;
        break;

    case Selection:
        break;

    case Moving:
        GLProjector projector;
        VertexList vList = mesh.Vertices();
        int diffX = x - lastX;
        int diffY = y - lastY;
        if (diffX == 0 && diffY == 0)
            break;
        Eigen::Vector3d offset1 = projector.UnProject(x, -y, 0);
        Eigen::Vector3d offset2 = projector.UnProject(lastX, -lastY, 0);
        Eigen::Vector3d offset = (offset1 - offset2) * 20;
        for (size_t i = 0; i < vList.size(); i++) {
            if (vList[i]->Flag() == selectedHandleIndex)
                vList[i]->SetPosition(vList[i]->Position() + offset);
        }
//         if (deformer) {
//             deformer->Deform();
//             mesh.ComputeVertexNormals();
//             mesh.ComputeVertexCurvatures();
//         }
        break;
    }

    lastX = x;
    lastY = y;
    glutPostRedisplay();
}

// select the ROI
void SelectVertexByRect() {
    // get the selection rectangle
    int x1 = downX, x2 = lastX, y1 = winHeight - downY, y2 = winHeight - lastY;
    if (x2 < x1) {
        int tmp = x1;
        x1 = x2;
        x2 = tmp;
    }
    if (y2 < y1) {
        int tmp = y1;
        y1 = y2;
        y2 = tmp;
    }
    GLProjector projector;

    VertexList vList = mesh.Vertices();
    if (shiftDown) // the case in which we wanna add more points to the mesh
    {
        for (size_t i = 0; i < vList.size(); i++) {
            Eigen::Vector3d v = projector.Project(vList[i]->Position());
            if ((v(0) >= 0 && v(0) < winWidth && v(1) >= 0 && v(1) < winHeight) &&
                (v(0) > x1 && v(0) < x2 && v(1) > y1 && v(1) < y2))
                vList[i]->SetFlag(1);
        }
    } else {
        for (size_t i = 0; i < vList.size(); i++) {
            Eigen::Vector3d v = projector.Project(vList[i]->Position());
            if ((v(0) >= 0 && v(0) < winWidth && v(1) >= 0 && v(1) < winHeight) &&
                (v(0) > x1 && v(0) < x2 && v(1) > y1 && v(1) < y2))
                vList[i]->SetFlag(1);
            else
                vList[i]->SetFlag(0);
        }
    }
}


// perform a single editing operation
int StartMove() {
	GLProjector projector;
	VertexList vList = mesh.Vertices();

	double minDis = 0;
	int minIndex = -1;

	// find the closest vertex to the mouse position ?
	for (size_t i = 0; i < vList.size(); i++) {
		if (vList[i]->Flag() == 0)
			continue;

		Eigen::Vector3d v = projector.Project(vList[i]->Position());
		double diffX = v(0) - lastX;
		double diffY = v(1) - (winHeight - lastY);
		double dis = diffX * diffX + diffY * diffY;
		if (dis < minDis || minIndex == -1) {
			minDis = dis;
			minIndex = (int)i;
		}
	}

	return (minIndex != -1) ? vList[minIndex]->Flag() : 0;
}


void Fix_Hole()
{
	displayMode = SMOOTHSHADED;

	mesh.Mesh_RBF_fix();

}



// main function
void main(int argc, char** argv) {
    glutInit(&argc, argv);
    InitGL();
    InitMenu();
	std::cout << "°´ÏÂÓÒ¼üÑ¡Ôñfixhole" << std::endl;
	if (argc >= 2)
	{
		mesh.LoadObjFile(argv[1]);
		//mesh.Least_squares_1();
		//mesh.Mesh_RBF_fix();
	}
	else
	{
		InitGeometry();
	}

    SetBoundaryBox(mesh.MinCoord(), mesh.MaxCoord());
    mesh.ComputeVertexNormals();
    mesh.ComputeVertexCurvatures();

    /************************************************************************/
    /* activate the following code if you finish the corresponding functions*/
    mesh.DisplayMeshInfo();
    /************************************************************************/
    
	

    glutMainLoop();
}
