/*
 *
 * Demonstrates how to load and display an Wavefront OBJ file. 
 * Using triangles and normals as static object. No texture mapping.
 *
 * OBJ files must be triangulated!!!
 * Non triangulated objects wont work!
 * You can use Blender to triangulate
 *http://openglsamples.sourceforge.net/files/glut_obj.cpp
 */

#include "ObjModel.hpp"

// for mac osx
#ifdef __APPLE__
#include <OpenGL/gl.h>
#include <OpenGL/glu.h>
#include <GLUT/glut.h>
#else
// only for windows
#ifdef _WIN32
#include <windows.h>
#endif
// for windows and linux
#include <GL/gl.h>
#include <GL/glu.h>
#include <GL/freeglut.h>
#endif

#include <iostream>
#include <fstream>
#include <vector>
#include <stdio.h>
#include <string.h>
#include <assert.h>

 
#define KEY_ESCAPE 27

#define DELTA_ANGLE_X	5
#define DELTA_ANGLE_Y	5
#define DELTA_DISTANCE	0.3
#define DISTANCE_MIN	0.0
 
using namespace std;
 

//  Window
typedef struct {
    int width;
	int height;
	char* title;
 
	float field_of_view_angle;
	float z_near;
	float z_far;
} glutWindow;
 
 
  
//************************************
// global variable containing the OBJ model
//************************************
ObjModel obj;


int angle_y=0;
int angle_x=0;
float camDistance = 5;
bool flat = true;
bool wireframe = true;
glutWindow win;

// define a material in terms of its components
void define_material (	GLfloat ar, GLfloat ag, GLfloat ab, // ambient
						GLfloat dr, GLfloat dg, GLfloat db, // diffuse
						GLfloat sr, GLfloat sg, GLfloat sb, // specular
						GLfloat sh							// shininess
						) 
{
	GLfloat mat_ambient[4];
	GLfloat mat_diffuse[4];
	GLfloat mat_specular[4];
	
	mat_ambient[0] = ar;
	mat_ambient[1] = ag;
	mat_ambient[2] = ab;
	mat_ambient[3] = 1.0;
	glMaterialfv (GL_FRONT, GL_AMBIENT, mat_ambient);
	
	mat_diffuse[0] = dr;
	mat_diffuse[1] = dg;
	mat_diffuse[2] = db;
	mat_diffuse[3] = 1.0;
	glMaterialfv (GL_FRONT, GL_DIFFUSE, mat_diffuse);	

	mat_specular[0] = sr;
	mat_specular[1] = sg;
	mat_specular[2] = sb;
	mat_specular[3] = 1.0;
	glMaterialfv (GL_FRONT, GL_SPECULAR, mat_specular);

	glMaterialf (GL_FRONT, GL_SHININESS, sh);	
}


void DrawAxis(float scale)
{
    glPushMatrix();
    glDisable(GL_LIGHTING);

    glScalef(scale, scale, scale);

    glBegin(GL_LINES);

		glColor3f(1.0, 0.0, 0.0);
		/*  X axis */
		glVertex3f(0.0, 0.0, 0.0);
		glVertex3f(1.0, 0.0, 0.0);	
		/*  Letter X */
		glVertex3f(.8f, 0.05f, 0.0); // "backslash""
		glVertex3f(1.0, 0.25f, 0.0);	
		glVertex3f(0.8f, .25f, 0.0); // "slash"
		glVertex3f(1.0, 0.05f, 0.0);

		/*  Y axis */
		glColor3f(0.0, 1.0, 0.0);
		glVertex3f(0.0, 0.0, 0.0);
		glVertex3f(0.0, 1.0, 0.0);	

		// z-axis
		glVertex3f(0.0, 0.0, 0.0);
		glVertex3f(0.0, 0.0, 1.0);	
		/*  Letter Z */
		glColor3f(0.0, 0.0, 1.0);
		glVertex3f(0.0, 0.05f, 0.8);  // bottom horizontal leg
		glVertex3f(0.0, 0.05f, 1.0);	
		glVertex3f(0.0, 0.05f, 1.0);   // slash
		glVertex3f(0.0, 0.25f, 0.8);
		glVertex3f(0.0, 0.25f, 0.8);  // upper horizontal leg
		glVertex3f(0.0, 0.25f, 1.0);	    
	
    glEnd();

	glEnable(GL_LIGHTING);

    glColor3f(1.0, 1.0, 1.0);
    glPopMatrix();
}

// place the light in x,y,z
void place_light (GLfloat x, GLfloat y, GLfloat z) 
{

	GLfloat light_position[4];
	GLfloat light_ambient[] = {0.2, 0.2, 0.2, 1.0};
	GLfloat light_diffuse[] = {1.0, 1.0, 1.0, 1.0};	
	GLfloat light_specular[] = {1.0, 1.0, 1.0, 1.0};


	light_position[0] = x;
	light_position[1] = y;
	light_position[2] = z;
	light_position[3] = 1.0;

	glLightfv (GL_LIGHT0, GL_AMBIENT, light_ambient);
	glLightfv (GL_LIGHT0, GL_DIFFUSE, light_diffuse);
	glLightfv (GL_LIGHT0, GL_SPECULAR, light_specular);
	glLightfv (GL_LIGHT0, GL_POSITION, light_position);

	glEnable (GL_LIGHT0);
	glEnable (GL_LIGHTING);

}


void display() 
{
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();

	glPushMatrix();
		place_light (0, 0, 4);
		glTranslatef (0, 0, -camDistance);
		glRotatef (angle_x, 1, 0, 0);
		glRotatef (angle_y, 0, 1, 0);

//		angle_y += 1;

		DrawAxis(1.0f);
		
		define_material (
			0.2, 0.2, 0.2,
			0.8, 0.8, 0.8,
			1.0, 0.8, 0.8,
			250	);
		//***********************************************
		// draw the model
		//***********************************************
        if(flat)
        {
//            cout << "flat";
            obj.draw();
//            obj.flatDraw();
        }
        else
        {
//            cout << "noflat";
            obj.indexDraw();
        }
        if(wireframe)
            obj.wireframetDraw();
		
	glPopMatrix();
	
	glutSwapBuffers();
}
 

// initialize the opengl
void initialize () 
{
    glMatrixMode(GL_MODELVIEW);
    glShadeModel( GL_SMOOTH );
    
	glEnable (GL_CULL_FACE);
	
    glEnable( GL_DEPTH_TEST );
    glDepthFunc( GL_LEQUAL );
    glHint( GL_PERSPECTIVE_CORRECTION_HINT, GL_NICEST );

	glEnable (GL_NORMALIZE);
	
    glMatrixMode(GL_PROJECTION);
	glViewport(0, 0, win.width, win.height);
	GLfloat aspect = (GLfloat) win.width / win.height;
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
	gluPerspective(win.field_of_view_angle, aspect, win.z_near, win.z_far);
    glMatrixMode(GL_MODELVIEW);
}



 
void keyboard ( unsigned char key, int x, int y ) 
{
  switch ( key ) {
    case KEY_ESCAPE:        
      exit ( 0 );   
      break;
  case 's':
      PRINTVAR(flat);
      flat = !flat;
    break;
  case 'w':
      PRINTVAR(wireframe);
      wireframe = !wireframe;
    break;
    default:      
      break;
  }
  glutPostRedisplay ();
}

void arrows (int key, int x, int y) {
	switch (key) 
	{
		case GLUT_KEY_UP:
			angle_x = (angle_x + DELTA_ANGLE_X) % 360;
		break;	
		case GLUT_KEY_DOWN:
			angle_x = (angle_x - DELTA_ANGLE_X) % 360;
		break;	
		case GLUT_KEY_LEFT:
			angle_y = (angle_y + DELTA_ANGLE_Y) % 360;
		break;	
		case GLUT_KEY_RIGHT:
			angle_y = (angle_y - DELTA_ANGLE_Y) % 360;
		break;
		case GLUT_KEY_PAGE_DOWN:
			camDistance += DELTA_DISTANCE;
		break;	
		case GLUT_KEY_PAGE_UP:
			camDistance -= (camDistance>DISTANCE_MIN)? DELTA_DISTANCE: 0.0;
		break;
		default: break;
	}
	glutPostRedisplay ();
}
 

int main(int argc, char **argv) 
{
	// set window values
	win.width = 640;
	win.height = 480;
	win.title = "OpenGL/GLUT OBJ Loader.";
	win.field_of_view_angle = 45;
    win.z_near = 0.25f;
	win.z_far = 500.0f;

//    vector<point3d> myvec;

//    myvec.push_back(point3d(1.f,2.f,3.f));
//    myvec.push_back(point3d(4.f,5.f,6.f));
//    myvec.push_back(point3d(7.f,8.f,9.f));

//    float* a = (float*)&myvec[0];

//    for(int i=0; i<9; ++i)
//    {
//        cout << a[i] << endl;
//    }



 
	// initialize and run program
	glutInit(&argc, argv);                                      
	glutInitDisplayMode(GLUT_RGB | GLUT_DOUBLE | GLUT_DEPTH ); 
	glutInitWindowSize(win.width,win.height);					
	glutCreateWindow(win.title);								
	glutDisplayFunc(display);	
//	glutIdleFunc( display );									// register Idle Function
    glutKeyboardFunc( keyboard );								
    glutSpecialFunc( arrows );
	initialize();

	//***********************************************
	// Load the obj model from file
	//***********************************************
	obj.load(argv[1]);
	
	//***********************************************
	// Make it unitary
	//***********************************************
	obj.unitizeModel();
	
	glutMainLoop();						
	
	
	return EXIT_SUCCESS;
}
