
#ifndef inc_muhaha_h
#define inc_muhaha_h
//}

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <sys/file.h>
#include <string.h>
#include <pthread.h>
#include <SDL/SDL.h>
#include <SDL/SDL_thread.h>
#include <SDL/SDL_audio.h>
#include <SDL/SDL_timer.h>
#include <linux/videodev2.h>
#include <sys/ioctl.h>
#include <sys/mman.h>
#include <errno.h>
#include <fcntl.h>
#include <time.h>
#include <sys/time.h>
#include <signal.h>
#include <X11/Xlib.h>
#include <SDL/SDL_syswm.h>
#include "v4l2uvc.h"
#include "gui.h"
#include "utils.h"
#include "color.h"

#include <math.h>


typedef unsigned char u08;
typedef unsigned short u16;
typedef unsigned int u32;


typedef unsigned long u00;
typedef signed long s00;

typedef u00 ui;
typedef s00 si;

typedef struct {
	si x, y;
}tV2si;

typedef struct {
	float x, y;
}tV2f;

typedef struct {
	float x, y, z;
}tV3f;

typedef struct {
	float x, y, z, w;
}tV4f;

typedef struct {
	double x, y;
}tV2d;


typedef struct {
	float a, b, c, d;
}tM2f;

typedef struct {
	union {
	struct {
		float x00, x01, x02;
		float x10, x11, x12;
		float x20, x21, x22;
	};
	float f[3][3];
	};
}tM3f;

typedef struct {
	union {
	struct {
		float x00, x01, x02, x03;
		float x10, x11, x12, x13;
		float x20, x21, x22, x23;
		float x30, x31, x32, x33;
	};
	float f[4][4];
	};
}tM4f;

#define CALIBRATIONPOINTS    9
typedef struct {
	//tV2f  calipoints[CALIBRATIONPOINTS];       //conversion from eye to scene calibration points
	tV2f  scenecalipoints[CALIBRATIONPOINTS];  //captured (with mouse) calibration points
	//tV2f  pucalipoints[CALIBRATIONPOINTS];     //captured eye points while looking at the calibration points in the scene
	//tV2f  crcalipoints[CALIBRATIONPOINTS];     //captured corneal reflection points while looking at the calibration points in the scene
	tV2f  vectors[CALIBRATIONPOINTS];          //differences between the corneal reflection and pupil center
	
	double map_matrix[3][3];
	
	float aa, bb, cc, dd, ee;                       //pupil X coefficients
	float ff, gg, hh, ii, jj;			//pupil Y coefficients
	
	float centx, centy;                             // translation to center pupil data after biquadratics
	int inx, iny;                                   // translation to center pupil data before biquadratics
	float cmx[4], cmy[4];                           // corner correctioncoefficients
}tHomo;

enum {
	eEye_Fit_Default,
	eEye_Fit_SFit,
	eEye_Fit_S0,
	eEye_Fit_S2Fit,
	eEye_Fit_S3Fit_START,
	eEye_Fit_S3Fit_Point,
	eEye_Fit_S3Fit_Eye,
	eEye_Fit_S3Fit_Eye3,
	eEye_Fit_S3Fit_END,
	eEye_Fit_S4Fit_START,
	eEye_Fit_S4Fit_Edge0,
	eEye_Fit_S4Fit_Edge1,
	eEye_Fit_S4Fit_Edge2,
	eEye_Fit_S4Fit_END,
};

typedef struct {
	tV2f P; tV2f V; // Position, Velocity
	float Ax, Ay, Aa;
	
	si Exp_R;
	
	u08 Fit, Fit_Edge;
	float Fit_Scale, Fit_Trans;
	float S2Fit_Scale, S2Fit_Trans;
	
	si Point_N, Point_Max;
	tV2f*	paPoint;
	
	tV2si LinView;
	
	tHomo Homo;
}tEye;


typedef struct {
	tEye DotC, DotL, DotR;
	tM3f M, MI;
	
	tM4f M4, M4I;
}tHead;


typedef struct {
	ui Full_W, Full_H;
	float Full_FOV;
	
	ui Image_W, Image_H;
	float Image_FOV;
	
	ui Zoom;
}tCam;

typedef struct {
	int Y_Level;
	
	u08 DeSat;
	
	tCam Cam;
	
	si View_W, View_H;
	float Proj_L, Proj_R, Proj_B, Proj_T;
	float Proj_W, Proj_H, Proj_N, Proj_F;
	tM4f Proj, World;
	
	
	tHead Head, HeadC;
//		tEye DotC, DotL, DotR;
//	}Cent;
//	tM3f HeadC, Head, HeadT;
	
	tEye Left, Right;
	tV2f L_Vec, R_Vec;
	
	SDL_mutex* pGaze_Mutex;
	tV2f GazeL, GazeR, Gaze;
	
	u08* pDst;
	
	SDL_Surface* pScreen;
	SDL_Overlay* pOverlay;
	
	
	struct {
		float x, y, z;
	}tmp;
}tM;

extern tM gM;

extern struct vdIn *videoIn;

void	muhaha_Init	();
void	muhaha	();


struct pt_data {
	SDL_Surface **ptscreen;
	SDL_Event *ptsdlevent;
	SDL_Rect *drect;
	struct vdIn *ptvideoIn;
	float frmrate;
	SDL_mutex *affmutex;
} ptdata;

int muhaha_eventThread(void *data);



typedef struct _dyn_config_entry dyn_config_entry;
typedef struct _dyn_config_mapping dyn_config_mapping;
typedef struct _dyn_config dyn_config;
typedef enum DYN_ENTRY_TYPE dyn_entry_t;

enum { dyn_config_entry_value_max = 64 };
struct _dyn_config_entry {
	dyn_config_entry *next;
	dyn_config_entry *child;
	char name[32];
	char value[dyn_config_entry_value_max];
//	float value;
};

struct _dyn_config {
	dyn_config_entry *entries;
	int count;
};


void dyn_config_read(dyn_config *dc, const char *f_name);
void dyn_config_watch(dyn_config *dc, const char *f_name);




void svd(int m, int n, double **a, double **p, double *d, double **q);
tV2d* normalize_point_set(tV2d* point_set, double *dis_scale, tV2d *nor_center, int num);

extern ui gM_edge_point_N;
extern tV2d gM_edge_point[204];
extern double pupil_param[5];//parameters of an ellipse {ellipse_a, ellipse_b, cx, cy, theta}; a & b is the major or minor axis; 







typedef struct TObject{
	int nbPts;
	double **objectPts;/* Array of coordinates of points */
	double **objectVects;/* Array of coordinates of vectors from reference pt to all pts */
	double **objectCopy;/* Copy of objectVects, used because SVD code destroys input data */
	double **objectMatrix;/* Pseudoinverse of objectVects */
}TObject;

typedef struct TImage{
	int nbPts;
	int imageCenter[2];
	int **imagePts;
	double **imageVects;/* scaled orthographic projections */
	double **oldImageVects;/* projections at previous iteration step */
	double *epsilon;/* Corrections to scaled orthographic projections at each iteration */
}TImage;

typedef struct TCamera{
	int focalLength;
	double rotation[3][3];/* Rotation of SCENE in camera reference, NOT other way around */
	double translation[3];/* Translation of SCENE in camera reference */
}TCamera;


double **InitDoubleArray(int nbRows, int nbCols);


//{
#endif

