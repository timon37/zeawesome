
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
	double x, y;
}tV2d;

typedef struct {
	float a, b, c, d;
}tM2f;

typedef struct {
	float x00, x01, x02;
	float x10, x11, x12;
	float x20, x21, x22;
}tM3f;

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
}tHead;

typedef struct {
	int Y_Level;
	
	u08 DeSat;
	
	tHead Head, HeadC;
//		tEye DotC, DotL, DotR;
//	}Cent;
//	tM3f HeadC, Head, HeadT;
	
	tEye Left, Right;
	tV2f L_Vec, R_Vec;
	
	tV2f GazeL, GazeR, Gaze;
	
	u08* pDst;
	
	SDL_Surface* pScreen;
	SDL_Overlay* pOverlay;
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

void svd(int m, int n, double **a, double **p, double *d, double **q);
//tV2d* normalize_point_set(tV2d* point_set, double *dis_scale, tV2d *nor_center, int num);

extern ui gM_edge_point_N;
extern tV2d gM_edge_point[204];
extern double pupil_param[5];//parameters of an ellipse {ellipse_a, ellipse_b, cx, cy, theta}; a & b is the major or minor axis; 

//{
#endif

