
#ifndef inc_muhaha_h
#define inc_muhaha_h
#if 0
}
#endif

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <sys/file.h>
#include <string.h>
#include <pthread.h>
#include <SDL.h>
#include <SDL_thread.h>
#include <SDL_audio.h>
#include <SDL_timer.h>
#include <linux/videodev2.h>
#include <sys/ioctl.h>
#include <sys/mman.h>
#include <errno.h>
#include <fcntl.h>
#include <time.h>
#include <sys/time.h>
#include <signal.h>
#include <X11/Xlib.h>
#include <SDL_syswm.h>
#include "v4l2uvc.h"
#include "gui.h"
#include "utils.h"
#include "color.h"

#include "uvcvideo.h"

#include <math.h>

#include <cv.h>

#include <X11/Xlib.h>
#include <X11/Xutil.h>
#include <X11/extensions/XInput2.h>
#include <X11/extensions/shape.h>
#include <X11/cursorfont.h>

#include <fann.h>

#include "gab.h"

//#include "xdo.h"

typedef unsigned char u08;
typedef unsigned short u16;
typedef unsigned int u32;

typedef signed char s08;

typedef unsigned long u00;
typedef signed long s00;

typedef u00 ui;
typedef s00 si;

typedef double	f00;

typedef struct {
	si x, y;
}tV2si;

typedef struct {
	f00 x, y;
}tV2f;

typedef struct {
	f00 x, y, z;
}tV3f;

typedef struct {
	f00 x, y, z, w;
}tV4f;

typedef struct {
	double x, y;
}tV2d;


typedef struct {
	f00 a, b, c, d;
}tM2f;

typedef struct {
	union {
	struct {
		f00 x00, x01, x02;
		f00 x10, x11, x12;
		f00 x20, x21, x22;
	};
	f00 f[3][3];
	};
}tM3f;

typedef struct {
	union {
	struct {
		f00 x00, x01, x02, x03;
		f00 x10, x11, x12, x13;
		f00 x20, x21, x22, x23;
		f00 x30, x31, x32, x33;
	};
	f00 f[4][4];
	f00 ff[16];
	};
}tM4f;

typedef struct {//YUV 4:2:2
	u08 Y;
	u08 U:4;
	u08 V:4;
}__attribute__((packed)) tPix;

typedef struct {
	u08 B, G, R, A;
}__attribute__((packed)) tRGBA;


#define CALIBRATIONPOINTS    9
typedef struct {
	u08 aPoint_State[CALIBRATIONPOINTS];  //0 not set, 1 set
	//tV2f  calipoints[CALIBRATIONPOINTS];       //conversion from eye to scene calibration points
	tV2f  scenecalipoints[CALIBRATIONPOINTS];  //captured (with mouse) calibration points
	//tV2f  pucalipoints[CALIBRATIONPOINTS];     //captured eye points while looking at the calibration points in the scene
	//tV2f  crcalipoints[CALIBRATIONPOINTS];     //captured corneal reflection points while looking at the calibration points in the scene
	tV2f  vectors[CALIBRATIONPOINTS];          //differences between the corneal reflection and pupil center
	
	double map_matrix[3][3];
	
	f00 aa, bb, cc, dd, ee;                       //pupil X coefficients
	f00 ff, gg, hh, ii, jj;			//pupil Y coefficients
	
	f00 centx, centy;                             // translation to center pupil data after biquadratics
	int inx, iny;                                   // translation to center pupil data before biquadratics
	f00 cmx[4], cmy[4];                           // corner correctioncoefficients
}tHomo;

enum {
	eEye_Fit_Default,
	eEye_Fit_SFit,
	eEye_Fit_S0,
	eEye_Fit_S2Fit,
	eEye_Fit_S3Fit_START,
	eEye_Fit_S3Fit_Point,
	eEye_Fit_S3Fit_Point4,
	eEye_Fit_S3Fit_Eye,
	eEye_Fit_S3Fit_Eye3,
	eEye_Fit_S3Fit_END,
	eEye_Fit_S4Fit_START,
	eEye_Fit_S4Fit_Edge0,
	eEye_Fit_S4Fit_Edge1,
	eEye_Fit_S4Fit_Edge2,
	eEye_Fit_S4Fit_END,
	eEye_Fit_S5_START,
	eEye_Fit_S5_Diff,
	eEye_Fit_S5_END,
	eEye_Fit_FF_START,
	eEye_Fit_FF_Diode,
	eEye_Fit_FF_END,
	eEye_Fit_C,
};

enum {
	eEye_PFit_S1,
	eEye_PFit_S2,
	eEye_PFit_S3,
	eEye_PFit_S3const,
	eEye_PFit_RANSAC,
};


enum {
	dScreen_MAX = 6,
	dEye_Screen_Cal_NUM = 3,
	dEye_Screen_Cal_LAST = dEye_Screen_Cal_NUM-1,
	
	dHead_Point_NUM = 1,
};

enum {//Calibration point state
	dEye_Screen_Cal_NULL,
	dEye_Screen_Cal_Set,	//set by the user
	dEye_Screen_Cal_Interp,	//interpolated between other points set by the user
	dEye_Screen_Cal_Extra,	//guess
};

typedef struct _tEye
{
	tV4f P0, P1, Vec;	//global position (start of retina), and global gaze vector, P1 = P0+Vec
	tV4f PS;		//global position calculated from two cameras
	
	tV2f P, OP, V; // Position, OldPosition, Velocity
	f00 Ax, Ay, Aa;
	
	si Pix_Dark, Pix_Bright;
	si S4_Pix_Bright;
	
	f00 Exp_R;
	f00 Min_R, Max_R, PFit_R;
	
	u08 Fit, Fit_Edge, PFit;
	f00 AngRes;
	f00 Fit_Scale, Fit_Trans;
	f00 S2Fit_Scale, S2Fit_Trans;
	
	struct {
		si Min_N;
		f00 Min_R;
		f00 Break_Dist;
		si Pix_Dark, Pix_Bright;
		si Pix_Diff_Start;
		si Pix_Diff_Min;
		
		f00 Diff_Dist;
	}S5;
	
/*	struct {
		u08 tmpID;
		si tmpNum;
		tV2f tmpP;
		ui Max_R;
		f00 Perf_R, MaxDiff_R;
		si Search_R, GSearch_R;
		si Y, Y_Marg, GY;
		
		tV2si Mark_P;
		u08* paMark;
	}FF;*/
	
	si Point_N, Point_Max;
	tV2f*	paPoint;
	
	tV2f GTF;	//Glint Transform Fix, couldn't resist GTF sounds too cool
	tV2f G0, G1;	//glint points on screen
	tV2f GV;		//normalized glint-retina vector
	tV4f GP;	//Lense sphere position based on G0, G1
	SDL_mutex *GV_mutex;
	f00 LR;	//Lense sphere radius
	
	
	tV2si LinView, CirView;
	
	struct {
		tV2f P, OP, V; // Position, OldPosition, Velocity
		
		
		struct {
			u08 tmpID;
			si tmpNum;
			tV2f tmpP;
			ui Max_R;
			f00 Perf_R, MaxDiff_R;
			si Search_R, GSearch_R;
			si Y, Y_Marg, GY;
			
			tV2si Mark_P;
			u08* paMark;
		}FF;
		
	}aCam[2];
	
	tHomo Homo;
	struct {
		ui Line_N;
		struct {
			tV4f P, V, P0, P1;
		}aLine[128];
		
		tV4f P;
		f00 R;
	}InHead;
	
	struct {
		struct {
			u08 State;
			si SX, SY;
			tV4f P;
		}aaCal[dEye_Screen_Cal_NUM][dEye_Screen_Cal_NUM];
		
		struct {
			tV2si Top, Left, Front;
		}View;
	}aScreen[dScreen_MAX];
	
	struct {
		tV2f P[2], GV[2];
	}tmp;
}tEye;



typedef struct {
	tV2f P;
	si W, H;	//in Pixels
	
	struct fann* pANN;
	
	
}tTrack_Point;

enum {
	Head_Type_eDot,
	Head_Type_ePoint,
	
	Head_Point_MAX = 3,
};
typedef struct {
	u08 Type;
	
	tV4f P, N, R;	//Position, Normal, Rotation about axes
	
	//tTrack_Point aPoint[dHead_Point_NUM];
	
	tEye DotC, DotL, DotR;
	tM3f M, MI;
	
	tM4f M4, M4_T, M4_R;
	tM4f M4I, M4I_T, M4I_R;
	
	f00 R_X, R_Y, R_Z, SRX, SRY;
	
	si Point_N;
	tEye aPoint[Head_Point_MAX];
	
	struct {
		tV4f PC, PL, PR;
		tV4f TInc, RInc;
		
		struct {
			tV4f P, Psum;
		}aPoint[Head_Point_MAX];
		si Sample_N;
		
		CvPOSITObject* cvPOSIT;
	}Mod;
	
	struct {
		tV4f P, N, R;	//Position, Normal, Rotation about axes
		
		tM4f M4, M4_T, M4_R;
		tM4f M4I, M4I_T, M4I_R;
		
		f00 R_X, R_Y, R_Z, SRX, SRY;
	}aCam[2];
}tHead;


typedef struct {
	ui Idx;
	u08 bGood;
	
	ui PixW, PixH;
	
	struct {
		tV4f C; //center
		f00 W, H;
	};/**/
	struct {
		tV4f C, Rot; //center
		f00 W, H;
	}Real;/**/
	
	Window Win;
	tV2si	Off;
	
	XWindowAttributes Attr;
	
	struct {
		ui sx, sy;
		si ix, iy;
	}Cal;
}tScreen;


typedef struct {
	ui Idx;
	struct vdIn* UVC;
	
	SDL_Window* SDL_Win;
	SDL_Surface* SDL_Surf;
	SDL_Thread* SDL_Thread;
	
	ui Full_W, Full_H;
	f00 Full_FOV;
	
	ui Image_W, Image_H;
	f00 Image_Zoom, Image_FOV, Image_FOV_W, Image_FOV_H;
	
	si Focus, Exposure, Zoom;
	
	CvMat* cvCam, *cvCamI, *cvDist;
	
	
	si View_W, View_H;
	f00 Proj_L, Proj_R, Proj_B, Proj_T;
	f00 Proj_W, Proj_H, Proj_N, Proj_F;
	
	tM4f Proj;
	tM4f  World, WorldI;
	tM4f  World_R, WorldI_R;
	
	tV4f T, R;
}tCam;


typedef struct {
	tV4f P;
	
}tLight;


typedef struct {
	u08 Type, Dirty;
	ui Win_W;
	Window Win;
	u32 Col;
}tMarker;


typedef struct _tDbg
{
	SDL_Window* SDL_Win;
	SDL_Surface* SDL_Surf;
	
	tM4f Proj, World;
	si View_X, View_Y, View_W, View_H;	//x y is the center
	
	f00 Scale, T_X, T_Y, R_X, R_Y;
	tV2si Off;
}tDbg;

enum {
	dM_X_Queue_M_Down,
	dM_X_Queue_M_Up,
	
	dM_X_Queue_NUM = 8,
	
	M_Dbg_NUM = 3
};

//actions should be triggered by 0 - off, 1 - X11 key grab, 2 - kernel zeawesome module
#define dM_Actions_Mode	2

typedef struct {
	u08 Eye_Line_Ray, Eye_GlintMode, GazeMode, GazeAvg, Pointer_Mode;
	
	u08 bGazeHold, bHead_Eye_LineDraw, bEye_S4_EdgeMark3_Micro;
	
	SDL_mutex *mutMainIsOut, *mutMainIsIn;
	
	int Y_Level;
	
	u08 DeSat;
	
	u08 Cam_N;
	tCam aCam[2];
	//si Cam_TCal, Cam_RCal;
	//si Head_Snap;
	
	SDL_sem* sWaitForCams, *sWaitForUpdate;
	
	tLight aLight[2];
	
	si Draw_X, Draw_Y, Draw_W, Draw_H;
	
	tHead Head, HeadC;
	
	tEye Left, Right;
	tV2f L_Vec, R_Vec;
	
	ui Screen_N;
	tScreen aScreen[dScreen_MAX];
	
	si Screen_CalIdx;
	u08 Cal_bRecord;
	si Cal_Record_Num;
	
	
	SDL_mutex* pGaze_Mutex;
	tV2f GazeL, GazeR, Gaze;
	
	
	tMarker CrossHair, CalPoint;
	
	tMarker ScreenPoint[dEye_Screen_Cal_NUM][dEye_Screen_Cal_NUM];
	
	u08* pDst;
	
	SDL_Surface* pScreen;
	//SDL_Overlay* pOverlay;
	
	struct {
		Display* pDisp;
		
	//	Window RootWin;
	//	XWindowAttributes RootAttr;
		
		XIDeviceInfo Pointer, Keyboard;
	//	xdo_t* pDo;
		
		struct {
			u08 Act;
			
			union {
				u08 Button;
			};
		}aQueue[dM_X_Queue_NUM];
		si Queue_C, Queue_N;	//current idx, number of remaning
	}X;
	
	struct {
		#define dact(name,press_stuff) name,
		#define dact2(name,press_stuff,release_stuff) name,
		#define dactD(name)	name,
		int
		#include "actions.h"
		NONE;
		#undef dact
		#undef dact2
		#undef dactD
	}Action;
	
	struct {
		#define dact(name,press_stuff)
		#define dact2(name,press_stuff,release_stuff)
		#define dactD(name)	name,
		int
		#include "actions.h"
		NONE;
		#undef dact
		#undef dact2
		#undef dactD
	}ActionD;
	
	f00 GazeAvg_3_MinAlpha, GazeAvg_3_MinDist, GazeAvg_3_Dist;
	struct {
		u08 State;
		tM4f Head_CR;
		tV2f Gaze;
		f00 SX, SY;	//scale
		f00 R_X, R_Y;
	}Micro;
	
	tDbg aDbg[M_Dbg_NUM];
	
	struct {
		f00 scale;
		f00 x, y, z;
	}tmp;
}tM;

extern tM gM;
//extern tPix gCol;
extern tRGBA gCol;
extern u32 gColARGB;
extern tCam* pcam;

extern si ay, au, av, n;

//extern struct vdIn *videoIn;

//For blocking main thread either at start or end of "muhahaha", but not within or outside of it
#define dSafe_Main_S()	\
	do {	\
		SDL_mutexP (gM.mutMainIsIn);	\
		SDL_mutexP (gM.mutMainIsOut);	\
	}while(0)

#define dSafe_Main_E()	\
	do {	\
		SDL_mutexV (gM.mutMainIsIn);	\
		SDL_mutexV (gM.mutMainIsOut);	\
	}while(0)




void	muhaha_Init		();
void	muhaha_DeInit	();
void	muhaha		();


void	Ehh_Draw_Line_2d	(tCam* pcam, si x0, si y0, si x1, si y1);



void	Eye_Init	(tEye* peye);
void	Eye_Conf	(tEye* peye);

void	Eye_CalcAYUV	(tEye* peye, tCam* pcam, si dd);

void	Eye_V_Pre		(tEye* peye, tCam* pcam);
void	Eye_V_Post		(tEye* peye, tCam* pcam);

void	Eye_Fit		(tEye* peye, tCam* pcam);

void	Eye_Draw		(tEye* peye, tCam* pcam);

tV2f	Eye_map_point	(tEye* peye, tV2f p);


void	Head_Eye_CalcP	(tHead* p, tEye* peye);
void	Screen_Eye_PreCal	(tScreen* p, tEye* peye);


int muhaha_eventThread(void *data);



static inline int		SDL_SemPostN	(SDL_sem* ps, si num)
{
	while (num--) {
		int ret = SDL_SemPost(ps);
		if (ret) {
			printf ("SDL_SemPostN: %s\n", SDL_GetError());
			abort();
		}
	}
	return 0;
}

static inline int		SDL_SemWaitN	(SDL_sem* ps, si num)
{
	while (num--) {
		int ret = SDL_SemWait(ps);
		if (ret) {
			printf ("SDL_SemWaitN: %s\n", SDL_GetError());
			abort();
		}
	}
	return 0;
}





struct pt_data {
	SDL_Surface **ptscreen;
	SDL_Event *ptsdlevent;
	SDL_Rect *drect;
	struct vdIn *ptvideoIn;
	f00 frmrate;
	SDL_mutex *affmutex;
} ptdata;








typedef struct _dyn_config_entry dyn_config_entry;
typedef struct _dyn_config_mapping dyn_config_mapping;
typedef struct _dyn_config dyn_config;
typedef enum DYN_ENTRY_TYPE dyn_entry_t;

enum { dyn_config_entry_name_max = 64 };
enum { dyn_config_entry_value_max = 64 };
struct _dyn_config_entry {
	dyn_config_entry *next;
	dyn_config_entry *child;
	char name[dyn_config_entry_name_max];
	char value[dyn_config_entry_value_max];
//	f00 value;
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

