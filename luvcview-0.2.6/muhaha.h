
#ifndef inc_muhaha_h
#define inc_muhaha_h

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

#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_linalg.h>

typedef unsigned char u08;
typedef unsigned short u16;


typedef unsigned long u00;
typedef signed long s00;

typedef u00 ui;
typedef s00 si;


typedef struct {
	int Y_Level;
	
	u08 DeSat;
	
	float X, Y;
	float Ax, Ay, Aa, Axy, Ar;
	
	u08* pDst;
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

#endif

