
#include "muhaha.h"

tM gM =
{
	.Y_Level = 128,
	.DeSat = 0,
	.X = 100,
	.Y = 100,
};

/*
typedef struct {//YUV 4:2:2
	u16 Y1:4;
	u16 U:2;
	u16 Y2:4;
	u16 V:2;
}tPix __attribute__((packed));
typedef struct {//YUV 4:2:2
	u16 U:2;
	u16 Y1:4;
	u16 V:2;
	u16 Y2:4;
}tPix __attribute__((packed));
*/
typedef struct {//YUV 4:2:2
	u08 Y;
	u08 U:4;
	u08 V:4;
}__attribute__((packed)) tPix;



#define  TPL_WIDTH       12      /* template width       */
#define  TPL_HEIGHT      12      /* template height      */
#define  WINDOW_WIDTH    24      /* search window width  */
#define  WINDOW_HEIGHT   24      /* search window height */
#define  THRESHOLD       0.3

#if 0

IplImage *frame, *tpl, *tm;

	/* create template image */
	tpl = cvCreateImage( cvSize( TPL_WIDTH, TPL_HEIGHT ),
	                     frame->depth, frame->nChannels );
	
	/* create image for template matching result */
	tm = cvCreateImage( cvSize( WINDOW_WIDTH  - TPL_WIDTH  + 1,
	                            WINDOW_HEIGHT - TPL_HEIGHT + 1 ),
	                    IPL_DEPTH_32F, 1 );
if( event == CV_EVENT_LBUTTONUP ) {
		object_x0 = x - ( TPL_WIDTH  / 2 );
		object_y0 = y - ( TPL_HEIGHT / 2 );
		
		cvSetImageROI( frame,
		               cvRect( object_x0,
		                       object_y0,
		                       TPL_WIDTH,
		                       TPL_HEIGHT ) );
		
		cvCopy( frame, tpl, NULL );
		
		cvResetImageROI( frame );
		
		/* template is available, start tracking! */
		fprintf( stdout, "Template selected. Start tracking... \n" );
		
		is_tracking = 1;
	}


void trackObject()
{
	CvPoint minloc, maxloc;
	double  minval, maxval;
	
	/* setup position of search window */
	
	int win_x0 = object_x0 - ( ( WINDOW_WIDTH  - TPL_WIDTH  ) / 2 );
	int win_y0 = object_y0 - ( ( WINDOW_HEIGHT - TPL_HEIGHT ) / 2 );
	
	/*
	 * Ooops, some bugs here.
	 * If the search window exceed the frame boundaries,
	 * it will trigger errors.
	 *
	 * Add some code to make sure that the search window
	 * is still within the frame.
	 */
	
	/* search object in search window */
	cvSetImageROI( frame,
	               cvRect( win_x0,
	                       win_y0,
	                       WINDOW_WIDTH,
	                       WINDOW_HEIGHT ) );
	
	cvMatchTemplate( frame, tpl, tm, CV_TM_SQDIFF_NORMED );
	cvMinMaxLoc( tm, &minval, &maxval, &minloc, &maxloc, 0 );
	cvResetImageROI( frame );
	
	/* if object found... */
	if( minval <= THRESHOLD ) {
		/* save object's current location */
		object_x0 = win_x0 + minloc.x;
		object_y0 = win_y0 + minloc.y;
		
		/* and draw a box there */
		cvRectangle( frame,
		             cvPoint( object_x0, object_y0 ),
		             cvPoint( object_x0 + TPL_WIDTH,
		                      object_y0 + TPL_HEIGHT ),
		             cvScalar( 0, 0, 255, 0 ), 1, 0, 0 );
	} else {
		/* if not found... */
		fprintf( stdout, "Lost object.\n" );
		is_tracking = 0;
	}
}

#endif

void	muhaha_Init	()
{
	
}

#define ddist2(x0,y0,x1,y1)	(((x1)-(x0))*((x1)-(x0)) + ((y1)-(y0))*((y1)-(y0)))

#define ddist(x0,y0,x1,y1) sqrt(ddist2(x0,y0,x1,y1))

#define dopix(_x,_y) ((tPix*)videoIn->framebuffer + ((_x) + (_y)*videoIn->width))
#define dnpix(_x,_y) ((tPix*)gM.pDst + ((_x) + (_y)*videoIn->width))
#define dpix(_x,_y) dopix(_x,_y)

si ay = 0, au = 0, av = 0, n = 0;

void	crap_ff	(si x, si y)
{
	if (x < 0 || y < 0 || x >= videoIn->width || y >= videoIn->height)
		return;
	
	if (	dnpix(x,y)->Y == 0xFF
		&& dnpix(x,y)->U == 0
		&& dnpix(x,y)->V == 0
	)
		return;
	
	if (ddist2(x,y,gM.X,gM.Y) >= 600)
		return;
	
	if (	dopix(x,y)->Y <= ay + 16
		//abs(dopix(x,y)->Y - ay) <= 16
		//	&& abs(dopix(x,y)->U - au) <= 12
		//	&& abs(dopix(x,y)->V - av) <= 12
	) {
		dnpix(x,y)->Y = 0xFF;
		dnpix(x,y)->U = 0x0;
		dnpix(x,y)->V = 0x0;
		crap_ff (x-1, y);
		crap_ff (x+1, y);
		crap_ff (x, y-1);
		crap_ff (x, y+1);
	}
};

void	muhaha	()
{
	tPix* pin = videoIn->framebuffer, *pin1;
	int x, y;
/*	for (y = 0; y < videoIn->height; ++y) {
		for (x = 0; x < videoIn->width; ++x) {
		//	pin->Y = 0;
		//	pin->Y1 = 0;
		//	pin->UV = 0;
			pin->U = 0;
			pin->V = 0;
			pin++;
		}
	}*/
	
/*	while (pin < (u00)videoIn->framebuffer + videoIn->width * (videoIn->height) * 2) {
	//	pin->Y = 0xFF;
		pin->U = 0x8;
		pin->V = 0x8;
		pin++;
	}/**/
	if (gM.DeSat) {
		for (y = 0; y < videoIn->height; ++y) {
			for (x = 0; x < videoIn->width-1; ++x) {
				pin = videoIn->framebuffer + (y*videoIn->width + x) * 2;
			/*	pin1 = pin+1;
			//	pin->Y = 0;
				si y = pin->Y, y1 = pin1->Y;
			//	pin->Y = abs(y1 - y);
				
				if (y < gM.Y_Level) {
			//		pin->Y = 0x00;
				}else {
			//		pin->Y = 0xFF;
				}*/
				
			//	pin->U = ((si)pin1->U - ((si)pin->U)) +7;
			//	pin->V = ((si)pin1->V - ((si)pin->V)) +8;
			//	pin->V = abs((si)pin1->V - (si)pin->V);
			//	pin->U = 0x7;
			//	pin->V = 0x8;
				dpix(x,y)->U = 0x7;
				dpix(x,y)->V = 0x8;
			}
		}/**/
	}
	
	memcpy(gM.pDst, videoIn->framebuffer, videoIn->width * (videoIn->height) * 2);
	
//	si ay = 0, au = 0, av = 0, n = 0;
	ay = 0;
	au = 0;
	av = 0;
	n = 0;
	for (y = gM.Y-8; y < gM.Y+8; ++y) {
		for (x = gM.X-8; x < gM.X+8; ++x) {
			ay += dpix(x,y)->Y;
			au += dpix(x,y)->U;
			av += dpix(x,y)->V;
			++n;
		}
	}
	ay /= n;
	au /= n;
	av /= n;/**/
	
//	ay = dnpix(gM.X,gM.Y)->Y;
//	au = dnpix(gM.X,gM.Y)->U;
//	av = dnpix(gM.X,gM.Y)->V;
	
//	printf ("n %d\n", n);
	#define dd 2
	for (y = gM.Y-dd; y <= gM.Y+dd; ++y) {
		for (x = gM.X-dd; x <= gM.X+dd; ++x) {
			crap_ff(x,y);
		}
	}
	#undef dd
	
//	memcpy(gM.pDst, videoIn->framebuffer, videoIn->width * (videoIn->height) * 2);
	
/*	for (y = 0; y < videoIn->height; ++y) {
		for (x = 0; x < videoIn->width-1; ++x) {
			if (abs(dpix(x,y)->Y - ay) <= 24
			//	&& abs(dpix(x,y)->U - au) <= 4
			//	&& abs(dpix(x,y)->V - av) <= 4
			) {
			//	dpix(x,y)->Y = 0xFF;
				dpix(x,y)->U = 0x0;
				dpix(x,y)->V = 0x0;
			}
		}
	}/**/
//	if (videoIn->formatIn != V4L2_PIX_FMT_YUYV)
//		printf ("videoIn->formatIn %d\n", videoIn->formatIn);
	
//	memcpy(gM.pDst, videoIn->framebuffer, videoIn->width * (videoIn->height) * 2);
	dnpix(gM.X,gM.Y)->Y = 0xFF;
	dnpix(gM.X,gM.Y)->U = 0xF;
	dnpix(gM.X,gM.Y)->V = 0xF;
}


int muhaha_eventThread(void *data)
{
	struct pt_data *gdata = (struct pt_data *) data;
	struct v4l2_control control;
	SDL_Surface *pscreen = *gdata->ptscreen;
	struct vdIn *videoIn = gdata->ptvideoIn;
	SDL_Event *sdlevent = gdata->ptsdlevent;
	SDL_Rect *drect = gdata->drect;
	SDL_mutex *affmutex = gdata->affmutex;
	
	int x, y;
	int mouseon = 0;
	int value = 0;
	int len = 0;
//	short incpantilt = INCPANTILT;
	int boucle = 0;
//	action_gui curr_action = A_VIDEO;
	while (videoIn->signalquit) {
		SDL_LockMutex(affmutex);
		
		float frmrate = gdata->frmrate;
		
		while (SDL_PollEvent(sdlevent)) {	//scan the event queue
			switch (sdlevent->type) {
			case SDL_KEYUP:
			case SDL_MOUSEBUTTONUP:
				mouseon = 0;
				boucle = 0;
				break;
			case SDL_MOUSEBUTTONDOWN:
				SDL_GetMouseState(&x, &y);
				gM.X = x; 
				gM.Y = y;
				break;
			case SDL_MOUSEMOTION:
				SDL_GetMouseState(&x, &y);
			//	curr_action = GUI_whichbutton(x, y, pscreen, videoIn);
				break;
			case SDL_VIDEORESIZE:
			/*	pscreen =
				      SDL_SetVideoMode(sdlevent->resize.w,
				                       sdlevent->resize.h, 0,
				                       SDL_VIDEO_Flags);
				drect->w = sdlevent->resize.w;
				drect->h = sdlevent->resize.h;*/
				break;
			case SDL_KEYDOWN:
				printf("Key down\n");
				if (sdlevent->key.keysym.sym == SDLK_s) {
					gM.DeSat = !gM.DeSat;
				}
			//	curr_action = GUI_keytoaction(sdlevent->key.keysym.sym);
			//	if (curr_action != A_VIDEO)
			//		mouseon = 1;
				break;
			case SDL_QUIT:
				printf("\nQuit signal received.\n");
				videoIn->signalquit = 0;
				break;
			}
		}			//end if poll
		SDL_UnlockMutex(affmutex);
		/* traiter les actions */
		value = 0;
		
	//	len = 100 * sizeof(char);	// as allocated in init_videoIn
	//	snprintf(videoIn->status, len, "%s, %.1f fps", title_act[curr_action].title, frmrate);
		
		SDL_Delay(50);
		//printf("fp/s %d\n",frmrate);
	}				//end main loop
	
	/* Close the stream capture file */
	if (videoIn->captureFile) {
		fclose(videoIn->captureFile);
		printf("Stopped raw stream capturing to stream.raw. %u bytes written for %u frames.\n",
		       videoIn->bytesWritten, videoIn->framesWritten);
	}
	/* Display stats for raw frame stream capturing */
	if (videoIn->rawFrameCapture == 2) {
		printf("Stopped raw frame stream capturing. %u bytes written for %u frames.\n",
		       videoIn->rfsBytesWritten, videoIn->rfsFramesWritten);
	}
}

