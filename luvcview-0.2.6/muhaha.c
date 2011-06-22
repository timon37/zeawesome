
#include "muhaha.h"

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


tPix gCol;

tM gM =
{
	.Y_Level = 128,
	.DeSat = 0,
	.Head = {
		.DotC = {
			.P = {320, 50},
			.Ax = 20,
			.Ay = 20,
			.Aa = 0,
		//	.Aa = M_PI/4,
			
			.Exp_R = 15,
			
			.Fit_Scale = 0.2f,
			.Fit_Trans = 0.2f,
			.S2Fit_Scale = 0.8f,
			.S2Fit_Trans = 0.7f,
			
			.Proc_N = 128,
			.paProc = 0,
		//	.Axy = -4,
		//	.Ar = 100,
		},
		.DotL = {
			.P = {280, 100},
			.Ax = 20,
			.Ay = 20,
			.Aa = 0,
		//	.Aa = M_PI/4,
			
			.Exp_R = 15,
			
			.Fit_Scale = 0.2f,
			.Fit_Trans = 0.2f,
			.S2Fit_Scale = 0.8f,
			.S2Fit_Trans = 0.7f,
			
			.Proc_N = 128,
			.paProc = 0,
		//	.Axy = -4,
		//	.Ar = 100,
		},
		.DotR = {
			.P = {360, 100},
			.Ax = 20,
			.Ay = 20,
			.Aa = 0,
		//	.Aa = M_PI/4,
			
			.Exp_R = 15,
			
			.Fit_Scale = 0.2f,
			.Fit_Trans = 0.2f,
			.S2Fit_Scale = 0.8f,
			.S2Fit_Trans = 0.7f,
			
			.Proc_N = 128,
			.paProc = 0,
		//	.Axy = -4,
		//	.Ar = 100,
		},
	},
	.Left = {
		.P = {280, 240},
		.Ax = 30,
		.Ay = 30,
		.Aa = 0,
	//	.Aa = M_PI/4,
		
		.Exp_R = 20,
		
		.Fit_Scale = 0.1f,
		.Fit_Trans = 0.1f,
		.S2Fit_Scale = 0.8f,
		.S2Fit_Trans = 0.7f,
		
		.Proc_N = 128,
		.paProc = 0,
	//	.Axy = -4,
	//	.Ar = 100,
	},
	.Right = {
		.P = {360, 240},
		.Ax = 30,
		.Ay = 30,
		.Aa = 0,
	//	.Aa = M_PI/4,
		
		.Exp_R = 20,
		
		.Fit_Scale = 0.1f,
		.Fit_Trans = 0.1f,
		.S2Fit_Scale = 0.8f,
		.S2Fit_Trans = 0.7f,
		
		.Proc_N = 128,
		.paProc = 0,
	//	.Axy = -4,
	//	.Ar = 100,
	},
};

#define ddist2(x0,y0,x1,y1)	(((x1)-(x0))*((x1)-(x0)) + ((y1)-(y0))*((y1)-(y0)))

#define ddist(x0,y0,x1,y1) sqrt(ddist2(x0,y0,x1,y1))

#define dopix(_x,_y) ((tPix*)videoIn->framebuffer + ((si)(_x) + (si)(_y)*videoIn->width))
#define dnpix(_x,_y) ((tPix*)gM.pDst + ((si)(_x) + (si)(_y)*videoIn->width))
#define dpix(_x,_y) dnpix(_x,_y)


#define dspix(_x,_y) ((u32*)gM.pScreen->pixels + ((si)(_x) + (si)(_y)*gM.pScreen->pitch/4))

#define dpixout(_x,_y) ((_x) < 0 || (_y) < 0 || (_x) >= videoIn->width || (_y) >= videoIn->height)

float	V2f_dist	(tV2f* pv)
{
	return sqrtf(pv->x*pv->x + pv->y*pv->y);
}

void	V2f_sub_V2f	(tV2f* pv0, tV2f* pv1)
{
	float x, y;
	x = pv0->x - pv1->x;
	y = pv0->y - pv1->y;
	pv0->x = x;
	pv0->y = y;
}

void	V2f_mul_M2f		(tV2f* pv0, tM2f* pm0)
{
	float x, y;
	x = pv0->x*pm0->a + pv0->y*pm0->b;
	y = pv0->x*pm0->c + pv0->y*pm0->d;
	pv0->x = x;
	pv0->y = y;
}

void	V2f_mul_M3f		(tV2f* pv0, tM3f* pm0)
{
	float x, y;
	x = pv0->x*pm0->x00 + pv0->y*pm0->x01 + pm0->x02;
	y = pv0->x*pm0->x10 + pv0->y*pm0->x11 + pm0->x12;
	pv0->x = x;
	pv0->y = y;
}


void	M3f_Iden		(tM3f* pm0)
{
	pm0->x00 = 1;	pm0->x01 = 0;	pm0->x02 = 0;
	pm0->x10 = 0;	pm0->x11 = 1;	pm0->x12 = 0;
	pm0->x20 = 0;	pm0->x21 = 0;	pm0->x22 = 1;
}

void	M3f_mul_M3f		(tM3f* pm0, tM3f* pm1)
{
	tM3f m;
	m.x00 = pm0->x00*pm1->x00	+ pm0->x01*pm1->x10	+ pm0->x02*pm1->x20;
	m.x01 = pm0->x00*pm1->x01	+ pm0->x01*pm1->x11	+ pm0->x02*pm1->x21;
	m.x02 = pm0->x00*pm1->x02	+ pm0->x01*pm1->x12	+ pm0->x02*pm1->x22;
	
	m.x10 = pm0->x10*pm1->x00	+ pm0->x11*pm1->x10	+ pm0->x12*pm1->x20;
	m.x11 = pm0->x10*pm1->x01	+ pm0->x11*pm1->x11	+ pm0->x12*pm1->x21;
	m.x12 = pm0->x10*pm1->x02	+ pm0->x11*pm1->x12	+ pm0->x12*pm1->x22;
	
	m.x20 = pm0->x20*pm1->x00	+ pm0->x21*pm1->x10	+ pm0->x22*pm1->x20;
	m.x21 = pm0->x20*pm1->x01	+ pm0->x21*pm1->x11	+ pm0->x22*pm1->x21;
	m.x22 = pm0->x20*pm1->x02	+ pm0->x21*pm1->x12	+ pm0->x22*pm1->x22;
	
	*pm0 = m;
}


void	Ehh_Draw_Line_2d	(si x0, si y0, si x1, si y1);

void	Vec_Draw		(si x0, si y0, si x1, si y1)
{
	Ehh_Draw_Line_2d (x0, y0, x1, y1);
	
	si x, y;
	#define dd 2
	for (y = y1-dd; y <= y1+dd; ++y) {
		for (x = x1-dd; x <= x1+dd; ++x) {
			if (dpixout(x,y))
				continue;
			*dnpix(x,y) = gCol;
		}
	}
	#undef dd
}

void	V2f_DrawPosVec	(tV2f* ppos, tV2f* pvec)
{
	Vec_Draw (ppos->x, ppos->y, ppos->x+pvec->x, ppos->y+pvec->y);
}
void	V2f_DrawPosPos	(tV2f* ppos0, tV2f* ppos1)
{
	Vec_Draw (ppos0->x, ppos0->y, ppos1->x, ppos1->y);
}

void	Eye_Xset	(tEye* peye, float x)
{
	if (x < 20)
		x = 20;
	else if (x > videoIn->width-20)
		x = videoIn->width-20;
	peye->P.x = x;
}
void	Eye_Yset	(tEye* peye, float y)
{
	if (y < 20)
		y = 20;
	else if (y > videoIn->height-20)
		y = videoIn->height-20;
	peye->P.y = y;
}
void	Eye_XYset	(tEye* peye, float x, float y)
{
	Eye_Xset (peye, x);
	Eye_Yset (peye, y);
}


#if 0
//}

float a, b, c, d, e;                            //temporary storage of coefficients
//float aa, bb, cc, dd, ee;                       //pupil X coefficients
//float ff, gg, hh, ii, jj;			//pupil Y coefficients
 
//float centx, centy;                             // translation to center pupil data after biquadratics
//int inx, iny;                                   // translation to center pupil data before biquadratics
//float cmx[4], cmy[4];                           // corner correctioncoefficients

void dqfit(	tHomo* phomo,
		float x1, float y1, 
		float x2, float y2, 
		float x3, float y3, 
		float x4, float y4, 
		float x5, float y5,
		float X1, float X2, float X3, float X4, float X5 )
{
	float den;
	float x22,x32,x42,x52;    // squared terms 
	float y22,y32,y42,y52;
	
	phomo->inx = (int)x1;            // record eye tracker centering constants 
	phomo->iny = (int)y1;
	a = X1;                    // first coefficient 
	X2 -= X1;  X3 -= X1;       // center screen points 
	X4 -= X1;  X5 -= X1;
	x2 -= x1;  x3 -= x1;       // center eye tracker points 
	x4 -= x1;  x5 -= x1;
	y2 -= y1;  y3 -= y1;  
	y4 -= y1;  y5 -= y1;
	x22 = x2*x2; x32 = x3*x3;   // squared terms of biquadratic 
	x42 = x4*x4; x52 = x5*x5;
	y22 = y2*y2; y32 = y3*y3;
	y42 = y4*y4; y52 = y5*y5;
	
	//Cramer's rule solution of 4x4 matrix */
	den = -x2*y3*x52*y42-x22*y3*x4*y52+x22*y5*x4*y32-y22*x42*y3*x5-
		x32*y22*x4*y5-x42*x2*y5*y32+x32*x2*y5*y42-y2*x52*x4*y32+
		x52*x2*y4*y32+y22*x52*y3*x4+y2*x42*x5*y32+x22*y3*x5*y42-
		x32*x2*y4*y52-x3*y22*x52*y4+x32*y22*x5*y4-x32*y2*x5*y42+
		x3*y22*x42*y5+x3*y2*x52*y42+x32*y2*x4*y52+x42*x2*y3*y52-
		x3*y2*x42*y52+x3*x22*y4*y52-x22*y4*x5*y32-x3*x22*y5*y42;
	
	b =  (-y32*y2*x52*X4-X2*y3*x52*y42-x22*y3*X4*y52+x22*y3*y42*X5+
		y32*y2*x42*X5-y22*x42*y3*X5+y22*y3*x52*X4+X2*x42*y3*y52+
		X3*y2*x52*y42-X3*y2*x42*y52-X2*x42*y5*y32+x32*y42*y5*X2+
		X2*x52*y4*y32-x32*y4*X2*y52-x32*y2*y42*X5+x32*y2*X4*y52+
		X4*x22*y5*y32-y42*x22*y5*X3-x22*y4*y32*X5+x22*y4*X3*y52+
		y22*x42*y5*X3+x32*y22*y4*X5-y22*x52*y4*X3-x32*y22*y5*X4)/den;
	
	c =  (-x32*x4*y22*X5+x32*x5*y22*X4-x32*y42*x5*X2+x32*X2*x4*y52+
		x32*x2*y42*X5-x32*x2*X4*y52-x3*y22*x52*X4+x3*y22*x42*X5+
		x3*x22*X4*y52-x3*X2*x42*y52+x3*X2*x52*y42-x3*x22*y42*X5-
		y22*x42*x5*X3+y22*x52*x4*X3+x22*y42*x5*X3-x22*x4*X3*y52-
		x2*y32*x42*X5+X2*x42*x5*y32+x2*X3*x42*y52+x2*y32*x52*X4+
		x22*x4*y32*X5-x22*X4*x5*y32-X2*x52*x4*y32-x2*X3*x52*y42)/den;
	
	d = -(-x4*y22*y3*X5+x4*y22*y5*X3-x4*y2*X3*y52+x4*y2*y32*X5-
		x4*y32*y5*X2+x4*y3*X2*y52-x3*y22*y5*X4+x3*y22*y4*X5+
		x3*y2*X4*y52-x3*y2*y42*X5+x3*y42*y5*X2-x3*y4*X2*y52-
		y22*y4*x5*X3+y22*X4*y3*x5-y2*X4*x5*y32+y2*y42*x5*X3+
		x2*y3*y42*X5-y42*y3*x5*X2+X4*x2*y5*y32+y4*X2*x5*y32-
		y42*x2*y5*X3-x2*y4*y32*X5+x2*y4*X3*y52-x2*y3*X4*y52)/den;
	
	e = -(-x3*y2*x52*X4+x22*y3*x4*X5+x22*y4*x5*X3-x3*x42*y5*X2-
		x42*x2*y3*X5+x42*x2*y5*X3+x42*y3*x5*X2-y2*x42*x5*X3+
		x32*x2*y4*X5-x22*y3*x5*X4+x32*y2*x5*X4-x22*y5*x4*X3+
		x2*y3*x52*X4-x52*x2*y4*X3-x52*y3*x4*X2-x32*y2*x4*X5+
		x3*x22*y5*X4+x3*y2*x42*X5+y2*x52*x4*X3-x32*x5*y4*X2-
		x32*x2*y5*X4+x3*x52*y4*X2+x32*x4*y5*X2-x3*x22*y4*X5)/den;
}


int CalculateCalibration(tHomo* phomo)
{
	int i, j;
	float x, y, wx[9], wy[9];	//work data points
	int calx[10], caly[10];		//scene coordinate interpolation variables
	int eye_x[10], eye_y[10];	//scene coordinate interpolation variables
	
	// Place scene coordinates into calx and caly
	for(i = 0; i<9;i++) {
		calx[i] = phomo->scenecalipoints[i].x;  caly[i] = phomo->scenecalipoints[i].y;
	}
	
	// Set the last "tenth"  point
	calx[9] = phomo->scenecalipoints[0].x;  caly[9] = phomo->scenecalipoints[0].y;
	
	// Store pupil into eye_x and eye_y
	for(i = 0; i < 9; i++) {
		eye_x[i] = phomo->vectors[i].x;
		eye_y[i] = phomo->vectors[i].y;
	}
	
	// Solve X biquadratic
	dqfit(phomo,
		(float)eye_x[0],(float)eye_y[0],(float)eye_x[1],(float)eye_y[1],(float)eye_x[2],   
		(float)eye_y[2],(float)eye_x[3],(float)eye_y[3],(float)eye_x[4],(float)eye_y[4],
		(float)calx[0],(float)calx[1],(float)calx[2],(float)calx[3],(float)calx[4]);
	phomo->aa = a; phomo->bb = b; phomo->cc = c; phomo->dd = d; phomo->ee = e;
	
	// Solve Y biquadratic
	dqfit(phomo,
		(float)eye_x[0],(float)eye_y[0],(float)eye_x[1],(float)eye_y[1],(float)eye_x[2],
		(float)eye_y[2],(float)eye_x[3],(float)eye_y[3],(float)eye_x[4],(float)eye_y[4],
		(float)caly[0],(float)caly[1],(float)caly[2],(float)caly[3],(float)caly[4]);
	phomo->ff = a; phomo->gg = b; phomo->hh = c; phomo->ii = d; phomo->jj = e;
	
	// Biquadratic mapping of points
	for(i = 0; i < 9; i++) {
		x = (float)(eye_x[i] - phomo->inx);
		y = (float)(eye_y[i] - phomo->iny);
		wx[i] = phomo->aa+phomo->bb*x+phomo->cc*y+phomo->dd*x*x+phomo->ee*y*y;
		wy[i] = phomo->ff+phomo->gg*x+phomo->hh*y+phomo->ii*x*x+phomo->jj*y*y;
	}
	
	// Shift screen points to center for quadrant compute
	phomo->centx = wx[0];      
	phomo->centy = wy[0];
	
	// Normalize to center:
	for(i = 0; i < 9; i++) {
		wx[i] -= phomo->centx;
		wy[i] -= phomo->centy;
	}
	
	// Compute coefficents for each quadrant
	for(i = 0; i < 4; i++) {
		j = i + 5;
		phomo->cmx[i] = (calx[j]-wx[j]-phomo->centx)/(wx[j]*wy[j]);
		phomo->cmy[i] = (caly[j]-wy[j]-phomo->centy)/(wx[j]*wy[j]);
	}
	
	return 0;
}

tV2f map_point (tHomo* phomo, tV2f p)
{
	tV2f p2;
	int quad=0;
	float x1,y1,xx,yy;
	
	// correct eye position by recentering offset:
	x1 = (float) p.x;
	y1 = (float) p.y;
	
	// translate before biquadratic:
	x1 -= phomo->inx;
	y1 -= phomo->iny;
	
	// biquadratic mapping:
	xx = phomo->aa+phomo->bb*x1+phomo->cc*y1+phomo->dd*x1*x1+phomo->ee*y1*y1;
	yy = phomo->ff+phomo->gg*x1+phomo->hh*y1+phomo->ii*x1*x1+phomo->jj*y1*y1;
	
	// translate after biquadratic:
	x1 = xx - phomo->centx;
	y1 = yy - phomo->centy;
	
	// determine quadrant of point:
	if      (( x1<0 )&&( y1<0 )) quad = 0;
	else if (( x1>0 )&&( y1<0 )) quad = 1;
	else if (( x1<0 )&&( y1>0 )) quad = 2;
	else if (( x1>0 )&&( y1>0 )) quad = 3;
	
	// fix up by quadrant:
	p2.x = (int)(xx + x1*y1*phomo->cmx[quad]);
	p2.y = (int)(yy + x1*y1*phomo->cmy[quad]);
	
	return p2;
}

//{
#endif


#if 1
//}

void svd(int m, int n, double **a, double **p, double *d, double **q);
tV2d* normalize_point_set(tV2d* point_set, double *dis_scale, tV2d *nor_center, int num);
static double   radius(double u, double v);
/*
#define CALIBRATIONPOINTS    9
//tV2f  calipoints[CALIBRATIONPOINTS];       //conversion from eye to scene calibration points
tV2f  phomo->scenecalipoints[CALIBRATIONPOINTS];  //captured (with mouse) calibration points
//tV2f  pucalipoints[CALIBRATIONPOINTS];     //captured eye points while looking at the calibration points in the scene
//tV2f  crcalipoints[CALIBRATIONPOINTS];     //captured corneal reflection points while looking at the calibration points in the scene
tV2f  phomo->vectors[CALIBRATIONPOINTS];          //differences between the corneal reflection and pupil center

double map_matrix[3][3];
*/
//------------ map pupil coordinates to screen coordinates ---------/
tV2f map_point (tHomo* phomo, tV2f p)
{
  tV2f p2;
  double z = phomo->map_matrix[2][0]*p.x + phomo->map_matrix[2][1]*p.y + phomo->map_matrix[2][2];
  p2.x = (int)((phomo->map_matrix[0][0]*p.x + phomo->map_matrix[0][1]*p.y + phomo->map_matrix[0][2])/z);
  p2.y = (int)((phomo->map_matrix[1][0]*p.x + phomo->map_matrix[1][1]*p.y + phomo->map_matrix[1][2])/z);
  return p2;
}


// r is result matrix
void affine_matrix_inverse(double a[][3], double r[][3])
{
  double det22 = a[0][0]*a[1][1] - a[0][1]*a[1][0];
  r[0][0] = a[1][1]/det22;
  r[0][1] = -a[0][1]/det22;
  r[1][0] = -a[1][0]/det22;
  r[1][1] = a[0][0]/det22;

  r[2][0] = r[2][1] = 0;
  r[2][2] = 1/a[2][2];

  r[0][2] = -r[2][2] * (r[0][0]*a[0][2] + r[0][1]*a[1][2]);
  r[1][2] = -r[2][2] * (r[1][0]*a[0][2] + r[1][1]*a[1][2]);
}

// r is result matrix
void matrix_multiply33(double a[][3], double b[][3], double r[][3])
{
  int i, j;
  double result[9];
  double v = 0;
  for (j = 0; j < 3; j++) {
    for (i = 0; i < 3; i++) {
      v = a[j][0]*b[0][i];
      v += a[j][1]*b[1][i];
      v += a[j][2]*b[2][i];
      result[j*3+i] = v;
    }
  }
  for (i = 0; i < 3; i++) {
    r[i][0] = result[i*3];
    r[i][1] = result[i*3+1];
    r[i][2] = result[i*3+2];
  }
}

void CalculateCalibration(tHomo* phomo)
{
  int i, j;
  tV2d cal_scene[9], cal_eye[9];
  tV2d scene_center, eye_center, *eye_nor, *scene_nor;
  double dis_scale_scene, dis_scale_eye;  

  for (i = 0; i < 9; i++) {
    cal_scene[i].x = phomo->scenecalipoints[i].x;  
    cal_scene[i].y = phomo->scenecalipoints[i].y;
    cal_eye[i].x = phomo->vectors[i].x;
    cal_eye[i].y = phomo->vectors[i].y;
  }

  scene_nor = normalize_point_set(cal_scene, &dis_scale_scene, &scene_center, CALIBRATIONPOINTS);
  eye_nor = normalize_point_set(cal_eye, &dis_scale_eye, &eye_center, CALIBRATIONPOINTS);

  printf("normalize_point_set end\n");
  printf("scene scale:%lf  center (%lf, %lf)\n", dis_scale_scene, scene_center.x, scene_center.y);
  printf("eye scale:%lf  center (%lf, %lf)\n", dis_scale_eye, eye_center.x, eye_center.y);

  const int homo_row=18, homo_col=9;
  double A[homo_row][homo_col];
  int M = homo_row, N = homo_col; //M is row; N is column
  double **ppa = (double**)malloc(sizeof(double*)*M);
  double **ppu = (double**)malloc(sizeof(double*)*M);
  double **ppv = (double**)malloc(sizeof(double*)*N);
  double pd[homo_col];
  for (i = 0; i < M; i++) {
    ppa[i] = A[i];
    ppu[i] = (double*)malloc(sizeof(double)*N);
  }
  for (i = 0; i < N; i++) {
    ppv[i] = (double*)malloc(sizeof(double)*N);
  }

  for (j = 0;  j< M; j++) {
    if (j%2 == 0) {
      A[j][0] = A[j][1] = A[j][2] = 0;
      A[j][3] = -eye_nor[j/2].x;
      A[j][4] = -eye_nor[j/2].y;
      A[j][5] = -1;
      A[j][6] = scene_nor[j/2].y * eye_nor[j/2].x;
      A[j][7] = scene_nor[j/2].y * eye_nor[j/2].y;
      A[j][8] = scene_nor[j/2].y;
    } else {
      A[j][0] = eye_nor[j/2].x;
      A[j][1] = eye_nor[j/2].y;
      A[j][2] = 1;
      A[j][3] = A[j][4] = A[j][5] = 0;
      A[j][6] = -scene_nor[j/2].x * eye_nor[j/2].x;
      A[j][7] = -scene_nor[j/2].x * eye_nor[j/2].y;
      A[j][8] = -scene_nor[j/2].x;
    }
  }

  printf("normalize_point_set end\n");

  svd(M, N, ppa, ppu, pd, ppv);
  int min_d_index = 0;
  for (i = 1; i < N; i++) {
    if (pd[i] < pd[min_d_index])
      min_d_index = i;
  }

  for (i = 0; i < N; i++) {
      phomo->map_matrix[i/3][i%3] = ppv[i][min_d_index];  //the column of v that corresponds to the smallest singular value,
                                                //which is the solution of the equations
  }

  double T[3][3] = {0}, T1[3][3] = {0};
  printf("\nT1: \n");
  for (j = 0; j < 3; j++) {
    for (i = 0; i < 3; i++) {
      printf("%8lf ", T1[j][i]);
    }
    printf("\n");
  }  

  T[0][0] = T[1][1] = dis_scale_eye;
  T[0][2] = -dis_scale_eye*eye_center.x;
  T[1][2] = -dis_scale_eye*eye_center.y;
  T[2][2] = 1;

  printf("\nmap_matrix: \n");
  for (j = 0; j < 3; j++) {
    for (i = 0; i < 3; i++) {
      printf("%8lf ", phomo->map_matrix[j][i]);
    }
    printf("\n");
  }   
  printf("\nT: \n");
  for (j = 0; j < 3; j++) {
    for (i = 0; i < 3; i++) {
      printf("%8lf ", T[j][i]);
    }
    printf("\n");
  }  

  matrix_multiply33(phomo->map_matrix, T, phomo->map_matrix); 

  T[0][0] = T[1][1] = dis_scale_scene;
  T[0][2] = -dis_scale_scene*scene_center.x;
  T[1][2] = -dis_scale_scene*scene_center.y;
  T[2][2] = 1;

  printf("\nmap_matrix: \n");
  for (j = 0; j < 3; j++) {
    for (i = 0; i < 3; i++) {
      printf("%8lf ", phomo->map_matrix[j][i]);
    }
    printf("\n");
  } 
  printf("\nT: \n");
  for (j = 0; j < 3; j++) {
    for (i = 0; i < 3; i++) {
      printf("%8lf ", T[j][i]);
    }
    printf("\n");
  }   

  affine_matrix_inverse(T, T1);
  matrix_multiply33(T1, phomo->map_matrix, phomo->map_matrix);

  printf("\nmap_matrix: \n");
  for (j = 0; j < 3; j++) {
    for (i = 0; i < 3; i++) {
      printf("%8lf ", phomo->map_matrix[j][i]);
    }
    printf("\n");
  }   

  for (i = 0; i < M; i++) {
    free(ppu[i]);
  }
  for (i = 0; i < N; i++) {
    free(ppv[i]);
  }
  free(ppu);
  free(ppv);
  free(ppa);
      
  free(eye_nor);
  free(scene_nor);
  printf("\nfinish calculate calibration\n");
}

//{
#endif


tV2f	point_clip	(tV2f* ppoint)
{
	if (ppoint->x < 0)
		ppoint->x = 0;
	else if (ppoint->x >= videoIn->width)
		ppoint->x = videoIn->width-1;
		
	if (ppoint->y < 0)
		ppoint->y = 0;
	else if (ppoint->y >= videoIn->height)
		ppoint->y = videoIn->height-1;
}

si ay = 0, au = 0, av = 0, n = 0;

void	Eye_crap_ff		(tEye* peye, si x, si y)
{
	if (x < 0 || y < 0 || x >= videoIn->width || y >= videoIn->height)
		return;
	
	if (	1//dnpix(x,y)->Y == 0xFF
		&& dnpix(x,y)->U == 0
		&& dnpix(x,y)->V == 0
	)
		return;
	
	if (ddist2(x,y,peye->P.x,peye->P.y) >= peye->Exp_R*peye->Exp_R*2)
		return;
	
	if (	dopix(x,y)->Y <= ay + 2
		//abs(dopix(x,y)->Y - ay) <= 16
		//	&& abs(dopix(x,y)->U - au) <= 12
		//	&& abs(dopix(x,y)->V - av) <= 12
	) {
	//	dnpix(x,y)->Y = 0xFF;
		dnpix(x,y)->U = 0x0;
		dnpix(x,y)->V = 0x0;
		Eye_crap_ff (peye, x-1, y);
		Eye_crap_ff (peye, x+1, y);
		Eye_crap_ff (peye, x, y-1);
		Eye_crap_ff (peye, x, y+1);
	}
};

void	Eye_FF		(tEye* peye)
{
	si x, y;
//	si ay = 0, au = 0, av = 0, n = 0;
	ay = 0;
	au = 0;
	av = 0;
	n = 0;
	#define dd 4
	for (y = peye->P.y-dd; y <= peye->P.y+dd; ++y) {
		for (x = peye->P.x-dd; x <= peye->P.x+dd; ++x) {
			ay += dpix(x,y)->Y;
			au += dpix(x,y)->U;
			av += dpix(x,y)->V;
			++n;
		}
	}
	#undef dd
	ay /= n;
	au /= n;
	av /= n;/**/
	
//	ay = dnpix(peye->P.x,peye->P.y)->Y;
//	au = dnpix(peye->P.x,peye->P.y)->U;
//	av = dnpix(peye->P.x,peye->P.y)->V;
	
//	printf ("n %d\n", n);
	#define dd 2
	for (y = peye->P.y-dd; y <= peye->P.y+dd; ++y) {
		for (x = peye->P.x-dd; x <= peye->P.x+dd; ++x) {
			Eye_crap_ff (peye, x, y);
		}
	}
	#undef dd
	
}

void	Eye_CalcAYUV	(tEye* peye)
{
	si x, y;
//	si ay = 0, au = 0, av = 0, n = 0;
	ay = 0;
	au = 0;
	av = 0;
	n = 0;
	#define dd 8
	for (y = peye->P.y-dd; y <= peye->P.y+dd; ++y) {
		for (x = peye->P.x-dd; x <= peye->P.x+dd; ++x) {
			ay += dpix(x,y)->Y;
			au += dpix(x,y)->U;
			av += dpix(x,y)->V;
			++n;
		}
	}
	#undef dd
	ay /= n;
	au /= n;
	av /= n;
}

void	Eye_EdgeMark	(tEye* peye)
{
	si x, y;
	for (y = 0; y < videoIn->height; ++y) {
		for (x = 0; x < videoIn->width-1; ++x) {
			if (	1//dnpix(x,y)->Y == 0xFF
				&& dnpix(x,y)->U == 0
				&& dnpix(x,y)->V == 0
			) {
				float t, dx, dy, xx, yy;
				dx = peye->P.x - x;
				dy = peye->P.y - y;
				t = atan2(dy,dx) - peye->Aa;
				
				xx = (peye->Ax * cos(t) * cos(peye->Aa) - peye->Ay * sin(t) * sin(peye->Aa));
				yy = (peye->Ax * cos(t) * sin(peye->Aa) + peye->Ay * sin(t) * cos(peye->Aa));
				
				if (dx*dx+dy*dy > xx*xx+yy*yy) {
				//	printf ("f1\n");
				//	dnpix(x,y)->U = 0xF;
				//	dnpix(x,y)->V = 0xF;
				}//else
				{
					#if 1
					{	float nx = x, ny = y, d = sqrt(dx*dx+dy*dy);
					//	printf ("d %f\n", d);
					/*	if (t >= -M_PI_4 && t < M_PI_4) {
							nx = x-1;
							ny = y;
						}else if (t >= M_PI_4 && t < 3*M_PI_4) {
							nx = x;
							ny = y-1;
						}else if (t >= -3*M_PI_4 && t < -M_PI_4) {
							nx = x;
							ny = y+1;
						}else {
							nx = x+1;
							ny = y;
						}/**/
						while (1) {
							float tx = x + d*cos(t);
							float ty = y + d*sin(t);
						//	printf ("tx %f ty %f\n", tx, ty);
							if (	1//dnpix(nx,ny)->Y == 0xFF
								&& dnpix(tx,ty)->U == 0
								&& dnpix(tx,ty)->V == 0
							//	0
							//	|| (dnpix(tx,ny)->U != 0 && dnpix(tx,ty)->U != 0xF)
							//	|| (dnpix(tx,ny)->V != 0 && dnpix(tx,ty)->V != 0xF)
							) {
							//	printf ("ehhh\n");
								nx = tx+cos(t);
								ny = ty+sin(t);
							}
							d += 1;
							if (d >= 100) {
								break;
							}
						}
					//	dnpix(nx,ny)->U = 0xF;
					//	dnpix(nx,ny)->V = 0xF;
						if (	0//dnpix(nx,ny)->Y == 0xFF
							|| (dnpix(nx,ny)->U != 0 && dnpix(nx,ny)->U != 0xF)
							|| (dnpix(nx,ny)->V != 0 && dnpix(nx,ny)->V != 0xF)
						) {
							if (dpixout(nx,ny))
								continue;
						//	printf ("f2\n");
						//	printf ("nx %f ny %f\n", nx, ny);
							dnpix(nx,ny)->Y = 0xFF;
							dnpix(nx,ny)->U = 0xF;
							dnpix(nx,ny)->V = 0xF;
						//	printf ("Ax %f Ay %f\n", cos(t), sin(t));
						}
					}
					#endif
				}
			}
		}
	}/**/
}


void	Eye_SFit		(tEye* peye)
{
	si x, y;
/*	for (y = 0; y < videoIn->height; ++y) {
		for (x = 0; x < videoIn->width; ++x) {
			if (	1//dnpix(x,y)->Y == 0xFF
				&& dnpix(x,y)->U == 0xF
				&& dnpix(x,y)->V == 0xF
			) {
				float t, dx, dy, xx, yy;
				dx = peye->P.x - x;
				dy = peye->P.y - y;
				t = atan2(dy,dx) - peye->Aa;
				
				xx = (peye->Ax * cos(t) * cos(peye->Aa) - peye->Ay * sin(t) * sin(peye->Aa));
				yy = (peye->Ax * cos(t) * sin(peye->Aa) + peye->Ay * sin(t) * cos(peye->Aa));
				
				float diff = sqrt(dx*dx+dy*dy) - sqrt(xx*xx+yy*yy);
				if (fabs(diff) >= 0.01f) {
				//	printf ("f1\n");
					peye->Ax += peye->Fit_Scale*diff*fabs(cos(t));
					peye->Ay += peye->Fit_Scale*diff*fabs(sin(t));
					
					Eye_Xset (peye, peye->P.x - peye->Fit_Trans*diff*cos(t));
					Eye_Yset (peye, peye->P.y - peye->Fit_Trans*diff*sin(t));
				}
			}
		}
	}/**/
	for (y = 0; y < videoIn->height; ++y) {
		for (x = 0; x < videoIn->width; ++x) {
			if (	1//dnpix(x,y)->Y == 0xFF
				&& dnpix(x,y)->U == 0xF
				&& dnpix(x,y)->V == 0xF
			) {
				float t, dx, dy, xx, yy;
				dx = x - peye->P.x;
				dy = y - peye->P.y;
				t = atan2(dy,dx)/* - peye->Aa*/;
				
				xx = (peye->Ax * cos(t) * cos(peye->Aa) - peye->Ay * sin(t) * sin(peye->Aa));
				yy = (peye->Ax * cos(t) * sin(peye->Aa) + peye->Ay * sin(t) * cos(peye->Aa));
				
				float diff = sqrt(dx*dx+dy*dy) - sqrt(xx*xx+yy*yy);
				if (fabs(diff) >= 0.1f) {
				//	printf ("f1\n");
					peye->Ax += peye->Fit_Scale*diff*fabs(cos(t));
					peye->Ay += peye->Fit_Scale*diff*fabs(sin(t));
					
					Eye_Xset (peye, peye->P.x + peye->Fit_Trans*(diff )*cos(t));
					Eye_Yset (peye, peye->P.y + peye->Fit_Trans*(diff )*sin(t));
				}
			}
		}
	}/**/
}

#if 0
void	Eye_CFit_EdgeGet	(tEye* peye, float reqt, si* px, si* py)
{
	si border = peye->Exp_R*peye->Exp_R*3;
	struct {
		float t;
		si x, y;
	}min;
	min.t = reqt + 100;
	si x, y;
	for (y = peye->P.y - border; y < peye->P.y + border; ++y) {
		for (x = peye->P.x - border; x < peye->P.x + border; ++x) {
		/*	if (	1//dnpix(x,y)->Y == 0xFF
				&& dnpix(x,y)->U == 0xF
				&& dnpix(x,y)->V == 0xF
			) {
				float t, dx, dy, xx, yy;
				dx = peye->P.x - x;
				dy = peye->P.y - y;
				t = atan2(dy,dx) - peye->Aa;
				if (fabsf(t - reqt) < fabsf(min.t - reqt)) {
					min.t = t;
					min.x = x;
					min.y = y;
				}
			}*/
		}
	}/**/
	*px = min.x;
	*py = min.y;
}
#endif
si	Eye_S2Fit_EdgeGet	(tEye* peye, float t, si* px, si* py)
{
//	si border = peye->Exp_R*peye->Exp_R*3;
	float x = peye->P.x, y = peye->P.y;
	
	for (; ddist2(x,y,peye->P.x,peye->P.y) < peye->Exp_R*peye->Exp_R*3; ) {
		if (	1//dnpix(x,y)->Y == 0xFF
			&& dnpix(x,y)->U == 0xF
			&& dnpix(x,y)->V == 0xF
		) {
			*px = x;
			*py = y;
			return 1;
		}
		x += cos(t);
		y += sin(t);
	}/**/
	return 0;
}
void	Eye_S2Fit		(tEye* peye)
{
//	Eye_SFit (peye);	return;
	si border = 4*peye->Exp_R;
/*	float a = 0;
	for (a = 0; a < M_PI; a += M_PI/180.0f) {
		si x0, y0;
		Eye_CFit_EdgeGet (peye, a, &x0, &y0);
	//	dnpix(x0,y0)->Y = 0x0;
	//	dnpix(x0,y0)->U = 0xF;
	//	dnpix(x0,y0)->V = 0xF;
	/*	float t, dx, dy, xx, yy;
		dx = peye->P.x - x;
		dy = peye->P.y - y;
		t = atan2(dy,dx) - peye->Aa;
		
		xx = (peye->Ax * cos(t) * cos(peye->Aa) - peye->Ay * sin(t) * sin(peye->Aa));
		yy = (peye->Ax * cos(t) * sin(peye->Aa) + peye->Ay * sin(t) * cos(peye->Aa));
		
		float diff = sqrt(dx*dx+dy*dy) - sqrt(xx*xx+yy*yy);
		if (fabs(diff) >= 0.01f) {
		//	printf ("f1\n");
			peye->Ax += peye->Fit_Scale*diff*fabs(cos(t));
			peye->Ay += peye->Fit_Scale*diff*fabs(sin(t));
			
			Eye_Xset (peye, peye->P.x - peye->Fit_Trans*diff*cos(t));
			Eye_Yset (peye, peye->P.y - peye->Fit_Trans*diff*sin(t));
		}*/
//	}
	float ax = 0, ay = 0;
	si n = 0;
	
	si x0, y0;
	for (y0 = peye->P.y - border; y0 < peye->P.y + border; ++y0) {		//will have to stop this crap soon
//	for (y0 = peye->P.y + border; y0 > peye->P.y - border; --y0) {
		for (x0 = peye->P.x - border; x0 < peye->P.x + border; ++x0) {
			if (	1//dnpix(x0,y0)->Y == 0xFF
				&& dnpix(x0,y0)->U == 0xF
				&& dnpix(x0,y0)->V == 0xF
			) {
				float t0, dx0, dy0;
				dx0 = x0 - peye->P.x;
				dy0 = y0 - peye->P.y;
				t0 = atan2(dy0,dx0)/* - peye->Aa*/;
				float t1 = t0 - M_PI;
				si x1, y1;
				
			//	printf ("Heee? eye %f %f %f %f\n", peye->P.x, peye->P.y, dx0, dy0);
				
				if (Eye_S2Fit_EdgeGet (peye, t1, &x1, &y1)) {
				//	printf ("Heee? eye %ld %ld   %ld %ld\n", x0 - (si)peye->P.x, y0 - (si)peye->P.y, x1 - (si)peye->P.x, y1 - (si)peye->P.y);
				/*	dnpix(x0,y0)->Y = 0x0;
					dnpix(x0,y0)->U = 0x0;
					dnpix(x0,y0)->V = 0x0;
					
					dnpix(x1,y1)->Y = 0xFF;
					dnpix(x1,y1)->U = 0xF;
					dnpix(x1,y1)->V = 0xF;/**/
					if (0) {
						float dx1, dy1;
						dx1 = x1 - peye->P.x;
						dy1 = y1 - peye->P.y;
						
						float ox0 = (peye->Ax * cos(t0) * cos(peye->Aa) - peye->Ay * sin(t0) * sin(peye->Aa));
						float ox1 = (peye->Ax * cos(t1) * cos(peye->Aa) - peye->Ay * sin(t1) * sin(peye->Aa));
						
						float oy0 = (peye->Ax * cos(t0) * sin(peye->Aa) + peye->Ay * sin(t0) * cos(peye->Aa));
						float oy1 = (peye->Ax * cos(t1) * sin(peye->Aa) + peye->Ay * sin(t1) * cos(peye->Aa));
						
					/*	float diff = sqrt(dx0*dx0+dy0*dy0) - sqrt(ox0*ox0+oy0*oy0);
						if (fabs(diff) >= 0.01f) {
						//	printf ("f1\n");
							peye->Ax += peye->Fit_Scale*diff*fabs(cos(t0));
							peye->Ay += peye->Fit_Scale*diff*fabs(sin(t0));
							
							Eye_Xset (peye, peye->P.x + peye->Fit_Trans*(diff )*cos(t0));
							Eye_Yset (peye, peye->P.y + peye->Fit_Trans*(diff )*sin(t0));
						}/**/
					/*	float diff = sqrt(dx1*dx1+dy1*dy1) - sqrt(ox1*ox1+oy1*oy1);
						if (fabs(diff) >= 0.01f) {
						//	printf ("f1\n");
							peye->Ax += peye->Fit_Scale*diff*fabs(cos(t1));
							peye->Ay += peye->Fit_Scale*diff*fabs(sin(t1));
							
							Eye_Xset (peye, peye->P.x + peye->Fit_Trans*(diff )*cos(t1));
							Eye_Yset (peye, peye->P.y + peye->Fit_Trans*(diff )*sin(t1));
						}/**/
					}else {
						float ox0 = (peye->Ax * cos(t0) * cos(peye->Aa) - peye->Ay * sin(t0) * sin(peye->Aa));
						float ox1 = (peye->Ax * cos(t1) * cos(peye->Aa) - peye->Ay * sin(t1) * sin(peye->Aa));
						
						float oy0 = (peye->Ax * cos(t0) * sin(peye->Aa) + peye->Ay * sin(t0) * cos(peye->Aa));
						float oy1 = (peye->Ax * cos(t1) * sin(peye->Aa) + peye->Ay * sin(t1) * cos(peye->Aa));
						
					//	printf ("Heee? pic %f %f   %f %f\n", x0 - peye->P.x, y0 - peye->P.y, x1 - peye->P.x, y1 - peye->P.y);
					//	printf ("Heee? eli %f %f   %f %f\n", ox0, oy0, ox1, oy1);
						
						float sdx = fabsf(x1 - x0) - fabsf(ox1 - ox0);
						float sdy = fabsf(y1 - y0) - fabsf(oy1 - oy0);
						
						ax += fabsf((x1 - x0)*cos(t0));
						ay += fabsf((y1 - y0)*sin(t0));
						++n;
						
						float pdx = 0.2f + (x0+x1)*0.5f - (2*peye->P.x + ox0+ox1)*0.5f;
						float pdy = -0.1f + (y0+y1)*0.5f - (2*peye->P.y + oy0+oy1)*0.5f;
						
						Eye_Xset (peye, peye->P.x + peye->S2Fit_Trans * pdx );// *cos(t));
						Eye_Yset (peye, peye->P.y + peye->S2Fit_Trans * pdy );// *sin(t));
						
					//	peye->Ax += sdx;
					//	peye->Ay += sdy;
						
					//	peye->Ax += 0.5f*sdx*fabs(cos(t0));
					//	peye->Ay += 0.5f*sdy*fabs(sin(t0));
						
					//	peye->Ax += peye->Fit_Scale * sdx;
					//	peye->Ay += peye->Fit_Scale * sdy;
						
					//	peye->Ax += peye->Fit_Scale * ( sdx*fabs(cos(t0)) );
					//	peye->Ay += peye->Fit_Scale * ( sdy*fabs(sin(t0)) );
						
					//	peye->Ax += peye->Fit_Scale * ( sdx* (fabs(cos(t0))*fabs(cos(t0))) );
					//	peye->Ay += peye->Fit_Scale * ( sdy* (fabs(sin(t0))*fabs(sin(t0))) );
						
					//	peye->Ax += peye->Fit_Scale * ( sdx* (fabs(cos(t0))+fabs(sin(t0))) );
					//	peye->Ay += peye->Fit_Scale * ( sdy* (fabs(cos(t0))+fabs(sin(t0))) );
						/*
						float pdx = (x0+x1)*0.5f - (ox0+ox1)*0.5f;
						float pdy = (y0+y1)*0.5f - (oy0+oy1)*0.5f;
						
						Eye_Xset (peye, peye->P.x + peye->Fit_Trans * pdx );// *cos(t));
						Eye_Yset (peye, peye->P.y + peye->Fit_Trans * pdy );// *sin(t));
						/**/
						
					//	Eye_Xset (peye, peye->P.x + peye->Fit_Trans * ( (x0+x1)*0.5f - peye->P.x-(ox0+ox1)*0.5f ) );// *cos(t));
					//	Eye_Yset (peye, peye->P.y + peye->Fit_Trans * ( (y0+y1)*0.5f - peye->P.y-(oy0+oy1)*0.5f ) );// *sin(t));
						
					/*	dnpix(ox0, oy0)->Y = 0x0;
						dnpix(ox0, oy0)->U = 0x0;
						dnpix(ox0, oy0)->V = 0x0;
						
						dnpix(ox1, oy1)->Y = 0xFF;
						dnpix(ox1, oy1)->U = 0xF;
						dnpix(ox1, oy1)->V = 0xF;/**/
					}
				}else {
				//	printf ("Heee? eye %f %f %ld %ld\n", peye->P.x, peye->P.y, x0, y0);
				//	printf ("Heee? %f %f\n", dx0, dy0);
				}
			}
		}
	}/**/
	if (n > 0) {
		peye->Ax -= peye->S2Fit_Scale*(peye->Ax - (ax / (float)n));
		peye->Ay -= peye->S2Fit_Scale*(peye->Ay - (ay / (float)n));
	}
}

si	Eye_S3Fit_EdgeGet	(tEye* peye, float t, float* px, float* py)
{
//	si border = peye->Exp_R*peye->Exp_R*3;
	float x = peye->P.x, y = peye->P.y;
	
	for (; ddist2(x,y,peye->P.x,peye->P.y) < peye->Exp_R*peye->Exp_R*3; ) {
		if (	1//dnpix(x,y)->Y == 0xFF
			&& dnpix(x,y)->U == 0xF
			&& dnpix(x,y)->V == 0xF
		) {
			*px = x+0.000000000001f;
			*py = y+0.000000000001f;
			return 1;
		}
		x += cos(t);
		y += sin(t);
	}/**/
	return 0;
}
si	Eye_S3Fit_EdgeMark	(tEye* peye, si i, float t, float* px, float* py)
{
//	si border = peye->Exp_R*peye->Exp_R*3;
	u08 wrote = 0;
	float x = peye->P.x, y = peye->P.y;
	si n = 0;
	
	for (; ddist2(x,y,peye->P.x,peye->P.y) < peye->Exp_R*peye->Exp_R; ) {
		*dspix(videoIn->width+n,i) = dnpix(x,y)->Y | dnpix(x,y)->Y<<8 | dnpix(x,y)->Y<<16;
		if (	dopix(x,y)->Y > ay + 18
			&& ddist2(x,y,peye->P.x,peye->P.y) >= 14*14
	//	if (	1//dnpix(x,y)->Y == 0xFF
	//		&& dnpix(x,y)->U == 0xF
	//		&& dnpix(x,y)->V == 0xF
		) {
		//	return 1;
			float x1 = x, y1 = y;
			for (; ddist2(x1,y1,peye->P.x,peye->P.y) < peye->Exp_R*peye->Exp_R*4; ) {
				if (	dopix(x1,y1)->Y <= ay + 8) {
					break;
				}
				x1 += cos(t);
				y1 += sin(t);
			}
			if (ddist2(x1,y1,x,y) >= 2*2) {
				if (!wrote) {
					wrote = 1;
					*px = x+0.000000000001f;
					*py = y+0.000000000001f;
					*dspix(videoIn->width+n,i) = 0xFF<<8;
				}
				*dspix(videoIn->width+n,i) = 0xFF<<8;
			}
		}/**/
		
		x += cos(t);
		y += sin(t);
		++n;
	}/**/
	return wrote;
	return 0;
}
void	Eye_S3Fit		(tEye* peye)
{
	SDL_LockSurface (gM.pScreen);
//	Eye_SFit (peye);	return;
	si border = 4*peye->Exp_R;
//	float ax = 0, ay = 0;
//	si n = 0;
	si i = 0;
	
	float ao, a;
//	for (a = 0; a < M_PI; a += M_PI/100.0f) {
//	for (a = -M_PI_2; a < M_PI_2; a += M_PI/100.0f) {
//	for (a = -M_PI_4; a < M_PI_4; a += M_PI/100.0f) {
	for (ao = 0; ao < M_PI_2; ao += M_PI/100.0f) {
	for (a = ao; a <= ao+M_PI_2+M_PI_4; a += M_PI_2) {
//	{
	//	printf ("mino ");
		float x0, y0;
		float t0 = a;
	//	u08 skip = Eye_S3Fit_EdgeMark (peye, i, t0, &x0, &y0);
		u08 skip = Eye_S3Fit_EdgeMark (peye, a/(M_PI/100.0f), t0, &x0, &y0);
		
		float t1 = t0 - M_PI;
		float x1, y1;
	//	skip &= Eye_S3Fit_EdgeMark (peye, 110 + i, t1, &x1, &y1);
		skip &= Eye_S3Fit_EdgeMark (peye, 110 + a/(M_PI/100.0f), t1, &x1, &y1);
		++i;
		
	//	printf ("Heee? eye %f %f %f %f\n", peye->P.x, peye->P.y, dx0, dy0);
		
		if (!skip)
			continue;
		
		dnpix(x0,y0)->Y = 0xFF;
		dnpix(x0,y0)->U = 0x0;
		dnpix(x0,y0)->V = 0x0;
		
		dnpix(x1,y1)->Y = 0xFF;
		dnpix(x1,y1)->U = 0x0;
		dnpix(x1,y1)->V = 0x0;
	//	printf ("Heee? eye %ld %ld   %ld %ld\n", x0 - (si)peye->P.x, y0 - (si)peye->P.y, x1 - (si)peye->P.x, y1 - (si)peye->P.y);
	/*	dnpix(x0,y0)->Y = 0x0;
		dnpix(x0,y0)->U = 0x0;
		dnpix(x0,y0)->V = 0x0;
		
		dnpix(x1,y1)->Y = 0xFF;
		dnpix(x1,y1)->U = 0xF;
		dnpix(x1,y1)->V = 0xF;/**/
		if (0) {
			float dx1, dy1;
			dx1 = x1 - peye->P.x;
			dy1 = y1 - peye->P.y;
			
			float ox0 = (peye->Ax * cos(t0) * cos(peye->Aa) - peye->Ay * sin(t0) * sin(peye->Aa));
			float oy0 = (peye->Ax * cos(t0) * sin(peye->Aa) + peye->Ay * sin(t0) * cos(peye->Aa));
			
			float ox1 = (peye->Ax * cos(t1) * cos(peye->Aa) - peye->Ay * sin(t1) * sin(peye->Aa));
			float oy1 = (peye->Ax * cos(t1) * sin(peye->Aa) + peye->Ay * sin(t1) * cos(peye->Aa));
			
		/*	float diff = sqrt(dx0*dx0+dy0*dy0) - sqrt(ox0*ox0+oy0*oy0);
			if (fabs(diff) >= 0.01f) {
			//	printf ("f1\n");
				peye->Ax += peye->Fit_Scale*diff*fabs(cos(t0));
				peye->Ay += peye->Fit_Scale*diff*fabs(sin(t0));
				
				Eye_Xset (peye, peye->P.x + peye->Fit_Trans*(diff )*cos(t0));
				Eye_Yset (peye, peye->P.y + peye->Fit_Trans*(diff )*sin(t0));
			}/**/
		/*	float diff = sqrt(dx1*dx1+dy1*dy1) - sqrt(ox1*ox1+oy1*oy1);
			if (fabs(diff) >= 0.01f) {
			//	printf ("f1\n");
				peye->Ax += peye->Fit_Scale*diff*fabs(cos(t1));
				peye->Ay += peye->Fit_Scale*diff*fabs(sin(t1));
				
				Eye_Xset (peye, peye->P.x + peye->Fit_Trans*(diff )*cos(t1));
				Eye_Yset (peye, peye->P.y + peye->Fit_Trans*(diff )*sin(t1));
			}/**/
		}else {
			float ox0 = (peye->Ax * cos(t0) * cos(peye->Aa) - peye->Ay * sin(t0) * sin(peye->Aa));
			float ox1 = (peye->Ax * cos(t1) * cos(peye->Aa) - peye->Ay * sin(t1) * sin(peye->Aa));
			
			float oy0 = (peye->Ax * cos(t0) * sin(peye->Aa) + peye->Ay * sin(t0) * cos(peye->Aa));
			float oy1 = (peye->Ax * cos(t1) * sin(peye->Aa) + peye->Ay * sin(t1) * cos(peye->Aa));
			
		//	printf ("Heee? pic %f %f   %f %f\n", x0 - peye->P.x, y0 - peye->P.y, x1 - peye->P.x, y1 - peye->P.y);
		//	printf ("Heee? eli %f %f   %f %f\n", ox0, oy0, ox1, oy1);
			
			float sdx = fabsf(x1 - x0) - fabsf(ox1 - ox0);
			float sdy = fabsf(y1 - y0) - fabsf(oy1 - oy0);
			
		//	ax += fabsf((x1 - x0)*cos(t0));
		//	ay += fabsf((y1 - y0)*sin(t0));
		//	++n;
			
			float pdx = (x0+x1)*0.5f - (2.0f*peye->P.x + ox0+ox1)*0.5f;
			float pdy = (y0+y1)*0.5f - (2.0f*peye->P.y + oy0+oy1)*0.5f;
			
			Eye_Xset (peye, peye->P.x + 0.10f * pdx );// *cos(t));
			Eye_Yset (peye, peye->P.y + 0.10f * pdy );// *sin(t));
			
			peye->Ax += 0.10f*sdx*fabs(cos(t0));
			peye->Ay += 0.10f*sdy*fabs(sin(t0));
			
		//	peye->Ax += 0.5f*sdx*fabs(cos(t0));
		//	peye->Ay += 0.5f*sdy*fabs(sin(t0));
			
		//	peye->Ax += peye->Fit_Scale * sdx;
		//	peye->Ay += peye->Fit_Scale * sdy;
			
		//	peye->Ax += peye->Fit_Scale * ( sdx*fabs(cos(t0)) );
		//	peye->Ay += peye->Fit_Scale * ( sdy*fabs(sin(t0)) );
			
		//	peye->Ax += peye->Fit_Scale * ( sdx* (fabs(cos(t0))*fabs(cos(t0))) );
		//	peye->Ay += peye->Fit_Scale * ( sdy* (fabs(sin(t0))*fabs(sin(t0))) );
			
		//	peye->Ax += peye->Fit_Scale * ( sdx* (fabs(cos(t0))+fabs(sin(t0))) );
		//	peye->Ay += peye->Fit_Scale * ( sdy* (fabs(cos(t0))+fabs(sin(t0))) );
			/*
			float pdx = (x0+x1)*0.5f - (ox0+ox1)*0.5f;
			float pdy = (y0+y1)*0.5f - (oy0+oy1)*0.5f;
			
			Eye_Xset (peye, peye->P.x + peye->Fit_Trans * pdx );// *cos(t));
			Eye_Yset (peye, peye->P.y + peye->Fit_Trans * pdy );// *sin(t));
			/**/
			
		//	Eye_Xset (peye, peye->P.x + peye->Fit_Trans * ( (x0+x1)*0.5f - peye->P.x-(ox0+ox1)*0.5f ) );// *cos(t));
		//	Eye_Yset (peye, peye->P.y + peye->Fit_Trans * ( (y0+y1)*0.5f - peye->P.y-(oy0+oy1)*0.5f ) );// *sin(t));
			
		/*	dnpix(ox0, oy0)->Y = 0x0;
			dnpix(ox0, oy0)->U = 0x0;
			dnpix(ox0, oy0)->V = 0x0;
			
			dnpix(ox1, oy1)->Y = 0xFF;
			dnpix(ox1, oy1)->U = 0xF;
			dnpix(ox1, oy1)->V = 0xF;/**/
		}
	}
	}
/*	if (n > 0) {
		peye->Ax -= peye->S2Fit_Scale*(peye->Ax - (ax / (float)n));
		peye->Ay -= peye->S2Fit_Scale*(peye->Ay - (ay / (float)n));
	}/**/
	SDL_UnlockSurface (gM.pScreen);
	return;
}

float	Eye_CFit_AvgDiff	(tEye* peye)
{
	si border = 4*peye->Exp_R;
	float avgdiff = 0, diff_num = 0;
	si x, y;
	for (y = peye->P.y - border; y < peye->P.y + border; ++y) {
		for (x = peye->P.x - border; x < peye->P.x + border; ++x) {
			if (	1//dnpix(x,y)->Y == 0xFF
				&& dnpix(x,y)->U == 0xF
				&& dnpix(x,y)->V == 0xF
			) {
				float t, dx, dy, xx, yy;
				dx = x - peye->P.x;
				dy = y - peye->P.y;
				t = atan2(dy,dx)/* - peye->Aa*/;
				
				xx = (peye->Ax * cos(t) * cos(peye->Aa) - peye->Ay * sin(t) * sin(peye->Aa));
				yy = (peye->Ax * cos(t) * sin(peye->Aa) + peye->Ay * sin(t) * cos(peye->Aa));
				
				float diff = sqrt(dx*dx+dy*dy) - sqrt(xx*xx+yy*yy);
				avgdiff += fabsf(diff);
				++diff_num;
			}
		}
	}/**/
	return avgdiff /= diff_num;
}

si	Eye_CFit_Remove	(tEye* peye, float remdiff)
{
	si border = 4*peye->Exp_R;
	si ret = 0, x, y;
	for (y = peye->P.y - border; y < peye->P.y + border; ++y) {
		for (x = peye->P.x - border; x < peye->P.x + border; ++x) {
			if (	1//dnpix(x,y)->Y == 0xFF
				&& dnpix(x,y)->U == 0xF
				&& dnpix(x,y)->V == 0xF
			) {
				float t, dx, dy, xx, yy;
				dx = x - peye->P.x;
				dy = y - peye->P.y;
				t = atan2(dy,dx)/* - peye->Aa*/;
				
				xx = (peye->Ax * cos(t) * cos(peye->Aa) - peye->Ay * sin(t) * sin(peye->Aa));
				yy = (peye->Ax * cos(t) * sin(peye->Aa) + peye->Ay * sin(t) * cos(peye->Aa));
				
				float diff = sqrt(dx*dx+dy*dy) - sqrt(xx*xx+yy*yy);
				if (fabsf(diff) >= remdiff) {
					dnpix(x,y)->Y = 0x00;
					dnpix(x,y)->U = 0x0;
					dnpix(x,y)->V = 0x0;
					++ret;
				}
			}
		}
	}/**/
	return ret;
}

void	Eye_CFit		(tEye* peye)
{
//	Eye_SFit (peye);	return;
	Eye_S3Fit (peye);	return;
	
	Eye_SFit (peye);
	printf ("%ld  ", Eye_CFit_Remove (peye, Eye_CFit_AvgDiff (peye)*4.0f));
	Eye_SFit (peye);
	printf ("%ld  ", Eye_CFit_Remove (peye, Eye_CFit_AvgDiff (peye)*3.3f));
	Eye_SFit (peye);
	printf ("%ld  ", Eye_CFit_Remove (peye, Eye_CFit_AvgDiff (peye)*2.2f));
	Eye_SFit (peye);
	printf ("%ld  ", Eye_CFit_Remove (peye, Eye_CFit_AvgDiff (peye)*1.2f));
	Eye_SFit (peye);
	
	printf ("\n");
}

void	Eye_Clip		(tEye* peye)
{
	si x, y;
	if (peye->Ax > 100)
		peye->Ax = 100;
	if (peye->Ay > 100)
		peye->Ay = 100;
	if (peye->Ax < 5)
		peye->Ax = 5;
	if (peye->Ay < 5)
		peye->Ay = 5;
}

void	Eye_Draw		(tEye* peye)
{
	si x, y;
	float t;
	for (t = 0; t < 2*M_PI; t += M_PI/180) {
	//	peye->Ax*x*x + peye->Ay*y*y + peye->Axy*x*y = peye->Ar
		x = (peye->Ax * cos(t) * cos(peye->Aa) - peye->Ay * sin(t) * sin(peye->Aa));
		y = (peye->Ax * cos(t) * sin(peye->Aa) + peye->Ay * sin(t) * cos(peye->Aa));
	//	printf ("x %4d  y %4d\n", x, y);
	//	y = peye->P.y;
		if (peye->P.x+x < 0 || peye->P.y+y < 0 || peye->P.x+x >= videoIn->width || peye->P.y+y >= videoIn->height)
			continue;
		dnpix(peye->P.x+x,peye->P.y+y)->Y = 0xFF;
		dnpix(peye->P.x+x,peye->P.y+y)->U = 0xF;
		dnpix(peye->P.x+x,peye->P.y+y)->V = 0xF;
	}/**/
	dnpix(gM.Left.P.x,gM.Left.P.y)->Y = 0xFF;
	dnpix(gM.Left.P.x,gM.Left.P.y)->U = 0xF;
	dnpix(gM.Left.P.x,gM.Left.P.y)->V = 0xF;
}



tV2f	Eye_map_point	(tEye* peye, tV2f p)
{
	return map_point (&peye->Homo, p);
}


void	Head_Init		(tHead* phead)
{
	M3f_Iden (&phead->M);
	M3f_Iden (&phead->MI);
	
}

float	Head_Calc_a3	(tHead* p, float* pa)
{
	pa[0] = atan2(p->DotR.P.y - p->DotL.P.y, p->DotR.P.x - p->DotL.P.x);
	pa[1] = atan2(p->DotR.P.y - p->DotL.P.y, p->DotR.P.x - p->DotL.P.x);
	pa[2] = atan2(p->DotR.P.y - p->DotL.P.y, p->DotR.P.x - p->DotL.P.x);
}

float	Head_Calc_scalex	(tHead* p, float* ps)
{
	tV2f dlr = p->DotR.P;	V2f_sub_V2f (&dlr, &p->DotL.P);
	tV2f dlc = p->DotC.P;	V2f_sub_V2f (&dlc, &p->DotL.P);
	tV2f drc = p->DotC.P;	V2f_sub_V2f (&drc, &p->DotR.P);
	ps[0] = dlr.x;
	ps[1] = dlc.x;
	ps[2] = drc.x;
}
float	Head_Calc_scaley	(tHead* p, float* ps)
{
//	tV2f dlr = p->DotR.P;	V2f_sub_V2f (&dlr, &p->DotL.P);
	tV2f dlc = p->DotC.P;	V2f_sub_V2f (&dlc, &p->DotL.P);
	tV2f drc = p->DotC.P;	V2f_sub_V2f (&drc, &p->DotR.P);
	ps[0] = dlc.y;
	ps[1] = dlc.y;
	ps[2] = drc.y;
}


float	Head_Calc_scale	(tHead* p, float* ps)
{
	tV2f dlr = p->DotR.P;	V2f_sub_V2f (&dlr, &p->DotL.P);
	tV2f dlc = p->DotC.P;	V2f_sub_V2f (&dlc, &p->DotL.P);
	tV2f drc = p->DotC.P;	V2f_sub_V2f (&drc, &p->DotR.P);
	ps[0] = V2f_dist (&dlr);
	ps[1] = V2f_dist (&dlc);
	ps[2] = V2f_dist (&drc);
	return V2f_dist (&dlc) + V2f_dist (&drc) + V2f_dist (&dlr);
}
float	Head_Calc_sheerx	(tHead* p, float* ps)
{
	tV2f dlr = p->DotR.P;	V2f_sub_V2f (&dlr, &p->DotL.P);
	tV2f dlc = p->DotC.P;	V2f_sub_V2f (&dlc, &p->DotL.P);
	tV2f drc = p->DotC.P;	V2f_sub_V2f (&drc, &p->DotR.P);
	
	return V2f_dist (&dlc) / V2f_dist (&drc);
}

float	Head_Calc_sheery	(tHead* p, float* ps)
{
	tV2f dlr = p->DotR.P;	V2f_sub_V2f (&dlr, &p->DotL.P);
	tV2f dlc = p->DotC.P;	V2f_sub_V2f (&dlc, &p->DotL.P);
	tV2f drc = p->DotC.P;	V2f_sub_V2f (&drc, &p->DotR.P);
	
	return (V2f_dist (&dlc) + V2f_dist (&drc)) / V2f_dist (&dlr);
}

float	avg_diff	(int n, float* pv0, float* pv1) {
	float a = 0.0f;
	int i;
	for (i = 0; i < n; ++i)
		a += pv0[i] - pv1[i];
	a /= n;
	return a;
}
float	diff_avg	(int n, float* pv0, float* pv1) {
	float a0 = 0, a1 = 0;
	int i;
	for (i = 0; i < n; ++i) {
		a0 += pv0[i];
		a1 += pv1[i];
	}
	a0 /= n;
	a1 /= n;
	return a0-a1;
}


float	avg_div	(int n, float* pv0, float* pv1) {
	float a = 0.0f;
	int i;
	for (i = 0; i < n; ++i)
		a += pv0[i] / pv1[i];
	a /= n;
	return a;
}


void	Head_Calc_M_Rel	(tHead* p, tHead* pc)
{
	//	L		6.2			R
	//		3.3			4
	//			C
	//	0 0					6.2 0
	//			
	//			2.7 1.9
	int i;
	tV2f dlr = p->DotR.P;	V2f_sub_V2f (&dlr, &p->DotL.P);
	tV2f dlc = p->DotC.P;	V2f_sub_V2f (&dlr, &p->DotL.P);
	tV2f drc = p->DotC.P;	V2f_sub_V2f (&dlr, &p->DotR.P);
	
	float a = 0.0f, aa[3], ac[3];
	Head_Calc_a3 (p, aa);
	Head_Calc_a3 (pc, ac);
	a = avg_diff(3, aa, ac);
	
//	p->x02 = -(p->DotL.P.x + p->DotR.P.x) / 2.0f;
//	p->x12 = -(p->DotL.P.y + p->DotR.P.y) / 2.0f;
	
//	gM.Head = gM.Head;
	tM3f rot;	M3f_Iden (&rot);
	tM3f trans;	M3f_Iden (&trans);
	tM3f scale;	M3f_Iden (&scale);
	tM3f sheer;	M3f_Iden (&sheer);
	tM3f roti = rot;
	tM3f transi = rot;
	tM3f scalei = rot;
	tM3f sheeri = rot;
	
	M3f_Iden (&p->M);
	M3f_Iden (&p->MI);
	
	
	rot.x00 = cos(a);		rot.x01 = -sin(a);
	rot.x10 = sin(a);		rot.x11 = cos(a);
	
	roti.x00 = rot.x00;	roti.x01 = -rot.x01;
	roti.x10 = -rot.x10;	roti.x11 = rot.x11;
	
	
	trans.x02 = (p->DotC.P.x + p->DotL.P.x + p->DotR.P.x) / 3.0f;
	trans.x12 = (p->DotC.P.y + p->DotL.P.y + p->DotR.P.y) / 3.0f;
	
	transi.x02 = -trans.x02;
	transi.x12 = -trans.x12;
	
	
//	M3f_mul_M3f (&p->MI, &roti);
	{
		float s = 1.0f;
		float ssx[3], scx[3];	Head_Calc_scale (p, ssx);	Head_Calc_scale (pc, scx);
		s = avg_div(3, ssx, scx);
	//	s = Head_Calc_scale(p) / Head_Calc_scale(pc);
	//	printf ("scale %f\n", s);
		
		scale.x00 = s;		scale.x01 = 0.0f;
		scale.x10 = 0.0f;		scale.x11 = s;
		
		scalei.x00 = 1.0f/scale.x00;	scalei.x01 = scale.x01;
		scalei.x10 = scale.x10;		scalei.x11 = 1.0f/scale.x11;
	}
/*	{
		float sx = 1.0f, sy = 1.0f;
		float ssx[3], scx[3];	Head_Calc_scalex (p, ssx);	Head_Calc_scalex (pc, scx);
	//	printf ("p  %f %f %f\n", ssx[0], ssx[1], ssx[2]);
	//	printf ("pc %f %f %f\n", scx[0], scx[1], scx[2]);
		sx = avg_div(3, ssx, scx);
	//	printf ("scale x %f\n", sx);
		
		float ssy[3], scy[3];	Head_Calc_scaley (p, ssy);	Head_Calc_scaley (pc, scy);
		sy = avg_div(3, ssy, scy);
		
	//	printf ("scale y %f\n", sy);
		
		scale.x00 = sx;		scale.x01 = 0.0f;
		scale.x10 = 0.0f;		scale.x11 = sy;
		
		scalei.x00 = 1.0f/scale.x00;	scalei.x01 = scale.x01;
		scalei.x10 = scale.x10;		scalei.x11 = 1.0f/scale.x11;
	}/**/
	{
		float sx = 0.0f, sy = 0.0f;
		float ssx[3], scx[3];	Head_Calc_sheerx (p, ssx);	Head_Calc_sheerx (pc, scx);
	//	printf ("p  %f %f %f\n", ssx[0], ssx[1], ssx[2]);
	//	printf ("pc %f %f %f\n", scx[0], scx[1], scx[2]);
		
	//	sx = diff_avg(3, ssx, scx);
	//	sx = Head_Calc_sheerx(p, ssx) - Head_Calc_sheerx(pc, scx);
		sx = Head_Calc_sheerx(pc, scx) - Head_Calc_sheerx(p, ssx);
	//	printf ("sheer x %f\n", sx);
		
	//	sy = Head_Calc_sheery(p, ssx) - Head_Calc_sheery(pc, scx);
		sy = Head_Calc_sheery(pc, scx) - Head_Calc_sheery(p, ssx);
	//	printf ("sheer y %f\n", sy);
		
	//	float ssy[3], scy[3];	Head_Calc_scaley (p, ssy);	Head_Calc_scaley (pc, scy);
	//	sy = avg_div(3, ssy, scy);
		
		sheeri.x00 = 1.0f;		sheeri.x01 = sx;
		sheeri.x10 = -sy;			sheeri.x11 = 1.0f;
		
	//	printf ("sheer x %f\n", sheeri.x01);
		
	}/**/
/*	{
		float sx = 0.0f, ssx[3], scx[3];	Head_Calc_scalex (p, ssx);	Head_Calc_scalex (pc, scx);
		sx = avg_div(3, scx, ssx);
		sx = 1.0f - sx;
		printf ("sheer x %f\n", sx);
		
		sheeri.x00 = 1.0f;		sheeri.x01 = sx;
		sheeri.x10 = 0.0f;		sheeri.x11 = 1.0f;
	}/**/
	M3f_mul_M3f (&p->M, &trans);
	M3f_mul_M3f (&p->M, &rot);
	M3f_mul_M3f (&p->M, &scale);
	
//	M3f_mul_M3f (&p->MI, &transi);
	M3f_mul_M3f (&p->MI, &roti);
	M3f_mul_M3f (&p->MI, &scalei);
	M3f_mul_M3f (&p->MI, &sheeri);
}

void	muhaha_Init	()
{
	Head_Init (&gM.Head);
	
	
}

void	muhaha	()
{
	tPix* pin = (tPix*)videoIn->framebuffer, *pin1;
	si x, y;
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
	/*	for (y = 0; y < videoIn->height; ++y) {
			for (x = 0; x < videoIn->width/2; ++x) {
				tPix tmp = *dopix(x,y);
				*dopix(x,y) = *dopix(videoIn->width-x-1,y);
				*dopix(videoIn->width-x-1,y) = tmp;
			}
		}/**/
		for (y = 0; y < videoIn->height; ++y) {
			for (x = 0; x < videoIn->width-1; ++x) {
				pin = (tPix*)videoIn->framebuffer + (y*videoIn->width + x);
				pin1 = pin+1;
			//	pin->Y = 0;
				si y = pin->Y, y1 = pin1->Y;
			//	pin->Y = abs(y1 - y);
				
			//	pin->U = ((si)pin1->U - ((si)pin->U)) +7;
			//	pin->V = ((si)pin1->V - ((si)pin->V)) +8;
			//	pin->V = abs((si)pin1->V - (si)pin->V);
				pin->U = 0x7;
				pin->V = 0x8;
			//	dpix(x,y)->U = 0x7;
			//	dpix(x,y)->V = 0x8;
			}
		}/**/
	}
	
	
	memcpy(gM.pDst, videoIn->framebuffer, videoIn->width * (videoIn->height) * 2);
	
//	Eye_CalcAYUV (&gM.Head.DotC);
//	Eye_S3Fit (&gM.Head.DotC);
	Eye_FF (&gM.Head.DotC);
	Eye_EdgeMark (&gM.Head.DotC);
	Eye_SFit (&gM.Head.DotC);
	Eye_Clip (&gM.Head.DotC);
	
	memcpy(gM.pDst, videoIn->framebuffer, videoIn->width * (videoIn->height) * 2);
	
//	Eye_CalcAYUV (&gM.Head.DotL);
//	Eye_S3Fit (&gM.Head.DotL);
	Eye_FF (&gM.Head.DotL);
	Eye_EdgeMark (&gM.Head.DotL);
	Eye_SFit (&gM.Head.DotL);
	Eye_Clip (&gM.Head.DotL);
	
	memcpy(gM.pDst, videoIn->framebuffer, videoIn->width * (videoIn->height) * 2);
	
//	Eye_CalcAYUV (&gM.Head.DotR);
//	Eye_S3Fit (&gM.Head.DotR);
	Eye_FF (&gM.Head.DotR);
	Eye_EdgeMark (&gM.Head.DotR);
	Eye_SFit (&gM.Head.DotR);
	Eye_Clip (&gM.Head.DotR);
	
	memcpy(gM.pDst, videoIn->framebuffer, videoIn->width * (videoIn->height) * 2);
	
	Eye_CalcAYUV (&gM.Right);
	Eye_S3Fit (&gM.Right);
//	Eye_FF (&gM.Right);
//	Eye_EdgeMark (&gM.Right);
//	Eye_SFit (&gM.Right);
	
	Eye_Clip (&gM.Right);
	
	
	//hehehe evil
	memcpy(gM.pDst, videoIn->framebuffer, videoIn->width * (videoIn->height) * 2);
	
	Eye_CalcAYUV (&gM.Left);
	Eye_S3Fit (&gM.Left);
//	Eye_FF (&gM.Left);
//	Eye_EdgeMark (&gM.Left);
//	Eye_CFit (&gM.Left);
	Eye_Clip (&gM.Left);
	
	
	
	Eye_Draw (&gM.Left);
	Eye_Draw (&gM.Right);
	
	Eye_Draw (&gM.Head.DotC);
	Eye_Draw (&gM.Head.DotL);
	Eye_Draw (&gM.Head.DotR);
	
	gCol.Y = 0xFF;
	gCol.U = 0x0;
	gCol.V = 0x0;
//	Vec_Draw (gM.Head.DotL.P.x, gM.Head.DotL.P.y, gM.Head.DotR.P.x, gM.Head.DotR.P.y);
	
	gCol.U = 0xF;
	gCol.V = 0xF;
//	Vec_Draw (gM.Head.DotL.P.x, gM.Head.DotL.P.y, gM.Head.DotL.P.x+gM.L_Vec.x, gM.Head.DotL.P.y+gM.L_Vec.y);
//	Vec_Draw (gM.Head.DotR.P.x, gM.Head.DotR.P.y, gM.Head.DotR.P.x+gM.R_Vec.x, gM.Head.DotR.P.y+gM.R_Vec.y);
	
	Vec_Draw ((640/4*1)-gM.L_Vec.x, (480/2)-gM.L_Vec.y, (640/4*1), (480/2));
	Vec_Draw ((640/4*3)-gM.R_Vec.x, (480/2)-gM.R_Vec.y, (640/4*3), (480/2));
	
	
	gCol.Y = 0xFF;
	gCol.U = 0x0;
	gCol.V = 0x0;
	
	tV2f pc = {0, 0};
	tV2f lx = {40, 0};
	tV2f ly = {0, 20};
	
	V2f_mul_M3f (&pc, &gM.Head.M);
	V2f_mul_M3f (&lx, &gM.Head.M);
	V2f_mul_M3f (&ly, &gM.Head.M);
	
	V2f_DrawPosPos (&pc, &lx);
	V2f_DrawPosPos (&pc, &ly);
	
//	if (videoIn->formatIn != V4L2_PIX_FMT_YUYV)
//		printf ("videoIn->formatIn %d\n", videoIn->formatIn);
	
//	memcpy(gM.pDst, videoIn->framebuffer, videoIn->width * (videoIn->height) * 2);
	
	#define dd 1
	for (y = gM.Gaze.y-dd; y <= gM.Gaze.y+dd; ++y) {
		for (x = gM.Gaze.x-dd; x <= gM.Gaze.x+dd; ++x) {
			if (dpixout(x,y))
				continue;
			dnpix(x,y)->Y = 0xFF;
			dnpix(x,y)->U = 0xF;
			dnpix(x,y)->V = 0xF;
		}
	}
	#undef dd
	#define dd 4
	x = gM.Gaze.x;
	for (y = gM.Gaze.y-dd; y <= gM.Gaze.y+dd; ++y) {
		if (dpixout(x,y))
			continue;
		dnpix(x,y)->Y = 0xFF;
		dnpix(x,y)->U = 0xF;
		dnpix(x,y)->V = 0xF;
	}
	y = gM.Gaze.y;
	for (x = gM.Gaze.x-dd; x <= gM.Gaze.x+dd; ++x) {
		if (dpixout(x,y))
			continue;
		dnpix(x,y)->Y = 0xFF;
		dnpix(x,y)->U = 0xF;
		dnpix(x,y)->V = 0xF;
	}
	#undef dd
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
	
	static si ii = 0;
	
//	action_gui curr_action = A_VIDEO;
	while (videoIn->signalquit) {
		SDL_LockMutex(affmutex);
		
		float frmrate = gdata->frmrate;
		
		while (SDL_PollEvent(sdlevent)) {	//scan the event queue
			switch (sdlevent->type) {
			case SDL_KEYDOWN:
				printf("Key down\n");
				switch (sdlevent->key.keysym.sym) {
				case SDLK_BACKSPACE:
					printf("Reset points\n");
					ii = 0;
					break;
				case SDLK_s:
					printf("gM.DeSat\n");
					gM.DeSat = !gM.DeSat;
					break;
				case SDLK_o:
					gM.Head.DotC.P.x = x;
					gM.Head.DotC.P.y = y;
					break;
				case SDLK_a:
					gM.Head.DotL.P.x = x;
					gM.Head.DotL.P.y = y;
					break;
				case SDLK_e:
					gM.Head.DotR.P.x = x;
					gM.Head.DotR.P.y = y;
					break;
				case SDLK_SEMICOLON:
					gM.Left.P.x = x;
					gM.Left.P.y = y;
					break;
				case SDLK_j:
					gM.Right.P.x = x;
					gM.Right.P.y = y;
					break;
				case SDLK_SPACE:
					gM.HeadC = gM.Head;
					break;
				}
			//	curr_action = GUI_keytoaction(sdlevent->key.keysym.sym);
			//	if (curr_action != A_VIDEO)
			//		mouseon = 1;
				break;
			case SDL_KEYUP:
			case SDL_MOUSEBUTTONUP:
				mouseon = 0;
				boucle = 0;
				break;
			case SDL_MOUSEBUTTONDOWN:
				printf ("mouse down\n");
				SDL_GetMouseState(&x, &y);
				switch (sdlevent->button.button) {
				case SDL_BUTTON_LEFT: {
					gM.Left.Homo.scenecalipoints[ii].x = x;
					gM.Left.Homo.scenecalipoints[ii].y = y;
					gM.Right.Homo.scenecalipoints[ii].x = x;
					gM.Right.Homo.scenecalipoints[ii].y = y;
					
				/*	gM.Left.Homo.vectors[ii].x = gM.Left.P.x;
					gM.Left.Homo.vectors[ii].y = gM.Left.P.y;
					
					gM.Right.Homo.vectors[ii].x = gM.Right.P.x;
					gM.Right.Homo.vectors[ii].y = gM.Right.P.y;/**/
					
				/*	gM.Left.Homo.vectors[ii].x = gM.Left.P.x - gM.Head.DotL.P.x;
					gM.Left.Homo.vectors[ii].y = gM.Left.P.y - gM.Head.DotL.P.y;
					
					gM.Right.Homo.vectors[ii].x = gM.Right.P.x - gM.Head.DotR.P.x;
					gM.Right.Homo.vectors[ii].y = gM.Right.P.y - gM.Head.DotR.P.y;/**/
					
					gM.Left.Homo.vectors[ii].x = gM.L_Vec.x;
					gM.Left.Homo.vectors[ii].y = gM.L_Vec.y;
					
					gM.Right.Homo.vectors[ii].x = gM.R_Vec.x;
					gM.Right.Homo.vectors[ii].y = gM.R_Vec.y;/**/
					printf ("got point %ld\n", ii);
					
					if (++ii >= CALIBRATIONPOINTS) {
						ii = 0;
						CalculateCalibration(&gM.Left.Homo);
						CalculateCalibration(&gM.Right.Homo);
					}
					break;
				}
				case SDL_BUTTON_RIGHT: {
					gM.Left.P.x = x;
					gM.Left.P.y = y;
					break;
				}
				}
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
		/*
		{
			tV2f pl, pr;
			tM2f rot;
			{	float a = atan2(gM.Head.DotR.P.y - gM.Head.DotL.P.y, gM.Head.DotR.P.x - gM.Head.DotL.P.x);
				rot.a = cos(a);	rot.b = sin(a);
				rot.c = -sin(a);	rot.d = cos(a);
			}
			
		//	pl.x = gM.Left.P.x;	pl.y = gM.Left.P.y;
		//	pr.x = gM.Right.P.x;	pr.y = gM.Right.P.y;
			
			pl.x = gM.Left.P.x - gM.Head.DotL.P.x;	pl.y = gM.Left.P.y - gM.Head.DotL.P.y;
			pr.x = gM.Right.P.x - gM.Head.DotR.P.x;	pr.y = gM.Right.P.y - gM.Head.DotR.P.y;
			
			V2f_mul_M2f (&pl, &rot);
			V2f_mul_M2f (&pr, &rot);
			
			gM.L_Vec = pl;
			gM.R_Vec = pr;
			
		//	printf ("eye %f %f", p.x, p.y);
			pl = Eye_map_point(&gM.Left, pl);
			pr = Eye_map_point(&gM.Right, pr);
			point_clip (&pl);
			point_clip (&pr);
			
			
			gM.Gaze.x = pl.x;
			gM.Gaze.y = pl.y;
			
			gM.Gaze.x = pr.x;
			gM.Gaze.y = pr.y;
			
			gM.Gaze.x = (pl.x+pr.x)/2.0f;
			gM.Gaze.y = (pl.y+pr.y)/2.0f;
		}/**/
		
	/*	{	float a = atan2(gM.Head.DotR.P.y - gM.Head.DotL.P.y, gM.Head.DotR.P.x - gM.Head.DotL.P.x);
			gM.Head.x00 = cos(a);		gM.Head.x01 = -sin(a);
			gM.Head.x10 = sin(a);		gM.Head.x11 = cos(a);
			
			gM.Head.x02 = (gM.Head.DotL.P.x + gM.Head.DotR.P.x) / 2.0f;
			gM.Head.x12 = (gM.Head.DotL.P.y + gM.Head.DotR.P.y) / 2.0f;
		}/**/
		
		{	
			Head_Calc_M_Rel (&gM.Head, &gM.HeadC);
			
			tV2f pel, per;
			
			pel = gM.Left.P;		V2f_sub_V2f (&pel, &gM.Head.DotL.P);
			per = gM.Right.P;		V2f_sub_V2f (&per, &gM.Head.DotR.P);
			
		//	pel = gM.Head.DotL.P;	V2f_sub_V2f (&pel, &gM.Left.P);
		//	per = gM.Head.DotR.P;	V2f_sub_V2f (&per, &gM.Right.P);
			
			
		/*	pel = gM.Left.P;
			per = gM.Right.P;
			pel.x -= gM.Head.M.x02;	pel.y -= gM.Head.M.x12;
			per.x -= gM.Head.M.x02;	per.y -= gM.Head.M.x12;
		//	pel.x = gM.Left.P.x - gM.Head.DotL.P.x;	pel.y = gM.Left.P.y - gM.Head.DotL.P.y;
		//	per.x = gM.Right.P.x - gM.Head.DotR.P.x;	per.y = gM.Right.P.y - gM.Head.DotR.P.y;
		
		//	tV2f delr = gM.Right.P;	V2f_sub_V2f (&delr, &gM.Left.P);
		/*	{	float x, y;
				x = (gM.Left.P.x + gM.Right.P.x) / 2.0f;
				y = (gM.Left.P.y + gM.Right.P.y) / 2.0f;
				pel.x -= x;		pel.y -= y;
				per.x -= x;		per.y -= y;
			}*/
			
			V2f_mul_M3f (&pel, &gM.Head.MI);
			V2f_mul_M3f (&per, &gM.Head.MI);
			
			gM.L_Vec = pel;
			gM.R_Vec = per;
			
			pel = Eye_map_point(&gM.Left, pel);
			per = Eye_map_point(&gM.Right, per);
			
			gM.Gaze.x = (pel.x+per.x)/2.0f;
			gM.Gaze.y = (pel.y+per.y)/2.0f;
		}
		
		point_clip (&gM.Gaze);
		
	//	printf (" to %f %f\n", p.x, p.y);
	//	printf("fp/s %d\n",frmrate);
	}
	
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




tV2d* normalize_point_set(tV2d* point_set, double *dis_scale, tV2d *nor_center, int num)
{
	double sumx = 0, sumy = 0;
	double sumdis = 0;
	tV2d *edge = point_set;
	int i;
	for (i = 0; i < num; i++) {
		sumx += edge->x;
		sumy += edge->y;
		sumdis += sqrt((double)(edge->x*edge->x + edge->y*edge->y));
		edge++;
	}
	
	*dis_scale = sqrt((double)2)*num/sumdis;
	nor_center->x = sumx*1.0/num;
	nor_center->y = sumy*1.0/num;
	tV2d *edge_point_nor = (tV2d*)malloc(sizeof(tV2d)*num);
	edge = point_set;
	for (i = 0; i < num; i++) {
		edge_point_nor[i].x = (edge->x - nor_center->x)*(*dis_scale);
		edge_point_nor[i].y = (edge->y - nor_center->y)*(*dis_scale);
		edge++;
	}
	return edge_point_nor;
}


#define SIGN(u, v)     ( (v)>=0.0 ? fabs(u) : -fabs(u) )
#define MAX(x, y)     ( (x) >= (y) ? (x) : (y) )  

static double   radius(double u, double v)
{
	double          w;
	u = fabs(u);
	v = fabs(v);
	if (u > v) {
		w = v / u;
		return (u * sqrt(1. + w * w));
    } else {
        if (v) {
	         w = u / v;
             return (v * sqrt(1. + w * w));
		} else
			return 0.0;
	}
}

/*
 Given matrix a[m][n], m>=n, using svd decomposition a = p d q' to get
 p[m][n], diag d[n] and q[n][n].
*/
void svd(int m, int n, double **a, double **p, double *d, double **q)
{
        int             flag, i, its, j, jj, k, l, nm, nm1 = n - 1, mm1 = m - 1;
        double          c, f, h, s, x, y, z;
        double          anorm = 0, g = 0, scale = 0;
        //double         *r = tvector_alloc(0, n, double);
		double			*r = (double*)malloc(sizeof(double)*n);

        for (i = 0; i < m; i++)
                for (j = 0; j < n; j++)
                        p[i][j] = a[i][j];
        //for (i = m; i < n; i++)
        //                p[i][j] = 0;

        /* Householder reduction to bidigonal form */
        for (i = 0; i < n; i++)
        {
                l = i + 1;
                r[i] = scale * g;
                g = s = scale = 0.0;
                if (i < m)
                {
                        for (k = i; k < m; k++)
                                scale += fabs(p[k][i]);
                        if (scale)
                        {
                                for (k = i; k < m; k++)
                                {
                                        p[k][i] /= scale;
                                        s += p[k][i] * p[k][i];
                                }
                                f = p[i][i];
                                g = -SIGN(sqrt(s), f);
                                h = f * g - s;
                                p[i][i] = f - g;
                                if (i != nm1)
                                {
                                        for (j = l; j < n; j++)
                                        {
                                                for (s = 0.0, k = i; k < m; k++)
                                                        s += p[k][i] * p[k][j];
                                                f = s / h;
                                                for (k = i; k < m; k++)
                                                        p[k][j] += f * p[k][i];
                                        }
                                }
                                for (k = i; k < m; k++)
                                        p[k][i] *= scale;
                        }
                }
                d[i] = scale * g;
                g = s = scale = 0.0;
                if (i < m && i != nm1)
                {
                        for (k = l; k < n; k++)
                                scale += fabs(p[i][k]);
                        if (scale)
                        {
                                for (k = l; k < n; k++)
                                {
                                        p[i][k] /= scale;
                                        s += p[i][k] * p[i][k];
                                }
                                f = p[i][l];
                                g = -SIGN(sqrt(s), f);
                                h = f * g - s;
                                p[i][l] = f - g;
                                for (k = l; k < n; k++)
                                        r[k] = p[i][k] / h;
                                if (i != mm1)
                                {
                                        for (j = l; j < m; j++)
                                        {
                                                for (s = 0.0, k = l; k < n; k++)
                                                        s += p[j][k] * p[i][k];
                                                for (k = l; k < n; k++)
                                                        p[j][k] += s * r[k];
                                        }
                                }
                                for (k = l; k < n; k++)
                                        p[i][k] *= scale;
                        }
                }
                anorm = MAX(anorm, fabs(d[i]) + fabs(r[i]));
        }

        /* Accumulation of right-hand transformations */
        for (i = n - 1; i >= 0; i--)
        {
                if (i < nm1)
                {
                        if (g)
                        {
                                for (j = l; j < n; j++)
                                        q[j][i] = (p[i][j] / p[i][l]) / g;
                                for (j = l; j < n; j++)
                                {
                                        for (s = 0.0, k = l; k < n; k++)
                                                s += p[i][k] * q[k][j];
                                        for (k = l; k < n; k++)
                                                q[k][j] += s * q[k][i];
                                }
                        }
                        for (j = l; j < n; j++)
                                q[i][j] = q[j][i] = 0.0;
                }
                q[i][i] = 1.0;
                g = r[i];
                l = i;
        }
        /* Accumulation of left-hand transformations */
        for (i = n - 1; i >= 0; i--)
        {
                l = i + 1;
                g = d[i];
                if (i < nm1)
                        for (j = l; j < n; j++)
                                p[i][j] = 0.0;
                if (g)
                {
                        g = 1.0 / g;
                        if (i != nm1)
                        {
                                for (j = l; j < n; j++)
                                {
                                        for (s = 0.0, k = l; k < m; k++)
                                                s += p[k][i] * p[k][j];
                                        f = (s / p[i][i]) * g;
                                        for (k = i; k < m; k++)
                                                p[k][j] += f * p[k][i];
                                }
                        }
                        for (j = i; j < m; j++)
                                p[j][i] *= g;
                } else
                        for (j = i; j < m; j++)
                                p[j][i] = 0.0;
                ++p[i][i];
        }
        /* diagonalization of the bidigonal form */
        for (k = n - 1; k >= 0; k--)
        {                       /* loop over singlar values */
                for (its = 0; its < 30; its++)
                {               /* loop over allowed iterations */
                        flag = 1;
                        for (l = k; l >= 0; l--)
                        {       /* test for splitting */
                                nm = l - 1;     /* note that r[l] is always
                                                 * zero */
                                if (fabs(r[l]) + anorm == anorm)
                                {
                                        flag = 0;
                                        break;
                                }
                                if (fabs(d[nm]) + anorm == anorm)
                                        break;
                        }
                        if (flag)
                        {
                                c = 0.0;        /* cancellation of r[l], if
                                                 * l>1 */
                                s = 1.0;
                                for (i = l; i <= k; i++)
                                {
                                        f = s * r[i];
                                        if (fabs(f) + anorm != anorm)
                                        {
                                                g = d[i];
                                                h = radius(f, g);
                                                d[i] = h;
                                                h = 1.0 / h;
                                                c = g * h;
                                                s = (-f * h);
                                                for (j = 0; j < m; j++)
                                                {
                                                        y = p[j][nm];
                                                        z = p[j][i];
                                                        p[j][nm] = y * c + z * s;
                                                        p[j][i] = z * c - y * s;
                                                }
                                        }
                                }
                        }
                        z = d[k];
                        if (l == k)
                        {       /* convergence */
                                if (z < 0.0)
                                {
                                        d[k] = -z;
                                        for (j = 0; j < n; j++)
                                                q[j][k] = (-q[j][k]);
                                }
                                break;
                        }
                        if (its == 30)
                        {
                                //error("svd: No convergence in 30 svd iterations", non_fatal);
                                return;
                        }
                        x = d[l];       /* shift from bottom 2-by-2 minor */
                        nm = k - 1;
                        y = d[nm];
                        g = r[nm];
                        h = r[k];
                        f = ((y - z) * (y + z) + (g - h) * (g + h)) / (2.0 * h * y);
                        g = radius(f, 1.0);
                        /* next QR transformation */
                        f = ((x - z) * (x + z) + h * ((y / (f + SIGN(g, f))) - h)) / x;
                        c = s = 1.0;
                        for (j = l; j <= nm; j++)
                        {
                                i = j + 1;
                                g = r[i];
                                y = d[i];
                                h = s * g;
                                g = c * g;
                                z = radius(f, h);
                                r[j] = z;
                                c = f / z;
                                s = h / z;
                                f = x * c + g * s;
                                g = g * c - x * s;
                                h = y * s;
                                y = y * c;
                                for (jj = 0; jj < n; jj++)
                                {
                                        x = q[jj][j];
                                        z = q[jj][i];
                                        q[jj][j] = x * c + z * s;
                                        q[jj][i] = z * c - x * s;
                                }
                                z = radius(f, h);
                                d[j] = z;       /* rotation can be arbitrary
                                                 * id z=0 */
                                if (z)
                                {
                                        z = 1.0 / z;
                                        c = f * z;
                                        s = h * z;
                                }
                                f = (c * g) + (s * y);
                                x = (c * y) - (s * g);
                                for (jj = 0; jj < m; jj++)
                                {
                                        y = p[jj][j];
                                        z = p[jj][i];
                                        p[jj][j] = y * c + z * s;
                                        p[jj][i] = z * c - y * s;
                                }
                        }
                        r[l] = 0.0;
                        r[k] = f;
                        d[k] = x;
                }
        }
        free(r);

		// dhli add: the original code does not sort the eigen value
		// should do that and change the eigen vector accordingly

}




void	Ehh_Draw_Line_2d	(si x0, si y0, si x1, si y1)
{
	si dx = x1 - x0;
	si dy = y1 - y0;
	
	ui dst_w = videoIn->width;
	ui dst_h = videoIn->height;
	tPix *pdst = (tPix*)gM.pDst;
	
	tPix col = gCol;
	if (dx < 0) {
		dx = -dx;
		x0 = x1;
		x1 = x0 + dx;
		dy = -dy;
		y0 = y1;
		y1 = y0 + dy;
	}
	si x, y, num;
	if (abs(dx) >= abs(dy)) {
		if (dy >= 0) {
			y = y0;
			num = dx / 2;
			for (x = x0; x <= x1; x++) {
				if ( (x >= 0) && (x < dst_w) && (y >= 0) && (y < dst_h) ) //just so it works for now
					*(pdst + x + (y * dst_w)) = col;
				num += dy;
				if (num >= dx) {
					num -= dx;
					y++;
				}
			}
		}else {
			y = y0;
			num = dx / 2;
			for (x = x0; x <= x1; x++) {
				if ( (x >= 0) && (x < dst_w) && (y >= 0) && (y < dst_h) ) //just so it works for now
					*(pdst + x + (y * dst_w)) = col;
				num += -dy;
				if (num >= dx) {
					num -= dx;
					y--;
				}
			}
		}
	}else {
		if (dy >= 0) {
			x = x0;
			num = dy / 2;
			for (y = y0; y <= y1; y++) {
				if ( (x >= 0) && (x < dst_w) && (y >= 0) && (y < dst_h) ) //just so it works for now
					*(pdst + x + (y * dst_w)) = col;
				num += dx;
				if (num >= dy) {
					num -= dy;
					x++;
				}
			}
		}else {
			x = x0;
			num = -dy / 2;
			for (y = y0; y >= y1; y--) {
				if ( (x >= 0) && (x < dst_w) && (y >= 0) && (y < dst_h) ) //just so it works for now
					*(pdst + x + (y * dst_w)) = col;
				num += dx;
				if (num >= -dy) {
					num -= -dy;
					x++;
				}
			}
		}
	}
	return;
}
