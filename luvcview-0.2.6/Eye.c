
#include "muhaha.h"

#include "priv_macros.h"


#define dcam	peye->aCam[pcam->Idx]


void	Homo_PreCalPrint	(tHomo* phomo)
{
	si ix, iy, i;
	
	printf ("ScreenPoints: \n");
	for (iy = 0; iy < dEye_Screen_Cal_NUM; ++iy) {
		for (ix = 0; ix < dEye_Screen_Cal_NUM; ++ix) {
			i = ix + iy * dEye_Screen_Cal_NUM;
			printf ("%f\t%f\t\t| ", phomo->scenecalipoints[i].x, phomo->scenecalipoints[i].y);
		}
		printf ("\n");
	}
	printf ("Vectors: \n");
	for (iy = 0; iy < dEye_Screen_Cal_NUM; ++iy) {
		for (ix = 0; ix < dEye_Screen_Cal_NUM; ++ix) {
			i = ix + iy * dEye_Screen_Cal_NUM;
			printf ("%f\t%f\t\t| ", phomo->vectors[i].x, phomo->vectors[i].y);
		}
		printf ("\n");
	}
	
	printf ("\n");
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
	Homo_PreCalPrint (phomo);
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
	Homo_PreCalPrint (phomo);
  int i, j;
  tV2d cal_scene[CALIBRATIONPOINTS], cal_eye[CALIBRATIONPOINTS];
  tV2d scene_center, eye_center, *eye_nor, *scene_nor;
  double dis_scale_scene, dis_scale_eye;  

  for (i = 0; i < CALIBRATIONPOINTS; i++) {
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



void	Eye_Init	(tEye* peye)
{
	printf ("INIT  %lx\n", peye);
	peye->Point_N = 0;
	peye->Point_Max = 0;
	peye->paPoint = 0;
	
	peye->InHead.Line_N = 0;
	memset (peye->InHead.aLine, 0, sizeof(peye->InHead.aLine));
	
	for (int i = 0; i < gM.Cam_N; ++i) {
		tCam* pcam = &gM.aCam[i];
		dcam.FF.paMark = 0;
	}
	
	peye->GV_mutex = SDL_CreateMutex();
	
	if (!peye->GV_mutex){
		fprintf(stderr, "Couldn't create mutex\n");
		exit(-1);
	}
/*	if(SDL_mutexP(peye->GV_mutex)==-1){
		fprintf(stderr, "Couldn't lock mutex\n");
		exit(-1);
	}
	/**/
}

void	Eye_Conf	(tEye* peye)
{
	printf ("CONF  %lx\n", peye);
	for (int i = 0; i < gM.Cam_N; ++i) {
		tCam* pcam = &gM.aCam[i];
		dcam.FF.paMark = realloc (dcam.FF.paMark, dpow2(dcam.FF.Max_R) * sizeof(dcam.FF.paMark[0]));
	}
	
	
	
}

void	Eye_CopyParam	(tEye* pdst, tEye* psrc)
{
	pdst->aCam[pcam->Idx].P = psrc->aCam[pcam->Idx].P;
	pdst->Ax = psrc->Ax;
	pdst->Ay = psrc->Ay;
	pdst->Aa = psrc->Aa;
}


float	Eye_Ellipse_xydist2		(tEye* peye, float t, float* px, float* py)
{
	float x, y;
	x = (peye->Ax * cos(t) * cos(peye->Aa) - peye->Ay * sin(t) * sin(peye->Aa));
	y = (peye->Ax * cos(t) * sin(peye->Aa) + peye->Ay * sin(t) * cos(peye->Aa));
	
	if (px) {
		*px = x;
		*py = y;
	}
	return ddist2(0,0, x,y);
}


void	Eye_Points_InsOrd		(tEye* peye, float a, tV2f* pp)
{
	if (peye->Point_N >= peye->Point_Max) {
		printf ("shit not enough points\n");
		return;
	}
	si i;
	for (i = 0; i < peye->Point_N; ++i) {
		float x1 = peye->paPoint[i].x, y1 = peye->paPoint[i].y;
		float t1 = atan2(y1 - dcam.P.y, x1 - dcam.P.x);
		
		if (a < t1) {
			memmove (
				peye->paPoint + i+1,
				peye->paPoint + i,
				(peye->Point_N - i) * sizeof(peye->paPoint[0])
			);
			break;
		}
	}
	peye->paPoint[i] = *pp;
	peye->Point_N++;
}
void	Eye_Points_Ins		(tEye* peye, float a, tV2f* pp)
{
	if (peye->Point_N >= peye->Point_Max) {
		printf ("shit not enough points\n");
		return;
	}
	peye->paPoint[peye->Point_N] = *pp;
	peye->Point_N++;
}

void	Eye_Points_Sort		(tEye* peye)
{
	si i;
	for(i = 1; i < peye->Point_N; ++i)
	{
		float x0 = peye->paPoint[i].x, y0 = peye->paPoint[i].y;
		float t0 = atan2(y0 - dcam.P.y, x0 - dcam.P.x);
		
		si j = i - 1;
		while (1) {
			float x1 = peye->paPoint[j].x, y1 = peye->paPoint[j].y;
			float t1 = atan2(y1 - dcam.P.y, x1 - dcam.P.x);
			
			if (t0 >= t1 || j < 0)
				break;
			
			peye->paPoint[j+1] = peye->paPoint[j];
			--j;
		}
		peye->paPoint[j+1].x = x0;
		peye->paPoint[j+1].y = y0;
	}
}

si	Eye_Points_FindRot	(tEye* peye, si idx, float ang)
{
	#if 0
	float x0 = peye->paPoint[idx].x, y0 = peye->paPoint[idx].y;
	float t0 = atan2(y0 - dcam.P.y, x0 - dcam.P.x);
	
	struct {
		si idx;
		float diff_t;
	}min = {-1, 4*deg2rad};
//	printf ("Eye_S4_Fit_FindRot\n");
	si i;
	for (i = 0; i < peye->Point_N; ++i) {
		float x1 = peye->paPoint[i].x, y1 = peye->paPoint[i].y;
		float t1 = atan2(y1 - dcam.P.y, x1 - dcam.P.x);
	//	printf ("%ld\t%f\n", i, t1);
		
		float diff_t = fabsf(angle_diff_norm_pi_pi(t1, t0) - ang);
		if (diff_t < min.diff_t) {
			min.idx = i;
			min.diff_t = diff_t;
		}
	}
	return min.idx;
	#else
	si i;
	float x0 = peye->paPoint[idx].x, y0 = peye->paPoint[idx].y;
	float t0 = atan2(y0 - dcam.P.y, x0 - dcam.P.x);
	float tgt = angle_norm_pi_pi(t0 + ang);
	
	struct {
		si idx;
		float diff_t;
	}min = {-1, 4*deg2rad};
	
//	printf ("Eye_Points_FindRot tgt %f\n", tgt);
	
	float x1, y1, t1;
	si l = 0, r = peye->Point_N-1;
	while (r-l > 1) {
		i = (l+r) / 2;
	//	printf ("l %ld\ti %ld\tr %ld\n", l, i, r);
		
		x1 = peye->paPoint[i].x;
		y1 = peye->paPoint[i].y;
		t1 = atan2(y1 - dcam.P.y, x1 - dcam.P.x);
		
		if (tgt > t1) {
			l = i;
		}else {
			r = i;
		}
	}
	i = l;
	x1 = peye->paPoint[i].x, y1 = peye->paPoint[i].y;
	t1 = atan2(y1 - dcam.P.y, x1 - dcam.P.x);
//	printf ("got l %ld\t%f\n", l, t1);
	float diff_t = fabsf(angle_diff_norm_pi_pi(t1, t0) - ang);
	if (diff_t < min.diff_t) {
		min.idx = i;
		min.diff_t = diff_t;
	}
	i = r;
	x1 = peye->paPoint[i].x, y1 = peye->paPoint[i].y;
	t1 = atan2(y1 - dcam.P.y, x1 - dcam.P.x);
//	printf ("got r %ld\t%f\n", r, t1);
	diff_t = fabsf(angle_diff_norm_pi_pi(t1, t0) - ang);
	if (diff_t < min.diff_t) {
		min.idx = i;
		min.diff_t = diff_t;
	}
//	printf ("finish %ld\n", min.idx);
	return min.idx;
	#endif
}

float	Eye_Points_Fit2		(tEye* peye)
{
	si i;
	for (i = 0; i < peye->Point_N; ++i) {
		float x0 = peye->paPoint[i].x, y0 = peye->paPoint[i].y;
		float t0, dx, dy, ox0, oy0;
		dx = x0 - dcam.P.x;
		dy = y0 - dcam.P.y;
		t0 = atan2(dy,dx);
		
		ox0 = (peye->Ax * cos(t0) * cos(peye->Aa) - peye->Ay * sin(t0) * sin(peye->Aa));
		oy0 = (peye->Ax * cos(t0) * sin(peye->Aa) + peye->Ay * sin(t0) * cos(peye->Aa));
		
		{
			si idx = Eye_Points_FindRot (peye, i, M_PI);
			if (idx >= 0) {
				float x1 = peye->paPoint[idx].x, y1 = peye->paPoint[idx].y;
				float t1 = atan2(y1 - dcam.P.y, x1 - dcam.P.x);
				
				dnpix(x0,y0)->Y = 0xFF;
				dnpix(x0,y0)->U = 0x0;
				dnpix(x0,y0)->V = 0x0;
			//	dnpix(x1,y1)->Y = 0xFF;
			//	dnpix(x1,y1)->U = 0x0;
			//	dnpix(x1,y1)->V = 0x0;
				
				
			//	float ox0 = (peye->Ax * cos(t0) * cos(peye->Aa) - peye->Ay * sin(t0) * sin(peye->Aa));
			//	float oy0 = (peye->Ax * cos(t0) * sin(peye->Aa) + peye->Ay * sin(t0) * cos(peye->Aa));
				
				float oy1 = (peye->Ax * cos(t1) * sin(peye->Aa) + peye->Ay * sin(t1) * cos(peye->Aa));
				float ox1 = (peye->Ax * cos(t1) * cos(peye->Aa) - peye->Ay * sin(t1) * sin(peye->Aa));
				
			//	printf ("Heee? pic %f %f   %f %f\n", x0 - dcam.P.x, y0 - dcam.P.y, x1 - dcam.P.x, y1 - dcam.P.y);
			//	printf ("Heee? eli %f %f   %f %f\n", ox0, oy0, ox1, oy1);
				
				float sdx = fabsf(x1 - x0) - fabsf(ox1 - ox0);
				float sdy = fabsf(y1 - y0) - fabsf(oy1 - oy0);
				
			//	ax += fabsf((x1 - x0)*cos(t0));
			//	ay += fabsf((y1 - y0)*sin(t0));
			//	++n;
				
				float pdx = (x0+x1)*0.5f - (2.0f*dcam.P.x + ox0+ox1)*0.5f;
				float pdy = (y0+y1)*0.5f - (2.0f*dcam.P.y + oy0+oy1)*0.5f;
				
				Eye_Xset (peye, dcam.P.x + 0.10f * pdx );// *cos(t));
				Eye_Yset (peye, dcam.P.y + 0.10f * pdy );// *sin(t));
				
				peye->Ax += 0.10f*sdx*fabs(cos(t0));
				peye->Ay += 0.10f*sdy*fabs(sin(t0));
			}else {
				dnpix(x0,y0)->Y = 0xFF;
				dnpix(x0,y0)->U = 0xF;
				dnpix(x0,y0)->V = 0xF;
				float diff = sqrt(dx*dx+dy*dy) - sqrt(ox0*ox0+oy0*oy0);
				if (fabs(diff) >= 0.01f) {
				//	printf ("f1\n");
					peye->Ax += peye->Fit_Scale*diff*fabs(cos(t0));
					peye->Ay += peye->Fit_Scale*diff*fabs(sin(t0));
					
					Eye_Xset (peye, dcam.P.x + peye->Fit_Trans*(diff )*cos(t0));
					Eye_Yset (peye, dcam.P.y + peye->Fit_Trans*(diff )*sin(t0));
				}
			}
		}
	}/**/
	float avgerr = 0;
	for (i = 0; i < peye->Point_N; ++i) {
		float x = peye->paPoint[i].x, y = peye->paPoint[i].y;
		float t, dx, dy, xx, yy;
		dx = x - dcam.P.x;
		dy = y - dcam.P.y;
		t = atan2(dy,dx);
		
		xx = (peye->Ax * cos(t) * cos(peye->Aa) - peye->Ay * sin(t) * sin(peye->Aa));
		yy = (peye->Ax * cos(t) * sin(peye->Aa) + peye->Ay * sin(t) * cos(peye->Aa));
		
		float diff = sqrt(dx*dx+dy*dy) - sqrt(xx*xx+yy*yy);
	//	avgerr += fabsf(diff);
		avgerr += dpow2(diff);
	}/**/
	return avgerr / peye->Point_N;
}

float	Eye_Points_Fit_Tri	(tEye* peye)
{
	si i;
	for (i = 0; i < peye->Point_N; ++i) {
		float x0 = peye->paPoint[i].x, y0 = peye->paPoint[i].y;
		float t0, dx, dy, ox0, oy0;
		dx = x0 - dcam.P.x;
		dy = y0 - dcam.P.y;
		t0 = atan2(dy,dx);
		
		ox0 = (peye->Ax * cos(t0) * cos(peye->Aa) - peye->Ay * sin(t0) * sin(peye->Aa));
		oy0 = (peye->Ax * cos(t0) * sin(peye->Aa) + peye->Ay * sin(t0) * cos(peye->Aa));
		
		tV2f p0 = {x0, y0}, op0 = {ox0, oy0};
		{
			si i1 = Eye_Points_FindRot (peye, i, M_PI_4/2);
			si i2 = Eye_Points_FindRot (peye, i, -M_PI_4/2);
			
			if (i1 >= 0 && i2 >= 0) {
				float x1 = peye->paPoint[i1].x, y1 = peye->paPoint[i1].y;
				float t1 = atan2(y1 - dcam.P.y, x1 - dcam.P.x);
				
				float x2 = peye->paPoint[i2].x, y2 = peye->paPoint[i2].y;
				float t2 = atan2(y2 - dcam.P.y, x2 - dcam.P.x);
				
				
				float oy1 = (peye->Ax * cos(t1) * sin(peye->Aa) + peye->Ay * sin(t1) * cos(peye->Aa));
				float ox1 = (peye->Ax * cos(t1) * cos(peye->Aa) - peye->Ay * sin(t1) * sin(peye->Aa));
				
				float oy2 = (peye->Ax * cos(t2) * sin(peye->Aa) + peye->Ay * sin(t2) * cos(peye->Aa));
				float ox2 = (peye->Ax * cos(t2) * cos(peye->Aa) - peye->Ay * sin(t2) * sin(peye->Aa));
				
				tV2f p1 = {x1, y1}, op1 = {ox1, oy1};
				tV2f p2 = {x2, y2}, op2 = {ox2, oy2};
				
				tV2f v1 = p1;	V2f_sub_V2f (&v1, &p0);
				tV2f v2 = p2;	V2f_sub_V2f (&v2, &p0);
				
				tV2f ov1 = op1;	V2f_sub_V2f (&ov1, &op0);
				tV2f ov2 = op2;	V2f_sub_V2f (&ov2, &op0);
				
				float dot = V2f_dot_V2f (&v1, &v2);
				float odot = V2f_dot_V2f (&ov1, &ov2);
				
				float diff = odot-dot;
				peye->Ax += +0.01f*diff*fabsf(sin(t0))	-0.01f*diff*fabsf(cos(t0));
				peye->Ay += -0.01f*diff*fabsf(sin(t0))	+0.01f*diff*fabsf(cos(t0));
				
				
			/*	float sdx = fabsf(x1 - x0) - fabsf(ox1 - ox0);
				float sdy = fabsf(y1 - y0) - fabsf(oy1 - oy0);
				
				float pdx = (x0+x1)*0.5f - (2.0f*dcam.P.x + ox0+ox1)*0.5f;
				float pdy = (y0+y1)*0.5f - (2.0f*dcam.P.y + oy0+oy1)*0.5f;
				
				Eye_Xset (peye, dcam.P.x + 0.10f * pdx );// *cos(t));
				Eye_Yset (peye, dcam.P.y + 0.10f * pdy );// *sin(t));
				
				peye->Ax += 0.10f*sdx*fabs(cos(t0));
				peye->Ay += 0.10f*sdy*fabs(sin(t0));/**/
				
				dnpix(x0,y0)->Y = 0xFF;
				dnpix(x0,y0)->U = 0x0;
				dnpix(x0,y0)->V = 0x0;
			}else {
				dnpix(x0,y0)->Y = 0xFF;
				dnpix(x0,y0)->U = 0xF;
				dnpix(x0,y0)->V = 0xF;
			}
			if (1) {
				si idx = Eye_Points_FindRot (peye, i, M_PI);
				if (idx >= 0) {
					float x1 = peye->paPoint[idx].x, y1 = peye->paPoint[idx].y;
					float t1 = atan2(y1 - dcam.P.y, x1 - dcam.P.x);
					
				//	dnpix(x0,y0)->Y = 0xFF;
				//	dnpix(x0,y0)->U = 0x0;
				//	dnpix(x0,y0)->V = 0x0;
				//	dnpix(x1,y1)->Y = 0xFF;
				//	dnpix(x1,y1)->U = 0x0;
				//	dnpix(x1,y1)->V = 0x0;
					
					
				//	float ox0 = (peye->Ax * cos(t0) * cos(peye->Aa) - peye->Ay * sin(t0) * sin(peye->Aa));
				//	float oy0 = (peye->Ax * cos(t0) * sin(peye->Aa) + peye->Ay * sin(t0) * cos(peye->Aa));
					
					float oy1 = (peye->Ax * cos(t1) * sin(peye->Aa) + peye->Ay * sin(t1) * cos(peye->Aa));
					float ox1 = (peye->Ax * cos(t1) * cos(peye->Aa) - peye->Ay * sin(t1) * sin(peye->Aa));
					
				//	printf ("Heee? pic %f %f   %f %f\n", x0 - dcam.P.x, y0 - dcam.P.y, x1 - dcam.P.x, y1 - dcam.P.y);
				//	printf ("Heee? eli %f %f   %f %f\n", ox0, oy0, ox1, oy1);
					
					float sdx = fabsf(x1 - x0) - fabsf(ox1 - ox0);
					float sdy = fabsf(y1 - y0) - fabsf(oy1 - oy0);
					
				//	ax += fabsf((x1 - x0)*cos(t0));
				//	ay += fabsf((y1 - y0)*sin(t0));
				//	++n;
					
					float pdx = (x0+x1)*0.5f - (2.0f*dcam.P.x + ox0+ox1)*0.5f;
					float pdy = (y0+y1)*0.5f - (2.0f*dcam.P.y + oy0+oy1)*0.5f;
					
					Eye_Xset (peye, dcam.P.x + 0.10f * pdx );// *cos(t));
					Eye_Yset (peye, dcam.P.y + 0.10f * pdy );// *sin(t));
					
					peye->Ax += 0.10f*sdx*fabs(cos(t0));
					peye->Ay += 0.10f*sdy*fabs(sin(t0));
					continue;
				}
			}
			if (1) {
				float diff = sqrt(dx*dx+dy*dy) - sqrt(ox0*ox0+oy0*oy0);
				if (fabs(diff) >= 0.01f) {
				//	printf ("f1\n");
					peye->Ax += 0.1f*diff*fabs(cos(t0));
					peye->Ay += 0.1f*diff*fabs(sin(t0));
					
					Eye_Xset (peye, dcam.P.x + 0.1f*(diff)*cos(t0));
					Eye_Yset (peye, dcam.P.y + 0.1f*(diff)*sin(t0));
				}
			}
		}
	}/**/
	return 0;
}


float	Eye_Points_Fit_Const	(tEye* peye)
{
	peye->Ax = peye->PFit_R/2;
	peye->Ay = peye->PFit_R/2;
	si i;
	for (i = 0; i < peye->Point_N; ++i) {
		float x0 = peye->paPoint[i].x, y0 = peye->paPoint[i].y;
		float t0, dx, dy, ox0, oy0;
		dx = x0 - dcam.P.x;
		dy = y0 - dcam.P.y;
		t0 = atan2(dy,dx);
		
		ox0 = (peye->Ax * cos(t0) * cos(peye->Aa) - peye->Ay * sin(t0) * sin(peye->Aa));
		oy0 = (peye->Ax * cos(t0) * sin(peye->Aa) + peye->Ay * sin(t0) * cos(peye->Aa));
		
		tV2f p0 = {x0, y0}, op0 = {ox0, oy0};
		{
			if (1) {
				si idx = Eye_Points_FindRot (peye, i, M_PI);
				if (idx >= 0) {
					float x1 = peye->paPoint[idx].x, y1 = peye->paPoint[idx].y;
					float t1 = atan2(y1 - dcam.P.y, x1 - dcam.P.x);
					
					dnpix(x0,y0)->Y = 0xFF;
					dnpix(x0,y0)->U = 0x0;
					dnpix(x0,y0)->V = 0x0;
				//	dnpix(x1,y1)->Y = 0xFF;
				//	dnpix(x1,y1)->U = 0x0;
				//	dnpix(x1,y1)->V = 0x0;
					
					
				//	float ox0 = (peye->Ax * cos(t0) * cos(peye->Aa) - peye->Ay * sin(t0) * sin(peye->Aa));
				//	float oy0 = (peye->Ax * cos(t0) * sin(peye->Aa) + peye->Ay * sin(t0) * cos(peye->Aa));
					
					float oy1 = (peye->Ax * cos(t1) * sin(peye->Aa) + peye->Ay * sin(t1) * cos(peye->Aa));
					float ox1 = (peye->Ax * cos(t1) * cos(peye->Aa) - peye->Ay * sin(t1) * sin(peye->Aa));
					
				//	printf ("Heee? pic %f %f   %f %f\n", x0 - dcam.P.x, y0 - dcam.P.y, x1 - dcam.P.x, y1 - dcam.P.y);
				//	printf ("Heee? eli %f %f   %f %f\n", ox0, oy0, ox1, oy1);
					
					float sdx = fabsf(x1 - x0) - fabsf(ox1 - ox0);
					float sdy = fabsf(y1 - y0) - fabsf(oy1 - oy0);
					
				//	ax += fabsf((x1 - x0)*cos(t0));
				//	ay += fabsf((y1 - y0)*sin(t0));
				//	++n;
					
					float pdx = (x0+x1)*0.5f - (2.0f*dcam.P.x + ox0+ox1)*0.5f;
					float pdy = (y0+y1)*0.5f - (2.0f*dcam.P.y + oy0+oy1)*0.5f;
					
					Eye_Xset (peye, dcam.P.x + 0.10f * pdx );// *cos(t));
					Eye_Yset (peye, dcam.P.y + 0.10f * pdy );// *sin(t));
					
				//	peye->Ax += 0.10f*sdx*fabs(cos(t0));
				//	peye->Ay += 0.10f*sdy*fabs(sin(t0));
					continue;
				}
			}
			if (1) {
				float diff = sqrt(dx*dx+dy*dy) - sqrt(ox0*ox0+oy0*oy0);
				if (fabs(diff) >= 0.01f) {
				//	printf ("f1\n");
				//	peye->Ax += 0.1f*diff*fabs(cos(t0));
				//	peye->Ay += 0.1f*diff*fabs(sin(t0));
					
					Eye_Xset (peye, dcam.P.x + 0.1f*(diff)*cos(t0));
					Eye_Yset (peye, dcam.P.y + 0.1f*(diff)*sin(t0));
				}
				dnpix(x0,y0)->Y = 0xFF;
				dnpix(x0,y0)->U = 0xF;
				dnpix(x0,y0)->V = 0xF;
			}
		}
	}/**/
	return 0;
}

float	Eye_Points_Fit		(tEye* peye)
{
	switch (peye->PFit) {
	case eEye_PFit_S1:
	case eEye_PFit_S2:
		Eye_Points_Fit2 (peye);
		Eye_Points_Fit2 (peye);
		break;
//	case eEye_PFit_RANSAC:
//		Eye_S4_Fit_Ransac (peye);
//		break;
	case eEye_PFit_S3:
		Eye_Points_Fit_Tri (peye);
		Eye_Points_Fit_Tri (peye);
		break;
	case eEye_PFit_S3const:
		Eye_Points_Fit_Const (peye);
		Eye_Points_Fit_Const (peye);
		break;
	}
}


si ay = 0, au = 0, av = 0, n = 0;

void	Eye_crap_ff		(tEye* peye, si x, si y)
{
	if (x < 0 || y < 0 || x >= pcam->Image_W || y >= pcam->Image_H)
		return;
	
	if (	1//dnpix(x,y)->Y == 0xFF
		&& dnpix(x,y)->U == 0
		&& dnpix(x,y)->V == 0
	)
		return;
	
	if (ddist2(x,y,dcam.P.x,dcam.P.y) >= peye->Exp_R*peye->Exp_R*2)
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

void	Eye_OldFF		(tEye* peye)
{
	si x, y;
//	si ay = 0, au = 0, av = 0, n = 0;
	ay = 0;
	au = 0;
	av = 0;
	n = 0;
	#define dd 4
	for (y = dcam.P.y-dd; y <= dcam.P.y+dd; ++y) {
		for (x = dcam.P.x-dd; x <= dcam.P.x+dd; ++x) {
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
	
//	ay = dnpix(dcam.P.x,dcam.P.y)->Y;
//	au = dnpix(dcam.P.x,dcam.P.y)->U;
//	av = dnpix(dcam.P.x,dcam.P.y)->V;
	
//	printf ("n %d\n", n);
	#define dd 2
	for (y = dcam.P.y-dd; y <= dcam.P.y+dd; ++y) {
		for (x = dcam.P.x-dd; x <= dcam.P.x+dd; ++x) {
			Eye_crap_ff (peye, x, y);
		}
	}
	#undef dd
	
}

void	Eye_CalcAYUV	(tEye* peye, tCam* pcam, si dd)
{
	si x, y;
//	si ay = 0, au = 0, av = 0, n = 0;
	ay = 0;
	au = 0;
	av = 0;
	n = 0;
	
	for (y = dcam.P.y-dd; y <= dcam.P.y+dd; ++y) {
		for (x = dcam.P.x-dd; x <= dcam.P.x+dd; ++x) {
			ay += dopix(x,y)->Y;
			au += dopix(x,y)->U;
			av += dopix(x,y)->V;
			++n;
		}
	}
	
	ay /= n;
	au /= n;
	av /= n;
}

void	Eye_V_Pre		(tEye* peye, tCam* pcam)	//motion compensation
{
	return;
	dcam.OP = dcam.P;
	dcam.P.x += dcam.V.x;
	dcam.P.y += dcam.V.y;
}
void	Eye_V_Post		(tEye* peye, tCam* pcam)
{
	return;
//	dcam.V.x = dcam.P.x - oldpos.x;
//	dcam.V.y = dcam.P.y - oldpos.y;
	dcam.V.x = (dcam.V.x + dcam.P.x - dcam.OP.x) * 0.5f;
	dcam.V.y = (dcam.V.y + dcam.P.y - dcam.OP.y) * 0.5f;
	
	dcam.V.x *= 0.9f;
	dcam.V.y *= 0.9f;
}


void	Eye_GV_Calc		(tEye* peye)
{
	SDL_mutexP (peye->GV_mutex);
	if (gM.Eye_GlintMode == 2) {
		tV2f v0 = dcam.P;	V2f_sub_V2f (&v0, &peye->G0);
		tV2f v1 = dcam.P;	V2f_sub_V2f (&v1, &peye->G1);
		tV2f v01 = peye->G1;	V2f_sub_V2f (&v01, &peye->G0);
		
		V2f_add_V2f (&v0, &v1);
		V2f_mul_S (&v0, 1/V2f_dist (&v01));
		peye->GV = v0;
	}else if (gM.Eye_GlintMode == 1) {
		tV2f v0 = dcam.P;	V2f_sub_V2f (&v0, &peye->G0);
		
		tV2f gtf = peye->G0;
		gtf.x -= gM.aCam[0].Image_W/2;
		gtf.y -= gM.aCam[0].Image_H/2;
		
		gtf.x *= peye->GTF.x;
		gtf.y *= peye->GTF.y;
		
		V2f_sub_V2f (&v0, &gtf);
		peye->GV = v0;
	}
	SDL_mutexV (peye->GV_mutex);
}

void	Eye_GV_Get	(tEye* peye, tV2f* pret)
{
	SDL_mutexP (peye->GV_mutex);
	*pret = peye->GV;
	SDL_mutexV (peye->GV_mutex);
}


void	Eye_GP_Calc		(tEye* peye)
{
	if (gM.Eye_GlintMode != 2)
		return;
	
	SDL_mutexP (peye->GV_mutex);
	
	tV4f g0_v, g1_v;
	Cam_Pos2Ray (peye->G0, NULL, NULL, &g0_v);
	Cam_Pos2Ray (peye->G1, NULL, NULL, &g1_v);
	V4f_norm (&g0_v);
	V4f_norm (&g1_v);
	float gd = V4f_dot_V4f (&g0_v, &g1_v);
	
	
	tV2f ag = peye->G0;	V2f_add_V2f (&ag, &peye->G1);
	ag.x /= 2;	ag.y /= 2;
	tV4f sr_p0, sr_norm, sr_p1;
	Cam_Pos2Ray (ag, &sr_p0, &sr_p1, &sr_norm);
	V4f_norm (&sr_norm);
	
	float z = 16;
	si i;
	for (i = 0; i < 100; ++i) {
		peye->GP = sr_norm;
		V4f_mul_S (&peye->GP, z);
		
		tV4f l0_v, l1_v;
		Sphere_Reflect (&l0_v, &peye->GP, peye->LR, &gM.aLight[0].P);
		Sphere_Reflect (&l1_v, &peye->GP, peye->LR, &gM.aLight[1].P);
		V4f_norm (&l0_v);
		V4f_norm (&l1_v);
		
		float ld = V4f_dot_V4f (&l0_v, &l1_v);
		
		z += 1000*(gd - ld);
		
	//	printf ("z %f	gd ld	%f %f  diff %e\n", z, gd, ld, gd-ld);
	//	printf ("g0_v "); V4f_Print (&g0_v);
	//	printf ("l0_v "); V4f_Print (&l0_v);
	//	printf ("g1_v "); V4f_Print (&g1_v);
	//	printf ("l1_v "); V4f_Print (&l1_v);
	}
	{
	//	printf ("start\n");
		float g0_ry = V4f_ry (&g0_v);
		float g1_ry = V4f_ry (&g1_v);
		float rot = 0.0f;
		for (i = 0; i < 10; ++i) {
			peye->GP = sr_norm;
			V4f_mul_S (&peye->GP, z);
			V4f_roty (&peye->GP, rot);
			
			tV4f l0_v, l1_v;
			Sphere_Reflect (&l0_v, &peye->GP, peye->LR, &gM.aLight[0].P);
			Sphere_Reflect (&l1_v, &peye->GP, peye->LR, &gM.aLight[1].P);
			V4f_norm (&l0_v);
			V4f_norm (&l1_v);
			
			float avg = V4f_ry(&l0_v)-g0_ry + V4f_ry(&l1_v)-g1_ry;
			
			rot += avg/2;
			
		//	printf ("rot %f	g0 l0	%f %f  diff %e\n", rot, g0_ry, V4f_ry(&l0_v), g0_ry-V4f_ry(&l0_v));
		//	printf ("rot %f	g1 l1	%f %f  diff %e\n", z, g1_ry, V4f_ry(&l1_v), g1_ry-V4f_ry(&l1_v));
		}
	}/**/
	
	SDL_mutexV (peye->GV_mutex);
}

void	Eye_GP_GazeVector	(tEye* peye, tV4f* pr0, tV4f* pr1, tV4f* prv)
{
	tV4f p0, p1, v, l0, l1;
	p0 = peye->GP;
	
	Cam_Retina_Ray (peye, &l0, &l1, 0);
	
	V4f_Intersect_Sphere0R_Line01	(
		&p1,
		&p0, peye->LR,
		&l0, &l1
	);
	if (pr0)
		*pr0 = p0;
	if (pr1)
		*pr1 = p1;
	if (prv) {
		*prv = p1;	V4f_sub_V4f (prv, &p0);
	}
	return;
}



void	Eye_EdgeMark	(tEye* peye)
{
	si x, y;
	for (y = 0; y < pcam->Image_H; ++y) {
		for (x = 0; x < pcam->Image_W-1; ++x) {
			if (	1//dnpix(x,y)->Y == 0xFF
				&& dnpix(x,y)->U == 0
				&& dnpix(x,y)->V == 0
			) {
				float t, dx, dy, xx, yy;
				dx = dcam.P.x - x;
				dy = dcam.P.y - y;
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
				dx = dcam.P.x - x;
				dy = dcam.P.y - y;
				t = atan2(dy,dx) - peye->Aa;
				
				xx = (peye->Ax * cos(t) * cos(peye->Aa) - peye->Ay * sin(t) * sin(peye->Aa));
				yy = (peye->Ax * cos(t) * sin(peye->Aa) + peye->Ay * sin(t) * cos(peye->Aa));
				
				float diff = sqrt(dx*dx+dy*dy) - sqrt(xx*xx+yy*yy);
				if (fabs(diff) >= 0.01f) {
				//	printf ("f1\n");
					peye->Ax += peye->Fit_Scale*diff*fabs(cos(t));
					peye->Ay += peye->Fit_Scale*diff*fabs(sin(t));
					
					Eye_Xset (peye, dcam.P.x - peye->Fit_Trans*diff*cos(t));
					Eye_Yset (peye, dcam.P.y - peye->Fit_Trans*diff*sin(t));
				}
			}
		}
	}/**/
	for (y = 0; y < pcam->Image_H; ++y) {
		for (x = 0; x < pcam->Image_W; ++x) {
			if (	1//dnpix(x,y)->Y == 0xFF
				&& dnpix(x,y)->U == 0xF
				&& dnpix(x,y)->V == 0xF
			) {
				float t, dx, dy, xx, yy;
				dx = x - dcam.P.x;
				dy = y - dcam.P.y;
				t = atan2(dy,dx)/* - peye->Aa*/;
				
				xx = (peye->Ax * cos(t) * cos(peye->Aa) - peye->Ay * sin(t) * sin(peye->Aa));
				yy = (peye->Ax * cos(t) * sin(peye->Aa) + peye->Ay * sin(t) * cos(peye->Aa));
				
				float diff = sqrt(dx*dx+dy*dy) - sqrt(xx*xx+yy*yy);
				if (fabs(diff) >= 0.1f) {
				//	printf ("f1\n");
					peye->Ax += peye->Fit_Scale*diff*fabs(cos(t));
					peye->Ay += peye->Fit_Scale*diff*fabs(sin(t));
					
					Eye_Xset (peye, dcam.P.x + peye->Fit_Trans*(diff )*cos(t));
					Eye_Yset (peye, dcam.P.y + peye->Fit_Trans*(diff )*sin(t));
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
	for (y = dcam.P.y - border; y < dcam.P.y + border; ++y) {
		for (x = dcam.P.x - border; x < dcam.P.x + border; ++x) {
		/*	if (	1//dnpix(x,y)->Y == 0xFF
				&& dnpix(x,y)->U == 0xF
				&& dnpix(x,y)->V == 0xF
			) {
				float t, dx, dy, xx, yy;
				dx = dcam.P.x - x;
				dy = dcam.P.y - y;
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
	float x = dcam.P.x, y = dcam.P.y;
	
	for (; ddist2(x,y,dcam.P.x,dcam.P.y) < peye->Exp_R*peye->Exp_R*3; ) {
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
		dx = dcam.P.x - x;
		dy = dcam.P.y - y;
		t = atan2(dy,dx) - peye->Aa;
		
		xx = (peye->Ax * cos(t) * cos(peye->Aa) - peye->Ay * sin(t) * sin(peye->Aa));
		yy = (peye->Ax * cos(t) * sin(peye->Aa) + peye->Ay * sin(t) * cos(peye->Aa));
		
		float diff = sqrt(dx*dx+dy*dy) - sqrt(xx*xx+yy*yy);
		if (fabs(diff) >= 0.01f) {
		//	printf ("f1\n");
			peye->Ax += peye->Fit_Scale*diff*fabs(cos(t));
			peye->Ay += peye->Fit_Scale*diff*fabs(sin(t));
			
			Eye_Xset (peye, dcam.P.x - peye->Fit_Trans*diff*cos(t));
			Eye_Yset (peye, dcam.P.y - peye->Fit_Trans*diff*sin(t));
		}*/
//	}
	float ax = 0, ay = 0;
	si n = 0;
	
	si x0, y0;
	for (y0 = dcam.P.y - border; y0 < dcam.P.y + border; ++y0) {		//will have to stop this crap soon
//	for (y0 = dcam.P.y + border; y0 > dcam.P.y - border; --y0) {
		for (x0 = dcam.P.x - border; x0 < dcam.P.x + border; ++x0) {
			if (	1//dnpix(x0,y0)->Y == 0xFF
				&& dnpix(x0,y0)->U == 0xF
				&& dnpix(x0,y0)->V == 0xF
			) {
				float t0, dx0, dy0;
				dx0 = x0 - dcam.P.x;
				dy0 = y0 - dcam.P.y;
				t0 = atan2(dy0,dx0)/* - peye->Aa*/;
				float t1 = t0 - M_PI;
				si x1, y1;
				
			//	printf ("Heee? eye %f %f %f %f\n", dcam.P.x, dcam.P.y, dx0, dy0);
				
				if (Eye_S2Fit_EdgeGet (peye, t1, &x1, &y1)) {
				//	printf ("Heee? eye %ld %ld   %ld %ld\n", x0 - (si)dcam.P.x, y0 - (si)dcam.P.y, x1 - (si)dcam.P.x, y1 - (si)dcam.P.y);
				/*	dnpix(x0,y0)->Y = 0x0;
					dnpix(x0,y0)->U = 0x0;
					dnpix(x0,y0)->V = 0x0;
					
					dnpix(x1,y1)->Y = 0xFF;
					dnpix(x1,y1)->U = 0xF;
					dnpix(x1,y1)->V = 0xF;/**/
					if (0) {
						float dx1, dy1;
						dx1 = x1 - dcam.P.x;
						dy1 = y1 - dcam.P.y;
						
						float ox0 = (peye->Ax * cos(t0) * cos(peye->Aa) - peye->Ay * sin(t0) * sin(peye->Aa));
						float ox1 = (peye->Ax * cos(t1) * cos(peye->Aa) - peye->Ay * sin(t1) * sin(peye->Aa));
						
						float oy0 = (peye->Ax * cos(t0) * sin(peye->Aa) + peye->Ay * sin(t0) * cos(peye->Aa));
						float oy1 = (peye->Ax * cos(t1) * sin(peye->Aa) + peye->Ay * sin(t1) * cos(peye->Aa));
						
					/*	float diff = sqrt(dx0*dx0+dy0*dy0) - sqrt(ox0*ox0+oy0*oy0);
						if (fabs(diff) >= 0.01f) {
						//	printf ("f1\n");
							peye->Ax += peye->Fit_Scale*diff*fabs(cos(t0));
							peye->Ay += peye->Fit_Scale*diff*fabs(sin(t0));
							
							Eye_Xset (peye, dcam.P.x + peye->Fit_Trans*(diff )*cos(t0));
							Eye_Yset (peye, dcam.P.y + peye->Fit_Trans*(diff )*sin(t0));
						}/**/
					/*	float diff = sqrt(dx1*dx1+dy1*dy1) - sqrt(ox1*ox1+oy1*oy1);
						if (fabs(diff) >= 0.01f) {
						//	printf ("f1\n");
							peye->Ax += peye->Fit_Scale*diff*fabs(cos(t1));
							peye->Ay += peye->Fit_Scale*diff*fabs(sin(t1));
							
							Eye_Xset (peye, dcam.P.x + peye->Fit_Trans*(diff )*cos(t1));
							Eye_Yset (peye, dcam.P.y + peye->Fit_Trans*(diff )*sin(t1));
						}/**/
					}else {
						float ox0 = (peye->Ax * cos(t0) * cos(peye->Aa) - peye->Ay * sin(t0) * sin(peye->Aa));
						float ox1 = (peye->Ax * cos(t1) * cos(peye->Aa) - peye->Ay * sin(t1) * sin(peye->Aa));
						
						float oy0 = (peye->Ax * cos(t0) * sin(peye->Aa) + peye->Ay * sin(t0) * cos(peye->Aa));
						float oy1 = (peye->Ax * cos(t1) * sin(peye->Aa) + peye->Ay * sin(t1) * cos(peye->Aa));
						
					//	printf ("Heee? pic %f %f   %f %f\n", x0 - dcam.P.x, y0 - dcam.P.y, x1 - dcam.P.x, y1 - dcam.P.y);
					//	printf ("Heee? eli %f %f   %f %f\n", ox0, oy0, ox1, oy1);
						
						float sdx = fabsf(x1 - x0) - fabsf(ox1 - ox0);
						float sdy = fabsf(y1 - y0) - fabsf(oy1 - oy0);
						
						ax += fabsf((x1 - x0)*cos(t0));
						ay += fabsf((y1 - y0)*sin(t0));
						++n;
						
						float pdx = 0.2f + (x0+x1)*0.5f - (2*dcam.P.x + ox0+ox1)*0.5f;
						float pdy = -0.1f + (y0+y1)*0.5f - (2*dcam.P.y + oy0+oy1)*0.5f;
						
						Eye_Xset (peye, dcam.P.x + peye->S2Fit_Trans * pdx );// *cos(t));
						Eye_Yset (peye, dcam.P.y + peye->S2Fit_Trans * pdy );// *sin(t));
						
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
						
						Eye_Xset (peye, dcam.P.x + peye->Fit_Trans * pdx );// *cos(t));
						Eye_Yset (peye, dcam.P.y + peye->Fit_Trans * pdy );// *sin(t));
						/**/
						
					//	Eye_Xset (peye, dcam.P.x + peye->Fit_Trans * ( (x0+x1)*0.5f - dcam.P.x-(ox0+ox1)*0.5f ) );// *cos(t));
					//	Eye_Yset (peye, dcam.P.y + peye->Fit_Trans * ( (y0+y1)*0.5f - dcam.P.y-(oy0+oy1)*0.5f ) );// *sin(t));
						
					/*	dnpix(ox0, oy0)->Y = 0x0;
						dnpix(ox0, oy0)->U = 0x0;
						dnpix(ox0, oy0)->V = 0x0;
						
						dnpix(ox1, oy1)->Y = 0xFF;
						dnpix(ox1, oy1)->U = 0xF;
						dnpix(ox1, oy1)->V = 0xF;/**/
					}
				}else {
				//	printf ("Heee? eye %f %f %ld %ld\n", dcam.P.x, dcam.P.y, x0, y0);
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



si	Eye_S0_EdgeMark	(tEye* peye, si i, float t, float* px, float* py)
{
//	si border = peye->Exp_R*peye->Exp_R*3;
	u08 wrote = 0;
	float x = dcam.P.x, y = dcam.P.y;
	si n = 0;
	
	for (; ddist2(x,y,dcam.P.x,dcam.P.y) < peye->Exp_R*peye->Exp_R; ) {
		if (peye->LinView.x != 0)
			dspix(peye->LinView.x+n, peye->LinView.y+i) = dnpix(x,y)->Y | dnpix(x,y)->Y<<8 | dnpix(x,y)->Y<<16;
		
		if (	dopix(x,y)->Y > ay + 50
		) {
		//	return 1;
			float x1 = x, y1 = y;
			for (; ddist2(x1,y1,dcam.P.x,dcam.P.y) < peye->Exp_R*peye->Exp_R*4; ) {
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
					dspix(peye->LinView.x+n, peye->LinView.y+i) = 0xFF<<8;
				}
				dspix(peye->LinView.x+n, peye->LinView.y+i) = 0xFF<<8;
			}
		}/**/
		
		x += cos(t);
		y += sin(t);
		++n;
	}/**/
	return wrote;
	return 0;
}

void	Eye_S0_Fit		(tEye* peye)
{
	SDL_LockSurface (gM.pScreen);
//	Eye_SFit (peye);	return;
	si border = 4*peye->Exp_R;
//	float ax = 0, ay = 0;
//	si n = 0;
	si i = 0;
	
	si	(*edge)	(tEye* peye, si i, float t, float* px, float* py);
	edge = Eye_S0_EdgeMark;
	
	float ao, a;
//	for (a = 0; a < M_PI; a += M_PI/100.0f) {
//	for (a = -M_PI_2; a < M_PI_2; a += M_PI/100.0f) {
//	for (a = -M_PI_4; a < M_PI_4; a += M_PI/100.0f) {
	for (ao = 0; ao < M_PI_2; ao += M_PI/100.0f) {
	for (a = ao; a <= ao+M_PI_2+M_PI_4; a += M_PI_2) {
//	{
	//	printf ("mino ");
		u08 skip;
		float x0, y0, x1, y1;
		float t0 = a, t1 = t0 - M_PI;
		
	//	skip = edge (peye, i, t0, &x0, &y0);
		skip = edge (peye, a/(M_PI/100.0f), t0, &x0, &y0);
		
		if (skip == 0) {
			t0 = t1;
		//	skip = edge (peye, 110 + i, t0, &x0, &y0);
			skip = edge (peye, 110 + a/(M_PI/100.0f), t0, &x0, &y0);
		}else {
		//	skip += edge (peye, 110 + i, t1, &x1, &y1);
			skip += edge (peye, 110 + a/(M_PI/100.0f), t1, &x1, &y1);
		}
		++i;
		
		if (skip == 0)
			continue;
		
		if (skip == 1) {
		//	printf ("fuuu\n");
			if (!dpixout(x0,y0)) {
				dnpix(x0,y0)->Y = 0xFF;
				dnpix(x0,y0)->U = 0xF;
				dnpix(x0,y0)->V = 0xF;
			}
			float t, dx, dy, xx, yy;
			dx = x0 - dcam.P.x;
			dy = y0 - dcam.P.y;
			t = atan2(dy,dx)/* - peye->Aa*/;
			
			xx = (peye->Ax * cos(t) * cos(peye->Aa) - peye->Ay * sin(t) * sin(peye->Aa));
			yy = (peye->Ax * cos(t) * sin(peye->Aa) + peye->Ay * sin(t) * cos(peye->Aa));
			
			float diff = sqrt(dx*dx+dy*dy) - sqrt(xx*xx+yy*yy);
			if (fabs(diff) >= 0.1f) {
			//	printf ("f1\n");
			//	peye->Ax += peye->Fit_Scale*diff*fabs(cos(t));
			//	peye->Ay += peye->Fit_Scale*diff*fabs(sin(t));
				
				Eye_Xset (peye, dcam.P.x + 1*peye->Fit_Trans*(diff )*cos(t));
				Eye_Yset (peye, dcam.P.y + 1*peye->Fit_Trans*(diff )*sin(t));
			}
			continue;
		}
		if (!dpixout(x0,y0)) {
			dnpix(x0,y0)->Y = 0xFF;
			dnpix(x0,y0)->U = 0x0;
			dnpix(x0,y0)->V = 0x0;
		}
		if (!dpixout(x1,y1)) {
			dnpix(x1,y1)->Y = 0xFF;
			dnpix(x1,y1)->U = 0x0;
			dnpix(x1,y1)->V = 0x0;
		}
	//	printf ("Heee? eye %ld %ld   %ld %ld\n", x0 - (si)dcam.P.x, y0 - (si)dcam.P.y, x1 - (si)dcam.P.x, y1 - (si)dcam.P.y);
	/*	dnpix(x0,y0)->Y = 0x0;
		dnpix(x0,y0)->U = 0x0;
		dnpix(x0,y0)->V = 0x0;
		
		dnpix(x1,y1)->Y = 0xFF;
		dnpix(x1,y1)->U = 0xF;
		dnpix(x1,y1)->V = 0xF;/**/
		if (0) {
			float dx1, dy1;
			dx1 = x1 - dcam.P.x;
			dy1 = y1 - dcam.P.y;
			
			float ox0 = (peye->Ax * cos(t0) * cos(peye->Aa) - peye->Ay * sin(t0) * sin(peye->Aa));
			float oy0 = (peye->Ax * cos(t0) * sin(peye->Aa) + peye->Ay * sin(t0) * cos(peye->Aa));
			
			float ox1 = (peye->Ax * cos(t1) * cos(peye->Aa) - peye->Ay * sin(t1) * sin(peye->Aa));
			float oy1 = (peye->Ax * cos(t1) * sin(peye->Aa) + peye->Ay * sin(t1) * cos(peye->Aa));
			
		/*	float diff = sqrt(dx0*dx0+dy0*dy0) - sqrt(ox0*ox0+oy0*oy0);
			if (fabs(diff) >= 0.01f) {
			//	printf ("f1\n");
				peye->Ax += peye->Fit_Scale*diff*fabs(cos(t0));
				peye->Ay += peye->Fit_Scale*diff*fabs(sin(t0));
				
				Eye_Xset (peye, dcam.P.x + peye->Fit_Trans*(diff )*cos(t0));
				Eye_Yset (peye, dcam.P.y + peye->Fit_Trans*(diff )*sin(t0));
			}/**/
		/*	float diff = sqrt(dx1*dx1+dy1*dy1) - sqrt(ox1*ox1+oy1*oy1);
			if (fabs(diff) >= 0.01f) {
			//	printf ("f1\n");
				peye->Ax += peye->Fit_Scale*diff*fabs(cos(t1));
				peye->Ay += peye->Fit_Scale*diff*fabs(sin(t1));
				
				Eye_Xset (peye, dcam.P.x + peye->Fit_Trans*(diff )*cos(t1));
				Eye_Yset (peye, dcam.P.y + peye->Fit_Trans*(diff )*sin(t1));
			}/**/
		}else {
			float ox0 = (peye->Ax * cos(t0) * cos(peye->Aa) - peye->Ay * sin(t0) * sin(peye->Aa));
			float ox1 = (peye->Ax * cos(t1) * cos(peye->Aa) - peye->Ay * sin(t1) * sin(peye->Aa));
			
			float oy0 = (peye->Ax * cos(t0) * sin(peye->Aa) + peye->Ay * sin(t0) * cos(peye->Aa));
			float oy1 = (peye->Ax * cos(t1) * sin(peye->Aa) + peye->Ay * sin(t1) * cos(peye->Aa));
			
		//	printf ("Heee? pic %f %f   %f %f\n", x0 - dcam.P.x, y0 - dcam.P.y, x1 - dcam.P.x, y1 - dcam.P.y);
		//	printf ("Heee? eli %f %f   %f %f\n", ox0, oy0, ox1, oy1);
			
			float sdx = fabsf(x1 - x0) - fabsf(ox1 - ox0);
			float sdy = fabsf(y1 - y0) - fabsf(oy1 - oy0);
			
		//	ax += fabsf((x1 - x0)*cos(t0));
		//	ay += fabsf((y1 - y0)*sin(t0));
		//	++n;
			
			float pdx = (x0+x1)*0.5f - (2.0f*dcam.P.x + ox0+ox1)*0.5f;
			float pdy = (y0+y1)*0.5f - (2.0f*dcam.P.y + oy0+oy1)*0.5f;
			
			Eye_Xset (peye, dcam.P.x + 0.10f * pdx );// *cos(t));
			Eye_Yset (peye, dcam.P.y + 0.10f * pdy );// *sin(t));
			
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
			
			Eye_Xset (peye, dcam.P.x + peye->Fit_Trans * pdx );// *cos(t));
			Eye_Yset (peye, dcam.P.y + peye->Fit_Trans * pdy );// *sin(t));
			/**/
			
		//	Eye_Xset (peye, dcam.P.x + peye->Fit_Trans * ( (x0+x1)*0.5f - dcam.P.x-(ox0+ox1)*0.5f ) );// *cos(t));
		//	Eye_Yset (peye, dcam.P.y + peye->Fit_Trans * ( (y0+y1)*0.5f - dcam.P.y-(oy0+oy1)*0.5f ) );// *sin(t));
			
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

void	Eye_S0		(tEye* peye)
{
//	Eye_V_Pre (peye);
	
	Eye_S0_Fit (peye);
	
//	Eye_V_Post (peye);
}


#define dmap_pix_ad(a,d)	dspix(peye->LinView.x+ 2*(d), peye->LinView.y+ (si)(angle_norm_0_2pi(a)/(peye->AngRes*deg2rad)))

void	Eye_Ellipse2LinDraw	(tEye* peye)
{
//	return;
	if (peye->LinView.x == 0)
		return;
//	si border = peye->Exp_R*peye->Exp_R*3;
	si i = 0;
	
	float a;
	for (a = 0; a < 2*M_PI; a += peye->AngRes*deg2rad) {
		si n = 0;
		float x = dcam.P.x, y = dcam.P.y;
		tPix prevpix = *dopix(x,y);
		for (; ddist2(x,y,dcam.P.x,dcam.P.y) < dpow2(peye->Exp_R*2); ) {
			dspix(peye->LinView.x+n, peye->LinView.y+i) = dmono2rgb(dopix(x,y)->Y);
			dspix(peye->LinView.x+peye->Exp_R*2*2+n, peye->LinView.y+i) = dmono2rgb(2*abs(dopix(x,y)->Y-prevpix.Y));
			++n;
			dspix(peye->LinView.x+n, peye->LinView.y+i) = dmono2rgb(dopix(x,y)->Y);
			dspix(peye->LinView.x+peye->Exp_R*2*2+n, peye->LinView.y+i) = dmono2rgb(2*abs(dopix(x,y)->Y-prevpix.Y));
			++n;
			prevpix = *dopix(x,y);
			
			x += cos(a);
			y += sin(a);
		}
	//	dspix(peye->LinView.x+(2*peye->Exp_R), peye->LinView.y+i) = 0xFF<<16;
	//	dspix(peye->LinView.x+(2*peye->Exp_R+1), peye->LinView.y+i) = 0xFF<<16;
		++i;
	}/**/
}
void	Eye_Ellipse2LinDraw_Pix_ad	(tEye* peye, float a, float d, u32 col)
{
//	return;
	dmap_pix_ad(a, d) = col;
//	dmap_pix_ad(a, peye->Exp_R*2.0f + d) = col;
}

void	Eye_CirView_Point_xy		(tEye* peye, float x, float y, u32 col)
{
	float sx = 2.0f, sy = sx;
	if (peye->CirView.x == 0)
		return;
//	printf ("CirView out	%f	%f\n", x*sx, y*sy);
	si xx = peye->CirView.x + x*sx, yy = peye->CirView.y + y*sy;
	if (!dsout(xx,yy))
		dspix(xx,yy) = col;
}


si	Eye_S3Fit_EdgeMark	(tEye* peye, si i, float t, float* px, float* py)
{
//	si border = peye->Exp_R*peye->Exp_R*3;
	u08 wrote = 0;
	float x = dcam.P.x, y = dcam.P.y;
	si n = 0;
	
	for (; ddist2(x,y,dcam.P.x,dcam.P.y) < peye->Exp_R*peye->Exp_R; ) {
		if (peye->LinView.x != 0)
			dspix(peye->LinView.x+n, peye->LinView.y+i) = dnpix(x,y)->Y | dnpix(x,y)->Y<<8 | dnpix(x,y)->Y<<16;
		
		if (	dopix(x,y)->Y > peye->Pix_Bright
		) {
		//	return 1;
			float x1 = x, y1 = y;
			for (; ddist2(x1,y1,dcam.P.x,dcam.P.y) < peye->Exp_R*peye->Exp_R*4; ) {
				if (	dopix(x1,y1)->Y <= peye->Pix_Bright + 8) {
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
					dspix(peye->LinView.x+n, peye->LinView.y+i) = 0xFF<<8;
				}
				dspix(peye->LinView.x+n, peye->LinView.y+i) = 0xFF<<8;
			}
		}/**/
		
		x += cos(t);
		y += sin(t);
		++n;
	}/**/
	return wrote;
	return 0;
}
si	Eye_S3Fit_EdgeMark2	(tEye* peye, si i, float t, float* px, float* py)
{
	float border_s = peye->Min_R;
	float border_e = peye->Max_R;
	
	u08 wrote = 0;
	float x = dcam.P.x, y = dcam.P.y;
	si n = 0;
	
	for (; ddist2(x,y,dcam.P.x,dcam.P.y) < dpow2(border_e); ) {
		if (peye->LinView.x != 0)
			dspix(peye->LinView.x+n, peye->LinView.y+i) = dnpix(x,y)->Y | dnpix(x,y)->Y<<8 | dnpix(x,y)->Y<<16;
		
		if (
		//	dopix(x,y)->Y > ay + 30
			dopix(x,y)->Y >= peye->Pix_Bright
		//	&& ddist2(x,y,dcam.P.x,dcam.P.y) >= 14*14
			&& ddist2(x,y,dcam.P.x,dcam.P.y) >= dpow2(border_s)
		) {
		//	return 1;
			float x1 = x, y1 = y;
			for (; ddist2(x1,y1,dcam.P.x,dcam.P.y) < dpow2(peye->Exp_R*2); ) {
				if (	dopix(x1,y1)->Y <= peye->Pix_Bright + 8) {
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
					dspix(peye->LinView.x+n, peye->LinView.y+i) = 0xFF<<8;
				}
				dspix(peye->LinView.x+n, peye->LinView.y+i) = 0xFF<<8;
			}
		}/**/
		
		x += cos(t);
		y += sin(t);
		++n;
	}/**/
	return wrote;
	return 0;
}
si	Eye_S3Fit_EdgeMark3	(tEye* peye, si i, float t, float* px, float* py)
{
//	si border = peye->Exp_R*peye->Exp_R*3;
	u08 wrote = 0;
	float x = dcam.P.x, y = dcam.P.y, mindiff = peye->Exp_R*peye->Exp_R;
	si n = 0;
	
	for (; ddist2(x,y,dcam.P.x,dcam.P.y) < peye->Exp_R*peye->Exp_R; ) {
		if (peye->LinView.x != 0)
			dspix(peye->LinView.x+n, peye->LinView.y+i) = dnpix(x,y)->Y | dnpix(x,y)->Y<<8 | dnpix(x,y)->Y<<16;
		
		if (	dopix(x,y)->Y > ay + 18
		//	&& ddist2(x,y,dcam.P.x,dcam.P.y) >= 14*14
		) {
		//	return 1;
			float x1 = x, y1 = y;
			for (; ddist2(x1,y1,dcam.P.x,dcam.P.y) < peye->Exp_R*peye->Exp_R*4; ) {
				if (	dopix(x1,y1)->Y <= ay + 8) {
					break;
				}
				x1 += cos(t);
				y1 += sin(t);
			}
			if (ddist2(x1,y1,x,y) >= 2*2) {
			//	if (!wrote) {
					float diff = ddist2(x,y,dcam.P.x,dcam.P.y) - peye->Exp_R*peye->Exp_R;
					if (fabsf(diff) < mindiff) {
						wrote = 1;
						mindiff = fabsf(diff);
						*px = x+0.000000000001f;
						*py = y+0.000000000001f;
					}
					dspix(peye->LinView.x+n, peye->LinView.y+i) = 0xFF<<8;
			//	}
				dspix(peye->LinView.x+n, peye->LinView.y+i) = 0xFF<<8;
				x = x1;
				y = y1;
				continue;
			}
		}/**/
		
		x += cos(t);
		y += sin(t);
		++n;
	}/**/
	return wrote;
	return 0;
}

si	Eye_S3Fit_EdgeMark4	(tEye* peye, si i, float t, float* px, float* py)
{
	u08 wrote = 0;
	float x = dcam.P.x, y = dcam.P.y;
	si n = 0;
	
	for (; ddist2(x,y,dcam.P.x,dcam.P.y) < dpow2(peye->Max_R); ) {
		if (peye->LinView.x != 0)
			dspix(peye->LinView.x+n, peye->LinView.y+i) = dnpix(x,y)->Y | dnpix(x,y)->Y<<8 | dnpix(x,y)->Y<<16;
		
		if (	dopix(x,y)->Y > peye->Pix_Bright
		) {
			if (!wrote) {
				wrote = 1;
				float x0 = x, y0 = y, x1 = x, y1 = y;
				x0 -= cos(t);	y0 -= sin(t);
				x1 += cos(t);	y1 += sin(t);
				/*
				float a = (float)(dopix(x1,y1)->Y - dopix(x,y)->Y) / 180;
				
				printf ("a %f\n", a);
				if (a < 0)		a = 0;
				else if (a > 1)	a = 1;
				a = 1.0f-a;
				*px = x*a + x1*(1.0f-a) + 0.000000000001f;
				*py = y*a + y1*(1.0f-a) + 0.000000000001f;/**/
				
				float a = (float)(dopix(x,y)->Y - dopix(x0,y0)->Y) / 180;
			//	printf ("a %f\n", a);
				
				if (a < 0)		a = 0;
				else if (a > 1)	a = 1;
				a = 1.0f-a;
				*px = x*a + x0*(1.0f-a) + 0.000000000001f;
				*py = y*a + y0*(1.0f-a) + 0.000000000001f;/**/
				
			//	*px = x+0.000000000001f;	*py = y+0.000000000001f;
				dspix(peye->LinView.x+n, peye->LinView.y+i) = 0xFF<<8;
			//	return 1;
			}
			dspix(peye->LinView.x+n, peye->LinView.y+i) = 0xFF<<8;
			
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
//	si border = 4*peye->Exp_R;
	
//	float ax = 0, ay = 0;
//	si n = 0;
	si i = 0;
	
	si	(*edge)	(tEye* peye, si i, float t, float* px, float* py);
	switch (peye->Fit) {
	case eEye_Fit_S3Fit_Point:	edge = Eye_S3Fit_EdgeMark;	break;
	default:
	case eEye_Fit_S3Fit_Eye:	edge = Eye_S3Fit_EdgeMark2;	break;
	case eEye_Fit_S3Fit_Eye3:	edge = Eye_S3Fit_EdgeMark3;	break;
	
	case eEye_Fit_S3Fit_Point4:	edge = Eye_S3Fit_EdgeMark4;	break;	//soft thresholding
	}
	
	float inc = peye->AngRes*deg2rad;
	
	float ao, a;
//	for (a = 0; a < M_PI; a += inc) {
//	for (a = -M_PI_2; a < M_PI_2; a += inc) {
//	for (a = -M_PI_4; a < M_PI_4; a += inc) {
	for (ao = 0; ao < M_PI_2; ao += inc) {
	for (a = ao; a <= ao+M_PI_2+M_PI_4; a += M_PI_2) {
//	{
	//	printf ("mino ");
		u08 skip;
		float x0, y0, x1, y1;
		float t0 = a, t1 = t0 - M_PI;
		
	//	skip = edge (peye, i, t0, &x0, &y0);
		skip = edge (peye, a/(inc), t0, &x0, &y0);
		
		if (skip == 0) {
			t0 = t1;
		//	skip = edge (peye, 110 + i, t0, &x0, &y0);
			skip = edge (peye, 110 + a/(inc), t0, &x0, &y0);
		}else {
		//	skip += edge (peye, 110 + i, t1, &x1, &y1);
			skip += edge (peye, 110 + a/(inc), t1, &x1, &y1);
		}
		++i;
		
		if (skip == 0)
			continue;
		
		if (skip == 1) {
			dnpix(x0,y0)->Y = 0xFF;
			dnpix(x0,y0)->U = 0xF;
			dnpix(x0,y0)->V = 0xF;
			
			float t, dx, dy, xx, yy;
			dx = x0 - dcam.P.x;
			dy = y0 - dcam.P.y;
			t = atan2(dy,dx)/* - peye->Aa*/;
			
			xx = (peye->Ax * cos(t) * cos(peye->Aa) - peye->Ay * sin(t) * sin(peye->Aa));
			yy = (peye->Ax * cos(t) * sin(peye->Aa) + peye->Ay * sin(t) * cos(peye->Aa));
			
			float diff = sqrt(dx*dx+dy*dy) - sqrt(xx*xx+yy*yy);
			if (fabs(diff) >= 0.1f) {
			//	printf ("f1\n");
				peye->Ax += peye->Fit_Scale*diff*fabs(cos(t));
				peye->Ay += peye->Fit_Scale*diff*fabs(sin(t));
				
				Eye_Xset (peye, dcam.P.x + peye->Fit_Trans*(diff )*cos(t));
				Eye_Yset (peye, dcam.P.y + peye->Fit_Trans*(diff )*sin(t));
			}
			continue;
		}
		dnpix(x0,y0)->Y = 0xFF;
		dnpix(x0,y0)->U = 0x0;
		dnpix(x0,y0)->V = 0x0;
		
		dnpix(x1,y1)->Y = 0xFF;
		dnpix(x1,y1)->U = 0x0;
		dnpix(x1,y1)->V = 0x0;
		
		Eye_CirView_Point_xy (peye, dcam.P.x - x0, dcam.P.y - y0, 0x00FF00);
		Eye_CirView_Point_xy (peye, dcam.P.x - x1, dcam.P.y - y1, 0x00FF00);
	//	printf ("Heee? eye %ld %ld   %ld %ld\n", x0 - (si)dcam.P.x, y0 - (si)dcam.P.y, x1 - (si)dcam.P.x, y1 - (si)dcam.P.y);
	/*	dnpix(x0,y0)->Y = 0x0;
		dnpix(x0,y0)->U = 0x0;
		dnpix(x0,y0)->V = 0x0;
		
		dnpix(x1,y1)->Y = 0xFF;
		dnpix(x1,y1)->U = 0xF;
		dnpix(x1,y1)->V = 0xF;/**/
		if (0) {
			float dx1, dy1;
			dx1 = x1 - dcam.P.x;
			dy1 = y1 - dcam.P.y;
			
			float ox0 = (peye->Ax * cos(t0) * cos(peye->Aa) - peye->Ay * sin(t0) * sin(peye->Aa));
			float oy0 = (peye->Ax * cos(t0) * sin(peye->Aa) + peye->Ay * sin(t0) * cos(peye->Aa));
			
			float ox1 = (peye->Ax * cos(t1) * cos(peye->Aa) - peye->Ay * sin(t1) * sin(peye->Aa));
			float oy1 = (peye->Ax * cos(t1) * sin(peye->Aa) + peye->Ay * sin(t1) * cos(peye->Aa));
			
		/*	float diff = sqrt(dx0*dx0+dy0*dy0) - sqrt(ox0*ox0+oy0*oy0);
			if (fabs(diff) >= 0.01f) {
			//	printf ("f1\n");
				peye->Ax += peye->Fit_Scale*diff*fabs(cos(t0));
				peye->Ay += peye->Fit_Scale*diff*fabs(sin(t0));
				
				Eye_Xset (peye, dcam.P.x + peye->Fit_Trans*(diff )*cos(t0));
				Eye_Yset (peye, dcam.P.y + peye->Fit_Trans*(diff )*sin(t0));
			}/**/
		/*	float diff = sqrt(dx1*dx1+dy1*dy1) - sqrt(ox1*ox1+oy1*oy1);
			if (fabs(diff) >= 0.01f) {
			//	printf ("f1\n");
				peye->Ax += peye->Fit_Scale*diff*fabs(cos(t1));
				peye->Ay += peye->Fit_Scale*diff*fabs(sin(t1));
				
				Eye_Xset (peye, dcam.P.x + peye->Fit_Trans*(diff )*cos(t1));
				Eye_Yset (peye, dcam.P.y + peye->Fit_Trans*(diff )*sin(t1));
			}/**/
		}else {
			float ox0 = (peye->Ax * cos(t0) * cos(peye->Aa) - peye->Ay * sin(t0) * sin(peye->Aa));
			float ox1 = (peye->Ax * cos(t1) * cos(peye->Aa) - peye->Ay * sin(t1) * sin(peye->Aa));
			
			float oy0 = (peye->Ax * cos(t0) * sin(peye->Aa) + peye->Ay * sin(t0) * cos(peye->Aa));
			float oy1 = (peye->Ax * cos(t1) * sin(peye->Aa) + peye->Ay * sin(t1) * cos(peye->Aa));
			
		//	printf ("Heee? pic %f %f   %f %f\n", x0 - dcam.P.x, y0 - dcam.P.y, x1 - dcam.P.x, y1 - dcam.P.y);
		//	printf ("Heee? eli %f %f   %f %f\n", ox0, oy0, ox1, oy1);
			
			float sdx = fabsf(x1 - x0) - fabsf(ox1 - ox0);
			float sdy = fabsf(y1 - y0) - fabsf(oy1 - oy0);
			
		//	ax += fabsf((x1 - x0)*cos(t0));
		//	ay += fabsf((y1 - y0)*sin(t0));
		//	++n;
			
			float pdx = (x0+x1)*0.5f - (2.0f*dcam.P.x + ox0+ox1)*0.5f;
			float pdy = (y0+y1)*0.5f - (2.0f*dcam.P.y + oy0+oy1)*0.5f;
			
			Eye_Xset (peye, dcam.P.x + 0.10f * pdx );// *cos(t));
			Eye_Yset (peye, dcam.P.y + 0.10f * pdy );// *sin(t));
			
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
			
			Eye_Xset (peye, dcam.P.x + peye->Fit_Trans * pdx );// *cos(t));
			Eye_Yset (peye, dcam.P.y + peye->Fit_Trans * pdy );// *sin(t));
			/**/
			
		//	Eye_Xset (peye, dcam.P.x + peye->Fit_Trans * ( (x0+x1)*0.5f - dcam.P.x-(ox0+ox1)*0.5f ) );// *cos(t));
		//	Eye_Yset (peye, dcam.P.y + peye->Fit_Trans * ( (y0+y1)*0.5f - dcam.P.y-(oy0+oy1)*0.5f ) );// *sin(t));
			
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


si	Eye_S4_EdgeMark00	(tEye* peye, float t, float* px, float* py)
{
	float border_s = peye->Min_R;
	float border_e = peye->Max_R;
	
	u08 wrote = 0;
	float x = dcam.P.x, y = dcam.P.y;
	si n = 0;
	
	for (; ddist2(x,y,dcam.P.x,dcam.P.y) < dpow2(border_e); ) {
	//	if (peye->LinView.x != 0)
	//		dspix(peye->LinView.x+n, peye->LinView.y+i) = dnpix(x,y)->Y | dnpix(x,y)->Y<<8 | dnpix(x,y)->Y<<16;
		
		if (
		//	dopix(x,y)->Y > ay + 30
			dopix(x,y)->Y >= peye->S4_Pix_Bright
			&& ddist2(x,y,dcam.P.x,dcam.P.y) >= dpow2(border_s)
		) {
		//	return 1;
			float x1 = x, y1 = y;
			for (; ddist2(x1,y1,dcam.P.x,dcam.P.y) < dpow2(peye->Exp_R*2); ) {
			//	if (	dopix(x1,y1)->Y <= ay + 8) {
				if (	dopix(x1,y1)->Y <= peye->S4_Pix_Bright + 8) {
					break;
				}
				x1 += cos(t);
				y1 += sin(t);
			}/**/
		//	if (ddist2(x1,y1,x,y) >= dpow2(2)) {
				if (!wrote) {
					wrote = 1;
					*px = x+0.000000000001f;
					*py = y+0.000000000001f;
				//	dspix(peye->LinView.x+n, peye->LinView.y+i) = 0xFF<<8;
				}
			//	dspix(peye->LinView.x+n, peye->LinView.y+i) = 0xFF<<8;
		//	}
		}/**/
		
		x += cos(t);
		y += sin(t);
		++n;
	}/**/
	return wrote;
}
si	Eye_S4_EdgeMark0	(tEye* peye, float t, float* px, float* py)
{
//	si border = peye->Exp_R*peye->Exp_R*3;
	u08 wrote = 0;
	float x = dcam.P.x, y = dcam.P.y, mindiff = peye->Exp_R*peye->Exp_R;
	si n = 0, minn = 0;
	
	for (; ddist2(x,y,dcam.P.x,dcam.P.y) < peye->Exp_R*peye->Exp_R; ) {
		
		if (	dopix(x,y)->Y > ay + 36
		//	&& ddist2(x,y,dcam.P.x,dcam.P.y) >= 14*14
		) {
		//	return 1;
			float x1 = x, y1 = y;
			for (; ddist2(x1,y1,dcam.P.x,dcam.P.y) < peye->Exp_R*peye->Exp_R*4; ) {
				if (	dopix(x1,y1)->Y <= ay + 32) {
					break;
				}
				x1 += cos(t);
				y1 += sin(t);
			}
			if (ddist2(x1,y1,x,y) >= 3*2) {
			//	if (!wrote) {
					float diff = ddist2(x,y,dcam.P.x,dcam.P.y) - peye->Exp_R*peye->Exp_R;
					if (fabsf(diff) < mindiff) {
						wrote = 1;
						mindiff = fabsf(diff);
						minn = n;
						*px = x+0.000000000001f;
						*py = y+0.000000000001f;
					}
				//	dspix(peye->LinView.x+n, peye->LinView.y+i) = 0xFF<<8;
			//	}
			//	dspix(peye->LinView.x+n, peye->LinView.y+i) = 0xFF<<8;
				x = x1;
				y = y1;
				continue;
			}
		}/**/
		x += cos(t);
		y += sin(t);
		++n;
	}/**/
	if (wrote) {
	//	dspix(peye->LinView.x+peye->Exp_R, peye->LinView.y+ (si)(angle_norm_0_2pi(t)/(M_PI/100.0f))) = 0xFF<<16;
		Eye_Ellipse2LinDraw_Pix_ad (peye, t,peye->Exp_R, 0xFF<<16);
	//	dspix(peye->LinView.x+minn, peye->LinView.y+ (si)(angle_norm_0_2pi(t)/(M_PI/100.0f))) = 0xFF<<8;
		Eye_Ellipse2LinDraw_Pix_ad (peye, t,minn, 0xFF<<8);
	}
	return wrote;
	return 0;
}
si	Eye_S4_EdgeMark1	(tEye* peye, float t, float* px, float* py)
{
//	si border = peye->Exp_R*peye->Exp_R*3;
	u08 wrote = 0;
	float x = dcam.P.x, y = dcam.P.y, mindiff = peye->Exp_R*peye->Exp_R;
	si n = 0, minn = 0;
	
	tPix prev = *dopix(x,y);
	
	for (; ddist2(x,y,dcam.P.x,dcam.P.y) < peye->Exp_R*peye->Exp_R; )
	{
		if (	dopix(x,y)->Y - prev.Y >= 0x20
			&& ddist(x,y,dcam.P.x,dcam.P.y) > peye->Exp_R*0.6f
		) {
		//	Eye_Ellipse2LinDraw_Pix_ad (peye, t, n, 0xFF<<0);
			
			float diff = ddist2(x,y,dcam.P.x,dcam.P.y) - peye->Exp_R*peye->Exp_R;
			if (fabsf(diff) < mindiff) {
				wrote = 1;
				mindiff = fabsf(diff);
				minn = n;
				*px = x+0.000000000001f;
				*py = y+0.000000000001f;
			}
		}
		prev = *dopix(x,y);
		x += cos(t);
		y += sin(t);
		++n;
	}/**/
	if (wrote) {
		Eye_Ellipse2LinDraw_Pix_ad (peye, t,peye->Exp_R*0.6f, 0xFF<<16);
		Eye_Ellipse2LinDraw_Pix_ad (peye, t,peye->Exp_R, 0xFF<<16);
		Eye_Ellipse2LinDraw_Pix_ad (peye, t,minn, 0xFF<<8);
	}
	return wrote;
	return 0;
}
si	Eye_S4_EdgeMark2	(tEye* peye, float t, float* px, float* py)
{
	float border_s = peye->Min_R;
	float border_e = peye->Max_R;
	
	u08 wrote = 0;
	float x = dcam.P.x, y = dcam.P.y;
	si n = 0;
	struct {
		float x, y;
		tPix pix;
	}min, max;
	
	tPix prev = *dopix(x,y);
	
	n = border_s;
	x += border_s * cos(t);
	y += border_s * sin(t);
	
	min.x = max.x = x;
	min.y = max.y = y;
	min.pix = max.pix = *dopix(x,y);
	
	for (; ddist2(x,y,dcam.P.x,dcam.P.y) < border_e*border_e; )
	{
		if (	dopix(x,y)->Y < min.pix.Y
		//	&& ddist(x,y,dcam.P.x,dcam.P.y) > peye->Exp_R*0.6f
		) {
			min.x = x;
			min.y = y;
			min.pix = *dopix(x,y);
		}
		if (	dopix(x,y)->Y > max.pix.Y
		//	&& ddist(x,y,dcam.P.x,dcam.P.y) > peye->Exp_R*0.6f
		) {
			max.x = x;
			max.y = y;
			max.pix = *dopix(x,y);
		}
		prev = *dopix(x,y);
		x += cos(t);
		y += sin(t);
		++n;
	}/**/
	Eye_Ellipse2LinDraw_Pix_ad (peye, t,border_s, 0xFF<<0);
	Eye_Ellipse2LinDraw_Pix_ad (peye, t,peye->Exp_R, 0xFF<<16);
	Eye_Ellipse2LinDraw_Pix_ad (peye, t,border_e, 0xFF<<0);
	
	if (ddist(min.x,min.y, dcam.P.x,dcam.P.y) < ddist(max.x,max.y, dcam.P.x,dcam.P.y)) {
		if (max.pix.Y - min.pix.Y < 0x30)
			return 0;
		
		*px = min.x + max.x;	*px /= 2.0f;
		*py = min.y + max.y;	*py /= 2.0f;
		Eye_Ellipse2LinDraw_Pix_ad (peye, t, ddist(*px,*py, dcam.P.x,dcam.P.y), 0xFF<<8);
		return 1;
	}
//	if (wrote) {
		
	//	Eye_Ellipse2LinDraw_Pix_ad (peye, t, ddist(min.x,min.y, dcam.P.x,dcam.P.y), 0xFF<<0);
	//	Eye_Ellipse2LinDraw_Pix_ad (peye, t, ddist(max.x,max.y, dcam.P.x,dcam.P.y), 0xFF<<8);
//	}
//	return wrote;
	return 0;
}
si	Eye_S4_EdgeMark3	(tEye* peye, float t, float* px, float* py)
{
	float border_s = peye->Min_R;
	float border_e = peye->Max_R;
	
	si ret = Eye_S4_EdgeMark00 (peye, t, px, py);
	if (!ret) {
		return ret;
	}
	
	Eye_Ellipse2LinDraw_Pix_ad (peye, t,border_s, 0xFF<<0);
//	Eye_Ellipse2LinDraw_Pix_ad (peye, t,peye->Exp_R, 0xFF<<16);
	Eye_Ellipse2LinDraw_Pix_ad (peye, t,border_e, 0xFF<<0);
	
	if (!gM.bEye_S4_EdgeMark3_Micro) {
		Eye_Ellipse2LinDraw_Pix_ad (peye, t, ddist(*px,*py, dcam.P.x,dcam.P.y), 0xFF<<8);
		return ret;
	}
	float x0 = *px, y0 = *py;
	float x2 = *px, y2 = *py;
	
	x0 -= 1*cos(t);	y0 -= 1*sin(t);
	x2 += 1*cos(t);	y2 += 1*sin(t);
	if (*px == x0 && *py == y0) {
		printf ("crap 0\n");
	}
	if (*px == x2 && *py == y2) {
		printf ("crap 2\n");
	}
	
	float d10 = dopix(*px,*py)->Y - dopix(x0,y0)->Y;
	float d21 = dopix(x2,y2)->Y - dopix(*px,*py)->Y;
	
	if (d10 < 0)
		d10 = 0;
	if (d21 < 0)
		d21 = 0;
	
	float sd = (d21 - d10) / 100.0f;
//	printf ("sd %f\n", sd);
	*px += sd * cos(t);
	*py += sd * sin(t);
	
	Eye_Ellipse2LinDraw_Pix_ad (peye, t, ddist(*px,*py, dcam.P.x,dcam.P.y), 0xFF<<8);
	return ret;
	
	
	u08 wrote = 0;
	float x = dcam.P.x, y = dcam.P.y;
	si n = 0;
	struct {
		float x, y;
		tPix pix;
	}min, max;
	
	tPix prev = *dopix(x,y);
	
	n = border_s;
	x += border_s * cos(t);
	y += border_s * sin(t);
	
	min.x = max.x = x;
	min.y = max.y = y;
	min.pix = max.pix = *dopix(x,y);
	
	for (; ddist2(x,y,dcam.P.x,dcam.P.y) < border_e*border_e; )
	{
		if (	dopix(x,y)->Y < min.pix.Y
		//	&& ddist(x,y,dcam.P.x,dcam.P.y) > peye->Exp_R*0.6f
		) {
			min.x = x;
			min.y = y;
			min.pix = *dopix(x,y);
		}
		if (	dopix(x,y)->Y > max.pix.Y
		//	&& ddist(x,y,dcam.P.x,dcam.P.y) > peye->Exp_R*0.6f
		) {
			max.x = x;
			max.y = y;
			max.pix = *dopix(x,y);
		}
		prev = *dopix(x,y);
		x += cos(t);
		y += sin(t);
		++n;
	}/**/
	Eye_Ellipse2LinDraw_Pix_ad (peye, t,border_s, 0xFF<<0);
	Eye_Ellipse2LinDraw_Pix_ad (peye, t,peye->Exp_R, 0xFF<<16);
	Eye_Ellipse2LinDraw_Pix_ad (peye, t,border_e, 0xFF<<0);
	
	if (ddist(min.x,min.y, dcam.P.x,dcam.P.y) < ddist(max.x,max.y, dcam.P.x,dcam.P.y)) {
		if (max.pix.Y - min.pix.Y < 0x10)
			return 0;
		
		*px = min.x + max.x;	*px /= 2.0f;
		*py = min.y + max.y;	*py /= 2.0f;
		Eye_Ellipse2LinDraw_Pix_ad (peye, t, ddist(*px,*py, dcam.P.x,dcam.P.y), 0xFF<<8);
		return 1;
	}
//	if (wrote) {
		
	//	Eye_Ellipse2LinDraw_Pix_ad (peye, t, ddist(min.x,min.y, dcam.P.x,dcam.P.y), 0xFF<<0);
	//	Eye_Ellipse2LinDraw_Pix_ad (peye, t, ddist(max.x,max.y, dcam.P.x,dcam.P.y), 0xFF<<8);
//	}
//	return wrote;
	return 0;
}

#if 0
si	Eye_S4_Edge_Line	(tEye* peye, float as, float inc, float* pae, float* pr)
{
	si	(*edgemark)	(tEye* peye, float t, float* px, float* py);
	
	switch (peye->Fit) {
	default:
	case eEye_Fit_S4Fit_Edge0:	edgemark = Eye_S4_EdgeMark0;	break;
	case eEye_Fit_S4Fit_Edge1:	edgemark = Eye_S4_EdgeMark1;	break;
	case eEye_Fit_S4Fit_Edge2:	edgemark = Eye_S4_EdgeMark2;	break;
	}
	
	
	float a;//, inc = M_PI/100.0f;
	*pr = 0;
	si num = 0;
	for (a = as; ; a += inc) {
		if (inc > 0) {
			if (a > as+2*M_PI)
				break;
		}else if (a < as-2*M_PI)
			break;
		
		u08 skip;
		float x0, y0, x1, y1;
		float t0 = a, t1 = t0 + inc;
		
		if (!edgemark (peye, t0, &x0, &y0))
			continue;
		
		si tries = 11;
		do {
			skip = edgemark (peye, t1, &x1, &y1);
		}while (--tries && (!skip || x1 != x0 && y1 != y0));
	//	if (n == 0)
	//		break;
		if (!skip)
			continue;
		
		if (ddist2(x0,y0,x1,y1) < 2*2) {
		//	dnpix(x0,y0)->Y = 0xFF;
		//	dnpix(x0,y0)->U = 0x0;
		//	dnpix(x0,y0)->V = 0x0;
			
			if (peye->Point_N < peye->Point_Max) {
				peye->paPoint[peye->Point_N].x = x0;
				peye->paPoint[peye->Point_N].y = y0;
				peye->Point_N++;
				*pr += ddist(dcam.P.x, dcam.P.y, x0, y0);
				++num;
			}else {
				printf ("Well shit... i guess it's not enough\n");
			}
			
		//	dnpix(x1,y1)->Y = 0xFF;
		//	dnpix(x1,y1)->U = 0x0;
		//	dnpix(x1,y1)->V = 0x0;
		}else
			break;
	}
	*pae = a;
	*pr /= (float)num;
	return num;
}
#endif
#if 1
si	Eye_S4_Edge_Line	(tEye* peye, float as, float inc, float* pae, float* pr)
{
	si	(*edgemark)	(tEye* peye, float t, float* px, float* py);
	
	switch (peye->Fit) {
	case eEye_Fit_S4Fit_Edge0:	edgemark = Eye_S4_EdgeMark00;	break;
	case eEye_Fit_S4Fit_Edge1:	edgemark = Eye_S4_EdgeMark1;	break;
//	case eEye_Fit_S4Fit_Edge2:	edgemark = Eye_S4_EdgeMark2;	break;
	default:
	case eEye_Fit_S4Fit_Edge2:	edgemark = Eye_S4_EdgeMark3;	break;
	}
	
	tV2f prev = {0,0}, p;
	float a;
	*pr = 0;
	si num = 0;
	for (a = as; ; a += inc) {
		if (inc > 0) {
			if (a > as+2*M_PI)
				break;
			if (a > *pae)
				break;
		}else {
			if (a < as-2*M_PI) 
				break;
			if (a < *pae)
				break;
		}
	//	printf ("iter\n");
		
		if (!edgemark (peye, a, &p.x, &p.y))
			continue;
		if (prev.x == 0) {
			prev = p;
			continue;
		}
		if (p.x == prev.x && p.y == prev.y) {
			continue;
		}
		
		{
			si du = dopix(prev.x,prev.y)->U - dopix(p.x,p.y)->U;
			si dv = dopix(prev.x,prev.y)->V - dopix(p.x,p.y)->V;
		//	printf ("U %d %d\tV %d %d\n", dopix(p.x,p.y)->U, dopix(prev.x,prev.y)->U, dopix(p.x,p.y)->V, dopix(prev.x,prev.y)->V);
		//	if (abs(du) > 14 || abs(dv) > 14)
		//		continue;
		//	if (abs(du) > 8 && abs(dv) > 8)
		//		continue;
		}/**/
		if (fabsf(ddist(p.x,p.y, dcam.P.x, dcam.P.y) - ddist(prev.x,prev.y, dcam.P.x, dcam.P.y)) > (2.5f))
			continue;
		if (ddist2(p.x,p.y,prev.x,prev.y) > dpow2(3.5f))
			continue;
		
		//	dnpix(x0,y0)->Y = 0xFF;
		//	dnpix(x0,y0)->U = 0x0;
		//	dnpix(x0,y0)->V = 0x0;
			
			if (peye->Point_N < peye->Point_Max) {
				peye->paPoint[peye->Point_N].x = p.x;
				peye->paPoint[peye->Point_N].y = p.y;
				peye->Point_N++;
				*pr += ddist(dcam.P.x, dcam.P.y, p.x, p.y);
				++num;
			}else {
				printf ("Well shit... i guess it's not enough\n");
			}
			
			prev = p;
		//	dnpix(x1,y1)->Y = 0xFF;
		//	dnpix(x1,y1)->U = 0x0;
		//	dnpix(x1,y1)->V = 0x0;
	//	}else
	//		break;
	}
	*pae = a;
	*pr /= (float)num;
	return num;
}
#endif

void	Eye_S4_Edge		(tEye* peye)
{
	float ae;
	ui iter_start[4];
	ui iter_n[4];
	float r[4];
	
	peye->Point_N = 0;
//	printf ("Eye_S4_Edge:\n");
	si i = 0;
	float inc = peye->AngRes*deg2rad;//REMEMBER Y IS DOWN SO THE ANGLE GOES CLOCKWISE
	
/*	iter_start[i] = peye->Point_N;	ae = 0 + M_PI_2;		iter_n[i] = Eye_S4_Edge_Line (peye, 0, inc, &ae, &r[i]);		++i;
	iter_start[i] = peye->Point_N;	ae = 0 - M_PI_4;		iter_n[i] = Eye_S4_Edge_Line (peye, 0, -inc, &ae, &r[i]);		++i;
	
	iter_start[i] = peye->Point_N;	ae = M_PI - M_PI_2;	iter_n[i] = Eye_S4_Edge_Line (peye, M_PI, -inc, &ae, &r[i]);	++i;
	iter_start[i] = peye->Point_N;	ae = M_PI + M_PI_4;	iter_n[i] = Eye_S4_Edge_Line (peye, M_PI, inc, &ae, &r[i]);		++i;
	/**/
	iter_start[i] = peye->Point_N;	ae = 0 + M_PI_2;		iter_n[i] = Eye_S4_Edge_Line (peye, M_PI_4/2, inc, &ae, &r[i]);		++i;
	iter_start[i] = peye->Point_N;	ae = 0 - M_PI_4;		iter_n[i] = Eye_S4_Edge_Line (peye, M_PI_4/2, -inc, &ae, &r[i]);		++i;
	
	iter_start[i] = peye->Point_N;	ae = M_PI - M_PI_2;	iter_n[i] = Eye_S4_Edge_Line (peye, M_PI-M_PI_4/2, -inc, &ae, &r[i]);	++i;
	iter_start[i] = peye->Point_N;	ae = M_PI + M_PI_4;	iter_n[i] = Eye_S4_Edge_Line (peye, M_PI-M_PI_4/2, inc, &ae, &r[i]);	++i;
	/**/
	
//	iter_start[i] = peye->Point_N;	iter_n[i] = Eye_S4_Edge_Line (peye, -M_PI_2, inc, &ae, &r[i]);		++i;
//	iter_start[i] = peye->Point_N;	iter_n[i] = Eye_S4_Edge_Line (peye, -M_PI_2, -inc, &ae, &r[i]);		++i;
	
//	iter_start[i] = peye->Point_N;	iter_n[i] = Eye_S4_Edge_Line (peye, M_PI_2, inc, &ae, &r[i]);		++i;
//	iter_start[i] = peye->Point_N;	iter_n[i] = Eye_S4_Edge_Line (peye, M_PI_2, -inc, &ae, &r[i]);		++i;
	
	si num = i;
	si n = 0;
	float ar = 0;
	for (i = 0; i < num; ++i) {
	//	printf ("\t%ld: r %f\n", iter_n[i], r[i]);
		if (iter_n[i]) {
			ar += iter_n[i] * r[i];
			n += iter_n[i];
		}
	}
	ar /= n;
	
//	printf ("Eye_S4_Edge: remove from %ld:\n", peye->Point_N);
/*	for (i = num-1; i >= 0; --i) {
		if (iter_n[i]) {
			if (fabsf(r[i] - ar) > 0.1f*ar
				|| iter_n[i] < 20
			) {
				printf ("\t%ld: %ld r %f\n", i, iter_n[i], r[i]);
				memmove (	peye->paPoint + iter_start[i],
						peye->paPoint + iter_start[i] + iter_n[i],
						(peye->Point_N - (iter_start[i] + iter_n[i])) * sizeof(float)
				);
				peye->Point_N -= iter_n[i];
			}
			++n;
		}
	}/**/
	
//	Eye_S4_Edge_Line (peye, M_PI, -M_PI/100.0f, &ae, &r);
	
/*	for (a = ao; a <= ao+M_PI_2+M_PI_4; a += M_PI_2) {
//	{
	//	printf ("mino ");
		u08 skip;
		float x0, y0, x1, y1;
		float t0 = a, t1 = t0 - M_PI;
		
	//	skip = edge (peye, i, t0, &x0, &y0);
		skip = edge (peye, a/(M_PI/100.0f), t0, &x0, &y0);
		
		if (skip == 0) {
			t0 = t1;
		//	skip = edge (peye, 110 + i, t0, &x0, &y0);
			skip = edge (peye, 110 + a/(M_PI/100.0f), t0, &x0, &y0);
		}else {
		//	skip += edge (peye, 110 + i, t1, &x1, &y1);
			skip += edge (peye, 110 + a/(M_PI/100.0f), t1, &x1, &y1);
		}
		++i;
		
		if (skip == 0)
			continue;
		
	}
	}/**/
	
	return;
}


si	Eye_S4_Fit_FindOposite	(tEye* peye, si idx)
{
	float x0 = peye->paPoint[idx].x, y0 = peye->paPoint[idx].y;
	float t0 = atan2(y0 - dcam.P.y, x0 - dcam.P.x);
	
	struct {
		si idx;
		float diff_t;
	}min = {-1, deg2rad};
	si i;
	for (i = 0; i < peye->Point_N; ++i) {
		float x1 = peye->paPoint[i].x, y1 = peye->paPoint[i].y;
		float t1 = atan2(y1 - dcam.P.y, x1 - dcam.P.x);
		
		float diff_t = fabsf(fabsf(angle_norm_0_2pi(t1) - angle_norm_0_2pi(t0)) - M_PI);
		if (diff_t < min.diff_t) {
			min.idx = i;
			min.diff_t = diff_t;
		}
	}
	return min.idx;
}

si	Eye_S4_Fit_FindRot	(tEye* peye, si idx, float ang)
{
	float x0 = peye->paPoint[idx].x, y0 = peye->paPoint[idx].y;
	float t0 = atan2(y0 - dcam.P.y, x0 - dcam.P.x);
	
	struct {
		si idx;
		float diff_t;
	}min = {-1, 4*deg2rad};
//	printf ("Eye_S4_Fit_FindRot\n");
	si i;
	for (i = 0; i < peye->Point_N; ++i) {
		float x1 = peye->paPoint[i].x, y1 = peye->paPoint[i].y;
		float t1 = atan2(y1 - dcam.P.y, x1 - dcam.P.x);
	//	printf ("%ld\t%f\n", i, t1);
		
		float diff_t = fabsf(angle_diff_norm_pi_pi(t1, t0) - ang);
		if (diff_t < min.diff_t) {
			min.idx = i;
			min.diff_t = diff_t;
		}
	}
	return min.idx;
}


float	Eye_S4_Fit		(tEye* peye)
{
	si i;
/*	for (i = 0; i < peye->Point_N; ++i) {
		float x = peye->paPoint[i].x, y = peye->paPoint[i].y;
		float t, dx, dy, xx, yy;
		dx = x - dcam.P.x;
		dy = y - dcam.P.y;
		t = atan2(dy,dx);
		
		xx = (peye->Ax * cos(t) * cos(peye->Aa) - peye->Ay * sin(t) * sin(peye->Aa));
		yy = (peye->Ax * cos(t) * sin(peye->Aa) + peye->Ay * sin(t) * cos(peye->Aa));
		
		float diff = sqrt(dx*dx+dy*dy) - sqrt(xx*xx+yy*yy);
		if (fabs(diff) >= 0.01f) {
		//	printf ("f1\n");
			peye->Ax += peye->Fit_Scale*diff*fabs(cos(t));
			peye->Ay += peye->Fit_Scale*diff*fabs(sin(t));
			
			Eye_Xset (peye, dcam.P.x + peye->Fit_Trans*(diff )*cos(t));
			Eye_Yset (peye, dcam.P.y + peye->Fit_Trans*(diff )*sin(t));
		}
		dnpix(x,y)->Y = 0xFF;
		dnpix(x,y)->U = 0x0;
		dnpix(x,y)->V = 0x0;
	}/**/
/*	for (i = peye->Point_N-1; i >= 0; --i) {
		float x = peye->paPoint[i].x, y = peye->paPoint[i].y;
		float t, dx, dy, xx, yy;
		dx = x - dcam.P.x;
		dy = y - dcam.P.y;
		t = atan2(dy,dx);
		
		xx = (peye->Ax * cos(t) * cos(peye->Aa) - peye->Ay * sin(t) * sin(peye->Aa));
		yy = (peye->Ax * cos(t) * sin(peye->Aa) + peye->Ay * sin(t) * cos(peye->Aa));
		
		float diff = sqrt(dx*dx+dy*dy) - sqrt(xx*xx+yy*yy);
		if (fabs(diff) >= 0.01f) {
		//	printf ("f1\n");
			peye->Ax += peye->Fit_Scale*diff*fabs(cos(t));
			peye->Ay += peye->Fit_Scale*diff*fabs(sin(t));
			
			Eye_Xset (peye, dcam.P.x + peye->Fit_Trans*(diff )*cos(t));
			Eye_Yset (peye, dcam.P.y + peye->Fit_Trans*(diff )*sin(t));
		}
		dnpix(x,y)->Y = 0xFF;
		dnpix(x,y)->U = 0x0;
		dnpix(x,y)->V = 0x0;
	}/**/
	for (i = 0; i < peye->Point_N; ++i) {
		float x0 = peye->paPoint[i].x, y0 = peye->paPoint[i].y;
		float t0, dx, dy, ox0, oy0;
		dx = x0 - dcam.P.x;
		dy = y0 - dcam.P.y;
		t0 = atan2(dy,dx);
		
		ox0 = (peye->Ax * cos(t0) * cos(peye->Aa) - peye->Ay * sin(t0) * sin(peye->Aa));
		oy0 = (peye->Ax * cos(t0) * sin(peye->Aa) + peye->Ay * sin(t0) * cos(peye->Aa));
		
		{
			si idx = Eye_S4_Fit_FindOposite (peye, i);
			if (idx >= 0) {
				float x1 = peye->paPoint[idx].x, y1 = peye->paPoint[idx].y;
				float t1 = atan2(y1 - dcam.P.y, x1 - dcam.P.x);
				
				dnpix(x0,y0)->Y = 0xFF;
				dnpix(x0,y0)->U = 0x0;
				dnpix(x0,y0)->V = 0x0;
			//	dnpix(x1,y1)->Y = 0xFF;
			//	dnpix(x1,y1)->U = 0x0;
			//	dnpix(x1,y1)->V = 0x0;
				
				
			//	float ox0 = (peye->Ax * cos(t0) * cos(peye->Aa) - peye->Ay * sin(t0) * sin(peye->Aa));
			//	float oy0 = (peye->Ax * cos(t0) * sin(peye->Aa) + peye->Ay * sin(t0) * cos(peye->Aa));
				
				float oy1 = (peye->Ax * cos(t1) * sin(peye->Aa) + peye->Ay * sin(t1) * cos(peye->Aa));
				float ox1 = (peye->Ax * cos(t1) * cos(peye->Aa) - peye->Ay * sin(t1) * sin(peye->Aa));
				
			//	printf ("Heee? pic %f %f   %f %f\n", x0 - dcam.P.x, y0 - dcam.P.y, x1 - dcam.P.x, y1 - dcam.P.y);
			//	printf ("Heee? eli %f %f   %f %f\n", ox0, oy0, ox1, oy1);
				
				float sdx = fabsf(x1 - x0) - fabsf(ox1 - ox0);
				float sdy = fabsf(y1 - y0) - fabsf(oy1 - oy0);
				
			//	ax += fabsf((x1 - x0)*cos(t0));
			//	ay += fabsf((y1 - y0)*sin(t0));
			//	++n;
				
				float pdx = (x0+x1)*0.5f - (2.0f*dcam.P.x + ox0+ox1)*0.5f;
				float pdy = (y0+y1)*0.5f - (2.0f*dcam.P.y + oy0+oy1)*0.5f;
				
				Eye_Xset (peye, dcam.P.x + 0.10f * pdx );// *cos(t));
				Eye_Yset (peye, dcam.P.y + 0.10f * pdy );// *sin(t));
				
				peye->Ax += 0.10f*sdx*fabs(cos(t0));
				peye->Ay += 0.10f*sdy*fabs(sin(t0));
			}else {
				dnpix(x0,y0)->Y = 0xFF;
				dnpix(x0,y0)->U = 0xF;
				dnpix(x0,y0)->V = 0xF;
				float diff = sqrt(dx*dx+dy*dy) - sqrt(ox0*ox0+oy0*oy0);
				if (fabs(diff) >= 0.01f) {
				//	printf ("f1\n");
					peye->Ax += peye->Fit_Scale*diff*fabs(cos(t0));
					peye->Ay += peye->Fit_Scale*diff*fabs(sin(t0));
					
					Eye_Xset (peye, dcam.P.x + peye->Fit_Trans*(diff )*cos(t0));
					Eye_Yset (peye, dcam.P.y + peye->Fit_Trans*(diff )*sin(t0));
				}
			}
		}
	}/**/
	float avgerr = 0;
	for (i = 0; i < peye->Point_N; ++i) {
		float x = peye->paPoint[i].x, y = peye->paPoint[i].y;
		float t, dx, dy, xx, yy;
		dx = x - dcam.P.x;
		dy = y - dcam.P.y;
		t = atan2(dy,dx);
		
		xx = (peye->Ax * cos(t) * cos(peye->Aa) - peye->Ay * sin(t) * sin(peye->Aa));
		yy = (peye->Ax * cos(t) * sin(peye->Aa) + peye->Ay * sin(t) * cos(peye->Aa));
		
		float diff = sqrt(dx*dx+dy*dy) - sqrt(xx*xx+yy*yy);
	//	avgerr += fabsf(diff);
		avgerr += dpow2(diff);
	}/**/
	return avgerr / peye->Point_N;
}

float	Eye_S4_Fit_Ransac		(tEye* peye)
{
	si i;
	gM_edge_point_N = 0;
	for (i = 0; i < peye->Point_N; ++i) {
		dnpix(peye->paPoint[i].x,peye->paPoint[i].y)->Y = 0xFF;
		dnpix(peye->paPoint[i].x,peye->paPoint[i].y)->U = 0x0;
		dnpix(peye->paPoint[i].x,peye->paPoint[i].y)->V = 0x0;
		gM_edge_point[i].x = peye->paPoint[i].x;
		gM_edge_point[i].y = peye->paPoint[i].y;
	}
	gM_edge_point_N = peye->Point_N;
	
	int max_inliers_num;
	int* ret = pupil_fitting_inliers (pcam->Image_W, pcam->Image_H, &max_inliers_num);
//	int* ret = simple_fit (videoIn->width, videoIn->height, &max_inliers_num);
	
	//double pupil_param[5];//parameters of an ellipse {ellipse_a, ellipse_b, cx, cy, theta}; a & b is the major or minor axis; 
	
//	printf ("%d: x %f y %f  Ax %f Ay %f   Aa %f\n", max_inliers_num, pupil_param[2], pupil_param[3], pupil_param[0], pupil_param[1], pupil_param[4]);
	
	float x = dcam.P.x, y = dcam.P.y, ax = peye->Ax, ay = peye->Ay, aa = peye->Aa;
	
	dcam.P.x = pupil_param[2];
	dcam.P.y = pupil_param[3];
	peye->Ax = pupil_param[0];
	peye->Ay = pupil_param[1];
	peye->Aa = pupil_param[4];
	
//	Eye_Draw_Ellipse (peye, 0.0f, M_PI*2);
	
	if (ddist(dcam.P.x, dcam.P.y, x, y) >= dpow2(6)) {
		dcam.P.x = x;
		dcam.P.y = y;
		peye->Ax = ax;
		peye->Ay = ay;
		
	}/**/
	return 0;
}

float	Eye_S4_Fit_Tri		(tEye* peye)
{
	si i;
	for (i = 0; i < peye->Point_N; ++i) {
		float x0 = peye->paPoint[i].x, y0 = peye->paPoint[i].y;
		float t0, dx, dy, ox0, oy0;
		dx = x0 - dcam.P.x;
		dy = y0 - dcam.P.y;
		t0 = atan2(dy,dx);
		
		ox0 = (peye->Ax * cos(t0) * cos(peye->Aa) - peye->Ay * sin(t0) * sin(peye->Aa));
		oy0 = (peye->Ax * cos(t0) * sin(peye->Aa) + peye->Ay * sin(t0) * cos(peye->Aa));
		
		tV2f p0 = {x0, y0}, op0 = {ox0, oy0};
		{
			si i1 = Eye_S4_Fit_FindRot (peye, i, M_PI_4/2);
			si i2 = Eye_S4_Fit_FindRot (peye, i, -M_PI_4/2);
			
		//	si i1 = Eye_Points_FindRot (peye, i, M_PI_4/2);
		//	si i2 = Eye_Points_FindRot (peye, i, -M_PI_4/2);
			
			if (i1 >= 0 && i2 >= 0) {
				float x1 = peye->paPoint[i1].x, y1 = peye->paPoint[i1].y;
				float t1 = atan2(y1 - dcam.P.y, x1 - dcam.P.x);
				
				float x2 = peye->paPoint[i2].x, y2 = peye->paPoint[i2].y;
				float t2 = atan2(y2 - dcam.P.y, x2 - dcam.P.x);
				
				
				float oy1 = (peye->Ax * cos(t1) * sin(peye->Aa) + peye->Ay * sin(t1) * cos(peye->Aa));
				float ox1 = (peye->Ax * cos(t1) * cos(peye->Aa) - peye->Ay * sin(t1) * sin(peye->Aa));
				
				float oy2 = (peye->Ax * cos(t2) * sin(peye->Aa) + peye->Ay * sin(t2) * cos(peye->Aa));
				float ox2 = (peye->Ax * cos(t2) * cos(peye->Aa) - peye->Ay * sin(t2) * sin(peye->Aa));
				
				tV2f p1 = {x1, y1}, op1 = {ox1, oy1};
				tV2f p2 = {x2, y2}, op2 = {ox2, oy2};
				
				tV2f v1 = p1;	V2f_sub_V2f (&v1, &p0);
				tV2f v2 = p2;	V2f_sub_V2f (&v2, &p0);
				
				tV2f ov1 = op1;	V2f_sub_V2f (&ov1, &op0);
				tV2f ov2 = op2;	V2f_sub_V2f (&ov2, &op0);
				
				float dot = V2f_dot_V2f (&v1, &v2);
				float odot = V2f_dot_V2f (&ov1, &ov2);
				
				float diff = odot-dot;
				peye->Ax += +0.01f*diff*fabsf(sin(t0))	-0.01f*diff*fabsf(cos(t0));
				peye->Ay += -0.01f*diff*fabsf(sin(t0))	+0.01f*diff*fabsf(cos(t0));
				
				
			/*	float sdx = fabsf(x1 - x0) - fabsf(ox1 - ox0);
				float sdy = fabsf(y1 - y0) - fabsf(oy1 - oy0);
				
				float pdx = (x0+x1)*0.5f - (2.0f*dcam.P.x + ox0+ox1)*0.5f;
				float pdy = (y0+y1)*0.5f - (2.0f*dcam.P.y + oy0+oy1)*0.5f;
				
				Eye_Xset (peye, dcam.P.x + 0.10f * pdx );// *cos(t));
				Eye_Yset (peye, dcam.P.y + 0.10f * pdy );// *sin(t));
				
				peye->Ax += 0.10f*sdx*fabs(cos(t0));
				peye->Ay += 0.10f*sdy*fabs(sin(t0));/**/
				
				dnpix(x0,y0)->Y = 0xFF;
				dnpix(x0,y0)->U = 0x0;
				dnpix(x0,y0)->V = 0x0;
			}else {
				dnpix(x0,y0)->Y = 0xFF;
				dnpix(x0,y0)->U = 0xF;
				dnpix(x0,y0)->V = 0xF;
			}
			if (1) {
				si idx = Eye_S4_Fit_FindOposite (peye, i);
				if (idx >= 0) {
					float x1 = peye->paPoint[idx].x, y1 = peye->paPoint[idx].y;
					float t1 = atan2(y1 - dcam.P.y, x1 - dcam.P.x);
					
				//	dnpix(x0,y0)->Y = 0xFF;
				//	dnpix(x0,y0)->U = 0x0;
				//	dnpix(x0,y0)->V = 0x0;
				//	dnpix(x1,y1)->Y = 0xFF;
				//	dnpix(x1,y1)->U = 0x0;
				//	dnpix(x1,y1)->V = 0x0;
					
					
				//	float ox0 = (peye->Ax * cos(t0) * cos(peye->Aa) - peye->Ay * sin(t0) * sin(peye->Aa));
				//	float oy0 = (peye->Ax * cos(t0) * sin(peye->Aa) + peye->Ay * sin(t0) * cos(peye->Aa));
					
					float oy1 = (peye->Ax * cos(t1) * sin(peye->Aa) + peye->Ay * sin(t1) * cos(peye->Aa));
					float ox1 = (peye->Ax * cos(t1) * cos(peye->Aa) - peye->Ay * sin(t1) * sin(peye->Aa));
					
				//	printf ("Heee? pic %f %f   %f %f\n", x0 - dcam.P.x, y0 - dcam.P.y, x1 - dcam.P.x, y1 - dcam.P.y);
				//	printf ("Heee? eli %f %f   %f %f\n", ox0, oy0, ox1, oy1);
					
					float sdx = fabsf(x1 - x0) - fabsf(ox1 - ox0);
					float sdy = fabsf(y1 - y0) - fabsf(oy1 - oy0);
					
				//	ax += fabsf((x1 - x0)*cos(t0));
				//	ay += fabsf((y1 - y0)*sin(t0));
				//	++n;
					
					float pdx = (x0+x1)*0.5f - (2.0f*dcam.P.x + ox0+ox1)*0.5f;
					float pdy = (y0+y1)*0.5f - (2.0f*dcam.P.y + oy0+oy1)*0.5f;
					
					Eye_Xset (peye, dcam.P.x + 0.10f * pdx );// *cos(t));
					Eye_Yset (peye, dcam.P.y + 0.10f * pdy );// *sin(t));
					
					peye->Ax += 0.10f*sdx*fabs(cos(t0));
					peye->Ay += 0.10f*sdy*fabs(sin(t0));
					continue;
				}
			}
			if (1) {
				float diff = sqrt(dx*dx+dy*dy) - sqrt(ox0*ox0+oy0*oy0);
				if (fabs(diff) >= 0.01f) {
				//	printf ("f1\n");
					peye->Ax += 0.1f*diff*fabs(cos(t0));
					peye->Ay += 0.1f*diff*fabs(sin(t0));
					
					Eye_Xset (peye, dcam.P.x + 0.1f*(diff)*cos(t0));
					Eye_Yset (peye, dcam.P.y + 0.1f*(diff)*sin(t0));
				}
			}
		}
	}/**/
	return 0;
}

void	Eye_S4		(tEye* peye)
{
	Eye_Ellipse2LinDraw (peye);
	
	if (!peye->Point_Max) {
		peye->Point_Max = 2048;
		peye->paPoint = malloc (peye->Point_Max * sizeof(tV2f));
	}
	
//	Eye_S4_Edge (peye);	Eye_S4_Fit (peye);		return;
//	Eye_S4_Edge (peye);	Eye_S4_Fit_Ransac (peye);	return;
	Eye_S4_Edge (peye);	Eye_S4_Fit_Tri (peye);		Eye_S4_Fit_Tri (peye);	return;
	
	tEye mineye, teye = *peye;
	
	Eye_S4_Edge (&teye);
	float minerr = Eye_S4_Fit (&teye);	//printf ("Error %f %f  e %f\n", teye.P.x, teye.P.y, minerr);
	mineye = teye;
	
	float dd = 2;
	float x, y, ox = dcam.P.x, oy = dcam.P.y;
	for (y = oy - dd; y <= oy + dd; y += 1) {
		for (x = ox - dd; x <= ox + dd; x += 1) {
			teye = *peye;
			teye.aCam[pcam->Idx].P.x = x;
			teye.aCam[pcam->Idx].P.y = y;
			Eye_S4_Edge (&teye);
			float avgerr = Eye_S4_Fit (&teye);	//printf ("Error %f %f  e %f\n", teye.P.x, teye.P.y, avgerr);
			if (teye.Point_N >= mineye.Point_N && avgerr < minerr) {
			//	printf ("teye.Point_N %ld  mineye %ld\n", teye.Point_N, mineye.Point_N);
				minerr = avgerr;
				mineye = teye;
			}
		}
	}/**/
	Eye_CopyParam (peye, &mineye);
	
	return;
}




si	Eye_S5_Edge_Mark00	(tEye* peye, float t, si* pnum, tV2f* papoint)
{
	float x = dcam.P.x, y = dcam.P.y;
	si n = 0;
	u08 got_dark = 0;
	
	for (; ddist2(x,y,dcam.P.x,dcam.P.y) < dpow2(peye->Exp_R*2); ) {
		if (!got_dark) {
			if (
			//	dopix(x,y)->Y < ay + 8
				dopix(x,y)->Y <= peye->S5.Pix_Dark
			) {
				Eye_Ellipse2LinDraw_Pix_ad (peye, t, ddist(x,y, dcam.P.x,dcam.P.y), 0x80<<8);
				got_dark = 1;
			}
		}else if (
		//	dopix(x,y)->Y > ay + 40
			dopix(x,y)->Y >= peye->S5.Pix_Bright
		//	&& ddist2(x,y,dcam.P.x,dcam.P.y) >= dpow2(peye->Min_R)
		) {
		//	return 1;
			float x1 = x, y1 = y;
			for (; ddist2(x1,y1,dcam.P.x,dcam.P.y) < dpow2(peye->Exp_R*2); ) {
			//	if (	dopix(x1,y1)->Y <= ay + 8) {
				if (	dopix(x1,y1)->Y <= peye->S5.Pix_Bright-2) {
					break;
				}
				x1 += cos(t);
				y1 += sin(t);
			}/**/
			if (
				n < *pnum
			//	&& ddist2(x,y, x1,y1) >= dpow2(6)
			) {
			//	if (ddist2(x1,y1,x,y) >= dpow2(2)) {
			//	Eye_Ellipse2LinDraw_Pix_ad (peye, t, ddist(x,y, dcam.P.x,dcam.P.y), 0xFF<<8);
				papoint[n].x = x+0.000000000001f;
				papoint[n].y = y+0.000000000001f;
				
				if (gM.bEye_S4_EdgeMark3_Micro) {
					float x0 = papoint[n].x, y0 = papoint[n].y;
					float x2 = papoint[n].x, y2 = papoint[n].y;
					
					x0 -= 1*cos(t);	y0 -= 1*sin(t);
					x2 += 1*cos(t);	y2 += 1*sin(t);
					if (papoint[n].x == x0 && papoint[n].y == y0) {
						printf ("crap 0\n");
					}
					if (papoint[n].x == x2 && papoint[n].y == y2) {
						printf ("crap 2\n");
					}
					
					float d10 = dopix(papoint[n].x,papoint[n].y)->Y - dopix(x0,y0)->Y;
					float d21 = dopix(x2,y2)->Y - dopix(papoint[n].x,papoint[n].y)->Y;
					
					if (d10 < 0)
						d10 = 0;
					if (d21 < 0)
						d21 = 0;
					
					float sd = (d21 - d10) / 100.0f;
				//	printf ("sd %f\n", sd);
					papoint[n].x += sd * cos(t);
					papoint[n].y += sd * sin(t);
					
					Eye_Ellipse2LinDraw_Pix_ad (peye, t, ddist(papoint[n].x,papoint[n].y, dcam.P.x,dcam.P.y), 0xFF<<8);
				}else
					Eye_Ellipse2LinDraw_Pix_ad (peye, t, ddist(papoint[n].x,papoint[n].y, dcam.P.x,dcam.P.y), 0xFF<<8);
				++n;
			}
			x = x1;
			y = y1;
			got_dark = 0;
		}/**/
		
		x += cos(t);
		y += sin(t);
	}/**/
	*pnum = n;
	return 0;
}


si	Eye_S5_Edge_Mark1		(tEye* peye, float t, si* pnum, tV2f* papoint)
{
	float x0 = dcam.P.x, y0 = dcam.P.y;
	si n = 0;
	u08 got_fail = 0;
	
	for (; ddist2(x0,y0,dcam.P.x,dcam.P.y) < dpow2(peye->Exp_R*2 - peye->S5.Diff_Dist); )
	{
		float x1 = x0 + cos(t)*peye->S5.Diff_Dist, y1 = y0 + sin(t)*peye->S5.Diff_Dist;
		
		if (dopix(x1,y1)->Y - dopix(x0,y0)->Y >= peye->S5.Pix_Diff_Start) {
			float x,y;
			x = x0; y = y0;
			while (1) {
				x0 = x; y0 = y;
				x += cos(t);
				y += sin(t);
				if (dopix(x,y)->Y - dopix(x0,y0)->Y >= peye->S5.Pix_Diff_Min) {
					if (got_fail < 1)
						++got_fail;
					else
						break;
				}
			}
			x = x1; y = y1;
			while (1) {
				x1 = x; y1 = y;
				x += cos(t);
				y += sin(t);
				if (dopix(x,y)->Y - dopix(x1,y1)->Y < peye->S5.Pix_Diff_Min) {
					if (got_fail < 1)
						++got_fail;
					else
						break;
				}
			}
			x = x0+x1;	x /= 2;
			y = y0+y1;	y /= 2;
			papoint[n].x = x+0.000000000001f;
			papoint[n].y = y+0.000000000001f;
			
			Eye_Ellipse2LinDraw_Pix_ad (peye, t, ddist(x0,y0, dcam.P.x,dcam.P.y), 0x80<<8);
			Eye_Ellipse2LinDraw_Pix_ad (peye, t, ddist(x1,y1, dcam.P.x,dcam.P.y), 0xFF<<8);
			Eye_Ellipse2LinDraw_Pix_ad (peye, t, ddist(papoint[n].x,papoint[n].y, dcam.P.x,dcam.P.y), 0xFF<<16);
			++n;
			if (n >= *pnum)
				break;
			
			x0 = x1;
			y0 = y1;
		}
		x0 += cos(t);
		y0 += sin(t);
		
	}/**/
	*pnum = n;
	return 0;
}

si	Eye_S5_Edge_Mark		(tEye* peye, float t, si* pnum, tV2f* papoint)
{
	switch (peye->Fit) {
	case eEye_Fit_S5_Diff:	return Eye_S5_Edge_Mark1 (peye, t, pnum, papoint);
	default:			return Eye_S5_Edge_Mark00 (peye, t, pnum, papoint);
	}
}

void	Eye_S5_Edge_Scan	(tEye* peye)
{
	float a = 0;
	for (a = 0; a < 2*M_PI; a += peye->AngRes*deg2rad) {
		si n = 20;
		tV2f point[20];
		
		Eye_S5_Edge_Mark (peye, a, &n, point);
		
	//	Eye_Ellipse2LinDraw_Pix_ad (peye, a,border_s, 0xFF<<0);
	//	Eye_Ellipse2LinDraw_Pix_ad (peye, a,peye->Exp_R, 0xFF<<16);
	//	Eye_Ellipse2LinDraw_Pix_ad (peye, a,border_e, 0xFF<<0);
		si i;
		for (i = 0; i < n; ++i) {
			dnpix(point[i].x,point[i].y)->Y = 0xFF;
			dnpix(point[i].x,point[i].y)->U = 0xF;
			dnpix(point[i].x,point[i].y)->V = 0xF;
		}
		
	}
}

void	Eye_S5_Edge_Trace	(tEye* peye, float as, si ns, float inc)
{
	tV2f prev_point;
	float a = as;
	{
		si n = 20;	tV2f point[20];
		Eye_S5_Edge_Mark (peye, a, &n, point);
		prev_point = point[ns];
		Eye_Points_Ins (peye, a, &prev_point);
	//	Eye_Ellipse2LinDraw_Pix_ad (peye, a, ddist(point[ns].x,point[ns].y, dcam.P.x,dcam.P.y), col);
	}
	
	for (; ; a += inc) {
		if (inc >= 0) {
			if (a >= as + M_PI)
				break;
		}else {
			if (a <= as - M_PI)
				break;
		}
		si n = 20;
		tV2f point[20];
		
		Eye_S5_Edge_Mark (peye, a, &n, point);
		
		float mind = dpow2(peye->S5.Break_Dist);
		si i, mini = -1;
		for (i = 0; i < n; ++i) {
			if (ddist2(point[i].x,point[i].y, prev_point.x,prev_point.y) < mind) {
				mind = ddist2(point[i].x,point[i].y, prev_point.x,prev_point.y);
				mini = i;
			}
		}
		if (mini >= 0) {
		//	Eye_Ellipse2LinDraw_Pix_ad (peye, a, peye->Exp_R, 0xFF<<16);
		//	Eye_Ellipse2LinDraw_Pix_ad (peye, a, ddist(point[mini].x,point[mini].y, dcam.P.x,dcam.P.y), 0xFF<<16);
			prev_point = point[mini];
			
			Eye_Points_Ins (peye, a, &prev_point);
		}
	}/**/
}

void	Eye_S5_Edges	(tEye* peye)
{
	si i, edge_n = 0, edge_max = 10;
	struct {
		ui start;
		ui num;
		float r;
	}edge[10];
	peye->Point_N = 0;
	
	if (0) {
		si n = 20;
		tV2f point[20];
	//	float a = 3*M_PI_4/2;
		float a = M_PI_2;
		Eye_Ellipse2LinDraw_Pix_ad (peye, a, peye->Exp_R, 0xFF<<8);
		Eye_Ellipse2LinDraw_Pix_ad (peye, a, peye->Exp_R+1, 0xFF<<8);
		Eye_Ellipse2LinDraw_Pix_ad (peye, a, peye->Exp_R+2, 0xFF<<8);
		Eye_Ellipse2LinDraw_Pix_ad (peye, a, peye->Exp_R+3, 0xFF<<8);
		
		Eye_S5_Edge_Mark (peye, a, &n, point);
		return;
	}
//	Eye_S5_Edge_Scan (peye);	return;
	
//	Eye_S5_Edge_Trace (peye, 0, 0);
	{
		si n = 20;
		tV2f point[20];
	//	float a = 3*M_PI_4/2;
		float a = M_PI_2;
		Eye_Ellipse2LinDraw_Pix_ad (peye, a, peye->Exp_R, 0xFF<<8);
		Eye_Ellipse2LinDraw_Pix_ad (peye, a, peye->Exp_R+1, 0xFF<<8);
		Eye_Ellipse2LinDraw_Pix_ad (peye, a, peye->Exp_R+2, 0xFF<<8);
		Eye_Ellipse2LinDraw_Pix_ad (peye, a, peye->Exp_R+3, 0xFF<<8);
		
		Eye_S5_Edge_Mark (peye, a, &n, point);
		
		for (i = 0; i < n; ++i) {
			if (edge_n >= edge_max)
				break;
		//	printf ("hello i %d\n", i);
			edge[edge_n].start = peye->Point_N;
			Eye_S5_Edge_Trace (peye, a, i, peye->AngRes*deg2rad);
			Eye_S5_Edge_Trace (peye, a, i, -peye->AngRes*deg2rad);
			edge[edge_n].num = peye->Point_N - edge[edge_n].start;
			edge_n++;
		}
	}
	if (0) {
		si n = 20;
		tV2f point[20];
		float a = -M_PI_2;
		Eye_Ellipse2LinDraw_Pix_ad (peye, a, peye->Exp_R, 0xFF<<8);
		Eye_Ellipse2LinDraw_Pix_ad (peye, a, peye->Exp_R+1, 0xFF<<8);
		Eye_Ellipse2LinDraw_Pix_ad (peye, a, peye->Exp_R+2, 0xFF<<8);
		Eye_Ellipse2LinDraw_Pix_ad (peye, a, peye->Exp_R+3, 0xFF<<8);
		
		Eye_S5_Edge_Mark (peye, a, &n, point);
		
		for (i = 0; i < n; ++i) {
			if (edge_n >= edge_max)
				break;
		//	printf ("hello i %d\n", i);
			edge[edge_n].start = peye->Point_N;
			Eye_S5_Edge_Trace (peye, a, i, peye->AngRes*deg2rad);
			Eye_S5_Edge_Trace (peye, a, i, -peye->AngRes*deg2rad);
			edge[edge_n].num = peye->Point_N - edge[edge_n].start;
			edge_n++;
		}
	}
/*	for (i = 0; i < edge_n; ++i) {
	//	printf ("\t%ld: r %f\n", edge[i].num, r[i]);
		si ip;
		for (ip = edge[i].start; ip < edge[i].start+edge[i].num; ++ip) {
			u32 col;
			switch (i) {
			case 0:	col = 0xFF<<16;				break;
			case 1:	col = 0xFF<<0;				break;
			case 2:	col = 0xFF<<16 | 0xFF<<8;		break;
			}
			Eye_Ellipse2LinDraw_Pix_ad (peye,
				angle_norm_0_2pi(atan2(peye->paPoint[ip].y-dcam.P.y,peye->paPoint[ip].x-dcam.P.x)),
				ddist(peye->paPoint[ip].x,peye->paPoint[ip].y, dcam.P.x,dcam.P.y),
				col);
			
		}
	}/**/
	
	for (i = 0; i < edge_n; ++i) {
	//	printf ("\t%ld: r %f\n", edge[i].num, r[i]);
		edge[i].r = 0;
		si ip;
		for (ip = edge[i].start; ip < edge[i].start+edge[i].num; ++ip) {
			edge[i].r += ddist(peye->paPoint[ip].x,peye->paPoint[ip].y, dcam.P.x,dcam.P.y);
		}
		edge[i].r /= (float)edge[i].num;
	}/**/
//	printf ("Eye_S5_Edges: remove from %ld:\n", peye->Point_N);
	for (i = edge_n-1; i >= 0; --i) {
		if (
			edge[i].num < peye->S5.Min_N
			|| edge[i].r < peye->S5.Min_R
		) {
		//	printf ("\t%ld: %ld\n", i, edge[i].num);
			memmove (	peye->paPoint + edge[i].start,
					peye->paPoint + edge[i].start + edge[i].num,
					(peye->Point_N - (edge[i].start + edge[i].num)) * sizeof(peye->paPoint[0])
			);
			peye->Point_N -= edge[i].num;
		}
		++n;
	}/**/
	
//	Eye_S4_Edge_Line (peye, M_PI, -M_PI/100.0f, &ae, &r);
	
/*	for (a = ao; a <= ao+M_PI_2+M_PI_4; a += M_PI_2) {
//	{
	//	printf ("mino ");
		u08 skip;
		float x0, y0, x1, y1;
		float t0 = a, t1 = t0 - M_PI;
		
	//	skip = edge (peye, i, t0, &x0, &y0);
		skip = edge (peye, a/(M_PI/100.0f), t0, &x0, &y0);
		
		if (skip == 0) {
			t0 = t1;
		//	skip = edge (peye, 110 + i, t0, &x0, &y0);
			skip = edge (peye, 110 + a/(M_PI/100.0f), t0, &x0, &y0);
		}else {
		//	skip += edge (peye, 110 + i, t1, &x1, &y1);
			skip += edge (peye, 110 + a/(M_PI/100.0f), t1, &x1, &y1);
		}
		++i;
		
		if (skip == 0)
			continue;
		
	}
	}/**/
	
//	printf ("Eye_S5_Edges: got %ld points\n", peye->Point_N);
	for (i = 0; i < peye->Point_N; ++i) {
		Eye_Ellipse2LinDraw_Pix_ad (peye,
			angle_norm_0_2pi(atan2(peye->paPoint[i].y-dcam.P.y,peye->paPoint[i].x-dcam.P.x)),
			ddist(peye->paPoint[i].x,peye->paPoint[i].y, dcam.P.x,dcam.P.y),
			0xFF<<16);
	}/**/
	return;
}


void	Eye_S5		(tEye* peye)
{
	Eye_Ellipse2LinDraw (peye);
	
	if (!peye->Point_Max) {
		peye->Point_Max = 2048;
		peye->paPoint = malloc (peye->Point_Max * sizeof(tV2f));
	}
	
	Eye_S5_Edges (peye);
	
	Eye_Points_Sort (peye);
	Eye_Points_Fit (peye);
	
	return;
}


//haha back to the roots

#define dGlint_FID 16
#define dGlint_Num 4

#define dMarkRawOut(_x,_y)	((_x) < 0 || (_x) >= dcam.FF.Max_R || (_y) < 0 || (_y) >= dcam.FF.Max_R)
#define dMarkOut(_x,_y)	((_x) < dcam.FF.Mark_P.x || (_x) >= dcam.FF.Mark_P.x+dcam.FF.Max_R || (_y) < dcam.FF.Mark_P.y || (_y) >= dcam.FF.Mark_P.y+dcam.FF.Max_R)

#define dMarkRaw(_x,_y)	(*(dcam.FF.paMark + (si)(_x) + (si)(_y)*(dcam.FF.Max_R)))
#define dMark(_x,_y)	dMarkRaw((si)(_x)-dcam.FF.Mark_P.x, (si)(_y)-dcam.FF.Mark_P.y)


void	Eye_FF_Mark_Crap		(tEye* peye, tCam* pcam, si x, si y)
{
	//if (dpixout(x,y))
		//return;
	if (dMarkOut(x,y))
		return;
	if (dMark(x,y))
		return;
	
	//printf("dopix(x,y)->Y %d\n", dopix(x,y)->Y);
	if (dopix(x,y)->Y < dcam.FF.Y) {
		dMark(x,y) = dcam.FF.tmpID;
		
		dset_c0(x,y);
		//dnpix(x,y)->U = 0;
		//dnpix(x,y)->V = 0;
		
		++dcam.FF.tmpNum;
		dcam.FF.tmpP.x += x;
		dcam.FF.tmpP.y += y;
		//Eye_CirView_Point_xy (peye, pcam, x-dcam.FF.Mark_P.x, y-dcam.FF.Mark_P.y, gColARGB);
		Eye_FF_Mark_Crap (peye, pcam, x-1, y);
		Eye_FF_Mark_Crap (peye, pcam, x+1, y);
		Eye_FF_Mark_Crap (peye, pcam, x, y-1);
		Eye_FF_Mark_Crap (peye, pcam, x, y+1);
	}
}

si	Eye_FF_G_Mark_Crap	(tEye* peye, tCam* pcam, si x, si y)
{
	//if (dpixout(x,y))
		//return;
	if (dMarkOut(x,y))
		return;
	if (dMark(x,y))
		return;
	
	if (dopix(x,y)->Y > dcam.FF.GY) {
		dMark(x,y) = dcam.FF.tmpID;
		
		//dnpix(x,y)->U = 0xF;
		//dnpix(x,y)->V = 0xF;
		
		++dcam.FF.tmpNum;
		dcam.FF.tmpP.x += x;
		dcam.FF.tmpP.y += y;
		//Eye_CirView_Point_xy (peye, pcam, x-dcam.FF.Mark_P.x, y-dcam.FF.Mark_P.y, gColARGB);
		Eye_FF_G_Mark_Crap (peye, pcam, x-1, y);
		Eye_FF_G_Mark_Crap (peye, pcam, x+1, y);
		Eye_FF_G_Mark_Crap (peye, pcam, x, y-1);
		Eye_FF_G_Mark_Crap (peye, pcam, x, y+1);
	}
}

si	Eye_FF_Mark2Pos		(tEye* peye, tCam* pcam, u08 id, tV2f *pret)
{
	si ax = 0, ay = 0;
	si num = 0;
	si x, y;
	for (y = 0; y < dcam.FF.Max_R; ++y) {
		for (x = 0; x < dcam.FF.Max_R; ++x) {
			if (	dMarkRaw(x,y) == id) {
				ax += x;
				ay += y;
				//Eye_CirView_Point_xy (peye, pcam, x, y, gColARGB);
				++num;
			}
		}
	}/**/
	pret->x = (float)ax / num;
	pret->y = (float)ay / num;
	return num;
}


si	Eye_FF_EdgeGet		(tEye* peye, tCam* pcam, u08 id, si __i, float t, float* px, float* py)
{
//	si border = peye->Exp_R*peye->Exp_R*3;
	u08 wrote = 0;
	float x = dcam.P.x, y = dcam.P.y;
	si n = 0;
//	printf ("Eye_FF_EdgeGet  id %d  t %f\n", id, t);
	for (; !dMarkOut(x,y); ) {
		if (	dMark(x,y) != id
		) {
		//	return 1;
			float x1 = x, y1 = y;
			for (; !dMarkOut(x1,y1); ) {
				if (	dMark(x1,y1) == id
				) {
					break;
				}
				x1 += cos(t);
				y1 += sin(t);
			}/**/
			if (ddist2(x1,y1,x,y) >= dpow2(16)) {
				if (!wrote) {
					wrote = 1;
					*px = x+0.000000000001f;
					*py = y+0.000000000001f;
					
					//Eye_CirView_Point_xy (peye, pcam, x-dcam.FF.Mark_P.x, y-dcam.FF.Mark_P.y, 0xFFFFFF);
					return 1;
				}
			}
			x = x1;
			y = y1;
			continue;
		}/**/
		
		x += cos(t);
		y += sin(t);
		++n;
	}/**/
	return wrote;
	return 0;
}

void	Eye_FF_Fit			(tEye* peye, tCam* pcam, u08 id)
{
//	Eye_SFit (peye);	return;
//	si border = 4*peye->Exp_R;
	
//	float ax = 0, ay = 0;
//	si n = 0;
	
	si i = 0;
	
	si	(*edge)	(tEye* peye, tCam* pcam, u08 id, si i, float t, float* px, float* py);
	edge = Eye_FF_EdgeGet;
	
	float inc = peye->AngRes*deg2rad;
	
//	if (peye == &gM.Left)
//		printf ("Eye_FF_Fit  %d  %f %f    inc %f\n", id, dcam.P.x, dcam.P.y,     inc);
	
	float ao, a;
	for (a = 0; a < M_PI; a += inc) {
//	for (a = -M_PI_2; a < M_PI_2; a += inc) {
//	for (a = -M_PI_4; a < M_PI_4; a += inc) {
//	for (ao = 0; ao < M_PI_2; ao += inc) {
//	for (a = ao; a <= ao+M_PI_2+M_PI_4; a += M_PI_2) {
//	{
	//	printf ("mino ");
		u08 skip;
		float x0, y0, x1, y1;
		float t0 = a, t1 = t0 - M_PI;
	//	printf ("Eye_FF_Fit  %f %f\n", t0, t1);
		
	//	skip = edge (peye, pcam, id, i, t0, &x0, &y0);
		skip = edge (peye, pcam, id, a/(inc), t0, &x0, &y0);
		
		if (skip == 0) {
			t0 = t1;
		//	skip = edge (peye, pcam, id, 110 + i, t0, &x0, &y0);
			skip = edge (peye, pcam, id, 110 + a/(inc), t0, &x0, &y0);
		}else {
		//	skip += edge (peye, pcam, id, 110 + i, t1, &x1, &y1);
			skip += edge (peye, pcam, id, 110 + a/(inc), t1, &x1, &y1);
		}
		++i;
		
		if (skip == 0)
			continue;
		
		if (skip == 1) {
		//	dnpix(x0,y0)->Y = 0xFF;
		//	dnpix(x0,y0)->U = 0xF;
		//	dnpix(x0,y0)->V = 0xF;
			
			float t, dx, dy, xx, yy;
			dx = x0 - dcam.P.x;
			dy = y0 - dcam.P.y;
			t = atan2(dy,dx)/* - peye->Aa*/;
			
			xx = (peye->Ax * cos(t) * cos(peye->Aa) - peye->Ay * sin(t) * sin(peye->Aa));
			yy = (peye->Ax * cos(t) * sin(peye->Aa) + peye->Ay * sin(t) * cos(peye->Aa));
			
			float diff = sqrt(dx*dx+dy*dy) - sqrt(xx*xx+yy*yy);
			if (fabs(diff) >= 0.1f) {
			//	printf ("f1\n");
				peye->Ax += peye->Fit_Scale*diff*fabs(cos(t));
				peye->Ay += peye->Fit_Scale*diff*fabs(sin(t));
				
				EyeC_Xset (peye, pcam, dcam.P.x + peye->Fit_Trans*(diff )*cos(t));
				EyeC_Yset (peye, pcam, dcam.P.y + peye->Fit_Trans*(diff )*sin(t));
			}
			continue;
		}
	/*	dnpix(x0,y0)->Y = 0xFF;
		dnpix(x0,y0)->U = 0x0;
		dnpix(x0,y0)->V = 0x0;
		
		dnpix(x1,y1)->Y = 0xFF;
		dnpix(x1,y1)->U = 0x0;
		dnpix(x1,y1)->V = 0x0;/**/
		
	//	Eye_CirView_Point_xy (peye, pcam, dcam.FF.Mark_P.x - x0, dcam.FF.Mark_P.y - y0, 0x00FF00);
	//	Eye_CirView_Point_xy (peye, pcam, dcam.FF.Mark_P.x - x1, dcam.FF.Mark_P.y - y1, 0x00FF00);
		{
			float ox0 = (peye->Ax * cos(t0) * cos(peye->Aa) - peye->Ay * sin(t0) * sin(peye->Aa));
			float ox1 = (peye->Ax * cos(t1) * cos(peye->Aa) - peye->Ay * sin(t1) * sin(peye->Aa));
			
			float oy0 = (peye->Ax * cos(t0) * sin(peye->Aa) + peye->Ay * sin(t0) * cos(peye->Aa));
			float oy1 = (peye->Ax * cos(t1) * sin(peye->Aa) + peye->Ay * sin(t1) * cos(peye->Aa));
			
		//	printf ("Heee? pic %f %f   %f %f\n", x0 - dcam.P.x, y0 - dcam.P.y, x1 - dcam.P.x, y1 - dcam.P.y);
		//	printf ("Heee? eli %f %f   %f %f\n", ox0, oy0, ox1, oy1);
			
			float sdx = fabsf(x1 - x0) - fabsf(ox1 - ox0);
			float sdy = fabsf(y1 - y0) - fabsf(oy1 - oy0);
			
		//	ax += fabsf((x1 - x0)*cos(t0));
		//	ay += fabsf((y1 - y0)*sin(t0));
		//	++n;
			
			float pdx = (x0+x1)*0.5f - (2.0f*dcam.P.x + ox0+ox1)*0.5f;
			float pdy = (y0+y1)*0.5f - (2.0f*dcam.P.y + oy0+oy1)*0.5f;
			
			EyeC_Xset (peye, pcam, dcam.P.x + 0.10f * pdx );// *cos(t));
			EyeC_Yset (peye, pcam, dcam.P.y + 0.10f * pdy );// *sin(t));
			
			peye->Ax += 0.10f*sdx*fabs(cos(t0));
			peye->Ay += 0.10f*sdy*fabs(sin(t0));
		}
	}
//	}
	return;
}

float	Eye_FF_FitConf		(tEye* peye, tCam* pcam, u08 id)
{
	float mse = 0;
	
	si i = 0;
	si	(*edge)	(tEye* peye, tCam* pcam, u08 id, si i, float t, float* px, float* py);
	edge = Eye_FF_EdgeGet;
	
	float inc = peye->AngRes*deg2rad;
	
	float ao, a;
//	for (a = 0; a < M_PI; a += inc) {
//	for (a = -M_PI_2; a < M_PI_2; a += inc) {
//	for (a = -M_PI_4; a < M_PI_4; a += inc) {
	for (ao = 0; ao < M_PI_2; ao += inc) {
	for (a = ao; a <= ao+M_PI_2+M_PI_4; a += M_PI_2) {
//	{
		u08 skip;
		float x0, y0, x1, y1;
		float t0 = a, t1 = t0 - M_PI;
		
	//	skip = edge (peye, pcam, id, i, t0, &x0, &y0);
		skip = edge (peye, pcam, id, a/(inc), t0, &x0, &y0);
		
		if (skip == 0) {
			t0 = t1;
		//	skip = edge (peye, pcam, id, 110 + i, t0, &x0, &y0);
			skip = edge (peye, pcam, id, 110 + a/(inc), t0, &x0, &y0);
		}else {
		//	skip += edge (peye, pcam, id, 110 + i, t1, &x1, &y1);
			skip += edge (peye, pcam, id, 110 + a/(inc), t1, &x1, &y1);
		}
		++i;
		
		if (skip == 0)
			continue;
		
		if (skip == 1) {
			dnpix(x0,y0)->Y = 0xFF;
			dnpix(x0,y0)->U = 0xF;
			dnpix(x0,y0)->V = 0xF;
			
			float t, dx, dy, xx, yy;
			dx = x0 - dcam.P.x;
			dy = y0 - dcam.P.y;
			t = atan2(dy,dx)/* - peye->Aa*/;
			
			xx = (peye->Ax * cos(t) * cos(peye->Aa) - peye->Ay * sin(t) * sin(peye->Aa));
			yy = (peye->Ax * cos(t) * sin(peye->Aa) + peye->Ay * sin(t) * cos(peye->Aa));
			
			float diff = sqrt(dx*dx+dy*dy) - sqrt(xx*xx+yy*yy);
			mse += dpow2(diff);
			continue;
		}
		{
			float ox0 = (peye->Ax * cos(t0) * cos(peye->Aa) - peye->Ay * sin(t0) * sin(peye->Aa));
			float ox1 = (peye->Ax * cos(t1) * cos(peye->Aa) - peye->Ay * sin(t1) * sin(peye->Aa));
			
			float oy0 = (peye->Ax * cos(t0) * sin(peye->Aa) + peye->Ay * sin(t0) * cos(peye->Aa));
			float oy1 = (peye->Ax * cos(t1) * sin(peye->Aa) + peye->Ay * sin(t1) * cos(peye->Aa));
			
		//	printf ("Heee? pic %f %f   %f %f\n", x0 - dcam.P.x, y0 - dcam.P.y, x1 - dcam.P.x, y1 - dcam.P.y);
		//	printf ("Heee? eli %f %f   %f %f\n", ox0, oy0, ox1, oy1);
			
			float sdx = fabsf(x1 - x0) - fabsf(ox1 - ox0);
			float sdy = fabsf(y1 - y0) - fabsf(oy1 - oy0);
			
			mse += dpow2(sdx*fabs(cos(t0)));
			mse += dpow2(sdy*fabs(sin(t0)));
		//	ax += fabsf((x1 - x0)*cos(t0));
		//	ay += fabsf((y1 - y0)*sin(t0));
		//	++n;
			
			float pdx = (x0+x1)*0.5f - (2.0f*dcam.P.x + ox0+ox1)*0.5f;
			float pdy = (y0+y1)*0.5f - (2.0f*dcam.P.y + oy0+oy1)*0.5f;
			
			mse += dpow2(pdx);
			mse += dpow2(pdy);
			
		}
	}
	}
	return mse;
}

float	Eye_FF_Mark2AvgDist	(tEye* peye, tCam* pcam, float ox, float oy, u08 id)
{
	float adist = 0;
	si num = 0;
	si x, y;
	for (y = 0; y < dcam.FF.Max_R; ++y) {
		for (x = 0; x < dcam.FF.Max_R; ++x) {
			if (dMarkRaw(x,y) == id) {
				adist += ddist(x,y,ox,oy);
				++num;
			}
		}
	}/**/
	return adist / num;
}


float	Eye_FF_Mark2Confidence	(tEye* peye, tCam* pcam, float ox, float oy, u08 id, float perfr)
{
	float adist = 0;
	si num = 0;
	si x, y;
	perfr = dpow2(perfr);
	for (y = 0; y < dcam.FF.Max_R; ++y) {
		for (x = 0; x < dcam.FF.Max_R; ++x) {
			if (dMarkRaw(x,y) == id) {
				if (ddist2(x,y,ox,oy) <= perfr) {
					++adist;
				}
				++num;
			}else {
				if (ddist2(x,y,ox,oy) <= perfr) {
					--adist;
				}
			}
		}
	}/**/
	return adist / num;
}

float	Eye_FF_Mark2Conf		(tEye* peye, tCam* pcam, u08 id, float ox, float oy)
{
	float mse = 0;
	si num = 0;
	si x, y;
	
	for (y = oy-peye->Ay; y <  oy+peye->Ay; ++y) {
		if (dMarkRawOut(0,y))
			continue;
		for (x = ox-peye->Ax; x < ox+peye->Ax; ++x) {
			if (peye == &gM.Left)
			if (dMarkRawOut(x,0))
				continue;
			float diff = ddist2(x,y,ox,oy) - Eye_Ellipse_xydist2(peye,atan2f(y-oy,x-ox), 0,0);
		//	float diff = ddist2(x,y,ox,oy) - dpow2(peye->Ax);
			if (dMarkRaw(x,y) == id) {
				if (diff > 0) {
					mse += dpow2(diff);
					//Eye_CirView_Point_xy (peye, pcam, x, y, 0xFF0000);
				}
			}else {
				if (diff < 0) {
					mse += dpow2(diff);
					//Eye_CirView_Point_xy (peye, pcam, x, y, 0xFF0000);
				}
			}
			++num;
		}
	}/**/
	return mse / num / (peye->Ax+peye->Ay);
}


void	Eye_FF			(tEye* peye, tCam* pcam)
{
	ay = dcam.FF.Y;
	au = 0;
	av = 0;
	
	memset (dcam.FF.paMark, 0, dpow2(dcam.FF.Max_R) * sizeof(dcam.FF.paMark[0]));
	
	dcam.FF.Mark_P.x = (dcam.P.x + 0.5f) - dcam.FF.Max_R/2;
	dcam.FF.Mark_P.y = (dcam.P.y + 0.5f) - dcam.FF.Max_R/2;
	
	gColARGB = 0x00FF00;
	
	dcam.FF.tmpID = 1;
	
/*	if (peye->Fit == eEye_Fit_FF_Diode) {
		dcam.FF.tmpNum = 0;
		dcam.FF.tmpP.x = 0;
		dcam.FF.tmpP.y = 0;
		Eye_FF_G_Mark_Crap (peye, pcam, dcam.P.x, dcam.P.y);
		
		if (dcam.FF.tmpNum >= 50) {
			dcam.FF.tmpP.x /= dcam.FF.tmpNum;
			dcam.FF.tmpP.y /= dcam.FF.tmpNum;
			
			dcam.P.x = dcam.FF.tmpP.x;
			dcam.P.y = dcam.FF.tmpP.y;
		}
		
		peye->Ax = dcam.FF.Perf_R;
		peye->Ay = dcam.FF.Perf_R;
		return;
	}/**/
	si	(*pfmark)	(tEye* peye, tCam* pcam, si x, si y);
	if (peye->Fit == eEye_Fit_FF_Diode)
		pfmark = Eye_FF_G_Mark_Crap;
	else
		pfmark = Eye_FF_Mark_Crap;
	
	if (dcam.FF.Perf_R != 0) {
		u08 best_id = 0;
		float best_ar = NAN;
	//	printf ("Search for main %f\n", dcam.FF.Perf_R);
		float best_x = dcam.P.x, best_y = dcam.P.y;
		
		si ox = dcam.P.x, oy = dcam.P.y, x, y;
		si search_r = dcam.FF.Search_R;
		for (y = oy - search_r; y < oy + search_r; y += 10) {
			for (x = ox - search_r; x < ox + search_r; x += 10) {
				dcam.FF.tmpNum = 0;
				dcam.FF.tmpP.x = 0;
				dcam.FF.tmpP.y = 0;
				pfmark (peye, pcam, x, y);
				
				if (dcam.FF.tmpNum >= 100) {
					dcam.FF.tmpP.x /= dcam.FF.tmpNum;
					dcam.FF.tmpP.y /= dcam.FF.tmpNum;
					#if 0
					float ar = Eye_FF_Mark2AvgDist (peye, pcam, dcam.FF.tmpP.x-dcam.FF.Mark_P.x, dcam.FF.tmpP.y-dcam.FF.Mark_P.y, dcam.FF.tmpID);
					
				//	printf ("AvgDist	%d	%f\n", dcam.FF.tmpID, ar);
					if (isfinite(ar)
						&& fabsf(dcam.FF.Perf_R - ar) <= dcam.FF.MaxDiff_R
						&& (!isfinite(best_ar) 
							|| fabsf(dcam.FF.Perf_R - ar) < fabsf(dcam.FF.Perf_R - best_ar))
					) {
						best_id = dcam.FF.tmpID;
						best_ar = ar;
						best_x = dcam.FF.tmpP.x;
						best_y = dcam.FF.tmpP.y;
					}
					#elif 1
					float ar = Eye_FF_Mark2Confidence (peye, pcam, dcam.FF.tmpP.x-dcam.FF.Mark_P.x, dcam.FF.tmpP.y-dcam.FF.Mark_P.y, dcam.FF.tmpID, dcam.FF.Perf_R);
					
				//	printf ("Confidence	%d	%f\n", dcam.FF.tmpID, ar);
					if (ar >= dcam.FF.MaxDiff_R
						&& (!isfinite(best_ar) 
							|| (ar > best_ar))
					) {
						best_id = dcam.FF.tmpID;
						best_ar = ar;
						best_x = dcam.FF.tmpP.x;
						best_y = dcam.FF.tmpP.y;
					}
					#elif 1
					peye->P = dcam.FF.tmpP;
					Eye_FF_Fit (peye, pcam, dcam.FF.tmpID);
					
					float ar = Eye_FF_FitConf (peye, pcam, dcam.FF.tmpID);
				//	float ar = Eye_FF_Mark2Conf (peye, pcam, dcam.FF.tmpID, dcam.P.x-dcam.FF.Mark_P.x, dcam.P.y-dcam.FF.Mark_P.y);
					
					if (peye == &gM.Left)
						printf ("Confidence	Ax %f	%d	%f\n", peye->Ax, dcam.FF.tmpID, ar);
					
					if (fabsf(peye->Ax - dcam.FF.Perf_R) <= 15
						&& ar >= 0.5
						&& (!isfinite(best_ar) 
							|| (ar < best_ar))
					) {
						best_id = dcam.FF.tmpID;
						best_ar = ar;
						best_x = dcam.FF.tmpP.x;
						best_y = dcam.FF.tmpP.y;
					}
					#else
					peye->P = dcam.FF.tmpP;
					Eye_FF_Fit (peye, pcam, dcam.FF.tmpID);
					
					float ar = Eye_FF_Mark2Confidence (peye, pcam, dcam.FF.tmpP.x-dcam.FF.Mark_P.x, dcam.FF.tmpP.y-dcam.FF.Mark_P.y, dcam.FF.tmpID, (peye->Ax+peye->Ay)/2.0f);
					
				//	printf ("Confidence	%d	%f\n", dcam.FF.tmpID, ar);
					if (ar >= dcam.FF.MaxDiff_R
						&& (!isfinite(best_ar) 
							|| (ar > best_ar))
					) {
						best_id = dcam.FF.tmpID;
						best_ar = ar;
						best_x = dcam.FF.tmpP.x;
						best_y = dcam.FF.tmpP.y;
					}
					#endif
					dcam.FF.tmpID++;
				}
			}
		}
		if (peye == &gM.Left)
			printf ("best id %d ar %f\n", best_id, best_ar);
		dcam.P.x = best_x;
		dcam.P.y = best_y;
		
	/*	if (best_id) {
			dcam.FF.tmpP = dcam.P;
			Eye_FF_Fit (peye, pcam, best_id);
			peye->P = dcam.FF.tmpP;
		}/**/
		
		peye->Ax = dcam.FF.Perf_R;
		peye->Ay = dcam.FF.Perf_R;
	}else {
		dcam.FF.tmpNum = 0;
		dcam.FF.tmpP.x = 0;
		dcam.FF.tmpP.y = 0;
		pfmark (peye, pcam, dcam.P.x, dcam.P.y);
		
		if (dcam.FF.tmpNum > 15) {
			dcam.P.x = dcam.FF.tmpP.x / dcam.FF.tmpNum;
			dcam.P.y = dcam.FF.tmpP.y / dcam.FF.tmpNum;
		}/**/
	/*	tV2f pos;
		if (Eye_FF_Mark2Pos (peye, pcam, 1, &pos) > 20) {
			dcam.P.x = pos.x + dcam.FF.Mark_P.x;
			dcam.P.y = pos.y + dcam.FF.Mark_P.y;
		}/**/
	}
	
	if (1) {
		if (peye == &gM.Left)
		/*	printf ("Glint get Left\n")/**/;
		else if (peye == &gM.Right)
		/*	printf ("Glint get Right\n")/**/;
		else
			return;
		
		struct {
			si num;
			tV2f pos;
			u32 col;
		}g[dGlint_Num];
		
		g[0].col = 0xFFFF00;
		g[1].col = 0xFF00FF;
		g[2].col = 0x00FFFF;
		g[3].col = 0xFF0000;
		
		si glint_idx = 0;
		si search_r = dcam.FF.GSearch_R;
		gColARGB = 0xFF;
		si x, y;
		for (y = dcam.P.y; y < dcam.P.y + search_r; ++y) {
			for (x = dcam.P.x - search_r; x < dcam.P.x + search_r; ++x) {
				gColARGB = g[glint_idx].col;
				dcam.FF.tmpNum = 0;
				dcam.FF.tmpP.x = 0;
				dcam.FF.tmpP.y = 0;
				dcam.FF.tmpID = dGlint_FID+glint_idx;
				Eye_FF_G_Mark_Crap (peye, pcam, x, y);
				g[glint_idx].num = dcam.FF.tmpNum;
				if (g[glint_idx].num) {
					g[glint_idx].pos.x = dcam.FF.tmpP.x / dcam.FF.tmpNum;
					g[glint_idx].pos.y = dcam.FF.tmpP.y / dcam.FF.tmpNum;
					
				//	printf ("Got Glint %d num pixs %d,   xy %f %f\n", glint_idx, g[glint_idx].num, g[glint_idx].pos.x, g[glint_idx].pos.y);
					
					++glint_idx;
				}/**/
				if (glint_idx >= dGlint_Num) {
					goto Got_Glints;
				}
			}
		}/**/
Got_Glints:	if (gM.Eye_GlintMode == 2) {
			si min0 = -1, min1 = -1;
			float mind = 100;
			si i, j;
			for (i = 0; i < glint_idx; ++i) {
				for (j = i+1; j < glint_idx; ++j) {
					if (fabsf(g[i].pos.y - g[j].pos.y) <= mind) {
						mind = fabsf(g[i].pos.y - g[j].pos.y);
						min0 = i;
						min1 = j;
					}
				}
			}
			if (min0 != -1) {
				if (g[min0].pos.x <= g[min1].pos.x) {
					peye->G0 = g[min0].pos;
					peye->G1 = g[min1].pos;
				}else  {
					peye->G0 = g[min1].pos;
					peye->G1 = g[min0].pos;
				}
				Eye_GV_Calc (peye);
			//	printf ("Eye GV %f %f\n", peye->GV.x, peye->GV.y);
				Eye_GP_Calc (peye);
			}
		}else if (gM.Eye_GlintMode == 1) {
			si idx = -1;
			si max_num = 0;
			si i, j;
			for (i = 0; i < glint_idx; ++i) {
				if (g[i].num >= max_num) {
					max_num = g[i].num;
					idx = i;
				}
			}
			if (idx != -1) {
				peye->G0 = g[idx].pos;
				Eye_GV_Calc (peye);
			//	printf ("Eye GV %f %f\n", peye->GV.x, peye->GV.y);
			}
		}/**/
		
	}
}

#undef dGlint_FID
#undef dGlint_Num
#undef dMarkOut
#undef dMarkRaw
#undef dMark


float	Eye_CFit_AvgDiff	(tEye* peye)
{
	si border = 4*peye->Exp_R;
	float avgdiff = 0, diff_num = 0;
	si x, y;
	for (y = dcam.P.y - border; y < dcam.P.y + border; ++y) {
		for (x = dcam.P.x - border; x < dcam.P.x + border; ++x) {
			if (	1//dnpix(x,y)->Y == 0xFF
				&& dnpix(x,y)->U == 0xF
				&& dnpix(x,y)->V == 0xF
			) {
				float t, dx, dy, xx, yy;
				dx = x - dcam.P.x;
				dy = y - dcam.P.y;
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
	for (y = dcam.P.y - border; y < dcam.P.y + border; ++y) {
		for (x = dcam.P.x - border; x < dcam.P.x + border; ++x) {
			if (	1//dnpix(x,y)->Y == 0xFF
				&& dnpix(x,y)->U == 0xF
				&& dnpix(x,y)->V == 0xF
			) {
				float t, dx, dy, xx, yy;
				dx = x - dcam.P.x;
				dy = y - dcam.P.y;
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
/*	Eye_S3Fit (peye);	return;
	
	Eye_SFit (peye);
	printf ("%ld  ", Eye_CFit_Remove (peye, Eye_CFit_AvgDiff (peye)*4.0f));
	Eye_SFit (peye);
	printf ("%ld  ", Eye_CFit_Remove (peye, Eye_CFit_AvgDiff (peye)*3.3f));
	Eye_SFit (peye);
	printf ("%ld  ", Eye_CFit_Remove (peye, Eye_CFit_AvgDiff (peye)*2.2f));
	Eye_SFit (peye);
	printf ("%ld  ", Eye_CFit_Remove (peye, Eye_CFit_AvgDiff (peye)*1.2f));
	Eye_SFit (peye);
	
	printf ("\n");/**/
	
	Eye_CalcAYUV (peye, pcam, 15);
	Eye_S5 (peye);
	Eye_CalcAYUV (peye, pcam, 15);
	Eye_S3Fit (peye);
	
}

void	Eye_Fit		(tEye* peye, tCam* pcam)
{
	switch (peye->Fit) {
	case eEye_Fit_S0:
		Eye_CalcAYUV (peye, pcam, 12);
		Eye_S0 (peye);
		break;
	case eEye_Fit_S2Fit:
		Eye_OldFF (peye);
		Eye_EdgeMark (peye);
		Eye_S2Fit (peye);
		break;
	case eEye_Fit_S3Fit_START ... eEye_Fit_S3Fit_END:
	//	Eye_CalcAYUV (peye, pcam, 12);
		ay = peye->Pix_Bright;
		Eye_S3Fit (peye);
		break;
	case eEye_Fit_S4Fit_START ... eEye_Fit_S4Fit_END:
		Eye_CalcAYUV (peye, pcam, 12);
		Eye_S4 (peye);
		break;
	case eEye_Fit_S5_START ... eEye_Fit_S5_END:
		Eye_CalcAYUV (peye, pcam, 12);
		Eye_S5 (peye);
		break;
	case eEye_Fit_FF_START ... eEye_Fit_FF_END:
		Eye_FF (peye, pcam);
		break;
	case eEye_Fit_C:
		Eye_CFit (peye);
		break;
	case eEye_Fit_SFit:
	default:
		abort();
		Eye_OldFF (peye);
		Eye_EdgeMark (peye);
		Eye_SFit (peye);
		break;
	}
}

void	Eye_Clip		(tEye* peye)
{
	si x, y;
	#define dmarg 10
	if (dcam.P.x < dmarg)
		dcam.P.x = dmarg;
	if (dcam.P.x > pcam->Image_W-dmarg)
		dcam.P.x = pcam->Image_W-dmarg;
	if (dcam.P.y < dmarg)
		dcam.P.y = dmarg;
	if (dcam.P.y > pcam->Image_H-dmarg)
		dcam.P.y = pcam->Image_H-dmarg;
	#undef dmarg
	if (peye->Ax > 100)
		peye->Ax = 100;
	if (peye->Ay > 100)
		peye->Ay = 100;
	if (peye->Ax < 5)
		peye->Ax = 5;
	if (peye->Ay < 5)
		peye->Ay = 5;
	
//	if (isnan(peye->Ax))
//		peye->Ax = peye->Exp_R*2;
//	if (isnan(peye->Ay))
//		peye->Ay = peye->Exp_R*2;
	
}



void	Eye_Draw_Ellipse	(tEye* peye, float tb, float te)
{
	si x, y;
	float t;
	for (t = tb; t < te; t += M_PI/180) {
	//	peye->Ax*x*x + peye->Ay*y*y + peye->Axy*x*y = peye->Ar
		x = (peye->Ax * cos(t) * cos(peye->Aa) - peye->Ay * sin(t) * sin(peye->Aa));
		y = (peye->Ax * cos(t) * sin(peye->Aa) + peye->Ay * sin(t) * cos(peye->Aa));
	//	printf ("x %4d  y %4d\n", x, y);
	//	y = dcam.P.y;
		if (dcam.P.x+x < 0 || dcam.P.y+y < 0 || dcam.P.x+x >= pcam->Image_W || dcam.P.y+y >= pcam->Image_H)
			continue;
		dnpix(dcam.P.x+x,dcam.P.y+y)->Y = 0xFF;
		dnpix(dcam.P.x+x,dcam.P.y+y)->U = 0xF;
		dnpix(dcam.P.x+x,dcam.P.y+y)->V = 0xF;
	}/**/
}
void	Eye_Draw		(tEye* peye, tCam* pcam)
{
	si x, y;
	dset_c1(dcam.P.x,dcam.P.y);
//	Eye_Draw_Ellipse (peye, 0, M_PI*2);	return;
/*	Eye_Draw_Ellipse (peye, 0 - 10*M_PI/180,		0 + 10*M_PI/180);
	Eye_Draw_Ellipse (peye, M_PI_2 - 10*M_PI/180,	M_PI_2 + 10*M_PI/180);
	Eye_Draw_Ellipse (peye, M_PI - 10*M_PI/180,	M_PI + 10*M_PI/180);
	Eye_Draw_Ellipse (peye, -M_PI_2 - 10*M_PI/180,	-M_PI_2 + 10*M_PI/180);
	return;
	/**/
	dset_c1(dcam.P.x - peye->Ax,		dcam.P.y - peye->Ay);
	dset_c1(dcam.P.x - peye->Ax+1,	dcam.P.y - peye->Ay);
	dset_c1(dcam.P.x - peye->Ax+2,	dcam.P.y - peye->Ay);
	dset_c1(dcam.P.x - peye->Ax,		dcam.P.y - peye->Ay+1);
	dset_c1(dcam.P.x - peye->Ax,		dcam.P.y - peye->Ay+2);
	
	dset_c1(dcam.P.x + peye->Ax,		dcam.P.y - peye->Ay);
	dset_c1(dcam.P.x + peye->Ax-1,	dcam.P.y - peye->Ay);
	dset_c1(dcam.P.x + peye->Ax-2,	dcam.P.y - peye->Ay);
	dset_c1(dcam.P.x + peye->Ax,		dcam.P.y - peye->Ay+1);
	dset_c1(dcam.P.x + peye->Ax,		dcam.P.y - peye->Ay+2);
	
	dset_c1(dcam.P.x - peye->Ax,		dcam.P.y + peye->Ay);
	dset_c1(dcam.P.x - peye->Ax+1,	dcam.P.y + peye->Ay);
	dset_c1(dcam.P.x - peye->Ax+2,	dcam.P.y + peye->Ay);
	dset_c1(dcam.P.x - peye->Ax,		dcam.P.y + peye->Ay-1);
	dset_c1(dcam.P.x - peye->Ax,		dcam.P.y + peye->Ay-2);
	
	dset_c1(dcam.P.x + peye->Ax,		dcam.P.y + peye->Ay);
	dset_c1(dcam.P.x + peye->Ax-1,	dcam.P.y + peye->Ay);
	dset_c1(dcam.P.x + peye->Ax-2,	dcam.P.y + peye->Ay);
	dset_c1(dcam.P.x + peye->Ax,		dcam.P.y + peye->Ay-1);
	dset_c1(dcam.P.x + peye->Ax,		dcam.P.y + peye->Ay-2);
	/**/
	
}



tV2f	Eye_map_point	(tEye* peye, tV2f p)
{
	return map_point (&peye->Homo, p);
}


