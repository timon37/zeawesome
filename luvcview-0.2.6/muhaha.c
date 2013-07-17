
#include "muhaha.h"


void	S_Draw_Line_2d	(SDL_Surface* surf, si x0, si y0, si x1, si y1);

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


tRGBA gCol;
//tPix gCol;
u32 gColARGB;

tCam* pcam;

tM gM =
{
	.Y_Level = 128,
	.DeSat = 0,
	.Head = {
		.DotC = {
			.Ax = 20,
			.Ay = 20,
			.Aa = 0,
		//	.Aa = M_PI/4,
			
			.Exp_R = 15,
			
			.Fit = eEye_Fit_S3Fit_Point,
			.Fit_Scale = 0.2f,
			.Fit_Trans = 0.2f,
			.S2Fit_Scale = 0.8f,
			.S2Fit_Trans = 0.7f,
			
			.LinView = {950, 0},
		},
		.DotL = {
			.Ax = 20,
			.Ay = 20,
			.Aa = 0,
		//	.Aa = M_PI/4,
			
			.Exp_R = 15,
			
			.Fit = eEye_Fit_S3Fit_Point,
			.Fit_Scale = 0.2f,
			.Fit_Trans = 0.2f,
			.S2Fit_Scale = 0.8f,
			.S2Fit_Trans = 0.7f,
			
			.LinView = {1050, 0},
		},
		.DotR = {
			.Ax = 20,
			.Ay = 20,
			.Aa = 0,
		//	.Aa = M_PI/4,
			
			.Exp_R = 15,
			
			.Fit = eEye_Fit_S3Fit_Point,
			.Fit_Scale = 0.2f,
			.Fit_Trans = 0.2f,
			.S2Fit_Scale = 0.8f,
			.S2Fit_Trans = 0.7f,
			
			.LinView = {1150, 0},
		},
	},
	.Left = {
		.Ax = 30,
		.Ay = 30,
		.Aa = 0,
	//	.Aa = M_PI/4,
		
		.Exp_R = 20,
		.Exp_R = 16,
		
		.Fit = eEye_Fit_S3Fit_Eye,
	//	.Fit = eEye_Fit_S4Fit_Edge0,
	//	.Fit = eEye_Fit_S4Fit_Edge1,
		.Fit = eEye_Fit_S4Fit_Edge2,
		.Fit_Scale = 0.1f,
		.Fit_Trans = 0.1f,
		.S2Fit_Scale = 0.8f,
		.S2Fit_Trans = 0.7f,
		
		.InHead.P.x = 0,
		.InHead.P.y = 0,
		
		.LinView = {640, 0},
	},
	.Right = {
		.Ax = 30,
		.Ay = 30,
		.Aa = 0,
	//	.Aa = M_PI/4,
		
		.Exp_R = 20,
	//	.Exp_R = 16,
		
		.Fit = eEye_Fit_S3Fit_Eye,
	//	.Fit = eEye_Fit_S4Fit_Edge2,
		
		.Fit_Scale = 0.1f,
		.Fit_Trans = 0.1f,
		.S2Fit_Scale = 0.8f,
		.S2Fit_Trans = 0.7f,
		
		.InHead.P.x = 0,
		.InHead.P.y = 0,
		
		.LinView = {640, 240},
	},
};


#include "priv_macros.h"


void	V4f_Print	(tV4f* p)
{
	printf ("%f %f %f %f\n", p->x, p->y, p->z, p->w);
}


void	M3f_Print	(tM3f* p)
{
	si y;
	for (y = 0; y < 3; ++y) {
		printf ("%f\t%f\t%f\n", p->f[y][0], p->f[y][1], p->f[y][2]);
	}
}

void	M4f_Print	(tM4f* p)
{
	si y;
	for (y = 0; y < 4; ++y) {
		printf ("%f\t%f\t%f\t%f\n", p->f[y][0], p->f[y][1], p->f[y][2], p->f[y][3]);
	}
}



void	V2f_add_V2f		(tV2f* pv0, tV2f* pv1)
{
	float x, y;
	x = pv0->x + pv1->x;
	y = pv0->y + pv1->y;
	pv0->x = x;
	pv0->y = y;
}
void	V2f_sub_V2f		(tV2f* pv0, tV2f* pv1)
{
	float x, y;
	x = pv0->x - pv1->x;
	y = pv0->y - pv1->y;
	pv0->x = x;
	pv0->y = y;
}

void	V2f_mul_V2f		(tV2f* pv0, tV2f* pv1)
{
	float x, y;
	x = pv0->x * pv1->x;
	y = pv0->y * pv1->y;
	pv0->x = x;
	pv0->y = y;
}
void	V2f_div_V2f		(tV2f* pv0, tV2f* pv1)
{
	float x, y;
	x = pv0->x / pv1->x;
	y = pv0->y / pv1->y;
	pv0->x = x;
	pv0->y = y;
}

void	V2f_mul_S		(tV2f* pv0, float s)
{
	pv0->x *= s;
	pv0->y *= s;
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



tV4f*	V4f_lneg		(tV4f* pv)
{
	pv->x = -pv->x;
	pv->y = -pv->y;
	pv->z = -pv->z;
	return pv;
}
tV4f	V4f_rneg		(tV4f* pv)
{
	tV4f ret = *pv;
	ret.x = -ret.x;
	ret.y = -ret.y;
	ret.z = -ret.z;
	return ret;
}

void	V4f_add_V4f		(tV4f* pv0, tV4f* pv1)
{
	float x, y, z;
	pv0->x += pv1->x;
	pv0->y += pv1->y;
	pv0->z += pv1->z;
}

void	V4f_sub_V4f		(tV4f* pv0, tV4f* pv1)
{
	float x, y, z;
	pv0->x -= pv1->x;
	pv0->y -= pv1->y;
	pv0->z -= pv1->z;
}


/*void	V4f_mul_M4f		(tV4f* pv0, tM4f* pm1)
{
	tV4f m = {0,0,0,0};
	
	m.x = pv0->x*pm1->x00	+ pv0->y*pm1->x10	+ pv0->z*pm1->x20	+ pv0->w*pm1->x30;
	m.y = pv0->x*pm1->x01	+ pv0->y*pm1->x11	+ pv0->z*pm1->x21	+ pv0->w*pm1->x31;
	m.z = pv0->x*pm1->x02	+ pv0->y*pm1->x12	+ pv0->z*pm1->x22	+ pv0->w*pm1->x32;
	m.w = pv0->x*pm1->x03	+ pv0->y*pm1->x13	+ pv0->z*pm1->x23	+ pv0->w*pm1->x33;
	*pv0 = m;
}/**/

void	M4f_mul_V4f		(tM4f* pm, tV4f* pv)
{
	tV4f m;
/*	m.x = pv->x*pm->x00	+ pv->y*pm->x10	+ pv->z*pm->x20	+ pv->w*pm->x30;
	m.y = pv->x*pm->x01	+ pv->y*pm->x11	+ pv->z*pm->x21	+ pv->w*pm->x31;
	m.z = pv->x*pm->x02	+ pv->y*pm->x12	+ pv->z*pm->x22	+ pv->w*pm->x32;
	m.w = pv->x*pm->x03	+ pv->y*pm->x13	+ pv->z*pm->x23	+ pv->w*pm->x33;/**/
	m.x = pv->x*pm->x00	+ pv->y*pm->x01	+ pv->z*pm->x02	+ pv->w*pm->x03;
	m.y = pv->x*pm->x10	+ pv->y*pm->x11	+ pv->z*pm->x12	+ pv->w*pm->x13;
	m.z = pv->x*pm->x20	+ pv->y*pm->x21	+ pv->z*pm->x22	+ pv->w*pm->x23;
	m.w = pv->x*pm->x30	+ pv->y*pm->x31	+ pv->z*pm->x32	+ pv->w*pm->x33;/**/
	*pv = m;
}


void	M3f_Iden		(tM3f* pm0)
{
	pm0->x00 = 1;	pm0->x01 = 0;	pm0->x02 = 0;
	pm0->x10 = 0;	pm0->x11 = 1;	pm0->x12 = 0;
	pm0->x20 = 0;	pm0->x21 = 0;	pm0->x22 = 1;
}

void	M4f_Iden		(tM4f* pm0)
{
	pm0->x00 = 1;	pm0->x01 = 0;	pm0->x02 = 0;	pm0->x03 = 0;
	pm0->x10 = 0;	pm0->x11 = 1;	pm0->x12 = 0;	pm0->x13 = 0;
	pm0->x20 = 0;	pm0->x21 = 0;	pm0->x22 = 1;	pm0->x23 = 0;
	pm0->x30 = 0;	pm0->x31 = 0;	pm0->x32 = 0;	pm0->x33 = 1;
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

void	M4f_add_M4f		(tM4f* pm0, tM4f* pm1)
{
	int x,y;
	for (y = 0; y < 4; ++y)
		for (x = 0; x < 4; ++x)
			pm0->f[y][x] += pm1->f[y][x];
}

void	M4f_mul_M4f		(tM4f* pm0, tM4f* pm1)
{
	tM4f m;
	m.x00 = pm0->x00*pm1->x00	+ pm0->x01*pm1->x10	+ pm0->x02*pm1->x20	+ pm0->x03*pm1->x30;
	m.x01 = pm0->x00*pm1->x01	+ pm0->x01*pm1->x11	+ pm0->x02*pm1->x21	+ pm0->x03*pm1->x31;
	m.x02 = pm0->x00*pm1->x02	+ pm0->x01*pm1->x12	+ pm0->x02*pm1->x22	+ pm0->x03*pm1->x32;
	m.x03 = pm0->x00*pm1->x03	+ pm0->x01*pm1->x13	+ pm0->x02*pm1->x23	+ pm0->x03*pm1->x33;
	
	m.x10 = pm0->x10*pm1->x00	+ pm0->x11*pm1->x10	+ pm0->x12*pm1->x20	+ pm0->x13*pm1->x30;
	m.x11 = pm0->x10*pm1->x01	+ pm0->x11*pm1->x11	+ pm0->x12*pm1->x21	+ pm0->x13*pm1->x31;
	m.x12 = pm0->x10*pm1->x02	+ pm0->x11*pm1->x12	+ pm0->x12*pm1->x22	+ pm0->x13*pm1->x32;
	m.x13 = pm0->x10*pm1->x03	+ pm0->x11*pm1->x13	+ pm0->x12*pm1->x23	+ pm0->x13*pm1->x33;
	
	m.x20 = pm0->x20*pm1->x00	+ pm0->x21*pm1->x10	+ pm0->x22*pm1->x20	+ pm0->x23*pm1->x30;
	m.x21 = pm0->x20*pm1->x01	+ pm0->x21*pm1->x11	+ pm0->x22*pm1->x21	+ pm0->x23*pm1->x31;
	m.x22 = pm0->x20*pm1->x02	+ pm0->x21*pm1->x12	+ pm0->x22*pm1->x22	+ pm0->x23*pm1->x32;
	m.x23 = pm0->x20*pm1->x03	+ pm0->x21*pm1->x13	+ pm0->x22*pm1->x23	+ pm0->x23*pm1->x33;
	
	m.x30 = pm0->x30*pm1->x00	+ pm0->x31*pm1->x10	+ pm0->x32*pm1->x20	+ pm0->x33*pm1->x30;
	m.x31 = pm0->x30*pm1->x01	+ pm0->x31*pm1->x11	+ pm0->x32*pm1->x21	+ pm0->x33*pm1->x31;
	m.x32 = pm0->x30*pm1->x02	+ pm0->x31*pm1->x12	+ pm0->x32*pm1->x22	+ pm0->x33*pm1->x32;
	m.x33 = pm0->x30*pm1->x03	+ pm0->x31*pm1->x13	+ pm0->x32*pm1->x23	+ pm0->x33*pm1->x33;
	
	*pm0 = m;
}

void	M4f_trans		(tM4f* p, float x, float y, float z)
{
//	tM4f op = {1,0,0,0, 0,1,0,0, 0,0,1,0, x,y,z,1};
	tM4f op = {1,0,0,x, 0,1,0,y, 0,0,1,z, 0,0,0,1};
	M4f_mul_M4f (p, &op);
}
void	M4f_rotx	(tM4f* p, float a)
{
	tM4f op = {1,0,0,0, 0,cosf(a),sinf(a),0, 0,-sinf(a),cosf(a),0, 0,0,0,1};
	M4f_mul_M4f (p, &op);
}
void	M4f_roty	(tM4f* p, float a)
{
	tM4f op = {cosf(a),0,sinf(a),0, 0,1,0,0, -sinf(a),0,cosf(a),0, 0,0,0,1};
	M4f_mul_M4f (p, &op);
}
void	M4f_rotz	(tM4f* p, float a)
{
	tM4f op = {cosf(a),sinf(a),0,0, -sinf(a),cosf(a),0,0, 0,0,1,0, 0,0,0,1};
	M4f_mul_M4f (p, &op);
}

void	V4f_rotx	(tV4f* p, float a)
{
	tM4f op = {1,0,0,0, 0,cosf(a),sinf(a),0, 0,-sinf(a),cosf(a),0, 0,0,0,1};
	M4f_mul_V4f (&op, p);
}
void	V4f_roty	(tV4f* p, float a)
{
	tM4f op = {cosf(a),0,sinf(a),0, 0,1,0,0, -sinf(a),0,cosf(a),0, 0,0,0,1};
	M4f_mul_V4f (&op, p);
}
void	V4f_rotz	(tV4f* p, float a)
{
	tM4f op = {cosf(a),sinf(a),0,0, -sinf(a),cosf(a),0,0, 0,0,1,0, 0,0,0,1};
	M4f_mul_V4f (&op, p);
}



void	M4f_Frustrum		(tM4f* pproj, float w, float h, float n, float f)
{
	pproj->x00 = n / (w/2);
	pproj->x11 = n / (h/2);
	
	pproj->x22 = -(f+n) / (f-n);
	
	pproj->x23 = (-2*f*n) / (f-n);
	pproj->x32 = -1;
	
	printf ("\nM4f_Frustrum\n"); M4f_Print (pproj);
	printf ("\n");
}
void	M4f_Ortho			(tM4f* pproj, float w, float h, float n, float f)
{
	pproj->x00 = 1 / (w/2);
	pproj->x11 = 1 / (h/2);
	
	pproj->x22 = -2 / (f-n);
	pproj->x23 = -(f+n) / (f-n);
	
	pproj->x33 = 1;
	
}


void angle_test ()
{
	#define test(a0,a1)	printf ("test   " #a0 "\t\t" #a1 "\t\t%f\n", rad2deg*angle_diff_norm_pi_pi (a0*deg2rad, a1*deg2rad));
	
	test (0, 10);
	test (10, 0);
	
	test (-10, 0);
	test (0, -10);
	
	test (-10, 10);
	test (10, -10);
	
	test (0, 110);
	test (110, 0);
	
	test (180, 190);
	test (190, 180);
	
	test (170, 180);
	test (180, 170);
	
	test (170, 190);
	test (190, 170);
	
	#undef test
	exit (1);
}

void	Vec_Draw		(tCam* pcam, si x0, si y0, si x1, si y1)
{
	Ehh_Draw_Line_2d (pcam, x0, y0, x1, y1);
	
	si x, y;
	#define dd 2
	for (y = y1-dd; y <= y1+dd; ++y) {
		for (x = x1-dd; x <= x1+dd; ++x) {
			if (dpixout(x,y))
				continue;
			*(tRGBA*)dnpix(x,y) = gCol;
		}
	}
	#undef dd
}

void	V2f_DrawPosVec	(tCam* pcam, tV2f* ppos, tV2f* pvec)
{
	Vec_Draw (pcam, ppos->x, ppos->y, ppos->x+pvec->x, ppos->y+pvec->y);
}
void	V2f_DrawPosPos	(tCam* pcam, tV2f* ppos0, tV2f* ppos1)
{
	Vec_Draw (pcam, ppos0->x, ppos0->y, ppos1->x, ppos1->y);
}

void	V4f_DrawPosPos	(tCam* pcam, tV4f* ppos0, tV4f* ppos1)
{
	tV4f p0 = *ppos0, p1 = *ppos1;
	
//	V4f_mul_M4f (&p0, &pcam->World);	V4f_mul_M4f (&p1, &pcam->World);
//	V4f_mul_M4f (&p0, &pcam->Proj);	V4f_mul_M4f (&p1, &pcam->Proj);
	
	M4f_mul_V4f (&pcam->World, &p0);		M4f_mul_V4f (&pcam->World, &p1);
	M4f_mul_V4f (&pcam->Proj, &p0);		M4f_mul_V4f (&pcam->Proj, &p1);
	
	p0.x /= p0.w;
	p0.y /= p0.w;
	p0.z /= p0.w;
	
	p1.x /= p1.w;
	p1.y /= p1.w;
	p1.z /= p1.w;
	
	p0.x *= pcam->View_W/2;	p0.x += pcam->View_W/2;
	p0.y *= pcam->View_H/2;	p0.y += pcam->View_H/2;
	p1.x *= pcam->View_W/2;	p1.x += pcam->View_W/2;
	p1.y *= pcam->View_H/2;	p1.y += pcam->View_H/2;
	
//	printf ("p0 %f %f  p1 %f %f\n", p0.x, p0.y, p1.x, p1.y);
	Vec_Draw (pcam, p0.x, p0.y, p1.x, p1.y);
	
/*	printf ("p0 %f %f\n", p0.x, p0.y);
	
	CvMat* cm0 = cvCreateMatHeader(3, 1, CV_32FC1);
	cvCreateData(cm0);
	cvSet2D (cm0, 0, 0, cvScalarAll(ppos0->x));
	cvSet2D (cm0, 1, 0, cvScalarAll(ppos0->y));
	cvSet2D (cm0, 2, 0, cvScalarAll(ppos0->z));
	
	CvMat* cm_t = cvCreateMatHeader(3, 1, CV_32FC1);
	cvCreateData(cm_t);
	cvSet (cm_t, cvScalarAll(0), 0);
	CvMat* cm_r = cvCreateMatHeader(3, 1, CV_32FC1);
	cvCreateData(cm_r);
	cvSet (cm_r, cvScalarAll(0), 0);
	
	CvMat* cmret = cvCreateMatHeader(2, 1, CV_32FC1);
	cvCreateData(cmret);
	
	PrintMat (cm0);
	
	cvProjectPoints2(cm0, cm_r, cm_t, gM.Cam.cvCam, gM.Cam.cvDist, cmret, NULL, NULL, NULL, NULL, NULL, 1);
	
	printf ("c0 %f %f\n", cvGet2D(cmret, 0, 0).val[0], cvGet2D(cmret, 1, 0).val[0]);
	/**/
}


void	V4f_ScreenPosNorm	(tCam* pcam, tV4f* ppos0, tV2f* pret)
{
	tV4f p0 = *ppos0;
	if (gM.Eye_Line_Ray) {
	/*	float ax = M_PI_2 - gM.Cam.Image_FOV_W/2 + (peye->P.x/gM.Cam.Image_W)*gM.Cam.Image_FOV_W;
		float ay = M_PI_2 - gM.Cam.Image_FOV_H/2 + (peye->P.y/gM.Cam.Image_H)*gM.Cam.Image_FOV_H;
		
		ax = M_PI_2 - ax;
		ay = M_PI_2 - ay;
		
	//	printf ("Image_FOV   %f\n", gM.Cam.Image_FOV);
	//	printf ("Image_FOV_W %f\n", gM.Cam.Image_FOV_W*rad2deg);
	//	printf ("axy %f %f\n", ax*rad2deg, ay*rad2deg);
		V4f_roty (&p1, ax);
		V4f_rotx (&p1, ay);
		/**/
		
	//	p0.x = atan2()
		
	//	printf ("V4f_ScreenPosNorm: "); V4f_Print (ppos0);
		
		float lx = M_PI_2/* + gM.Cam.Image_FOV_W/2*/;
		float ly = M_PI_2/* + gM.Cam.Image_FOV_H/2*/;
		
		float ax = atan2(-ppos0->z, ppos0->x);
		float ay = atan2(-ppos0->z, ppos0->y);
		
	//	printf ("V4f_ScreenPosNorm: a %f %f   l %f %f\n", ax, ay, lx, ly); 
		
		p0.x = -2*(ax - lx) / pcam->Image_FOV_W;
		p0.y = -2*(ay - ly) / pcam->Image_FOV_H;
		
	//	printf ("V4f_ScreenPosNorm: return %f %f\n", p0.x, p0.y); 
	//	static int count = 5;
	//	if (!--count)
	//		exit(1);
	}else {
		M4f_mul_V4f (&pcam->World, &p0);
		M4f_mul_V4f (&pcam->Proj, &p0);
		
		p0.x /= p0.w;
		p0.y /= p0.w;
		p0.z /= p0.w;
	}
	pret->x = p0.x;	pret->y = p0.y;
}
void	V4f_ScreenPosf	(tCam* pcam, tV4f* ppos0, tV2f* pret)
{
	V4f_ScreenPosNorm (pcam, ppos0, pret);
	
	pret->x *= pcam->Image_W/2;	pret->x += pcam->Image_W/2;
	pret->y *= pcam->Image_H/2;	pret->y += pcam->Image_H/2;
}
void	V4f_ScreenPos	(tCam* pcam, tV4f* ppos0, tV2si* pret)
{
	tV4f p0 = *ppos0;
	
	M4f_mul_V4f (&pcam->World, &p0);
	M4f_mul_V4f (&pcam->Proj, &p0);
	
	p0.x /= p0.w;
	p0.y /= p0.w;
	p0.z /= p0.w;
	
	p0.x *= pcam->Image_W/2;	p0.x += pcam->Image_W/2;
	p0.y *= pcam->Image_H/2;	p0.y += pcam->Image_H/2;
	
	pret->x = p0.x;	pret->y = p0.y;
}




void	Dbg_V4f_DrawPosPos	(tDbg* pdbg, tV4f* ppos0, tV4f* ppos1)
{
	tV4f p0 = *ppos0, p1 = *ppos1;
	
	if (p0.w != 1 || p1.w != 1)
		return;
	assert (p0.w == 1);
	
	//printf ("4f p0 ");	V4f_Print (&p0);
	
	M4f_mul_V4f (&pdbg->World, &p0);		M4f_mul_V4f (&pdbg->World, &p1);
//	printf ("4f p0 ");	V4f_Print (&p0);
	
	M4f_mul_V4f (&pdbg->Proj, &p0);		M4f_mul_V4f (&pdbg->Proj, &p1);
//	printf ("4f p0 ");	V4f_Print (&p0);
	
	p0.x /= p0.w;
	p0.y /= p0.w;
	p0.z /= p0.w;
	
	p1.x /= p1.w;
	p1.y /= p1.w;
	p1.z /= p1.w;
	
//	printf ("p0 %f %f  p1 %f %f\n", p0.x, p0.y, p1.x, p1.y);
	
	p0.x *= pdbg->View_W/2;	p0.x += pdbg->View_X;
	p0.y *= pdbg->View_H/2;	p0.y += pdbg->View_Y;
	p1.x *= pdbg->View_W/2;	p1.x += pdbg->View_X;
	p1.y *= pdbg->View_H/2;	p1.y += pdbg->View_Y;
	
//	printf ("p0 %f %f  p1 %f %f\n", p0.x, p0.y, p1.x, p1.y);
	S_Draw_Line_2d (pdbg->SDL_Surf, p0.x, p0.y, p1.x, p1.y);
}

void	Dbg_V4f_DrawPoint		(tDbg* pdbg, tV4f* ppos0)
{
	tV4f p0 = *ppos0;
	
	if (p0.w != 1)
		return;
	assert (p0.w == 1);
	
//	printf ("4f p0 ");	V4f_Print (&p0);
	
	M4f_mul_V4f (&pdbg->World, &p0);
//	printf ("4f p0 ");	V4f_Print (&p0);
	
	M4f_mul_V4f (&pdbg->Proj, &p0);
//	printf ("4f p0 ");	V4f_Print (&p0);
	
	p0.x /= p0.w;
	p0.y /= p0.w;
	p0.z /= p0.w;
	
	p0.x *= pdbg->View_W/2;	p0.x += pdbg->View_X;
	p0.y *= pdbg->View_H/2;	p0.y += pdbg->View_Y;
	
//	printf ("p0 %f %f  p1 %f %f\n", p0.x, p0.y, p1.x, p1.y);
//	S_Draw_Line_2d (p0.x, p0.y, p1.x, p1.y);
	
	si x, y;
	#define dd 1
	for (y = p0.y-dd; y <= p0.y+dd; ++y) {
		for (x = p0.x-dd; x <= p0.x+dd; ++x) {
			if (x < 0 || x >= gM.pScreen->w || y < 0 || y >= gM.pScreen->h)
				continue;
			//dspix(x,y) = gColARGB;
		}
	}
	#undef dd
}


void	Dbg_M4f_Conf		(tDbg* pdbg)
{
	pdbg->View_X = pdbg->SDL_Surf->w/2;
	pdbg->View_Y = pdbg->SDL_Surf->w/2;
	pdbg->View_W = pdbg->SDL_Surf->w;
	pdbg->View_H = pdbg->SDL_Surf->h;
	
	M4f_Iden (&pdbg->Proj);
	
	M4f_Ortho (&pdbg->Proj, pdbg->Scale*pdbg->View_W, pdbg->Scale*pdbg->View_H, 0.1f, 100);
//	M4f_Frustrum (&pdbg->Proj, pdbg->View_W, pdbg->View_H, 0.1f, 100);
	
	M4f_Iden (&pdbg->World);
	
	M4f_trans (&pdbg->World, pdbg->T_X, pdbg->T_Y, -30);
	
	M4f_rotx (&pdbg->World, pdbg->R_X*deg2rad);
	M4f_roty (&pdbg->World, pdbg->R_Y*deg2rad);
	
//	M4f_rotx (&pdbg->World, 15*deg2rad);
//	M4f_trans (&pdbg->World, 0, 0, -15);
//	M4f_rotx (&pdbg->World, 15*deg2rad);
}



void	Dbg_V4f_ADrawPosPos	(tV4f* ppos0, tV4f* ppos1)
{
	for (int i = 0; i < M_Dbg_NUM; ++i) {
		Dbg_V4f_DrawPosPos (&gM.aDbg[i], ppos0, ppos1);
	}
}

void	Dbg_V4f_ADrawPoint	(tV4f* ppos0)
{
	for (int i = 0; i < M_Dbg_NUM; ++i) {
		Dbg_V4f_DrawPoint (&gM.aDbg[i], ppos0);
	}
}


float	M3f_Det	(tM3f* pm)
{
	return	pm->f[0][0]*pm->f[1][1]*pm->f[2][2]
			+ pm->f[0][1]*pm->f[1][2]*pm->f[2][0]
			+ pm->f[0][2]*pm->f[1][0]*pm->f[2][1]
			- pm->f[0][2]*pm->f[1][1]*pm->f[2][0]
			- pm->f[1][2]*pm->f[2][1]*pm->f[0][0]
			- pm->f[2][2]*pm->f[0][1]*pm->f[1][0];
}

void	M4f_Arg	(tM4f* pm, tM3f* pret, ui x, ui y)
{
	int i0,i1,j0,j1;
	for(i0=0, i1=0; i1 < 3; i0++,i1++)
	{
		if(i0 == y)
			++i0;
		for(j0=0,j1=0; j1 < 3; j0++,j1++)
		{
			if(j0 == x)
				++j0;
			pret->f[i1][j1] = pm->f[i0][j0];
		}
	}
}

float	M4f_Det	(tM4f* pm)
{
	float sum = 0;
	int i,j;
	for(i=0,j=0; j < 4; j++)
	{
	//	*n=m;
	//	arg(p,&d[0][0],n,i,j);
		tM3f arg;
		M4f_Arg (pm, &arg, j, i);
	//	sum=sum+*(p+10*i+j)*pow(-1,(i+j))*det(&d[0][0],n);
		sum += pm->f[i][j] * pow(-1, (i+j)) * M3f_Det (&arg);
	}
	return sum;
}

void	M4f_Inv	(tM4f* pm, tM4f* pret)
{
	int i,j;
	float d = M4f_Det (pm);
	if (d > -0.000001 && d < 0.000001) {
		printf("INVERSE DOES NOT EXIST");
	}
	
	for(i = 0; i < 4; i++) {
		for(j = 0; j < 4; j++)
		{
			tM3f arg;
			M4f_Arg (pm, &arg, j, i);
			pret->f[j][i] = pow(-1,(i+j)) * M3f_Det(&arg) / d;
		}
	}
	
}

/*
void	V4f_Intersect_Line01_Plane012	(
	tV4f* pret,
	tV4f l0, tV4f l1,
	tV4f p0, tV4f p1, tV4f p2
) {
	
}/**/

void	V4f_Intersect_Line01_Plane0N	(
	tV4f* pret,
	tV4f* pl0, tV4f* pl1,
	tV4f* pp0, tV4f* ppn
) {
	tV4f negw = *pp0;	V4f_sub_V4f (&negw, pl0);
	tV4f u = *pl1;	V4f_sub_V4f (&u, pl0);
	
	float s = V4f_dot_V4f (ppn, &negw) / V4f_dot_V4f (ppn, &u);
	
	u.x *= s;
	u.y *= s;
	u.z *= s;
	
	*pret = *pl0;
	V4f_add_V4f (pret, &u);
}


u08	V4f_Intersect_Line01_Tri012	(
	tV4f* ppos, float* ps, float* pt,
	tV4f* pl0, tV4f* pl1,
	tV4f* pt0, tV4f* pt1, tV4f* pt2
) {
	tV4f p0 = *pl0;
	tV4f p1 = *pl1;
	
	tV4f v0 = *pt0;
	tV4f v1 = *pt1;
	tV4f v2 = *pt2;
	
	tV4f u = v1;	V4f_sub_V4f (&u, &v0);
	tV4f v = v2;	V4f_sub_V4f (&v, &v0);
	tV4f w;
	
	{
		tV4f n = {0, 0, 1, 1};
		V4f_cross (&n, &u, &v);
		V4f_norm (&n);
		
		V4f_Intersect_Line01_Plane0N (
			&w,
			&p0, &p1,
			&v0, &n
		);
	}
	V4f_sub_V4f (&w, &v0);
	
/*	printf ("u ");	V4f_Print (&u);
	printf ("v ");	V4f_Print (&v);
	printf ("w ");	V4f_Print (&w);/**/
	
	float den = V4f_dot_V4f(&u, &v)*V4f_dot_V4f(&u, &v) - V4f_dot_V4f(&u, &u)*V4f_dot_V4f(&v, &v);
	
	float s = V4f_dot_V4f(&u, &v)*V4f_dot_V4f(&w, &v) - V4f_dot_V4f(&v, &v)*V4f_dot_V4f(&w, &u);
	s /= den;
	
	float t = V4f_dot_V4f(&u, &v)*V4f_dot_V4f(&w, &u) - V4f_dot_V4f(&u, &u)*V4f_dot_V4f(&w, &v);
	t /= den;
	
//	printf ("st %f\t%f\n", s, t);
	if (ps)
		*ps = s;
	if (pt)
		*pt = t;
	if (ppos) {
		ppos->x = v0.x + s*u.x + t*v.x;
		ppos->y = v0.y + s*u.y + t*v.y;
		ppos->z = v0.z + s*u.z + t*v.z;
		ppos->w = 1;
	}
}

//float	dist3D_Line_to_Line( Line L1, Line L2)
void	V4f_Intersect_Line01_Line01	(
	tV4f* ppos0, tV4f* ppos1,
	tV4f* pl00, tV4f* pl01,
	tV4f* pl10, tV4f* pl11
) {
//	Vector   u = L1.P1 - L1.P0;
//	Vector   v = L2.P1 - L2.P0;
//	Vector   w = L1.P0 - L2.P0;
	tV4f u = *pl01;	V4f_sub_V4f (&u, pl00);
	tV4f v = *pl11;	V4f_sub_V4f (&v, pl10);
	tV4f w = *pl00;	V4f_sub_V4f (&w, pl10);
	float a = V4f_dot_V4f (&u,&u);        // always >= 0
	float b = V4f_dot_V4f (&u,&v);
	float c = V4f_dot_V4f (&v,&v);        // always >= 0
	float d = V4f_dot_V4f (&u,&w);
	float e = V4f_dot_V4f (&v,&w);
	float D = a*c - b*b;       // always >= 0
	float sc, tc;
	
	// compute the line parameters of the two closest points
	if (D < dSMALL_NUM) {         // the lines are almost parallel
		sc = 0.0;
		tc = (b>c ? d/b : e/c);   // use the largest denominator
	}
	else {
		sc = (b*e - c*d) / D;
		tc = (a*e - b*d) / D;
	}
	if (ppos0) {
		*ppos0 = u;
		V4f_mul_S (ppos0, sc);
		V4f_add_V4f (ppos0, pl00);
	}
	if (ppos1) {
		*ppos1 = v;
		V4f_mul_S (ppos1, tc);
		V4f_add_V4f (ppos1, pl10);
	}
	// get the difference of the two closest points
//	Vector dP = w + (sc * u) - (tc * v);  // = L1(sc) - L2(tc)
	
//	return norm(dP);   // return the closest distance
}


void	V4f_Intersect_Sphere0R_Line01	(
	tV4f* pret,
	tV4f* ppos, float r,
	tV4f* pl0, tV4f* pl1
) {
	tV4f vec;
	vec = *pl1; V4f_add_V4f (&vec, pl0);
	
	float a = V4f_dist2 (&vec);
	float b = 2*( (pl1->x-pl0->x)*(pl0->x-ppos->x)
				+ (pl1->y-pl0->y)*(pl0->y-ppos->y)
				+ (pl1->z-pl0->z)*(pl0->z-ppos->z)
	);
	float c = V4f_dist2 (ppos) + V4f_dist2 (pl0) - 2*(ppos->x*pl0->x + ppos->y*pl0->y + ppos->z*pl0->z) - dpow2(r);
	
	float d = b*b - 4*a*c;
	if (d < 0) {
		pret->x = 0;
		pret->y = 0;
		pret->z = 0;
		pret->w = 1;
		return;
	}
//	printf ("d %f\n", d);
	float u;
	
	u = (-b - sqrt(d)) / (2*a);
	pret->x = pl0->x + u*vec.x;
	pret->y = pl0->y + u*vec.y;
	pret->z = pl0->z + u*vec.z;
	pret->w = 1;
//	printf ("pret0 "); V4f_Print (pret);
	
	u = (-b + sqrt(d)) / (2*a);
	pret->x = pl0->x + u*vec.x;
	pret->y = pl0->y + u*vec.y;
	pret->z = pl0->z + u*vec.z;
	pret->w = 1;
//	printf ("pret1 "); V4f_Print (pret);
	
	u = (-b - sqrt(d)) / (2*a);
	pret->x = pl0->x + u*vec.x;
	pret->y = pl0->y + u*vec.y;
	pret->z = pl0->z + u*vec.z;
	pret->w = 1;
//	printf ("pret0 "); V4f_Print (pret);/**/
}


struct {
	int state;
}Mouse;

void mouse_coords (int *x, int *y)
{
	XEvent event;
	XQueryPointer (gM.X.pDisp, DefaultRootWindow (gM.X.pDisp),
	&event.xbutton.root, &event.xbutton.window,
	&event.xbutton.x_root, &event.xbutton.y_root,
	&event.xbutton.x, &event.xbutton.y,
	&event.xbutton.state);
	*x = event.xbutton.x;
	*y = event.xbutton.y;
}

void mouse_move (int dx, int dy)
{
	XWarpPointer (gM.X.pDisp, None, None, 0,0,0,0, dx, dy);
	XFlush (gM.X.pDisp);
	usleep (1);
}

void setup_event(XEvent * event, int button)
{
	memset(event, 0, sizeof (event));
	
	event->xbutton.button = button;
	event->xbutton.same_screen = True;
	event->xbutton.subwindow = DefaultRootWindow(gM.X.pDisp);
	
	while (event->xbutton.subwindow)
	{
		event->xbutton.window = event->xbutton.subwindow;
		XQueryPointer (gM.X.pDisp, event->xbutton.window,
			&event->xbutton.root, &event->xbutton.subwindow,
			&event->xbutton.x_root, &event->xbutton.y_root,
			&event->xbutton.x, &event->xbutton.y,
			&event->xbutton.state
		);
	}
}

void mouse_press	(int button)
{
	XEvent event;
	int mask = 0;
	switch (button) {
	case 1:	//mask = Button1Mask;	break;
			system("xdotool mousedown 1");	return;
	
	case 2:	//mask = Button2Mask;	break;
			system("xdotool mousedown 2");	return;
	case 3:	//mask = Button3Mask;	break;
			system("xdotool mousedown 3");	return;
	}
//	printf ("mouse press %x %x\n", Mouse.state, mask);
	if (Mouse.state & mask)
		return;
	
	setup_event(&event, button);
	event.type = ButtonPress;
	XSetInputFocus (gM.X.pDisp, event.xbutton.window, RevertToParent, CurrentTime);
	
	event.xbutton.state = Mouse.state;
	
	XSendEvent (gM.X.pDisp, event.xbutton.window, True, ButtonPressMask, &event);
	XFlush (gM.X.pDisp);
	
	Mouse.state |= mask;
	
	usleep (1);
}

void mouse_release	(int button)
{
	XEvent event;
	int mask = 0;
	switch (button) {
	case 1:	//mask = Button1Mask;	break;
			system("xdotool mouseup 1");	return;
	
	case 2:	//mask = Button2Mask;	break;
			system("xdotool mouseup 2");	return;
	case 3:	//mask = Button3Mask;	break;
			system("xdotool mouseup 3");	return;
	}
//	printf ("mouse_release %x %x\n", Mouse.state, mask);
	if (Mouse.state & mask == 0)
		return;
	
	setup_event(&event, button);
	event.type = ButtonRelease;
	XSetInputFocus(gM.X.pDisp, event.xbutton.window, RevertToParent, CurrentTime);
	
	event.xbutton.state = Mouse.state;
	
	XSendEvent (gM.X.pDisp, event.xbutton.window, True, ButtonReleaseMask, &event);
	XFlush (gM.X.pDisp);
	
	Mouse.state &= ~mask;
	usleep (1);
}


s08	Win_PropSet_Cardinal32	(Window win, Atom prop, int val)
{
	if (XChangeProperty (gM.X.pDisp, win,
			prop,		//Atom		// property
			XA_CARDINAL,			//Atom		// type
			32,				//int		// format
			PropModeReplace,		//int		// mode
			(unsigned char *)&val,		//_Xconst unsigned char*	/* data
			1				//int			/* nelements
		) != Success
	) {
		return False;
	}
	return True;
}
s08	Win_PropGet_Cardinal32	(Window win, Atom prop, int* pval)
{
	Atom actual_type_return;
	int actual_format_return;
	unsigned long nitems_return, bytes_after_return;
	
	unsigned char *prop_return;
	
	if (XGetWindowProperty (gM.X.pDisp, win,
			prop,			//Atom property;
			0, 1,			//long long_offset, long_length;
			False,			//Bool delete;
			XA_CARDINAL,		//Atom req_type;
			&actual_type_return,	//Atom *actual_type_return;
			&actual_format_return,	//int *actual_format_return;
			&nitems_return,		//unsigned long *nitems_return;
			&bytes_after_return,	//unsigned long *bytes_after_return;
			&prop_return		//unsigned char **prop_return;
		) != Success
	) {
		return False;
	}
	if (actual_type_return != XA_CARDINAL) {
		return False;
	}
	*pval = *(int*)prop_return;
	return True;
}





void	Marker_Show		(tMarker* p)
{
	XMapWindow(gM.X.pDisp, p->Win);
	XFlush(gM.X.pDisp);
	p->Dirty = 1;
}
void	Marker_Hide		(tMarker* p)
{
	XUnmapWindow(gM.X.pDisp, p->Win);
	XFlush(gM.X.pDisp);
}

void	Marker_Move		(tMarker* p, si x, si y)
{
	p->Dirty = 1;
	XMoveWindow (gM.X.pDisp, p->Win,
		x - p->Win_W/2,
		y - p->Win_W/2
	);
	XFlush(gM.X.pDisp);
	XClearArea (gM.X.pDisp, p->Win, 0, 0, 0, 0, True);
}

void	Marker_Paint	(tMarker* p)
{
//	if (!p->Dirty)
//		return;
	p->Dirty = 0;
	
	int blackColor = BlackPixel (gM.X.pDisp, DefaultScreen(gM.X.pDisp));
	int whiteColor = WhitePixel (gM.X.pDisp, DefaultScreen(gM.X.pDisp));
	
//	XSelectInput(gM.X.pDisp, p->Win, StructureNotifyMask);
	// "Map" the window (that is, make it appear on the screen)
//	XMapWindow(gM.X.pDisp, w);
	
/*	for(;;) {
	XEvent e;
	XNextEvent(gM.X.pDisp, &e);
	if (e.type == MapNotify)
	break;
	}/**/
	
	GC gc = XCreateGC(gM.X.pDisp, p->Win, 0, NULL);
	
	switch (p->Type) {
	case 0: {
		XSetForeground(gM.X.pDisp, gc, p->Col);
		
		XDrawRectangle (gM.X.pDisp, p->Win, gc, 1, 1, p->Win_W-3, p->Win_W-3);
		
	//	XDrawLine (gM.X.pDisp, p->Win, gc, 1, 1, p->Win_W-2, 1);
	//	XDrawLine (gM.X.pDisp, p->Win, gc, 1, 1, 1, p->Win_W-2);
	//	XDrawLine (gM.X.pDisp, p->Win, gc, 1, p->Win_W-2, p->Win_W-2, p->Win_W-2);
	//	XDrawLine (gM.X.pDisp, p->Win, gc, p->Win_W-2, 1, p->Win_W-2, p->Win_W-2);
		break;
	}
	case 1: {
		XSetForeground(gM.X.pDisp, gc, p->Col);
	//	XSetBackground(gM.X.pDisp, gc, blackColor);
		
		XFillRectangle (gM.X.pDisp, p->Win, gc, 0, 0, p->Win_W, p->Win_W);
		
		XSetForeground(gM.X.pDisp, gc, blackColor);
		
		XDrawPoint (gM.X.pDisp, p->Win, gc, p->Win_W/2, p->Win_W/2);
		
		XDrawRectangle (gM.X.pDisp, p->Win, gc, 1, 1, p->Win_W-3, p->Win_W-3);
		break;
	}
	}
	XFlush(gM.X.pDisp);
	
	XFreeGC (gM.X.pDisp, gc);
}

void	Marker_ColorSet	(tMarker* p, u32 col)
{
	if (col == p->Col)
		return;
	p->Dirty = 1;
	p->Col = col;
	XClearArea (gM.X.pDisp, p->Win, 0, 0, 0, 0, True);
}

void	Marker_Init	(tMarker* p, u08 type)
{
	p->Type = type;
	p->Col = 0xFFFFFFFF;
	switch (p->Type) {
	case 0:	p->Win_W = 31;	break;
	case 1:	p->Win_W = 21;	break;
	}
	
/*	p->Win = XCreateWindow (
			gM.X.pDisp,
			gM.aScreen[0].Win,
			0, 0, p->Win_W, p->Win_W,
			0,
			CopyFromParent,
			InputOutput,
			CopyFromParent,
			0,
			NULL
	);/**/
	
	p->Win = XCreateSimpleWindow (
			gM.X.pDisp,
			gM.aScreen[0].Win,
			1300, 900, p->Win_W, p->Win_W,
			0, 0,
			0
	);/**/
	
	switch (p->Type) {
	case 0: {
		XRectangle arect[5];
		arect[0].x = 0;	arect[0].y = 0;	arect[0].width = p->Win_W;	arect[0].height = 2;
		arect[1].x = 0;	arect[1].y = 0;	arect[1].width = 2;		arect[1].height = p->Win_W;
		
		arect[2].x = 0;			arect[2].y = p->Win_W-2;	arect[2].width = p->Win_W;	arect[2].height = 2;
		arect[3].x = p->Win_W-2;	arect[3].y = 0;			arect[3].width = 2;		arect[3].height = p->Win_W;
		
		arect[4].x = p->Win_W/2;	arect[4].y = p->Win_W/2;	arect[4].width = 1;		arect[4].height = 1;
		XShapeCombineRectangles (gM.X.pDisp, p->Win,
			ShapeBounding,
		//	ShapeClip,
		//	ShapeInput,
			0, 0,	//offset
			arect,
			5,
			ShapeSet,
			0
		);
		break;
	}
	case 1: {
		break;
	}
	}
	XFlush(gM.X.pDisp);
	
	{
		XSetWindowAttributes attr;
		attr.event_mask = ExposureMask;
		
		XChangeWindowAttributes (gM.X.pDisp, p->Win,
			CWEventMask,
			&attr
		);
	}/**/
//	XStoreName (gM.X.pDisp, p->Win, "fuuuuu");
	{
		XClassHint class_hint;
		class_hint.res_name = "luvcview";
		class_hint.res_class = class_hint.res_name;
		XSetClassHint(gM.X.pDisp, p->Win, &class_hint);
	}/**/
	{
		Atom prop = XInternAtom (gM.X.pDisp, "_NET_WM_STATE", True);
		Atom data = XInternAtom (gM.X.pDisp, "_NET_WM_STATE_ABOVE", True);
		
		XChangeProperty (gM.X.pDisp, p->Win,
			prop,				//Atom		// property
			XA_ATOM,			//Atom		// type
			32,				//int			// format
			PropModeReplace,		//int			// mode
			(unsigned char *)&data,			//_Xconst unsigned char*	/* data
			1				//int			/* nelements
		);
	}/**/
	{
		Atom prop = XInternAtom (gM.X.pDisp, "_NET_WM_WINDOW_TYPE", True);
		Atom data = XInternAtom (gM.X.pDisp, "_NET_WM_WINDOW_TYPE_DOCK", True);
		
		XChangeProperty (gM.X.pDisp, p->Win,
			prop,				//Atom		// property
			XA_ATOM,			//Atom		// type
			32,				//int			// format
			PropModeReplace,		//int			// mode
			(unsigned char *)&data,			//_Xconst unsigned char*	/* data
			1				//int			/* nelements
		);
	}/**/
	Win_PropSet_Cardinal32 (p->Win, XInternAtom (gM.X.pDisp, "_NET_WM_DESKTOP", True), -1);
	{
		typedef struct {//WHAT THE FUUUUUCK IS WROOOONG WITH YOUR FAAAAACEEEEEE
			unsigned long   flags;
			unsigned long   functions;
			unsigned long   decorations;
			long            inputMode;
			unsigned long   status;
		} Hints;
		Hints   hints;
		memset (&hints, 0, sizeof(hints));
		hints.flags = 2;        // Specify that we're changing the window decorations.
		hints.decorations = 0;  // 0 (false) means that window decorations should go bye-bye.
		
		Atom prop = XInternAtom(gM.X.pDisp, "_MOTIF_WM_HINTS", True);
		
		XChangeProperty (gM.X.pDisp, p->Win,
			prop,				//Atom		// property
			prop,				//Atom		// type
			32,				//int			// format
			PropModeReplace,		//int			// mode
			(unsigned char *)&hints,			//_Xconst unsigned char*	/* data
			5				//int			/* nelements
		);
	}/**/
	
//	XMapWindow(gM.X.pDisp, p->Win);
	
	XFlush(gM.X.pDisp);
//	XSync(gM.X.pDisp,False);
	
	p->Dirty = 1;
}



tV2f	point_clip	(tV2f* ppoint)
{
	if (ppoint->x < 0)
		ppoint->x = 0;
	else if (ppoint->x >= pcam->Image_W)
		ppoint->x = pcam->Image_W-1;
		
	if (ppoint->y < 0)
		ppoint->y = 0;
	else if (ppoint->y >= pcam->Image_W)
		ppoint->y = pcam->Image_W-1;
}



typedef struct {
	double A, B, C, D, E;
	
	ui Root_N;
	double aRoot[4];
}tQuart;

double	Quart_F	(tQuart* p, double x)
{
	return p->A*x*x*x*x
		+ p->B*x*x*x
		+ p->C*x*x
		+ p->D*x
		+ p->E;
}

void		Quart_NumSolve	(tQuart* p, double err)	//room not hot enough
{
	p->Root_N = 0;
	double b = -1, e = 1, d = 0.01;
	double x0, x1;
	for (x0 = b; x0 < e;) {
		x1 = x0 + d;
		if (signbit(Quart_F(p, x0)) != signbit(Quart_F(p, x1))) {
			double xn = x1;
			while (x1 - x0 >= err) {
				double xi = (x0+x1)/2;
				if (signbit(Quart_F(p, x0)) == signbit(Quart_F(p, xi)))
					x0 = xi;
				else
					x1 = xi;
			}
			p->aRoot[p->Root_N++] = (x0+x1)/2;
			if (p->Root_N == 4)
				return;
			x0 = xn;
		}else {
			x0 = x1;
		}
	}
	return;
}
void		Quart_AlgSolve	(tQuart* p)	//screw this crap
{
	double A = p->A, B = p->B, C = p->C, D = p->D, E = p->E;
	
	printf ("A %f\n", A);
	printf ("B %f\n", B);
	printf ("C %f\n", C);
	printf ("D %f\n", D);
	printf ("E %f\n", E);
	
	double alpha = -(3*B*B) / (8*A*A)
				+ C/A;
	double beta = B*B*B / (8*A*A*A)
				- B*C / (2*A*A)
				+ D/A;
	double gamma = -3*B*B*B*B / (256*A*A*A*A)
				+ C*B*B / (16*A*A*A)
				- B*D / (4*A*A)
				+ E/A;
	
	printf ("alpha %f\n", alpha);
	printf ("beta %f\n", beta);
	printf ("gamma %f\n", gamma);
	if (fabs(beta) < dSMALL_NUM) {
		p->aRoot[0] = -B/(4*A) + sqrt(-alpha + sqrt(alpha*alpha - 4*gamma) / 2);
		p->aRoot[1] = -B/(4*A) + sqrt(-alpha - sqrt(alpha*alpha - 4*gamma) / 2);
		p->aRoot[2] = -B/(4*A) - sqrt(-alpha + sqrt(alpha*alpha - 4*gamma) / 2);
		p->aRoot[3] = -B/(4*A) - sqrt(-alpha - sqrt(alpha*alpha - 4*gamma) / 2);
		return;
	}
	double P = -alpha*alpha/12 - gamma;
	double Q = -alpha*alpha*alpha/108 + alpha*gamma/3 - beta*beta/8;
	double R = -Q/2 + sqrtf(Q*Q/4 + P*P*P/27);
	printf ("Q*Q/4 + P*P*P/27 %f\n", Q*Q/4 + P*P*P/27);
	if (fabs(Q*Q/4 + P*P*P/27) <= dSMALL_NUM) {
		printf ("zero\n");
		R = -Q/2;
	}
	double U = cbrt(R);
	
	printf ("P %f\n", P);
	printf ("Q %f\n", Q);
	printf ("R %f\n", R);
	printf ("U %f\n", U);
	
	double yy;
	if (fabs(U) <= dSMALL_NUM) {
	}else {
		yy = -5.0/6.0 * alpha + U - P/(3*U);
	}
	double W = sqrt(alpha + 2*yy);
	
	p->aRoot[0] = -B/(4*A) + (W + sqrt(-(3*alpha + 2*yy + 2*beta/W))) / 2;
	p->aRoot[1] = -B/(4*A) + (W - sqrt(-(3*alpha + 2*yy + 2*beta/W))) / 2;
	p->aRoot[2] = -B/(4*A) + (-W + sqrt(-(3*alpha + 2*yy - 2*beta/W))) / 2;
	p->aRoot[3] = -B/(4*A) + (-W - sqrt(-(3*alpha + 2*yy - 2*beta/W))) / 2;
	
	return;
}


void	Sphere_Reflect	(tV4f* pret, tV4f* ppos, float r, tV4f* pl)	//Not yet quite what it should be
{
	si i;
	tV2f N, S = {-ppos->x, -ppos->z}, L = {pl->x-ppos->x, pl->z-ppos->z};
//	printf ("S %f %f\n", S.x, S.y);
//	printf ("L %f %f\n", L.x, L.y);
	
	V2f_mul_S (&S, 1.0f/r);
	V2f_mul_S (&L, 1.0f/r);
	{
		double a = V2f_dot_V2f (&S, &S);
		double b = V2f_dot_V2f (&S, &L);
		double c = V2f_dot_V2f (&L, &L);
		
	//	printf ("abc %f %f %f\n", a, b, c);
		
	//	printf ("a*c - b*b = %f\n", a*c - b*b); 
		tQuart poly;
		poly.A = 4*c*(a*c - b*b);
		poly.B = -4*(a*c - b*b);
		poly.C = (a + 2*b + c - 4*a*c);
		poly.D = 2*(a-b);
		poly.E = a - 1;
		
	//	Quart_AlgSolve (&poly);
		Quart_NumSolve (&poly, 0.000001);
		
		for (i = 0; i < 4; ++i) {
		//	printf ("poly.aRoot[%d] = %f", i, poly.aRoot[i]);
			
		//	printf ("\treally? %f\t", Quart_F(&poly, poly.aRoot[i]));
			if (poly.aRoot[i] > 0) {
				double ss = (-2*c*poly.aRoot[i]*poly.aRoot[i] + poly.aRoot[i] + 1) / (2*b*poly.aRoot[i] + 1);
		//		printf ("\tss = %f", ss);
				if (ss >= 0) {
					N.x = ss*S.x + poly.aRoot[i]*L.x;
					N.y = ss*S.y + poly.aRoot[i]*L.y;
					
				}
			}
		//	printf ("\n");
		}
	}
//	printf ("N %f %f\n", N.x, N.y);
	
	pret->x = ppos->x + N.x * r;
	pret->y = 0;
	pret->z = ppos->z + N.y * r;
	pret->w = 1;
	
	V4f_rotx (pret, -atan2(ppos->y, -ppos->z));
	
}


float	NN_Sphere	(tV2si* ppos, float ir, float or, u08 dark, u08 bright)
{
	si x, y;
	float iay = 0, oay = 0;
	ui in = 0, on = 0;
	for (y = ppos->y-or; y <= ppos->y+or; ++y) {
		for (x = ppos->x-or; x <= ppos->x+or; ++x) {
			if (ddist2(x,y,ppos->x,ppos->y) <= dpow2(ir)) {
				iay += dopix(x,y)->Y;
				++in;
			}else if (ddist2(x,y,ppos->x,ppos->y) <= dpow2(or)) {
				oay += dopix(x,y)->Y;
				++on;
			}
		}
	}
	iay /= in;
	oay /= on;
	
	return oay - iay;
}



void	Cam_Pos2Ray	(tV2f pos, tV4f* pp0, tV4f* pp1, tV4f* pvec)
{
	abort();
	#if 0
	tV4f p0 = {0, 0, 0, 1};
	tV4f p1 = {0, 0, -1, 1};
	
	if (gM.Eye_Line_Ray) {
		float ax = M_PI_2 - gM.Cam.Image_FOV_W/2 + (pos.x/gM.Cam.Image_W)*gM.Cam.Image_FOV_W;
		float ay = M_PI_2 - gM.Cam.Image_FOV_H/2 + (pos.y/gM.Cam.Image_H)*gM.Cam.Image_FOV_H;
		
		ax = M_PI_2 - ax;
		ay = M_PI_2 - ay;
		
	//	printf ("Image_FOV   %f\n", gM.Cam.Image_FOV);
	//	printf ("Image_FOV_W %f\n", gM.Cam.Image_FOV_W*rad2deg);
	//	printf ("axy %f %f\n", ax*rad2deg, ay*rad2deg);
		V4f_roty (&p1, ax);
		V4f_rotx (&p1, ay);
	}else {
		p1.x = pcam->Proj_L + (pos.x/gM.Cam.Image_W)*pcam->Proj_W;
		p1.y = pcam->Proj_B + (pos.y/gM.Cam.Image_H)*pcam->Proj_H;
	}
	V4f_mul_S (&p1, 1.0/V4f_dist(&p1));
//	printf ("Cam_Pos2Ray:	");	V4f_Print (&p1);
	if (pp0)
		*pp0 = p0;
	if (pp1)
		*pp1 = p1;
	if (pvec) {
		*pvec = p1;
	}
	#endif
}


void	Cam_Retina_Ray	(tEye* peye, tV4f* pp0, tV4f* pp1, tV4f* pvec)
{
	abort();
	//Cam_Pos2Ray (peye->P, pp0, pp1, pvec);
}


void	Cam_Param_Set	(tCam* pcam)
{
	si value;
	#define dcam_set(_id,_val)		\
		do {		\
			struct v4l2_control control;		\
			control.id    = _id;		\
			control.value = _val;		\
			if ((value = ioctl(pcam->UVC->fd, VIDIOC_S_CTRL, &control)) < 0)		\
				printf("Set " #_id " error\n");		\
		/*	else	printf(#_id " set to %d\n", control.value);/**/		\
		}while (0)
	
	dcam_set (V4L2_CID_AUTO_WHITE_BALANCE, 1);
	dcam_set (V4L2_CID_AUTOGAIN, 1);
	
	if (pcam->Focus == -1)
		dcam_set(V4L2_CID_FOCUS_AUTO, 1);
	else {
		dcam_set(V4L2_CID_FOCUS_AUTO, 0);
		
		if ((value = v4l2SetControl(pcam->UVC, V4L2_CID_FOCUS_ABSOLUTE, 1)) < 0)
			printf("Set CT_FOCUS_ABSOLUTE_CONTROL to %ld error\n", value);
		if ((value = v4l2SetControl(pcam->UVC, V4L2_CID_FOCUS_ABSOLUTE, pcam->Focus)) < 0)
			printf("Set CT_FOCUS_ABSOLUTE_CONTROL to %ld error\n", value);
	}
	printf ("V4L2_CID_FOCUS_ABSOLUTE %d\n", v4l2GetControl(pcam->UVC, V4L2_CID_FOCUS_ABSOLUTE));
	
	if (pcam->Exposure == -1) {
		//V4L2_EXPOSURE_SHUTTER_PRIORITY	V4L2_EXPOSURE_APERTURE_PRIORITY
		dcam_set(V4L2_CID_EXPOSURE_AUTO, V4L2_EXPOSURE_AUTO);
		dcam_set(V4L2_CID_EXPOSURE_AUTO_PRIORITY, V4L2_EXPOSURE_APERTURE_PRIORITY);
		dcam_set(V4L2_CID_EXPOSURE_AUTO, V4L2_EXPOSURE_AUTO);
	//	if ((value = v4l2SetControl(videoIn, V4L2_CID_EXPOSURE_AUTO, V4L2_EXPOSURE_AUTO)) < 0)
	//		printf("Set V4L2_CID_EXPOSURE_AUTO to %ld error\n", value);
	}else {
		if ((value = v4l2SetControl(pcam->UVC, V4L2_CID_EXPOSURE_AUTO, V4L2_EXPOSURE_MANUAL)) < 0)
			printf("Set V4L2_CID_EXPOSURE_AUTO to %ld error\n", value);
	//	dcam_set(V4L2_CID_EXPOSURE_AUTO, V4L2_EXPOSURE_MANUAL);
		
	//	if ((value = v4l2SetControl(videoIn, V4L2_CID_EXPOSURE_ABSOLUTE, pcam->Exposure-30)) < 0)
	//		printf("Set V4L2_CID_EXPOSURE_ABSOLUTE to %ld error\n", value);
		if ((value = v4l2SetControl(pcam->UVC, V4L2_CID_EXPOSURE_ABSOLUTE, pcam->Exposure)) < 0)
			printf("Set V4L2_CID_EXPOSURE_ABSOLUTE to %ld error\n", value);
	}
	if ((value = v4l2SetControl(pcam->UVC, V4L2_CID_ZOOM_ABSOLUTE, pcam->Zoom)) < 0)
		printf("Set V4L2_CID_ZOOM_ABSOLUTE to %ld error\n", value);
	#undef dcam_set
}




#if 0
void	Track_Point_Init		(tTrack_Point* p)
{
	p->P.x = 0;
	p->P.y = 0;
	
	p->W = p->H = 64;
	
	p->pANN = 0;
}

void	Track_Point_Conf		(tTrack_Point* p)
{
	p->P.x = gM.Cam.Image_W/2;
	p->P.y = gM.Cam.Image_H/2;
	
}



#if 0
TrainUtil_Example*	Track_Point_ExampleCreate	(tTrack_Point* p, float ox, float oy, float decision)
{
	TrainUtil_Example* example = (TrainUtil_Example*) malloc(sizeof *example);
	
	example->length = p->W*p->H;
	example->decision = decision;
	example->attrs = (float*)malloc(sizeof(float)*example->length);
	
	si idx = 0, yi, xi;
	for (yi = oy - p->H/2; yi < oy + p->H/2; yi++) {
		for (xi = ox - p->W/2; xi < ox + p->W/2; xi++) {
		//	int image_idx = yi*image->widthStep + xi;
		//	example->attrs[idx] = (unsigned char) image->imageData[image_idx];
			example->attrs[idx] = dopix(xi,yi)->Y;
			idx++;
		}
	}
	
	assert(example->length==idx);
	assert(example->length==64*64);
	
	return example;
}



void	Track_Point_Train		(tTrack_Point* p)
{
	
	TrainUtil_ExampleList *training_set = NULL;
	
//	for (int xi=-1; xi <= 1; xi++) {
//		for (int yi=-1; yi <= 1; yi++) {
			TrainUtil_Example *example = Track_Point_ExampleCreate (p, p->P.x, p->P.y, 1.0);
			TrainUtil_ExampleList_Add (&training_set, example);
//		}
//	}
	
	si n_points = 0;
	si margin = p->W / 2;
	si dead_zone = 6;
	
	while (n_points < 16) {
		int xi = dead_zone + rand()%margin;
		xi *= -1 + 2 * (rand()%2);
		
		int yi = dead_zone + rand()%margin;
		yi *= -1 + 2 * (rand()%2);
		
		//check if point lies outside dead zone
		if (xi > -dead_zone && xi < dead_zone && xi > -dead_zone && yi < dead_zone)
			continue;
		
	//	TrainUtil_ExampleParams mod_params = params;
	//	mod_params.x1 += xi;
	//	mod_params.x2 += xi;
	//	mod_params.y1 += yi;
	//	mod_params.y2 += yi;
		
	//	IplImage *image_sample = TrainUtil_SampleImage(gray_image, &mod_params);
	//	TrainUtil_Example *example = TrainUtil_CreateExample(image_sample, -1.0);
		TrainUtil_Example *example = Track_Point_ExampleCreate (p, p->P.x + xi, p->P.y + yi, -1.0);
		TrainUtil_ExampleList_Add(&training_set, example);
		n_points++;
	}
	
	//to sort we have to convert list to array
	int X_size = TrainUtil_ExampleList_Length(training_set);
	TrainUtil_Example **training_arr = (TrainUtil_Example**) malloc(sizeof(TrainUtil_Example*)*X_size);
	
	int idx = 0;
	while (training_set) {
		training_arr[idx] = training_set->example;
		training_set = training_set->next;
		idx++;
	}
	
	TrainUtil_Classifier *out_cls;
	
	TrainUtil_GentleBoost(
		//input:
		training_arr,
		X_size,
		training_arr,
		X_size,
		66,
		//output:
		&out_cls);
	/**/
	
//	si i;
//	for (i = 0; i < 66; ++i)
//		printf ("Track_Point_Train	index %d	error %f	a %f	b %f	th %f\n", out_cls[i].index, out_cls[i].error, out_cls[i].a, out_cls[i].b, out_cls[i].th);
	
	/*
	while (fgets(buf, 999, f)) {
		int len = (int) strlen(buf);
		while (len > 0 && isspace(buf[len-1]))
			len--;
		buf[len] = '\0';
		
		IplImage *image = cvLoadImage( buf, 1 );
		
		if (image) {
			//convert_image to gray
			IplImage* gray_image = cvCreateImage(cvSize(image->width,image->height), 8, 1);
			cvCvtColor( image, gray_image, CV_BGR2GRAY );
			cvEqualizeHist(gray_image, gray_image);
			//free rgb image
			cvReleaseImage(&image);
			
			TrainUtil_ExampleParams params;
			
			//open image in new window in order to select desired feature
			cvNamedWindow("original", 1);
			cvSetMouseCallback("original", TrainUtil_MouseHandler, &params );
			
			cvShowImage("original", gray_image);
			if (cvWaitKey(0) == 'e')
				exit(1);
			
			cvDestroyWindow("original");
			
			cvNamedWindow("test", 1);
			//create sample from mouse-selected area
			//mark 9 positive examples
			for (int xi=-1; xi <= 1; xi++) {
				for (int yi=-1; yi <= 1; yi++) {
					TrainUtil_ExampleParams mod_params = params;
					mod_params.x1 += xi;
					mod_params.x2 += xi;
					mod_params.y1 += yi;
					mod_params.y2 += yi;
					
					IplImage *image_sample = TrainUtil_SampleImage(gray_image, &mod_params);
					cvShowImage("test", image_sample);
					cvWaitKey(0);
				//	system("sleep 5");
					TrainUtil_Example *example = TrainUtil_CreateExample(image_sample, 1.0);
					TrainUtil_ExampleList_Add(&training_set, example);
					cvReleaseImage(&image_sample);
				}
			}
			cvDestroyWindow("test");
			
			//mark 16 negative points
			int n_points = 0;
			int margin = 32;
			int dead_zone = 6;
			
			while (n_points < 16) {
				int xi = dead_zone + rand() % margin;
				xi *= -1 + 2 * (rand() % 2);
				
				int yi = dead_zone + rand() % margin;
				yi *= -1 + 2 * (rand() % 2);
				//check if point lies outside dead zone
				if (xi > -dead_zone && xi < dead_zone
				    && xi > -dead_zone && yi < dead_zone)
					continue;
					
				TrainUtil_ExampleParams mod_params = params;
				mod_params.x1 += xi;
				mod_params.x2 += xi;
				mod_params.y1 += yi;
				mod_params.y2 += yi;
				
				IplImage *image_sample = TrainUtil_SampleImage(gray_image, &mod_params);
				TrainUtil_Example *example = TrainUtil_CreateExample(image_sample, -1.0);
				TrainUtil_ExampleList_Add(&training_set, example);
				cvReleaseImage(&image_sample);
				n_points++;
			}
			
			cvReleaseImage(&gray_image);
			n_lines++;
		}
	}

	TrainUtil_Classifier *out_cls;
	
	//to sort we have to convert list to array
	int X_size = TrainUtil_ExampleList_Length(training_set);
	TrainUtil_Example **training_arr = (TrainUtil_Example**) malloc(sizeof(TrainUtil_Example*)*X_size);
	
	int idx = 0;
	while (training_set) {
		training_arr[idx] = training_set->example;
		training_set = training_set->next;
		idx++;
	}
	
	TrainUtil_GentleBoost(
		//input:
		training_arr,
		X_size,
		training_arr,
		X_size,
		66,
		//output:
		&out_cls);
	/**/
	
	
	printf ("Track_Point_Train	END\n");
}
#endif


//#if 0

void	Track_Point_DataInWrite	(tTrack_Point* p, si ox, si oy, fann_type* ain)
{
	ox -= p->W/2;	oy -= p->H/2;
	si x , y;
	for (y = 0; y < p->H; y++) {
		for (x = 0; x < p->W; x++) {
			ain[x + y * p->W] = dopix(ox+x,oy+y)->Y;
			ain[x + y * p->W] /= 255;
		//	ain[x + y * p->W] -= dopix(ox+x+1,oy+y)->Y;
		//	ain[x + y * p->W] /= 128;
		}
	}
}

tTrack_Point* gTrack_Point__p = 0;
void	Track_Point_cbData	(unsigned int num, unsigned int num_in, unsigned int num_out, fann_type* ain, fann_type* aout)
{
	tTrack_Point* p = gTrack_Point__p;
	
/*	si ix = num%3;
	si iy = num/3;
	
	ix -= 1;	iy -= 1;
	ix *= 5;	iy *= 5;
	
	Track_Point_DataInWrite (p, p->P.x + ix, p->P.y + iy, ain);
	
	aout[0] = -1;
	if (ix == 0 && iy == 0)
		aout[0] = 1;/**/
	
	if (num < 9) {
		si ix = num%3;
		si iy = num/3;
		ix -= 1;	iy -= 1;
		Track_Point_DataInWrite (p, p->P.x, p->P.y, ain);
		if (ix == 0 && iy == 0)
			aout[0] = 1;
		else
			aout[0] = 0.8;
	}else {
		si ox, oy;
		do {
			ox = gM.Cam.Image_W - p->W;
			oy = gM.Cam.Image_H - p->H;
			ox = rand() % ox;
			oy = rand() % oy;
			ox += p->W/2;
			oy += p->H/2;
		}while (ddist2(ox, oy, p->P.x,p->P.y) <= dpow2(3));
		
		printf ("Track_Point_cbData	oxy	%d	%d\n", ox, oy);
		
		Track_Point_DataInWrite (p, ox, oy, ain);
		
		aout[0] = -1;
	}/**/
}

void	Track_Point_Train		(tTrack_Point* p)
{
	gTrack_Point__p = p;
	
	const unsigned int num_input = p->W*p->H;
	const unsigned int num_output = 1;
	const unsigned int num_layers = 3;
	const unsigned int num_neurons_hidden = 100;
	const float desired_error = (const float) 0.001;
	const unsigned int max_epochs = 100;
	const unsigned int epochs_between_reports = 10;
	
	p->pANN = fann_create_standard (num_layers, num_input, num_neurons_hidden, num_output);
	
	fann_set_activation_function_hidden (p->pANN, FANN_SIGMOID_SYMMETRIC);
	fann_set_activation_function_output (p->pANN, FANN_SIGMOID_SYMMETRIC);
	
//	fann_train_on_file (ann, "xor.data", max_epochs, epochs_between_reports, desired_error);
	struct fann_train_data *data = fann_create_train_from_callback (
		9 + 32,
		num_input,
		num_output,
		Track_Point_cbData
	);
	
	fann_train_on_data (p->pANN, data, max_epochs, epochs_between_reports, desired_error);
	
	
//	fann_save(ann, "xor_float.net");
	
//	fann_destroy(p->pANN);
	
	return 0;
}

//#endif


void	Track_Point_Step		(tTrack_Point* p)
{
	if (!p->pANN)
		return;
	
	gTrack_Point__p = p;
	fann_type best_out = 0.0f, *calc, input[p->W*p->H];
	tV2si best_pos = {p->P.x, p->P.y};
	
	si x, y, dd = 4;
	for (y = -dd; y <= dd; y++) {
		for (x = -dd; x <= dd; x++) {
			Track_Point_DataInWrite (p, p->P.x+x, p->P.y+y, input);
			
			calc = fann_run (p->pANN, input);
			if (calc[0] > best_out) {
				best_out = calc[0];
				best_pos.x = p->P.x+x;
				best_pos.y = p->P.y+y;
			}
		}
	}
	
	printf ("Track_Point_Step best_out %f\n", best_out);
	if (best_out >= 0.1f) {
		p->P.x = dclip_lh (best_pos.x, p->W, gM.Cam.Image_W - p->W);
		p->P.y = dclip_lh (best_pos.y, p->H, gM.Cam.Image_H - p->H);
	}
	
}

void	Track_Point_Draw		(tTrack_Point* p)
{
	
	dset_c0(p->P.x,p->P.y);
	
	dset_c0(p->P.x - (p->W/2),		p->P.y - (p->H/2));
	dset_c0(p->P.x - (p->W/2)+1,		p->P.y - (p->H/2));
	dset_c0(p->P.x - (p->W/2)+2,		p->P.y - (p->H/2));
	dset_c0(p->P.x - (p->W/2),		p->P.y - (p->H/2)+1);
	dset_c0(p->P.x - (p->W/2),		p->P.y - (p->H/2)+2);
	
	dset_c0(p->P.x + (p->W/2),		p->P.y - (p->H/2));
	dset_c0(p->P.x + (p->W/2)-1,		p->P.y - (p->H/2));
	dset_c0(p->P.x + (p->W/2)-2,		p->P.y - (p->H/2));
	dset_c0(p->P.x + (p->W/2),		p->P.y - (p->H/2)+1);
	dset_c0(p->P.x + (p->W/2),		p->P.y - (p->H/2)+2);
	
	dset_c0(p->P.x - (p->W/2),		p->P.y + (p->H/2));
	dset_c0(p->P.x - (p->W/2)+1,		p->P.y + (p->H/2));
	dset_c0(p->P.x - (p->W/2)+2,		p->P.y + (p->H/2));
	dset_c0(p->P.x - (p->W/2),		p->P.y + (p->H/2)-1);
	dset_c0(p->P.x - (p->W/2),		p->P.y + (p->H/2)-2);
	
	dset_c0(p->P.x + (p->W/2),		p->P.y + (p->H/2));
	dset_c0(p->P.x + (p->W/2)-1,		p->P.y + (p->H/2));
	dset_c0(p->P.x + (p->W/2)-2,		p->P.y + (p->H/2));
	dset_c0(p->P.x + (p->W/2),		p->P.y + (p->H/2)-1);
	dset_c0(p->P.x + (p->W/2),		p->P.y + (p->H/2)-2);
	
}

#endif


void	Head_Point_Train	(tHead* p)
{
	si i = 0;
	//for (; i < dHead_Point_NUM; ++i)
		//Track_Point_Train (p->aPoint + i);
}

void	Head_Init		(tHead* p)
{
	si i = 0;
	//for (; i < dHead_Point_NUM; ++i)
		//Track_Point_Init (p->aPoint + i);
	
	Eye_Init (&p->DotC);
	Eye_Init (&p->DotL);
	Eye_Init (&p->DotR);
	
	M3f_Iden (&p->M);
	M3f_Iden (&p->MI);
	
	p->R_X = p->R_Y = p->R_Z = 0;
	
	M4f_Iden (&p->M4);	M4f_trans (&p->M4, 0, 0, -30);
	M4f_Iden (&p->M4_R);
	p->M4_T = p->M4;
	M4f_Inv (&p->M4, &p->M4I);
}

void	Head_Conf		(tHead* p)
{
	Eye_Conf (&p->DotC);
	Eye_Conf (&p->DotL);
	Eye_Conf (&p->DotR);
	
	si i = 0;
	//for (; i < dHead_Point_NUM; ++i)
		//Track_Point_Conf (p->aPoint + i);
}


void	HeadC_Draw		(tHead* p, tCam* pcam)
{
	
	#define dline3(x0,y0,z0,x1,y1,z1)	\
		do {	\
			tV4f __p0 = {x0, y0, z0, 1}, __p1 = {x1, y1, z1, 1};	\
			V4f_mul_M4f (&__p0, &p->M4);	\
			V4f_mul_M4f (&__p1, &p->M4);	\
			V4f_DrawPosPos (pcam, &__p0, &__p1);	\
		}while(0)
	
	#undef dpoint3
	#define dpoint3(name,x0,y0,z0)	\
		tV4f name = {x0,y0,z0,1}; M4f_mul_V4f (&p->M4, &name);
	
	gCol.R = 0x0;
	gCol.G = 0xFF;
	gCol.B = 0x0;
	
	si i;
	//for (i = 0; i < dHead_Point_NUM; ++i)
		//Track_Point_Draw (p->aPoint + i);
	
	{
		tV4f p0 = p->P; V4f_add_V4f (&p0, &p->N);
		V4f_DrawPosPos (pcam, &p->P, &p0);
	}
	
	if (0) {
		float d = 5;
		tV4f c = {0, 0, 0, 1};
		tV4f dx = {d, 0, 0, 1};
		tV4f dy = {0, d, 0, 1};
		tV4f dz = {0, 0, d, 1};
		
	//	V4f_mul_M4f (&c, &p->M4);
		M4f_mul_V4f (&p->M4, &c);
	//	printf (" c: ");	V4f_Print (&c);
	//	printf ("dx: ");	V4f_Print (&dx);
		
	//	V4f_mul_M4f (&dx, &p->M4);
	//	V4f_mul_M4f (&dy, &p->M4);
	//	V4f_mul_M4f (&dz, &p->M4);
		M4f_mul_V4f (&p->M4, &dx);
		M4f_mul_V4f (&p->M4, &dy);
		M4f_mul_V4f (&p->M4, &dz);
		
	//	printf ("closer\n");
		V4f_DrawPosPos (pcam, &c, &dx);
		V4f_DrawPosPos (pcam, &c, &dy);
		V4f_DrawPosPos (pcam, &c, &dz);
	}
	if (1) {
	/*	tV4f pc = {0,	-1,	3.5,	1};
		tV4f pl = {-7.5,	0,	0,	1};
		tV4f pr = {7.5,	0,	0,	1};/**/
	/*	tV4f pc = {0,	-1.5,	0,	1};
		tV4f pl = {-2.5,	0,	0,	1};
		tV4f pr = {2.5,	0,	0,	1};/**/
	/*	tV4f pc = {0,	2,	0,	1};
		tV4f pl = {-2.7,	0,	0,	1};
		tV4f pr = {3.5,	0,	0,	1};/**/
		tV4f pc = p->Mod.PC;
		tV4f pl = p->Mod.PL;
		tV4f pr = p->Mod.PR;/**/
	/*	printf ("pc: ");	V4f_Print (&pc);
		printf ("pl: ");	V4f_Print (&pl);
		printf ("pr: ");	V4f_Print (&pr);/**/
		
	/*	V4f_mul_M4f (&pc, &p->M4);
		V4f_mul_M4f (&pl, &p->M4);
		V4f_mul_M4f (&pr, &p->M4);/**/
		M4f_mul_V4f (&p->M4, &pc);
		M4f_mul_V4f (&p->M4, &pl);
		M4f_mul_V4f (&p->M4, &pr);
		
		V4f_DrawPosPos (pcam, &pl, &pc);
		V4f_DrawPosPos (pcam, &pr, &pc);
		V4f_DrawPosPos (pcam, &pl, &pr);
		if (0) {
			tV2f ret;
			V4f_ScreenPosNorm (pcam, &pl, &ret);	printf ("x %f y %f", ret.x, ret.y);
			V4f_ScreenPosNorm (pcam, &pr, &ret);	printf ("    x %f y %f\n", ret.x, ret.y);
		}
		if (0) {
			tV2si ret;
			V4f_ScreenPos (pcam, &pl, &ret);	printf ("x %ld y %ld", ret.x, ret.y);
			V4f_ScreenPos (pcam, &pr, &ret);	printf ("    x %ld y %ld\n", ret.x, ret.y);
		}
	}
	if (1) {
		float xdd = 10, ydd = 5;
	//	float xdd = 3.25, ydd = 2.25;
		tV4f p0 = {-xdd,	-ydd,	0,	1};
		tV4f p1 = {xdd,	-ydd,	0,	1};
		tV4f p2 = {xdd,	ydd,	0,	1};
		tV4f p3 = {-xdd,	ydd,	0,	1};
		
		M4f_mul_V4f (&p->M4, &p0);
		M4f_mul_V4f (&p->M4, &p1);
		M4f_mul_V4f (&p->M4, &p2);
		M4f_mul_V4f (&p->M4, &p3);
		
		V4f_DrawPosPos (pcam, &p0, &p1);
		V4f_DrawPosPos (pcam, &p1, &p2);
		V4f_DrawPosPos (pcam, &p2, &p3);
		V4f_DrawPosPos (pcam, &p3, &p0);
	}
	
	
	if (gM.bHead_Eye_LineDraw) {
		Head_Eye_LineDraw (p, &gM.Left);
		Head_Eye_LineDraw (p, &gM.Right);
	}
	
	if (0) {
		tV4f eye, ret, vec, re, rr, rv;
		Head_Eye_Vector (p, &gM.Left, &ret);
		Head_Eye_Vector (p, &gM.Right, &rr);
		
		eye = gM.Left.InHead.P;
		vec = ret;	V4f_sub_V4f (&vec, &gM.Left.InHead.P);
		
		re = gM.Right.InHead.P;
		rv = rr;	V4f_sub_V4f (&rv, &gM.Right.InHead.P);
		
	//	printf ("hoho2 ");	V4f_Print (&ret);
		
		M4f_mul_V4f (&p->M4, &eye);
		M4f_mul_V4f (&p->M4, &ret);
		
		M4f_mul_V4f (&p->M4, &re);
		M4f_mul_V4f (&p->M4, &rr);
		
	//	printf ("hoho1 ");	V4f_Print (&ret);
		
		gCol.R = 0x0;
		gCol.G = 0xFF;
		gCol.B = 0x0;
	///	V4f_DrawPosPos (pcam, &eye, &ret);
		
	//	printf ("vec ");	V4f_Print (&vec);
		
		ret = vec;
		
		float ax, ay;
		{
			M4f_mul_V4f (&p->M4_R, &vec);
			ax = atan2 (vec.z, vec.x);
			ay = atan2 (vec.z, vec.y);
			
			gM.L_Vec.x = ax;		gM.L_Vec.y = ay;
		//	printf ("axy  %f\t%f\n", ax, ay);
		}
		
		vec = ret;
		{
			M4f_mul_V4f (&p->M4I_R, &vec);
			ax = atan2 (vec.z, vec.x);
			ay = atan2 (vec.z, vec.y);
			
		//	gM.L_Vec.x = ax;		gM.L_Vec.y = ay;
		//	printf ("axyI %f\t%f\n", ax, ay);
		}
	//	gM.GazeL = Eye_map_point(&gM.Left, gM.L_Vec);
		
		
	}
	
	if (gM.GazeMode == 0) {
		Head_Dbg_Print (p);
	}else if (gM.GazeMode == 1) {
		Glint_Dbg_Print ();
	}else if (gM.GazeMode == 2) {
		Glint_Dbg_Print ();
	}
	
	
	if (0) {//hihihi couldn't resist
		dpoint3(pbl, -1.8,	3.7,		3.0);
		dpoint3(pbr, 1.8,		3.7,		3.0);
		
		dpoint3(pt, 0,		-1.0,		2.2);
		
		dpoint3(pc, 0,		3.0,		6.0);
		
		V4f_DrawPosPos (pcam, &pbl, &pbr);
		V4f_DrawPosPos (pcam, &pbl, &pt);
		V4f_DrawPosPos (pcam, &pbr, &pt);
		
		V4f_DrawPosPos (pcam, &pbl, &pc);
		V4f_DrawPosPos (pcam, &pbr, &pc);
		V4f_DrawPosPos (pcam, &pt, &pc);
	}
	
	if (1) {//ehhh gotta check this
		#define drect(www,hhh,ddd)	\
			do {	\
				float ww = www/2.0f;	\
				float hh = hhh/2.0f;	\
				float dd = ddd/2.0f;	\
				dpoint3(p00, pos.x-ww,	pos.y+hh,		pos.z-dd);	\
				dpoint3(p01, pos.x+ww,	pos.y+hh,		pos.z-dd);	\
				dpoint3(p02, pos.x+ww,	pos.y-hh,		pos.z-dd);	\
				dpoint3(p03, pos.x-ww,	pos.y-hh,		pos.z-dd);	\
					\
				dpoint3(p10, pos.x-ww,	pos.y+hh,		pos.z+dd);	\
				dpoint3(p11, pos.x+ww,	pos.y+hh,		pos.z+dd);	\
				dpoint3(p12, pos.x+ww,	pos.y-hh,		pos.z+dd);	\
				dpoint3(p13, pos.x-ww,	pos.y-hh,		pos.z+dd);	\
					\
				V4f_DrawPosPos (pcam, &p00, &p01);	\
				V4f_DrawPosPos (pcam, &p01, &p02);	\
				V4f_DrawPosPos (pcam, &p02, &p03);	\
				V4f_DrawPosPos (pcam, &p03, &p00);	\
					\
				V4f_DrawPosPos (pcam, &p10, &p11);	\
				V4f_DrawPosPos (pcam, &p11, &p12);	\
				V4f_DrawPosPos (pcam, &p12, &p13);	\
				V4f_DrawPosPos (pcam, &p13, &p10);	\
					\
				V4f_DrawPosPos (pcam, &p00, &p10);	\
				V4f_DrawPosPos (pcam, &p01, &p11);	\
				V4f_DrawPosPos (pcam, &p02, &p12);	\
				V4f_DrawPosPos (pcam, &p03, &p13);	\
			}while(0)
		
		float ehh_y = -1.6, ehh_z = 3.2;
		tV4f pos;
		pos.x = 0; pos.y = ehh_y+0.15;	pos.z = ehh_z-0.75;
		drect (16, 0.3f, 1.5f);
		
		pos.x = 0; pos.y = ehh_y + 0.3f + 0.9f/2.0f;	pos.z = ehh_z - 0.75f/2;
		drect (8, 0.9f, 0.75f);
		
		pos.x = 0; pos.y = ehh_y + 1.2f + 0.9f/2.0f;	pos.z = ehh_z - 0.75f/2;
		drect (1.6f, 0.9f, 0.75f);
		
		pos.x = -8; pos.y = ehh_y+0.3+0.15;	pos.z = ehh_z-1.5f;
		drect (1.5f, 0.3f, 4.8f);
		
		pos.x = 8;
		drect (1.5f, 0.3f, 4.8f);
		
	/*	tV4f pos;
		pos.x = 0; pos.y = p->Mod.PC.y+0.15;	pos.z = p->Mod.PC.z-0.75;
		drect (16, 0.3f, 1.5f);
		
		pos.x = 0; pos.y = p->Mod.PC.y + 0.3f + 0.9f/2.0f;	pos.z = p->Mod.PC.z - 0.75f/2;
		drect (8, 0.9f, 0.75f);
		
		pos.x = 0; pos.y = p->Mod.PC.y + 1.2f + 0.9f/2.0f;	pos.z = p->Mod.PC.z - 0.75f/2;
		drect (1.6f, 0.9f, 0.75f);
		
		pos.x = p->Mod.PL.x; pos.y = p->Mod.PC.y+0.3+0.15;	pos.z = p->Mod.PC.z-1.5f;
		drect (1.5f, 0.3f, 4.8f);
		
		pos.x = p->Mod.PR.x;
		drect (1.5f, 0.3f, 4.8f);/**/
		
		#undef drect
	}
}

#if 1
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

void	HeadC_Calc_M4_Rel	(tHead* p, tCam* pcam, tHead* pcen)
{
	if (0) {
	//	return;
		si i, j;
	//	int fLength = 700;/* Focal Length in PIXELS */
	//	int imageCenterX = 0;/* In this test, origin of image coordinates is at image center */
	//	int imageCenterY = 0;/* More often, imageCenterX=256, imageCenterY=256 */
		int isMaxRank;
		TObject obj;
		TImage img;
		TCamera testCamera;
		
		img.imageCenter[0] = 400;
		img.imageCenter[1] = 300;
		
		/* Init structures and read object and image data:*/	
	//	GetObjectAndImage("TestCube", "CubeImage", &obj, &testImage);
		
		obj.nbPts = 3;
		obj.objectPts = InitDoubleArray(obj.nbPts, 3);
		
		obj.objectPts[0][0] = 0;	obj.objectPts[0][1] = 0;	obj.objectPts[0][2] = 3.5;
		obj.objectPts[1][0] = -7.5;	obj.objectPts[1][1] = -1;	obj.objectPts[1][2] = 0;
		obj.objectPts[2][0] = 7.5;	obj.objectPts[2][1] = -1;	obj.objectPts[2][2] = 0;/**/
		
	/*	obj.objectPts[0][0] = 0;	obj.objectPts[0][1] = -1.5;	obj.objectPts[0][2] = 1;
		obj.objectPts[1][0] = -2.5;	obj.objectPts[1][1] = 0;	obj.objectPts[1][2] = 0;
		obj.objectPts[2][0] = 2.5;	obj.objectPts[2][1] = 0;	obj.objectPts[2][2] = 0;/**/
		
		obj.objectVects = InitDoubleArray(obj.nbPts, 3);
		obj.objectCopy = InitDoubleArray(obj.nbPts, 3);
			
		for (i=0;i<obj.nbPts;i++){
			for(j=0;j<3;j++){
				obj.objectVects[i][j] = obj.objectCopy[i][j] = obj.objectPts[i][j] - obj.objectPts[0][j];
			}
		}
		obj.objectMatrix = InitDoubleArray(3, obj.nbPts);
		
		/* Image data */
		
		img.nbPts = obj.nbPts;
		img.imagePts = (int**)InitDoubleArray(img.nbPts, 2);
	//	img.imagePts = ReadImage(imageName, img.imageCenter, img.nbPts);
		
		img.imagePts[0][0] = p->DotC.P.x;	img.imagePts[0][1] = p->DotC.P.y;
		img.imagePts[1][0] = p->DotL.P.x;	img.imagePts[1][1] = p->DotL.P.y;
		img.imagePts[2][0] = p->DotR.P.x;	img.imagePts[2][1] = p->DotR.P.y;
		
		img.imageVects = InitDoubleArray(img.nbPts, 2);
		img.oldImageVects = InitDoubleArray(img.nbPts, 2);
		img.epsilon = (double *)malloc(img.nbPts * sizeof(double));
		
		testCamera.focalLength = 800;
	//	testCamera.focalLength = 800/84;
	//	double rotation[3][3];/* Rotation of SCENE in camera reference, NOT other way around */
	//	double translation[3];/* Translation of SCENE in camera reference */
		testCamera.rotation[0][0] = 1;	testCamera.rotation[0][1] = 0;	testCamera.rotation[0][2] = 0;
		testCamera.rotation[1][0] = 0;	testCamera.rotation[1][1] = 1;	testCamera.rotation[1][2] = 0;
		testCamera.rotation[2][0] = 0;	testCamera.rotation[2][1] = 0;	testCamera.rotation[2][2] = 1;
		testCamera.translation[0] = 0;
		testCamera.translation[1] = 0;
		testCamera.translation[2] = 50;
		
		/* Find object matrix as pseudoInverse of matrix of object vector coordinates */
		isMaxRank = PseudoInverse(obj.objectCopy, obj.nbPts, obj.objectMatrix);
		/* Check if matrix rank for object matrix is less than 3: */
		if(!isMaxRank){
			nrerror("object is too flat; another method is required ");/* exit */
		}else {
			/* Find object pose by POSIT, Pose from Orthography, Scaling, and ITerations */
			POSIT(obj, img, &testCamera);
			
			/* Print results */
			printf("\n");
			printf("Rotation matrix of scene in camera reference frame\n");
			PrintRotation(testCamera.rotation);
			printf("\n");
			printf("Translation vector of scene in camera reference frame\n" 
				"from camera projecction center to FIRST point of scene (object) file\n");
			PrintVector(testCamera.translation, 3);
			printf("\n");
			
			tM4f rot;	M4f_Iden (&rot);
			tM4f trans;	M4f_Iden (&trans);
			
			rot.x00 = testCamera.rotation[0][0];	rot.x01 = testCamera.rotation[0][1];	rot.x02 = testCamera.rotation[0][2];
			rot.x10 = testCamera.rotation[1][0];	rot.x11 = testCamera.rotation[1][1];	rot.x12 = testCamera.rotation[1][2];
			rot.x20 = testCamera.rotation[2][0];	rot.x21 = testCamera.rotation[2][1];	rot.x22 = testCamera.rotation[2][2];/**/
		/*	rot.x00 = testCamera.rotation[0][0];	rot.x01 = testCamera.rotation[1][0];	rot.x02 = testCamera.rotation[2][0];
			rot.x10 = testCamera.rotation[0][1];	rot.x11 = testCamera.rotation[1][1];	rot.x12 = testCamera.rotation[2][1];
			rot.x20 = testCamera.rotation[0][2];	rot.x21 = testCamera.rotation[1][2];	rot.x22 = testCamera.rotation[2][2];/**/
			
			trans.x30 = testCamera.translation[0];
			trans.x31 = testCamera.translation[1];
			trans.x32 = testCamera.translation[2];
			M4f_Iden (&p->M4);
			M4f_mul_M4f (&p->M4, &rot);
			M4f_mul_M4f (&p->M4, &trans);
		}
		return;
	}
	if (1) {
		tV2f osl = p->DotL.aCam[pcam->Idx].P;
		tV2f osr = p->DotR.aCam[pcam->Idx].P;
		tV2f osc = p->DotC.aCam[pcam->Idx].P;
		
		tV2f os_vrl = osl;	V2f_sub_V2f (&os_vrl, &osr);
		tV2f os_vlc = osc;	V2f_sub_V2f (&os_vlc, &osl);
		tV2f os_vrc = osc;	V2f_sub_V2f (&os_vrc, &osr);
		
		float osdlr = ddist (osl.x,osl.y,osr.x,osr.y);
		float osdlc = ddist (osl.x,osl.y,osc.x,osc.y);
		float osdrc = ddist (osr.x,osr.y,osc.x,osc.y);
		
		
		static tM4f rot = {1,0,0,0, 0,1,0,0, 0,0,1,0, 0,0,0,1};
		static tM4f trans = {1,0,0,0, 0,1,0,0, 0,0,1,0, 0,0,0,1};
		
		if (1) {
			p->dcam.P.x = p->dcam.P.y = 0;
			p->dcam.P.z = -30;
			p->R_X = p->R_Y = p->R_Z = 0;
			M4f_Iden (&rot);	M4f_Iden (&trans);	M4f_trans (&trans, 0, 0, -30);
		}
	/*	if (!finite(p->dcam.M4.x03)) {
			p->dcam.P.x = p->dcam.P.y = 0;
			p->dcam.P.z = -30;
			p->R_X = p->R_Y = p->R_Z = 0;
		}/**/
		
		float err = 0;
		static si i = 0;
		for (i = 0; i < 100; ++i) {
			if (1) {
			//	M4f_Iden (&trans);
			//	M4f_trans (&trans, p->dcam.P.x, p->dcam.P.y, p->dcam.P.z);
				
				M4f_Iden (&rot);
				M4f_rotx (&rot, p->R_X);
				M4f_roty (&rot, p->R_Y);
				M4f_rotz (&rot, p->R_Z);
			}
			M4f_Iden (&p->dcam.M4);
			M4f_mul_M4f (&p->dcam.M4, &rot);
			M4f_mul_M4f (&p->dcam.M4, &trans);
			
		/*	tV4f pc = {0,	-1,	3.5,	1};
			tV4f pl = {-7.5,	0,	0,	1};
			tV4f pr = {7.5,	0,	0,	1};/**/
		/*	tV4f pc = {0,	-1.5,	0,	1};
			tV4f pl = {-2.5,	0,	0,	1};
			tV4f pr = {2.5,	0,	0,	1};/**/
		/*	tV4f pc = {0,	2,	0,	1};
			tV4f pl = {-2.7,	0,	0,	1};
			tV4f pr = {3.5,	0,	0,	1};/**/
			tV4f pc = p->Mod.PC;
			tV4f pl = p->Mod.PL;
			tV4f pr = p->Mod.PR;/**/
			
		/*	V4f_mul_M4f (&pc, &p->dcam.M4);
			V4f_mul_M4f (&pl, &p->dcam.M4);
			V4f_mul_M4f (&pr, &p->dcam.M4);/**/
			
			M4f_mul_V4f (&p->dcam.M4, &pc);
			M4f_mul_V4f (&p->dcam.M4, &pl);
			M4f_mul_V4f (&p->dcam.M4, &pr);
			
			tV2f sc, sl, sr;
			V4f_ScreenPosf (pcam, &pl, &sl);
			V4f_ScreenPosf (pcam, &pr, &sr);
			V4f_ScreenPosf (pcam, &pc, &sc);
			
			tV2f s_vrl = sl;	V2f_sub_V2f (&s_vrl, &sr);
			tV2f s_vlc = sc;	V2f_sub_V2f (&s_vlc, &sl);
			tV2f s_vrc = sc;	V2f_sub_V2f (&s_vrc, &sr);
			
			float sdlr = ddist (sl.x,sl.y,sr.x,sr.y);
			float sdlc = ddist (sl.x,sl.y,sc.x,sc.y);
			float sdrc = ddist (sr.x,sr.y,sc.x,sc.y);
			
			err = 0;
			err += sc.x - osc.x;			err += sc.y - osc.y;
			err += sl.x - osl.x;			err += sl.y - osl.y;
			err += sr.x - osr.x;			err += sr.y - osr.y;
			
			float dx, dy, dz;
			dx = osc.x - sc.x;	dx += osl.x - sl.x;	dx += osr.x - sr.x;	dx /= 3.0f;
			dy = osc.y - sc.y;	dy += osl.y - sl.y;	dy += osr.y - sr.y;	dy /= 3.0f;
			dz = osdlr - sdlr;	dz += osdlc - sdlc;	dz += osdrc - sdrc;	dz /= 3.0f;
			
		/*	p->dcam.P.x += p->Mod.TInc.x*dx;
			p->dcam.P.y += p->Mod.TInc.y*dy;
			p->dcam.P.z += p->Mod.TInc.z*dz;
			p->dcam.P.w = 1;/**/
			M4f_trans (&trans, p->Mod.TInc.x*dx, p->Mod.TInc.y*dy, p->Mod.TInc.z*dz);
			
		/*	printf ("cross os lc rc %f   rl rc %f   rl lc %f\n",
				V2f_cross (&os_vlc, &os_vrc) / (osdlc*osdrc),
				V2f_cross (&os_vrl, &os_vrc) / (osdlr*osdrc),
				V2f_cross (&os_vrl, &os_vlc) / (osdlr*osdlc)
			);/**/
			
			if (1) {
			/*	printf ("atan2  rl %f   lc %f   rc %f\n",
					atan2(os_vrl.y,os_vrl.x),
					atan2(os_vlc.y,os_vlc.x),
					atan2(os_vrc.y,os_vrc.x)
				);/**/
				dz = angle_norm_pi_pi(atan2(s_vrl.y,s_vrl.x) - atan2(os_vrl.y,os_vrl.x));
				dz += angle_norm_pi_pi(atan2(s_vlc.y,s_vlc.x) - atan2(os_vlc.y,os_vlc.x));
				dz += angle_norm_pi_pi(atan2(s_vrc.y,s_vrc.x) - atan2(os_vrc.y,os_vrc.x));
			//	dz = -dz;
			//	printf ("rotz %f\n", dz);
				p->R_Z += p->Mod.RInc.z*dz;
			//	M4f_rotz (&rot, p->Mod.RInc.z*dz);
			}/**/
			
			if (1) {
				float ocross = V2f_cross (&os_vlc, &os_vrc) / (osdlc*osdrc);
				float cross = V2f_cross (&s_vlc, &s_vrc) / (sdlc*sdrc);
				dx = -ocross + cross;
			//	printf ("rotx %f\n", dx);
				p->R_X += -p->Mod.RInc.x*dx;
			//	M4f_rotx (&rot, -p->Mod.RInc.x*dx);
			}/**/
			if (1) {
				float ocross = V2f_cross (&os_vrl, &os_vlc) / (osdlr*osdlc) - V2f_cross (&os_vrl, &os_vrc) / (osdlr*osdrc);
				float cross = V2f_cross (&s_vrl, &s_vlc) / (sdlr*sdlc) - V2f_cross (&s_vrl, &s_vrc) / (sdlr*sdrc);
			//	float ocross = V2f_cross (&os_vrl, &os_vlc) - V2f_cross (&os_vrl, &os_vrc);
			//	float cross = V2f_cross (&s_vrl, &s_vlc) - V2f_cross (&s_vrl, &s_vrc);
				
				if (V2f_cross (&os_vrc, &os_vlc) < 0)
					ocross = -ocross;
				if (V2f_cross (&s_vrc, &s_vlc) < 0)
					cross = -cross;
				
			//	printf ("roty ocross %f\n", ocross);
				if (fabsf(ocross) >= 0.05f) {
					dy = ocross - cross;
				}else {
					dy = (sdlc-sdrc)/sdlr - (osdlc-osdrc)/osdlr;
				//	dy *= 0.1f;
			//		printf ("roty %f\n", dy);
				}
				p->R_Y += -p->Mod.RInc.y*dy;
			//	printf ("roty %f\n", dy);
			//	M4f_roty (&rot, -p->Mod.RInc.y*dy);
			}/**/
			
		/*	if (++i >= 100) {
				i = 0;
				M4f_Iden (&rot);
				M4f_Iden (&trans);
				M4f_trans (&trans, 0, 0, -80);
			}/**/
		}
		if (1) {
		//	M4f_Iden (&trans);
		//	M4f_trans (&trans, p->dcam.P.x, p->dcam.P.y, p->dcam.P.z);
			
			M4f_Iden (&rot);
			M4f_rotx (&rot, p->R_X);
			M4f_roty (&rot, p->R_Y);
			M4f_rotz (&rot, p->R_Z);
		}
	/*	if (fabsf(err) > 5 || !finite(trans.x03)) {
			M4f_Iden (&rot);	M4f_Iden (&trans);	M4f_trans (&trans, 0, 0, -40);
			p->dcam.P.x = p->dcam.P.y = 0;
			p->dcam.P.z = -40;
			p->R_X = p->R_Y = p->R_Z = 0;
		}/**/
		p->dcam.M4_T = trans;
		
		if (1) {
			M4f_Iden (&rot);
			M4f_rotx (&rot, p->R_X);
			M4f_roty (&rot, p->R_Y);
			M4f_rotz (&rot, p->R_Z);
		}
		if (0) {
			M4f_rotx (&rot, p->SRX*p->R_X);
			M4f_roty (&rot, p->SRY*p->R_Y);
		}
		p->dcam.M4_R = rot;
		
		M4f_Iden (&p->dcam.M4);
		M4f_mul_M4f (&p->dcam.M4, &p->dcam.M4_R);
		M4f_mul_M4f (&p->dcam.M4, &p->dcam.M4_T);
		
		p->dcam.P.x = p->dcam.M4.x03;
		p->dcam.P.y = p->dcam.M4.x13;
		p->dcam.P.z = p->dcam.M4.x23;
		p->dcam.P.w = 1;
		
		p->dcam.N.x = 0;
		p->dcam.N.y = 0;
		p->dcam.N.z = 1;
		p->dcam.N.w = 1;
		M4f_mul_V4f (&p->dcam.M4_R, &p->N);
		
	//	printf ("p->P: ");	V4f_Print (&p->P);
	//	printf ("p->N: ");	V4f_Print (&p->N);
	//	V4f_sub_V4f (&p->P, &p->N);
	//	V4f_sub_V4f (&p->P, &p->N);
	//	V4f_sub_V4f (&p->P, &p->N);
	//	V4f_sub_V4f (&p->P, &p->N);
		
		M4f_Inv (&p->dcam.M4, &p->dcam.M4I);
		M4f_Inv (&p->dcam.M4_R, &p->dcam.M4I_R);
		M4f_Inv (&p->dcam.M4_T, &p->dcam.M4I_T);
		
	}
	printf ("c:\t%f\t%f\t%f\n", p->dcam.M4.x03, p->dcam.M4.x13, p->dcam.M4.x23);
	printf ("c:\t%f\t%f\t%f\n", p->dcam.M4_T.x03, p->dcam.M4_T.x13, p->dcam.M4_T.x23);
	
}


tV2f	Head_Diff_Eye	(tHead* p, tEye* peye)
{
	tV2f ret = peye->P;
	
//	V2f_sub_V2f (&ret, &gM.Head.DotC.P);	return ret;
	
	if (1) {
		ret.x -= p->M.x02;
		ret.y -= p->M.x12;
		return ret;
	}
	
	if (peye->P.x > p->DotC.P.x)
		V2f_sub_V2f (&ret, &p->DotR.P);
	else
		V2f_sub_V2f (&ret, &p->DotL.P);
	return ret;
}

tV4f	Head_Diff4_Eye4	(tHead* p, tEye* peye)
{
	tV4f p0 = {0, 0, 0, 0};
	tV4f p1 = {0, 0, -1, 0};
	
/*	if (gM.Eye_Line_Ray) {
		float ax = M_PI_2 - gM.Cam.Image_FOV_W/2 + (peye->P.x/gM.Cam.Image_W)*gM.Cam.Image_FOV_W;
		float ay = M_PI_2 - gM.Cam.Image_FOV_H/2 + (peye->P.y/gM.Cam.Image_H)*gM.Cam.Image_FOV_H;
		
		ax = M_PI_2 - ax;
		ay = M_PI_2 - ay;
		
	//	printf ("Image_FOV   %f\n", gM.Cam.Image_FOV);
	//	printf ("Image_FOV_W %f\n", gM.Cam.Image_FOV_W*rad2deg);
	//	printf ("axy %f %f\n", ax*rad2deg, ay*rad2deg);
		V4f_roty (&p1, ax);
		V4f_rotx (&p1, ay);
	}else {
		p1.x = pcam->Proj_L + (peye->P.x/gM.Cam.Image_W)*pcam->Proj_W;
		p1.y = pcam->Proj_B + (peye->P.y/gM.Cam.Image_H)*pcam->Proj_H;
	}*/
	Cam_Retina_Ray (peye, &p0, &p1, 0);
//	printf ("p1 %f %f %f\n", p1.x, p1.y, p1.z);
	
	
	tV4f n = p->N, negn = V4f_rneg (&p->N);
	tV4f negw = p->P;
	tV4f u = p1; //V4f_sub_V4f (&u, &p0);
//	printf ("n*u %f\n", V4f_dot_V4f (&n, &u));
	float s = V4f_dot_V4f (&n, &negw) / V4f_dot_V4f (&n, &u);
	
	u.x *= s;
	u.y *= s;
	u.z *= s;
	return u;
}

void	Head_Eye_Line	(tHead* p, tEye* peye, tV4f* ppos, tV4f* pvec)	//line from camera to retina in head_space
{
	tV4f p0 = {0, 0, 0, 1};
	tV4f p1 = {0, 0, -1, 1};
	/*
	if (gM.Eye_Line_Ray) {
		float ax = M_PI_2 - gM.Cam.Image_FOV_W/2 + (peye->P.x/gM.Cam.Image_W)*gM.Cam.Image_FOV_W;
		float ay = M_PI_2 - gM.Cam.Image_FOV_H/2 + (peye->P.y/gM.Cam.Image_H)*gM.Cam.Image_FOV_H;
		
		ax = M_PI_2 - ax;
		ay = M_PI_2 - ay;
		
	//	printf ("Image_FOV   %f\n", gM.Cam.Image_FOV);
	//	printf ("Image_FOV_W %f\n", gM.Cam.Image_FOV_W*rad2deg);
	//	printf ("axy %f %f\n", ax*rad2deg, ay*rad2deg);
		V4f_roty (&p1, ax);
		V4f_rotx (&p1, ay);
	}else {
		p1.x = pcam->Proj_L + (peye->P.x/gM.Cam.Image_W)*pcam->Proj_W;
		p1.y = pcam->Proj_B + (peye->P.y/gM.Cam.Image_H)*pcam->Proj_H;
	}/**/
	Cam_Retina_Ray (peye, &p0, &p1, 0);
//	printf ("p0: ");	V4f_Print (&p0);
//	printf ("p1: ");	V4f_Print (&p1);
	
	tV4f u;
	V4f_Intersect_Line01_Plane0N (
		&u,
		&p0, &p1,
		&p->P, &p->N
	);
	
	/*
	tV4f n = p->N, negn = V4f_rneg (&p->N);
	tV4f negw = p->P;
	tV4f u = p1; //V4f_sub_V4f (&u, &p0);
//	printf ("n*u %f\n", V4f_dot_V4f (&n, &u));
	float s = V4f_dot_V4f (&n, &negw) / V4f_dot_V4f (&n, &u);
	
	u.x *= s;
	u.y *= s;
	u.z *= s;
	/**/
	*ppos = u;		V4f_sub_V4f (ppos, &p->P);
	
	*pvec = u;
	pvec->x *= 0.2f;	pvec->y *= 0.2f;	pvec->z *= 0.2f;
//	V4f_add_V4f (ppos, pvec);
	V4f_lneg (pvec);
//	pvec->x *= 2.0f;	pvec->y *= 2.0f;	pvec->z *= 2.0f;
	
//	printf ("onplane0 %f %f %f\n", onplane.x, onplane.y, onplane.z);
	M4f_mul_V4f (&gM.Head.M4I_R, ppos);
	M4f_mul_V4f (&gM.Head.M4I_R, pvec);
//	printf ("onplane1 %f %f %f\n", onplane.x, onplane.y, onplane.z);
	
}

void	Head_Eye_LineAdd	(tHead* p, tEye* peye)
{
	tV4f pos, vec;
	Head_Eye_Line (p, peye, &pos, &vec);
	
//	printf ("pos: ");	V4f_Print (&pos);
//	printf ("vec: ");	V4f_Print (&vec);
	
	peye->InHead.aLine[peye->InHead.Line_N].P = pos;
	peye->InHead.aLine[peye->InHead.Line_N].V = vec;
	
	peye->InHead.aLine[peye->InHead.Line_N].P0 = pos;
	V4f_sub_V4f (&peye->InHead.aLine[peye->InHead.Line_N].P0, &vec);
	peye->InHead.aLine[peye->InHead.Line_N].P1 = pos;
	V4f_add_V4f (&peye->InHead.aLine[peye->InHead.Line_N].P1, &vec);
	peye->InHead.Line_N++;
	
	if (peye->InHead.Line_N >= 3)
		Head_Eye_CalcP (p, peye);
}

void	Head_Eye_LineDraw	(tHead* p, tEye* peye)
{
	if (1) {
		gCol.R = 0x0;
		gCol.G = 0xFF;
		gCol.B = 0x0;
		si i;
		for (i = 0; i < peye->InHead.Line_N; ++i) {
			tV4f p0 = peye->InHead.aLine[i].P0, p1 = peye->InHead.aLine[i].P1;
			
		//	p1 = p0;
		//	V4f_add_V4f (&p1, &p->N);
		//	printf ("p0 %f %f %f\n", p0.x, p0.y, p0.z);
		//	printf ("p1 %f %f %f\n", p1.x, p1.y, p1.z);
			
			M4f_mul_V4f (&p->M4, &p0);
			M4f_mul_V4f (&p->M4, &p1);
			
		//	printf ("p0 %f %f %f\n", p0.x, p0.y, p0.z);
		//	printf ("p1 %f %f %f\n", p1.x, p1.y, p1.z);
			
			V4f_DrawPosPos (pcam, &p0, &p1);
			
			Dbg_V4f_ADrawPosPos (&p0, &p1);
		}
	}
	gCol.R = 0xFF;
	gCol.G = 0xFF;
	gCol.B = 0x0;
	tV4f pos = peye->InHead.P;
	M4f_mul_V4f (&p->M4, &pos);
	V4f_DrawPosPos (pcam, &pos, &pos);
}

void	Head_Eye_Vector	(tHead* p, tEye* peye, tV4f* pret)	//Point of retina in head_space
{
	float r = peye->InHead.R;
	tV4f p0, vec, p1;
	Head_Eye_Line (p, peye, &p0, &vec);
	
	p1 = p0;	V4f_add_V4f (&p1, &vec);
	
	float a = V4f_dist2 (&vec);
	float b = 2*( (p1.x-p0.x)*(p0.x-peye->InHead.P.x)
				+ (p1.y-p0.y)*(p0.y-peye->InHead.P.y)
				+ (p1.z-p0.z)*(p0.z-peye->InHead.P.z)
	);
	float c = V4f_dist2 (&peye->InHead.P) + V4f_dist2 (&p0) - 2*(peye->InHead.P.x*p0.x + peye->InHead.P.y*p0.y + peye->InHead.P.z*p0.z) - dpow2(r);
	
	float d = b*b - 4*a*c;
	if (d < 0) {
		pret->x = 0;
		pret->y = 0;
		pret->z = 0;
		pret->w = 1;
		
		tV4f pos = peye->InHead.P;
		M4f_mul_V4f (&gM.Head.M4, &pos);
		
		tV2si screen_pos;
		abort();//V4f_ScreenPos (&pos, &screen_pos);
		
		peye->P.x = screen_pos.x;
		peye->P.y = screen_pos.y;
	//	printf ("eye failed"
		return;
	}
//	printf ("d %f\n", d);
	float u;
	
	u = (-b - sqrt(d)) / (2*a);
	pret->x = p0.x + u*vec.x;
	pret->y = p0.y + u*vec.y;
	pret->z = p0.z + u*vec.z;
	pret->w = 1;
//	printf ("pret0 "); V4f_Print (pret);
	
	u = (-b + sqrt(d)) / (2*a);
	pret->x = p0.x + u*vec.x;
	pret->y = p0.y + u*vec.y;
	pret->z = p0.z + u*vec.z;
	pret->w = 1;
//	printf ("pret1 "); V4f_Print (pret);
	
/*	u = (-b - sqrt(d)) / (2*a);
	pret->x = p0.x + u*vec.x;
	pret->y = p0.y + u*vec.y;
	pret->z = p0.z + u*vec.z;
	pret->w = 1;
	printf ("pret0 "); V4f_Print (pret);/**/
}

void	Head_Eye_VectorGlob	(tHead* p, tEye* peye, tV4f* pr0, tV4f* pr1, tV4f* prv)	//gaze vector of retina in world_space
{
	if (gM.GazeMode == 2) {
		Eye_GP_GazeVector (peye, pr0, pr1, prv);
		return;
	}
	tV4f eye, ret, vec;
	Head_Eye_Vector (p, peye, &ret);
	
	eye = peye->InHead.P;
	
	M4f_mul_V4f (&p->M4, &eye);
	M4f_mul_V4f (&p->M4, &ret);
	
	if (pr0)
		*pr0 = eye;
	if (pr1)
		*pr1 = ret;
	if (prv) {
		*prv = ret;	V4f_sub_V4f (prv, &eye);
	}
}

#endif

void PrintMat(CvMat *A)
{
	int i, j;
	for (i = 0; i < A->rows; i++)
	{
	printf("\n"); 
	switch (CV_MAT_DEPTH(A->type))
	{
	case CV_32F:
	case CV_64F:
	for (j = 0; j < A->cols; j++)
	printf ("%8.3f ", (float)cvGetReal2D(A, i, j));
	break;
	case CV_8U:
	case CV_16U:
	for(j = 0; j < A->cols; j++)
	printf ("%6d",(int)cvGetReal2D(A, i, j));
	break;
	default:
	break;
	}
	}
	printf("\n");
}



float	Head_Eye_CalcP_PosEr	(tHead* p, tEye* peye, tV4f* pos)
{
	float mse = 0;
	si i;
	for (i = 0; i < peye->InHead.Line_N; ++i) {
		tV4f p0 = peye->InHead.aLine[i].P0;	V4f_sub_V4f (&p0, pos);
		tV4f p1 = peye->InHead.aLine[i].P1;	V4f_sub_V4f (&p1, pos);
		
		tV4f cross;
		V4f_cross (&cross, &p1, &p0);
		mse += V4f_dot_V4f (&cross, &cross) / V4f_dot_V4f(&p1, &p1);
	}
	return mse / (float)peye->InHead.Line_N;
}

void	Head_Eye_CalcP_Iter	(tHead* p, tEye* peye)
{
	printf ("Head_Eye_CalcP_Iter error %f\n", Head_Eye_CalcP_PosEr (p, peye, &peye->InHead.P));
	
	float dx = 0.1f, dy = dx, dz = dx;
	si i;
	for (i = 0; i < 100; ++i) {
		si fix = 0;
		
		#define dfix(dim,val)	\
			do {			\
				tV4f tmp = peye->InHead.P;			\
				float omse = Head_Eye_CalcP_PosEr (p, peye, &tmp);			\
				tmp.dim += val;			\
				float nmse = Head_Eye_CalcP_PosEr (p, peye, &tmp);			\
				if (nmse < omse) {			\
					peye->InHead.P = tmp;			\
					++fix;			\
				}else {			\
					tmp.dim -= 2*val;			\
					nmse = Head_Eye_CalcP_PosEr (p, peye, &tmp);			\
					if (nmse < omse) {			\
						peye->InHead.P = tmp;			\
						++fix;			\
					}			\
				}			\
			}while(0)
		
		dfix (x,dx);
		dfix (y,dy);
		dfix (z,dz);
		
		if (fix == 0) {
			dx /= 2;
			dy /= 2;
			dz /= 2;
		}
		
		printf ("Head_Eye_CalcP_Iter error %f\n", Head_Eye_CalcP_PosEr (p, peye, &peye->InHead.P));
	}
	
}

void	Head_Eye_CalcP	(tHead* p, tEye* peye)
{
	CvMat* ma = cvCreateMatHeader(peye->InHead.Line_N, 3, CV_32FC1);
	cvCreateData(ma);
	CvMat* mb = cvCreateMatHeader(peye->InHead.Line_N, 1, CV_32FC1);
	cvCreateData(mb);
	si i;
	for (i = 0; i < peye->InHead.Line_N; ++i) {
		float a, b, c, d;
		a = atan2(peye->InHead.aLine[i].V.x, peye->InHead.aLine[i].V.z);
		b = atan2(peye->InHead.aLine[i].V.y, peye->InHead.aLine[i].V.z);
		c = 1;
		
		cvSet2D (ma, i, 0, cvScalarAll(a*10));
		cvSet2D (ma, i, 1, cvScalarAll(b*10));
		cvSet2D (ma, i, 2, cvScalarAll(c*10));
		
		d = a*peye->InHead.aLine[i].P.x + b*peye->InHead.aLine[i].P.y;
		cvSet2D (mb, i, 0, cvScalarAll(d*10));
	}
	
	CvMat* mx = cvCreateMatHeader(3, 1, CV_32FC1);
	cvCreateData(mx);
	cvSet (mx, cvScalarAll(0), 0);
	
	printf ("got %d\n", cvSolve(ma, mb, mx, CV_SVD));
	
	PrintMat (ma);
	PrintMat (mb);
	PrintMat (mx);
	
	peye->InHead.P.x = cvGet2D (mx, 0, 0).val[0];
	peye->InHead.P.y = cvGet2D (mx, 1, 0).val[0];
	peye->InHead.P.z = cvGet2D (mx, 2, 0).val[0];
	peye->InHead.P.w = 1;
	
	Head_Eye_CalcP_Iter (p, peye);
	
	printf ("xyz %f %f %f\n", peye->InHead.P.x, peye->InHead.P.y, peye->InHead.P.z);
}

#define dcal(_ix,_iy)	(peye->aScreen[p->Idx].aaCal[_iy][_ix])

void	Screen_Eye_PreCal	(tScreen* p, tEye* peye)
{
	peye->aScreen[p->Idx].aaCal[0][0].SX = 1280;
	peye->aScreen[p->Idx].aaCal[0][0].SY = 0;
	peye->aScreen[p->Idx].aaCal[0][0].P = (tV4f){gM.aScreen[p->Idx].C.x - gM.aScreen[p->Idx].W/2, gM.aScreen[p->Idx].C.y - gM.aScreen[p->Idx].H/2, gM.aScreen[p->Idx].C.z, 1};
	
	peye->aScreen[p->Idx].aaCal[0][1].SX = 0;
	peye->aScreen[p->Idx].aaCal[0][1].SY = 0;
	peye->aScreen[p->Idx].aaCal[0][1].P = (tV4f){gM.aScreen[p->Idx].C.x + gM.aScreen[p->Idx].W/2, gM.aScreen[p->Idx].C.y - gM.aScreen[p->Idx].H/2, gM.aScreen[p->Idx].C.z, 1};
	
	peye->aScreen[p->Idx].aaCal[1][0].SX = 1280;
	peye->aScreen[p->Idx].aaCal[1][0].SY = 1024;
	peye->aScreen[p->Idx].aaCal[1][0].P = (tV4f){gM.aScreen[p->Idx].C.x - gM.aScreen[p->Idx].W/2, gM.aScreen[p->Idx].C.y + gM.aScreen[p->Idx].H/2, gM.aScreen[p->Idx].C.z, 1};
	
	peye->aScreen[p->Idx].aaCal[1][1].SX = 0;
	peye->aScreen[p->Idx].aaCal[1][1].SY = 1024;
	peye->aScreen[p->Idx].aaCal[1][1].P = (tV4f){gM.aScreen[p->Idx].C.x + gM.aScreen[p->Idx].W/2, gM.aScreen[p->Idx].C.y + gM.aScreen[p->Idx].H/2, gM.aScreen[p->Idx].C.z, 1};
	
	p->Cal.ix = 0;	p->Cal.iy = 0;	printf ("%d %d P: ", p->Cal.ix, p->Cal.iy);	V4f_Print (&peye->aScreen[p->Idx].aaCal[p->Cal.iy][p->Cal.ix].P);
	p->Cal.ix = 1;	p->Cal.iy = 0;	printf ("%d %d P: ", p->Cal.ix, p->Cal.iy);	V4f_Print (&peye->aScreen[p->Idx].aaCal[p->Cal.iy][p->Cal.ix].P);
	p->Cal.ix = 0;	p->Cal.iy = 1;	printf ("%d %d P: ", p->Cal.ix, p->Cal.iy);	V4f_Print (&peye->aScreen[p->Idx].aaCal[p->Cal.iy][p->Cal.ix].P);
	p->Cal.ix = 1;	p->Cal.iy = 1;	printf ("%d %d P: ", p->Cal.ix, p->Cal.iy);	V4f_Print (&peye->aScreen[p->Idx].aaCal[p->Cal.iy][p->Cal.ix].P);
}



void	Screen_Cal_Eye_ExtraXY	(tScreen* p, tEye* peye, si ix, si iy)
{
	si tmp, ix0 = -1, ix1 = -1;
	for (tmp = ix-1; tmp >= 0; --tmp) {
		if (ix1 < 0 && dcal(tmp,iy).State == dEye_Screen_Cal_Set) {
			ix1 = tmp;
			continue;
		}
		if (dcal(tmp,iy).State == dEye_Screen_Cal_Set) {
			ix0 = tmp;
			break;
		}
	}
	if (ix0 < 0 || ix1 < 0) {
		ix0 = -1; ix1 = -1;
		for (tmp = ix+1; tmp < dEye_Screen_Cal_NUM; ++tmp) {
			if (ix0 < 0 && dcal(tmp,iy).State == dEye_Screen_Cal_Set) {
				ix0 = tmp;
				continue;
			}
			if (dcal(tmp,iy).State == dEye_Screen_Cal_Set) {
				ix1 = tmp;
				break;
			}
		}
	}
	dcal(ix,iy).P.x = 0;
	dcal(ix,iy).P.y = 0;
	dcal(ix,iy).P.z = 0;
	dcal(ix,iy).P.w = 0;
	
	if (ix0 >= 0 && ix1 >= 0) {
		dcal(ix,iy).P.x += dcal(ix0,iy).P.x + (ix-ix0) * (dcal(ix1,iy).P.x-dcal(ix0,iy).P.x) / (ix1-ix0);
		dcal(ix,iy).P.y += dcal(ix0,iy).P.y + (ix-ix0) * (dcal(ix1,iy).P.y-dcal(ix0,iy).P.y) / (ix1-ix0);
		dcal(ix,iy).P.z += dcal(ix0,iy).P.z + (ix-ix0) * (dcal(ix1,iy).P.z-dcal(ix0,iy).P.z) / (ix1-ix0);
		dcal(ix,iy).P.w += 1;
	}
	si iy0 = -1, iy1 = -1;
	for (tmp = iy-1; tmp >= 0; --tmp) {
		if (iy1 < 0 && dcal(ix,tmp).State == dEye_Screen_Cal_Set) {
			iy1 = tmp;
			continue;
		}
		if (dcal(ix,tmp).State == dEye_Screen_Cal_Set) {
			iy0 = tmp;
			break;
		}
	}
	if (iy0 < 0 || iy1 < 0) {
		iy0 = -1; iy1 = -1;
		for (tmp = iy+1; tmp < dEye_Screen_Cal_NUM; ++tmp) {
			if (iy0 < 0 && dcal(ix,tmp).State == dEye_Screen_Cal_Set) {
				iy0 = tmp;
				continue;
			}
			if (dcal(ix,tmp).State == dEye_Screen_Cal_Set) {
				iy1 = tmp;
				break;
			}
		}
	}
	if (iy0 >= 0 && iy1 >= 0) {
		dcal(ix,iy).P.x += dcal(ix,iy0).P.x + (iy-iy0) * (dcal(ix,iy1).P.x - dcal(ix,iy0).P.x) / (iy1-iy0);
		dcal(ix,iy).P.y += dcal(ix,iy0).P.y + (iy-iy0) * (dcal(ix,iy1).P.y - dcal(ix,iy0).P.y) / (iy1-iy0);
		dcal(ix,iy).P.z += dcal(ix,iy0).P.z + (iy-iy0) * (dcal(ix,iy1).P.z - dcal(ix,iy0).P.z) / (iy1-iy0);
		dcal(ix,iy).P.w += 1;
	}
	
	if (dcal(ix,iy).P.w == 0) {
		return;
	}
	dcal(ix,iy).P.x /= dcal(ix,iy).P.w;
	dcal(ix,iy).P.y /= dcal(ix,iy).P.w;
	dcal(ix,iy).P.z /= dcal(ix,iy).P.w;
	dcal(ix,iy).P.w = 1;
	dcal(ix,iy).State = dEye_Screen_Cal_Extra;
	return;
}

void	Screen_Cal_Eye_InterXY	(tScreen* p, tEye* peye, si ix, si iy)
{
	si tmp, ix0 = -1, ix1 = -1;
	for (tmp = ix-1; tmp >= 0; --tmp) {
		if (dcal(tmp,iy).State == dEye_Screen_Cal_Set) {
			ix0 = tmp;
			break;
		}
		if (ix0 < 0 && dcal(tmp,iy).State == dEye_Screen_Cal_Interp)
			ix0 = tmp;
	}
	for (tmp = ix+1; tmp < dEye_Screen_Cal_NUM; ++tmp) {
		if (dcal(tmp,iy).State == dEye_Screen_Cal_Set) {
			ix1 = tmp;
			break;
		}
		if (ix1 < 0 && dcal(tmp,iy).State == dEye_Screen_Cal_Interp)
			ix1 = tmp;
	}
	si iy0 = -1, iy1 = -1;
	for (tmp = iy-1; tmp >= 0; --tmp) {
		if (dcal(ix,tmp).State == dEye_Screen_Cal_Set) {
			iy0 = tmp;
			break;
		}
		if (iy0 < 0 && dcal(ix,tmp).State == dEye_Screen_Cal_Interp)
			iy0 = tmp;
	}
	for (tmp = iy+1; tmp < dEye_Screen_Cal_NUM; ++tmp) {
		if (dcal(ix,tmp).State == dEye_Screen_Cal_Set) {
			iy1 = tmp;
			break;
		}
		if (iy1 < 0 && dcal(ix,tmp).State == dEye_Screen_Cal_Interp) {
			iy1 = tmp;
		}
	}
	dcal(ix,iy).P.x = 0;
	dcal(ix,iy).P.y = 0;
	dcal(ix,iy).P.z = 0;
	dcal(ix,iy).P.w = 0;
	
	if (dcal(ix0,iy).State == dEye_Screen_Cal_Set && dcal(ix1,iy).State == dEye_Screen_Cal_Set
		&& (dcal(ix,iy0).State != dEye_Screen_Cal_Set || dcal(ix,iy1).State != dEye_Screen_Cal_Set)
	)
		iy0 = -1;
	else if (dcal(ix,iy0).State == dEye_Screen_Cal_Set && dcal(ix,iy1).State == dEye_Screen_Cal_Set
		&& (dcal(ix0,iy).State != dEye_Screen_Cal_Set || dcal(ix1,iy).State != dEye_Screen_Cal_Set)
	)
		ix0 = -1;
	
	if (ix0 >= 0 && ix1 >= 0) {
		dcal(ix,iy).P.x += dcal(ix0,iy).P.x + (ix-ix0) * (dcal(ix1,iy).P.x-dcal(ix0,iy).P.x) / (ix1-ix0);
		dcal(ix,iy).P.y += dcal(ix0,iy).P.y + (ix-ix0) * (dcal(ix1,iy).P.y-dcal(ix0,iy).P.y) / (ix1-ix0);
		dcal(ix,iy).P.z += dcal(ix0,iy).P.z + (ix-ix0) * (dcal(ix1,iy).P.z-dcal(ix0,iy).P.z) / (ix1-ix0);
	//	dcal(ix,iy).P.x = dcal(ix0,iy).P.x + dcal(ix1,iy).P.x;		dcal(ix,iy).P.x /= 2;
	//	dcal(ix,iy).P.y = dcal(ix0,iy).P.y + dcal(ix1,iy).P.y;		dcal(ix,iy).P.y /= 2;
	//	dcal(ix,iy).P.z = dcal(ix0,iy).P.z + dcal(ix1,iy).P.z;		dcal(ix,iy).P.z /= 2;
		dcal(ix,iy).P.w += 1;
	}
	if (iy0 >= 0 && iy1 >= 0) {
		dcal(ix,iy).P.x += dcal(ix,iy0).P.x + (iy-iy0) * (dcal(ix,iy1).P.x - dcal(ix,iy0).P.x) / (iy1-iy0);
		dcal(ix,iy).P.y += dcal(ix,iy0).P.y + (iy-iy0) * (dcal(ix,iy1).P.y - dcal(ix,iy0).P.y) / (iy1-iy0);
		dcal(ix,iy).P.z += dcal(ix,iy0).P.z + (iy-iy0) * (dcal(ix,iy1).P.z - dcal(ix,iy0).P.z) / (iy1-iy0);
	//	dcal(ix,iy).P.x = dcal(ix,iy0).P.x + dcal(ix,iy1).P.x;		dcal(ix,iy).P.x /= 2;
	//	dcal(ix,iy).P.y = dcal(ix,iy0).P.y + dcal(ix,iy1).P.y;		dcal(ix,iy).P.y /= 2;
	//	dcal(ix,iy).P.z = dcal(ix,iy0).P.z + dcal(ix,iy1).P.z;		dcal(ix,iy).P.z /= 2;
		dcal(ix,iy).P.w += 1;
	}
	if (dcal(ix,iy).P.w == 0) {
		return;
	}
	dcal(ix,iy).P.x /= dcal(ix,iy).P.w;
	dcal(ix,iy).P.y /= dcal(ix,iy).P.w;
	dcal(ix,iy).P.z /= dcal(ix,iy).P.w;
	dcal(ix,iy).P.w = 1;
	dcal(ix,iy).State = dEye_Screen_Cal_Interp;
	return;
}

void	Screen_Eye_InterClear	(tScreen* p)
{
	si ix, iy;
	for (iy = 0; iy < dEye_Screen_Cal_NUM; ++iy) {
		for (ix = 0; ix < dEye_Screen_Cal_NUM; ++ix) {
			if (gM.Left.aScreen[p->Idx].aaCal[iy][ix].State == dEye_Screen_Cal_Interp
				|| gM.Left.aScreen[p->Idx].aaCal[iy][ix].State == dEye_Screen_Cal_Extra
			) {
				gM.Left.aScreen[p->Idx].aaCal[iy][ix].State = dEye_Screen_Cal_NULL;
			}
			if (gM.Right.aScreen[p->Idx].aaCal[iy][ix].State == dEye_Screen_Cal_Interp
				|| gM.Right.aScreen[p->Idx].aaCal[iy][ix].State == dEye_Screen_Cal_Extra
			) {
				gM.Right.aScreen[p->Idx].aaCal[iy][ix].State = dEye_Screen_Cal_NULL;
			}
		}
	}
}
void	Screen_Eye_InterAll	(tScreen* p)
{
	si ix, iy;
	for (iy = 0; iy < dEye_Screen_Cal_NUM; ++iy) {
		for (ix = 0; ix < dEye_Screen_Cal_NUM; ++ix) {
			if (	gM.Left.aScreen[p->Idx].aaCal[iy][ix].State == dEye_Screen_Cal_NULL
				|| gM.Left.aScreen[p->Idx].aaCal[iy][ix].State == dEye_Screen_Cal_Interp
			) {
				Screen_Cal_Eye_InterXY (p, &gM.Left, ix, iy);
			}
			if (	gM.Right.aScreen[p->Idx].aaCal[iy][ix].State == dEye_Screen_Cal_NULL
				|| gM.Right.aScreen[p->Idx].aaCal[iy][ix].State == dEye_Screen_Cal_Interp
			) {
				Screen_Cal_Eye_InterXY (p, &gM.Right, ix, iy);
			}
		}
	}
}
void	Screen_Eye_ExtraAll	(tScreen* p)
{
	si ix, iy;
	for (iy = 0; iy < dEye_Screen_Cal_NUM; ++iy) {
		for (ix = 0; ix < dEye_Screen_Cal_NUM; ++ix) {
			if (	gM.Left.aScreen[p->Idx].aaCal[iy][ix].State == dEye_Screen_Cal_NULL
				|| gM.Left.aScreen[p->Idx].aaCal[iy][ix].State == dEye_Screen_Cal_Extra
			) {
				Screen_Cal_Eye_ExtraXY (p, &gM.Left, ix, iy);
			}
			if (	gM.Right.aScreen[p->Idx].aaCal[iy][ix].State == dEye_Screen_Cal_NULL
				|| gM.Right.aScreen[p->Idx].aaCal[iy][ix].State == dEye_Screen_Cal_Extra
			) {
				Screen_Cal_Eye_ExtraXY (p, &gM.Right, ix, iy);
			}
		}
	}
}


void	Screen_MarkersUpdate	(tScreen* p)
{
	tEye* peye = &gM.Left;
	si ix, iy;
	for (iy = 0; iy < dEye_Screen_Cal_NUM; ++iy) {
		for (ix = 0; ix < dEye_Screen_Cal_NUM; ++ix) {
		//	Marker_Hide (&gM.ScreenPoint[iy][ix]);
			
			switch (dcal(ix,iy).State) {
			case dEye_Screen_Cal_NULL:	Marker_ColorSet (&gM.ScreenPoint[iy][ix], 0xFFC00000);	break;
			case dEye_Screen_Cal_Set:	Marker_ColorSet (&gM.ScreenPoint[iy][ix], 0xFF00C000);	break;
			case dEye_Screen_Cal_Interp:	Marker_ColorSet (&gM.ScreenPoint[iy][ix], 0xFF808080);	break;
			case dEye_Screen_Cal_Extra:	Marker_ColorSet (&gM.ScreenPoint[iy][ix], 0xFF404040);	break;
			}
			if (iy == p->Cal.iy && ix == p->Cal.ix)
				Marker_ColorSet (&gM.ScreenPoint[iy][ix], 0xFFFFFFFF);
			
		/*	Marker_Show (&gM.ScreenPoint[iy][ix]);
			
			Marker_Move (&gM.ScreenPoint[iy][ix],
				p->Off.x + dcal(ix,iy).SX,
				p->Off.y + dcal(ix,iy).SY
			);*/
		//	XClearWindow (gM.X.pDisp, gM.ScreenPoint[iy][ix].Win);
		}
	}
}

void	Screen_MarkersShow	(tScreen* p)
{
	tEye* peye = &gM.Left;
	si ix, iy;
	for (iy = 0; iy < dEye_Screen_Cal_NUM; ++iy) {
		for (ix = 0; ix < dEye_Screen_Cal_NUM; ++ix) {
		//	Marker_Hide (&gM.ScreenPoint[iy][ix]);
			
			Marker_Show (&gM.ScreenPoint[iy][ix]);
			
			Marker_Move (&gM.ScreenPoint[iy][ix],
				p->Off.x + dcal(ix,iy).SX,
				p->Off.y + dcal(ix,iy).SY
			);
		//	XClearWindow (gM.X.pDisp, gM.ScreenPoint[iy][ix].Win);
		//	XClearArea (gM.X.pDisp, gM.ScreenPoint[iy][ix].Win, 0, 0, 0, 0, True);
		}
	}
}
void	Screen_MarkersHide	(tScreen* p)
{
	tEye* peye = &gM.Left;
	si ix, iy;
	for (iy = 0; iy < dEye_Screen_Cal_NUM; ++iy) {
		for (ix = 0; ix < dEye_Screen_Cal_NUM; ++ix) {
			Marker_Hide (&gM.ScreenPoint[iy][ix]);
		}
	}
}


void	Screen_Cal_Prep		(tScreen* p)
{
	p->Cal.sx = gM.Left.aScreen[p->Idx].aaCal[p->Cal.iy][p->Cal.ix].SX;
	p->Cal.sy = gM.Left.aScreen[p->Idx].aaCal[p->Cal.iy][p->Cal.ix].SY;
//	printf ("Screen Cal Prep %ld %ld\n", p->Cal.sx, p->Cal.sy);
	Marker_Move (&gM.CalPoint,
		p->Off.x + p->Cal.sx,
		p->Off.y + p->Cal.sy
	);
	
	Screen_MarkersUpdate (p);
	p->bGood = 1;
}

void	Screen_Cal_Next		(tScreen* p)
{
	++p->Cal.ix;
	if (p->Cal.ix >= dEye_Screen_Cal_NUM) {
		p->Cal.ix = 0;
		++p->Cal.iy;
		if (p->Cal.iy >= dEye_Screen_Cal_NUM) {
			p->Cal.iy = 0;
			p->bGood = 1;
			Marker_Hide (&gM.CalPoint);
		}
	}
}
void	Screen_Cal_DoThis		(tScreen* p, si ix, si iy)
{
	p->Cal.ix = ix;
	p->Cal.iy = iy;
	Screen_Cal_Prep (p);
}
void	Screen_Cal_DoLeft		(tScreen* p)
{
	if (p->Cal.ix > 0)
		--p->Cal.ix;
	Screen_Cal_Prep (p);
}
void	Screen_Cal_DoRight	(tScreen* p)
{
	if (p->Cal.ix < dEye_Screen_Cal_LAST)
		++p->Cal.ix;
	Screen_Cal_Prep (p);
}
void	Screen_Cal_DoUp		(tScreen* p)
{
	if (p->Cal.iy > 0)
		--p->Cal.iy;
	Screen_Cal_Prep (p);
}
void	Screen_Cal_DoDown		(tScreen* p)
{
	if (p->Cal.iy < dEye_Screen_Cal_LAST)
		++p->Cal.iy;
	Screen_Cal_Prep (p);
}
void	Screen_Cal_DoPoint	(tScreen* p, si x, si y)
{
	si gridx = p->PixW / (dEye_Screen_Cal_NUM-1);
	si gridy = p->PixH / (dEye_Screen_Cal_NUM-1);
	
	x -= p->Off.x - gridx/2;
	y -= p->Off.y - gridy/2;
	
	si ix = x / gridx;
	si iy = y / gridy;
	
	if (ix < 0)
		p->Cal.ix = 0;
	if (ix > dEye_Screen_Cal_LAST)
		p->Cal.ix = dEye_Screen_Cal_LAST;
	else
		p->Cal.ix = ix;
	
	if (iy < 0)
		p->Cal.iy = 0;
	if (iy > dEye_Screen_Cal_LAST)
		p->Cal.iy = dEye_Screen_Cal_LAST;
	else
		p->Cal.iy = iy;
	Screen_Cal_Prep (p);
}



void	Screen_GCal_Save_Eye	(tScreen* p, tEye* peye)
{
	si ii = p->Cal.ix + p->Cal.iy * dEye_Screen_Cal_NUM;
	peye->Homo.aPoint_State[ii] = 1;
	peye->Homo.scenecalipoints[ii].x = dcal(p->Cal.ix, p->Cal.iy).SX;
	peye->Homo.scenecalipoints[ii].y = dcal(p->Cal.ix, p->Cal.iy).SY;
	
	SDL_mutexP (peye->GV_mutex);
	peye->Homo.vectors[ii].x = peye->GV.x;
	peye->Homo.vectors[ii].y = peye->GV.y;
	SDL_mutexV (peye->GV_mutex);
	
	ii = 0;
	for (; ii < CALIBRATIONPOINTS; ++ii) {
		if (!peye->Homo.aPoint_State[ii])
			return;
	}
	CalculateCalibration(&peye->Homo);
	return;
}
void	Screen_GCal_Save		(tScreen* p)
{
//	printf ("Screen_Cal_Save Idx %d\n", p->Idx);
	
//	printf ("hoho1 ");	V4f_Print (&ret);
//	printf ("vec ");	V4f_Print (&vec);
	
	Screen_GCal_Save_Eye (p, &gM.Left);
	Screen_GCal_Save_Eye (p, &gM.Right);
	
}


void	Screen_ReInterp		(tScreen* p)
{
	Screen_Eye_InterClear (p);
	Screen_Eye_InterAll (p);
	Screen_Eye_InterAll (p);
	Screen_Eye_InterAll (p);
	
	Screen_Eye_ExtraAll (p);
}

void	Screen_Cal_Save_Eye	(tScreen* p, tEye* peye, float z)
{
	printf ("Screen_Cal_Save_Eye Z %f\n", z);
//	peye->aScreen[p->Idx].aaCal[p->Cal.iy][p->Cal.ix].SX = p->Cal.sx;
//	peye->aScreen[p->Idx].aaCal[p->Cal.iy][p->Cal.ix].SY = p->Cal.sy;
	
	tV4f leye, lret, lvec;
/*	leye = peye->InHead.P;
	Head_Eye_Vector (&gM.Head, peye, &lret);
	lvec = lret;	V4f_sub_V4f (&lvec, &peye->InHead.P);
	
	M4f_mul_V4f (&gM.Head.M4, &leye);	M4f_mul_V4f (&gM.Head.M4, &lret);	M4f_mul_V4f (&gM.Head.M4, &lvec);/**/
	Head_Eye_VectorGlob (&gM.Head, peye, &leye, &lret, &lvec);
	
	tV4f pos = (tV4f){0, 0, z, 1};
	tV4f norm = (tV4f){0, 0, 1, 1};
	
	V4f_Intersect_Line01_Plane0N (
		&peye->aScreen[p->Idx].aaCal[p->Cal.iy][p->Cal.ix].P,
		&leye, &lret,
		&pos, &norm
	);
	
//	printf ("%d %d P: ", p->Cal.ix, p->Cal.iy);	V4f_Print (&peye->aScreen[p->Idx].aaCal[p->Cal.iy][p->Cal.ix].P);
//	peye->aScreen[p->Idx].aaCal[p->Cal.iy][p->Cal.ix].P = (tV4f){gM.aScreen[p->Idx].C.x - gM.aScreen[p->Idx].W/2, gM.aScreen[p->Idx].C.y - gM.aScreen[p->Idx].H/2, gM.aScreen[p->Idx].C.z, 1};
	
	peye->aScreen[p->Idx].aaCal[p->Cal.iy][p->Cal.ix].State = dEye_Screen_Cal_Set;
}
void	Screen_Cal_Save		(tScreen* p)
{
	if (dEye_Screen_Cal_NUM == 3) {
		Screen_GCal_Save (p);
	}
	float z = gM.aScreen[p->Idx].C.z;
	if (0) {
		tV4f e00, e01, e10, e11;
		tV4f p0, p1;
		
		Head_Eye_VectorGlob (&gM.Head, &gM.Left, &e00, &e01, NULL);
		Head_Eye_VectorGlob (&gM.Head, &gM.Right, &e10, &e11, NULL);
		
		V4f_Intersect_Line01_Line01 (
			&p0, &p1,
			&e00, &e01,
			&e10, &e11
		);
		z = (p0.z + p1.z) / 2;
		printf ("Screen_Cal_Save Z %f\n", z);
	}
	
//	return;
//	printf ("Screen_Cal_Save Idx %d\n", p->Idx);
	
//	printf ("hoho1 ");	V4f_Print (&ret);
//	printf ("vec ");	V4f_Print (&vec);
	
	Screen_Cal_Save_Eye (p, &gM.Left, z);
	Screen_Cal_Save_Eye (p, &gM.Right, z);
	
	Screen_ReInterp (p);
	
	Screen_Cal_Next (p);
	Screen_Cal_Prep (p);
}



void	Screen_Cal_Record		(tScreen* p, tEye* peye, float z)
{
//	printf ("Screen_Cal_Record Z %f\n", z);
	
	tV4f leye, lret, lvec;
	
	Head_Eye_VectorGlob (&gM.Head, peye, &leye, &lret, &lvec);
	
	tV4f pos = (tV4f){0, 0, z, 1};
	tV4f norm = (tV4f){0, 0, 1, 1};		V4f_rotx (&norm, -30*deg2rad);
	tV4f ret;
	
	
	V4f_Intersect_Line01_Plane0N (
		&ret,
		&leye, &lret,
		&pos, &norm
	);
//	printf ("Screen_Cal_Record Z %f	ret\n", ret.z);
	
	V4f_add_V4f (&peye->aScreen[p->Idx].aaCal[p->Cal.iy][p->Cal.ix].P, &ret);
}
void	Screen_Cal_RecordStart	(tScreen* p)
{
	gM.Cal_Record_Num = 0;
	
	gM.Left.aScreen[p->Idx].aaCal[p->Cal.iy][p->Cal.ix].P.x = 0;
	gM.Left.aScreen[p->Idx].aaCal[p->Cal.iy][p->Cal.ix].P.y = 0;
	gM.Left.aScreen[p->Idx].aaCal[p->Cal.iy][p->Cal.ix].P.z = 0;
	gM.Left.aScreen[p->Idx].aaCal[p->Cal.iy][p->Cal.ix].P.w = 1;
	
	gM.Right.aScreen[p->Idx].aaCal[p->Cal.iy][p->Cal.ix].P.x = 0;
	gM.Right.aScreen[p->Idx].aaCal[p->Cal.iy][p->Cal.ix].P.y = 0;
	gM.Right.aScreen[p->Idx].aaCal[p->Cal.iy][p->Cal.ix].P.z = 0;
	gM.Right.aScreen[p->Idx].aaCal[p->Cal.iy][p->Cal.ix].P.w = 1;
	
}
void	Screen_Cal_RecordEnd	(tScreen* p)
{
	gM.Left.aScreen[p->Idx].aaCal[p->Cal.iy][p->Cal.ix].P.x /= (float)gM.Cal_Record_Num;
	gM.Left.aScreen[p->Idx].aaCal[p->Cal.iy][p->Cal.ix].P.y /= (float)gM.Cal_Record_Num;
	gM.Left.aScreen[p->Idx].aaCal[p->Cal.iy][p->Cal.ix].P.z /= (float)gM.Cal_Record_Num;
	gM.Left.aScreen[p->Idx].aaCal[p->Cal.iy][p->Cal.ix].P.w = 1;
	
	gM.Right.aScreen[p->Idx].aaCal[p->Cal.iy][p->Cal.ix].P.x /= (float)gM.Cal_Record_Num;
	gM.Right.aScreen[p->Idx].aaCal[p->Cal.iy][p->Cal.ix].P.y /= (float)gM.Cal_Record_Num;
	gM.Right.aScreen[p->Idx].aaCal[p->Cal.iy][p->Cal.ix].P.z /= (float)gM.Cal_Record_Num;
	gM.Right.aScreen[p->Idx].aaCal[p->Cal.iy][p->Cal.ix].P.w = 1;
	
	gM.Left.aScreen[p->Idx].aaCal[p->Cal.iy][p->Cal.ix].State = dEye_Screen_Cal_Set;
	gM.Right.aScreen[p->Idx].aaCal[p->Cal.iy][p->Cal.ix].State = dEye_Screen_Cal_Set;
	
	Screen_ReInterp (p);
	
	Screen_Cal_Next (p);
	Screen_Cal_Prep (p);
}




void	Screen_Cal_PointDel	(tScreen* p)
{
	
	gM.Left.aScreen[p->Idx].aaCal[p->Cal.iy][p->Cal.ix].State = dEye_Screen_Cal_NULL;
	gM.Right.aScreen[p->Idx].aaCal[p->Cal.iy][p->Cal.ix].State = dEye_Screen_Cal_NULL;
	
	
	Screen_ReInterp (p);
	
//	Screen_Cal_Next (p);
	Screen_Cal_Prep (p);
}





void	Screen_Cal_Init		(tScreen* p)
{
	p->Cal.ix = 0;
	p->Cal.iy = 0;
	
//	Marker_Show (&gM.CalPoint);
	Screen_Cal_Prep (p);
	Screen_MarkersUpdate (p);
	Screen_MarkersShow (p);
}


void	Screen_Eye_Init		(tScreen* p, tEye* peye)
{
	si ix, iy;
	for (iy = 0; iy < dEye_Screen_Cal_NUM; ++iy) {
		for (ix = 0; ix < dEye_Screen_Cal_NUM; ++ix) {
			dcal(ix,iy).SX = ix * p->PixW/dEye_Screen_Cal_LAST;
			switch (ix) {
			case 0:				dcal(ix,iy).SX += gM.CalPoint.Win_W/2;	break;
			case dEye_Screen_Cal_LAST:	dcal(ix,iy).SX -= gM.CalPoint.Win_W/2;	break;
			}
			
			dcal(ix,iy).SY = iy * p->PixH/dEye_Screen_Cal_LAST;
			switch (iy) {
			case 0:				dcal(ix,iy).SY += gM.CalPoint.Win_W/2;	break;
			case dEye_Screen_Cal_LAST:	dcal(ix,iy).SY -= gM.CalPoint.Win_W/2;	break;
			}
		}
	}
}


u08	Screen_Eye_Point	(tEye* peye, tV4f* ppos, float* ps, float *pt)
{
	tV4f p0, p1, t0 = {0, 0, 0, 1}, t1, t2;
/*	Head_Eye_Vector (&gM.Head, peye, &p1);
	
	p0 = peye->InHead.P;
	
	M4f_mul_V4f (&gM.Head.M4, &p0);
	M4f_mul_V4f (&gM.Head.M4, &p1);/**/
	Head_Eye_VectorGlob (&gM.Head, peye, &p0, &p1, NULL);
	
	t0.z = gM.aScreen[0].C.z;
	t2 = t1 = t0;
	t1.x = 1;
	t2.y = 1;
	return V4f_Intersect_Line01_Tri012 (
		ppos, ps, pt,
		&p0, &p1,
	//	&peye->aScreen[0].aaCal[0][1].P, &peye->aScreen[0].aaCal[0][0].P, &peye->aScreen[0].aaCal[1][1].P
		&t0, &t1, &t2
	);
}

s08	Screen_Eye_XY_T	(tScreen* p, tEye* peye, si ix, si iy, float *px, float *py,  float* ps, float* pt)
{
	tV4f p0, p1;
/*	Head_Eye_Vector (&gM.Head, peye, &p1);
	
	p0 = peye->InHead.P;
	
	M4f_mul_V4f (&gM.Head.M4, &p0);
	M4f_mul_V4f (&gM.Head.M4, &p1);/**/
	
	Head_Eye_VectorGlob (&gM.Head, peye, &p0, &p1, NULL);
	
	float s, t;
	V4f_Intersect_Line01_Tri012 (
		0, &s, &t,
		&p0, &p1,
		&peye->aScreen[p->Idx].aaCal[iy+0][ix+0].P, &peye->aScreen[p->Idx].aaCal[iy+0][ix+1].P, &peye->aScreen[p->Idx].aaCal[iy+1][ix+0].P
	);
//	printf ("Screen_Eye_XY_T	st %f\t%f\n", s, t);
/*	printf ("t0		");	V4f_Print (&dcal(ix+0,iy+0).P);
	printf ("t1		");	V4f_Print (&dcal(ix+1,iy+0).P);
	printf ("t2		");	V4f_Print (&dcal(ix+0,iy+1).P);/**/
	
	*px = peye->aScreen[p->Idx].aaCal[iy+0][ix+0].SX + s * (peye->aScreen[p->Idx].aaCal[iy+0][ix+1].SX - peye->aScreen[p->Idx].aaCal[iy+0][ix+0].SX);
	*py = peye->aScreen[p->Idx].aaCal[iy+0][ix+0].SY + t * (peye->aScreen[p->Idx].aaCal[iy+1][ix+0].SY - peye->aScreen[p->Idx].aaCal[iy+0][ix+0].SY);
	*ps = s;
	*pt = t;
	if (s < -0.0001)
		return -1;
	if (t < -0.0001)
		return -1;
	if (!isfinite(s+t) || s+t > 1)
		return 1;
	return 0;
}

s08	Screen_Eye_XY_B	(tScreen* p, tEye* peye, si ix, si iy, float *px, float *py, float* ps, float* pt)
{
	tV4f p0, p1;
/*	Head_Eye_Vector (&gM.Head, peye, &p1);
	
	p0 = peye->InHead.P;
	
	M4f_mul_V4f (&gM.Head.M4, &p0);
	M4f_mul_V4f (&gM.Head.M4, &p1);/**/
	
	Head_Eye_VectorGlob (&gM.Head, peye, &p0, &p1, NULL);
	
	float s, t;
	V4f_Intersect_Line01_Tri012 (
		0, &s, &t,
		&p0, &p1,
		&peye->aScreen[p->Idx].aaCal[iy+1][ix+1].P, &peye->aScreen[p->Idx].aaCal[iy+1][ix+0].P, &peye->aScreen[p->Idx].aaCal[iy+0][ix+1].P
	);
	
//	printf ("Screen_Eye_XY_B	st %f\t%f\n", s, t);
/*	printf ("t0		");	V4f_Print (&dcal(ix+1,iy+1).P);
	printf ("t1		");	V4f_Print (&dcal(ix+0,iy+1).P);
	printf ("t2		");	V4f_Print (&dcal(ix+1,iy+0).P);/**/
	
	*px = peye->aScreen[p->Idx].aaCal[iy+1][ix+1].SX + s * (peye->aScreen[p->Idx].aaCal[iy+1][ix+0].SX - peye->aScreen[p->Idx].aaCal[iy+1][ix+1].SX);
	*py = peye->aScreen[p->Idx].aaCal[iy+1][ix+1].SY + t * (peye->aScreen[p->Idx].aaCal[iy+0][ix+1].SY - peye->aScreen[p->Idx].aaCal[iy+1][ix+1].SY);
	*ps = s;
	*pt = t;
	
	if (s < -0.0001)
		return -1;
	if (t < -0.0001)
		return -1;
	if (!isfinite(s+t) || s+t > 1)
		return 1;
	return 0;
}

s08	Screen_Eye_XY	(tScreen* p, tEye* peye, float *px, float *py)
{
//	printf ("Screen_Eye_XY\n");
	float ax = 0, ay = 0;
	si num = 0;
	si ix, iy;
	for (iy = 0; iy < dEye_Screen_Cal_NUM-1; ++iy) {
		for (ix = 0; ix < dEye_Screen_Cal_NUM-1; ++ix) {
			float x0, y0, s, t, tx, ty;
			s08 ret;
			ret = Screen_Eye_XY_T (p, peye, ix, iy, &x0, &y0, &s, &t);
			if (ret == 0) {
				*px = x0;
				*py = y0;
				return 1;
			}else if (ret == -1 && (ix == 0 || iy == 0)) {
				if (s < 0 && ix != 0)
					continue;
				if (t < 0 && iy != 0)
					continue;
				if (s > 1 && ix != dEye_Screen_Cal_LAST-1)
					continue;
				if (t > 1 && iy != dEye_Screen_Cal_LAST-1)
					continue;
				*px = x0;
				*py = y0;
				num = -1;	//TODO add proper picking closest outside tri
			}else {
				ret = Screen_Eye_XY_B (p, peye, ix, iy, &x0, &y0, &s, &t);
				if (ret == 0) {
					*px = x0;
					*py = y0;
					return 1;
				}else if (ret == -1 && (ix == dEye_Screen_Cal_LAST-1 || iy == dEye_Screen_Cal_LAST-1)) {
				//	printf ("on ix %d iy %d got %f %f\n", ix, iy, s, t);
					if (s < 0 && ix != dEye_Screen_Cal_LAST-1)
						continue;
					if (t < 0 && iy != dEye_Screen_Cal_LAST-1)
						continue;
					if (s > 1 && ix != 0)
						continue;
					if (t > 1 && iy != 0)
						continue;
					*px = x0;
					*py = y0;
					num = -1;
				}else {
					continue;
				}
			}
		}
	}
//	printf ("num %ld\n", num);
	return num;
}



void	Screen_Eye_PrintSub	(tScreen* p, tEye* peye)
{
	if (1) {
		ui ix = 0, iy = 0;
		for (iy = 1; iy < dEye_Screen_Cal_NUM; ++iy) {
			for (ix = 0; ix < dEye_Screen_Cal_NUM; ++ix) {
				if (dcal(ix,iy-1).State == dEye_Screen_Cal_NULL
					|| dcal(ix,iy).State == dEye_Screen_Cal_NULL
				)
					continue;
				tV4f p0 = dcal(ix,iy-1).P;
				tV4f p1 = dcal(ix,iy).P;
				Dbg_V4f_ADrawPosPos (&p0, &p1);
			}
		}
		for (iy = 0; iy < dEye_Screen_Cal_NUM; ++iy) {
			for (ix = 1; ix < dEye_Screen_Cal_NUM; ++ix) {
				if (dcal(ix-1,iy).State == dEye_Screen_Cal_NULL
					|| dcal(ix,iy).State == dEye_Screen_Cal_NULL
				)
					continue;
				tV4f p0 = dcal(ix-1,iy).P;
				tV4f p1 = dcal(ix,iy).P;
				Dbg_V4f_ADrawPosPos (&p0, &p1);
			}
		}
		for (iy = 0; iy < dEye_Screen_Cal_NUM-1; ++iy) {
			for (ix = 0; ix < dEye_Screen_Cal_NUM-1; ++ix) {
				if (dcal(ix,iy).State == dEye_Screen_Cal_NULL
					|| dcal(ix+1,iy+1).State == dEye_Screen_Cal_NULL
				)
					continue;
				tV4f p0 = dcal(ix,iy).P;
				tV4f p1 = dcal(ix+1,iy+1).P;
				Dbg_V4f_ADrawPosPos (&p0, &p1);
			}
		}
	}
	if (1) {
		float s, t;
		tV4f p0;
		Screen_Eye_Point (peye, &p0, &s, &t);
		
		Dbg_V4f_ADrawPoint (&p0);
	}/**/
}

void	Screen_Eye_Print	(tScreen* p, tEye* peye)
{
//	printf ("Screen_Eye_Print Idx %d\n", p->Idx);
	Col_EyeSet (peye);
	
	if (1) {
		Screen_Eye_PrintSub (p, peye);
	}
	if (1) {
		ui ix = 0, iy = 0;
		for (iy = 0; iy < dEye_Screen_Cal_NUM-1; ++iy) {
			for (ix = 0; ix < dEye_Screen_Cal_NUM-1; ++ix) {
				tV4f p0 = dcal(ix,iy-1).P;
				tV4f p1 = dcal(ix,iy).P;
				Dbg_V4f_ADrawPosPos (&p0, &p1);
			}
		}
	}
}

void	Screen_Print	(tScreen* p)
{
	Screen_Eye_Print (p, &gM.Left);
	Screen_Eye_Print (p, &gM.Right);
	
	if (p == gM.aScreen + 0) {
		tEye* peye = &gM.Left;
		printf ("width top %f\n", V4f_dist_V4f (&dcal(0,0).P, &dcal(dEye_Screen_Cal_LAST,0).P));
		printf ("width bot %f\n", V4f_dist_V4f (&dcal(0,dEye_Screen_Cal_LAST).P, &dcal(dEye_Screen_Cal_LAST,dEye_Screen_Cal_LAST).P));
		
		printf ("heigh lef %f\n", V4f_dist_V4f (&dcal(0,0).P, &dcal(0,dEye_Screen_Cal_LAST).P));
		printf ("heigh rig %f\n", V4f_dist_V4f (&dcal(dEye_Screen_Cal_LAST,0).P, &dcal(dEye_Screen_Cal_LAST,dEye_Screen_Cal_LAST).P));
		
	}
	
}

#undef dcal




void	Head_Eye_Print	(tHead* p, tEye* peye)	//print the eye vectors
{
	Col_EyeSet (peye);
	
	tV4f eye, ret, vec;
	if (1) {
		Head_Eye_VectorGlob (p, peye, &eye, NULL, &vec);
		V4f_mul_S (&vec, 60);
		V4f_add_V4f (&vec, &eye);
		
		Dbg_V4f_ADrawPosPos (&eye, &vec);
		Dbg_V4f_ADrawPoint (&eye);
		
		if (1) {
			Head_Eye_VectorGlob (p, peye, &eye, NULL, &vec);
			V4f_mul_S (&vec, peye->InHead.R/V4f_dist(&vec));
			V4f_add_V4f (&vec, &eye);
			
			Dbg_V4f_ADrawPoint (&vec);
			
			gCol.R = 0x0;
			gCol.G = 0xFF;
			gCol.B = 0x0;
			V4f_DrawPosPos (pcam, &vec, &vec);
		}
	}
	
	if (1) {
	//	float r = peye->InHead.R;
		tV4f p0, vec, p1;
	//	Head_Eye_Line (p, peye, &p0, &vec);
		
		Cam_Retina_Ray (peye, &p0, &p1, &vec);
		
		V4f_mul_S (&p1, 20);
		
		Dbg_V4f_ADrawPosPos (&p0, &p1);
		Dbg_V4f_ADrawPoint (&p0);
	}
}

void	Head_Dbg_Print	(tHead* p)
{
	if (1 && isfinite(gM.Head.P.x)) {
		tV4f p0 = gM.Head.P; V4f_add_V4f (&p0, &gM.Head.N);
		tV4f pc = gM.Head.Mod.PC;
		tV4f pl = gM.Head.Mod.PL;
		tV4f pr = gM.Head.Mod.PR;/**/
		
		M4f_mul_V4f (&gM.Head.M4, &pc);
		M4f_mul_V4f (&gM.Head.M4, &pl);
		M4f_mul_V4f (&gM.Head.M4, &pr);
		
		
		Dbg_V4f_ADrawPosPos (&gM.Head.P, &p0);
		Dbg_V4f_ADrawPosPos (&pl, &pc);
		Dbg_V4f_ADrawPosPos (&pr, &pc);
		Dbg_V4f_ADrawPosPos (&pl, &pr);
		
	}
	//Head_Eye_Print (p, &gM.Left);
	//Head_Eye_Print (p, &gM.Right);
	
	if (0) {
		tV4f e00, e01, e10, e11;
		tV4f p0, p1;
		
		Head_Eye_VectorGlob (p, &gM.Left, &e00, &e01, NULL);
		Head_Eye_VectorGlob (p, &gM.Right, &e10, &e11, NULL);
		
		V4f_Intersect_Line01_Line01 (
			&p0, &p1,
			&e00, &e01,
			&e10, &e11
		);
		
		gColARGB = 0xFFFFFF;
	//	Dbg_V4f_ADrawPosPos (&p0, &p1);
		Dbg_V4f_ADrawPoint (&p0);
		Dbg_V4f_ADrawPoint (&p1);
	}
	/**/
}

void	Head_Eye_Gaze	(tHead* p, tEye* peye)
{
	tV4f eye, ret, vec;
	Head_Eye_Vector (p, peye, &ret);
	
	eye = peye->InHead.P;
	vec = ret;	V4f_sub_V4f (&vec, &peye->InHead.P);
	
	M4f_mul_V4f (&p->M4, &eye);
	M4f_mul_V4f (&p->M4, &ret);
	
//	printf ("Head_Eye_Gaze\n");
//	s08 good[gM.Screen_N];
//	tV2f gaze[gM.Screen_N];
	si i, an = 0, mind = 600;	s08 got = 0;
	float ax = 0, ay = 0;
	
	for (i = 0; i < gM.Screen_N; ++i) {
		if (!gM.aScreen[i].bGood)
			continue;
		float x, y;
		s08 ret = Screen_Eye_XY (gM.aScreen + i, peye, &x, &y);
	//	printf ("Screen_Eye_XY ret %d\n", ret);
		if (ret) {
			si d = 0;
			if (x < 0) {
				d += -x;
				x = 0;
			}else if (x >= gM.aScreen[i].PixW) {
				d += x-gM.aScreen[i].PixW+1;
				x = gM.aScreen[i].PixW-1;
			}
			if (y < 0) {
				d += -y;
				y = 0;
			}else if (y >= gM.aScreen[i].PixH) {
				d += y-gM.aScreen[i].PixH+1;
				y = gM.aScreen[i].PixH-1;
			}
			if (got != 1) {
				if (!d) {
					got = 1;
					ax = gM.aScreen[i].Off.x+x;
					ay = gM.aScreen[i].Off.y+y;
					an = 1;
					break;
				}else if (d < mind) {
					got = -1;
					mind = d;
					ax = gM.aScreen[i].Off.x+x;
					ay = gM.aScreen[i].Off.y+y;
					an = 1;
				}
			}
		}
	}
//	printf ("an %ld\n", an);
	ax /= an;
	ay /= an;
	
	if (peye == &gM.Left) {
	//	printf ("Left\t%f\t%f\n", ax, ay);
		gM.GazeL.x = ax;
		gM.GazeL.y = ay;
	}else {
	//	printf ("Righ\t%f\t%f\n", ax, ay);
		gM.GazeR.x = ax;
		gM.GazeR.y = ay;
	}
}



void	Glint_Sphere_Dbg_Print	(tV4f* ppos, float r)
{
	tV4f cam = {0, 0, 0, 1}, cp = *ppos;
	
	gColARGB = 0xFFFF00;
	
	Dbg_V4f_ADrawPoint (&cp);
	
	tV4f l0, l1;
	
	Sphere_Reflect (&l0, &cp, r, &gM.aLight[0].P);
	Sphere_Reflect (&l1, &cp, r, &gM.aLight[1].P);
	
//	printf ("l0 rad %f\n", V4f_dist_V4f(&cp, &l0));
//	printf ("l1 rad %f\n", V4f_dist_V4f(&cp, &l1));
	
	gColARGB = 0xFFFF00;
	Dbg_V4f_ADrawPosPos (&cp, &l0);
	Dbg_V4f_ADrawPosPos (&cp, &l1);
	
	gColARGB = 0xFF00FF;
	
	Dbg_V4f_ADrawPosPos (&gM.aLight[0].P, &l0);
	Dbg_V4f_ADrawPosPos (&gM.aLight[1].P, &l1);
	Dbg_V4f_ADrawPosPos (&cam, &l0);
	Dbg_V4f_ADrawPosPos (&cam, &l1);
	Dbg_V4f_ADrawPoint (&l0);
	Dbg_V4f_ADrawPoint (&l1);
	
}

void	Eye_Glint_Dbg_Print	(tEye* peye)
{
	abort();
	/*
	Glint_Sphere_Dbg_Print (&peye->GP, peye->LR);
	
	Col_EyeSet (peye);
	tV4f p0, vec, p1;
	
	Cam_Pos2Ray (peye->P, &p0, &p1, &vec);
	
	V4f_mul_S (&p1, 20);
	
//	Dbg_V4f_ADrawPosPos (&p0, &p1);
//	Dbg_V4f_ADrawPoint (&p0);
	
	if (gM.Eye_GlintMode == 2) {
		Cam_Pos2Ray (peye->G0, &p0, &p1, &vec);
		V4f_mul_S (&p1, 20);
		Dbg_V4f_ADrawPosPos (&p0, &p1);
		Cam_Pos2Ray (peye->G1, &p0, &p1, &vec);
		V4f_mul_S (&p1, 20);
		Dbg_V4f_ADrawPosPos (&p0, &p1);
		
		
	//	Dbg_V4f_ADrawPoint (&peye->GP);
		
		if (1) {
			tV4f eye, ret, vec;
			Eye_GP_GazeVector (peye, &eye, NULL, &vec);
			V4f_mul_S (&vec, 60);
			V4f_add_V4f (&vec, &eye);
			
			Dbg_Ortho_Front ();
			Dbg_V4f_ADrawPosPos (&eye, &vec);
			Dbg_V4f_DrawPoint (&eye);
			
			Dbg_Ortho_Top ();
			Dbg_V4f_ADrawPosPos (&eye, &vec);
			Dbg_V4f_DrawPoint (&eye);
			
			Dbg_Ortho_Left ();
			Dbg_V4f_ADrawPosPos (&eye, &vec);
			Dbg_V4f_DrawPoint (&eye);
		}
	}*/
}

void	Glint_Dbg_Print	()
{
	gColARGB = 0xFFFFFF;
	Dbg_V4f_ADrawPoint (&gM.aLight[0].P);
	Dbg_V4f_ADrawPoint (&gM.aLight[1].P);
//	Cam_Pos2Ray (, &p0, &p1, &vec);
	
	if (0) {
		si i;
		tV4f cam = {0, 0, 0, 1}, cp = {gM.tmp.x, gM.tmp.y, gM.tmp.z, 1};
		float cr = 1.5f;
		
		gColARGB = 0xFFFF00;
		
		Dbg_V4f_ADrawPoint (&cp);
		
		tV4f l0, l1;
		
		Sphere_Reflect (&l0, &cp, cr, &gM.aLight[0].P);
		Sphere_Reflect (&l1, &cp, cr, &gM.aLight[1].P);
		
	//	printf ("l0 rad %f\n", V4f_dist_V4f(&cp, &l0));
	//	printf ("l1 rad %f\n", V4f_dist_V4f(&cp, &l1));
		
		gColARGB = 0xFFFF00;
		Dbg_V4f_ADrawPosPos (&cp, &l0);
		Dbg_V4f_ADrawPosPos (&cp, &l1);
		
		gColARGB = 0xFF00FF;
		
		Dbg_V4f_ADrawPosPos (&gM.aLight[0].P, &l0);
		Dbg_V4f_ADrawPosPos (&gM.aLight[1].P, &l1);
		Dbg_V4f_ADrawPosPos (&cam, &l0);
		Dbg_V4f_ADrawPosPos (&cam, &l1);
		Dbg_V4f_ADrawPoint (&l0);
		Dbg_V4f_ADrawPoint (&l1);
	}
	
	Eye_Glint_Dbg_Print (&gM.Left);
	Eye_Glint_Dbg_Print (&gM.Right);
	
	if (1) {
		tV4f e00, e01, e10, e11;
		tV4f p0, p1;
		
		Head_Eye_VectorGlob (&gM.Head, &gM.Left, &e00, &e01, NULL);
		Head_Eye_VectorGlob (&gM.Head, &gM.Right, &e10, &e11, NULL);
		
		V4f_Intersect_Line01_Line01 (
			&p0, &p1,
			&e00, &e01,
			&e10, &e11
		);
		
		gColARGB = 0xFFFFFF;
	//	Dbg_V4f_ADrawPosPos (&p0, &p1);
		Dbg_V4f_ADrawPoint (&p0);
		Dbg_V4f_ADrawPoint (&p1);
	}
}


void	Eye_Glint_Gaze	(tEye* peye)
{
/*	tV4f eye, ret, vec;
	Head_Eye_Vector (p, peye, &ret);
	
	eye = peye->InHead.P;
	vec = ret;	V4f_sub_V4f (&vec, &peye->InHead.P);
	
	M4f_mul_V4f (&p->M4, &eye);
	M4f_mul_V4f (&p->M4, &ret);
	
//	s08 good[gM.Screen_N];
//	tV2f gaze[gM.Screen_N];
	si i, an = 0;
	float ax = 0, ay = 0;
	for (i = 0; i < gM.Screen_N; ++i) {
		if (!gM.aScreen[i].bGood)
			continue;
		float x, y;
		if (Screen_Eye_XY (gM.aScreen + i, peye, &x, &y)) {
			ax += gM.aScreen[i].Off.x+x;
			ay += gM.aScreen[i].Off.y+y;
			++an;
		}
	}
//	printf ("an %ld\n", an);
	ax /= an;
	ay /= an;
	
	if (peye == &gM.Left) {
	//	printf ("Left\t%f\t%f\n", ax, ay);
		gM.GazeL.x = ax;
		gM.GazeL.y = ay;
	}else {
	//	printf ("Righ\t%f\t%f\n", ax, ay);
		gM.GazeR.x = ax;
		gM.GazeR.y = ay;
	}/**/
}



void	Cam_Proj_Conf		(tCam* pcam)
{
	if (0) {
		float fd = dpow2(pcam->Full_W) + dpow2(pcam->Full_H);
		float id = dpow2(pcam->Image_Zoom*pcam->Image_W) + dpow2(pcam->Image_Zoom*pcam->Image_H);
	//	float id = ddist(pcam->Image_W, pcam->Image_H);
		pcam->Image_FOV = id/fd * pcam->Full_FOV;
	}
	if (1) {
		float fd = pcam->Full_W;
		float id = pcam->Image_Zoom*pcam->Image_W;
	//	float id = ddist(pcam->Image_W, pcam->Image_H);
		pcam->Image_FOV = id/fd * pcam->Full_FOV;
	}
	
	float fov = pcam->Image_FOV*deg2rad;
	float a = atan2(pcam->Image_H, pcam->Image_W);
	
	if (pcam->Image_FOV_W == 0)
		pcam->Image_FOV_W = fov*cos(a);
	if (pcam->Image_FOV_H == 0)
		pcam->Image_FOV_H = fov*sin(a);
	
	printf ("Image_FOV %f  W %f H %f	a %f\n", pcam->Image_FOV, pcam->Image_FOV_W*rad2deg, pcam->Image_FOV_H*rad2deg, a);
	
	pcam->Proj_W = pcam->Proj_N*tan(pcam->Image_FOV_W/2.0f);
	pcam->Proj_H = pcam->Proj_N*tan(pcam->Image_FOV_H/2.0f);
	
	pcam->Proj_L = -pcam->Proj_W;
	pcam->Proj_R = pcam->Proj_W;
	pcam->Proj_B = -pcam->Proj_H;
	pcam->Proj_T = pcam->Proj_H;
	
	pcam->Proj_W *= 2;
	pcam->Proj_H *= 2;
	
/*	pcam->Proj.x00 = pcam->Proj_N / (pcam->Proj_W/2);
	pcam->Proj.x11 = pcam->Proj_N / (pcam->Proj_H/2);
	
	pcam->Proj.x22 = -(pcam->Proj_F+pcam->Proj_N) / (pcam->Proj_F - pcam->Proj_N);
	
	pcam->Proj.x23 = (-2*pcam->Proj_F*pcam->Proj_N) / (pcam->Proj_F - pcam->Proj_N);
	pcam->Proj.x32 = -1;/**/
	M4f_Frustrum (&pcam->Proj, pcam->Proj_W, pcam->Proj_H, pcam->Proj_N, pcam->Proj_F);
	
	
	pcam->cvCam = cvCreateMatHeader(3, 3, CV_32FC1);
	cvCreateData(pcam->cvCam);
	cvSet2D (pcam->cvCam, 0, 0, cvScalarAll(2058.4890867164763));
	cvSet2D (pcam->cvCam, 0, 1, cvScalarAll(0));
	cvSet2D (pcam->cvCam, 0, 2, cvScalarAll(959.50000000000000));
	cvSet2D (pcam->cvCam, 1, 0, cvScalarAll(0));
	cvSet2D (pcam->cvCam, 1, 1, cvScalarAll(2058.4890867164763));
	cvSet2D (pcam->cvCam, 1, 2, cvScalarAll(539.50000000000000));
	cvSet2D (pcam->cvCam, 2, 0, cvScalarAll(0));
	cvSet2D (pcam->cvCam, 2, 1, cvScalarAll(0));
	cvSet2D (pcam->cvCam, 2, 2, cvScalarAll(1));
	
	pcam->cvDist = cvCreateMatHeader(5, 1, CV_32FC1);
	cvCreateData(pcam->cvDist);
	cvSet2D (pcam->cvDist, 0, 0, cvScalarAll(-0.023321480666943440));
	cvSet2D (pcam->cvDist, 1, 0, cvScalarAll(0.69231202113335466));
	cvSet2D (pcam->cvDist, 2, 0, cvScalarAll(0));
	cvSet2D (pcam->cvDist, 3, 0, cvScalarAll(0));
	cvSet2D (pcam->cvDist, 4, 0, cvScalarAll(-3.5279699461953733));
	
	
	PrintMat (pcam->cvCam);
	PrintMat (pcam->cvDist);
}

int	Cam_Loop		(tCam* pcam)
{
	while (gM.aCam[0].UVC->signalquit) {
		SDL_SemWait (gM.sWaitForUpdate);
		
	/*	for (int y = 0; y < pcam->SDL_Surf->h; ++y) {
			for (int x = 0; x < pcam->SDL_Surf->w; ++x) {
				tPix* pix = (tPix*)pcam->UVC->framebuffer + x + y* gM.apCam[i]->UVC->width;
				uint32_t* out = (uint32_t*)pcam->SDL_Surf->pixels + x + y* pcam->SDL_Surf->w;
				*out = (pix->Y << 8) | (pix->Y << 0) | (pix->Y << 16);
			}
		}/**/
		//printf("cam %d write start\n", pcam->Idx);
		SDL_LockSurface (pcam->SDL_Surf);
		for (int y = 0; y < pcam->SDL_Surf->h; ++y) {
			for (int x = 0; x < pcam->SDL_Surf->w; ++x) {
				tPix* pix = (tPix*)pcam->UVC->framebuffer + x*2 + y*2 * pcam->UVC->width;
				uint32_t* out = (uint32_t*)pcam->SDL_Surf->pixels + x + y* pcam->SDL_Surf->w;
				uint16_t val = (pix[0].Y + pix[1].Y) / 2;
				*out = (val << 8) | (val << 0) | (val << 16);
			}
		}
		
		Eye_V_Pre (&gM.Head.DotC, pcam);
		Eye_Fit (&gM.Head.DotC, pcam);
		Eye_V_Post (&gM.Head.DotC, pcam);
		Eye_Clip (&gM.Head.DotC);
		
	//	memcpy(gM.pDst, videoIn->framebuffer, videoIn->width * (videoIn->height) * 2);
		
		Eye_V_Pre (&gM.Head.DotL, pcam);
		Eye_Fit (&gM.Head.DotL, pcam);
		Eye_V_Post (&gM.Head.DotL, pcam);
		Eye_Clip (&gM.Head.DotL);
		
	//	memcpy(gM.pDst, videoIn->framebuffer, videoIn->width * (videoIn->height) * 2);
		
		Eye_V_Pre (&gM.Head.DotR, pcam);
		Eye_Fit (&gM.Head.DotR, pcam);
		Eye_V_Post (&gM.Head.DotR, pcam);
		Eye_Clip (&gM.Head.DotR);
		
	//	memcpy(gM.pDst, videoIn->framebuffer, videoIn->width * (videoIn->height) * 2);
		/*
	//	Eye_V_Pre (&gM.Right);
		Eye_Fit (&gM.Right);
	//	Eye_V_Post (&gM.Right);
		Eye_Clip (&gM.Right);
		
	//	memcpy(gM.pDst, videoIn->framebuffer, videoIn->width * (videoIn->height) * 2);
		
	//	Eye_V_Pre (&gM.Left);
		Eye_Fit (&gM.Left);
	//	Eye_V_Post (&gM.Left);
		Eye_Clip (&gM.Left);
		
		Eye_Draw (&gM.Left);
		Eye_Draw (&gM.Right);
	//	printf ("Eye pos Left %f %f\n", gM.Left.P.x, gM.Left.P.y);
	//	printf ("Eye pos Righ %f %f\n", gM.Right.P.x, gM.Right.P.y);
		*/
		if (pcam->Idx == 0) {
			//printf ("DotC aCam[0].P %f %f\n", gM.Head.DotC.aCam[0].P.x, gM.Head.DotC.aCam[0].P.y);
			
			HeadC_Calc_M4_Rel (&gM.Head, pcam, &gM.HeadC);
			
			#define dcopy(_name)	gM.Head._name = gM.Head.aCam[pcam->Idx]._name
			
			dcopy(P);		dcopy(N);		dcopy(R);
			dcopy(M4);		dcopy(M4_T);	dcopy(M4_R);
			dcopy(M4I);		dcopy(M4I_T);	dcopy(M4I_R);
			dcopy(R_X);		dcopy(R_Y);		dcopy(R_Z);
			dcopy(SRX);		dcopy(SRY);
			#undef dcopy
			
			Head_Dbg_Print(&gM.Head);
		}
		HeadC_Draw (&gM.Head, pcam);
		
		#define dprojline3(x0,y0,z0,x1,y1,z1)	\
			do {	\
				tV4f __p0 = {x0, y0, z0, 1}, __p1 = {x1, y1, z1, 1};	\
				V4f_DrawPosPos (pcam, &__p0, &__p1);	\
			}while(0)
		
		dprojline3(0, 0, -40, 10, 0, -40);
		
		dprojline3(0, 0, -40, 0, 10, -40);
		
		
		#undef dprojline3
		
		Eye_Draw (&gM.Head.DotC, pcam);
		Eye_Draw (&gM.Head.DotL, pcam);
		Eye_Draw (&gM.Head.DotR, pcam);
		
		SDL_UnlockSurface (pcam->SDL_Surf);
		//printf("cam %d write start end\n", pcam->Idx);
		//muhaha();
		
		SDL_SemPost (gM.sWaitForCams);
	}
}

dyn_config gM_DC;
int	muhaha_Config_Thread	(void *data)
{
	dyn_config_watch (&gM_DC, "config.yaml");
	return 0;
}


#if dM_Actions_Mode == 2
int	muhaha_chrdev_Thread	(void *data)
{
	int fd = open ("/dev/zeawesome", O_RDONLY);
	if (!fd)
		exit(-1);
	
	while (1) {
		char buf[32];
		ssize_t size = read(fd, &buf, 1);
		printf ("YEAH got character 0x%hhx 0x%hhu\n", buf[0], buf[0] & 0x7F);
		do {
			if (buf[0] & 0x80) {
				buf[0] &= 0x7F;
				#define dact(name,press_stuff)	\
					if (buf[0] == gM.Action.name) {	\
						press_stuff	\
						break;	\
					}
				#define dact2(name,press_stuff,release_stuff)	dact(name,press_stuff)
				
				#include "actions.h"
				#undef dact
				#undef dact2
			}else {
				#define dact(name,press_stuff)	
				#define dact2(name,press_stuff,release_stuff)	\
					if (buf[0] == gM.Action.name) {	\
						release_stuff	\
						break;	\
					}
				
				#include "actions.h"
				#undef dact
				#undef dact2
			}
		}while (0);
	}
	
	close (fd);
	return 0;
}
#endif

int	muhaha_XEvent_Thread	(void *data)
{
	printf ("first\n");
	while (pcam->UVC->signalquit) {
		XEvent e;
		XNextEvent(gM.X.pDisp, &e);
		switch (e.type) {
		case Expose:
			Marker_Paint (&gM.CrossHair);
			Marker_Paint (&gM.CalPoint);
			
			si ix, iy;
			for (iy = 0; iy < dEye_Screen_Cal_NUM; ++iy) {
				for (ix = 0; ix < dEye_Screen_Cal_NUM; ++ix) {
					Marker_Paint (&gM.ScreenPoint[iy][ix]);
				}
			}
			break;
		#if dM_Actions_Mode == 1
		case KeyPress: {
			KeySym key;
			char buffer[20];
			int bufsize = 20;
			int charcount = XLookupString(&e, buffer, bufsize, &key, 0);
			
			#define dact(name,press_stuff)	\
				if (key == gM.Action.name) {	\
					press_stuff	\
					break;	\
				}
			#define dact2(name,press_stuff,release_stuff)	dact(name,press_stuff)
			
			#include "actions.h"
			#undef dact
			#undef dact2
			break;
		}
		case KeyRelease: {
			if (XEventsQueued (gM.X.pDisp, QueuedAfterReading))
			{
				XEvent nev;
				XPeekEvent(gM.X.pDisp, &nev);
				
				if (nev.type == KeyPress && nev.xkey.time == e.xkey.time &&
					nev.xkey.keycode == e.xkey.keycode)
				{
				//	fprintf (stdout, "key #%ld was retriggered.\n", (long) XLookupKeysym (&nev.xkey, 0));
					// delete retriggered KeyPress event
					XNextEvent (gM.X.pDisp, &nev);
					break;
				}
			}
			KeySym key;
			char buffer[20];
			int bufsize = 20;
			int charcount = XLookupString(&e, buffer, bufsize, &key, 0);
			
			#define dact(name,press_stuff)	
			#define dact2(name,press_stuff,release_stuff)	\
				if (key == gM.Action.name) {	\
					release_stuff	\
					break;	\
				}
			
			#include "actions.h"
			#undef dact
			#undef dact2
			break;
		}
		#endif
		}
		XSync (gM.X.pDisp, False);
		
	/*	while (gM.X.Queue_N) {
			si idx = gM.X.Queue_C;
			
			switch (gM.X.aQueue[idx].Act) {
			case dM_X_Queue_M_Down:
				xdo_mousedown(gM.X.pDo, CURRENTWINDOW, gM.X.aQueue[idx].Button);
				break;
			case dM_X_Queue_M_Up:
				xdo_mouseup(gM.X.pDo, CURRENTWINDOW, gM.X.aQueue[idx].Button);
				break;
			}
			--gM.X.Queue_N;
			++gM.X.Queue_C;
			gM.X.Queue_C %= dM_X_Queue_NUM;
		}/**/
		
		
	//	printf ("hahahahha\n");
	//	XFlush(gM.X.pDisp);/**/
	}
}



void	muhaha_Keyboard_Init	()
{
	
	#if dM_Actions_Mode == 1
	u32 mods[] = {
		Mod4Mask,
		Mod4Mask | Mod2Mask,
	};
	#define dgrab(_keycode)	\
		do {		\
			int j;		\
			for (j = 0; j < sizeof(mods)/sizeof(mods[0]); j++) {		\
				printf ("Key grab code %d mods 0x%x\n", _keycode, mods[j]);/**/	\
				XGrabKey(gM.X.pDisp, _keycode, mods[j], gM.aScreen[0].Win, True, GrabModeAsync, GrabModeAsync);		\
			}		\
		}while(0)
	#define dgrab(_keysym)	\
		do {		\
			int j;		\
			for (j = 0; j < sizeof(mods)/sizeof(mods[0]); j++) {		\
				int ret = XGrabKey(gM.X.pDisp, XKeysymToKeycode(gM.X.pDisp, _keysym), mods[j], gM.aScreen[0].Win, True, GrabModeAsync, GrabModeAsync);		\
				printf ("Key grab sym 0x%x code %d mods 0x%x ret %d\n", _keysym, XKeysymToKeycode(gM.X.pDisp, _keysym), mods[j], ret);/**/	\
				XSync(gM.X.pDisp, False);	\
			}		\
		}while(0)
	
	#define dact(name,press_stuff)	dgrab(gM.Action.name);
	#define dact2(name,press_stuff,release_stuff) dact(name,press_stuff)
	#include "actions.h"
	#undef dact
	#undef dact2
//	XGrabKey(gM.X.pDisp, gM.Action.Mod_Key, Mod2Mask, gM.aScreen[0].Win, True, GrabModeAsync, GrabModeAsync);		\
	
	#undef dgrab
	XSync(gM.X.pDisp, False);
	#endif
}
extern char** environ;
void	muhaha_Init	()
{
	if (!XInitThreads()) {
		printf("Xlib not thread safe\n");
		exit(1);
	}
	{
		char** p = environ;
		for (; ; ++p) {
			if (!*p) {
				printf("NO DISPLAY environment variable\n");
			}
			if (strncmp(*p, "DISPLAY=", 8) == 0) {
				gM.X.pDisp = XOpenDisplay((*p)+8);
			//	gM.X.pDo = xdo_new_with_opened_display (gM.X.pDisp, (*p)+8, 0);
				break;
			}
		}/**/
	//	gM.X.pDisp = XOpenDisplay(":0");
		if (gM.X.pDisp == NULL) {
			fprintf(stderr, "Couldn't open display\n");
			exit (1);
		}
	}
	
	gM.sWaitForCams = SDL_CreateSemaphore(0);
	gM.sWaitForUpdate = SDL_CreateSemaphore(0);
	
	gM.X.Queue_N = 0;
	gM.X.Queue_C = 0;
	
	gM.Draw_X = 0;
	gM.Draw_Y = 0;
	gM.Draw_W = 800;
	gM.Draw_H = 600;
	
	gM.Cal_bRecord = 0;
	
	gM.mutMainIsOut = SDL_CreateMutex();
	gM.mutMainIsIn = SDL_CreateMutex();
	
	Head_Init (&gM.Head);
	
	Eye_Init (&gM.Left);
	Eye_Init (&gM.Right);
	
	gM.pGaze_Mutex = SDL_CreateMutex();
	
	dyn_config_init(&gM_DC);
	dyn_config_read(&gM_DC, "config.yaml");
	/*mythread = */SDL_CreateThread(muhaha_Config_Thread, "muhaha_Config_Thread", (void *)NULL);
	
	
	for (int i = 0; i < gM.Cam_N; ++i) {
		tCam* pcam = &gM.aCam[i];
		Cam_Param_Set (pcam);
		//Cam_Param_Set (pcam);
		
		pcam->View_W = pcam->SDL_Surf->w;
		pcam->View_H = pcam->SDL_Surf->h;
		
		Cam_Proj_Conf (pcam);
		
		M4f_Iden (&pcam->World);
		if (i == 0)
			M4f_trans (&pcam->World, 0, 0, 0);
		else
			M4f_trans (&pcam->World, 0, -2.3, 0);
	//	M4f_rotx (&pcam->World, 15*deg2rad);
	//	M4f_trans (&pcam->World, 0, 0, -15);
	//	M4f_rotx (&pcam->World, 15*deg2rad);
		Cam_Param_Set (pcam);
	}
	
/*	for (gM.X.Screen_N = 0; gM.X.Screen_N < ScreenCount(gM.X.pDisp); ++gM.X.Screen_N) {
		gM.X.aScreen[gM.X.Screen_N].Win = RootWindow(gM.X.pDisp, gM.X.Screen_N);
		
		XGetWindowAttributes (gM.X.pDisp, gM.X.aScreen[gM.X.Screen_N].Win, &gM.X.aScreen[gM.X.Screen_N].Attr);
		
		printf ("Screen %d  xywh %d %d %d %d\n", gM.X.Screen_N,
			gM.X.aScreen[gM.X.Screen_N].Attr.x,
			gM.X.aScreen[gM.X.Screen_N].Attr.y,
			gM.X.aScreen[gM.X.Screen_N].Attr.width,
			gM.X.aScreen[gM.X.Screen_N].Attr.height
		);
	}/**/
//	Screen* screen = ScreenOfDisplay(display,0);
//	printf ("w %d\n", WidthOfScreen(screen));
	
	si i;
	for (i = 0; i < dScreen_MAX; ++i) {
		gM.aScreen[i].Idx = i;
		gM.aScreen[i].bGood = 0;
		
		gM.aScreen[i].Win = DefaultRootWindow(gM.X.pDisp);
		XGetWindowAttributes (gM.X.pDisp, gM.aScreen[i].Win, &gM.aScreen[i].Attr);
		
		gM.Left.aScreen[i].View.Top.x = 1100;
		gM.Left.aScreen[i].View.Top.y = 300;
		
		gM.Right.aScreen[i].View.Top.x = 1100;
		gM.Right.aScreen[i].View.Top.y = 300;
		
		
		gM.Left.aScreen[i].View.Left.x = 1300;
		gM.Left.aScreen[i].View.Left.y = 300;
		
		gM.Right.aScreen[i].View.Left.x = 1300;
		gM.Right.aScreen[i].View.Left.y = 300;
		
		
		gM.Left.aScreen[i].View.Front.x = 1100;
		gM.Left.aScreen[i].View.Front.y = 500;
		
		gM.Right.aScreen[i].View.Front.x = 1100;
		gM.Right.aScreen[i].View.Front.y = 500;
		
	}
	
	for (i = 0; i < M_Dbg_NUM; ++i) {
		tDbg* pdbg = &gM.aDbg[i];
		char buf[64];
		sprintf (buf, "Dbg %d", i);
		pdbg->SDL_Win = SDL_CreateWindow(buf,
			SDL_WINDOWPOS_UNDEFINED, SDL_WINDOWPOS_UNDEFINED,
			320, 240,
			SDL_WINDOW_SHOWN | SDL_WINDOW_OPENGL
		);
		if (!pdbg->SDL_Win)
			abort();
		
		pdbg->SDL_Surf = SDL_GetWindowSurface(pdbg->SDL_Win);
		if (!pdbg->SDL_Surf)
			abort();
		
		Dbg_M4f_Conf (pdbg);
	}
	
//	Window window = RootWindow(display, DefaultScreen(display));
//	Window window = DefaultRootWindow(gM.X.pDisp);
//	printf ("root %lx\n", window);
	
/*	XIAddMasterInfo add;
	add.type = XIAddMaster;
	add.name = "muhaha";
	add.send_core = False;
	add.enable = True;
	Status ret = XIChangeHierarchy (gM.X.pDisp, &add, 1);
	printf ("add master %d\n", ret);
	/**/
	
	XIDeviceInfo *devices;
	int device_num;
	
	devices = XIQueryDevice (gM.X.pDisp, XIAllDevices, &device_num);
	
	if ( device_num <= 0 ){
		return;
	}
	
	for (i = 0; i < device_num; ++i) {
		XIDeviceInfo *pdev;
	/*	typedef struct
		{
			int                 deviceid;
			char                *name;
			int                 use;
			int                 attachment;
			Bool                enabled;
			int                 num_classes;
			XIAnyClassInfo      **classes;
		} XIDeviceInfo;/**/
		pdev = devices + i;
		
		printf ("Dev %d  u %d %s\n", pdev->deviceid, pdev->use, pdev->name);
		if (strcmp ("muhaha pointer", pdev->name) == 0) {
			gM.X.Pointer = *pdev;
		}
		if (strcmp ("Virtual core pointer", pdev->name) == 0) {
			gM.X.Pointer = *pdev;
		}
		if (strcmp ("Virtual core keyboard", pdev->name) == 0) {
			gM.X.Keyboard = *pdev;
		}
	//	if (pdev->use == IsXPointer) {
	//	}
	}
	
//	Cursor cur = 0;
//	cur = XCreateFontCursor(gM.X.pDisp, XC_left_ptr);
//	printf ("cur %d\n", cur);
	
//	printf ("define cur %d\n", XIDefineCursor (gM.X.pDisp, gM.X.Pointer.deviceid, gM.aScreen[0].Win, cur));
//	printf ("define cur %d\n", XIDefineCursor (gM.X.pDisp, gM.X.Pointer.deviceid, 0x3400001, cur));
	
	XIFreeDeviceInfo (devices);
	
	XFlush(gM.X.pDisp);
	XSync(gM.X.pDisp, False);
//	XSync(gM.X.pDisp, True);
//	XFlush(gM.X.pDisp);
	
	
	si ix, iy;
	for (iy = 0; iy < dEye_Screen_Cal_NUM; ++iy) {
		for (ix = 0; ix < dEye_Screen_Cal_NUM; ++ix) {
			Marker_Init (&gM.ScreenPoint[iy][ix], 1);
		}
	}
	Marker_Init (&gM.CrossHair, 0);
	Marker_Init (&gM.CalPoint, 1);
	
	
	for (i = 0; i < dScreen_MAX; ++i) {
		Screen_Eye_Init (gM.aScreen + i, &gM.Left);
		Screen_Eye_Init (gM.aScreen + i, &gM.Right);
	}
	
	
	if (gM.Pointer_Mode == 1)
		Marker_Show (&gM.CrossHair);
	
	
	muhaha_Keyboard_Init ();
	
	
	/*mythread = */SDL_CreateThread(muhaha_XEvent_Thread, "muhaha_XEvent_Thread", (void *)NULL);
	
	#if dM_Actions_Mode == 2
	/*mythread = */SDL_CreateThread(muhaha_chrdev_Thread, "muhaha_chrdev_Thread", (void *)NULL);
	#endif
	
	/*mythread =*/ SDL_CreateThread(muhaha_eventThread, "muhaha_eventThread", (void *) &ptdata);
}


void	muhaha_Print	()
{
	printf ("EYE POINTS!!!!\n");
	printf ("Left:\n");
	printf ("    P:\n");
	printf ("      x: %f\n", gM.Left.InHead.P.x);
	printf ("      y: %f\n", gM.Left.InHead.P.y);
	printf ("      z: %f\n", gM.Left.InHead.P.z);
	printf ("Right:\n");
	printf ("    P:\n");
	printf ("      x: %f\n", gM.Right.InHead.P.x);
	printf ("      y: %f\n", gM.Right.InHead.P.y);
	printf ("      z: %f\n", gM.Right.InHead.P.z);
}

void	muhaha_DeInit	()
{
	Status ret;
	
/*	XIRemoveMasterInfo rem;
	rem.type = XIRemoveMaster;
	rem.deviceid = gM.X.Pointer.deviceid;
	rem.return_mode = XIFloating;
	rem.return_pointer = 0;
	rem.return_keyboard = 0;
	
	ret = XIChangeHierarchy (gM.X.pDisp, &rem, 1);
	printf ("rem master %d\n", ret);
	XFlush(gM.X.pDisp);/**/
//	xdo_free (gM.X.pDo);
	
	muhaha_Print ();
}

void	muhaha_Cross	(tV2f* ppoint, tPix col)
{
	float x, y;
	#define dd 1
	for (y = ppoint->y-dd; y <= ppoint->y+dd; ++y) {
		for (x = ppoint->x-dd; x <= ppoint->x+dd; ++x) {
			if (dpixout(x,y))
				continue;
			*dnpix(x,y) = col;
		}
	}
	#undef dd
	#define dd 4
	x = ppoint->x;
	for (y = ppoint->y-dd; y <= ppoint->y+dd; ++y) {
		if (dpixout(x,y))
			continue;
		*dnpix(x,y) = col;
	}
	y = ppoint->y;
	for (x = ppoint->x-dd; x <= ppoint->x+dd; ++x) {
		if (dpixout(x,y))
			continue;
		*dnpix(x,y) = col;
	}
	#undef dd
}

si fps = 0, tfps = 0;
void	muhaha	()
{
	SDL_mutexV (gM.mutMainIsIn);
	
	SDL_mutexP (gM.mutMainIsOut);
	
	++tfps;
	
	tPix* pin = (tPix*)pcam->UVC->framebuffer, *pin1;
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
		
		for (y = 0; y < pcam->UVC->height; ++y) {
			for (x = 0; x < pcam->UVC->width-1; ++x) {
				pin = (tPix*)pcam->UVC->framebuffer + (y*pcam->UVC->width + x);
				pin1 = pin+1;
			//	pin->Y = 0;
				si y = pin->Y, y1 = pin1->Y;
			//	pin->Y = abs(y1 - y);
				
			//	pin->Y = 0x80;
			//	pin->U = ((si)pin1->U - ((si)pin->U)) +7;
			//	pin->V = ((si)pin1->V - ((si)pin->V)) +8;
			//	pin->V = abs((si)pin1->V - (si)pin->V);
				
				pin->U = 0x7;
				pin->V = 0x7;
			//	dpix(x,y)->U = 0x7;
			//	dpix(x,y)->V = 0x8;
			}
		}
		memcpy(gM.pDst, pcam->UVC->framebuffer, pcam->UVC->width * (pcam->UVC->height) * 2);
		/**/
	/*	memcpy(gM.pDst, videoIn->framebuffer, videoIn->width * (videoIn->height) * 2);
		
		float or = 55;
		for (y = or+1; y < videoIn->height-or-1; y += 30) {
			for (x = or+1; x < videoIn->width-or-1; x += 30) {
				tV2si pos;
				pos.x = x;
				pos.y = y;
			//	printf ("%f   ", NN_Sphere (&pos, or-2, or, 0, 0));
				float res = 10*NN_Sphere (&pos, or-5, or, 0, 0);
				res = dclip_l (res, 0);
				res = dclip_h (res, 255);
				dnpix(x,y)->Y = res;
			//	dnpix(x,y)->Y = 255;
				dnpix(x,y)->U = 0;
				dnpix(x,y)->V = 0;
			}
		//	printf ("\n");
		}/**/
	}else {
		//memcpy(gM.pDst, pcam->UVC->framebuffer, pcam->UVC->width * (pcam->UVC->height) * 2);
	}
/*	printf ("World:\n");
	M4f_Print (&pcam->World);
	printf ("Proj:\n");
	M4f_Print (&pcam->Proj);/**/
	/*
	for (y = 0; y < 600; ++y)
		for (x = 800; x < 800+640; ++x)
			dspix(x, y) = 0;
	*/
	
	//Track_Point_Step (gM.Head.aPoint + 0);
	
	/*
	if (1) {
	//	memcpy(gM.pDst, videoIn->framebuffer, videoIn->width * (videoIn->height) * 2);
		
		Eye_V_Pre (&gM.Head.DotC);
		Eye_Fit (&gM.Head.DotC);
		Eye_V_Post (&gM.Head.DotC);
		Eye_Clip (&gM.Head.DotC);
		
	//	memcpy(gM.pDst, videoIn->framebuffer, videoIn->width * (videoIn->height) * 2);
		
		Eye_V_Pre (&gM.Head.DotL);
		Eye_Fit (&gM.Head.DotL);
		Eye_V_Post (&gM.Head.DotL);
		Eye_Clip (&gM.Head.DotL);
		
	//	memcpy(gM.pDst, videoIn->framebuffer, videoIn->width * (videoIn->height) * 2);
		
		Eye_V_Pre (&gM.Head.DotR);
		Eye_Fit (&gM.Head.DotR);
		Eye_V_Post (&gM.Head.DotR);
		Eye_Clip (&gM.Head.DotR);
		
	//	memcpy(gM.pDst, videoIn->framebuffer, videoIn->width * (videoIn->height) * 2);
		
	//	Eye_V_Pre (&gM.Right);
		Eye_Fit (&gM.Right);
	//	Eye_V_Post (&gM.Right);
		Eye_Clip (&gM.Right);
		
	//	memcpy(gM.pDst, videoIn->framebuffer, videoIn->width * (videoIn->height) * 2);
		
	//	Eye_V_Pre (&gM.Left);
		Eye_Fit (&gM.Left);
	//	Eye_V_Post (&gM.Left);
		Eye_Clip (&gM.Left);
		
		Eye_Draw (&gM.Left);
		Eye_Draw (&gM.Right);
	//	printf ("Eye pos Left %f %f\n", gM.Left.P.x, gM.Left.P.y);
	//	printf ("Eye pos Righ %f %f\n", gM.Right.P.x, gM.Right.P.y);
		
		Eye_Draw (&gM.Head.DotC);
		Eye_Draw (&gM.Head.DotL);
		Eye_Draw (&gM.Head.DotR);
	}
	*/
/*	if (0) {
		gCol.Y = 0xFF;
		gCol.U = 0x0;
		gCol.V = 0x0;
		
		#define dheadline(x0,y0,x1,y1)	\
			do {	\
				tV2f __p0 = {x0, y0}, __p1 = {x1, y1};	\
				V2f_mul_M3f (&__p0, &gM.Head.M);	V2f_mul_M3f (&__p1, &gM.Head.M);	\
				V2f_DrawPosPos (&__p0, &__p1);	\
			}while(0)
		
		tV2f pc = {0, 0};
		tV2f lx = {40, 0};
		tV2f ly = {0, 20};
		
		V2f_mul_M3f (&pc, &gM.Head.M);
		V2f_mul_M3f (&lx, &gM.Head.M);
		V2f_mul_M3f (&ly, &gM.Head.M);
		
		V2f_DrawPosPos (&pc, &lx);
		V2f_DrawPosPos (&pc, &ly);
	}*/
	
//	dheadline (20,30,40,50);
	
	#define dprojline3(x0,y0,z0,x1,y1,z1)	\
		do {	\
			tV4f __p0 = {x0, y0, z0, 1}, __p1 = {x1, y1, z1, 1};	\
			V4f_DrawPosPos (pcam, &__p0, &__p1);	\
		}while(0)
	
//	dprojline2(0,0,100,100);
	if (0) {
		gCol.R = 0x0;
		gCol.G = 0xFF;
		gCol.B = 0x0;
		static float ang = 0, zz = -30;
		
		float d = 7.0/2.0, z = 3*d;
		
	//	tM4f rot2 = {cosf(ang),sinf(ang),0,0, -sinf(ang),cosf(ang),0,0, 0,0,1,0, 0,0,0,1};
	/*	tM4f rot = {1,0,0,0, 0,cosf(ang),sinf(ang),0, 0,-sinf(ang),cosf(ang),0, 0,0,0,1};
		tM4f rot2 = {cosf(ang),0,sinf(ang),0, 0,1,0,0, -sinf(ang),0,cosf(ang),0, 0,0,0,1};
		M4f_mul_M4f (&rot, &rot2);
		
		tM4f trans = {1,0,0,0, 0,1,0,0, 0,0,1,0, 0,0,zz,1};*/
		
		tM4f trans, rot;
		M4f_Iden (&trans);
		M4f_trans (&trans, 0, -3, zz);
		M4f_roty (&trans, -45*deg2rad);
		M4f_rotx (&trans, 45*deg2rad);
		
		tV4f p00 = {-d, -d, -d, 1};
		tV4f p01 = {+d, -d, -d, 1};
		tV4f p02 = {+d, +d, -d, 1};
		tV4f p03 = {-d, +d, -d, 1};
		
		tV4f p10 = {-d, -d,  d, 1};
		tV4f p11 = {+d, -d,  d, 1};
		tV4f p12 = {+d, +d,  d, 1};
		tV4f p13 = {-d, +d,  d, 1};
		
	/*	V4f_mul_M4f (&p00, &rot);
		V4f_mul_M4f (&p01, &rot);
		V4f_mul_M4f (&p02, &rot);
		V4f_mul_M4f (&p03, &rot);
		V4f_mul_M4f (&p10, &rot);
		V4f_mul_M4f (&p11, &rot);
		V4f_mul_M4f (&p12, &rot);
		V4f_mul_M4f (&p13, &rot);/**/
		
	/*	V4f_mul_M4f (&p00, &trans);
		V4f_mul_M4f (&p01, &trans);
		V4f_mul_M4f (&p02, &trans);
		V4f_mul_M4f (&p03, &trans);
		V4f_mul_M4f (&p10, &trans);
		V4f_mul_M4f (&p11, &trans);
		V4f_mul_M4f (&p12, &trans);
		V4f_mul_M4f (&p13, &trans);/**/
		
	/*	M4f_mul_V4f (&rot, &p00);
		M4f_mul_V4f (&rot, &p01);
		M4f_mul_V4f (&rot, &p02);
		M4f_mul_V4f (&rot, &p03);
		M4f_mul_V4f (&rot, &p10);
		M4f_mul_V4f (&rot, &p11);
		M4f_mul_V4f (&rot, &p12);
		M4f_mul_V4f (&rot, &p13);/**/
		
		M4f_mul_V4f (&trans, &p00);
		M4f_mul_V4f (&trans, &p01);
		M4f_mul_V4f (&trans, &p02);
		M4f_mul_V4f (&trans, &p03);
		M4f_mul_V4f (&trans, &p10);
		M4f_mul_V4f (&trans, &p11);
		M4f_mul_V4f (&trans, &p12);
		M4f_mul_V4f (&trans, &p13);
		
	//	V4f_Print (&p00);
	//	printf ("closer\n");
		V4f_DrawPosPos (pcam, &p00, &p01);
		V4f_DrawPosPos (pcam, &p01, &p02);
		V4f_DrawPosPos (pcam, &p02, &p03);
		V4f_DrawPosPos (pcam, &p03, &p00);
		
	//	printf ("mid\n");
		V4f_DrawPosPos (pcam, &p00, &p10);
		V4f_DrawPosPos (pcam, &p01, &p11);
		V4f_DrawPosPos (pcam, &p02, &p12);
		V4f_DrawPosPos (pcam, &p03, &p13);
		
	//	printf ("far\n");
		V4f_DrawPosPos (pcam, &p10, &p11);
		V4f_DrawPosPos (pcam, &p11, &p12);
		V4f_DrawPosPos (pcam, &p12, &p13);
		V4f_DrawPosPos (pcam, &p13, &p10);
		
		ang += 10*M_PI/180.0f;
		//zz += 1.0f;
		//if (zz >= z*2)
			//zz = z/2;
	}/**/
	
	if (0) {
		gCol.R = 0x0;
		gCol.G = 0xFF;
		gCol.B = 0x0;
		
		float ang = 0, zz = 13*3;
		
		float d = 5.5/2, z = 0;
		
		
	//	tM4f rot2 = {cosf(ang),sinf(ang),0,0, -sinf(ang),cosf(ang),0,0, 0,0,1,0, 0,0,0,1};
	//	tM4f rot = {1,0,0,0, 0,cosf(ang),sinf(ang),0, 0,-sinf(ang),cosf(ang),0, 0,0,0,1};
	//	tM4f rot2 = {cosf(ang),0,sinf(ang),0, 0,1,0,0, -sinf(ang),0,cosf(ang),0, 0,0,0,1};
	//	M4f_mul_M4f (&rot, &rot2);
		
	//	tM4f trans = {1,0,0,0, 0,1,0,0, 0,0,1,0, 0,0,z,1};
	//	M4f_mul_M4f (&rot, &trans);
		
		tV4f p00 = {-d, -d, -d, 1};
		tV4f p01 = {+d, -d, -d, 1};
		tV4f p02 = {+d, +d, -d, 1};
		tV4f p03 = {-d, +d, -d, 1};
		
		tV4f p10 = {-d, -d,  d, 1};
		tV4f p11 = {+d, -d,  d, 1};
		tV4f p12 = {+d, +d,  d, 1};
		tV4f p13 = {-d, +d,  d, 1};
		
	/*	V4f_mul_M4f (&p00, &rot);
		V4f_mul_M4f (&p01, &rot);
		V4f_mul_M4f (&p02, &rot);
		V4f_mul_M4f (&p03, &rot);
		V4f_mul_M4f (&p10, &rot);
		V4f_mul_M4f (&p11, &rot);
		V4f_mul_M4f (&p12, &rot);
		V4f_mul_M4f (&p13, &rot);/**/
		
	//	V4f_Print (&p00);
	//	printf ("closer\n");
		V4f_DrawPosPos (pcam, &p00, &p01);
		V4f_DrawPosPos (pcam, &p01, &p02);
		V4f_DrawPosPos (pcam, &p02, &p03);
		V4f_DrawPosPos (pcam, &p03, &p00);
		
	//	printf ("mid\n");
		V4f_DrawPosPos (pcam, &p00, &p10);
		V4f_DrawPosPos (pcam, &p01, &p11);
		V4f_DrawPosPos (pcam, &p02, &p12);
		V4f_DrawPosPos (pcam, &p03, &p13);
		
	//	printf ("far\n");
		V4f_DrawPosPos (pcam, &p10, &p11);
		V4f_DrawPosPos (pcam, &p11, &p12);
		V4f_DrawPosPos (pcam, &p12, &p13);
		V4f_DrawPosPos (pcam, &p13, &p10);
		
	}/**/
	
	SDL_mutexP (gM.pGaze_Mutex);
	
/*	float ox, oy;
	ox = gM.Left.aScreen[0].View.Top.x; oy = gM.Left.aScreen[0].View.Top.y;
	for (y = oy-80; y < oy+100; ++y)
		for (x = ox-80; x < ox+100; ++x)
			dspix(x, y) = 0;
	
	ox = gM.Left.aScreen[0].View.Left.x; oy = gM.Left.aScreen[0].View.Left.y;
	for (y = oy-80; y < oy+100; ++y)
		for (x = ox-80; x < ox+80; ++x)
			dspix(x, y) = 0;
	
	ox = gM.Left.aScreen[0].View.Front.x; oy = gM.Left.aScreen[0].View.Front.y;
	for (y = oy-100; y < oy+100; ++y)
		for (x = ox-100; x < ox+100; ++x)
			dspix(x, y) = 0;
	/**/
	
	abort();//Head_Calc_M4_Rel (&gM.Head, &gM.HeadC);
	
/*	gM.GazeL = Head_Diff4_Eye (&gM.Head, &gM.Left);
	gM.GazeR = Head_Diff4_Eye (&gM.Head, &gM.Right);
	
	gM.L_Vec = gM.GazeL;
	gM.R_Vec = gM.GazeR;
	
	gM.GazeL = Eye_map_point(&gM.Left, gM.GazeL);
	gM.GazeR = Eye_map_point(&gM.Right, gM.GazeR);
	
	gM.Gaze.x = (gM.GazeL.x+gM.GazeR.x)/2.0f;
	gM.Gaze.y = (gM.GazeL.y+gM.GazeR.y)/2.0f;
	
	
	gM.GazeL = Head_Diff_Eye (&gM.Head, &gM.Left);
	gM.GazeR = Head_Diff_Eye (&gM.Head, &gM.Right);
	
	gM.Gaze.x = (gM.GazeL.x+gM.GazeR.x)/2.0f;
	gM.Gaze.y = (gM.GazeL.y+gM.GazeR.y)/2.0f;/**/
	
	
	Head_Draw (&gM.Head);
	
	if (0) {
		tM4f bproj = pcam->Proj;
		tM4f bworld = pcam->World;
		
		#undef dpoint3
		#define dpoint3(name,x0,y0,z0)	\
			tV4f name = {x0,y0,z0,1};//	V4f_mul_M4f (&name, &model)
		
		M4f_Iden (&pcam->Proj);
		M4f_Iden (&pcam->World);
		
	//	tV4f screen_c = {0, 10, 0, 0};
	//	tV4f head_c = {gM.Head.M4.x03, gM.Head.M4.x13, gM.Head.M4.x23, 0};
		tV4f head_c = {0, 10, 0, 0};
		tV4f screen_c = {gM.Head.M4.x03, gM.Head.M4.x13, gM.Head.M4.x23, 0};
		
		V4f_sub_V4f (&head_c, &screen_c);
		
		float ax = atan2(head_c.z,head_c.y), ay = atan2(head_c.z,head_c.x);
		
	//	M4f_rotx (&pcam->World, M_PI_2 - ax);
	//	M4f_roty (&pcam->World, -M_PI_2 + ay);
	//	head_c.x = 5;
	//	head_c.y = 2;
	//	pcam->World = gM.Head.M4I;
	//	pcam->World.x32 = 0;
	//	M4f_trans (&pcam->World, head_c.x/11, 0, -20);
		
	//	head_c.x = -head_c.x;
	//	head_c.y = -head_c.y;
	//	head_c.z = 1;
		
	//	head_c.x = gM.tmp.x;
	//	head_c.y = gM.tmp.y;
	//	head_c.z = gM.tmp.z;
		
		head_c.x = head_c.x/gM.tmp.x;
		head_c.y = -head_c.y/gM.tmp.y;
		head_c.z = head_c.z / gM.tmp.z;
		
		V4f_Print (&head_c);
		
		M4f_trans (&pcam->World, -head_c.x, -head_c.y, -head_c.z);
	//	model.x32 = gM.Head.M4.32;
	//	fov = atan2(V4f_dist(&head_c), 22/2);
	//	fov = atan2(head_c.z, 22/2);
		
		// Set up our view matrix. A view matrix can be defined given an eye point,
		// a point to lookat, and a direction for which way is up. Here, we set the
		// eye five units back along the z-axis and up three units, look at the
		// origin, and define "up" to be in the y-direction.
	/*	device.Transform.View = Matrix.LookAtLH(
			new Vector3(headX, headY, headDist),
			new Vector3(headX, headY, 0),
			new Vector3(0.0f, 1.0f, 0.0f));
		
		// For the projection matrix, we set up a perspective transform (which
		// transforms geometry from 3D view space to 2D viewport space, with
		// a perspective divide making objects smaller in the distance). To build
		// a perpsective transform, we need the field of view (1/4 pi is common),
		// the aspect ratio, and the near and far clipping planes (which define at
		// what distances geometry should be no longer be rendered).
		
		//compute the near plane so that the camera stays fixed to -.5f*screenAspect, .5f*screenAspect, -.5f,.5f
		//compting a closer plane rather than simply specifying xmin,xmax,ymin,ymax allows things to float in front of the display
		float nearPlane = .05f;
		device.Transform.Projection = Matrix.PerspectiveOffCenterLH(
				nearPlane*(-.5f * screenAspect + headX)/headDist, 
				nearPlane*(.5f * screenAspect + headX)/headDist, 
				nearPlane*(-.5f - headY)/headDist, 
				nearPlane*(.5f - headY)/headDist, 
				nearPlane, 100);
		*/
	//	head_c.x = 0;
	//	head_c.y = 0;
		pcam->Proj_L = -1 - head_c.x;
		pcam->Proj_R = 1 - head_c.x;
		pcam->Proj_B = -1 - head_c.y;
		pcam->Proj_T = 1 - head_c.y;/**/
		
		pcam->Proj_N = head_c.z;
	//	pcam->Proj_L /= -head_c.z;
	//	pcam->Proj_R /= -head_c.z;
	//	pcam->Proj_B /= -head_c.z;
	//	pcam->Proj_T /= -head_c.z;
		
	/*	double fov;
		fov = gM.Cam.Image_FOV;
		{
			double aspectRatio = (float)gM.Cam.Image_W/(float)gM.Cam.Image_H, front = 1, back = 10;
			aspectRatio = 1;
			double tangent = tan(fov/2 * deg2rad);   // tangent of half fovY
			double height = front * tangent;          // half height of near plane
			double width = height * aspectRatio;      // half width of near plane
			
			// params: left, right, bottom, top, near, far
			pcam->Proj_W = 2*width;
			pcam->Proj_H = 2*height;
		//	pcam->Proj_N = front;
		//	pcam->Proj_F = back;
		}/**/
	//	pcam->Proj_L = -pcam->Proj_W/2 - head_c.x/11;
	//	pcam->Proj_R = pcam->Proj_W/2 - head_c.x/11;
	//	pcam->Proj_B = -pcam->Proj_H/2 - head_c.y/11;
	//	pcam->Proj_T = pcam->Proj_H/2 - head_c.y/11;
		
	//	pcam->Proj_L /= head_c.z;
	//	pcam->Proj_R /= head_c.z;
	//	pcam->Proj_B /= head_c.z;
	//	pcam->Proj_T /= head_c.z;
		{
		//	float l = -pcam->Proj_W/4, r = 3*pcam->Proj_W/4;
		//	float t = pcam->Proj_W/2, b = -pcam->Proj_W/2;
			
			pcam->Proj.f[0][0] = 2*pcam->Proj_N / (pcam->Proj_R-pcam->Proj_L);
			pcam->Proj.f[1][1] = 2*pcam->Proj_N / (pcam->Proj_T-pcam->Proj_B);
			
			pcam->Proj.f[0][2] = (pcam->Proj_R+pcam->Proj_L) / (pcam->Proj_R-pcam->Proj_L);
			pcam->Proj.f[1][2] = (pcam->Proj_T+pcam->Proj_B) / (pcam->Proj_T-pcam->Proj_B);
			
			pcam->Proj.f[2][2] = -(pcam->Proj_F+pcam->Proj_N) / (pcam->Proj_F - pcam->Proj_N);
			
			pcam->Proj.f[2][3] = (-2*pcam->Proj_F*pcam->Proj_N) / (pcam->Proj_F - pcam->Proj_N);
			pcam->Proj.f[3][2] = -1;
		}/**/
		
	//	M4f_Print (&pcam->World);
	//	M4f_Print (&pcam->Proj);
		void draw (float x, float y, float z) {
			float dd = 0.2;
			dpoint3(pt, x+0,		y+-dd,	z);
			dpoint3(pb, x+0,		y+dd,		z);
			dpoint3(pl, x+-dd,	y+0,		z);
			dpoint3(pr, x+dd,		y+0,		z);
			V4f_DrawPosPos (pcam, &pl, &pt);
			V4f_DrawPosPos (pcam, &pt, &pr);
			V4f_DrawPosPos (pcam, &pr, &pb);
			V4f_DrawPosPos (pcam, &pb, &pl);
			
			dpoint3(pn, x,		y,		z);
			dpoint3(pf, x,		y,		-10);
			V4f_DrawPosPos (pcam, &pf, &pn);
		};
		
		draw (0,0,2);
		
		draw (0.5,0.5,-10);
		
		draw (-0.6,0.2,0);
		
		
		draw (-0.8,-0.8,-10);
		
		draw (0.5,-0.5,2);
		/**/
		si i;
		float z = 0;
		for (i = 0; i < 10; ++i, z -= 2) {
		//	printf ("aoaeuoaue\n");
			float xd = 1, dd = 1;
			dpoint3(pl, +xd,		-dd,		z);
			dpoint3(pt, +xd,		+dd,		z);
			dpoint3(pr, -xd,		+dd,		z);
			dpoint3(pb, -xd,		-dd,		z);
			
			V4f_DrawPosPos (pcam, &pl, &pt);
			V4f_DrawPosPos (pcam, &pt, &pr);
			V4f_DrawPosPos (pcam, &pr, &pb);
			V4f_DrawPosPos (pcam, &pb, &pl);
		};
		
		for (z = -1; z <= 1; z += 0.1f) {
			{
				dpoint3(pn, z,		1,		0);
				dpoint3(pf, z,		1,		-10);
				V4f_DrawPosPos (pcam, &pf, &pn);
			}
			{
				dpoint3(pn, z,		-1,		0);
				dpoint3(pf, z,		-1,		-10);
				V4f_DrawPosPos (pcam, &pf, &pn);
			}
			{
				dpoint3(pn, 1,		z,		0);
				dpoint3(pf, 1,		z,		-10);
				V4f_DrawPosPos (pcam, &pf, &pn);
			}
			{
				dpoint3(pn, -1,		z,		0);
				dpoint3(pf, -1,		z,		-10);
				V4f_DrawPosPos (pcam, &pf, &pn);
			}
		}/**/
		
		pcam->World = bworld;
		pcam->Proj = bproj;
	}
	
	
	Screen_Print (gM.aScreen + 0);
	Screen_Print (gM.aScreen + 1);
	Screen_Print (gM.aScreen + 2);
/*	Screen_Eye_Print (gM.aScreen + 1, &gM.Left);
	Screen_Eye_Print (gM.aScreen + 1, &gM.Right);
	Screen_Eye_Print (gM.aScreen + 2, &gM.Left);
	Screen_Eye_Print (gM.aScreen + 2, &gM.Right);/**/
	
/*	if (0) {
		tV2si pos;
		pos.x = gM.Head.DotL.P.x;
		pos.y = gM.Head.DotL.P.y;
		printf ("ehh pos	%f	%f\n", gM.Head.DotL.P.x, gM.Head.DotL.P.y);
		printf ("ehhh %f\n", NN_Sphere (&pos, 50, 55, 0, 0));
	}*/
	
//	dprojline3(0, 0, 1, 400.0, 300.0, 1);
	
	#if 0
	if (0) {
		Head_Calc_M_Rel (&gM.Head, &gM.HeadC);
		
	//	tV2f pel, per;
		
		gM.GazeL = Head_Diff_Eye (&gM.Head, &gM.Left);
		gM.GazeR = Head_Diff_Eye (&gM.Head, &gM.Right);
		
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
		
		V2f_mul_M3f (&gM.GazeL, &gM.Head.MI);
		V2f_mul_M3f (&gM.GazeR, &gM.Head.MI);
		
		gM.L_Vec = gM.GazeL;
		gM.R_Vec = gM.GazeR;
		
		gM.GazeL = Eye_map_point(&gM.Left, gM.GazeL);
		gM.GazeR = Eye_map_point(&gM.Right, gM.GazeR);
		
		gM.Gaze.x = (gM.GazeL.x+gM.GazeR.x)/2.0f;
		gM.Gaze.y = (gM.GazeL.y+gM.GazeR.y)/2.0f;
	}else {
		float x, y;
		
		
		if (gM.GazeMode == 0) {	//Normal head tracking
			Head_Eye_Gaze (&gM.Head, &gM.Left);
			Head_Eye_Gaze (&gM.Head, &gM.Right);
			
		}else if (gM.GazeMode == 1) {	//Glint normalization homography mapping
			gM.GazeL = Eye_map_point(&gM.Left, gM.Left.GV);
			
			gM.GazeL.x += gM.aScreen[0].Off.x;
			gM.GazeL.y += gM.aScreen[0].Off.y;
			
			gM.GazeR = Eye_map_point(&gM.Right, gM.Right.GV);
			gM.GazeR.x += gM.aScreen[0].Off.x;
			gM.GazeR.y += gM.aScreen[0].Off.y;
			
		//	gM.GazeR = gM.GazeL;
			
		//	printf ("Gaze						%f %f\n", gM.GazeL.x, gM.GazeL.y);
		}else if (gM.GazeMode == 2) {	//Glint sphere tracking
			Head_Eye_Gaze (&gM.Head, &gM.Left);
		//	Head_Eye_Gaze (&gM.Head, &gM.Right);
			gM.GazeR = gM.GazeL;
		//	Eye_Glint_Gaze (&gM.Left);
		//	Eye_Glint_Gaze (&gM.Right);
		}
		
		x = gM.GazeL.x;
		y = gM.GazeL.y;
		
		if (!isnan(x) && !isnan(y)) {
			if (!isnan(gM.GazeR.x) && !isnan(gM.GazeR.y)) {
				x = (gM.GazeL.x+gM.GazeR.x)/2.0f;
				y = (gM.GazeL.y+gM.GazeR.y)/2.0f;
			}
		}else {
			x = gM.GazeR.x;
			y = gM.GazeR.y;
		}
		
		if (!isnan(x) && !isnan(y)) {
			switch (gM.GazeAvg) {
			case 0:
				gM.Gaze.x = x;
				gM.Gaze.y = y;
				break;
			case 1:
				gM.Gaze.x = (gM.Gaze.x+x)/2.0f;
				gM.Gaze.y = (gM.Gaze.y+y)/2.0f;
				break;
			case 2: {
				#define dn (sizeof(weight)/sizeof(weight[0]))
				static float weight[] = {0.2, 0.3, 0.5};
				static tV2f prev[dn];
				si i;
				for (i = 0; i < dn-1; ++i) {
					prev[i] = prev[i+1];
				}
				prev[dn-1].x = x;
				prev[dn-1].y = y;
				x = 0;
				y = 0;
				for (i = 0; i < dn; ++i) {
					x += prev[i].x * weight[i];
					y += prev[i].y * weight[i];
				}
				prev[dn-1].x = x;
				prev[dn-1].y = y;
				gM.Gaze.x = x;
				gM.Gaze.y = y;
				#undef dn
				break;
			}
			case 3: {
				float a = sqrt(dpow2(gM.Gaze.x-x) + dpow2(gM.Gaze.y-y)) - gM.GazeAvg_3_MinDist;
				
				a /= gM.GazeAvg_3_Dist;
				if (a > 1)
					a = 1;
				if (a < gM.GazeAvg_3_MinAlpha)
					a = gM.GazeAvg_3_MinAlpha;
				
				gM.Gaze.x = (1-a)*gM.Gaze.x + a*x;
				gM.Gaze.y = (1-a)*gM.Gaze.y + a*y;
				break;
			}
			}
		}
	}
	#endif
//	point_clip (&gM.GazeL);
//	point_clip (&gM.GazeR);
//	point_clip (&gM.Gaze);
	
	gCol.R = 0xFF;
	gCol.G = 0xFF;
	gCol.B = 0x0;
//	Vec_Draw (gM.Head.DotL.P.x, gM.Head.DotL.P.y, gM.Head.DotL.P.x+gM.L_Vec.x, gM.Head.DotL.P.y+gM.L_Vec.y);
//	Vec_Draw (gM.Head.DotR.P.x, gM.Head.DotR.P.y, gM.Head.DotR.P.x+gM.R_Vec.x, gM.Head.DotR.P.y+gM.R_Vec.y);
	
	//Vec_Draw ((640/4*1)-gM.L_Vec.x, (480/2)-gM.L_Vec.y, (640/4*1), (480/2));
	//Vec_Draw ((640/4*3)-gM.R_Vec.x, (480/2)-gM.R_Vec.y, (640/4*3), (480/2));
	
	SDL_mutexV (gM.pGaze_Mutex);
	
//	muhaha_Cross (&gM.GazeL, (tPix){0xFF, 0x0, 0x0});		muhaha_Cross (&gM.GazeR, (tPix){0xFF, 0xF, 0xF});
//	muhaha_Cross (&gM.Gaze, (tPix){0xFF, 0xF, 0xF});
	
	
	
//	if (videoIn->formatIn != V4L2_PIX_FMT_YUYV)
//		printf ("videoIn->formatIn %d\n", videoIn->formatIn);
	
//	memcpy(gM.pDst, videoIn->framebuffer, videoIn->width * (videoIn->height) * 2);
	
	
	
//	XIWarpPointer (gM.X.pDisp, gM.X.Pointer.deviceid, None, None, 0, 0, 0, 0, 1, 1);
	
//	printf ("gM.Gaze %f\t%f\n", gM.Gaze.x, gM.Gaze.y);
	if (!gM.bGazeHold) {
		if (gM.Micro.State) {
			gM.Gaze = gM.Micro.Gaze;
			float dx, dy;
			dx = gM.Head.R_Y - gM.Micro.R_Y;
			dy = gM.Head.R_X - gM.Micro.R_X;
			dx *= gM.Micro.SX;
			dy *= gM.Micro.SY;
		//	printf ("micro dxy %f %f\n", dx, dy);
			gM.Gaze.x += dx;
			gM.Gaze.y += dy;
		}
		switch (gM.Pointer_Mode) {
		case 1:{
			Marker_Move (&gM.CrossHair, gM.Gaze.x, gM.Gaze.y);
		//	Marker_Move (&gM.CrossHair, gM.aScreen[0].Off.x + gM.Gaze.x, gM.aScreen[0].Off.y + gM.Gaze.y);
			break;
		}
		case 2: {
			XIWarpPointer (gM.X.pDisp, gM.X.Pointer.deviceid, 
				None, gM.aScreen[0].Win,
				0, 0, 0, 0,
				gM.Gaze.x,
				gM.Gaze.y
			);
		/*	XIWarpPointer (gM.X.pDisp, gM.X.Pointer.deviceid, 
				None, gM.aScreen[0].Win,
				0, 0, 0, 0,
				gM.aScreen[0].Off.x + gM.Gaze.x,
				gM.aScreen[0].Off.y + gM.Gaze.y
			);/**/
			XFlush(gM.X.pDisp);
			break;
		}
		}
	}
	
	if (gM.Cal_bRecord) {
		float z = gM.aScreen[gM.Screen_CalIdx].C.z;
		if (0) {
			tV4f e00, e01, e10, e11;
			tV4f p0, p1;
			
			Head_Eye_VectorGlob (&gM.Head, &gM.Left, &e00, &e01, NULL);
			Head_Eye_VectorGlob (&gM.Head, &gM.Right, &e10, &e11, NULL);
			
			V4f_Intersect_Line01_Line01 (
				&p0, &p1,
				&e00, &e01,
				&e10, &e11
			);
			z = (p0.z + p1.z) / 2;
			printf ("Screen_Cal_Save Z %f\n", z);
		}
		
		Screen_Cal_Record (gM.aScreen + gM.Screen_CalIdx, &gM.Left, z);
		Screen_Cal_Record (gM.aScreen + gM.Screen_CalIdx, &gM.Right, z);
		
		++gM.Cal_Record_Num;
	}
/*	{
		printf ("crap\n");
		XEvent e;
		XNextEvent(gM.X.pDisp, &e);
		printf ("yeah\n");
	}/**/
	SDL_mutexV (gM.mutMainIsOut);
	
	
//	SDL_mutexP (gM.mutMainIsIn);
//	printf ("muhahaha		mutMainIsIn	E\n");
	
	SDL_mutexP (gM.mutMainIsIn);
}

void	muhaha_Loop	()
{
	
	for (int i = 0; i < gM.Cam_N; ++i) {
		tCam* pcam = &gM.aCam[i];
		
		//SDL_LockSurface (pcam->SDL_Surf);
		
		pcam->SDL_Thread = SDL_CreateThread (Cam_Loop, "Cam Thread", (void *)pcam);
	}/**/
	
	for (si i = 0; i < M_Dbg_NUM; ++i) {
		tDbg* pdbg = &gM.aDbg[i];
		SDL_LockSurface (pdbg->SDL_Surf);
	}
	
	SDL_SemPostN (gM.sWaitForUpdate, gM.Cam_N);
	/* main big loop */
	while (gM.aCam[0].UVC->signalquit) {
		// Measure the frame rate every (fps/2) frames
	/*	if(loop_counter ++ % frmrate_update == 0) {
			currtime = SDL_GetTicks();	// [ms]
			if (currtime - lasttime > 0) {
				frmrate = frmrate_update * (1000.0 / (currtime - lasttime));
			}
			lasttime = currtime;
		}*/
		
	/*	for (int i = 0; i < gM.Cam_N; ++i) {
			tCam* pcam = &gM.aCam[i];
			
		}*/
		
		//SDL_LockMutex(affmutex);
		//ptdata.frmrate = frmrate;
		//SDL_WM_SetCaption(videoIn->status, NULL);
		//SDL_UnlockMutex(affmutex);
	//	SDL_Delay(10);
		
		SDL_SemWaitN (gM.sWaitForCams, gM.Cam_N);
		//printf ("Main update start  ");
		
		for (si i = 0; i < M_Dbg_NUM; ++i) {
			tDbg* pdbg = &gM.aDbg[i];
			SDL_UnlockSurface (pdbg->SDL_Surf);
			SDL_UpdateWindowSurface (pdbg->SDL_Win);
		}
		
		for (si i = 0; i < gM.Cam_N; ++i) {
			tCam* pcam = &gM.aCam[i];
			
			SDL_UpdateWindowSurface (pcam->SDL_Win);
		}
		for (si i = 0; i < gM.Cam_N; ++i) {
			tCam* pcam = &gM.aCam[i];
			if (uvcGrab(pcam->UVC) < 0) {
				printf("Error grabbing\n");
				abort();
			}
		}
		
		for (si i = 0; i < M_Dbg_NUM; ++i) {
			tDbg* pdbg = &gM.aDbg[i];
			SDL_LockSurface (pdbg->SDL_Surf);
			memset(pdbg->SDL_Surf->pixels, 0, pdbg->SDL_Surf->w * pdbg->SDL_Surf->h * sizeof(tRGBA));
		}
		
		//for (si i = 0; i < gM.Cam_N; ++i) {
			//tCam* pcam = &gM.aCam[i];
		//}
	
		//printf ("end\n");
		
		SDL_SemPostN (gM.sWaitForUpdate, gM.Cam_N);
		
	}
}

int muhaha_eventThread(void *data)
{
	struct pt_data *gdata = (struct pt_data *) data;
	struct v4l2_control control;
	SDL_Surface *pscreen = *gdata->ptscreen;
	//struct vdIn *videoIn = gdata->ptvideoIn;
	SDL_Event *sdlevent = gdata->ptsdlevent;
	SDL_Rect *drect = gdata->drect;
	SDL_mutex *affmutex = gdata->affmutex;
	
	static int x, y;
	static tCam* pcam = NULL;
	int mouseon = 0;
	int value = 0;
	int len = 0;
//	short incpantilt = INCPANTILT;
	int boucle = 0;
	
	static si ii = 0;
	
	si last_ms = SDL_GetTicks();
//	action_gui curr_action = A_VIDEO;
	while (gM.aCam[0].UVC->signalquit) {
		SDL_LockMutex(affmutex);
		
		tV2f pos;
		if (pcam) {
			pos.x = x*pcam->Image_W / pcam->SDL_Surf->w;
			pos.y = y*pcam->Image_H / pcam->SDL_Surf->h;
		}
		float frmrate = gdata->frmrate;
		
		while (SDL_PollEvent(sdlevent)) {	//scan the event queue
			switch (sdlevent->type) {
			case SDL_KEYDOWN:
				printf("Key down\n");
				switch (sdlevent->key.keysym.sym) {
			/*	case SDLK_BACKSPACE:
					printf("Reset points\n");
					ii = 0;
					break;/**/
				case SDLK_s:
					printf("gM.DeSat\n");
					gM.DeSat = !gM.DeSat;
					break;
				case SDLK_o:
					if (pcam) {
						gM.Head.DotC.aCam[pcam->Idx].P = pos;
						gM.Head.DotC.aCam[pcam->Idx].V.x = 0;
						gM.Head.DotC.aCam[pcam->Idx].V.y = 0;
					}
					break;
				case SDLK_a:
					if (pcam) {
						gM.Head.DotL.aCam[pcam->Idx].P = pos;
						gM.Head.DotL.aCam[pcam->Idx].V.x = 0;
						gM.Head.DotL.aCam[pcam->Idx].V.y = 0;
					}
					break;
				case SDLK_e:
					if (pcam) {
						gM.Head.DotR.aCam[pcam->Idx].P = pos;
						gM.Head.DotR.aCam[pcam->Idx].V.x = 0;
						gM.Head.DotR.aCam[pcam->Idx].V.y = 0;
					}
					break;
				case SDLK_SEMICOLON:
					gM.Left.P.x = x;
					gM.Left.P.y = y;
					break;
				case SDLK_j:
					gM.Right.P.x = x;
					gM.Right.P.y = y;
					break;
				case SDLK_y:
					gM.Head.aPoint[0].P.x = x;
					gM.Head.aPoint[0].P.y = y;
					break;/**/
			//	case SDLK_SPACE:
			//		gM.HeadC = gM.Head;
			//		break;
				
			/*	case SDLK_BACKSPACE:
					gM.Left.InHead.Line_N = 0;
					gM.Right.InHead.Line_N = 0;
					break;
				case SDLK_BACKSLASH: {
					Head_Eye_LineAdd (&gM.Head, &gM.Left);
					Head_Eye_LineAdd (&gM.Head, &gM.Right);
					break;
				}
				case SDLK_1:
					gM.Screen_CalIdx = 0;
					Screen_Cal_Init (gM.aScreen + gM.Screen_CalIdx);
					break;
				case SDLK_2:
					gM.Screen_CalIdx = 1;
					Screen_Cal_Init (gM.aScreen + gM.Screen_CalIdx);
					break;
				case SDLK_3:
					gM.Screen_CalIdx = 2;
					Screen_Cal_Init (gM.aScreen + gM.Screen_CalIdx);
					break;
				case SDLK_SPACE:
					Screen_Cal_Save (gM.aScreen + gM.Screen_CalIdx);
					break;
				
				case SDLK_RETURN: {
					switch (gM.Pointer_Mode) {
					case 0:
						++gM.Pointer_Mode;
						Marker_Show (&gM.CrossHair);
						break;
					case 1:
						++gM.Pointer_Mode;
					//	gM.Pointer_Mode = 0;
						Marker_Hide (&gM.CrossHair);
						break;
					case 2:
						gM.Pointer_Mode = 0;
						break;
					}
					break;
				}/**/
				
				case SDLK_ESCAPE: {
					pcam->UVC->signalquit = 0;
					break;
				}
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
			//	printf ("mouse down %d %d\n", x, y);
				switch (sdlevent->button.button) {
				case SDL_BUTTON_LEFT: {
					SDL_mutexP (gM.pGaze_Mutex);
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
					printf ("	left %f %f\n", gM.Left.Homo.vectors[ii].x, gM.Left.Homo.vectors[ii].y);
					printf ("	righ %f %f\n", gM.Right.Homo.vectors[ii].x, gM.Right.Homo.vectors[ii].y);
					
					if (++ii >= CALIBRATIONPOINTS) {
						ii = 0;
						CalculateCalibration(&gM.Left.Homo);
						CalculateCalibration(&gM.Right.Homo);
					}
					SDL_mutexV (gM.pGaze_Mutex);
					break;
				}
				case SDL_BUTTON_RIGHT: {
				//	gM.Left.P.x = x;
				//	gM.Left.P.y = y;
					ui iy, ix;
					ix = x / (gM.pScreen->w/dEye_Screen_Cal_NUM);
					iy = y / (gM.pScreen->h/dEye_Screen_Cal_NUM);
					
					gM.Left.aScreen[0].aaCal[iy][ix].SX = x;
					gM.Left.aScreen[0].aaCal[iy][ix].SY = y;
					
				//	gM.aScreen[0].C.x
				//	gM.Left.aScreen[0].aaCal[iy][ix].P = 
				//	gM.Left.aScreen[0].aaCal[iy][ix].Y = y;
					gM.Right.Homo.vectors[ii].x = gM.L_Vec.x;
					gM.Right.Homo.vectors[ii].y = gM.L_Vec.y;/**/
					printf ("got point %ld\n", ii);
					printf ("	left %f %f\n", gM.Left.Homo.vectors[ii].x, gM.Left.Homo.vectors[ii].y);
					printf ("	righ %f %f\n", gM.Right.Homo.vectors[ii].x, gM.Right.Homo.vectors[ii].y);
					
					if (++ii >= CALIBRATIONPOINTS) {
						ii = 0;
						CalculateCalibration(&gM.Left.Homo);
						CalculateCalibration(&gM.Right.Homo);
					}
				//	SDL_mutexV (gM.pGaze_Mutex);
					break;
				}
				}
				break;
			case SDL_MOUSEMOTION:
				//curr_action = GUI_whichbutton(x, y, pscreen, videoIn);
				x = sdlevent->motion.x;
				y = sdlevent->motion.y;
				pcam = NULL;
				for (int i = 0; i < gM.Cam_N; ++i) {
					if (sdlevent->motion.windowID == SDL_GetWindowID(gM.aCam[i].SDL_Win)) {
						pcam = &gM.aCam[i];
						break;
					}
				}
				break;
		/*	case SDL_VIDEORESIZE:
			/*	pscreen =
				      SDL_SetVideoMode(sdlevent->resize.w,
				                       sdlevent->resize.h, 0,
				                       SDL_VIDEO_Flags);
				drect->w = sdlevent->resize.w;
				drect->h = sdlevent->resize.h;*/
				break;
			case SDL_QUIT:
				printf("\nQuit signal received.\n");
				pcam->UVC->signalquit = 0;
				break;
			}
		}			//end if poll
		SDL_UnlockMutex(affmutex);
		/* traiter les actions */
		value = 0;
		
	//	len = 100 * sizeof(char);	// as allocated in init_videoIn
	//	snprintf(videoIn->status, len, "%s, %.1f fps", title_act[curr_action].title, frmrate);
		
		SDL_Delay(50);
		
		si ms = SDL_GetTicks();
		if (ms - last_ms >= 1000) {
			fps = tfps;
			printf ("FPS %ld\n", fps);
			tfps = 0;
			last_ms = SDL_GetTicks();
		}
	//	printf("fp/s %d\n",frmrate);
	}
	
	/* Close the stream capture file */
	if (pcam->UVC->captureFile) {
		fclose(pcam->UVC->captureFile);
		printf("Stopped raw stream capturing to stream.raw. %u bytes written for %u frames.\n",
		       pcam->UVC->bytesWritten, pcam->UVC->framesWritten);
	}
	/* Display stats for raw frame stream capturing */
	if (pcam->UVC->rawFrameCapture == 2) {
		printf("Stopped raw frame stream capturing. %u bytes written for %u frames.\n",
		       pcam->UVC->rfsBytesWritten, pcam->UVC->rfsFramesWritten);
	}
}



/*
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
/**/

#define SIGN(u, v)     ( (v)>=0.0 ? fabs(u) : -fabs(u) )
//#define MAX(x, y)     ( (x) >= (y) ? (x) : (y) )  

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




void	Ehh_Draw_Line_2d	(tCam* pcam, si x0, si y0, si x1, si y1)
{
	if (x0 < -2*pcam->SDL_Surf->w || x0 > 3*pcam->SDL_Surf->w
		|| x1  < -2*pcam->SDL_Surf->w || x1 > 3*pcam->SDL_Surf->w
		|| y0 < -2*pcam->SDL_Surf->h || y0 > 3*pcam->SDL_Surf->h
		|| y1  < -2*pcam->SDL_Surf->h || y1 > 3*pcam->SDL_Surf->h
	)
		return;
	//x0 /= 2;
	//y0 /= 2;
	//x1 /= 2;
	//y1 /= 2;
	si dx = x1 - x0;
	si dy = y1 - y0;
	
	ui dst_w = pcam->SDL_Surf->w;
	ui dst_h = pcam->SDL_Surf->h;
	tRGBA *pdst = (tRGBA*)pcam->SDL_Surf->pixels;
	
	tRGBA col = gCol;
	//tRGBA col;	col.R = 0xFF;	col.G = 0xFF;	col.B = 0xFF;	col.A = 0xFF;
	
	
	if (dx == 0 && dy == 0) {
		if ( (x0 >= 0) && (x0 < dst_w) && (y0 >= 0) && (y0 < dst_h) )
			*(pdst + x0 + (y0 * dst_w)) = col;
		return;
	}
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


void	S_Draw_Line_2d	(SDL_Surface* surf, si x0, si y0, si x1, si y1)
{
//	printf ("xy0 %ld\t%ld\t%ld\t%ld\n", x0, y0, x1, y1);
	si dx = x1 - x0;
	si dy = y1 - y0;
	
	ui dst_w = surf->pitch/4;
	ui dst_h = surf->h;
	tRGBA *pdst = surf->pixels;
	
	tRGBA col = gCol;
	
	if (dx == 0 && dy == 0) {
		if ( (x0 >= 0) && (x0 < dst_w) && (y0 >= 0) && (y0 < dst_h) )
			*(pdst + x0 + (y0 * dst_w)) = col;
		return;
	}
	//#define dspix(_x,_y) (*((u32*)surf->pixels + ((si)(_x) + (si)(_y)*surf->pitch/4)))
	
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


