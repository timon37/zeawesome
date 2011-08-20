
#include "muhaha.h"


void	S_Draw_Line_2d	(si x0, si y0, si x1, si y1);
void	Ehh_Draw_Line_2d	(si x0, si y0, si x1, si y1);

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
u32 gColARGB;

tM gM =
{
	.Y_Level = 128,
	.DeSat = 0,
	.Head = {
		.DotC = {
			.P = {400, 280},
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
			.P = {200, 300},
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
			.P = {600, 300},
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
		.P = {280, 240},
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
		.P = {360, 240},
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

#define deg2rad	(M_PI/180.0f)
#define rad2deg	(180.0f/M_PI)
#define dpow2(_num) ((_num)*(_num))

#define ddot(x0,y0,x1,y1) ((x0)*(x1) + (y0)*(y1))

#define ddist2(x0,y0,x1,y1)	(((x1)-(x0))*((x1)-(x0)) + ((y1)-(y0))*((y1)-(y0)))

#define ddist(x0,y0,x1,y1) sqrt(ddist2(x0,y0,x1,y1))

#define dopix(_x,_y) ((tPix*)videoIn->framebuffer + ((si)(_x) + (si)(_y)*videoIn->width))
#define dnpix(_x,_y) ((tPix*)gM.pDst + ((si)(_x) + (si)(_y)*videoIn->width))
#define dpix(_x,_y) dnpix(_x,_y)


#define dspix(_x,_y) (*((u32*)gM.pScreen->pixels + ((si)(_x) + (si)(_y)*gM.pScreen->pitch/4)))

#define dmono2rgb(_g) ((u32)(_g) | (u32)(_g)<<8 | (u32)(_g)<<16)


#define dpixout(_x,_y) ((_x) < 0 || (_y) < 0 || (_x) >= videoIn->width || (_y) >= videoIn->height)


#define dset_cyuv(_x,_y,y,u,v)	\
	do {	\
		dnpix(_x,_y)->Y = y;	\
		dnpix(_x,_y)->U = u;	\
		dnpix(_x,_y)->V = v;	\
	}while(0)
#define dset_c1(_x,_y) dset_cyuv (_x,_y,0xFF,0xF,0xF)

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


float	V2f_dot_V2f	(tV2f* pv0, tV2f* pv1)
{
	return pv0->x*pv1->x + pv0->y*pv1->y;
}
float	V4f_dot_V4f	(tV4f* pv0, tV4f* pv1)
{
	return pv0->x*pv1->x + pv0->y*pv1->y + pv0->z*pv1->z;
}

void	V4f_mul_S	(tV4f* pv0, float s)
{
	pv0->x *= s;
	pv0->y *= s;
	pv0->z *= s;
}


float	V2f_cross	(tV2f* pv0, tV2f* pv1)
{
	return pv0->x*pv1->y - pv0->y*pv1->x;
}

void	V4f_cross	(tV4f* pret, tV4f* pv0, tV4f* pv1)
{
	tV4f ret;
	ret.x = pv0->y*pv1->z - pv0->z*pv1->y;
	ret.y = pv0->z*pv1->x - pv0->x*pv1->z;
	ret.z = pv0->x*pv1->y - pv0->y*pv1->x;
	ret.w = 1;
	*pret = ret;
}

float	V2f_dist	(tV2f* pv)
{
	return sqrtf(pv->x*pv->x + pv->y*pv->y);
}

float	V4f_dist2	(tV4f* pv)
{
	return pv->x*pv->x + pv->y*pv->y + pv->z*pv->z;
}
float	V4f_dist	(tV4f* pv)
{
	return sqrtf(V4f_dist2(pv));
}


void	V4f_norm	(tV4f* pv0)
{
	V4f_mul_S (pv0, V4f_dist (pv0));
}

float	V4f_dist_V4f	(tV4f* pv0, tV4f* pv1)
{
	return sqrtf(dpow2(pv1->x-pv0->x) + dpow2(pv1->y-pv0->y) + dpow2(pv1->z-pv0->z));
}


void	V2f_sub_V2f		(tV2f* pv0, tV2f* pv1)
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


float	angle_norm_0_2pi	(float a)
{
	while (a < 0) {
		a += 2*M_PI;
	}
	while (a > 2*M_PI) {
		a -= 2*M_PI;
	}
	return a;
}

float	angle_norm_pi_pi	(float a)
{
	while (a < -M_PI) {
		a += 2*M_PI;
	}
	while (a > M_PI) {
		a -= 2*M_PI;
	}
	return a;
}

float	angle_diff_norm_pi_pi	(float a0, float a1)
{
	return angle_norm_pi_pi(a0 - a1);
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

void	V4f_DrawPosPos	(tV4f* ppos0, tV4f* ppos1)
{
	tV4f p0 = *ppos0, p1 = *ppos1;
	
//	V4f_mul_M4f (&p0, &gM.World);	V4f_mul_M4f (&p1, &gM.World);
//	V4f_mul_M4f (&p0, &gM.Proj);	V4f_mul_M4f (&p1, &gM.Proj);
	
	M4f_mul_V4f (&gM.World, &p0);		M4f_mul_V4f (&gM.World, &p1);
	M4f_mul_V4f (&gM.Proj, &p0);		M4f_mul_V4f (&gM.Proj, &p1);
	
	p0.x /= p0.w;
	p0.y /= p0.w;
	p0.z /= p0.w;
	
	p1.x /= p1.w;
	p1.y /= p1.w;
	p1.z /= p1.w;
	
	p0.x *= gM.View_W/2;	p0.x += gM.View_W/2;
	p0.y *= gM.View_H/2;	p0.y += gM.View_H/2;
	p1.x *= gM.View_W/2;	p1.x += gM.View_W/2;
	p1.y *= gM.View_H/2;	p1.y += gM.View_H/2;
	
//	printf ("p0 %f %f  p1 %f %f\n", p0.x, p0.y, p1.x, p1.y);
	Vec_Draw (p0.x, p0.y, p1.x, p1.y);
}

void	V4f_ScreenPosNorm	(tV4f* ppos0, tV2f* pret)
{
	tV4f p0 = *ppos0;
	
	M4f_mul_V4f (&gM.World, &p0);
	M4f_mul_V4f (&gM.Proj, &p0);
	
	p0.x /= p0.w;
	p0.y /= p0.w;
	p0.z /= p0.w;
	
	pret->x = p0.x;	pret->y = p0.y;
}
void	V4f_ScreenPosf	(tV4f* ppos0, tV2f* pret)
{
	V4f_ScreenPosNorm (ppos0, pret);
	
	pret->x *= gM.View_W/2;	pret->x += gM.View_W/2;
	pret->y *= gM.View_H/2;	pret->y += gM.View_H/2;
}
void	V4f_ScreenPos	(tV4f* ppos0, tV2si* pret)
{
	tV4f p0 = *ppos0;
	
	M4f_mul_V4f (&gM.World, &p0);
	M4f_mul_V4f (&gM.Proj, &p0);
	
	p0.x /= p0.w;
	p0.y /= p0.w;
	p0.z /= p0.w;
	
	p0.x *= gM.View_W/2;	p0.x += gM.View_W/2;
	p0.y *= gM.View_H/2;	p0.y += gM.View_H/2;
	
	pret->x = p0.x;	pret->y = p0.y;
}


void	Proj_Cam		()
{
	{
		float fd = dpow2(gM.Cam.Full_W) + dpow2(gM.Cam.Full_H);
		float id = dpow2(gM.Cam.Image_Zoom*gM.Cam.Image_W) + dpow2(gM.Cam.Image_Zoom*gM.Cam.Image_H);
	//	float id = ddist(gM.Cam.Image_W, gM.Cam.Image_H);
		gM.Cam.Image_FOV = id/fd * gM.Cam.Full_FOV;
		printf ("Image_FOV %f\n", gM.Cam.Image_FOV);
	}
	
	float fov = gM.Cam.Image_FOV*deg2rad;
	float a = atan2(gM.Cam.Image_H, gM.Cam.Image_W);
	
	gM.Cam.Image_FOV_W = fov*cos(a);
	gM.Cam.Image_FOV_H = fov*sin(a);
	
	gM.Proj_W = gM.Proj_N*tan(gM.Cam.Image_FOV_W/2.0f);
	gM.Proj_H = gM.Proj_N*tan(gM.Cam.Image_FOV_H/2.0f);
	
	gM.Proj_L = -gM.Proj_W;
	gM.Proj_R = gM.Proj_W;
	gM.Proj_B = -gM.Proj_H;
	gM.Proj_T = gM.Proj_H;
	
	gM.Proj_W *= 2;
	gM.Proj_H *= 2;
	
	gM.Proj.x00 = gM.Proj_N / (gM.Proj_W/2);
	gM.Proj.x11 = gM.Proj_N / (gM.Proj_H/2);
	
	gM.Proj.x22 = -(gM.Proj_F+gM.Proj_N) / (gM.Proj_F - gM.Proj_N);
	
	gM.Proj.x23 = (-2*gM.Proj_F*gM.Proj_N) / (gM.Proj_F - gM.Proj_N);
	gM.Proj.x32 = -1;
	
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
	}
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
//	XMapWindow(dpy, w);
	
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





void	Eye_Init	(tEye* peye)
{
	peye->Point_N = 0;
	peye->Point_Max = 0;
	peye->paPoint = 0;
	
	peye->InHead.Line_N = 0;
	memset (peye->InHead.aLine, 0, sizeof(peye->InHead.aLine));
}


void	Eye_CopyParam	(tEye* pdst, tEye* psrc)
{
	pdst->P = psrc->P;
	pdst->Ax = psrc->Ax;
	pdst->Ay = psrc->Ay;
	pdst->Aa = psrc->Aa;
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
		float t1 = atan2(y1 - peye->P.y, x1 - peye->P.x);
		
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
		float t0 = atan2(y0 - peye->P.y, x0 - peye->P.x);
		
		si j = i - 1;
		while (1) {
			float x1 = peye->paPoint[j].x, y1 = peye->paPoint[j].y;
			float t1 = atan2(y1 - peye->P.y, x1 - peye->P.x);
			
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
	float t0 = atan2(y0 - peye->P.y, x0 - peye->P.x);
	
	struct {
		si idx;
		float diff_t;
	}min = {-1, 4*deg2rad};
//	printf ("Eye_S4_Fit_FindRot\n");
	si i;
	for (i = 0; i < peye->Point_N; ++i) {
		float x1 = peye->paPoint[i].x, y1 = peye->paPoint[i].y;
		float t1 = atan2(y1 - peye->P.y, x1 - peye->P.x);
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
	float t0 = atan2(y0 - peye->P.y, x0 - peye->P.x);
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
		t1 = atan2(y1 - peye->P.y, x1 - peye->P.x);
		
		if (tgt > t1) {
			l = i;
		}else {
			r = i;
		}
	}
	i = l;
	x1 = peye->paPoint[i].x, y1 = peye->paPoint[i].y;
	t1 = atan2(y1 - peye->P.y, x1 - peye->P.x);
//	printf ("got l %ld\t%f\n", l, t1);
	float diff_t = fabsf(angle_diff_norm_pi_pi(t1, t0) - ang);
	if (diff_t < min.diff_t) {
		min.idx = i;
		min.diff_t = diff_t;
	}
	i = r;
	x1 = peye->paPoint[i].x, y1 = peye->paPoint[i].y;
	t1 = atan2(y1 - peye->P.y, x1 - peye->P.x);
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
		dx = x0 - peye->P.x;
		dy = y0 - peye->P.y;
		t0 = atan2(dy,dx);
		
		ox0 = (peye->Ax * cos(t0) * cos(peye->Aa) - peye->Ay * sin(t0) * sin(peye->Aa));
		oy0 = (peye->Ax * cos(t0) * sin(peye->Aa) + peye->Ay * sin(t0) * cos(peye->Aa));
		
		{
			si idx = Eye_Points_FindRot (peye, i, M_PI);
			if (idx >= 0) {
				float x1 = peye->paPoint[idx].x, y1 = peye->paPoint[idx].y;
				float t1 = atan2(y1 - peye->P.y, x1 - peye->P.x);
				
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
			}else {
				dnpix(x0,y0)->Y = 0xFF;
				dnpix(x0,y0)->U = 0xF;
				dnpix(x0,y0)->V = 0xF;
				float diff = sqrt(dx*dx+dy*dy) - sqrt(ox0*ox0+oy0*oy0);
				if (fabs(diff) >= 0.01f) {
				//	printf ("f1\n");
					peye->Ax += peye->Fit_Scale*diff*fabs(cos(t0));
					peye->Ay += peye->Fit_Scale*diff*fabs(sin(t0));
					
					Eye_Xset (peye, peye->P.x + peye->Fit_Trans*(diff )*cos(t0));
					Eye_Yset (peye, peye->P.y + peye->Fit_Trans*(diff )*sin(t0));
				}
			}
		}
	}/**/
	float avgerr = 0;
	for (i = 0; i < peye->Point_N; ++i) {
		float x = peye->paPoint[i].x, y = peye->paPoint[i].y;
		float t, dx, dy, xx, yy;
		dx = x - peye->P.x;
		dy = y - peye->P.y;
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
		dx = x0 - peye->P.x;
		dy = y0 - peye->P.y;
		t0 = atan2(dy,dx);
		
		ox0 = (peye->Ax * cos(t0) * cos(peye->Aa) - peye->Ay * sin(t0) * sin(peye->Aa));
		oy0 = (peye->Ax * cos(t0) * sin(peye->Aa) + peye->Ay * sin(t0) * cos(peye->Aa));
		
		tV2f p0 = {x0, y0}, op0 = {ox0, oy0};
		{
			si i1 = Eye_Points_FindRot (peye, i, M_PI_4/2);
			si i2 = Eye_Points_FindRot (peye, i, -M_PI_4/2);
			
			if (i1 >= 0 && i2 >= 0) {
				float x1 = peye->paPoint[i1].x, y1 = peye->paPoint[i1].y;
				float t1 = atan2(y1 - peye->P.y, x1 - peye->P.x);
				
				float x2 = peye->paPoint[i2].x, y2 = peye->paPoint[i2].y;
				float t2 = atan2(y2 - peye->P.y, x2 - peye->P.x);
				
				
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
				
				float pdx = (x0+x1)*0.5f - (2.0f*peye->P.x + ox0+ox1)*0.5f;
				float pdy = (y0+y1)*0.5f - (2.0f*peye->P.y + oy0+oy1)*0.5f;
				
				Eye_Xset (peye, peye->P.x + 0.10f * pdx );// *cos(t));
				Eye_Yset (peye, peye->P.y + 0.10f * pdy );// *sin(t));
				
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
					float t1 = atan2(y1 - peye->P.y, x1 - peye->P.x);
					
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
					continue;
				}
			}
			if (1) {
				float diff = sqrt(dx*dx+dy*dy) - sqrt(ox0*ox0+oy0*oy0);
				if (fabs(diff) >= 0.01f) {
				//	printf ("f1\n");
					peye->Ax += 0.1f*diff*fabs(cos(t0));
					peye->Ay += 0.1f*diff*fabs(sin(t0));
					
					Eye_Xset (peye, peye->P.x + 0.1f*(diff)*cos(t0));
					Eye_Yset (peye, peye->P.y + 0.1f*(diff)*sin(t0));
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
		dx = x0 - peye->P.x;
		dy = y0 - peye->P.y;
		t0 = atan2(dy,dx);
		
		ox0 = (peye->Ax * cos(t0) * cos(peye->Aa) - peye->Ay * sin(t0) * sin(peye->Aa));
		oy0 = (peye->Ax * cos(t0) * sin(peye->Aa) + peye->Ay * sin(t0) * cos(peye->Aa));
		
		tV2f p0 = {x0, y0}, op0 = {ox0, oy0};
		{
			if (1) {
				si idx = Eye_Points_FindRot (peye, i, M_PI);
				if (idx >= 0) {
					float x1 = peye->paPoint[idx].x, y1 = peye->paPoint[idx].y;
					float t1 = atan2(y1 - peye->P.y, x1 - peye->P.x);
					
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
					
					Eye_Xset (peye, peye->P.x + 0.1f*(diff)*cos(t0));
					Eye_Yset (peye, peye->P.y + 0.1f*(diff)*sin(t0));
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

void	Eye_OldFF		(tEye* peye)
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
	#define dd 12
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



si	Eye_S0_EdgeMark	(tEye* peye, si i, float t, float* px, float* py)
{
//	si border = peye->Exp_R*peye->Exp_R*3;
	u08 wrote = 0;
	float x = peye->P.x, y = peye->P.y;
	si n = 0;
	
	for (; ddist2(x,y,peye->P.x,peye->P.y) < peye->Exp_R*peye->Exp_R; ) {
		if (peye->LinView.x != 0)
			dspix(peye->LinView.x+n, peye->LinView.y+i) = dnpix(x,y)->Y | dnpix(x,y)->Y<<8 | dnpix(x,y)->Y<<16;
		
		if (	dopix(x,y)->Y > ay + 50
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
			dx = x0 - peye->P.x;
			dy = y0 - peye->P.y;
			t = atan2(dy,dx)/* - peye->Aa*/;
			
			xx = (peye->Ax * cos(t) * cos(peye->Aa) - peye->Ay * sin(t) * sin(peye->Aa));
			yy = (peye->Ax * cos(t) * sin(peye->Aa) + peye->Ay * sin(t) * cos(peye->Aa));
			
			float diff = sqrt(dx*dx+dy*dy) - sqrt(xx*xx+yy*yy);
			if (fabs(diff) >= 0.1f) {
			//	printf ("f1\n");
			//	peye->Ax += peye->Fit_Scale*diff*fabs(cos(t));
			//	peye->Ay += peye->Fit_Scale*diff*fabs(sin(t));
				
				Eye_Xset (peye, peye->P.x + 1*peye->Fit_Trans*(diff )*cos(t));
				Eye_Yset (peye, peye->P.y + 1*peye->Fit_Trans*(diff )*sin(t));
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

void	Eye_S0		(tEye* peye)
{
	tV2f oldpos = peye->P;
	peye->P.x += peye->V.x;
	peye->P.y += peye->V.y;
	
	Eye_S0_Fit (peye);
	
//	peye->V.x = peye->P.x - oldpos.x;
//	peye->V.y = peye->P.y - oldpos.y;
	peye->V.x = (peye->V.x + peye->P.x - oldpos.x) * 0.5f;
	peye->V.y = (peye->V.y + peye->P.y - oldpos.y) * 0.5f;
	
	peye->V.x *= 0.9f;
	peye->V.y *= 0.9f;
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
		float x = peye->P.x, y = peye->P.y;
		tPix prevpix = *dopix(x,y);
		for (; ddist2(x,y,peye->P.x,peye->P.y) < dpow2(peye->Exp_R*2); ) {
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

si	Eye_S3Fit_EdgeMark	(tEye* peye, si i, float t, float* px, float* py)
{
//	si border = peye->Exp_R*peye->Exp_R*3;
	u08 wrote = 0;
	float x = peye->P.x, y = peye->P.y;
	si n = 0;
	
	for (; ddist2(x,y,peye->P.x,peye->P.y) < peye->Exp_R*peye->Exp_R; ) {
		if (peye->LinView.x != 0)
			dspix(peye->LinView.x+n, peye->LinView.y+i) = dnpix(x,y)->Y | dnpix(x,y)->Y<<8 | dnpix(x,y)->Y<<16;
		
		if (	dopix(x,y)->Y > ay + 18
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
	float x = peye->P.x, y = peye->P.y;
	si n = 0;
	
	for (; ddist2(x,y,peye->P.x,peye->P.y) < dpow2(border_e); ) {
		if (peye->LinView.x != 0)
			dspix(peye->LinView.x+n, peye->LinView.y+i) = dnpix(x,y)->Y | dnpix(x,y)->Y<<8 | dnpix(x,y)->Y<<16;
		
		if (
		//	dopix(x,y)->Y > ay + 30
			dopix(x,y)->Y >= peye->Pix_Bright
		//	&& ddist2(x,y,peye->P.x,peye->P.y) >= 14*14
			&& ddist2(x,y,peye->P.x,peye->P.y) >= dpow2(border_s)
		) {
		//	return 1;
			float x1 = x, y1 = y;
			for (; ddist2(x1,y1,peye->P.x,peye->P.y) < dpow2(peye->Exp_R*2); ) {
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
si	Eye_S3Fit_EdgeMark3	(tEye* peye, si i, float t, float* px, float* py)
{
//	si border = peye->Exp_R*peye->Exp_R*3;
	u08 wrote = 0;
	float x = peye->P.x, y = peye->P.y, mindiff = peye->Exp_R*peye->Exp_R;
	si n = 0;
	
	for (; ddist2(x,y,peye->P.x,peye->P.y) < peye->Exp_R*peye->Exp_R; ) {
		if (peye->LinView.x != 0)
			dspix(peye->LinView.x+n, peye->LinView.y+i) = dnpix(x,y)->Y | dnpix(x,y)->Y<<8 | dnpix(x,y)->Y<<16;
		
		if (	dopix(x,y)->Y > ay + 18
		//	&& ddist2(x,y,peye->P.x,peye->P.y) >= 14*14
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
			//	if (!wrote) {
					float diff = ddist2(x,y,peye->P.x,peye->P.y) - peye->Exp_R*peye->Exp_R;
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
			dx = x0 - peye->P.x;
			dy = y0 - peye->P.y;
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
			continue;
		}
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


si	Eye_S4_EdgeMark00	(tEye* peye, float t, float* px, float* py)
{
	float border_s = peye->Min_R;
	float border_e = peye->Max_R;
	
	u08 wrote = 0;
	float x = peye->P.x, y = peye->P.y;
	si n = 0;
	
	for (; ddist2(x,y,peye->P.x,peye->P.y) < dpow2(border_e); ) {
	//	if (peye->LinView.x != 0)
	//		dspix(peye->LinView.x+n, peye->LinView.y+i) = dnpix(x,y)->Y | dnpix(x,y)->Y<<8 | dnpix(x,y)->Y<<16;
		
		if (
		//	dopix(x,y)->Y > ay + 30
			dopix(x,y)->Y >= peye->S4_Pix_Bright
			&& ddist2(x,y,peye->P.x,peye->P.y) >= dpow2(border_s)
		) {
		//	return 1;
			float x1 = x, y1 = y;
			for (; ddist2(x1,y1,peye->P.x,peye->P.y) < dpow2(peye->Exp_R*2); ) {
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
	float x = peye->P.x, y = peye->P.y, mindiff = peye->Exp_R*peye->Exp_R;
	si n = 0, minn = 0;
	
	for (; ddist2(x,y,peye->P.x,peye->P.y) < peye->Exp_R*peye->Exp_R; ) {
		
		if (	dopix(x,y)->Y > ay + 36
		//	&& ddist2(x,y,peye->P.x,peye->P.y) >= 14*14
		) {
		//	return 1;
			float x1 = x, y1 = y;
			for (; ddist2(x1,y1,peye->P.x,peye->P.y) < peye->Exp_R*peye->Exp_R*4; ) {
				if (	dopix(x1,y1)->Y <= ay + 32) {
					break;
				}
				x1 += cos(t);
				y1 += sin(t);
			}
			if (ddist2(x1,y1,x,y) >= 3*2) {
			//	if (!wrote) {
					float diff = ddist2(x,y,peye->P.x,peye->P.y) - peye->Exp_R*peye->Exp_R;
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
	float x = peye->P.x, y = peye->P.y, mindiff = peye->Exp_R*peye->Exp_R;
	si n = 0, minn = 0;
	
	tPix prev = *dopix(x,y);
	
	for (; ddist2(x,y,peye->P.x,peye->P.y) < peye->Exp_R*peye->Exp_R; )
	{
		if (	dopix(x,y)->Y - prev.Y >= 0x20
			&& ddist(x,y,peye->P.x,peye->P.y) > peye->Exp_R*0.6f
		) {
		//	Eye_Ellipse2LinDraw_Pix_ad (peye, t, n, 0xFF<<0);
			
			float diff = ddist2(x,y,peye->P.x,peye->P.y) - peye->Exp_R*peye->Exp_R;
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
	float x = peye->P.x, y = peye->P.y;
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
	
	for (; ddist2(x,y,peye->P.x,peye->P.y) < border_e*border_e; )
	{
		if (	dopix(x,y)->Y < min.pix.Y
		//	&& ddist(x,y,peye->P.x,peye->P.y) > peye->Exp_R*0.6f
		) {
			min.x = x;
			min.y = y;
			min.pix = *dopix(x,y);
		}
		if (	dopix(x,y)->Y > max.pix.Y
		//	&& ddist(x,y,peye->P.x,peye->P.y) > peye->Exp_R*0.6f
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
	
	if (ddist(min.x,min.y, peye->P.x,peye->P.y) < ddist(max.x,max.y, peye->P.x,peye->P.y)) {
		if (max.pix.Y - min.pix.Y < 0x30)
			return 0;
		
		*px = min.x + max.x;	*px /= 2.0f;
		*py = min.y + max.y;	*py /= 2.0f;
		Eye_Ellipse2LinDraw_Pix_ad (peye, t, ddist(*px,*py, peye->P.x,peye->P.y), 0xFF<<8);
		return 1;
	}
//	if (wrote) {
		
	//	Eye_Ellipse2LinDraw_Pix_ad (peye, t, ddist(min.x,min.y, peye->P.x,peye->P.y), 0xFF<<0);
	//	Eye_Ellipse2LinDraw_Pix_ad (peye, t, ddist(max.x,max.y, peye->P.x,peye->P.y), 0xFF<<8);
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
		Eye_Ellipse2LinDraw_Pix_ad (peye, t, ddist(*px,*py, peye->P.x,peye->P.y), 0xFF<<8);
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
	
	Eye_Ellipse2LinDraw_Pix_ad (peye, t, ddist(*px,*py, peye->P.x,peye->P.y), 0xFF<<8);
	return ret;
	
	
	u08 wrote = 0;
	float x = peye->P.x, y = peye->P.y;
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
	
	for (; ddist2(x,y,peye->P.x,peye->P.y) < border_e*border_e; )
	{
		if (	dopix(x,y)->Y < min.pix.Y
		//	&& ddist(x,y,peye->P.x,peye->P.y) > peye->Exp_R*0.6f
		) {
			min.x = x;
			min.y = y;
			min.pix = *dopix(x,y);
		}
		if (	dopix(x,y)->Y > max.pix.Y
		//	&& ddist(x,y,peye->P.x,peye->P.y) > peye->Exp_R*0.6f
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
	
	if (ddist(min.x,min.y, peye->P.x,peye->P.y) < ddist(max.x,max.y, peye->P.x,peye->P.y)) {
		if (max.pix.Y - min.pix.Y < 0x10)
			return 0;
		
		*px = min.x + max.x;	*px /= 2.0f;
		*py = min.y + max.y;	*py /= 2.0f;
		Eye_Ellipse2LinDraw_Pix_ad (peye, t, ddist(*px,*py, peye->P.x,peye->P.y), 0xFF<<8);
		return 1;
	}
//	if (wrote) {
		
	//	Eye_Ellipse2LinDraw_Pix_ad (peye, t, ddist(min.x,min.y, peye->P.x,peye->P.y), 0xFF<<0);
	//	Eye_Ellipse2LinDraw_Pix_ad (peye, t, ddist(max.x,max.y, peye->P.x,peye->P.y), 0xFF<<8);
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
				*pr += ddist(peye->P.x, peye->P.y, x0, y0);
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
		if (fabsf(ddist(p.x,p.y, peye->P.x, peye->P.y) - ddist(prev.x,prev.y, peye->P.x, peye->P.y)) > (2.5f))
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
				*pr += ddist(peye->P.x, peye->P.y, p.x, p.y);
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
	float t0 = atan2(y0 - peye->P.y, x0 - peye->P.x);
	
	struct {
		si idx;
		float diff_t;
	}min = {-1, deg2rad};
	si i;
	for (i = 0; i < peye->Point_N; ++i) {
		float x1 = peye->paPoint[i].x, y1 = peye->paPoint[i].y;
		float t1 = atan2(y1 - peye->P.y, x1 - peye->P.x);
		
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
	float t0 = atan2(y0 - peye->P.y, x0 - peye->P.x);
	
	struct {
		si idx;
		float diff_t;
	}min = {-1, 4*deg2rad};
//	printf ("Eye_S4_Fit_FindRot\n");
	si i;
	for (i = 0; i < peye->Point_N; ++i) {
		float x1 = peye->paPoint[i].x, y1 = peye->paPoint[i].y;
		float t1 = atan2(y1 - peye->P.y, x1 - peye->P.x);
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
		dx = x - peye->P.x;
		dy = y - peye->P.y;
		t = atan2(dy,dx);
		
		xx = (peye->Ax * cos(t) * cos(peye->Aa) - peye->Ay * sin(t) * sin(peye->Aa));
		yy = (peye->Ax * cos(t) * sin(peye->Aa) + peye->Ay * sin(t) * cos(peye->Aa));
		
		float diff = sqrt(dx*dx+dy*dy) - sqrt(xx*xx+yy*yy);
		if (fabs(diff) >= 0.01f) {
		//	printf ("f1\n");
			peye->Ax += peye->Fit_Scale*diff*fabs(cos(t));
			peye->Ay += peye->Fit_Scale*diff*fabs(sin(t));
			
			Eye_Xset (peye, peye->P.x + peye->Fit_Trans*(diff )*cos(t));
			Eye_Yset (peye, peye->P.y + peye->Fit_Trans*(diff )*sin(t));
		}
		dnpix(x,y)->Y = 0xFF;
		dnpix(x,y)->U = 0x0;
		dnpix(x,y)->V = 0x0;
	}/**/
/*	for (i = peye->Point_N-1; i >= 0; --i) {
		float x = peye->paPoint[i].x, y = peye->paPoint[i].y;
		float t, dx, dy, xx, yy;
		dx = x - peye->P.x;
		dy = y - peye->P.y;
		t = atan2(dy,dx);
		
		xx = (peye->Ax * cos(t) * cos(peye->Aa) - peye->Ay * sin(t) * sin(peye->Aa));
		yy = (peye->Ax * cos(t) * sin(peye->Aa) + peye->Ay * sin(t) * cos(peye->Aa));
		
		float diff = sqrt(dx*dx+dy*dy) - sqrt(xx*xx+yy*yy);
		if (fabs(diff) >= 0.01f) {
		//	printf ("f1\n");
			peye->Ax += peye->Fit_Scale*diff*fabs(cos(t));
			peye->Ay += peye->Fit_Scale*diff*fabs(sin(t));
			
			Eye_Xset (peye, peye->P.x + peye->Fit_Trans*(diff )*cos(t));
			Eye_Yset (peye, peye->P.y + peye->Fit_Trans*(diff )*sin(t));
		}
		dnpix(x,y)->Y = 0xFF;
		dnpix(x,y)->U = 0x0;
		dnpix(x,y)->V = 0x0;
	}/**/
	for (i = 0; i < peye->Point_N; ++i) {
		float x0 = peye->paPoint[i].x, y0 = peye->paPoint[i].y;
		float t0, dx, dy, ox0, oy0;
		dx = x0 - peye->P.x;
		dy = y0 - peye->P.y;
		t0 = atan2(dy,dx);
		
		ox0 = (peye->Ax * cos(t0) * cos(peye->Aa) - peye->Ay * sin(t0) * sin(peye->Aa));
		oy0 = (peye->Ax * cos(t0) * sin(peye->Aa) + peye->Ay * sin(t0) * cos(peye->Aa));
		
		{
			si idx = Eye_S4_Fit_FindOposite (peye, i);
			if (idx >= 0) {
				float x1 = peye->paPoint[idx].x, y1 = peye->paPoint[idx].y;
				float t1 = atan2(y1 - peye->P.y, x1 - peye->P.x);
				
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
			}else {
				dnpix(x0,y0)->Y = 0xFF;
				dnpix(x0,y0)->U = 0xF;
				dnpix(x0,y0)->V = 0xF;
				float diff = sqrt(dx*dx+dy*dy) - sqrt(ox0*ox0+oy0*oy0);
				if (fabs(diff) >= 0.01f) {
				//	printf ("f1\n");
					peye->Ax += peye->Fit_Scale*diff*fabs(cos(t0));
					peye->Ay += peye->Fit_Scale*diff*fabs(sin(t0));
					
					Eye_Xset (peye, peye->P.x + peye->Fit_Trans*(diff )*cos(t0));
					Eye_Yset (peye, peye->P.y + peye->Fit_Trans*(diff )*sin(t0));
				}
			}
		}
	}/**/
	float avgerr = 0;
	for (i = 0; i < peye->Point_N; ++i) {
		float x = peye->paPoint[i].x, y = peye->paPoint[i].y;
		float t, dx, dy, xx, yy;
		dx = x - peye->P.x;
		dy = y - peye->P.y;
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
	int* ret = pupil_fitting_inliers (videoIn->width, videoIn->height, &max_inliers_num);
//	int* ret = simple_fit (videoIn->width, videoIn->height, &max_inliers_num);
	
	//double pupil_param[5];//parameters of an ellipse {ellipse_a, ellipse_b, cx, cy, theta}; a & b is the major or minor axis; 
	
//	printf ("%d: x %f y %f  Ax %f Ay %f   Aa %f\n", max_inliers_num, pupil_param[2], pupil_param[3], pupil_param[0], pupil_param[1], pupil_param[4]);
	
	float x = peye->P.x, y = peye->P.y, ax = peye->Ax, ay = peye->Ay, aa = peye->Aa;
	
	peye->P.x = pupil_param[2];
	peye->P.y = pupil_param[3];
	peye->Ax = pupil_param[0];
	peye->Ay = pupil_param[1];
	peye->Aa = pupil_param[4];
	
//	Eye_Draw_Ellipse (peye, 0.0f, M_PI*2);
	
	if (ddist(peye->P.x, peye->P.y, x, y) >= dpow2(6)) {
		peye->P.x = x;
		peye->P.y = y;
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
		dx = x0 - peye->P.x;
		dy = y0 - peye->P.y;
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
				float t1 = atan2(y1 - peye->P.y, x1 - peye->P.x);
				
				float x2 = peye->paPoint[i2].x, y2 = peye->paPoint[i2].y;
				float t2 = atan2(y2 - peye->P.y, x2 - peye->P.x);
				
				
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
				
				float pdx = (x0+x1)*0.5f - (2.0f*peye->P.x + ox0+ox1)*0.5f;
				float pdy = (y0+y1)*0.5f - (2.0f*peye->P.y + oy0+oy1)*0.5f;
				
				Eye_Xset (peye, peye->P.x + 0.10f * pdx );// *cos(t));
				Eye_Yset (peye, peye->P.y + 0.10f * pdy );// *sin(t));
				
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
					float t1 = atan2(y1 - peye->P.y, x1 - peye->P.x);
					
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
					continue;
				}
			}
			if (1) {
				float diff = sqrt(dx*dx+dy*dy) - sqrt(ox0*ox0+oy0*oy0);
				if (fabs(diff) >= 0.01f) {
				//	printf ("f1\n");
					peye->Ax += 0.1f*diff*fabs(cos(t0));
					peye->Ay += 0.1f*diff*fabs(sin(t0));
					
					Eye_Xset (peye, peye->P.x + 0.1f*(diff)*cos(t0));
					Eye_Yset (peye, peye->P.y + 0.1f*(diff)*sin(t0));
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
	float x, y, ox = peye->P.x, oy = peye->P.y;
	for (y = oy - dd; y <= oy + dd; y += 1) {
		for (x = ox - dd; x <= ox + dd; x += 1) {
			teye = *peye;
			teye.P.x = x;
			teye.P.y = y;
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
	float x = peye->P.x, y = peye->P.y;
	si n = 0;
	u08 got_dark = 0;
	
	for (; ddist2(x,y,peye->P.x,peye->P.y) < dpow2(peye->Exp_R*2); ) {
		if (!got_dark) {
			if (
			//	dopix(x,y)->Y < ay + 8
				dopix(x,y)->Y <= peye->S5.Pix_Dark
			) {
				Eye_Ellipse2LinDraw_Pix_ad (peye, t, ddist(x,y, peye->P.x,peye->P.y), 0x80<<8);
				got_dark = 1;
			}
		}else if (
		//	dopix(x,y)->Y > ay + 40
			dopix(x,y)->Y >= peye->S5.Pix_Bright
		//	&& ddist2(x,y,peye->P.x,peye->P.y) >= dpow2(peye->Min_R)
		) {
		//	return 1;
			float x1 = x, y1 = y;
			for (; ddist2(x1,y1,peye->P.x,peye->P.y) < dpow2(peye->Exp_R*2); ) {
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
			//	Eye_Ellipse2LinDraw_Pix_ad (peye, t, ddist(x,y, peye->P.x,peye->P.y), 0xFF<<8);
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
					
					Eye_Ellipse2LinDraw_Pix_ad (peye, t, ddist(papoint[n].x,papoint[n].y, peye->P.x,peye->P.y), 0xFF<<8);
				}else
					Eye_Ellipse2LinDraw_Pix_ad (peye, t, ddist(papoint[n].x,papoint[n].y, peye->P.x,peye->P.y), 0xFF<<8);
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
	float x0 = peye->P.x, y0 = peye->P.y;
	si n = 0;
	u08 got_fail = 0;
	
	for (; ddist2(x0,y0,peye->P.x,peye->P.y) < dpow2(peye->Exp_R*2 - peye->S5.Diff_Dist); )
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
			
			Eye_Ellipse2LinDraw_Pix_ad (peye, t, ddist(x0,y0, peye->P.x,peye->P.y), 0x80<<8);
			Eye_Ellipse2LinDraw_Pix_ad (peye, t, ddist(x1,y1, peye->P.x,peye->P.y), 0xFF<<8);
			Eye_Ellipse2LinDraw_Pix_ad (peye, t, ddist(papoint[n].x,papoint[n].y, peye->P.x,peye->P.y), 0xFF<<16);
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
	//	Eye_Ellipse2LinDraw_Pix_ad (peye, a, ddist(point[ns].x,point[ns].y, peye->P.x,peye->P.y), col);
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
		//	Eye_Ellipse2LinDraw_Pix_ad (peye, a, ddist(point[mini].x,point[mini].y, peye->P.x,peye->P.y), 0xFF<<16);
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
				angle_norm_0_2pi(atan2(peye->paPoint[ip].y-peye->P.y,peye->paPoint[ip].x-peye->P.x)),
				ddist(peye->paPoint[ip].x,peye->paPoint[ip].y, peye->P.x,peye->P.y),
				col);
			
		}
	}/**/
	
	for (i = 0; i < edge_n; ++i) {
	//	printf ("\t%ld: r %f\n", edge[i].num, r[i]);
		edge[i].r = 0;
		si ip;
		for (ip = edge[i].start; ip < edge[i].start+edge[i].num; ++ip) {
			edge[i].r += ddist(peye->paPoint[ip].x,peye->paPoint[ip].y, peye->P.x,peye->P.y);
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
			angle_norm_0_2pi(atan2(peye->paPoint[i].y-peye->P.y,peye->paPoint[i].x-peye->P.x)),
			ddist(peye->paPoint[i].x,peye->paPoint[i].y, peye->P.x,peye->P.y),
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

void	Eye_FF_MarkCrap		(tEye* peye, si x, si y)
{
	if (x < 0 || y < 0 || x >= videoIn->width || y >= videoIn->height)
		return;
	
	if (	1//dnpix(x,y)->Y == 0xFF
		&& dnpix(x,y)->U == 0
		&& dnpix(x,y)->V == 0
	)
		return;
	
	if (ddist2(x,y,peye->P.x,peye->P.y) >= dpow2(peye->Max_R))
		return;
	
	if (	dopix(x,y)->Y < ay
		//abs(dopix(x,y)->Y - ay) <= 16
	//	&& abs(dopix(x,y)->U - au) <= 12
	//	&& abs(dopix(x,y)->V - av) <= 12
	) {
		++peye->FF.tmpNum;
	//	dnpix(x,y)->Y = 0xFF;
		dnpix(x,y)->U = 0x0;
		dnpix(x,y)->V = 0x0;
		Eye_crap_ff (peye, x-1, y);
		Eye_crap_ff (peye, x+1, y);
		Eye_crap_ff (peye, x, y-1);
		Eye_crap_ff (peye, x, y+1);
	}
}

void	Eye_FF_Mark2Pos	(tEye* peye)
{
	si ax = 0, ay = 0;
	si num = 0;
	si x, y;
	for (y = peye->P.y - peye->Max_R-1; y < peye->P.y + peye->Max_R+1; ++y) {
		for (x = peye->P.x - peye->Max_R-1; x < peye->P.x + peye->Max_R+1; ++x) {
			if (	1//dnpix(x,y)->Y == 0xFF
				&& dnpix(x,y)->U == 0x0
				&& dnpix(x,y)->V == 0x0
			) {
				ax += x;
				ay += y;
				++num;
			}
		}
	}/**/
	if (num > 20) {
		peye->P.x = (float)ax / num;
		peye->P.y = (float)ay / num;
	}
}


void	Eye_FF		(tEye* peye)
{
	
	ay = peye->FF.Y;
	au = 0;
	av = 0;
	
	peye->FF.tmpNum = 0;
	
//	Eye_FF_MarkCrap (peye, peye->P.x, peye->P.y);
	si search_r = peye->FF.Search_R;
	si x, y;
	for (y = peye->P.y - search_r; y < peye->P.y + search_r; ++y) {
		for (x = peye->P.x - search_r; x < peye->P.x + search_r; ++x) {
			Eye_FF_MarkCrap (peye, x, y);
		}
	}/**/
	
	Eye_FF_Mark2Pos (peye);
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
	
	Eye_CalcAYUV (peye);
	Eye_S5 (peye);
	Eye_CalcAYUV (peye);
	Eye_S3Fit (peye);
	
}

void	Eye_Fit		(tEye* peye)
{
	switch (peye->Fit) {
	case eEye_Fit_S0:
		Eye_CalcAYUV (peye);
		Eye_S0 (peye);
		break;
	case eEye_Fit_S2Fit:
		Eye_OldFF (peye);
		Eye_EdgeMark (peye);
		Eye_S2Fit (peye);
		break;
	case eEye_Fit_S3Fit_START ... eEye_Fit_S3Fit_END:
		Eye_CalcAYUV (peye);
		Eye_S3Fit (peye);
		break;
	case eEye_Fit_S4Fit_START ... eEye_Fit_S4Fit_END:
		Eye_CalcAYUV (peye);
		Eye_S4 (peye);
		break;
	case eEye_Fit_S5_START ... eEye_Fit_S5_END:
		Eye_CalcAYUV (peye);
		Eye_S5 (peye);
		break;
	case eEye_Fit_FF_START ... eEye_Fit_FF_END:
		Eye_FF (peye);
		break;
	case eEye_Fit_C:
		Eye_CFit (peye);
		break;
	case eEye_Fit_SFit:
	default:
		Eye_OldFF (peye);
		Eye_EdgeMark (peye);
		Eye_SFit (peye);
		break;
	}
}

void	Eye_Clip		(tEye* peye)
{
	si x, y;
	if (peye->P.y < 100)
		peye->P.y = 100;
	if (peye->P.y > videoIn->height-100)
		peye->P.y = videoIn->height-100;
	
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
	//	y = peye->P.y;
		if (peye->P.x+x < 0 || peye->P.y+y < 0 || peye->P.x+x >= videoIn->width || peye->P.y+y >= videoIn->height)
			continue;
		dnpix(peye->P.x+x,peye->P.y+y)->Y = 0xFF;
		dnpix(peye->P.x+x,peye->P.y+y)->U = 0xF;
		dnpix(peye->P.x+x,peye->P.y+y)->V = 0xF;
	}/**/
}
void	Eye_Draw		(tEye* peye)
{
	si x, y;
	dset_c1(peye->P.x,peye->P.y);
//	Eye_Draw_Ellipse (peye, 0, M_PI*2);	return;
/*	Eye_Draw_Ellipse (peye, 0 - 10*M_PI/180,		0 + 10*M_PI/180);
	Eye_Draw_Ellipse (peye, M_PI_2 - 10*M_PI/180,	M_PI_2 + 10*M_PI/180);
	Eye_Draw_Ellipse (peye, M_PI - 10*M_PI/180,	M_PI + 10*M_PI/180);
	Eye_Draw_Ellipse (peye, -M_PI_2 - 10*M_PI/180,	-M_PI_2 + 10*M_PI/180);
	return;
	/**/
	dset_c1(peye->P.x - peye->Ax,		peye->P.y - peye->Ay);
	dset_c1(peye->P.x - peye->Ax+1,	peye->P.y - peye->Ay);
	dset_c1(peye->P.x - peye->Ax+2,	peye->P.y - peye->Ay);
	dset_c1(peye->P.x - peye->Ax,		peye->P.y - peye->Ay+1);
	dset_c1(peye->P.x - peye->Ax,		peye->P.y - peye->Ay+2);
	
	dset_c1(peye->P.x + peye->Ax,		peye->P.y - peye->Ay);
	dset_c1(peye->P.x + peye->Ax-1,	peye->P.y - peye->Ay);
	dset_c1(peye->P.x + peye->Ax-2,	peye->P.y - peye->Ay);
	dset_c1(peye->P.x + peye->Ax,		peye->P.y - peye->Ay+1);
	dset_c1(peye->P.x + peye->Ax,		peye->P.y - peye->Ay+2);
	
	dset_c1(peye->P.x - peye->Ax,		peye->P.y + peye->Ay);
	dset_c1(peye->P.x - peye->Ax+1,	peye->P.y + peye->Ay);
	dset_c1(peye->P.x - peye->Ax+2,	peye->P.y + peye->Ay);
	dset_c1(peye->P.x - peye->Ax,		peye->P.y + peye->Ay-1);
	dset_c1(peye->P.x - peye->Ax,		peye->P.y + peye->Ay-2);
	
	dset_c1(peye->P.x + peye->Ax,		peye->P.y + peye->Ay);
	dset_c1(peye->P.x + peye->Ax-1,	peye->P.y + peye->Ay);
	dset_c1(peye->P.x + peye->Ax-2,	peye->P.y + peye->Ay);
	dset_c1(peye->P.x + peye->Ax,		peye->P.y + peye->Ay-1);
	dset_c1(peye->P.x + peye->Ax,		peye->P.y + peye->Ay-2);
	/**/
	
}


void	Cam_Retina_Ray	(tEye* peye, tV4f* pp0, tV4f* pp1, tV4f* pvec)
{
	tV4f p0 = {0, 0, 0, 1};
	tV4f p1 = {0, 0, -1, 1};
	
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
		p1.x = gM.Proj_L + (peye->P.x/gM.Cam.Image_W)*gM.Proj_W;
		p1.y = gM.Proj_B + (peye->P.y/gM.Cam.Image_H)*gM.Proj_H;
	}
	if (pp0)
		*pp0 = p0;
	if (pp1)
		*pp1 = p1;
	if (pvec) {
		*pvec = p1;
	}
}


void	Cam_Param_Set	(tCam* pcam)
{
	si value;
	#define dcam_set(_id,_val)		\
		do {		\
			struct v4l2_control control;		\
			control.id    = _id;		\
			control.value = _val;		\
			if ((value = ioctl(videoIn->fd, VIDIOC_S_CTRL, &control)) < 0)		\
				printf("Set " #_id " error\n");		\
		/*	else	printf(#_id " set to %d\n", control.value);/**/		\
		}while (0)
	
	dcam_set (V4L2_CID_AUTO_WHITE_BALANCE, 1);
	dcam_set (V4L2_CID_AUTOGAIN, 1);
	
	if (gM.Cam.Focus == -1)
		dcam_set(V4L2_CID_FOCUS_AUTO, 1);
	else {
		dcam_set(V4L2_CID_FOCUS_AUTO, 0);
		
	//	if ((value = v4l2SetControl(videoIn, V4L2_CID_FOCUS_ABSOLUTE, gM.Cam.Focus-30)) < 0)
	//		printf("Set CT_FOCUS_ABSOLUTE_CONTROL to %ld error\n", value);
		if ((value = v4l2SetControl(videoIn, V4L2_CID_FOCUS_ABSOLUTE, gM.Cam.Focus)) < 0)
			printf("Set CT_FOCUS_ABSOLUTE_CONTROL to %ld error\n", value);
	}
	if (gM.Cam.Exposure == -1) {
		//V4L2_EXPOSURE_SHUTTER_PRIORITY	V4L2_EXPOSURE_APERTURE_PRIORITY
		dcam_set(V4L2_CID_EXPOSURE_AUTO, V4L2_EXPOSURE_AUTO);
		dcam_set(V4L2_CID_EXPOSURE_AUTO_PRIORITY, V4L2_EXPOSURE_APERTURE_PRIORITY);
		dcam_set(V4L2_CID_EXPOSURE_AUTO, V4L2_EXPOSURE_AUTO);
	//	if ((value = v4l2SetControl(videoIn, V4L2_CID_EXPOSURE_AUTO, V4L2_EXPOSURE_AUTO)) < 0)
	//		printf("Set V4L2_CID_EXPOSURE_AUTO to %ld error\n", value);
	}else {
		if ((value = v4l2SetControl(videoIn, V4L2_CID_EXPOSURE_AUTO, V4L2_EXPOSURE_MANUAL)) < 0)
			printf("Set V4L2_CID_EXPOSURE_AUTO to %ld error\n", value);
	//	dcam_set(V4L2_CID_EXPOSURE_AUTO, V4L2_EXPOSURE_MANUAL);
		
	//	if ((value = v4l2SetControl(videoIn, V4L2_CID_EXPOSURE_ABSOLUTE, gM.Cam.Exposure-30)) < 0)
	//		printf("Set V4L2_CID_EXPOSURE_ABSOLUTE to %ld error\n", value);
		if ((value = v4l2SetControl(videoIn, V4L2_CID_EXPOSURE_ABSOLUTE, gM.Cam.Exposure)) < 0)
			printf("Set V4L2_CID_EXPOSURE_ABSOLUTE to %ld error\n", value);
	}
	if ((value = v4l2SetControl(videoIn, V4L2_CID_ZOOM_ABSOLUTE, gM.Cam.Zoom)) < 0)
		printf("Set V4L2_CID_ZOOM_ABSOLUTE to %ld error\n", value);
	#undef dcam_set
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

void	Head_Calc_M4_Rel	(tHead* p, tHead* pcen)
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
		tV2f osl = p->DotL.P;
		tV2f osr = p->DotR.P;
		tV2f osc = p->DotC.P;
		
		tV2f os_vrl = osl;	V2f_sub_V2f (&os_vrl, &osr);
		tV2f os_vlc = osc;	V2f_sub_V2f (&os_vlc, &osl);
		tV2f os_vrc = osc;	V2f_sub_V2f (&os_vrc, &osr);
		
		float osdlr = ddist (osl.x,osl.y,osr.x,osr.y);
		float osdlc = ddist (osl.x,osl.y,osc.x,osc.y);
		float osdrc = ddist (osr.x,osr.y,osc.x,osc.y);
		
		p->R_X = p->R_Y = p->R_Z = 0;
		static tM4f rot = {1,0,0,0, 0,1,0,0, 0,0,1,0, 0,0,0,1};
		static tM4f trans = {1,0,0,0, 0,1,0,0, 0,0,1,0, 0,0,0,1};
		M4f_Iden (&rot);	M4f_Iden (&trans);	M4f_trans (&trans, 0, 0, -30);
		
		float err = 0;
		static si i = 0;
		for (i = 0; i < 400; ++i) {
			M4f_Iden (&p->M4);
			M4f_mul_M4f (&p->M4, &rot);
			M4f_mul_M4f (&p->M4, &trans);
			
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
			
		/*	V4f_mul_M4f (&pc, &p->M4);
			V4f_mul_M4f (&pl, &p->M4);
			V4f_mul_M4f (&pr, &p->M4);/**/
			
			M4f_mul_V4f (&p->M4, &pc);
			M4f_mul_V4f (&p->M4, &pl);
			M4f_mul_V4f (&p->M4, &pr);
			
			tV2f sc, sl, sr;
			V4f_ScreenPosf (&pl, &sl);
			V4f_ScreenPosf (&pr, &sr);
			V4f_ScreenPosf (&pc, &sc);
			
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
				M4f_rotz (&rot, p->Mod.RInc.z*dz);
			}/**/
			
			if (1) {
				float ocross = V2f_cross (&os_vlc, &os_vrc) / (osdlc*osdrc);
				float cross = V2f_cross (&s_vlc, &s_vrc) / (sdlc*sdrc);
				dx = -ocross + cross;
			//	printf ("rotx %f\n", dx);
				p->R_X += -p->Mod.RInc.x*dx;
				M4f_rotx (&rot, -p->Mod.RInc.x*dx);
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
				M4f_roty (&rot, -p->Mod.RInc.y*dy);
			}/**/
			
		/*	if (++i >= 100) {
				i = 0;
				M4f_Iden (&rot);
				M4f_Iden (&trans);
				M4f_trans (&trans, 0, 0, -80);
			}/**/
		}
		if (fabsf(err) > 5) {
			M4f_Iden (&rot);	M4f_Iden (&trans);	M4f_trans (&trans, 0, 0, -40);
		}/**/
		p->M4_T = trans;
		p->M4_R = rot;
		
		M4f_Iden (&p->M4);
		M4f_mul_M4f (&p->M4, &p->M4_R);
		M4f_mul_M4f (&p->M4, &p->M4_T);
		
		p->P.x = p->M4.x03;
		p->P.y = p->M4.x13;
		p->P.z = p->M4.x23;
		
		p->N.x = 0;
		p->N.y = 0;
		p->N.z = 1;
		p->N.w = 1;
		M4f_mul_V4f (&p->M4_R, &p->N);
		
	//	printf ("p->P: ");	V4f_Print (&p->P);
	//	printf ("p->N: ");	V4f_Print (&p->N);
	//	V4f_sub_V4f (&p->P, &p->N);
	//	V4f_sub_V4f (&p->P, &p->N);
	//	V4f_sub_V4f (&p->P, &p->N);
	//	V4f_sub_V4f (&p->P, &p->N);
		
		M4f_Inv (&p->M4, &p->M4I);
		M4f_Inv (&p->M4_R, &p->M4I_R);
		M4f_Inv (&p->M4_T, &p->M4I_T);
		
		return;
	}
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
		p1.x = gM.Proj_L + (peye->P.x/gM.Cam.Image_W)*gM.Proj_W;
		p1.y = gM.Proj_B + (peye->P.y/gM.Cam.Image_H)*gM.Proj_H;
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

void	Head_Eye_Line	(tHead* p, tEye* peye, tV4f* ppos, tV4f* pvec)
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
		p1.x = gM.Proj_L + (peye->P.x/gM.Cam.Image_W)*gM.Proj_W;
		p1.y = gM.Proj_B + (peye->P.y/gM.Cam.Image_H)*gM.Proj_H;
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
		gCol.U = 0x0;
		gCol.V = 0x0;
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
			
			V4f_DrawPosPos (&p0, &p1);
		}
	}
	gCol.U = 0xF;
	gCol.V = 0xF;
	tV4f pos = peye->InHead.P;
	M4f_mul_V4f (&p->M4, &pos);
	V4f_DrawPosPos (&pos, &pos);
}

void	Head_Eye_Vector	(tHead* p, tEye* peye, tV4f* pret)	//Point of retina in relation to head
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
void	Screen_Cal_Save_Eye	(tScreen* p, tEye* peye)
{
//	peye->aScreen[p->Idx].aaCal[p->Cal.iy][p->Cal.ix].SX = p->Cal.sx;
//	peye->aScreen[p->Idx].aaCal[p->Cal.iy][p->Cal.ix].SY = p->Cal.sy;
	
	tV4f leye, lret, lvec;
	leye = peye->InHead.P;
	Head_Eye_Vector (&gM.Head, peye, &lret);
	lvec = lret;	V4f_sub_V4f (&lvec, &peye->InHead.P);
	
	M4f_mul_V4f (&gM.Head.M4, &leye);	M4f_mul_V4f (&gM.Head.M4, &lret);	M4f_mul_V4f (&gM.Head.M4, &lvec);
	
	tV4f pos = (tV4f){0, 0, gM.aScreen[p->Idx].C.z, 1};
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
//	printf ("Screen_Cal_Save Idx %d\n", p->Idx);
	
//	printf ("hoho1 ");	V4f_Print (&ret);
//	printf ("vec ");	V4f_Print (&vec);
	
	Screen_Cal_Save_Eye (p, &gM.Left);
	Screen_Cal_Save_Eye (p, &gM.Right);
	
	Screen_Eye_InterClear (p);
	Screen_Eye_InterAll (p);
	Screen_Eye_InterAll (p);
	Screen_Eye_InterAll (p);
	
	Screen_Eye_ExtraAll (p);
	
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
			case 0:				dcal(ix,iy).SX += gM.CalPoint.Win_W;	break;
			case dEye_Screen_Cal_LAST:	dcal(ix,iy).SX -= gM.CalPoint.Win_W;	break;
			}
			
			dcal(ix,iy).SY = iy * p->PixH/dEye_Screen_Cal_LAST;
			switch (iy) {
			case 0:				dcal(ix,iy).SY += gM.CalPoint.Win_W;	break;
			case dEye_Screen_Cal_LAST:	dcal(ix,iy).SY -= gM.CalPoint.Win_W;	break;
			}
		}
	}
}


u08	Screen_Eye_Point	(tEye* peye, tV4f* ppos, float* ps, float *pt)
{
	tV4f p0, p1, t0 = {0, 0, 0, 1}, t1, t2;
	Head_Eye_Vector (&gM.Head, peye, &p1);
	
	p0 = peye->InHead.P;
	
	M4f_mul_V4f (&gM.Head.M4, &p0);
	M4f_mul_V4f (&gM.Head.M4, &p1);
	
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

u08	Screen_Eye_XY_T	(tScreen* p, tEye* peye, si ix, si iy, float *px, float *py)
{
	tV4f p0, p1;
	Head_Eye_Vector (&gM.Head, peye, &p1);
	
	p0 = peye->InHead.P;
	
	M4f_mul_V4f (&gM.Head.M4, &p0);
	M4f_mul_V4f (&gM.Head.M4, &p1);
	
	float s, t;
	V4f_Intersect_Line01_Tri012 (
		0, &s, &t,
		&p0, &p1,
		&peye->aScreen[p->Idx].aaCal[iy+0][ix+0].P, &peye->aScreen[p->Idx].aaCal[iy+0][ix+1].P, &peye->aScreen[p->Idx].aaCal[iy+1][ix+0].P
	);
	
//	printf ("st %f\t%f\n", s, t);
	*px = peye->aScreen[p->Idx].aaCal[iy+0][ix+0].SX + s * (peye->aScreen[p->Idx].aaCal[iy+0][ix+1].SX - peye->aScreen[p->Idx].aaCal[iy+0][ix+0].SX);
	*py = peye->aScreen[p->Idx].aaCal[iy+0][ix+0].SY + t * (peye->aScreen[p->Idx].aaCal[iy+1][ix+0].SY - peye->aScreen[p->Idx].aaCal[iy+0][ix+0].SY);
	
	if (s < -0.001 || s > 1.001)
		return 0;
	if (t < -0.001 || t > 1.001)
		return 0;
	return 1;
}

void	Screen_Eye_XY_B	(tScreen* p, tEye* peye, si ix, si iy, float *px, float *py)
{
	tV4f p0, p1;
	Head_Eye_Vector (&gM.Head, peye, &p1);
	
	p0 = peye->InHead.P;
	
	M4f_mul_V4f (&gM.Head.M4, &p0);
	M4f_mul_V4f (&gM.Head.M4, &p1);
	
	float s, t;
	V4f_Intersect_Line01_Tri012 (
		0, &s, &t,
		&p0, &p1,
		&peye->aScreen[p->Idx].aaCal[iy+1][ix+1].P, &peye->aScreen[p->Idx].aaCal[iy+1][ix+0].P, &peye->aScreen[p->Idx].aaCal[iy+0][ix+1].P
	);
	
//	printf ("st %f\t%f\n", s, t);
	*px = peye->aScreen[p->Idx].aaCal[iy+1][ix+1].SX + s * (peye->aScreen[p->Idx].aaCal[iy+1][ix+0].SX - peye->aScreen[p->Idx].aaCal[iy+1][ix+1].SX);
	*py = peye->aScreen[p->Idx].aaCal[iy+1][ix+1].SY + t * (peye->aScreen[p->Idx].aaCal[iy+0][ix+1].SY - peye->aScreen[p->Idx].aaCal[iy+1][ix+1].SY);
}

s08	Screen_Eye_XY	(tScreen* p, tEye* peye, float *px, float *py)
{
	#if 0
	float x0, y0;
	float x1, y1;
	Screen_Eye_XY_TR (p, peye, &x0, &y0);
	if (0) {
		*px = x0;	*py = y0;
		return;
	}
	Screen_Eye_XY_BL (p, peye, &x1, &y1);
	
	float marg = 20;
	float a = (float)gM.aScreen[p->Idx].PixH / (float)gM.aScreen[p->Idx].PixW;
//	printf ("a %f   y0 %f   x0*a %f\n", a, y0, x0*a);
	if (y0+marg < x0*a) {
	//	printf ("tr\n");
		*px = x0;
		*py = y0;
	}else if (y0-marg > x0*a) {
	//	printf ("bl\n");
		*px = x1;
		*py = y1;
	}else {
		*px = (x0 + x1) / 2;
		*py = (y0 + y1) / 2;
	}
	float outmarg = 20;
	if (*px < -outmarg || *py < -outmarg) 
		return 0;
	if (*px > p->PixW+outmarg || *py > p->PixH+outmarg) 
		return 0;
	return 1;
	#else
	float ax = 0, ay = 0;
	si num = 0;
	si ix, iy;
	for (iy = 0; iy < dEye_Screen_Cal_NUM-1; ++iy) {
		for (ix = 0; ix < dEye_Screen_Cal_NUM-1; ++ix) {
			float x0, y0, x1, y1, tx, ty;
			u08 top = Screen_Eye_XY_T (p, peye, ix, iy, &x0, &y0);
			Screen_Eye_XY_B (p, peye, ix, iy, &x1, &y1);
		//	if (!top)
		//		continue;
			
			tx = x0;
			ty = y0;
			float marg = 20;
			float a = (float)(peye->aScreen[p->Idx].aaCal[iy+0][ix+1].SY-peye->aScreen[p->Idx].aaCal[iy+1][ix+0].SY)
					/ (float)(peye->aScreen[p->Idx].aaCal[iy+0][ix+1].SX - peye->aScreen[p->Idx].aaCal[iy+1][ix+0].SX);
			float b = peye->aScreen[p->Idx].aaCal[iy+1][ix+0].SY - a*peye->aScreen[p->Idx].aaCal[iy+1][ix+0].SX;
		//	printf ("a %f   y0 %f   x0*a %f\n", a, y0, x0*a);
			if (y0+marg < x0*a+b) {
			//	printf ("t\n");
				tx = x0;
				ty = y0;
			}else if (y0-marg > x0*a+b) {
			//	printf ("b\n");
				tx = x1;
				ty = y1;
			}else {
				tx = (x0 + x1) / 2;
				ty = (y0 + y1) / 2;
			}/**/
			float outmarg = 30;
			if (tx < peye->aScreen[p->Idx].aaCal[iy][ix].SX-outmarg || ty < peye->aScreen[p->Idx].aaCal[iy][ix].SY-outmarg) 
				continue;
			if (tx > peye->aScreen[p->Idx].aaCal[iy+1][ix+1].SX+outmarg || ty > peye->aScreen[p->Idx].aaCal[iy+1][ix+1].SY+outmarg) 
				continue;/**/
			
		//	printf ("adding %f %f\n", tx, ty);
			ax += tx;
			ay += ty;
			++num;
		}
	}
//	printf ("num %ld\n", num);
	if (num == 0)
		return 0;
	
	*px = ax/num;
	*py = ay/num;
	return 1;
	#endif
}

void	Screen_Eye_Print	(tScreen* p, tEye* peye)
{
//	printf ("Screen_Eye_Print Idx %d\n", p->Idx);
	
	if (peye == &gM.Left)
		gColARGB = 0x00FF00;
	else
		gColARGB = 0xFF0000;
	
	if (1) {
		float ox = peye->aScreen[p->Idx].View.Top.x, oy = peye->aScreen[p->Idx].View.Top.y;
		
		if (1) {
			ui ix, iy = 0;
			for (iy = 0; iy < dEye_Screen_Cal_NUM; ++iy) {
				for (ix = 1; ix < dEye_Screen_Cal_NUM; ++ix) {
					S_Draw_Line_2d (
						ox+peye->aScreen[p->Idx].aaCal[iy][ix-1].P.x,
						oy+peye->aScreen[p->Idx].aaCal[iy][ix-1].P.z,
						ox+peye->aScreen[p->Idx].aaCal[iy][ix].P.x,
						oy+peye->aScreen[p->Idx].aaCal[iy][ix].P.z
					);
				}
			}
		}
	/*	{
			float s, t;
			tV4f p0;
			Screen_Eye_Point (peye, &p0, &s, &t);
			if (ox+p0.x-1 >= 0 && ox+p0.x+1 < gM.pScreen->w
				&& oy+p0.z-1 >= 0 && oy+p0.z+1 < gM.pScreen->h
			) {
			//	printf ("st  %f\t%f\n", s, t);
			//	printf ("p0 ");	V4f_Print (&p0);
				dspix(ox+p0.x-1, oy+p0.z-1) = 0xFFFFFF;
				dspix(ox+p0.x+0, oy+p0.z-1) = 0xFFFFFF;
				dspix(ox+p0.x+1, oy+p0.z-1) = 0xFFFFFF;
				
				dspix(ox+p0.x-1, oy+p0.z+0) = 0xFFFFFF;
				dspix(ox+p0.x, oy+p0.z) = 0x0000FF;
				dspix(ox+p0.x+1, oy+p0.z+0) = 0xFFFFFF;
				
				dspix(ox+p0.x-1, oy+p0.z+1) = 0xFFFFFF;
				dspix(ox+p0.x+0, oy+p0.z+1) = 0xFFFFFF;
				dspix(ox+p0.x+1, oy+p0.z+1) = 0xFFFFFF;
			}
		}/**/
	}
	if (1) {
		float ox = peye->aScreen[p->Idx].View.Left.x, oy = peye->aScreen[p->Idx].View.Left.y;
		
		if (1) {
			ui ix = 0, iy = 0;
			for (iy = 1; iy < dEye_Screen_Cal_NUM; ++iy) {
				for (ix = 0; ix < dEye_Screen_Cal_NUM; ++ix) {
					S_Draw_Line_2d (
						ox+peye->aScreen[p->Idx].aaCal[iy-1][ix].P.z,
						oy+peye->aScreen[p->Idx].aaCal[iy-1][ix].P.y,
						ox+peye->aScreen[p->Idx].aaCal[iy][ix].P.z,
						oy+peye->aScreen[p->Idx].aaCal[iy][ix].P.y
					);
				}
			}
		}
	/*	{
			float s, t;
			tV4f p0;
			Screen_Eye_Point (peye, &p0, &s, &t);
			
			if (ox+p0.z-1 >= 0 && ox+p0.z+1 < gM.pScreen->w
				&& oy+p0.y-1 >= 0 && oy+p0.y+1 < gM.pScreen->h
			) {
				dspix(ox+p0.z-1, oy+p0.y-1) = 0xFFFFFF;
				dspix(ox+p0.z+0, oy+p0.y-1) = 0xFFFFFF;
				dspix(ox+p0.z+1, oy+p0.y-1) = 0xFFFFFF;
				
				dspix(ox+p0.z-1, oy+p0.y+0) = 0xFFFFFF;
				dspix(ox+p0.z, oy+p0.y) = 0x0000FF;
				dspix(ox+p0.z+1, oy+p0.y+0) = 0xFFFFFF;
				
				dspix(ox+p0.z-1, oy+p0.y+1) = 0xFFFFFF;
				dspix(ox+p0.z+0, oy+p0.y+1) = 0xFFFFFF;
				dspix(ox+p0.z+1, oy+p0.y+1) = 0xFFFFFF;
			}
		}/**/
	}
	
	if (1) {
		float scale = 6;
		float ox = peye->aScreen[p->Idx].View.Front.x, oy = peye->aScreen[p->Idx].View.Front.y;
		
		if (1) {
			ui ix = 0, iy = 0;
			for (iy = 1; iy < dEye_Screen_Cal_NUM; ++iy) {
				for (ix = 0; ix < dEye_Screen_Cal_NUM; ++ix) {
					S_Draw_Line_2d (
						ox+ scale*peye->aScreen[p->Idx].aaCal[iy-1][ix].P.x,
						oy+ scale*peye->aScreen[p->Idx].aaCal[iy-1][ix].P.y,
						ox+ scale*peye->aScreen[p->Idx].aaCal[iy][ix].P.x,
						oy+ scale*peye->aScreen[p->Idx].aaCal[iy][ix].P.y
					);
				}
			}
			for (iy = 0; iy < dEye_Screen_Cal_NUM; ++iy) {
				for (ix = 1; ix < dEye_Screen_Cal_NUM; ++ix) {
					S_Draw_Line_2d (
						ox+ scale*peye->aScreen[p->Idx].aaCal[iy][ix-1].P.x,
						oy+ scale*peye->aScreen[p->Idx].aaCal[iy][ix-1].P.y,
						ox+ scale*peye->aScreen[p->Idx].aaCal[iy][ix].P.x,
						oy+ scale*peye->aScreen[p->Idx].aaCal[iy][ix].P.y
					);
				}
			}
			for (iy = 0; iy < dEye_Screen_Cal_NUM-1; ++iy) {
				for (ix = 0; ix < dEye_Screen_Cal_NUM-1; ++ix) {
					S_Draw_Line_2d (
						ox+ scale*peye->aScreen[p->Idx].aaCal[iy][ix].P.x,
						oy+ scale*peye->aScreen[p->Idx].aaCal[iy][ix].P.y,
						ox+ scale*peye->aScreen[p->Idx].aaCal[iy+1][ix+1].P.x,
						oy+ scale*peye->aScreen[p->Idx].aaCal[iy+1][ix+1].P.y
					);
				}
			}
		}
		if (1) {
			float s, t;
			tV4f p0;
			Screen_Eye_Point (peye, &p0, &s, &t);
			
			if (ox+scale*p0.x-1 >= 0 && ox+scale*p0.x+1 < gM.pScreen->w
				&& oy+scale*p0.y-1 >= 0 && oy+scale*p0.y+1 < gM.pScreen->h
			) {
				dspix(ox+scale*p0.x-1, oy+scale*p0.y-1) = gColARGB;
				dspix(ox+scale*p0.x+0, oy+scale*p0.y-1) = gColARGB;
				dspix(ox+scale*p0.x+1, oy+scale*p0.y-1) = gColARGB;
				
				dspix(ox+scale*p0.x-1, oy+scale*p0.y+0) = gColARGB;
				dspix(ox+scale*p0.x,   oy+scale*p0.y) = 0x0000FF;
				dspix(ox+scale*p0.x+1, oy+scale*p0.y+0) = gColARGB;
				
				dspix(ox+scale*p0.x-1, oy+scale*p0.y+1) = gColARGB;
				dspix(ox+scale*p0.x+0, oy+scale*p0.y+1) = gColARGB;
				dspix(ox+scale*p0.x+1, oy+scale*p0.y+1) = gColARGB;
			}
		}/**/
	}
}

#undef dcal

void	Head_Eye_Print	(tHead* p, tEye* peye)
{
	tV4f eye, ret, vec;
	Head_Eye_Vector (p, peye, &ret);
	
	eye = peye->InHead.P;
	vec = ret;	V4f_sub_V4f (&vec, &peye->InHead.P);
	
	M4f_mul_V4f (&p->M4, &eye);
	M4f_mul_V4f (&p->M4, &ret);
	
	if (1) {
		float ox = peye->aScreen[0].View.Top.x, oy = peye->aScreen[0].View.Top.y;
		
		tV4f vv = ret;	V4f_sub_V4f (&vv, &eye);
		vv.x *= 40;
		vv.y *= 40;
		vv.z *= 40;
		
		if (peye == &gM.Left)
			gColARGB = 0x00FF00;
		else
			gColARGB = 0xFF0000;
		
		S_Draw_Line_2d (ox+eye.x, oy+eye.z, ox+eye.x+vv.x, oy+eye.z+vv.z);
		
	}
	if (1) {
		float ox = peye->aScreen[0].View.Left.x, oy = peye->aScreen[0].View.Left.y;
		
		tV4f vv = ret;	V4f_sub_V4f (&vv, &eye);
		vv.x *= 40;
		vv.y *= 40;
		vv.z *= 40;
		
		if (peye == &gM.Left)
			gColARGB = 0x00FF00;
		else
			gColARGB = 0xFF0000;
		
		S_Draw_Line_2d (ox+eye.z, oy+eye.y, ox+eye.z+vv.z, oy+eye.y+vv.y);
		
	}
}

void	Head_Eye_Gaze	(tHead* p, tEye* peye)
{
	tV4f eye, ret, vec;
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
	}
}

dyn_config gM_DC;
int	muhaha_Config_Thread	(void *data)
{
	dyn_config_watch (&gM_DC, "config.yaml");
	return 0;
}


int	muhaha_XEvent_Thread	(void *data)
{
	printf ("first\n");
	while (videoIn->signalquit) {
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
			#include "actions.h"
			#undef dact
			break;
		}
		}/**/
		
	//	printf ("hahahahha\n");
	//	XFlush(gM.X.pDisp);/**/
	}
}



void	muhaha_Keyboard_Init	()
{
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
				printf ("Key grab sym 0x%x code %d mods 0x%x\n", _keysym, XKeysymToKeycode(gM.X.pDisp, _keysym), mods[j]);/**/	\
				XGrabKey(gM.X.pDisp, XKeysymToKeycode(gM.X.pDisp, _keysym), mods[j], gM.aScreen[0].Win, True, GrabModeAsync, GrabModeAsync);		\
				XSync(gM.X.pDisp, False);	\
			}		\
		}while(0)
	
	#define dact(name,press_stuff)	dgrab(gM.Action.name);
	#include "actions.h"
	#undef dact
//	XGrabKey(gM.X.pDisp, gM.Action.Mod_Key, Mod2Mask, gM.aScreen[0].Win, True, GrabModeAsync, GrabModeAsync);		\
	
	#undef dgrab
	XSync(gM.X.pDisp, False);
}

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
				break;
			}
		}/**/
	//	gM.X.pDisp = XOpenDisplay(":0");
		if (gM.X.pDisp == NULL) {
			fprintf(stderr, "Couldn't open display\n");
			exit (1);
		}
	}
//	angle_test ();
	
	gM.Draw_X = 0;
	gM.Draw_Y = 0;
	gM.Draw_W = 800;
	gM.Draw_H = 600;
	
	
	Head_Init (&gM.Head);
	
	Eye_Init (&gM.Left);
	Eye_Init (&gM.Right);
	
	gM.pGaze_Mutex = SDL_CreateMutex();
	
	dyn_config_read(&gM_DC, "config.yaml");
	/*mythread = */SDL_CreateThread(muhaha_Config_Thread, (void *)NULL);
	
	
	M4f_Iden (&gM.World);
	
	
	M4f_trans (&gM.World, 0, 0, 0);
	
//	M4f_rotx (&gM.World, 15*deg2rad);
//	M4f_trans (&gM.World, 0, 0, -15);
//	M4f_rotx (&gM.World, 15*deg2rad);
	
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
	
	
	/*mythread = */SDL_CreateThread(muhaha_XEvent_Thread, (void *)NULL);
	
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
	++tfps;
	
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
				
			//	pin->Y = 0x80;
			//	pin->U = ((si)pin1->U - ((si)pin->U)) +7;
			//	pin->V = ((si)pin1->V - ((si)pin->V)) +8;
			//	pin->V = abs((si)pin1->V - (si)pin->V);
			//	pin->U = 0x7;
			//	pin->V = 0x8;
			//	dpix(x,y)->U = 0x7;
			//	dpix(x,y)->V = 0x8;
			}
		}/**/
	}
/*	printf ("World:\n");
	M4f_Print (&gM.World);
	printf ("Proj:\n");
	M4f_Print (&gM.Proj);/**/
	
	for (y = 0; y < 600; ++y)
		for (x = 800; x < 800+640; ++x)
			dspix(x, y) = 0;
	
	if (1) {
		memcpy(gM.pDst, videoIn->framebuffer, videoIn->width * (videoIn->height) * 2);
		
		Eye_Fit (&gM.Head.DotC);
		Eye_Clip (&gM.Head.DotC);
		
	//	memcpy(gM.pDst, videoIn->framebuffer, videoIn->width * (videoIn->height) * 2);
		
		Eye_Fit (&gM.Head.DotL);
		Eye_Clip (&gM.Head.DotL);
		
	//	memcpy(gM.pDst, videoIn->framebuffer, videoIn->width * (videoIn->height) * 2);
		
		Eye_Fit (&gM.Head.DotR);
		Eye_Clip (&gM.Head.DotR);
		
	//	memcpy(gM.pDst, videoIn->framebuffer, videoIn->width * (videoIn->height) * 2);
		
		Eye_Fit (&gM.Right);
		Eye_Clip (&gM.Right);
		
	//	memcpy(gM.pDst, videoIn->framebuffer, videoIn->width * (videoIn->height) * 2);
		
		Eye_Fit (&gM.Left);
		Eye_Clip (&gM.Left);
		
		Eye_Draw (&gM.Left);
		Eye_Draw (&gM.Right);
		printf ("Eye pos Left %f %f\n", gM.Left.P.x, gM.Left.P.y);
	//	printf ("Eye pos Righ %f %f\n", gM.Right.P.x, gM.Right.P.y);
		
		Eye_Draw (&gM.Head.DotC);
		Eye_Draw (&gM.Head.DotL);
		Eye_Draw (&gM.Head.DotR);
	}
	
	if (0) {
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
	}
	
//	dheadline (20,30,40,50);
	
	#define dprojline3(x0,y0,z0,x1,y1,z1)	\
		do {	\
			tV4f __p0 = {x0, y0, z0, 1}, __p1 = {x1, y1, z1, 1};	\
			V4f_DrawPosPos (&__p0, &__p1);	\
		}while(0)
	
//	dprojline2(0,0,100,100);
	if (0) {
		gCol.Y = 0xFF;
		gCol.U = 0x0;
		gCol.V = 0x0;
		static float ang = 0, zz = 13*3;
		
		float d = 13, z = 3*d;
		
	//	tM4f rot2 = {cosf(ang),sinf(ang),0,0, -sinf(ang),cosf(ang),0,0, 0,0,1,0, 0,0,0,1};
		tM4f rot = {1,0,0,0, 0,cosf(ang),sinf(ang),0, 0,-sinf(ang),cosf(ang),0, 0,0,0,1};
		tM4f rot2 = {cosf(ang),0,sinf(ang),0, 0,1,0,0, -sinf(ang),0,cosf(ang),0, 0,0,0,1};
		M4f_mul_M4f (&rot, &rot2);
		
		tM4f trans = {1,0,0,0, 0,1,0,0, 0,0,1,0, 0,0,zz,1};
		
		tV4f p00 = {-d, -d, -d, 1};
		tV4f p01 = {+d, -d, -d, 1};
		tV4f p02 = {+d, +d, -d, 1};
		tV4f p03 = {-d, +d, -d, 1};
		
		tV4f p10 = {-d, -d,  d, 1};
		tV4f p11 = {+d, -d,  d, 1};
		tV4f p12 = {+d, +d,  d, 1};
		tV4f p13 = {-d, +d,  d, 1};
		
		V4f_mul_M4f (&p00, &rot);
		V4f_mul_M4f (&p01, &rot);
		V4f_mul_M4f (&p02, &rot);
		V4f_mul_M4f (&p03, &rot);
		V4f_mul_M4f (&p10, &rot);
		V4f_mul_M4f (&p11, &rot);
		V4f_mul_M4f (&p12, &rot);
		V4f_mul_M4f (&p13, &rot);/**/
		
		V4f_mul_M4f (&p00, &trans);
		V4f_mul_M4f (&p01, &trans);
		V4f_mul_M4f (&p02, &trans);
		V4f_mul_M4f (&p03, &trans);
		V4f_mul_M4f (&p10, &trans);
		V4f_mul_M4f (&p11, &trans);
		V4f_mul_M4f (&p12, &trans);
		V4f_mul_M4f (&p13, &trans);
		
	//	V4f_Print (&p00);
	//	printf ("closer\n");
		V4f_DrawPosPos (&p00, &p01);
		V4f_DrawPosPos (&p01, &p02);
		V4f_DrawPosPos (&p02, &p03);
		V4f_DrawPosPos (&p03, &p00);
		
	//	printf ("mid\n");
		V4f_DrawPosPos (&p00, &p10);
		V4f_DrawPosPos (&p01, &p11);
		V4f_DrawPosPos (&p02, &p12);
		V4f_DrawPosPos (&p03, &p13);
		
	//	printf ("far\n");
		V4f_DrawPosPos (&p10, &p11);
		V4f_DrawPosPos (&p11, &p12);
		V4f_DrawPosPos (&p12, &p13);
		V4f_DrawPosPos (&p13, &p10);
		
		ang += 10*M_PI/180.0f;
		zz += 1.0f;
		if (zz >= z*2)
			zz = z/2;
	}/**/
	
	if (0) {
		gCol.Y = 0xFF;
		gCol.U = 0x0;
		gCol.V = 0x0;
		
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
		V4f_DrawPosPos (&p00, &p01);
		V4f_DrawPosPos (&p01, &p02);
		V4f_DrawPosPos (&p02, &p03);
		V4f_DrawPosPos (&p03, &p00);
		
	//	printf ("mid\n");
		V4f_DrawPosPos (&p00, &p10);
		V4f_DrawPosPos (&p01, &p11);
		V4f_DrawPosPos (&p02, &p12);
		V4f_DrawPosPos (&p03, &p13);
		
	//	printf ("far\n");
		V4f_DrawPosPos (&p10, &p11);
		V4f_DrawPosPos (&p11, &p12);
		V4f_DrawPosPos (&p12, &p13);
		V4f_DrawPosPos (&p13, &p10);
		
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
	
	
	if (1) {
		#define dline3(x0,y0,z0,x1,y1,z1)	\
			do {	\
				tV4f __p0 = {x0, y0, z0, 1}, __p1 = {x1, y1, z1, 1};	\
				V4f_mul_M4f (&__p0, &gM.Head.M4);	\
				V4f_mul_M4f (&__p1, &gM.Head.M4);	\
				V4f_DrawPosPos (&__p0, &__p1);	\
			}while(0)
		
		#undef dpoint3
		#define dpoint3(name,x0,y0,z0)	\
			tV4f name = {x0,y0,z0,1}; M4f_mul_V4f (&gM.Head.M4, &name);
		
		Head_Calc_M4_Rel (&gM.Head, &gM.HeadC);
		
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
		
		
		gCol.Y = 0xFF;
		gCol.U = 0x0;
		gCol.V = 0x0;
		
		{
			tV4f p0 = gM.Head.P; V4f_add_V4f (&p0, &gM.Head.N);
			V4f_DrawPosPos (&gM.Head.P, &p0);
		}
		
		printf ("c:\t%f\t%f\t%f\n", gM.Head.M4.x03, gM.Head.M4.x13, gM.Head.M4.x23);
	//	printf ("c:\t%f\t%f\t%f\n", gM.Head.M4_T.x03, gM.Head.M4_T.x13, gM.Head.M4_T.x23);
		
		if (0) {
			float d = 5;
			tV4f c = {0, 0, 0, 1};
			tV4f dx = {d, 0, 0, 1};
			tV4f dy = {0, d, 0, 1};
			tV4f dz = {0, 0, d, 1};
			
		//	V4f_mul_M4f (&c, &gM.Head.M4);
			M4f_mul_V4f (&gM.Head.M4, &c);
		//	printf (" c: ");	V4f_Print (&c);
		//	printf ("dx: ");	V4f_Print (&dx);
			
		//	V4f_mul_M4f (&dx, &gM.Head.M4);
		//	V4f_mul_M4f (&dy, &gM.Head.M4);
		//	V4f_mul_M4f (&dz, &gM.Head.M4);
			M4f_mul_V4f (&gM.Head.M4, &dx);
			M4f_mul_V4f (&gM.Head.M4, &dy);
			M4f_mul_V4f (&gM.Head.M4, &dz);
			
		//	printf ("closer\n");
			V4f_DrawPosPos (&c, &dx);
			V4f_DrawPosPos (&c, &dy);
			V4f_DrawPosPos (&c, &dz);
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
			tV4f pc = gM.Head.Mod.PC;
			tV4f pl = gM.Head.Mod.PL;
			tV4f pr = gM.Head.Mod.PR;/**/
		/*	printf ("pc: ");	V4f_Print (&pc);
			printf ("pl: ");	V4f_Print (&pl);
			printf ("pr: ");	V4f_Print (&pr);/**/
			
		/*	V4f_mul_M4f (&pc, &gM.Head.M4);
			V4f_mul_M4f (&pl, &gM.Head.M4);
			V4f_mul_M4f (&pr, &gM.Head.M4);/**/
			M4f_mul_V4f (&gM.Head.M4, &pc);
			M4f_mul_V4f (&gM.Head.M4, &pl);
			M4f_mul_V4f (&gM.Head.M4, &pr);
			
			V4f_DrawPosPos (&pl, &pc);
			V4f_DrawPosPos (&pr, &pc);
			V4f_DrawPosPos (&pl, &pr);
			if (0) {
				tV2f ret;
				V4f_ScreenPosNorm (&pl, &ret);	printf ("x %f y %f", ret.x, ret.y);
				V4f_ScreenPosNorm (&pr, &ret);	printf ("    x %f y %f\n", ret.x, ret.y);
			}
			if (0) {
				tV2si ret;
				V4f_ScreenPos (&pl, &ret);	printf ("x %ld y %ld", ret.x, ret.y);
				V4f_ScreenPos (&pr, &ret);	printf ("    x %ld y %ld\n", ret.x, ret.y);
			}
		}
		if (1) {
			float xdd = 10, ydd = 5;
		//	float xdd = 3.25, ydd = 2.25;
			tV4f p0 = {-xdd,	-ydd,	0,	1};
			tV4f p1 = {xdd,	-ydd,	0,	1};
			tV4f p2 = {xdd,	ydd,	0,	1};
			tV4f p3 = {-xdd,	ydd,	0,	1};
			
			M4f_mul_V4f (&gM.Head.M4, &p0);
			M4f_mul_V4f (&gM.Head.M4, &p1);
			M4f_mul_V4f (&gM.Head.M4, &p2);
			M4f_mul_V4f (&gM.Head.M4, &p3);
			
			V4f_DrawPosPos (&p0, &p1);
			V4f_DrawPosPos (&p1, &p2);
			V4f_DrawPosPos (&p2, &p3);
			V4f_DrawPosPos (&p3, &p0);
		}
		
		
		if (gM.bHead_Eye_LineDraw) {
			Head_Eye_LineDraw (&gM.Head, &gM.Left);
			Head_Eye_LineDraw (&gM.Head, &gM.Right);
		}
		
		if (0) {
			tV4f eye, ret, vec, re, rr, rv;
			Head_Eye_Vector (&gM.Head, &gM.Left, &ret);
			Head_Eye_Vector (&gM.Head, &gM.Right, &rr);
			
			eye = gM.Left.InHead.P;
			vec = ret;	V4f_sub_V4f (&vec, &gM.Left.InHead.P);
			
			re = gM.Right.InHead.P;
			rv = rr;	V4f_sub_V4f (&rv, &gM.Right.InHead.P);
			
		//	printf ("hoho2 ");	V4f_Print (&ret);
			
			M4f_mul_V4f (&gM.Head.M4, &eye);
			M4f_mul_V4f (&gM.Head.M4, &ret);
			
			M4f_mul_V4f (&gM.Head.M4, &re);
			M4f_mul_V4f (&gM.Head.M4, &rr);
			
		//	printf ("hoho1 ");	V4f_Print (&ret);
			
			gCol.U = 0x0;
			gCol.V = 0x0;
		///	V4f_DrawPosPos (&eye, &ret);
			
		//	printf ("vec ");	V4f_Print (&vec);
			
			ret = vec;
			
			float ax, ay;
			{
				M4f_mul_V4f (&gM.Head.M4_R, &vec);
				ax = atan2 (vec.z, vec.x);
				ay = atan2 (vec.z, vec.y);
				
				gM.L_Vec.x = ax;		gM.L_Vec.y = ay;
			//	printf ("axy  %f\t%f\n", ax, ay);
			}
			
			vec = ret;
			{
				M4f_mul_V4f (&gM.Head.M4I_R, &vec);
				ax = atan2 (vec.z, vec.x);
				ay = atan2 (vec.z, vec.y);
				
			//	gM.L_Vec.x = ax;		gM.L_Vec.y = ay;
			//	printf ("axyI %f\t%f\n", ax, ay);
			}
		//	gM.GazeL = Eye_map_point(&gM.Left, gM.L_Vec);
			
			
			if (1) {
				Head_Eye_Print (&gM.Head, &gM.Left);
				Head_Eye_Print (&gM.Head, &gM.Right);
			}
		}
		
		
		if (0) {//hihihi couldn't resist
			dpoint3(pbl, -1.8,	3.7,		3.0);
			dpoint3(pbr, 1.8,		3.7,		3.0);
			
			dpoint3(pt, 0,		-1.0,		2.2);
			
			dpoint3(pc, 0,		3.0,		6.0);
			
			V4f_DrawPosPos (&pbl, &pbr);
			V4f_DrawPosPos (&pbl, &pt);
			V4f_DrawPosPos (&pbr, &pt);
			
			V4f_DrawPosPos (&pbl, &pc);
			V4f_DrawPosPos (&pbr, &pc);
			V4f_DrawPosPos (&pt, &pc);
		}
		
	}
	
	if (0) {
		tM4f bproj = gM.Proj;
		tM4f bworld = gM.World;
		
		#undef dpoint3
		#define dpoint3(name,x0,y0,z0)	\
			tV4f name = {x0,y0,z0,1};//	V4f_mul_M4f (&name, &model)
		
		M4f_Iden (&gM.Proj);
		M4f_Iden (&gM.World);
		
	//	tV4f screen_c = {0, 10, 0, 0};
	//	tV4f head_c = {gM.Head.M4.x03, gM.Head.M4.x13, gM.Head.M4.x23, 0};
		tV4f head_c = {0, 10, 0, 0};
		tV4f screen_c = {gM.Head.M4.x03, gM.Head.M4.x13, gM.Head.M4.x23, 0};
		
		V4f_sub_V4f (&head_c, &screen_c);
		
		float ax = atan2(head_c.z,head_c.y), ay = atan2(head_c.z,head_c.x);
		
	//	M4f_rotx (&gM.World, M_PI_2 - ax);
	//	M4f_roty (&gM.World, -M_PI_2 + ay);
	//	head_c.x = 5;
	//	head_c.y = 2;
	//	gM.World = gM.Head.M4I;
	//	gM.World.x32 = 0;
	//	M4f_trans (&gM.World, head_c.x/11, 0, -20);
		
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
		
		M4f_trans (&gM.World, -head_c.x, -head_c.y, -head_c.z);
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
		gM.Proj_L = -1 - head_c.x;
		gM.Proj_R = 1 - head_c.x;
		gM.Proj_B = -1 - head_c.y;
		gM.Proj_T = 1 - head_c.y;/**/
		
		gM.Proj_N = head_c.z;
	//	gM.Proj_L /= -head_c.z;
	//	gM.Proj_R /= -head_c.z;
	//	gM.Proj_B /= -head_c.z;
	//	gM.Proj_T /= -head_c.z;
		
	/*	double fov;
		fov = gM.Cam.Image_FOV;
		{
			double aspectRatio = (float)gM.Cam.Image_W/(float)gM.Cam.Image_H, front = 1, back = 10;
			aspectRatio = 1;
			double tangent = tan(fov/2 * deg2rad);   // tangent of half fovY
			double height = front * tangent;          // half height of near plane
			double width = height * aspectRatio;      // half width of near plane
			
			// params: left, right, bottom, top, near, far
			gM.Proj_W = 2*width;
			gM.Proj_H = 2*height;
		//	gM.Proj_N = front;
		//	gM.Proj_F = back;
		}/**/
	//	gM.Proj_L = -gM.Proj_W/2 - head_c.x/11;
	//	gM.Proj_R = gM.Proj_W/2 - head_c.x/11;
	//	gM.Proj_B = -gM.Proj_H/2 - head_c.y/11;
	//	gM.Proj_T = gM.Proj_H/2 - head_c.y/11;
		
	//	gM.Proj_L /= head_c.z;
	//	gM.Proj_R /= head_c.z;
	//	gM.Proj_B /= head_c.z;
	//	gM.Proj_T /= head_c.z;
		{
		//	float l = -gM.Proj_W/4, r = 3*gM.Proj_W/4;
		//	float t = gM.Proj_W/2, b = -gM.Proj_W/2;
			
			gM.Proj.f[0][0] = 2*gM.Proj_N / (gM.Proj_R-gM.Proj_L);
			gM.Proj.f[1][1] = 2*gM.Proj_N / (gM.Proj_T-gM.Proj_B);
			
			gM.Proj.f[0][2] = (gM.Proj_R+gM.Proj_L) / (gM.Proj_R-gM.Proj_L);
			gM.Proj.f[1][2] = (gM.Proj_T+gM.Proj_B) / (gM.Proj_T-gM.Proj_B);
			
			gM.Proj.f[2][2] = -(gM.Proj_F+gM.Proj_N) / (gM.Proj_F - gM.Proj_N);
			
			gM.Proj.f[2][3] = (-2*gM.Proj_F*gM.Proj_N) / (gM.Proj_F - gM.Proj_N);
			gM.Proj.f[3][2] = -1;
		}/**/
		
	//	M4f_Print (&gM.World);
	//	M4f_Print (&gM.Proj);
		void draw (float x, float y, float z) {
			float dd = 0.2;
			dpoint3(pt, x+0,		y+-dd,	z);
			dpoint3(pb, x+0,		y+dd,		z);
			dpoint3(pl, x+-dd,	y+0,		z);
			dpoint3(pr, x+dd,		y+0,		z);
			V4f_DrawPosPos (&pl, &pt);
			V4f_DrawPosPos (&pt, &pr);
			V4f_DrawPosPos (&pr, &pb);
			V4f_DrawPosPos (&pb, &pl);
			
			dpoint3(pn, x,		y,		z);
			dpoint3(pf, x,		y,		-10);
			V4f_DrawPosPos (&pf, &pn);
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
			
			V4f_DrawPosPos (&pl, &pt);
			V4f_DrawPosPos (&pt, &pr);
			V4f_DrawPosPos (&pr, &pb);
			V4f_DrawPosPos (&pb, &pl);
		};
		
		for (z = -1; z <= 1; z += 0.1f) {
			{
				dpoint3(pn, z,		1,		0);
				dpoint3(pf, z,		1,		-10);
				V4f_DrawPosPos (&pf, &pn);
			}
			{
				dpoint3(pn, z,		-1,		0);
				dpoint3(pf, z,		-1,		-10);
				V4f_DrawPosPos (&pf, &pn);
			}
			{
				dpoint3(pn, 1,		z,		0);
				dpoint3(pf, 1,		z,		-10);
				V4f_DrawPosPos (&pf, &pn);
			}
			{
				dpoint3(pn, -1,		z,		0);
				dpoint3(pf, -1,		z,		-10);
				V4f_DrawPosPos (&pf, &pn);
			}
		}/**/
		
		gM.World = bworld;
		gM.Proj = bproj;
	}
	
	
	Screen_Eye_Print (gM.aScreen + 0, &gM.Left);
	Screen_Eye_Print (gM.aScreen + 0, &gM.Right);
	Screen_Eye_Print (gM.aScreen + 1, &gM.Left);
	Screen_Eye_Print (gM.aScreen + 1, &gM.Right);
	Screen_Eye_Print (gM.aScreen + 2, &gM.Left);
	Screen_Eye_Print (gM.aScreen + 2, &gM.Right);/**/
	
	Head_Eye_Gaze (&gM.Head, &gM.Left);
	Head_Eye_Gaze (&gM.Head, &gM.Right);
	
//	dprojline3(0, 0, 1, 400.0, 300.0, 1);
	
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
			}
		}
	}
	
//	point_clip (&gM.GazeL);
//	point_clip (&gM.GazeR);
//	point_clip (&gM.Gaze);
	
	gCol.U = 0xF;
	gCol.V = 0xF;
//	Vec_Draw (gM.Head.DotL.P.x, gM.Head.DotL.P.y, gM.Head.DotL.P.x+gM.L_Vec.x, gM.Head.DotL.P.y+gM.L_Vec.y);
//	Vec_Draw (gM.Head.DotR.P.x, gM.Head.DotR.P.y, gM.Head.DotR.P.x+gM.R_Vec.x, gM.Head.DotR.P.y+gM.R_Vec.y);
	
	Vec_Draw ((640/4*1)-gM.L_Vec.x, (480/2)-gM.L_Vec.y, (640/4*1), (480/2));
	Vec_Draw ((640/4*3)-gM.R_Vec.x, (480/2)-gM.R_Vec.y, (640/4*3), (480/2));
	
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
			printf ("micro dxy %f %f\n", dx, dy);
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
/*	{
		printf ("crap\n");
		XEvent e;
		XNextEvent(gM.X.pDisp, &e);
		printf ("yeah\n");
	}/**/
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
	
	si last_ms = SDL_GetTicks();
//	action_gui curr_action = A_VIDEO;
	while (videoIn->signalquit) {
		SDL_LockMutex(affmutex);
		
		SDL_GetMouseState(&x, &y);
		x = x*gM.pOverlay->w / 800;
		y = y*gM.pOverlay->h / 600;
		
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
					gM.Head.DotC.P.x = x;
					gM.Head.DotC.P.y = y;
					gM.Head.DotC.V.x = 0;
					gM.Head.DotC.V.y = 0;
					break;
				case SDLK_a:
					gM.Head.DotL.P.x = x;
					gM.Head.DotL.P.y = y;
					gM.Head.DotL.V.x = 0;
					gM.Head.DotL.V.y = 0;
					break;
				case SDLK_e:
					gM.Head.DotR.P.x = x;
					gM.Head.DotR.P.y = y;
					gM.Head.DotR.V.x = 0;
					gM.Head.DotR.V.y = 0;
					break;
				case SDLK_SEMICOLON:
					gM.Left.P.x = x;
					gM.Left.P.y = y;
					break;
				case SDLK_j:
					gM.Right.P.x = x;
					gM.Right.P.y = y;
					break;
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
					videoIn->signalquit = 0;
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
					SDL_mutexV (gM.pGaze_Mutex);
					break;
				}
				}
				break;
			case SDL_MOUSEMOTION:
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


void	S_Draw_Line_2d	(si x0, si y0, si x1, si y1)
{
//	printf ("xy0 %ld\t%ld\t%ld\t%ld\n", x0, y0, x1, y1);
	si dx = x1 - x0;
	si dy = y1 - y0;
	
	if (dx == 0 && dy == 0)
		return;
	
	ui dst_w = gM.pScreen->pitch/4;
	ui dst_h = gM.pScreen->h;
	u32 *pdst = gM.pScreen->pixels;
	
	//#define dspix(_x,_y) (*((u32*)gM.pScreen->pixels + ((si)(_x) + (si)(_y)*gM.pScreen->pitch/4)))
	
	u32 col = gColARGB;
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


