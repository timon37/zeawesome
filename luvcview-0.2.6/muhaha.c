
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

#define dSMALL_NUM  0.00000001

#define deg2rad	(M_PI/180.0f)
#define rad2deg	(180.0f/M_PI)
#define dpow2(_num) ((_num)*(_num))

#define dclip_l(val,min)	({typeof(val) ret; if ((val) < (min)) ret = (min); else ret = (val); ret;})
#define dclip_h(val,max)	({typeof(val) ret; if ((val) > (max)) ret = (max); else ret = (val); ret;})

#define ddot(x0,y0,x1,y1) ((x0)*(x1) + (y0)*(y1))

#define ddist2(x0,y0,x1,y1)	(((x1)-(x0))*((x1)-(x0)) + ((y1)-(y0))*((y1)-(y0)))

#define ddist(x0,y0,x1,y1) sqrt(ddist2(x0,y0,x1,y1))

#define dopix(_x,_y) ((tPix*)videoIn->framebuffer + ((si)(_x) + (si)(_y)*videoIn->width))
#define dnpix(_x,_y) ((tPix*)gM.pDst + ((si)(_x) + (si)(_y)*videoIn->width))
#define dpix(_x,_y) dnpix(_x,_y)


#define dsout(_x,_y) ((_x) < 0 || (_y) < 0 || (_x) >= gM.pScreen->w || (_y) >= gM.pScreen->h)
#define dspix(_x,_y) (*((u32*)gM.pScreen->pixels + ((si)(_x) + (si)(_y)*gM.pScreen->pitch/4)))

#define dmono2rgb(_g) ((u32)(_g) | (u32)(_g)<<8 | (u32)(_g)<<16)


#define dpixout(_x,_y) ((_x) < 0 || (_y) < 0 || (_x) >= videoIn->width || (_y) >= videoIn->height)


#define dset_cyuv(_x,_y,y,u,v)	\
	do {	\
		if (dpixout(_x,_y))	\
			break;		\
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
	V4f_mul_S (pv0, 1.0f/V4f_dist (pv0));
}

float	V4f_dist_V4f	(tV4f* pv0, tV4f* pv1)
{
	return sqrtf(dpow2(pv1->x-pv0->x) + dpow2(pv1->y-pv0->y) + dpow2(pv1->z-pv0->z));
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
		
		p0.x = -2*(ax - lx) / gM.Cam.Image_FOV_W;
		p0.y = -2*(ay - ly) / gM.Cam.Image_FOV_H;
		
	//	printf ("V4f_ScreenPosNorm: return %f %f\n", p0.x, p0.y); 
	//	static int count = 5;
	//	if (!--count)
	//		exit(1);
	}else {
		M4f_mul_V4f (&gM.World, &p0);
		M4f_mul_V4f (&gM.Proj, &p0);
		
		p0.x /= p0.w;
		p0.y /= p0.w;
		p0.z /= p0.w;
	}
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




void	Dbg_V4f_DrawPosPos	(tV4f* ppos0, tV4f* ppos1)
{
	tV4f p0 = *ppos0, p1 = *ppos1;
	
	if (p0.w != 1 || p1.w != 1)
		return;
	assert (p0.w == 1);
	
//	printf ("4f p0 ");	V4f_Print (&p0);
	
	M4f_mul_V4f (&gM.Dbg.World, &p0);		M4f_mul_V4f (&gM.Dbg.World, &p1);
//	printf ("4f p0 ");	V4f_Print (&p0);
	
	M4f_mul_V4f (&gM.Dbg.Proj, &p0);		M4f_mul_V4f (&gM.Dbg.Proj, &p1);
//	printf ("4f p0 ");	V4f_Print (&p0);
	
	p0.x /= p0.w;
	p0.y /= p0.w;
	p0.z /= p0.w;
	
	p1.x /= p1.w;
	p1.y /= p1.w;
	p1.z /= p1.w;
	
//	printf ("p0 %f %f  p1 %f %f\n", p0.x, p0.y, p1.x, p1.y);
	
	p0.x *= gM.Dbg.View_W/2;	p0.x += gM.Dbg.View_X;
	p0.y *= gM.Dbg.View_H/2;	p0.y += gM.Dbg.View_Y;
	p1.x *= gM.Dbg.View_W/2;	p1.x += gM.Dbg.View_X;
	p1.y *= gM.Dbg.View_H/2;	p1.y += gM.Dbg.View_Y;
	
//	printf ("p0 %f %f  p1 %f %f\n", p0.x, p0.y, p1.x, p1.y);
	S_Draw_Line_2d (p0.x, p0.y, p1.x, p1.y);
}

void	Dbg_V4f_DrawPoint		(tV4f* ppos0)
{
	tV4f p0 = *ppos0;
	
	if (p0.w != 1)
		return;
	assert (p0.w == 1);
	
//	printf ("4f p0 ");	V4f_Print (&p0);
	
	M4f_mul_V4f (&gM.Dbg.World, &p0);
//	printf ("4f p0 ");	V4f_Print (&p0);
	
	M4f_mul_V4f (&gM.Dbg.Proj, &p0);
//	printf ("4f p0 ");	V4f_Print (&p0);
	
	p0.x /= p0.w;
	p0.y /= p0.w;
	p0.z /= p0.w;
	
	p0.x *= gM.Dbg.View_W/2;	p0.x += gM.Dbg.View_X;
	p0.y *= gM.Dbg.View_H/2;	p0.y += gM.Dbg.View_Y;
	
//	printf ("p0 %f %f  p1 %f %f\n", p0.x, p0.y, p1.x, p1.y);
//	S_Draw_Line_2d (p0.x, p0.y, p1.x, p1.y);
	
	si x, y;
	#define dd 1
	for (y = p0.y-dd; y <= p0.y+dd; ++y) {
		for (x = p0.x-dd; x <= p0.x+dd; ++x) {
			if (x < 0 || x >= gM.pScreen->w || y < 0 || y >= gM.pScreen->h)
				continue;
			dspix(x,y) = gColARGB;
		}
	}
	#undef dd
}



void	M4f_Frustrum		(tM4f* pproj, float w, float h, float n, float f)
{
	pproj->x00 = n / (w/2);
	pproj->x11 = n / (h/2);
	
	pproj->x22 = -(f+n) / (f-n);
	
	pproj->x23 = (-2*f*n) / (f-n);
	pproj->x32 = -1;
	
}
void	M4f_Ortho		(tM4f* pproj, float w, float h, float n, float f)
{
	pproj->x00 = 1 / (w/2);
	pproj->x11 = 1 / (h/2);
	
	pproj->x22 = -2 / (f-n);
	pproj->x23 = -(f+n) / (f-n);
	
	pproj->x33 = 1;
	
}
void	Proj_Cam		()
{
	if (0) {
		float fd = dpow2(gM.Cam.Full_W) + dpow2(gM.Cam.Full_H);
		float id = dpow2(gM.Cam.Image_Zoom*gM.Cam.Image_W) + dpow2(gM.Cam.Image_Zoom*gM.Cam.Image_H);
	//	float id = ddist(gM.Cam.Image_W, gM.Cam.Image_H);
		gM.Cam.Image_FOV = id/fd * gM.Cam.Full_FOV;
	}
	if (1) {
		float fd = gM.Cam.Full_W;
		float id = gM.Cam.Image_Zoom*gM.Cam.Image_W;
	//	float id = ddist(gM.Cam.Image_W, gM.Cam.Image_H);
		gM.Cam.Image_FOV = id/fd * gM.Cam.Full_FOV;
	}
	
	float fov = gM.Cam.Image_FOV*deg2rad;
	float a = atan2(gM.Cam.Image_H, gM.Cam.Image_W);
	
	gM.Cam.Image_FOV_W = fov*cos(a);
	gM.Cam.Image_FOV_H = fov*sin(a);
	
	printf ("Image_FOV %f  W %f H %f\n", gM.Cam.Image_FOV, gM.Cam.Image_FOV_W*rad2deg, gM.Cam.Image_FOV_H*rad2deg);
	
	gM.Proj_W = gM.Proj_N*tan(gM.Cam.Image_FOV_W/2.0f);
	gM.Proj_H = gM.Proj_N*tan(gM.Cam.Image_FOV_H/2.0f);
	
	gM.Proj_L = -gM.Proj_W;
	gM.Proj_R = gM.Proj_W;
	gM.Proj_B = -gM.Proj_H;
	gM.Proj_T = gM.Proj_H;
	
	gM.Proj_W *= 2;
	gM.Proj_H *= 2;
	
/*	gM.Proj.x00 = gM.Proj_N / (gM.Proj_W/2);
	gM.Proj.x11 = gM.Proj_N / (gM.Proj_H/2);
	
	gM.Proj.x22 = -(gM.Proj_F+gM.Proj_N) / (gM.Proj_F - gM.Proj_N);
	
	gM.Proj.x23 = (-2*gM.Proj_F*gM.Proj_N) / (gM.Proj_F - gM.Proj_N);
	gM.Proj.x32 = -1;/**/
	M4f_Frustrum (&gM.Proj, gM.Proj_W, gM.Proj_H, gM.Proj_N, gM.Proj_F);
}

void	Dbg_Ortho	(struct sDbg_View* pv)
{
	gM.Dbg.View_X = pv->Off.x;
	gM.Dbg.View_Y = pv->Off.y;
	gM.Dbg.View_W = 100;
	gM.Dbg.View_H = 100;
	
	M4f_Iden (&gM.Dbg.Proj);
	
	M4f_Ortho (&gM.Dbg.Proj, pv->Scale*gM.Dbg.View_W, pv->Scale*gM.Dbg.View_H, 0.1f, 100);
//	M4f_Frustrum (&gM.Dbg.Proj, gM.Dbg.View_W, gM.Dbg.View_H, 0.1f, 100);
	
	M4f_Iden (&gM.Dbg.World);
	
	M4f_trans (&gM.Dbg.World, pv->T_X, pv->T_Y, -30);
	
	M4f_rotx (&gM.Dbg.World, pv->R_X*deg2rad);
	M4f_roty (&gM.Dbg.World, pv->R_Y*deg2rad);
	
//	M4f_rotx (&gM.Dbg.World, 15*deg2rad);
//	M4f_trans (&gM.Dbg.World, 0, 0, -15);
//	M4f_rotx (&gM.Dbg.World, 15*deg2rad);
}

void	Dbg_Ortho_Front	()
{
	Dbg_Ortho (&gM.Dbg.Front);
//	M4f_rotx (&gM.Dbg.World, 15*deg2rad);
//	M4f_trans (&gM.Dbg.World, 0, 0, -15);
//	M4f_rotx (&gM.Dbg.World, 15*deg2rad);
	
}

void	Dbg_Ortho_Top	()
{
	Dbg_Ortho (&gM.Dbg.Top);
	
//	M4f_rotx (&gM.Dbg.World, 90*deg2rad);
	
}
void	Dbg_Ortho_Left	()
{
	Dbg_Ortho (&gM.Dbg.Left);
	
//	M4f_roty (&gM.Dbg.World, 90*deg2rad);
	
}


void	Dbg_V4f_ADrawPosPos	(tV4f* ppos0, tV4f* ppos1)
{
	Dbg_Ortho_Front ();
	Dbg_V4f_DrawPosPos (ppos0, ppos1);
	Dbg_Ortho_Top ();
	Dbg_V4f_DrawPosPos (ppos0, ppos1);
	Dbg_Ortho_Left ();
	Dbg_V4f_DrawPosPos (ppos0, ppos1);
}

void	Dbg_V4f_ADrawPoint	(tV4f* ppos0)
{
	Dbg_Ortho_Front ();
	Dbg_V4f_DrawPoint (ppos0);
	Dbg_Ortho_Top ();
	Dbg_V4f_DrawPoint (ppos0);
	Dbg_Ortho_Left ();
	Dbg_V4f_DrawPoint (ppos0);
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



void	Col_EyeSet	(tEye* peye)
{
	if (peye == &gM.Left)
		gColARGB = 0x00FF00;
	else
		gColARGB = 0xFF0000;
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



void	Eye_Init	(tEye* peye)
{
	printf ("INIT  %lx\n", peye);
	peye->Point_N = 0;
	peye->Point_Max = 0;
	peye->paPoint = 0;
	
	peye->InHead.Line_N = 0;
	memset (peye->InHead.aLine, 0, sizeof(peye->InHead.aLine));
	
	
	peye->FF.paMark = 0;
	
	
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
	peye->FF.paMark = realloc (peye->FF.paMark, dpow2(peye->FF.Max_R) * sizeof(peye->FF.paMark[0]));
	
	
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

void	Eye_CalcAYUV	(tEye* peye, si dd)
{
	si x, y;
//	si ay = 0, au = 0, av = 0, n = 0;
	ay = 0;
	au = 0;
	av = 0;
	n = 0;
	
	for (y = peye->P.y-dd; y <= peye->P.y+dd; ++y) {
		for (x = peye->P.x-dd; x <= peye->P.x+dd; ++x) {
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

void	Eye_V_Pre	(tEye* peye)	//motion compensation
{
	return;
	peye->OP = peye->P;
	peye->P.x += peye->V.x;
	peye->P.y += peye->V.y;
}
void	Eye_V_Post	(tEye* peye)
{
	return;
//	peye->V.x = peye->P.x - oldpos.x;
//	peye->V.y = peye->P.y - oldpos.y;
	peye->V.x = (peye->V.x + peye->P.x - peye->OP.x) * 0.5f;
	peye->V.y = (peye->V.y + peye->P.y - peye->OP.y) * 0.5f;
	
	peye->V.x *= 0.9f;
	peye->V.y *= 0.9f;
}


void	Eye_GV_Calc		(tEye* peye)
{
	SDL_mutexP (peye->GV_mutex);
	if (gM.Eye_GlintMode == 2) {
		tV2f v0 = peye->P;	V2f_sub_V2f (&v0, &peye->G0);
		tV2f v1 = peye->P;	V2f_sub_V2f (&v1, &peye->G1);
		tV2f v01 = peye->G1;	V2f_sub_V2f (&v01, &peye->G0);
		
		V2f_add_V2f (&v0, &v1);
		V2f_mul_S (&v0, 1/V2f_dist (&v01));
		peye->GV = v0;
	}else if (gM.Eye_GlintMode == 1) {
		tV2f v0 = peye->P;	V2f_sub_V2f (&v0, &peye->G0);
		
		tV2f gtf = peye->G0;
		gtf.x -= gM.Cam.Image_W/2;
		gtf.y -= gM.Cam.Image_H/2;
		
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
	for (i = 0; i < 200; ++i) {
		peye->GP = sr_norm;
		V4f_mul_S (&peye->GP, z);
		
		tV4f l0_v, l1_v;
		Sphere_Reflect (&l0_v, &peye->GP, peye->LR, &gM.aLight[0].P);
		Sphere_Reflect (&l1_v, &peye->GP, peye->LR, &gM.aLight[1].P);
		V4f_norm (&l0_v);
		V4f_norm (&l1_v);
		
		float ld = V4f_dot_V4f (&l0_v, &l1_v);
	//	printf ("z %f	gd ld	%f %f\n", z, gd, ld);
		
		z += 1000*(gd - ld);
		
	/*	if (ld > gd) {
			z -= 0.1;
		}else {
			z += 0.1;
		}/**/
	//	printf ("g0_v "); V4f_Print (&g0_v);
	//	printf ("l0_v "); V4f_Print (&l0_v);
	//	printf ("g1_v "); V4f_Print (&g1_v);
	//	printf ("l1_v "); V4f_Print (&l1_v);
	}
	
	
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
	float x = peye->P.x, y = peye->P.y;
	si n = 0;
	
	for (; ddist2(x,y,peye->P.x,peye->P.y) < peye->Exp_R*peye->Exp_R; ) {
		if (peye->LinView.x != 0)
			dspix(peye->LinView.x+n, peye->LinView.y+i) = dnpix(x,y)->Y | dnpix(x,y)->Y<<8 | dnpix(x,y)->Y<<16;
		
		if (	dopix(x,y)->Y > peye->Pix_Bright
		) {
		//	return 1;
			float x1 = x, y1 = y;
			for (; ddist2(x1,y1,peye->P.x,peye->P.y) < peye->Exp_R*peye->Exp_R*4; ) {
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

si	Eye_S3Fit_EdgeMark4	(tEye* peye, si i, float t, float* px, float* py)
{
	u08 wrote = 0;
	float x = peye->P.x, y = peye->P.y;
	si n = 0;
	
	for (; ddist2(x,y,peye->P.x,peye->P.y) < dpow2(peye->Max_R); ) {
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
		
		Eye_CirView_Point_xy (peye, peye->P.x - x0, peye->P.y - y0, 0x00FF00);
		Eye_CirView_Point_xy (peye, peye->P.x - x1, peye->P.y - y1, 0x00FF00);
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

#define dGlint_FID 16
#define dGlint_Num 4

#define dMarkOut(_x,_y)	((_x) < peye->FF.Mark_P.x || (_x) >= peye->FF.Mark_P.x+peye->FF.Max_R || (_y) < peye->FF.Mark_P.y || (_y) >= peye->FF.Mark_P.y+peye->FF.Max_R)

#define dMarkRaw(_x,_y)	(*(peye->FF.paMark + (_x) + (_y)*(peye->FF.Max_R)))
#define dMark(_x,_y)	dMarkRaw((_x)-peye->FF.Mark_P.x, (_y)-peye->FF.Mark_P.y)


void	Eye_FF_Mark_Crap		(tEye* peye, si x, si y)
{
	if (dpixout(x,y))
		return;
	if (dMarkOut(x,y))
		return;
	if (dMark(x,y))
		return;
	
	if (dopix(x,y)->Y < peye->FF.Y) {
		dMark(x,y) = peye->FF.tmpID;
		dnpix(x,y)->U = 0;
		dnpix(x,y)->V = 0;
		++peye->FF.tmpNum;
		peye->FF.tmpP.x += x;
		peye->FF.tmpP.y += y;
		Eye_CirView_Point_xy (peye, x-peye->FF.Mark_P.x, y-peye->FF.Mark_P.y, gColARGB);
		Eye_FF_Mark_Crap (peye, x-1, y);
		Eye_FF_Mark_Crap (peye, x+1, y);
		Eye_FF_Mark_Crap (peye, x, y-1);
		Eye_FF_Mark_Crap (peye, x, y+1);
	}
}

si	Eye_FF_G_Mark_Crap	(tEye* peye, si x, si y)
{
	if (dpixout(x,y))
		return;
	if (dMarkOut(x,y))
		return;
	if (dMark(x,y))
		return;
	
	if (dopix(x,y)->Y > peye->FF.GY) {
		dMark(x,y) = peye->FF.tmpID;
		dnpix(x,y)->U = 0xF;
		dnpix(x,y)->V = 0xF;
		++peye->FF.tmpNum;
		peye->FF.tmpP.x += x;
		peye->FF.tmpP.y += y;
		Eye_CirView_Point_xy (peye, x-peye->FF.Mark_P.x, y-peye->FF.Mark_P.y, gColARGB);
		Eye_FF_G_Mark_Crap (peye, x-1, y);
		Eye_FF_G_Mark_Crap (peye, x+1, y);
		Eye_FF_G_Mark_Crap (peye, x, y-1);
		Eye_FF_G_Mark_Crap (peye, x, y+1);
	}
}

si	Eye_FF_Mark2Pos		(tEye* peye, u08 id, tV2f *pret)
{
	si ax = 0, ay = 0;
	si num = 0;
	si x, y;
	for (y = 0; y < peye->FF.Max_R; ++y) {
		for (x = 0; x < peye->FF.Max_R; ++x) {
			if (	dMarkRaw(x,y) == id) {
				ax += x;
				ay += y;
				Eye_CirView_Point_xy (peye, x, y, gColARGB);
				++num;
			}
		}
	}/**/
	pret->x = (float)ax / num;
	pret->y = (float)ay / num;
	return num;
}

float	Eye_FF_Mark2AvgDist	(tEye* peye, float ox, float oy, u08 id)
{
	float adist = 0;
	si num = 0;
	si x, y;
	for (y = 0; y < peye->FF.Max_R; ++y) {
		for (x = 0; x < peye->FF.Max_R; ++x) {
			if (dMarkRaw(x,y) == id) {
				adist += ddist(x,y,ox,oy);
				++num;
			}
		}
	}/**/
	return adist / num;
}

float	Eye_FF_Mark2Confidence	(tEye* peye, float ox, float oy, u08 id)
{
	float adist = 0;
	si num = 0;
	si x, y;
	for (y = 0; y < peye->FF.Max_R; ++y) {
		for (x = 0; x < peye->FF.Max_R; ++x) {
			if (dMarkRaw(x,y) == id) {
				if (ddist2(x,y,ox,oy) <= dpow2(peye->FF.Perf_R)) {
					++adist;
				}
				++num;
			}else {
				if (ddist2(x,y,ox,oy) <= dpow2(peye->FF.Perf_R)) {
					--adist;
				}
			}
		}
	}/**/
	return adist / num;
}


void	Eye_FF		(tEye* peye)
{
	ay = peye->FF.Y;
	au = 0;
	av = 0;
	
	memset (peye->FF.paMark, 0, dpow2(peye->FF.Max_R) * sizeof(peye->FF.paMark[0]));
	
	peye->FF.Mark_P.x = (peye->P.x + 0.5f) - peye->FF.Max_R/2;
	peye->FF.Mark_P.y = (peye->P.y + 0.5f) - peye->FF.Max_R/2;
	
	gColARGB = 0x00FF00;
	peye->FF.tmpID = 1;
	if (peye->FF.Perf_R != 0) {
		float best_ar = NAN;
	//	printf ("Search for main %f\n", peye->FF.Perf_R);
		float best_x = peye->P.x, best_y = peye->P.y;
		
		si ox = peye->P.x, oy = peye->P.y, x, y;
		si search_r = peye->FF.Search_R;
		for (y = oy - search_r; y < oy + search_r; y += 10) {
			for (x = ox - search_r; x < ox + search_r; x += 10) {
				peye->FF.tmpNum = 0;
				peye->FF.tmpP.x = 0;
				peye->FF.tmpP.y = 0;
				Eye_FF_Mark_Crap (peye, x, y);
				peye->FF.tmpP.x /= peye->FF.tmpNum;
				peye->FF.tmpP.y /= peye->FF.tmpNum;
				if (peye->FF.tmpNum) {
					#if 0
					float ar = Eye_FF_Mark2AvgDist (peye, peye->FF.tmpP.x-peye->FF.Mark_P.x, peye->FF.tmpP.y-peye->FF.Mark_P.y, peye->FF.tmpID);
					
				//	printf ("AvgDist	%d	%f\n", peye->FF.tmpID, ar);
					if (isfinite(ar)
						&& fabsf(peye->FF.Perf_R - ar) <= peye->FF.MaxDiff_R
						&& (!isfinite(best_ar) 
							|| fabsf(peye->FF.Perf_R - ar) < fabsf(peye->FF.Perf_R - best_ar))
					) {
						best_ar = ar;
						best_x = peye->FF.tmpP.x;
						best_y = peye->FF.tmpP.y;
					}
					#else
					float ar = Eye_FF_Mark2Confidence (peye, peye->FF.tmpP.x-peye->FF.Mark_P.x, peye->FF.tmpP.y-peye->FF.Mark_P.y, peye->FF.tmpID);
					
				//	printf ("Confidence	%d	%f\n", peye->FF.tmpID, ar);
					if (ar >= peye->FF.MaxDiff_R
						&& (!isfinite(best_ar) 
							|| (ar > best_ar))
					) {
						best_ar = ar;
						best_x = peye->FF.tmpP.x;
						best_y = peye->FF.tmpP.y;
					}
					
					#endif
					peye->FF.tmpID++;
				}
			}
		}
	//	printf ("best_ar %f\n", best_ar);
		peye->P.x = best_x;
		peye->P.y = best_y;
		
		peye->Ax = peye->FF.Perf_R;
		peye->Ay = peye->FF.Perf_R;
	}else {
		peye->FF.tmpNum = 0;
		peye->FF.tmpP.x = 0;
		peye->FF.tmpP.y = 0;
		Eye_FF_Mark_Crap (peye, peye->P.x, peye->P.y);
		
		if (peye->FF.tmpNum > 20) {
			peye->P.x = peye->FF.tmpP.x / peye->FF.tmpNum;
			peye->P.y = peye->FF.tmpP.y / peye->FF.tmpNum;
		}/**/
	/*	tV2f pos;
		if (Eye_FF_Mark2Pos (peye, 1, &pos) > 20) {
			peye->P.x = pos.x + peye->FF.Mark_P.x;
			peye->P.y = pos.y + peye->FF.Mark_P.y;
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
		si search_r = peye->FF.GSearch_R;
		gColARGB = 0xFF;
		si x, y;
		for (y = peye->P.y; y < peye->P.y + search_r; ++y) {
			for (x = peye->P.x - search_r; x < peye->P.x + search_r; ++x) {
				gColARGB = g[glint_idx].col;
				peye->FF.tmpNum = 0;
				peye->FF.tmpP.x = 0;
				peye->FF.tmpP.y = 0;
				peye->FF.tmpID = dGlint_FID+glint_idx;
				Eye_FF_G_Mark_Crap (peye, x, y);
				g[glint_idx].num = peye->FF.tmpNum;
				if (g[glint_idx].num) {
					g[glint_idx].pos.x = peye->FF.tmpP.x / peye->FF.tmpNum;
					g[glint_idx].pos.y = peye->FF.tmpP.y / peye->FF.tmpNum;
					
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
	
	Eye_CalcAYUV (peye, 15);
	Eye_S5 (peye);
	Eye_CalcAYUV (peye, 15);
	Eye_S3Fit (peye);
	
}

void	Eye_Fit		(tEye* peye)
{
	switch (peye->Fit) {
	case eEye_Fit_S0:
		Eye_CalcAYUV (peye, 12);
		Eye_S0 (peye);
		break;
	case eEye_Fit_S2Fit:
		Eye_OldFF (peye);
		Eye_EdgeMark (peye);
		Eye_S2Fit (peye);
		break;
	case eEye_Fit_S3Fit_START ... eEye_Fit_S3Fit_END:
	//	Eye_CalcAYUV (peye, 12);
		ay = peye->Pix_Bright;
		Eye_S3Fit (peye);
		break;
	case eEye_Fit_S4Fit_START ... eEye_Fit_S4Fit_END:
		Eye_CalcAYUV (peye, 12);
		Eye_S4 (peye);
		break;
	case eEye_Fit_S5_START ... eEye_Fit_S5_END:
		Eye_CalcAYUV (peye, 12);
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
	#define dmarg 10
	if (peye->P.x < dmarg)
		peye->P.x = dmarg;
	if (peye->P.x > videoIn->width-dmarg)
		peye->P.x = videoIn->width-dmarg;
	if (peye->P.y < dmarg)
		peye->P.y = dmarg;
	if (peye->P.y > videoIn->height-dmarg)
		peye->P.y = videoIn->height-dmarg;
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



tV2f	Eye_map_point	(tEye* peye, tV2f p)
{
	return map_point (&peye->Homo, p);
}



void	Cam_Pos2Ray	(tV2f pos, tV4f* pp0, tV4f* pp1, tV4f* pvec)
{
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
		p1.x = gM.Proj_L + (pos.x/gM.Cam.Image_W)*gM.Proj_W;
		p1.y = gM.Proj_B + (pos.y/gM.Cam.Image_H)*gM.Proj_H;
	}
	if (pp0)
		*pp0 = p0;
	if (pp1)
		*pp1 = p1;
	if (pvec) {
		*pvec = p1;
	}
}


void	Cam_Retina_Ray	(tEye* peye, tV4f* pp0, tV4f* pp1, tV4f* pvec)
{
	Cam_Pos2Ray (peye->P, pp0, pp1, pvec);
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
		
		if ((value = v4l2SetControl(videoIn, V4L2_CID_FOCUS_ABSOLUTE, 1)) < 0)
			printf("Set CT_FOCUS_ABSOLUTE_CONTROL to %ld error\n", value);
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



void	Head_Init		(tHead* p)
{
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
		
		
		static tM4f rot = {1,0,0,0, 0,1,0,0, 0,0,1,0, 0,0,0,1};
		static tM4f trans = {1,0,0,0, 0,1,0,0, 0,0,1,0, 0,0,0,1};
		
		if (1) {
			p->P.x = p->P.y = 0;
			p->P.z = -30;
			p->R_X = p->R_Y = p->R_Z = 0;
			M4f_Iden (&rot);	M4f_Iden (&trans);	M4f_trans (&trans, 0, 0, -30);
		}
	/*	if (!finite(p->M4.x03)) {
			p->P.x = p->P.y = 0;
			p->P.z = -30;
			p->R_X = p->R_Y = p->R_Z = 0;
		}/**/
		
		float err = 0;
		static si i = 0;
		for (i = 0; i < 100; ++i) {
			if (1) {
			//	M4f_Iden (&trans);
			//	M4f_trans (&trans, p->P.x, p->P.y, p->P.z);
				
				M4f_Iden (&rot);
				M4f_rotx (&rot, p->R_X);
				M4f_roty (&rot, p->R_Y);
				M4f_rotz (&rot, p->R_Z);
			}
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
			
		/*	p->P.x += p->Mod.TInc.x*dx;
			p->P.y += p->Mod.TInc.y*dy;
			p->P.z += p->Mod.TInc.z*dz;
			p->P.w = 1;/**/
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
		//	M4f_trans (&trans, p->P.x, p->P.y, p->P.z);
			
			M4f_Iden (&rot);
			M4f_rotx (&rot, p->R_X);
			M4f_roty (&rot, p->R_Y);
			M4f_rotz (&rot, p->R_Z);
		}
	/*	if (fabsf(err) > 5 || !finite(trans.x03)) {
			M4f_Iden (&rot);	M4f_Iden (&trans);	M4f_trans (&trans, 0, 0, -40);
			p->P.x = p->P.y = 0;
			p->P.z = -40;
			p->R_X = p->R_Y = p->R_Z = 0;
		}/**/
		p->M4_T = trans;
		
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
		p->M4_R = rot;
		
		M4f_Iden (&p->M4);
		M4f_mul_M4f (&p->M4, &p->M4_R);
		M4f_mul_M4f (&p->M4, &p->M4_T);
		
		p->P.x = p->M4.x03;
		p->P.y = p->M4.x13;
		p->P.z = p->M4.x23;
		p->P.w = 1;
		
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
			
			Dbg_V4f_ADrawPosPos (&p0, &p1);
		}
	}
	gCol.U = 0xF;
	gCol.V = 0xF;
	tV4f pos = peye->InHead.P;
	M4f_mul_V4f (&p->M4, &pos);
	V4f_DrawPosPos (&pos, &pos);
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
				Dbg_V4f_DrawPosPos (&p0, &p1);
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
				Dbg_V4f_DrawPosPos (&p0, &p1);
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
				Dbg_V4f_DrawPosPos (&p0, &p1);
			}
		}
	}
	if (1) {
		float s, t;
		tV4f p0;
		Screen_Eye_Point (peye, &p0, &s, &t);
		
		Dbg_V4f_DrawPoint (&p0);
	}/**/
}

void	Screen_Eye_Print	(tScreen* p, tEye* peye)
{
//	printf ("Screen_Eye_Print Idx %d\n", p->Idx);
	Col_EyeSet (peye);
	
	if (1) {
		Dbg_Ortho_Front ();
		Screen_Eye_PrintSub (p, peye);
		Dbg_Ortho_Top ();
		Screen_Eye_PrintSub (p, peye);
		Dbg_Ortho_Left ();
		Screen_Eye_PrintSub (p, peye);
	}
	if (1) {
		Dbg_Ortho_Front ();
		
		ui ix = 0, iy = 0;
		for (iy = 0; iy < dEye_Screen_Cal_NUM-1; ++iy) {
			for (ix = 0; ix < dEye_Screen_Cal_NUM-1; ++ix) {
				tV4f p0 = dcal(ix,iy-1).P;
				tV4f p1 = dcal(ix,iy).P;
				Dbg_V4f_DrawPosPos (&p0, &p1);
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
	//	printf ("width top %f\n", V4f_dist_V4f (&dcal(0,0).P, &dcal(dEye_Screen_Cal_LAST,0).P));
	//	printf ("width bot %f\n", V4f_dist_V4f (&dcal(0,dEye_Screen_Cal_LAST).P, &dcal(dEye_Screen_Cal_LAST,dEye_Screen_Cal_LAST).P));
		
		
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
		
		Dbg_Ortho_Front ();
		Dbg_V4f_DrawPosPos (&eye, &vec);
		Dbg_V4f_DrawPoint (&eye);
		
		Dbg_Ortho_Top ();
		Dbg_V4f_DrawPosPos (&eye, &vec);
		Dbg_V4f_DrawPoint (&eye);
		
		Dbg_Ortho_Left ();
		Dbg_V4f_DrawPosPos (&eye, &vec);
		Dbg_V4f_DrawPoint (&eye);
		
		if (1) {
			Head_Eye_VectorGlob (p, peye, &eye, NULL, &vec);
			V4f_mul_S (&vec, peye->InHead.R/V4f_dist(&vec));
			V4f_add_V4f (&vec, &eye);
			
			Dbg_V4f_ADrawPoint (&vec);
			
			gCol.U = 0x0;
			gCol.V = 0x0;
			V4f_DrawPosPos (&vec, &vec);
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
		
		
		Dbg_Ortho_Front ();
		Dbg_V4f_DrawPosPos (&gM.Head.P, &p0);
		Dbg_V4f_DrawPosPos (&pl, &pc);
		Dbg_V4f_DrawPosPos (&pr, &pc);
		Dbg_V4f_DrawPosPos (&pl, &pr);
		
		Dbg_Ortho_Top ();
		Dbg_V4f_DrawPosPos (&gM.Head.P, &p0);
		Dbg_V4f_DrawPosPos (&pl, &pc);
		Dbg_V4f_DrawPosPos (&pr, &pc);
		Dbg_V4f_DrawPosPos (&pl, &pr);
		
		Dbg_Ortho_Left ();
		Dbg_V4f_DrawPosPos (&gM.Head.P, &p0);
		Dbg_V4f_DrawPosPos (&pl, &pc);
		Dbg_V4f_DrawPosPos (&pr, &pc);
		Dbg_V4f_DrawPosPos (&pl, &pr);
		
	}
	Head_Eye_Print (p, &gM.Left);
	Head_Eye_Print (p, &gM.Right);
	
	if (1) {
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
			Dbg_V4f_DrawPosPos (&eye, &vec);
			Dbg_V4f_DrawPoint (&eye);
			
			Dbg_Ortho_Top ();
			Dbg_V4f_DrawPosPos (&eye, &vec);
			Dbg_V4f_DrawPoint (&eye);
			
			Dbg_Ortho_Left ();
			Dbg_V4f_DrawPosPos (&eye, &vec);
			Dbg_V4f_DrawPoint (&eye);
		}
	}
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
//	REvDev_Init (tREvDev* p);
	
//	XkbSetDetectableAutorepeat (gM.X.pDisp, True, NULL);
	
//	angle_test ();
	gM.X.Queue_N = 0;
	gM.X.Queue_C = 0;
	
	gM.Draw_X = 0;
	gM.Draw_Y = 0;
	gM.Draw_W = 800;
	gM.Draw_H = 600;
	
	gM.Main_mutex = SDL_CreateMutex();
	
	Head_Init (&gM.Head);
	
	Eye_Init (&gM.Left);
	Eye_Init (&gM.Right);
	
	gM.pGaze_Mutex = SDL_CreateMutex();
	
	dyn_config_init(&gM_DC);
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
	
	#if dM_Actions_Mode == 2
	/*mythread = */SDL_CreateThread(muhaha_chrdev_Thread, (void *)NULL);
	#endif
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
	SDL_mutexP (gM.Main_mutex);
	
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
				
				pin->U = 0x7;
				pin->V = 0x7;
			//	dpix(x,y)->U = 0x7;
			//	dpix(x,y)->V = 0x8;
			}
		}
		memcpy(gM.pDst, videoIn->framebuffer, videoIn->width * (videoIn->height) * 2);
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
		memcpy(gM.pDst, videoIn->framebuffer, videoIn->width * (videoIn->height) * 2);
	}
/*	printf ("World:\n");
	M4f_Print (&gM.World);
	printf ("Proj:\n");
	M4f_Print (&gM.Proj);/**/
	
	for (y = 0; y < 600; ++y)
		for (x = 800; x < 800+640; ++x)
			dspix(x, y) = 0;
	
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
		
	//	printf ("c:\t%f\t%f\t%f\n", gM.Head.M4.x03, gM.Head.M4.x13, gM.Head.M4.x23);
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
			
			
		}
		
		if (gM.GazeMode == 0) {
			Head_Dbg_Print (&gM.Head);
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
			
			V4f_DrawPosPos (&pbl, &pbr);
			V4f_DrawPosPos (&pbl, &pt);
			V4f_DrawPosPos (&pbr, &pt);
			
			V4f_DrawPosPos (&pbl, &pc);
			V4f_DrawPosPos (&pbr, &pc);
			V4f_DrawPosPos (&pt, &pc);
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
					V4f_DrawPosPos (&p00, &p01);	\
					V4f_DrawPosPos (&p01, &p02);	\
					V4f_DrawPosPos (&p02, &p03);	\
					V4f_DrawPosPos (&p03, &p00);	\
						\
					V4f_DrawPosPos (&p10, &p11);	\
					V4f_DrawPosPos (&p11, &p12);	\
					V4f_DrawPosPos (&p12, &p13);	\
					V4f_DrawPosPos (&p13, &p10);	\
						\
					V4f_DrawPosPos (&p00, &p10);	\
					V4f_DrawPosPos (&p01, &p11);	\
					V4f_DrawPosPos (&p02, &p12);	\
					V4f_DrawPosPos (&p03, &p13);	\
				}while(0)
			
			tV4f pos;
			pos.x = 0; pos.y = gM.Head.Mod.PC.y+0.15;	pos.z = gM.Head.Mod.PC.z-0.75;
			drect (16, 0.3f, 1.5f);
			
			pos.x = 0; pos.y = gM.Head.Mod.PC.y + 0.3f + 0.9f/2.0f;	pos.z = gM.Head.Mod.PC.z - 0.75f/2;
			drect (8, 0.9f, 0.75f);
			
			pos.x = 0; pos.y = gM.Head.Mod.PC.y + 1.2f + 0.9f/2.0f;	pos.z = gM.Head.Mod.PC.z - 0.75f/2;
			drect (1.6f, 0.9f, 0.75f);
			
			pos.x = gM.Head.Mod.PL.x; pos.y = gM.Head.Mod.PC.y+0.3+0.15;	pos.z = gM.Head.Mod.PC.z-1.5f;
			drect (1.5f, 0.3f, 4.8f);
			
			pos.x = gM.Head.Mod.PR.x;
			drect (1.5f, 0.3f, 4.8f);
			
			#undef drect
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
	
	
	Screen_Print (gM.aScreen + 0);
	Screen_Print (gM.aScreen + 1);
	Screen_Print (gM.aScreen + 2);
/*	Screen_Eye_Print (gM.aScreen + 1, &gM.Left);
	Screen_Eye_Print (gM.aScreen + 1, &gM.Right);
	Screen_Eye_Print (gM.aScreen + 2, &gM.Left);
	Screen_Eye_Print (gM.aScreen + 2, &gM.Right);/**/
	
	if (0) {
		tV2si pos;
		pos.x = gM.Head.DotL.P.x;
		pos.y = gM.Head.DotL.P.y;
		printf ("ehh pos	%f	%f\n", gM.Head.DotL.P.x, gM.Head.DotL.P.y);
		printf ("ehhh %f\n", NN_Sphere (&pos, 50, 55, 0, 0));
	}
	
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
/*	{
		printf ("crap\n");
		XEvent e;
		XNextEvent(gM.X.pDisp, &e);
		printf ("yeah\n");
	}/**/
	SDL_mutexV (gM.Main_mutex);
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
				//	SDL_mutexV (gM.pGaze_Mutex);
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


void	S_Draw_Line_2d	(si x0, si y0, si x1, si y1)
{
//	printf ("xy0 %ld\t%ld\t%ld\t%ld\n", x0, y0, x1, y1);
	si dx = x1 - x0;
	si dy = y1 - y0;
	
	ui dst_w = gM.pScreen->pitch/4;
	ui dst_h = gM.pScreen->h;
	u32 *pdst = gM.pScreen->pixels;
	
	u32 col = gColARGB;
	
	if (dx == 0 && dy == 0) {
		if ( (x0 >= 0) && (x0 < dst_w) && (y0 >= 0) && (y0 < dst_h) )
			*(pdst + x0 + (y0 * dst_w)) = col;
		return;
	}
	//#define dspix(_x,_y) (*((u32*)gM.pScreen->pixels + ((si)(_x) + (si)(_y)*gM.pScreen->pitch/4)))
	
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


