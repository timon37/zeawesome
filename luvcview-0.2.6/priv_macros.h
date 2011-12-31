
#ifndef inc_priv_macros_h
#define inc_priv_macros_h
//}

#define dSMALL_NUM  0.00000001

#define deg2rad	(M_PI/180.0f)
#define rad2deg	(180.0f/M_PI)
#define dpow2(_num) ((_num)*(_num))

#define dclip_l(val,min)	({typeof(val) ret; if ((val) < (min)) ret = (min); else ret = (val); ret;})
#define dclip_h(val,max)	({typeof(val) ret; if ((val) > (max)) ret = (max); else ret = (val); ret;})

#define dclip_lh(val,min,max)	({typeof(val) ret; if ((val) < (min)) ret = (min); else if ((val) > (max)) ret = (max); else ret = (val);; ret;})

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

#define dset_c0(_x,_y) dset_cyuv (_x,_y,0xFF,0x0,0x0)
#define dset_c1(_x,_y) dset_cyuv (_x,_y,0xFF,0xF,0xF)



static inline float	V2f_dot_V2f	(tV2f* pv0, tV2f* pv1)
{
	return pv0->x*pv1->x + pv0->y*pv1->y;
}
static inline float	V4f_dot_V4f	(tV4f* pv0, tV4f* pv1)
{
	return pv0->x*pv1->x + pv0->y*pv1->y + pv0->z*pv1->z;
}

static inline void	V4f_mul_S	(tV4f* pv0, float s)
{
	pv0->x *= s;
	pv0->y *= s;
	pv0->z *= s;
}


static inline float	V2f_cross	(tV2f* pv0, tV2f* pv1)
{
	return pv0->x*pv1->y - pv0->y*pv1->x;
}

static inline void	V4f_cross	(tV4f* pret, tV4f* pv0, tV4f* pv1)
{
	tV4f ret;
	ret.x = pv0->y*pv1->z - pv0->z*pv1->y;
	ret.y = pv0->z*pv1->x - pv0->x*pv1->z;
	ret.z = pv0->x*pv1->y - pv0->y*pv1->x;
	ret.w = 1;
	*pret = ret;
}

static inline float	V2f_dist	(tV2f* pv)
{
	return sqrtf(pv->x*pv->x + pv->y*pv->y);
}

static inline float	V4f_dist2	(tV4f* pv)
{
	return pv->x*pv->x + pv->y*pv->y + pv->z*pv->z;
}
static inline float	V4f_dist	(tV4f* pv)
{
	return sqrtf(V4f_dist2(pv));
}


static inline float	V4f_ry	(tV4f* pv)
{
	return  atan2 (pv->z, pv->x);
}



static inline void	V4f_norm	(tV4f* pv0)
{
	V4f_mul_S (pv0, 1.0f/V4f_dist (pv0));
}

static inline float	V4f_dist_V4f	(tV4f* pv0, tV4f* pv1)
{
	return sqrtf(dpow2(pv1->x-pv0->x) + dpow2(pv1->y-pv0->y) + dpow2(pv1->z-pv0->z));
}


static inline float	angle_norm_0_2pi	(float a)
{
	while (a < 0) {
		a += 2*M_PI;
	}
	while (a > 2*M_PI) {
		a -= 2*M_PI;
	}
	return a;
}

static inline float	angle_norm_pi_pi	(float a)
{
	while (a < -M_PI) {
		a += 2*M_PI;
	}
	while (a > M_PI) {
		a -= 2*M_PI;
	}
	return a;
}

static inline float	angle_diff_norm_pi_pi	(float a0, float a1)
{
	return angle_norm_pi_pi(a0 - a1);
}




static inline void	Col_EyeSet	(tEye* peye)
{
	if (peye == &gM.Left)
		gColARGB = 0x00FF00;
	else
		gColARGB = 0xFF0000;
}


static inline void	Eye_Xset	(tEye* peye, float x)
{
	if (x < 20)
		x = 20;
	else if (x > videoIn->width-20)
		x = videoIn->width-20;
	peye->P.x = x;
}
static inline void	Eye_Yset	(tEye* peye, float y)
{
	if (y < 20)
		y = 20;
	else if (y > videoIn->height-20)
		y = videoIn->height-20;
	peye->P.y = y;
}
static inline void	Eye_XYset	(tEye* peye, float x, float y)
{
	Eye_Xset (peye, x);
	Eye_Yset (peye, y);
}


#endif
