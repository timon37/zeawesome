
#ifndef inc_priv_macros_h
#define inc_priv_macros_h
//}

#define dSMALL_NUM  0.00000001

#define deg2rad	(M_PI/180.0f)
#define rad2deg	(180.0f/M_PI)
#define dpow2(_num) ((_num)*(_num))


#define dcam	aCam[pcam->Idx]

#define drgb(_r,_g,_b)	((tRGBA){.R = _r * 0xFF, .G = _g * 0xFF, .B = _b * 0xFF, .A = 0xFF})


#define dclip_l(val,min)	({typeof(val) ret; if ((val) < (min)) ret = (min); else ret = (val); ret;})
#define dclip_h(val,max)	({typeof(val) ret; if ((val) > (max)) ret = (max); else ret = (val); ret;})

#define dclip_lh(val,min,max)	({typeof(val) ret; if ((val) < (min)) ret = (min); else if ((val) > (max)) ret = (max); else ret = (val);; ret;})

#define ddot(x0,y0,x1,y1) ((x0)*(x1) + (y0)*(y1))

#define ddist2(x0,y0,x1,y1)	(((x1)-(x0))*((x1)-(x0)) + ((y1)-(y0))*((y1)-(y0)))

#define ddist(x0,y0,x1,y1) sqrt(ddist2(x0,y0,x1,y1))

#define dopix(_x,_y) ((tPix*)pcam->UVC->framebuffer + ((si)(_x) + (si)(_y)*pcam->Image_W))
#define dnpix(_x,_y) ((tPix*)((tRGBA*)pcam->SDL_Surf->pixels + ((si)(_x) + (si)(_y)*pcam->SDL_Surf->w)))
#define dpix(_x,_y) dnpix(_x,_y)


#define dsout(_x,_y) ((_x) < 0 || (_y) < 0 || (_x) >= pcam->SDL_Surf->w || (_y) >= pcam->SDL_Surf->h)
//#define dspix(_x,_y) (*((u32*)gM.pScreen->pixels + ((si)(_x) + (si)(_y)*gM.pScreen->pitch/4)))
//#define dspix(_x,_y) (*((u32*)pcam->SDL_Surf->pixels + ((si)(_x) + (si)(_y)*pcam->SDL_Surf->w)))
#define dspix(_x,_y) (*((u32*)pcam->SDL_Surf->pixels))

#define dmono2rgb(_g) ((u32)(_g) | (u32)(_g)<<8 | (u32)(_g)<<16)


#define dpixout(_x,_y) ((_x) < 0 || (_y) < 0 || (_x) >= pcam->Image_W || (_y) >= pcam->Image_H)

/*
#define dset_cyuv(_x,_y,y,u,v)	\
	do {	\
		if (dpixout(_x,_y))	\
			break;		\
		dnpix(_x,_y)->Y = y;	\
		dnpix(_x,_y)->U = u;	\
		dnpix(_x,_y)->V = v;	\
	}while(0)
*/

//#define dset_c0(_x,_y) dset_cyuv (_x,_y,0xFF,0x0,0x0)
//#define dset_c1(_x,_y) dset_cyuv (_x,_y,0xFF,0xF,0xF)

#define dset_crgb(_x,_y,r,g,b)	\
	do {	\
		if ((_x) < 0 || (_y) < 0 || (_x)/2 >= pcam->SDL_Surf->w || (_y)/2 >= pcam->SDL_Surf->h)	\
			break;		\
		tRGBA* pout = dnpix((_x)/2, (_y)/2);			\
		pout->R = r;			\
		pout->G = g;			\
		pout->B = b;			\
		pout->A = 0xFF;			\
	}while(0)

//#define dset_c0(_x,_y) ((*(tRGBA*)dnpix((_x)/2, (_y)/2)) = (tRGBA){ .R = 0xFF, .G = 0x00, .B = 0x00, .A = 0xFF })
//#define dset_c1(_x,_y) ((*(tRGBA*)dnpix((_x)/2, (_y)/2)) = (tRGBA){ .R = 0x00, .G = 0xFF, .B = 0x00, .A = 0xFF })
#define dset_c0(_x,_y) dset_crgb(_x,_y,0xFF,0x00,0x00)
#define dset_c1(_x,_y) dset_crgb(_x,_y,0x00,0xFF,0x00)



static inline f00	V2f_dot_V2f	(tV2f* pv0, tV2f* pv1)
{
	return pv0->x*pv1->x + pv0->y*pv1->y;
}
static inline f00	V4f_dot_V4f	(tV4f* pv0, tV4f* pv1)
{
	return pv0->x*pv1->x + pv0->y*pv1->y + pv0->z*pv1->z;
}

static inline void	V4f_set1	(tV4f* pv0, f00 s)
{
	pv0->x = s;
	pv0->y = s;
	pv0->z = s;
	pv0->w = 1;
}
static inline void	V4f_mul_S	(tV4f* pv0, f00 s)
{
	pv0->x *= s;
	pv0->y *= s;
	pv0->z *= s;
	pv0->w = 1;
}


static inline f00	V2f_cross	(tV2f* pv0, tV2f* pv1)
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

static inline f00	V2f_dist	(tV2f* pv)
{
	return sqrtf(pv->x*pv->x + pv->y*pv->y);
}

static inline f00	V2f_dist2_V2f	(tV2f* pv0, tV2f* pv1)
{
	return dpow2(pv1->x-pv0->x) + dpow2(pv1->y-pv0->y);
}
static inline f00	V2f_dist_V2f	(tV2f* pv0, tV2f* pv1)
{
	return sqrt(V2f_dist2_V2f(pv0, pv1));
}
static inline f00	V4f_dist2	(tV4f* pv)
{
	return pv->x*pv->x + pv->y*pv->y + pv->z*pv->z;
}
static inline f00	V4f_dist	(tV4f* pv)
{
	return sqrt(V4f_dist2(pv));
}


static inline f00	V4f_ry	(tV4f* pv)
{
	return  atan2 (pv->z, pv->x);
}



static inline void	V4f_norm	(tV4f* pv0)
{
	V4f_mul_S (pv0, 1.0f/V4f_dist (pv0));
}

static inline f00	V4f_dist_V4f	(tV4f* pv0, tV4f* pv1)
{
	return sqrt(dpow2(pv1->x-pv0->x) + dpow2(pv1->y-pv0->y) + dpow2(pv1->z-pv0->z));
}


static inline f00	angle_norm_0_2pi	(f00 a)
{
	while (a < 0) {
		a += 2*M_PI;
	}
	while (a > 2*M_PI) {
		a -= 2*M_PI;
	}
	return a;
}

static inline f00	angle_norm_pi_pi	(f00 a)
{
	while (a < -M_PI) {
		a += 2*M_PI;
	}
	while (a > M_PI) {
		a -= 2*M_PI;
	}
	return a;
}

static inline f00	angle_diff_norm_pi_pi	(f00 a0, f00 a1)
{
	return angle_norm_pi_pi(a0 - a1);
}




static inline void	Col_EyeSet	(tEye* peye)
{
	if (peye == &gM.Left)
		gCol = (tRGBA){.R = 0xFF, .G = 0x00, .B = 0x00, .A = 0x00};
	else
		gCol = (tRGBA){.R = 0x00, .G = 0xFF, .B = 0x00, .A = 0x00};
}

static inline void	Col_CamSet	(tCam* pcam)
{
	switch (pcam->Idx) {
	case 0:	gCol = (tRGBA){.R = 0xFF, .G = 0xFF, .B = 0x00, .A = 0x00};	break;
	case 1:	gCol = (tRGBA){.R = 0x00, .G = 0xFF, .B = 0xFF, .A = 0x00};	break;
	default:	gCol = (tRGBA){.R = 0xFF, .G = 0xFF, .B = 0xFF, .A = 0x00};	break;
	}
}


static inline void	Eye_Xset	(tEye* peye, f00 x)
{
	if (x < 20)
		x = 20;
	else if (x > pcam->Image_W-20)
		x = pcam->Image_W-20;
	peye->aCam[pcam->Idx].P.x = x;
}
static inline void	Eye_Yset	(tEye* peye, f00 y)
{
	if (y < 20)
		y = 20;
	else if (y > pcam->Image_H-20)
		y = pcam->Image_H-20;
	peye->aCam[pcam->Idx].P.y = y;
}
static inline void	Eye_XYset	(tEye* peye, f00 x, f00 y)
{
	Eye_Xset (peye, x);
	Eye_Yset (peye, y);
}


static inline void	EyeC_Xset	(tEye* peye, tCam* pcam, f00 x)
{
	if (x < 20)
		x = 20;
	else if (x > pcam->Image_W-20)
		x = pcam->Image_W-20;
	peye->aCam[pcam->Idx].P.x = x;
}
static inline void	EyeC_Yset	(tEye* peye, tCam* pcam, f00 y)
{
	if (y < 20)
		y = 20;
	else if (y > pcam->Image_H-20)
		y = pcam->Image_H-20;
	peye->aCam[pcam->Idx].P.y = y;
}
static inline void	EyeC_XYset	(tEye* peye, tCam* pcam, f00 x, f00 y)
{
	EyeC_Xset (peye, pcam, x);
	EyeC_Yset (peye, pcam, y);
}


static inline void	PrintMat	(CvMat *A)
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
	printf ("%8.3f ", (f00)cvGetReal2D(A, i, j));
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

static inline void	Iter_Optimize	(void* pclass, f00* apar, ui parn, f00 (*mse_cb)(void* pclass, f00* apar, ui parn))
{
	f00 omse = mse_cb(pclass, apar, parn);
	f00 err = 0;
	for (si i = 0; i < 500; ++i) {
		si fix = 0;
		for (si p = 0; p < parn; ++p) {
			f00 save = apar[p];
			f00 cor = apar[parn+p];
			apar[p] = save + cor;
			f00 nmse = mse_cb(pclass, apar, parn);
			if (nmse < omse) {
				//printf("fix %d %f\n", p, omse-nmse);
				omse = nmse;
				++fix;
				continue;
			}else {
				apar[p] = save - cor;
				nmse = mse_cb(pclass, apar, parn);
				//printf("omse %f	nmse1 %f\n", omse, nmse);
				if (nmse < omse) {
					//printf("fix %d %f\n", p, omse-nmse);
					omse = nmse;
					++fix;
					continue;
				}
				apar[p] = save;
			}
		}
		if (fix == 0) {
			for (si p = 0; p < parn; ++p) {
				apar[parn+p] *= 0.9;
			}
		}
	}
}
#endif
