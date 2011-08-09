#include <sys/inotify.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <stdio.h>
#include <yaml.h>
#include <assert.h>

#include "muhaha.h"

void dyn_config_init(dyn_config *dc) {
	dc->entries = NULL;
	dc->count = 0;
}

void traverse(int level, dyn_config_entry *de) {
	while (de) {
		if (de->child) {
		//	printf("M%d:%s ", level, de->name);
			traverse(level+1, de->child);
		} else {
		//	printf("%d:%s ", level, de->name);
		}
		de = de->next;
	} 
}

int	strpathcmp	(const char *name, const char *path)
{
	const char *dot = strchr(path, '.');
	if (!dot)
		return strcmp (name, path);
	
	if (name[dot-path] != '\0')
		return 1;
	
	return strncmp (name, path, dot - path);
}

dyn_config_entry* dyn_find_entry(dyn_config_entry *de, const char *path) {
	char *pch = strchr(path, '.');
	while (de) {
		if (strpathcmp(de->name, path) == 0)
			if (pch && de->child)
				return dyn_find_entry(de->child , pch+1);
			else {
				return de;
			}
		de = de->next;
	}
	return 0; 
}

int dyn_get_value_float(dyn_config_entry *de, const char *path, float *val) {
	de = dyn_find_entry(de, path);
	if (!de)
		return 0;
	sscanf(de->value, "%f", val);
	return 1; 
}
int dyn_get_value_si(dyn_config_entry *de, const char *path, si *val) {
	de = dyn_find_entry(de, path);
	if (!de)
		return 0;
	if (de->value[0] == '0' && de->value[1] == 'x')
		sscanf(de->value, "%lx", val);
	else
		sscanf(de->value, "%ld", val);
	return 1; 
}
int dyn_get_value_u08(dyn_config_entry *de, const char *path, u08 *val) {
	de = dyn_find_entry(de, path);
	if (!de)
		return 0;
	sscanf(de->value, "%hhu", val);
	return 1; 
}

int dyn_get_value_enum(dyn_config_entry *de, const char *path, int *val) {
	de = dyn_find_entry(de, path);
	if (!de)
		return 0;
	#define dhack(name)	if (strcmp (de->value, #name) == 0) { *val = name; return 1; }
	
	dhack(eEye_Fit_Default)
	dhack(eEye_Fit_SFit)
	dhack(eEye_Fit_S0)
	dhack(eEye_Fit_S2Fit)
	dhack(eEye_Fit_S3Fit_START)
	dhack(eEye_Fit_S3Fit_Point)
	dhack(eEye_Fit_S3Fit_Eye)
	dhack(eEye_Fit_S3Fit_Eye3)
	dhack(eEye_Fit_S3Fit_END)
	dhack(eEye_Fit_S4Fit_START)
	dhack(eEye_Fit_S4Fit_Edge0)
	dhack(eEye_Fit_S4Fit_Edge1)
	dhack(eEye_Fit_S4Fit_Edge2)
	dhack(eEye_Fit_S4Fit_END)
	dhack(eEye_Fit_S5_START)
	dhack(eEye_Fit_S5_Diff)
	dhack(eEye_Fit_S5_END)
	dhack(eEye_Fit_C)
	
	dhack(eEye_PFit_S1)
	dhack(eEye_PFit_S2)
	dhack(eEye_PFit_S3)
	dhack(eEye_PFit_S3const)
	dhack(eEye_PFit_RANSAC)
	
	#undef dhack
	return 0; 
}

dyn_config_entry* dyn_parse(yaml_parser_t *parser) {
	yaml_event_t event;
	int done = 0, error = 0;
	int stack_size = 0;
	int prev_event_type = YAML_SCALAR_EVENT;
	
	dyn_config_entry *stack[16];
	dyn_config_entry *var_stack[16];
	char label[32];
	
	dyn_config_entry *cur_entry = malloc(sizeof(dyn_config_entry));
	
	while (!done)
	{
		if (!yaml_parser_parse(parser, &event)) {
			error = 1; break;
		}
		
		if (event.type == YAML_MAPPING_START_EVENT) {
			cur_entry->next = NULL;
			cur_entry->child = NULL;

			var_stack[stack_size-1] = cur_entry;
			var_stack[stack_size] = stack[stack_size] = cur_entry;
			stack_size++;
			
			prev_event_type = YAML_MAPPING_START_EVENT;
		//	printf("MAPPING\n");
		}
		
		if (event.type == YAML_SCALAR_EVENT) {
			if (prev_event_type == YAML_SCALAR_EVENT) {
				cur_entry->child = NULL; 
				cur_entry->next = NULL;
			//	sscanf(event.data.scalar.value, "%f", &cur_entry->value);
			//	printf ("got %s\n", event.data.scalar.value);
				strncpy (cur_entry->value, event.data.scalar.value, dyn_config_entry_value_max);
				if (cur_entry->value[dyn_config_entry_value_max-1] != '\0') {
					fprintf (stderr, "OVERFLOWa\n");
					exit (1);
				}
				var_stack[stack_size-1] = cur_entry;
			//	printf("SCALAR VAL: %s\n", event.data.scalar.value);
				prev_event_type = YAML_NO_EVENT;
			} else if (prev_event_type == YAML_MAPPING_START_EVENT) { 
			//	printf("AFTER MAP SCALAR NAME: %s\n", event.data.scalar.value);
				dyn_config_entry *new_entry = malloc(sizeof(dyn_config_entry));
				var_stack[stack_size-1]->child = new_entry;
				cur_entry = new_entry;
				memcpy(cur_entry->name, event.data.scalar.value, event.data.scalar.length);
				prev_event_type = YAML_SCALAR_EVENT;
			} else {
			//	printf("SCALAR NAME: %s\n", event.data.scalar.value);
				dyn_config_entry *new_entry = malloc(sizeof(dyn_config_entry));
				var_stack[stack_size-1]->next = new_entry;
				cur_entry = new_entry;
				memcpy(cur_entry->name, event.data.scalar.value, event.data.scalar.length);
				prev_event_type = YAML_SCALAR_EVENT;
			}
		}
		
		if (event.type == YAML_MAPPING_END_EVENT) {
			prev_event_type = YAML_MAPPING_END_EVENT;
			stack_size--;
		//	printf("MAPPING END\n");
		}

		done = (event.type == YAML_STREAM_END_EVENT);
		yaml_event_delete(&event);
	}
	
	fflush(stdout);
	
	if (!stack[0])	//every once in a while it would segfault on this, maybe because some editors
		return 0;	//first truncate we get an event then they write and we get another event?
	
	dyn_config_entry *ret_val = stack[0]->child;
	free(stack[0]);
	traverse(0, ret_val);
	return ret_val;
}


void dyn_config_read(dyn_config *dc, const char *f_name) {
	FILE *file;
	yaml_parser_t parser;

	file = fopen(f_name, "rb");
	assert(file);
	assert(yaml_parser_initialize(&parser));
	yaml_parser_set_input_file(&parser, file);
	
	dyn_config_entry *dce = dyn_parse(&parser);
	
	if (!dce)
		return;
	
	int val_int;
	u08 val_u08;
	si val_si;
	float val_f;
	//test
/*	if (dyn_get_value_float(dce, strdup("ala"), &val))
		printf("\n%f\n", val);
	if (dyn_get_value_float(dce, strdup("a.g"), &val))
		printf("\n%f\n", val);
	if (dyn_get_value_float(dce, strdup("z.g"), &val))
		printf("\n%f\n", val);
	if (dyn_get_value_float(dce, strdup("lukasz.e.a"), &val))
		printf("\n%f\n", val);/**/
	//end test
	#define drw_full_f(name,tgt)		if (dyn_get_value_float(dce, name , &val_f)) { tgt = val_f; } else { printf ("fail to get " name "\n"); }
	#define drw_full_u08(name,tgt)	if (dyn_get_value_u08(dce, name , &val_u08)) { tgt = val_u08; } else { printf ("fail to get " name "\n"); }
	
	#define drw_u08(name)	drw_full_f(#name, gM.name)
	#define drw_si(name)	if (dyn_get_value_si(dce, #name , &val_si)) { /*printf(#name " = %ld\n", val_si);/**/ gM.name = val_si; } else { printf ("fail to get " #name "\n"); }
	#define drw_f(name)	if (dyn_get_value_float(dce, #name , &val_f)) { /*printf(#name " = %f\n", val_f);/**/ gM.name = val_f; } else { printf ("fail to get " #name "\n"); }
	#define drw_e(name)	if (dyn_get_value_enum(dce, #name , &val_int)) { /*printf(#name " = %d\n", val_int);/**/ gM.name = val_int; }
	
	#define drw_full_v2f(name,tgt)		\
		do {		\
			drw_full_f(name ".x", tgt.x);		\
			drw_full_f(name ".y", tgt.y);		\
		}while(0)
	
	#define drw_full_v4f(name,tgt)		\
		do {		\
			drw_full_f(name ".x", tgt.x);		\
			drw_full_f(name ".y", tgt.y);		\
			drw_full_f(name ".z", tgt.z);		\
			tgt.w = 1.0f;		\
		}while(0)
	#define drw_v4f(name)	drw_full_v4f(#name, gM.name)
	
	#define drw_a_f(name,idx,rest)	\
		do {		\
			drw_full_f (#name "[" #idx "]." #rest, gM.name[idx].rest);	\
		}while(0)
	#define drw_a_v2f(name,idx,rest)	\
		do {		\
			drw_full_v2f (#name "[" #idx "]." #rest, gM.name[idx].rest);	\
		}while(0)
	
	#define drw_aa_f(name,i0,i1,rest)	\
		do {		\
			drw_full_f (#name "[" #i0 "][" #i1 "]." #rest, gM.name[i0][i1].rest);	\
		}while(0)
	#define drw_aa_v4f(name,i0,i1,rest)	\
		do {		\
			drw_full_v4f (#name "[" #i0 "][" #i1 "]." #rest, gM.name[i0][i1].rest);	\
		}while(0)
	
	drw_u08 (Eye_Line_Ray)
	drw_u08 (GazeAvg)
	
	drw_u08 (Pointer_Mode)
	
	
	drw_u08 (bHead_Eye_LineDraw)
	
	drw_u08 (bEye_S4_EdgeMark3_Micro)
	
	drw_f (Cam.Full_W)
	drw_f (Cam.Full_H)
	drw_f (Cam.Full_FOV)
	
	drw_f (Cam.Image_W)
	drw_f (Cam.Image_H)
	drw_f (Cam.Image_Zoom)
	drw_f (Cam.Image_FOV)
	drw_f (Cam.Image_FOV_W)
	drw_f (Cam.Image_FOV_H)
	
	drw_f (Cam.Focus)
	drw_f (Cam.Exposure)
	drw_f (Cam.Zoom)
	
	Cam_Param_Set (&gM.Cam);
	
	drw_f (View_W)
	drw_f (View_H)
	
	drw_f (Proj_W)
	drw_f (Proj_H)
	drw_f (Proj_N)
	drw_f (Proj_F)
	
	drw_f (Proj_L)
	drw_f (Proj_R)
	drw_f (Proj_B)
	drw_f (Proj_T)
	//void makeFrustum(double fovY, double aspectRatio, double front, double back)
/*	{
		double fovY = gM.Cam.Full_FOV, aspectRatio = (float)gM.Cam.Image_W/(float)gM.Cam.Image_H, front = 1, back = 10;
		const double DEG2RAD = M_PI / 180;
		
		double tangent = tan(fovY/2 * DEG2RAD);   // tangent of half fovY
		double height = front * tangent;          // half height of near plane
		double width = height * aspectRatio;      // half width of near plane
		
		// params: left, right, bottom, top, near, far
		gM.Proj_W = 2*width;
		gM.Proj_H = 2*height;
	//	gM.Proj_N = front;
	//	gM.Proj_F = back;
	}
/*	{
		double fovY = gM.Cam.Full_FOV, aspectRatio = (float)gM.Cam.Image_W/(float)gM.Cam.Image_H, front = 1, back = 10;
		const double DEG2RAD = M_PI / 180;
		
		double tangent = tan(fovY/2 * DEG2RAD);   // tangent of half fovY
		double width = front * tangent;      // half width of near plane
		double height = width / aspectRatio;          // half height of near plane
		
		// params: left, right, bottom, top, near, far
		gM.Proj_W = 2*width;
		gM.Proj_H = 2*height;
	//	gM.Proj_N = front;
	//	gM.Proj_F = back;
	}/**/
	
/*	gM.Proj.x00 = gM.Proj_N / (gM.Proj_W/2);
	gM.Proj.x11 = gM.Proj_N / (gM.Proj_H/2);
	
	gM.Proj.x22 = -(gM.Proj_F+gM.Proj_N) / (gM.Proj_F - gM.Proj_N);
	
	gM.Proj.x23 = (-2*gM.Proj_F*gM.Proj_N) / (gM.Proj_F - gM.Proj_N);
	gM.Proj.x32 = -1;/**/
	
	Proj_Cam ();
	
	Cam_Param_Set (&gM.Cam);
	
	drw_v4f (Head.Mod.PC);
	drw_v4f (Head.Mod.PL);
	drw_v4f (Head.Mod.PR);
	
	drw_v4f (Head.Mod.TInc);
	drw_v4f (Head.Mod.RInc);
	
	drw_f (Head.DotC.AngRes)
	drw_f (Head.DotC.Exp_R)
	drw_f (Head.DotC.Min_R)
	drw_f (Head.DotC.Max_R)
	drw_f (Head.DotC.PFit_R)
	drw_e (Head.DotC.Fit)
	drw_e (Head.DotC.PFit)
	drw_f (Head.DotC.Fit_Scale)
	drw_f (Head.DotC.Fit_Trans)
	drw_f (Head.DotC.S2Fit_Scale)
	drw_f (Head.DotC.S2Fit_Trans)
	
	drw_f (Head.DotC.S4_Pix_Bright)
	
	drw_si (Head.DotC.S5.Min_N)
	drw_f (Head.DotC.S5.Min_R)
	drw_f (Head.DotC.S5.Pix_Dark)
	drw_f (Head.DotC.S5.Pix_Bright)
	drw_f (Head.DotC.S5.Diff_Dist)
	drw_si (Head.DotC.S5.Pix_Diff_Start)
	drw_si (Head.DotC.S5.Pix_Diff_Min)
	
	
	gM.Head.DotC.LinView.x = 0;
	gM.Head.DotC.LinView.y = 0;
	
	
	memcpy (&gM.Head.DotL.Ax, &gM.Head.DotC.Ax, sizeof(gM.Head.DotC) - ((si)&gM.Head.DotC.Ax - (si)&gM.Head.DotC));
	memcpy (&gM.Head.DotR.Ax, &gM.Head.DotC.Ax, sizeof(gM.Head.DotC) - ((si)&gM.Head.DotC.Ax - (si)&gM.Head.DotC));
	
	drw_f (Head.DotL.Exp_R)
	drw_e (Head.DotL.Fit)
	drw_f (Head.DotL.Fit_Scale)
	drw_f (Head.DotL.Fit_Trans)
	drw_f (Head.DotL.S2Fit_Scale)
	drw_f (Head.DotL.S2Fit_Trans)
	
	drw_f (Head.DotR.Exp_R)
	drw_e (Head.DotR.Fit)
	drw_f (Head.DotR.Fit_Scale)
	drw_f (Head.DotR.Fit_Trans)
	drw_f (Head.DotR.S2Fit_Scale)
	drw_f (Head.DotR.S2Fit_Trans)
	
	drw_f (Head.DotC.LinView.x)
	drw_f (Head.DotC.LinView.y)
	drw_f (Head.DotL.LinView.x)
	drw_f (Head.DotL.LinView.y)
	drw_f (Head.DotR.LinView.x)
	drw_f (Head.DotR.LinView.y)
	
	
	drw_f (Left.AngRes)
	drw_f (Left.Exp_R)
	drw_f (Left.Min_R)
	drw_f (Left.Max_R)
	drw_f (Left.PFit_R)
	drw_e (Left.Fit)
	drw_e (Left.PFit)
	drw_f (Left.Fit_Scale)
	drw_f (Left.Fit_Trans)
	drw_f (Left.S2Fit_Scale)
	drw_f (Left.S2Fit_Trans)
	
	drw_f (Left.Pix_Dark)
	drw_f (Left.Pix_Bright)
	drw_f (Left.S4_Pix_Bright)
	
	drw_si (Left.S5.Min_N)
	drw_f (Left.S5.Min_R)
	drw_f (Left.S5.Pix_Dark)
	drw_f (Left.S5.Pix_Bright)
	drw_f (Left.S5.Diff_Dist)
	drw_si (Left.S5.Pix_Diff_Start)
	drw_si (Left.S5.Pix_Diff_Min)
	
	drw_f (Left.InHead.R)
	drw_v4f (Left.InHead.P);
	
	drw_f (Left.LinView.x)
	drw_f (Left.LinView.y)
	
//	memcpy (&gM.Right, &gM.Left, sizeof(gM.Left));
	
	drw_f (Right.AngRes)
	drw_f (Right.Exp_R)
	drw_f (Right.Min_R)
	drw_f (Right.Max_R)
	drw_f (Right.PFit_R)
	drw_e (Right.Fit)
	drw_e (Right.PFit)
	drw_f (Right.Fit_Scale)
	drw_f (Right.Fit_Trans)
	drw_f (Right.S2Fit_Scale)
	drw_f (Right.S2Fit_Trans)
	
	drw_f (Right.Pix_Dark)
	drw_f (Right.Pix_Bright)
	drw_f (Right.S4_Pix_Bright)
	
	drw_si (Right.S5.Min_N)
	drw_f (Right.S5.Min_R)
	drw_f (Right.S5.Pix_Dark)
	drw_f (Right.S5.Pix_Bright)
	drw_f (Right.S5.Diff_Dist)
	drw_si (Right.S5.Pix_Diff_Start)
	drw_si (Right.S5.Pix_Diff_Min)
	
	drw_f (Right.InHead.R)
	drw_v4f (Right.InHead.P);
	
	drw_f (Right.LinView.x)
	drw_f (Right.LinView.y)
	
	
	drw_f (tmp.x)
	drw_f (tmp.y)
	drw_f (tmp.z)
	
	drw_full_v4f ("Screen.C", gM.aScreen[0].C);
	drw_full_f ("Screen.W", gM.aScreen[0].W);
	drw_full_f ("Screen.H", gM.aScreen[0].H);
	gM.aScreen[1].C = gM.aScreen[0].C;
	gM.aScreen[2].C = gM.aScreen[0].C;
	
	
	drw_f (Screen_N)
	
	drw_a_f (aScreen,0,PixW);
	drw_a_f (aScreen,0,PixH);
	drw_a_v2f (aScreen,0,Off);
	
	drw_a_f (aScreen,1,PixW);
	drw_a_f (aScreen,1,PixH);
	drw_a_v2f (aScreen,1,Off);
	
	drw_a_f (aScreen,2,PixW);
	drw_a_f (aScreen,2,PixH);
	drw_a_v2f (aScreen,2,Off);
	
	Screen_Eye_Init (gM.aScreen + 0, &gM.Left);
	Screen_Eye_Init (gM.aScreen + 0, &gM.Right);
	
	Screen_Eye_Init (gM.aScreen + 1, &gM.Left);
	Screen_Eye_Init (gM.aScreen + 1, &gM.Right);
	
	Screen_Eye_Init (gM.aScreen + 2, &gM.Left);
	Screen_Eye_Init (gM.aScreen + 2, &gM.Right);
	
//	printf ("aScreen[0].Off %d %d\n", gM.aScreen[0].Off.x, gM.aScreen[0].Off.y);
//	Screen_Eye_PreCal (gM.aScreen + 0, &gM.Left);
//	Screen_Eye_PreCal (gM.aScreen + 0, &gM.Right);
	
/*	drw_aa_f (Left.InHead.aaCal,0,0,SX);
	drw_aa_f (Left.InHead.aaCal,0,0,SY);
	drw_aa_v4f (Left.InHead.aaCal,0,0,P);
	
	drw_aa_f (Left.InHead.aaCal,0,1,SX);
	drw_aa_f (Left.InHead.aaCal,0,1,SY);
	drw_aa_v4f (Left.InHead.aaCal,0,1,P);
	
	drw_aa_f (Left.InHead.aaCal,1,0,SX);
	drw_aa_f (Left.InHead.aaCal,1,0,SY);
	drw_aa_v4f (Left.InHead.aaCal,1,0,P);
	
	drw_aa_f (Left.InHead.aaCal,1,1,SX);
	drw_aa_f (Left.InHead.aaCal,1,1,SY);
	drw_aa_v4f (Left.InHead.aaCal,1,1,P);
	/**/
	
/*	drw_aa_f (Right.InHead.aaCal,0,0,SX);
	drw_aa_f (Right.InHead.aaCal,0,0,SY);
	drw_aa_v4f (Right.InHead.aaCal,0,0,P);
	
	drw_aa_f (Right.InHead.aaCal,0,1,SX);
	drw_aa_f (Right.InHead.aaCal,0,1,SY);
	drw_aa_v4f (Right.InHead.aaCal,0,1,P);
	
	drw_aa_f (Right.InHead.aaCal,1,0,SX);
	drw_aa_f (Right.InHead.aaCal,1,0,SY);
	drw_aa_v4f (Right.InHead.aaCal,1,0,P);
	
	drw_aa_f (Right.InHead.aaCal,1,1,SX);
	drw_aa_f (Right.InHead.aaCal,1,1,SY);
	drw_aa_v4f (Right.InHead.aaCal,1,1,P);/**/
	
	#define dact(name,press_stuff) drw_si(Action.name)
	#include "actions.h"
	#undef dact
	
	
	drw_f (Micro.SX)
	drw_f (Micro.SY)
	
	
	#undef drw_v4f
	#undef drw_f
	#undef drw_e
	//TODO free the dyn_config_entries...
	
	printf ("W %f H %f", gM.Proj_W, gM.Proj_H);
	printf ("   N %f F %f\n", gM.Proj_N, gM.Proj_F);
	
	yaml_parser_delete(&parser);
	
	assert(!fclose(file));
}

void dyn_config_watch(dyn_config *dc, const char *f_name) {
	dc->entries = NULL;
	dc->count = 0;
	
	int fd = inotify_init();
	int wd = inotify_add_watch(fd, "config.yaml", IN_MODIFY | IN_CLOSE_WRITE);
	struct inotify_event evt;
	
	while (1) {
		ssize_t size = read(fd, &evt, sizeof(evt));
		dyn_config_read(dc, f_name);
		
	}
	
	
}

/*
int main(int argc, char **argv) {
	dyn_config dc;
	dyn_config_read(&dc, "config.yaml");
	return 0;
}
*/
