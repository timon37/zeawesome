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
	
	int val_int;
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
	#define drw_f(name)	if (dyn_get_value_float(dce, #name , &val_f)) { /*printf(#name " = %f\n", val_f);/**/ gM.name = val_f; }
	#define drw_e(name)	if (dyn_get_value_enum(dce, #name , &val_int)) { /*printf(#name " = %d\n", val_int);/**/ gM.name = val_int; }
//	.Y_Level = 128,
//	.DeSat = 0,
	/*
	drw_f (Proj.x00)
	drw_f (Proj.x01)
	drw_f (Proj.x02)
	drw_f (Proj.x03)
	drw_f (Proj.x10)
	drw_f (Proj.x11)
	drw_f (Proj.x12)
	drw_f (Proj.x13)
	drw_f (Proj.x20)
	drw_f (Proj.x21)
	drw_f (Proj.x22)
	drw_f (Proj.x23)
	drw_f (Proj.x30)
	drw_f (Proj.x31)
	drw_f (Proj.x32)
	drw_f (Proj.x33)*/
	
	drw_f (Cam.Full_W)
	drw_f (Cam.Full_H)
	drw_f (Cam.Full_FOV)
	
	drw_f (Cam.Image_W)
	drw_f (Cam.Image_H)
	drw_f (Cam.Image_FOV)
	
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
	{
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
	
	gM.Proj.x00 = gM.Proj_N / (gM.Proj_W/2);
	gM.Proj.x11 = gM.Proj_N / (gM.Proj_H/2);
	
	gM.Proj.x22 = -(gM.Proj_F+gM.Proj_N) / (gM.Proj_F - gM.Proj_N);
	
	gM.Proj.x23 = (-2*gM.Proj_F*gM.Proj_N) / (gM.Proj_F - gM.Proj_N);
	gM.Proj.x32 = -1;
	
	drw_f (Head.DotC.Exp_R)
	drw_e (Head.DotC.Fit)
	drw_f (Head.DotC.Fit_Scale)
	drw_f (Head.DotC.Fit_Trans)
	drw_f (Head.DotC.S2Fit_Scale)
	drw_f (Head.DotC.S2Fit_Trans)
	
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
	
	drw_f (Left.Exp_R)
	drw_e (Left.Fit)
	drw_f (Left.Fit_Scale)
	drw_f (Left.Fit_Trans)
	drw_f (Left.S2Fit_Scale)
	drw_f (Left.S2Fit_Trans)
	
	drw_f (Left.LinView.x)
	drw_f (Left.LinView.y)
	
	drw_f (Right.Exp_R)
	drw_e (Right.Fit)
	drw_f (Right.Fit_Scale)
	drw_f (Right.Fit_Trans)
	drw_f (Right.S2Fit_Scale)
	drw_f (Right.S2Fit_Trans)
	
	drw_f (Right.LinView.x)
	drw_f (Right.LinView.y)
	
	drw_f (tmp.x)
	drw_f (tmp.y)
	drw_f (tmp.z)
	
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
	int wd = inotify_add_watch(fd, "config.yaml", IN_MODIFY);
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
