#include <sys/inotify.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <stdio.h>
#include <yaml.h>
#include <assert.h>

typedef struct _dyn_config_entry dyn_config_entry;
typedef struct _dyn_config_mapping dyn_config_mapping;
typedef struct _dyn_config dyn_config;
typedef enum DYN_ENTRY_TYPE dyn_entry_t;

struct _dyn_config_entry {
	dyn_config_entry *next;
	dyn_config_entry *child;
	char name[32];
	float value;
};

struct _dyn_config {
	dyn_config_entry *entries;
	int count;
};

void dyn_config_init(dyn_config *dc) {
	dc->entries = NULL;
	dc->count = 0;
}

 traverse(int level, dyn_config_entry *de) {
	while (de) {
		if (de->child) {
			printf("M%d:%s ", level, de->name);
			traverse(level+1, de->child);
		} else {
			printf("%d:%s ", level, de->name);
		}
		de = de->next;
	} 
}

int dyn_get_value(dyn_config_entry *de, char *path, float *val) {
	char *pch = strchr(path, '.');
	if (pch) {
		int pos = pch-path;
		path[pos] = 0;
	}
	while (de) {
		if (strcmp(de->name, path) == 0)
			if (pch && de->child)
				return dyn_get_value(de->child , pch+1, val);
			else {
				*val = de->value;
				return 1;
			}
		de = de->next;
	}
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
			printf("MAPPING\n");
		}
		
		if (event.type == YAML_SCALAR_EVENT) {
			if (prev_event_type == YAML_SCALAR_EVENT) {
				cur_entry->child = NULL; 
				cur_entry->next = NULL;
				sscanf(event.data.scalar.value, "%f", &cur_entry->value);
				var_stack[stack_size-1] = cur_entry;
				printf("SCALAR VAL: %s\n", event.data.scalar.value);
				prev_event_type = YAML_NO_EVENT;
			} else if (prev_event_type == YAML_MAPPING_START_EVENT) { 
				printf("AFTER MAP SCALAR NAME: %s\n", event.data.scalar.value);
				dyn_config_entry *new_entry = malloc(sizeof(dyn_config_entry));
				var_stack[stack_size-1]->child = new_entry;
				cur_entry = new_entry;
				memcpy(cur_entry->name, event.data.scalar.value, event.data.scalar.length);
				prev_event_type = YAML_SCALAR_EVENT;
			} else {
				printf("SCALAR NAME: %s\n", event.data.scalar.value);
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
			printf("MAPPING END\n");
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
	
	
	//test
	float val;
	if (dyn_get_value(dce, strdup("ala"), &val))
		printf("\n%f\n", val);
	if (dyn_get_value(dce, strdup("a.g"), &val))
		printf("\n%f\n", val);
	if (dyn_get_value(dce, strdup("z.g"), &val))
		printf("\n%f\n", val);
	if (dyn_get_value(dce, strdup("lukasz.e.a"), &val))
		printf("\n%f\n", val);
	//end test
	
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

int main(int argc, char **argv) {
	dyn_config dc;
	dyn_config_read(&dc, "config.yaml");
	return 0;
}
