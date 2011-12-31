
#include <opencv/cv.h>
#include <opencv/highgui.h>

typedef struct TrainUtil_ExampleParams TrainUtil_ExampleParams;

struct TrainUtil_ExampleParams {
	int x1;
	int x2;
	int y1;
	int y2;
};

typedef struct TrainUtil_Example TrainUtil_Example;

struct TrainUtil_Example {
	float *attrs;
	int length;
	float decision;
	float w;
};

typedef struct TrainUtil_ExampleList TrainUtil_ExampleList;

struct TrainUtil_ExampleList {
	TrainUtil_ExampleList *next;
	TrainUtil_Example *example;
};

typedef struct TrainUtil_Classifier TrainUtil_Classifier;

struct TrainUtil_Classifier {
	double a, b, th, error;
	int index;
};

void TrainUtil_ExampleList_Add(TrainUtil_ExampleList **list, TrainUtil_Example *example);
int TrainUtil_ExampleList_Length(TrainUtil_ExampleList *list);
void    TrainUtil_GentleBoost    (
      //input:
      TrainUtil_Example **X,
      int X_size,
      TrainUtil_Example **T,
      int T_size,
      int n_rounds,
      //output:
      TrainUtil_Classifier **out_cls);
