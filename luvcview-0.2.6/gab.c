
#include <stdio.h>
#include <math.h>

#include "gab.h"

void    TrainUtil_ExampleList_Add  (TrainUtil_ExampleList **list, TrainUtil_Example *example)
{
	if (!*list) {
		*list = (TrainUtil_ExampleList*) malloc(sizeof(TrainUtil_ExampleList));
		(*list)->example = example;
		(*list)->next = NULL;
	} else {
		TrainUtil_ExampleList *first = *list;
		*list = (TrainUtil_ExampleList*) malloc(sizeof(TrainUtil_ExampleList));
		(*list)->example = example;
		(*list)->next = first;
	}
}

int TrainUtil_ExampleList_Length    (TrainUtil_ExampleList *list)
{
	int length = 0;
	while (list) {
		length++;
		list = list->next;
	}
	return length;
}

//dunno how to make it nicer, weak spot of C :]
int compare_examples_idx = 0;
int compare_examples (const TrainUtil_Example **a, const TrainUtil_Example **b)
{
	float a_val = (*a)->attrs[compare_examples_idx];
	float b_val = (*b)->attrs[compare_examples_idx];
	
	if (a_val < b_val)
		return -1;
	else
		return 1;
}

void TrainUtil_FitDecisionStump(TrainUtil_Example **X, int X_size, int f_index, TrainUtil_Classifier *out_cls)
{
	float *s_yw = (float*) malloc(sizeof(float)*X_size);
	float *s_w = (float*) malloc(sizeof(float)*X_size);
	float *a = (float*) malloc(sizeof(float)*X_size);
	float *b = (float*) malloc(sizeof(float)*X_size);
	float *error = (float*) malloc(sizeof(float)*X_size);
	
	float e_yw = 0.0f, e_w = 0.0f;
	
	for (int i = 0; i < X_size; i++) {
		e_w += X[i]->w;
		e_yw += X[i]->w * X[i]->decision;
		s_w[i] = e_w;
		s_yw[i] = e_yw;
	}
	//calculate a, b, th
	for (int i = 0; i < X_size; i++) {
		b[i] = s_yw[i] / s_w[i];           //check this shit... ;)
		a[i] = (e_yw - s_yw[i]) / ((s_w[i] == 1.0f) ? 1.0f : 1.0f - s_w[i]) - b[i];
		//estimate error for each classifier
		error[i] = 1.0f - 2.0f * a[i] * ( e_yw - s_yw[i] ) - 2.0f * b[i] * e_yw
		           + ( a[i] * a[i] + 2.0f * a[i] * b[i] ) * (1.0f - s_w[i]) + b[i] * b[i];
	}
	
	int min_index = 0;
	float min_error = error[min_index];
	
	for (int i = 0; i < X_size; i++) {
		if (error[i] < min_error) {
			min_error = error[i];
			min_index = i;
		}
	}
	
	out_cls->index = f_index;
	out_cls->a = a[min_index];
	out_cls->b = b[min_index];
	out_cls->error = error[min_index];
	
	if (min_index == X_size - 1)
		out_cls->th = X[min_index]->attrs[f_index];
	else
		out_cls->th = ( X[min_index]->attrs[f_index]
		                + X[min_index + 1]->attrs[f_index]) / 2.0f;
	
	free(s_yw);
	free(s_w);
	free(a);
	free(b);
	free(error);
}

void    TrainUtil_GentleBoost    (
      //input:
      TrainUtil_Example **X,
      int X_size,
      TrainUtil_Example **T,
      int T_size,
      int n_rounds,
      //output:
      TrainUtil_Classifier **out_cls)
{
	int total_dim = X[0]->length;
	printf("TOTAL_DIM %d\n", total_dim);
	printf("X_size %d	T_size %d\n", X_size, T_size);
	
	float *stump_out = (float*) malloc(sizeof(float)*X_size);
	float *class_out = (float*) malloc(sizeof(float)*X_size);
	
	//allocate output classifer list
	*out_cls = (TrainUtil_Classifier*) malloc(sizeof(TrainUtil_Classifier)*n_rounds);
	
	//creates sorted arrays for each feature
	TrainUtil_Example ***sorted_X = (TrainUtil_Example***) malloc(sizeof(TrainUtil_Example**)*total_dim);	
	for (int i = 0; i < total_dim; i++) {
		TrainUtil_Example **X_to_sort = (TrainUtil_Example**) malloc(sizeof(TrainUtil_Example*)*X_size);
		memcpy(X_to_sort, X, sizeof(TrainUtil_Example*)*X_size);
		//set attribute on which sorting is performed
		compare_examples_idx = i;
		qsort(X_to_sort, X_size, sizeof(TrainUtil_Example*), compare_examples);
		assert(X_to_sort[0]->attrs[i] <= X_to_sort[1]->attrs[i]);
		sorted_X[i] = X_to_sort;
	}
	
	//initialise weights
	for (int i = 0; i < X_size; i++) {
		X[i]->w = 1.0f/X_size;
	}
	
	//main boosting loop
	for (int r = 0; r < n_rounds; r++) {
		//weight normalization
		float sum_w = 0.0f;
		for (int i = 0; i < X_size; i++)
			sum_w += X[i]->w;
		for (int i = 0; i < X_size; i++)
			X[i]->w /= sum_w;
			
		//find regression stump with best parameters
		TrainUtil_Classifier best;
		TrainUtil_FitDecisionStump(sorted_X[0], X_size, 0, &best);
		for (int f = 1; f < total_dim; f++) {
			TrainUtil_Classifier stump;
			TrainUtil_FitDecisionStump(sorted_X[f], X_size, f, &stump);
			if (stump.error < best.error)
				best = stump;
		}
		
		int misclassified = 0;
		for (int i = 0; i < X_size; i++) {
			float dec = (X[i]->attrs[best.index] > best.th) ? 1.0f : 0.0f;
			stump_out[i] = best.a * dec + best.b;
			class_out[i] += stump_out[i];
			X[i]->w *= exp( -X[i]->decision * stump_out[i] );
			if (class_out[i] * X[i]->decision <= 0.0f )
				misclassified ++;
		}
		
		//store weak learner
		(*out_cls)[r] = best;
		
		printf("ROUND: %d, MISCLASSIFIED: %d\%\n", r, 100*misclassified/X_size);
	}
	free(stump_out);
	free(class_out);
}


void	TrainUtil_Classify	(TrainUtil_Classifier *p, TrainUtil_Example* pex)
{
	float stump_out, class_out;
	
	float dec = (pex->attrs[p->index] > p->th) ? 1.0f : 0.0f;
	stump_out = p->a * dec + p->b;
	class_out += stump_out;
	
//	pex->w *= exp(-pex->decision * stump_out);
	
	return class_out;
}



