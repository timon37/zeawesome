#include <stdio.h>
#include <ctype.h>
#include <math.h>
#include <assert.h>

#include <opencv/cv.h>
#include <opencv/highgui.h>

#include <glib.h>
#include <glib-object.h>

#include "gab.h"

void	TrainUtil_MouseHandler	(int event, int x, int y, int flags, void* _params)
{
	TrainUtil_ExampleParams *params = (TrainUtil_ExampleParams*) _params;
	
	if (flags & CV_EVENT_LBUTTONDBLCLK) {
		if (flags & CV_EVENT_FLAG_LBUTTON) {
			params->x1 = x;
			params->y1 = y;
		} else if ( flags & CV_EVENT_FLAG_RBUTTON) {
			params->x2 = x;
			params->y2 = y;
		}
	}
}

TrainUtil_Example*	TrainUtil_CreateExample	(IplImage *image, float decision)
{
	TrainUtil_Example* example = (TrainUtil_Example*) malloc(sizeof(TrainUtil_Example));
	
	example->length = image->width*image->height*image->nChannels;
	example->decision = decision;
	example->attrs = (float*) malloc(sizeof(float)*example->length);
	
	int idx = 0;
	
	for (int yi=0; yi<image->height; yi++) {
		for (int xi=0; xi<image->width; xi++) {
			int image_idx = yi*image->widthStep + xi;
			example->attrs[idx] = (unsigned char) image->imageData[image_idx];
			printf ("example->attrs[idx] %f\n", example->attrs[idx]);
			idx++;
			
		}
	}
	
	assert(example->length==idx);
	assert(example->length==64*64);
	
	return example;
}

//gray image as input
IplImage*	TrainUtil_SampleImage	(IplImage *image, TrainUtil_ExampleParams *params)
{
	int x_sample_size = 64, y_sample_size = 64;
	
	IplImage* sample = cvCreateImage(cvSize(x_sample_size, y_sample_size), IPL_DEPTH_8U, 1);
	
	float x_start = (float) params->x1;
	float x_end = (float) params->x2;
	float x_range = x_end - x_start;
	
	float y_start = (float) params->y1;
	float y_end = (float) params->y2;
	float y_range = y_end - y_start;
	
	
	for (int yi=0; yi<y_sample_size; yi++) {
		for (int xi=0; xi<x_sample_size; xi++) {
			float x_point = x_start + (xi * x_range / x_sample_size);
			
			int x_floor = floor(x_point);
			int x_ceil = x_floor + 1;
			float x_fraction = x_point - x_floor;
			
			float y_point = y_start + (yi * y_range / y_sample_size);
			
			int y_floor = floor(y_point);
			int y_ceil = y_floor + 1;
			float y_fraction = y_point - y_floor;
			
			int image_idx_tl = y_floor*image->widthStep + x_floor;
			int image_idx_tr = y_floor*image->widthStep + x_ceil;
			int image_idx_bl = y_ceil*image->widthStep + x_floor;
			int image_idx_br = y_ceil*image->widthStep + x_ceil;
			
			int sample_idx = yi*sample->widthStep + xi;
			
			/* Using openCV macros
			unsigned char* sample_pixel = &CV_IMAGE_ELEM(sample, unsigned char, yi, xi * sample->nChannels );
			*/
			
			float c_tl = (unsigned char) image->imageData[image_idx_tl];
			float c_tr = (unsigned char) image->imageData[image_idx_tr];
			float c_bl = (unsigned char) image->imageData[image_idx_bl];
			float c_br = (unsigned char) image->imageData[image_idx_br];
			
			float c_top = x_fraction * (c_tr - c_tl) + c_tl;
			float c_bottom = x_fraction * (c_br - c_bl) + c_bl;
			
			float c_final = y_fraction * (c_bottom - c_top) + c_top;
				
			sample->imageData[sample_idx] = (int) c_final;
		}
	}
	
	/*DEBUG
	cvNamedWindow("sample", 1);
	
	cvShowImage("sample", sample);
	if (cvWaitKey(0) == 'e')
		exit(1);
	cvDestroyWindow("sample");
	*/
	
	return sample;
}

int	TrainUtil_OpenList	(char *fileName)
{
	FILE* f = fopen(fileName, "rt");
	if (f) {
		int n_lines = 0;
		char buf[1000];
		
		TrainUtil_ExampleList *training_set = NULL;
		
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
		      
		      
		fclose(f);
		return 0;
	} else
		return 1;
}

int main	(int argc, char **argv)
{
	return TrainUtil_OpenList("list.txt");
}
