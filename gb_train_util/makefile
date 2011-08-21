main: main.c gab.c gab.h
	gcc -std=c99 -g main.c gab.c -lm -lopencv_core -lopencv_highgui -lopencv_imgproc `pkg-config --cflags --libs glib-2.0` -o main
clean:
	rm -f *.o
	rm -f main
astyle:
	astyle --style=linux -t6 -E *.c *.h
