##############################
# spcaview Makefile
##############################

INSTALLROOT=$(PWD)

CC=gcc
CPP=g++
INSTALL=install
APP_BINARY=luvcview
BIN=/usr/local/bin
MATH_LIB = -lm -lyaml
SDLLIBS = $(shell ../SDL2_bin/bin/sdl2-config --libs) 
SDLFLAGS = $(shell ../SDL2_bin/bin/sdl2-config --cflags)
V4L2LIBS = $(shell pkg-config --libs libv4l2)
V4L2FLAGS = $(shell pkg-config --cflags libv4l2)
#LIBX11FLAGS= -I/usr/X11R6/include -L/usr/X11R6/lib
VERSION = 0.2.6

#WARNINGS = -Werror
#WARNINGS = -Wall \
#           -Wundef -Wpointer-arith -Wbad-function-cast \
#           -Wcast-align -Wwrite-strings -Wstrict-prototypes \
#           -Wmissing-prototypes -Wmissing-declarations \
#           -Wnested-externs -Winline -Wcast-qual -W \
#           -Wno-unused
#           -Wunused

CFLAGS += -O2 -ggdb3 -std=gnu99 `pkg-config opencv --cflags --libs` `pkg-config fann --cflags --libs` -lX11 -lXi -lXtst  -DUSE_SDL  -DLINUX -DVERSION=\"$(VERSION)\" $(SDLFLAGS) $(V4L2FLAGS) $(WARNINGS) 
CPPFLAGS = $(CFLAGS)

OBJECTS = luvcview.o muhaha.o Eye.o muhaha_ransac_ellipse.o config.o gab.o Posit.o color.o utils.o v4l2uvc.o gui.o avilib.o


all:	luvcview

clean:
	@echo "Cleaning up directory."
	rm -f *.a *.o $(APP_BINARY) core *~ log errlog *.avi

$(OBJECTS):	muhaha.h
$(OBJECTS):	actions.h
$(OBJECTS):	priv_macros.h

# Applications:
luvcview:	$(OBJECTS)
	$(CC)	$(CFLAGS) $(OBJECTS) $(X11_LIB) $(XPM_LIB)\
		$(MATH_LIB) \
		$(SDLLIBS)\
		$(V4L2LIBS)\
		-o $(APP_BINARY)
#	$(CC)	-S $(CFLAGS) muhaha.c -o muhaha.S
	chmod 755 $(APP_BINARY)


install: luvcview
	$(INSTALL) -s -m 755 -g root -o root $(APP_BINARY) $(BIN) 
	rm -f $(BIN)/uvcview
