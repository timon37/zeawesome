# zeawesome
Code for gaze &amp; head tracking

Sorry for the rought state, I never found the time to clean it up.

This is a bit complex, if you wish to play around with it just email
timon37@gmail.com and I'll gladly help you.

Here's an old demo video https://www.youtube.com/watch?v=r9VmyHDFNjY
This is pretty close to the best that I got out of my setup.
Later on I made quite a few code improvements (e.g. stereo cameras,
model-less tracking) but unfortunately because I moved a bunch of times
I never had a properly setup environment.
I use IR moded logitech c920 webcams.

This is a hodgepodge of various bits, many of them aren't really used
they were there just as a quick experiment (e.g. ransac).
I hijacked luvcview (got a half-functional guvcview version around as well),
as the image source and wrote my code mostly in muhaha.c/h and Eye.c


Copyright and license
---------------------

Copyright (C) 2016 Tomasz Borowik
Copyright (C) 2005-2008 Michel Xhaard


This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.


Credits
-------

Original code from luvcview: Copyright (C) 2005-2008 Michel Xhaard
Original code from spcaview: Copyright (C) 2003-2006 Michel Xhaard
AVI writing code from Avilib: Copyright (C) 1999 Rainer Johanni
SDL (Simple DirectMedia Layer)

Contributions by:
	Laurent Pinchart, Linux UVC driver (http://linux-uvc.berlios.de/)
	Martin Rubli, Logitech (http://www.quickcamteam.net/)
	... and many others whose names are in the ChangeLog.
