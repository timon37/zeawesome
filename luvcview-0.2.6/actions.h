
dact (GazeHold,
	gM.bGazeHold = !gM.bGazeHold;
)
dact (Mode_Switch,
	switch (gM.Pointer_Mode) {
	case 0:
		++gM.Pointer_Mode;
		Marker_Show (&gM.CrossHair);
		break;
	case 1:
		++gM.Pointer_Mode;
	//	gM.Pointer_Mode = 0;
		Marker_Hide (&gM.CrossHair);
		break;
	case 2:
		gM.Pointer_Mode = 0;
		break;
	}
)

dact (EyeCent_CalNext,
	Head_Eye_LineAdd (&gM.Head, &gM.Left);
	Head_Eye_LineAdd (&gM.Head, &gM.Right);
)
dact (Back,
	if (gM.Left.InHead.Line_N)
		--gM.Left.InHead.Line_N;
	if (gM.Right.InHead.Line_N)
		--gM.Right.InHead.Line_N;
)

dact (Screen0,
	gM.Screen_CalIdx = 0;
	Screen_Cal_Init (gM.aScreen + gM.Screen_CalIdx);
)
dact (Screen1,
	gM.Screen_CalIdx = 1;
	Screen_Cal_Init (gM.aScreen + gM.Screen_CalIdx);
)
dact (Screen2,
	gM.Screen_CalIdx = 2;
	Screen_Cal_Init (gM.aScreen + gM.Screen_CalIdx);
)
dact (CalNext,
	Screen_Cal_Save (gM.aScreen + gM.Screen_CalIdx);
)
dact (Screen_Left,
	Screen_Cal_DoLeft (gM.aScreen + gM.Screen_CalIdx);
)
dact (Screen_Right,
	Screen_Cal_DoRight (gM.aScreen + gM.Screen_CalIdx);
)
dact (Screen_Down,
	Screen_Cal_DoDown (gM.aScreen + gM.Screen_CalIdx);
)
dact (Screen_Up,
	Screen_Cal_DoUp (gM.aScreen + gM.Screen_CalIdx);
)
#if 0
dact (Screen_TL,
	Screen_Cal_DoThis (gM.aScreen + gM.Screen_CalIdx, 0, 0);
)
dact (Screen_TC,
	Screen_Cal_DoThis (gM.aScreen + gM.Screen_CalIdx, 1, 0);
)
dact (Screen_TR,
	Screen_Cal_DoThis (gM.aScreen + gM.Screen_CalIdx, 2, 0);
)
dact (Screen_ML,
	Screen_Cal_DoThis (gM.aScreen + gM.Screen_CalIdx, 0, 1);
)
dact (Screen_MC,
	Screen_Cal_DoThis (gM.aScreen + gM.Screen_CalIdx, 1, 1);
)
dact (Screen_MR,
	Screen_Cal_DoThis (gM.aScreen + gM.Screen_CalIdx, 2, 1);
)
dact (Screen_BL,
	Screen_Cal_DoThis (gM.aScreen + gM.Screen_CalIdx, 0, 2);
)
dact (Screen_BC,
	Screen_Cal_DoThis (gM.aScreen + gM.Screen_CalIdx, 1, 2);
)
dact (Screen_BR,
	Screen_Cal_DoThis (gM.aScreen + gM.Screen_CalIdx, 2, 2);
)
#endif

dact (Micro_Tog,
	switch (gM.Micro.State) {
	case 0:
		gM.Micro.Gaze = gM.Gaze;
		gM.Micro.Head_CR = gM.Head.M4_R;
		
		gM.Micro.State = 1;
		break;
	case 1:
		gM.Micro.State = 0;
		break;
	}
)

dact (View_Zoom,
	static u08 zoomed = 0;
	if (zoomed) {
		zoomed = 0;
		gM.Draw_X = 0;
		gM.Draw_Y = 0;
		gM.Draw_W = 800;
		gM.Draw_H = 600;
	}else {
		zoomed = 1;
		gM.Draw_X = -gM.Cam.Image_W/2 + 400;
		gM.Draw_Y = -gM.Cam.Image_H/2 + 300;
		gM.Draw_W = gM.Cam.Image_W;
		gM.Draw_H = gM.Cam.Image_H;
	}
)



dact (Head_Eye_LineDraw,
	gM.bHead_Eye_LineDraw = !gM.bHead_Eye_LineDraw;
)
