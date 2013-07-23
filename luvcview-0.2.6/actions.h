

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

dact (ReThr,
	for (int i = 0; i < gM.Cam_N; ++i) {
		tCam* pcam = &gM.aCam[i];
		
		tEye* peye = &gM.Left;
		Eye_CalcAYUV (peye, pcam, 5);	peye->aCam[i].FF.Y = ay + peye->aCam[i].FF.Y_Marg;
		
		peye = &gM.Right;
		Eye_CalcAYUV (peye, pcam, 5);	peye->aCam[i].FF.Y = ay + peye->aCam[i].FF.Y_Marg;
		
		peye = &gM.Head.DotC;
		Eye_CalcAYUV (peye, pcam, 4);	peye->aCam[i].FF.Y = ay + peye->aCam[i].FF.Y_Marg;
		
		peye = &gM.Head.DotL;
		Eye_CalcAYUV (peye, pcam, 4);	peye->aCam[i].FF.Y = ay + peye->aCam[i].FF.Y_Marg;
		for (si n = 0; n < gM.Head.Point_N; ++n) {
			peye = gM.Head.aPoint+n;
			Eye_CalcAYUV (peye, pcam, 4);	peye->aCam[i].FF.Y = ay + peye->aCam[i].FF.Y_Marg;
		}
		
		peye = &gM.Head.DotR;
		Eye_CalcAYUV (peye, pcam, 4);	peye->aCam[i].FF.Y = ay + peye->aCam[i].FF.Y_Marg;
	}
	tHead* p = &gM.Head;
	p->P.x = p->P.y = 0;
	p->P.z = -40;
	p->R_X = p->R_Y = p->R_Z = 0;
)


dact (Back,
	if (gM.Left.InHead.Line_N)
		--gM.Left.InHead.Line_N;
	if (gM.Right.InHead.Line_N)
		--gM.Right.InHead.Line_N;
)

dactD (Cam_TCal)
dactD (Cam_RCal)

dactD (Head_Snap)

dactD (EyeCent_CalNext)

dact (Head_Point_Train,
	dSafe_Main_S ();
	
	Head_Point_Train (&gM.Head);
	
	dSafe_Main_E ();
)

dact (Screen_Hide,
	gM.Screen_CalIdx = 0;
	Screen_MarkersHide (gM.aScreen + 0);
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
//	Screen_Cal_Save (gM.aScreen + gM.Screen_CalIdx);
	if (gM.Cal_bRecord) {
		gM.Cal_bRecord = 0;
		Screen_Cal_RecordEnd (gM.aScreen + gM.Screen_CalIdx);
	}else {
		Screen_Cal_RecordStart (gM.aScreen + gM.Screen_CalIdx);
		gM.Cal_bRecord = 1;
	}/**/
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
dact (Screen_Point,
	Screen_Cal_DoPoint (gM.aScreen + gM.Screen_CalIdx, gM.Gaze.x, gM.Gaze.y);
)
dact (Screen_PointDel,
	Screen_Cal_PointDel (gM.aScreen + gM.Screen_CalIdx);
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

dact2 (Mouse_Click_L,
//	xdo_active_mods_t *active_mods = NULL;
//	active_mods = xdo_get_active_modifiers(context->xdo);
//	xdo_clear_active_modifiers(context->xdo, window, active_mods);
//	mouse_press (1);
//	Window win;
//	xdo_window_get_active (gM.X.pDo, &win);
//	xdo_window_focus (gM.X.pDo, win);
//	printf ("						win %d\n", win);
//	xdo_keysequence_up(gM.X.pDo, CURRENTWINDOW, 1);
	
	
//	u32 state = xdo_get_input_state (gM.X.pDo);
//	printf ("state = %x\n", state);
///	xdo_keysequence_up (gM.X.pDo, CURRENTWINDOW, "Super_R", 0);
//	printf ("state = %x\n", state);
	
//	_xdo_send_key (gM.X.pDo, CURRENTWINDOW, charcodemap_t *key,
//                          int modstate, int is_press, useconds_t delay) {
	
///	xdo_mousedown(gM.X.pDo, CURRENTWINDOW, 1);
	
///	xdo_keysequence_down (gM.X.pDo, CURRENTWINDOW, "Super_R", 0);
/*	if (gM.X.Queue_N >= dM_X_Queue_NUM)
		break;
	si idx = (gM.X.Queue_C+gM.X.Queue_N) % dM_X_Queue_NUM;
	gM.X.aQueue[idx].Act = dM_X_Queue_M_Down;
	gM.X.aQueue[idx].Button = 1;
	++gM.X.Queue_N;/**/
,
//	mouse_release (1);

///	xdo_keysequence_up (gM.X.pDo, CURRENTWINDOW, "Super_R", 0);
///	xdo_mouseup(gM.X.pDo, CURRENTWINDOW, 1);
///	xdo_keysequence_down (gM.X.pDo, CURRENTWINDOW, "Super_R", 0);
/*	if (gM.X.Queue_N >= dM_X_Queue_NUM)
		break;
	si idx = (gM.X.Queue_C+gM.X.Queue_N) % dM_X_Queue_NUM;
	gM.X.aQueue[idx].Act = dM_X_Queue_M_Up;
	gM.X.aQueue[idx].Button = 1;
	++gM.X.Queue_N;/**/
)
dact2 (Mouse_Click_R,
	printf ("hahahahahah 2\n");
	mouse_press (3);
,
	printf ("hahahahahah 2 relea\n");
	mouse_release (3);
)
dact2 (Mouse_Click_M,
	mouse_press (2);
,
	mouse_release (2);
)
/**/

dact (Micro_Tog,
	switch (gM.Micro.State) {
	case 0:
		gM.Micro.Gaze = gM.Gaze;
		gM.Micro.Head_CR = gM.Head.M4_R;
		
		gM.Micro.R_X = gM.Head.R_X;
		gM.Micro.R_Y = gM.Head.R_Y;
		
		gM.Micro.State = 1;
		break;
	case 1:
		gM.Micro.State = 0;
		break;
	}
)

dact (View_Zoom,
/*	static u08 zoomed = 0;
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
	}*/
)



dact (Head_Eye_LineDraw,
	gM.bHead_Eye_LineDraw = !gM.bHead_Eye_LineDraw;
)


dact (Num1,
	gM.Left.GTF.x = 0;
	gM.Left.GTF.y = 0;
	
	gM.Right.GTF.x = 0;
	gM.Right.GTF.y = 0;
)

/*
dact (Num2,
	tEye* peye = &gM.Left;
	
	peye->tmp.P[0] = peye->P;
	Eye_GV_Get (peye, &peye->tmp.GV[0]);
	
	peye = &gM.Right;
	
	peye->tmp.P[0] = peye->P;
	Eye_GV_Get (peye, &peye->tmp.GV[0]);
)

dact (Num3,
	tEye* peye = &gM.Left;
	
	peye->tmp.P[1] = peye->P;
	Eye_GV_Get (peye, &peye->tmp.GV[1]);
	
	peye = &gM.Right;
	
	peye->tmp.P[1] = peye->P;
	Eye_GV_Get (peye, &peye->tmp.GV[1]);
)
*/

dact (Num4,
	tEye* peye = &gM.Left;
	
	peye->GTF = peye->tmp.GV[0];		V2f_sub_V2f (&peye->GTF, &peye->tmp.GV[1]);
	tV2f p = peye->tmp.P[0];		V2f_sub_V2f (&p, &peye->tmp.P[1]);
	V2f_div_V2f (&peye->GTF, &p);
	
	peye = &gM.Right;
	
	peye->GTF = peye->tmp.GV[0];	V2f_sub_V2f (&peye->GTF, &peye->tmp.GV[1]);
	p = peye->tmp.P[0];		V2f_sub_V2f (&p, &peye->tmp.P[1]);
	V2f_div_V2f (&peye->GTF, &p);
	
)

