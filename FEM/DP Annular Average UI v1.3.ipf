#pragma rtGlobals=1		// Use modern global access method.
//#include "Annular Average2"

// Nice user interface for annular average of a diffraction pattern.  Creates
// a graph with some controls for finding the center, blocking the beam stop
// then computing the annular average.
//
// Start the routine with AnnularAverageGraph(<image wave>), or by selecting
// "Image Annular Average" from the "Macros" menu.
//
// v1.1  	added DP scale option wgs 8/29/05
// v1.2  	fixed bug in beam-stop block (Redimension/S not /R in AAGraphSetup)  pmv 9/8/05
//		added liberal wave name checking (and rejection)  pmv 9/8/05  
//		annular averaging excluding the corners is broken for some reason.  added a warning pmv 9/8/05
// v1.3	Added autocentering algorithm from Petrson T. C. Peterson et al. 
//		Ultramicoscopy 103, 275 (2005).  Added hand-centering "update" button

function AnnularAverageGraph(im)
	wave im
	
	AAGraphSetup(im)
	AAGraph(im)	
	
	showtools/a  //added WGS
	
end

function StartAnnularAverage()
	
	string wname
	Prompt wname,"Select an image wave", popup, WaveList("*", ";", "DIMS:2")
	DoPrompt "Select an image wave", wname

	if(cmpstr(wname, PossiblyQuoteName(wname)))
		DoAlert 0, "Liberal wave names (e.g. with spaces) are not supported.  Sorry."
		return 0
	endif
	
	if(!V_flag)
		wave im = $wname
		AnnularAverageGraph(im)
	endif
	
end

Menu "Macros"
	"Image Annular Average", StartAnnularAverage()
end

function AAGraphSetup(im)
	wave im

	string cfol = GetDataFolder(1)
	NewDataFolder/O/S root:Packages
	NewDataFolder/O/S dp

	// variables for auto-centering
	variable/g xcp_refined
	variable/g ycp_refined
	variable/g ac_size = 7, ac_radius = 300, ac_angle = 20
	variable/g xac_min, xac_max, yac_min, yac_max
	
	variable/g xcp
	variable/g ycp
	variable/g beam_blocked = 0
	variable/g dp_scale  //added 8/29/05 WGS
	string/g outname = "temp_aav"
	string/g imname = NameOfWave(im)
	string/g curfol = cfol

	Duplicate/O im dp_im
	Redimension/S dp_im	
	SetDataFolder $cfol
	
end

Function AnnularAverageButton(ctrlName) : ButtonControl
	String ctrlName

	SetDataFolder root:Packages:dp:
	wave im = dp_im

	variable xc, yc
	NVAR xcp_hand = xcp
	NVAR ycp_hand = ycp
	NVAR xcp_auto = xcp_refined
	NVAR ycp_auto = ycp_refined
	NVAR beam_blocked = beam_blocked
	
	NVAR dp_scale = dp_scale //added 8/29/05 WGS
	
	SVAR outname = outname
	SVAR imname = imname
	SVAR curfol = curfol
	UpdateCenters("")
	
	ControlInfo center_method
	if(V_value == 1)
		xc = xcp_hand
		yc = ycp_hand
	elseif(V_value == 2)
		xc = xcp_auto
		yc = ycp_auto
	else
		printf "Error in center method!  How did you get here?\r"
		return 0;
	endif
	
	ControlInfo cornerYN
	if(V_value)
		AnnularAverageCorners(im, xc, yc, 1)
	else
		DoAlert 1, "Excluding the corners is broken in v1.2.  Watch out!"
		AnnularAverage(im, xc, yc, 1)
	endif
	string n
	sprintf n, "Annular average of %s about (%g, %g)", imname, xc, yc
	if(beam_blocked)
		sprintf n, " excluding (%g, %g) to (%g, %g)", hcsr(A), vcsr(A), hcsr(B), hcsr(B)
	endif
	
	Note annular_av, n
	
	setscale/p x, 0,dp_scale, annular_av //added WGS 8/29/05
	
	duplicate/O annular_av $(curfol+outname)

	Killwaves annular_av
	SetDataFolder $curfol
	
	
End

Function BlockBeamStop(ctrlName) : ButtonControl
	String ctrlName

	NVAR beam_blocked = root:Packages:dp:beam_blocked
	wave im = root:Packages:dp:dp_im
	SVAR origname = root:Packages:dp:imname
	SVAR curfol = root:Packages:dp:curfol
	wave orig = $(curfol+origname)
	if(!waveexists(im))
		DoAlert 0, "Image wave not found."
	else
		im = orig
		im[pcsr(A),pcsr(B)][qcsr(A),qcsr(B)] = nan
	endif
	beam_blocked = 1
	SetDrawLayer/K ProgFront
	SetDrawEnv xcoord= top,ycoord= left, linethick=0
	DrawRect hcsr(A),vcsr(a),hcsr(B),vcsr(B)
	SetDrawLayer UserFront
End

function UpdateCenters(crtlName) : ButtonControl
	string crtlName

	NVAR xcp = root:Packages:dp:xcp
	NVAR ycp = root:Packages:dp:ycp

	nvar dp_scale = root:Packages:dp:dp_scale //added 8/29/05 WGS
	
	variable x0, x1, y0, y1
	string rec = WinRecreation("", 0)
	variable start = strsearch(rec, "DrawOval ", 0) + 9	
	variable stop = strsearch(rec, "\r", start)
	string pos = rec[start, stop]
	//print pos
	sscanf pos, "%g,%g,%g,%g", x0, y0, x1, y1
	xcp = round((x0 + x1)/2)
	ycp = round((y0 + y1)/2)
end

function UpdateAutocenterSearch()
	
	NVAR xac_min = root:Packages:dp:xac_min
	NVAR xac_max = root:Packages:dp:xac_max
	NVAR yac_min = root:Packages:dp:yac_min
	NVAR yac_max = root:Packages:dp:yac_max

	variable start, stop, x0, x1, y0, y1
	string rec = WinRecreation("", 0)
	start = strsearch(rec, "linefgc= (65280,16384,16384)", 0)
	start = strsearch(rec, "DrawRect", start) + 9
	stop = strsearch(rec, "\r", start)
	string pos = rec[start, stop]
	sscanf pos, "%g,%g,%g,%g", x0, y0, x1, y1
	xac_min = round(min(x0, x1))
	xac_max = round(max(x0, x1))
	yac_min = round(min(y0, y1))
	yac_max = round(max(y0 ,y1))

	//SetDrawEnv xcoord= top,ycoord= left,linethick= 2,linefgc= (65280,16384,16384),fillpat= 0
	//DrawRect 400,400,500,500	
	

end

Function AutocenterButton(ctrlName) : ButtonControl
	String ctrlName
	
	SetDataFolder root:Packages:dp:
	wave im = dp_im
	NVAR xcp_refined = xcp_refined
	NVAR ycp_refined = ycp_refined
	NVAR xmin = xac_min
	NVAR xmax = xac_max
	NVAR ymin = yac_min
	NVAR ymax = yac_max
	SVAR curfol = curfol
	NVAR size = ac_size
	NVAR radius = ac_radius
	NVAR angle = ac_angle

	string centers
	UpdateAutocenterSearch()
	printf "Searching area (%d, %d) to (%d, %d) for center.\r", xmin, ymin, xmax, ymax
	
	centers = AutoCenterFind(im, xmin, ymin, xmax, ymax, size, radius, angle)
	sscanf centers, "%d,%d", xcp_refined, ycp_refined

	variable i=1
	do
		sprintf centers, "centers_%d", i
		wave center_im = $centers
		if(!WaveExists(center_im))
			break
		endif
		Duplicate/O center_im $(curfol+centers)
		Killwaves center_im
		i+=1
	while(1)

	SetDataFolder $curfol
	i=1
	do
		sprintf centers, "centers_%d", i
		wave center_im = $centers
		if(!WaveExists(center_im))
			break
		endif
		AppendImage/T center_im
		i+=1
	while(1)
	
End

function AAGraph(im)
	wave im
	
	PauseUpdate; Silent 1		// building window...
	Display /W=(273,50.75,683.25,507.5) as "Compute Annular Average"
	AppendImage/T IM
	ModifyGraph margin(left)=14,margin(bottom)=126,margin(top)=14,margin(right)=14,width={Plan,1,top,left}
	ModifyGraph mirror=2
	ModifyGraph nticks=6
	ModifyGraph minor=1
	ModifyGraph fSize=8
	ModifyGraph standoff=0
	ModifyGraph tkLblRot(left)=90
	ModifyGraph btLen=3
	ModifyGraph tlOffset=-2
	Cursor/P/I/H=1/C=(65280,65280,0) A $NameOfWave(im), 300, 300; Cursor/P/I/H=1/C=(65280,65280,0) B $NameOfWave(im), 400, 400
	ShowTools
	Button ann_av,pos={19,583},size={90,20},proc=AnnularAverageButton,title="annular average"
	ValDisplay xcp_disp,pos={33,482},size={80,15},title="x center:"
	ValDisplay xcp_disp,labelBack=(65535,65535,65535),format="%d"
	ValDisplay xcp_disp,limits={0,0,0},barmisc={0,1000}
	ValDisplay xcp_disp,value= #"root:Packages:dp:xcp"
	ValDisplay ycp_disp,pos={126,482},size={80,15},title="y center"
	ValDisplay ycp_disp,labelBack=(65535,65535,65535),format="%d"
	ValDisplay ycp_disp,limits={0,0,0},barmisc={0,1000}
	ValDisplay ycp_disp,value= #"root:Packages:dp:ycp"
	Button blockstop,pos={19,556},size={90,20},proc=BlockBeamStop,title="block beam stop"
	Button blockstop,help={"Use the cursors to block a rectangular region containing the beam stop."}
	CheckBox cornerYN,pos={139,559},size={96,14},title="include corners?",value= 1
	SetVariable output_name,pos={153,585},size={300,16},title="output wavename"
	SetVariable output_name,labelBack=(65535,65535,65535)
	SetVariable output_name,value= root:Packages:dp:outname
	SetVariable dp_scale_online,pos={282,558},size={170,16},title="scaling (1/A / pixel)"
	SetVariable dp_scale_online,labelBack=(65535,65535,65535)
	SetVariable dp_scale_online,limits={0,inf,0},value= root:Packages:dp:dp_scale
	SetVariable ac_angle,pos={331,529},size={110,16},title="rotation angle"
	SetVariable ac_angle,labelBack=(65535,65535,65535)
	SetVariable ac_angle,limits={0,360,0},value= root:Packages:dp:ac_angle
	SetVariable ac_radius,pos={295,507},size={90,16},title="radius"
	SetVariable ac_radius,labelBack=(65535,65535,65535)
	SetVariable ac_radius,limits={0,inf,0},value= root:Packages:dp:ac_radius
	SetVariable ac_size,pos={397,507},size={60,16},title="size"
	SetVariable ac_size,labelBack=(65535,65535,65535)
	SetVariable ac_size,limits={1,inf,2},value= root:Packages:dp:ac_size
	ValDisplay xac,pos={295,482},size={80,15},title="x center:"
	ValDisplay xac,labelBack=(65535,65535,65535),format="%d"
	ValDisplay xac,limits={0,0,0},barmisc={0,1000}
	ValDisplay xac,value= #"root:Packages:dp:xcp_refined"
	ValDisplay yac,pos={397,482},size={80,15},title="y center:"
	ValDisplay yac,labelBack=(65535,65535,65535),format="%d"
	ValDisplay yac,limits={0,0,0},barmisc={0,1000}
	ValDisplay yac,value= #"root:Packages:dp:ycp_refined"
	Button autocenter,pos={397,455},size={75,20},proc=AutocenterButton,title="refine center"
	Button updatecenters,pos={126,455},size={85,20},proc=UpdateCenters,title="update center"
	PopupMenu center_method,pos={72,524},size={134,21},title="center method:"
	PopupMenu center_method,help={"Draw pattern center from by-hand or automatic method."}
	PopupMenu center_method,mode=2,popvalue="auto",value= #"\"hand;auto\""
	SetDrawLayer UserFront
	SetDrawEnv xcoord= top,ycoord= left,linethick= 2,linefgc= (65280,65280,0),fillpat= 0
	DrawOval 100,100,300,300
	SetDrawEnv fsize= 10
	DrawText 0.0269722297140329,1.07224465771627,"hand centering"
	SetDrawEnv fsize= 10
	DrawText 0.572785581854971,1.07460872391013,"auto centering"
	SetDrawEnv fillpat= 0
	DrawRect -0.00467289719626168,1.01240694789082,0.448598130841121,1.16873449131514
	SetDrawEnv fillpat= 0
	DrawRect 0.539007092198582,1.26315789473684,0.997635933806147,1.0125313283208
	SetDrawEnv xcoord= top,ycoord= left,linethick= 2,linefgc= (65280,16384,16384),fillpat= 0
	DrawRect 400,400,500,500
End
