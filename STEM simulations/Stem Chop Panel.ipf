#pragma rtGlobals=1		// Use modern global access method.
#include "Stem Chop"

// fixed SetPhonon bug 02-07-10
// added configurations for exectuable files and email addresses 03-09-10 pmv
// 01-03-11 added ability to use autostem's ability to generate output as a function of thickness from a single calculation
// 01-09-11 added support for autopacbed program

Menu "Macros"
		"STEM Chop", stem_chop_panel()
		help = {"Generate input files for a STEM image simulation."}
end

function stem_chop_panel() : Panel

	wave test = $"root:Packages:stem_chop:sim_p"
	if(!waveexists(test))
		MakeControlWaves()
	endif
	
	PauseUpdate; Silent 1		// building window...
	NewPanel /W=(296.25,62.75,910.25,520.25)/K=1 as "STEM Chop"
	DoWindow/C StemChopWin
	//ShowTools
	TabControl t,pos={4,4},size={604,423},proc=SwitchTabs
	TabControl t,tabLabel(0)="Microscope"
	TabControl t,tabLabel(1)="Simulation",tabLabel(2)="Output Image"
	TabControl t,tabLabel(3)="Detectors",tabLabel(4)="Thickness"
	TabControl t,tabLabel(5)="Output Files",value= 0
	
	ValDisplay memory_est,pos={6,435},size={160,18},title="Memory estimate: "
	ValDisplay memory_est,fSize=12,limits={0,0,0},barmisc={0,1000}
	ValDisplay memory_est,value= #"root:Packages:stem_chop:mem_estimate"
	ValDisplay njobs,pos={212,436},size={90,14},title="# of jobs"
	ValDisplay njobs,limits={0,0,0},barmisc={0,1000}
	ValDisplay njobs,value= #"root:Packages:stem_chop:njobs"
	Button generate_inputs,pos={517,433},size={90,20},proc=PanelGenerateInputs,title="Generate Inputs"
	mic_param()
	sim_param()
	SimParamLabels()
	detect()
	thickness()
	imag_param()
	ino_files()
	SwitchTabs("sim", 0)
End


function GlobalLabels()

	SetDrawLayer UserBack
	DrawText 173,453,"MB"
end

function SwitchTabs(name, tab)
	string name
	variable tab
	
	string tabstr
	switch (tab)
	case 0:
		tabstr = "mic"
		MicParamLabels()
		break
	case 1:
		tabstr = "sim"
		SimParamLabels()
		break

	case 2:
		tabstr = "ima"
		ImagParamLabels()
		break

	case 3:
		tabstr = "det"
		DetectLabels()
		break

	case 4:
		tabstr = "thi"
		ThicknessLabels()
		break

	case 5:
		tabstr = "ino"
		InoFilesLabels()
		break
	
	default:
		printf "Uncoded tab value in SwitchTabs.  Oops.\r"
		return 0

	endswitch
	GlobalLabels()
	
	// Get a list of all the controls in the window
	variable i = 0
	string all = ControlNameList( "StemChopWin" )
	string thisControl
	
	do
		thisControl = StringFromList( i, all )
		if( strlen( thisControl ) <= 0 )
			break
		endif
		
		// Found another control.  Does it start with two letters and an underscore?
		if( !CmpStr( thisControl[3], "_" ) )
			// If it matches the current tab, show it.  Otherwise, hide it
			if( !CmpStr( thisControl[0,2], tabStr ) )
				ShowControl( thisControl, 0 )
			else
				ShowControl( thisControl, 1)
			endif
		endif
		i += 1
	while( 1 )
	
end


Function PanelGenerateInputs(ctrlName) : ButtonControl
	String ctrlName

	ReadDetectors()
	ReadThickness()
	SetPhonons()
	SetAberrations()
	wave sim_p = $"root:Packages:stem_chop:sim_p"
	NVAR mem = $"root:Packages:stem_chop:mem_estimate"
	sim_p[8] = mem
	wave imageout_p = $"root:Packages:stem_chop:imageout_p"
	ControlInfo thi_yn
	imageout_p[8] = V_value	
	
	ControlInfo sim_program
	
	switch(V_Value)
	
	case 1:
		StemChopPanel("autostem")
		break
	case 2:
		StemChopPanel("autopacbed")
		break
	default:
		printf "How did you get here?  Error in PanelGenerateInputs.\r"
		return 0
	endswitch
	
End

function StemChopPanel(program)
	string program
	
	variable outnum, target
	
	wave/T names = $"root:Packages:stem_chop:names"
	wave stem_p = $"root:Packages:stem_chop:stem_p"
	wave aber = $"root:Packages:stem_chop:aber"
	wave sim_p = $"root:Packages:stem_chop:sim_p"
	wave detect_p = $"root:Packages:stem_chop:detect_p"
	wave/T detect_name = $"root:Packages:stem_chop:detect_name"
	wave imageout_p = $"root:Packages:stem_chop:imageout_p"
	wave/t sim_paths = $"root:Packages:stem_chop:sim_paths"
	wave thick_p = $"root:Packages:stem_chop:thick_p"
	outnum = StemChopInputsandReassemble(program, names[3], names[1],  names[0], stem_p, aber, sim_p, detect_p, detect_name, imageout_p, thick_p)

	if(!outnum)
		printf "Error writing the input or reassembly file.  Exiting.\r"
		return 0
	endif
	target = GetTarget()
	StemChopCmd(target, sim_paths, names[3], names[1], names[2], sim_p[8], outnum)
	StemChopDescription(program, names[3], names, stem_p, aber, sim_p, detect_p, detect_name, imageout_p, thick_p, outnum)
	SaveControlWaves(names[3])

end

function sim_param() : Panel
	PauseUpdate; Silent 1		// building window...
	//NewPanel /W=(513.75,64.25,1125,488) as "Simulation Parameters"
	SetDrawLayer UserBack
	SetVariable sim_tran_nx,pos={22,61},size={110,19},title="X pixels: ",fSize=12
	SetVariable sim_tran_nx,limits={2,Inf,1},value= root:Packages:stem_chop:sim_p[0]
	SetVariable sim_tran_ny,pos={22,86},size={110,19},title="Y pixels",fSize=12
	SetVariable sim_tran_ny,limits={2,Inf,1},value= root:Packages:stem_chop:sim_p[1]
	SetVariable sim_probe_nx,pos={22,147},size={110,19},title="X pixels:",fSize=12
	SetVariable sim_probe_nx,limits={2,Inf,1},value= root:Packages:stem_chop:sim_p[2]
	SetVariable sim_probe_ny,pos={22,172},size={110,19},title="Y pixels:",fSize=12
	SetVariable sim_probe_ny,limits={2,Inf,1},value= root:Packages:stem_chop:sim_p[3]
	SetVariable sim_slicet,pos={22,205},size={150,19},title="slice thickness: "
	SetVariable sim_slicet,fSize=12
	SetVariable sim_slicet,limits={1,Inf,0.1},value= root:Packages:stem_chop:sim_p[4]
	SetVariable sim_temp,pos={232,86},size={130,19},title="temperature",fSize=12
	SetVariable sim_temp,limits={0,Inf,1},value= root:Packages:stem_chop:sim_p[5]
	SetVariable sim_phon_n,pos={232,111},size={180,19},title="phonon configurations:"
	SetVariable sim_phon_n,fSize=12
	SetVariable sim_phon_n,limits={1,Inf,1},value= root:Packages:stem_chop:sim_p[6]
	SetVariable sim_source,pos={232,148},size={125,19},title="source size",fSize=12
	SetVariable sim_source,limits={0,Inf,1},value= root:Packages:stem_chop:sim_p[7]
	CheckBox sim_phononYN,pos={232,61},size={105,16},title="Add phonons?",fSize=12
	CheckBox sim_phononYN,value= 0
	PopupMenu sim_chop_target,pos={232,203},size={227,24},title="chop target"
	PopupMenu sim_chop_target,fSize=12
	PopupMenu sim_chop_target,mode=2,popvalue="cluster",value= #"\"cluster;condor\""
	PopupMenu sim_program,pos={232,175},size={227,24},title="simulation program"
	PopupMenu sim_program,fSize=12
	PopupMenu sim_program,mode=2,popvalue="autostem",value= #"\"autostem;autopacbed\""
	SetVariable sim_cluster_exec,pos={22,259},size={400,16},title="executable:"
	SetVariable sim_cluster_exec,value= root:Packages:stem_chop:sim_paths[0]
	SetVariable sim_email,pos={22,283},size={200,16},title="email:"
	SetVariable sim_email,value= root:Packages:stem_chop:sim_paths[1]
	SetVariable sim_condor_exec,pos={22,347},size={400,16},title="exectuable:"
	SetVariable sim_condor_exec,value= root:Packages:stem_chop:sim_paths[2]
	SetVariable sim_condor_work,pos={22,369},size={400,16},title="working directory:"
	SetVariable sim_condor_work,value= root:Packages:stem_chop:sim_paths[3]

End

function SimParamLabels()

	SetDrawLayer/K UserBack
	SetDrawEnv linebgc= (52224,52224,52224),fillpat= 0,fillfgc= (54016,52480,50432)
	DrawRect 16,32,190,112
	DrawText 22,51,"Transmission Wavefunction:"
	DrawText 22,137,"Probe Wavefunction:"
	DrawText 232,51,"Phonons:"
	SetDrawEnv fillpat= 0
	DrawRect 15,196,147,116
	SetDrawEnv fillpat= 0
	DrawRect 227,32,430,135
	DrawText 368,102,"K"
	DrawText 363,165,"Å"
	DrawText 180,220,"Å"
	DrawRect 16,319,430,394
	DrawRect 16,236,430,307
	DrawText 22,256,"Cluster:"
	DrawText 22,342,"Condor:"
end	

Function imag_param() : Panel
	PauseUpdate; Silent 1		// building window...
	//NewPanel /W=(405.75,65,1046.25,488.75) as "Image Parameters"
	//ShowTools
	SetVariable ima_xi,pos={27,58},size={90,19},title="start x",fSize=12
	SetVariable ima_xi,limits={0,Inf,1},value= root:Packages:stem_chop:imageout_p[0]
	SetVariable ima_xf,pos={148,58},size={90,19},title="end x",fSize=12
	SetVariable ima_xf,limits={0,Inf,1},value= root:Packages:stem_chop:imageout_p[1]
	SetVariable ima_yi,pos={27,85},size={90,19},title="start y",fSize=12
	SetVariable ima_yi,limits={0,Inf,1},value= root:Packages:stem_chop:imageout_p[2]
	SetVariable ima_yf,pos={148,85},size={90,19},title="end y",fSize=12
	SetVariable ima_yf,value= root:Packages:stem_chop:imageout_p[3]
	SetVariable ima_nx,pos={301,58},size={100,19},title="x pixels",fSize=12
	SetVariable ima_nx,limits={2,Inf,1},value= root:Packages:stem_chop:imageout_p[4]
	SetVariable ima_ny,pos={301,83},size={100,19},title="y pixels",fSize=12
	SetVariable ima_ny,limits={2,Inf,1},value= root:Packages:stem_chop:imageout_p[5]
	SetVariable ima_nchunksx,pos={444,58},size={120,19},title="chunks in x",fSize=12
	SetVariable ima_nchunksx,limits={1,Inf,1},value= root:Packages:stem_chop:imageout_p[6]
	SetVariable ima_nchunksy,pos={444,83},size={120,19},title="chunks in y",fSize=12
	SetVariable ima_nchunksy,limits={2,Inf,1},value= root:Packages:stem_chop:imageout_p[7]
	ValDisplay ima_xspatial,pos={27,148},size={152,18},title="x spatial sampling"
	ValDisplay ima_xspatial,fSize=12,format="%0.3g",limits={0,0,0},barmisc={0,1000}
	ValDisplay ima_xspatial,value= #"root:Packages:stem_chop:xspatial"
	ValDisplay ima_yspatial,pos={27,173},size={152,18},title="y spatial sampling"
	ValDisplay ima_yspatial,fSize=12,format="%0.3g",limits={0,0,0},barmisc={0,1000}
	ValDisplay ima_yspatial,value= #"root:Packages:stem_chop:yspatial"
	ValDisplay ima_xmainpix,pos={301,148},size={70,18},title="x main",fSize=12
	ValDisplay ima_xmainpix,limits={0,0,0},barmisc={0,1000}
	ValDisplay ima_xmainpix,value= #"root:Packages:stem_chop:xmainpix"
	ValDisplay ima_ymainpix,pos={444,148},size={70,18},title="y main",fSize=12
	ValDisplay ima_ymainpix,limits={0,0,0},barmisc={0,1000}
	ValDisplay ima_ymainpix,value= #"root:Packages:stem_chop:ymainpix"
	ValDisplay ima_xleft,pos={301,173},size={80,18},title="x left over",fSize=12
	ValDisplay ima_xleft,limits={0,0,0},barmisc={0,1000}
	ValDisplay ima_xleft,value= #"root:Packages:stem_chop:xpixleft"
	ValDisplay ima_yleft,pos={444,173},size={80,18},title="y leftover",fSize=12
	ValDisplay ima_yleft,limits={0,0,0},barmisc={0,1000}
	ValDisplay ima_yleft,value= #"root:Packages:stem_chop:ypixleft"
End

function ImagParamLabels()

	SetDrawLayer/K UserBack
	DrawText 27,50,"Physical Image Size:"
	DrawText 301,50,"Pixel Size:"
	SetDrawEnv fillpat= 0
	DrawRect 19,30,258,110
	SetDrawEnv fillpat= 0
	DrawRect 294,30,573,109
	DrawText 123,75,"Å"
	DrawText 244,101,"Å"
	DrawText 244,76,"Å"
	DrawText 124,101,"Å"
	DrawText 27,143,"Spatial Sampling:"
	DrawText 301,143,"Pixels per chunk:"
	SetDrawEnv fillpat= 0
	DrawRect 19,121,258,201
	SetDrawEnv fillpat= 0
	DrawRect 294,121,573,200
	DrawText 187,190,"Å"
	DrawText 187,165,"Å"
 end

function ino_files() : Panel
	PauseUpdate; Silent 1		// building window...
	//NewPanel /W=(501.75,62.75,1140.75,487.25) as "Input and Output Files"
	//ShowTools
	Button ino_model,pos={285,96},size={100,20},proc=SelectModelFile,title="Select Model File"
	SetVariable ino_basename,pos={17,39},size={250,19},title="Output basename"
	SetVariable ino_basename,fSize=12,value= root:Packages:stem_chop:names[1]
	SetVariable ino_stem_cmt,pos={17,68},size={400,19},title="comment:",fSize=12
	SetVariable ino_stem_cmt,value= root:Packages:stem_chop:names[2]
	SetVariable ino_modeldisp,pos={17,96},size={250,19},title="model file",fSize=12
	SetVariable ino_modeldisp,value= root:Packages:stem_chop:names[0],noedit= 1
	SetVariable ino_pathdisp,pos={17,124},size={580,19},title="folder",fSize=12
	SetVariable ino_pathdisp,value= root:Packages:stem_chop:names[3],noedit= 1
End


function InoFilesLabels()

	SetDrawLayer/K Userback
	
end

Function TitanDefaultAberrations(ba) : ButtonControl
	STRUCT WMButtonAction &ba

	switch( ba.eventCode )
		case 2: // mouse up
			wave/t aber_text = $"root:Packages:stem_chop:aber_text"
			wave aber_default = $"root:Packages:stem_chop:aber_default"
			aber_text[][1,2] = num2str(aber_default[p][q-1])
			break
	endswitch

	return 0
End

function SetAberrations()

	wave aber = $"root:Packages:stem_chop:aber"
	wave/t aber_text = $"root:Packages:stem_chop:aber_text"
	
	aber = str2num(aber_text[p][q+1])
	Note/K aber, "units=0"
end


// which program we're simulating form.  1 = autostem, 2 = autopacbed
//function GetProgram()
//	
//	ControlInfo/W=sim_param sim_program
//	return V_value
//	
//end

// which execution platform we're chopping for.  1 = cluster, 2 = condor
function GetTarget()

	ControlInfo/W=StemChopWin sim_chop_target
	return V_value
	
end

Function SelectModelFile(ctrlName) : ButtonControl
	String ctrlName
	
	wave/T names = $"root:Packages:stem_chop:names"
	
	variable f, nitems
	string model_path, model_name
	Open/R/D/T=".xyz" f  
	
	if(strlen(S_filename))
		nitems = ItemsInList(S_filename, ":")
		model_name = StringFromList(nitems-1, S_filename, ":")
		model_path = S_filename[0, (strlen(S_filename) - strlen(model_name)-1)]
		names[0] = model_name
		names[4] = model_path
	endif

	//printf "S_filename = %s\r", S_filename
	//printf "model_path = %s\r", model_path
	//printf "model_name = %s\r", model_name

End

function detect() : Panel

	string fol = GetDataFolder(1)
	SetDataFolder root:Packages:stem_chop:
	Make/O/T/N=(numpnts(detect_name), 3) detect_lb
	Make/O/B/U/N=(numpnts(detect_name), 3) detect_sw
	wave/t detect_name = $"detect_name"
	wave detect_p = $"detect_p"
	detect_lb[][0] = detect_name[p]
	detect_lb[][1] = num2str(detect_p[p][0])
	detect_lb[][2] = num2str(detect_p[p][1])
	detect_sw = 2
	
	SetDimLabel 1, 0, 'detector name', detect_lb
	SetDimlabel 1, 1, 'inner angle (mr)', detect_lb
	SetDimLabel 1, 2, 'outer angle (mr)', detect_lb	

	// draw the panel
	PauseUpdate; Silent 1		// building window...
	//NewPanel /W=(605.25,68,1218.75,491.75) as "Detectors"
	ListBox det_lb,pos={15,53},size={447,362}
	ListBox det_lb,listWave=root:Packages:stem_chop:detect_lb
	ListBox det_lb,selWave=root:Packages:stem_chop:detect_sw,mode= 7
	Button det_adddet,pos={474,54},size={70,20},proc=AddDetector,title="Add detector"
	Button det_removedet,pos={474,78},size={95,20},proc=RemoveDetector,title="Remove detector"
	Button det_defaults,pos={474,114},size={80,20},proc=DefaultDetectors,title="Go to default:"
	PopupMenu det_select_default,pos={474,149},size={116,21}
	PopupMenu det_select_default,mode=3,popvalue=" Generic STEM",value= #"\"Titan STEM; Titan EFSTEM; Generic STEM\""

	SetDataFolder fol

End

function DetectLabels()

	SetDrawLayer/K UserBack
	DrawText 15,45,"Define detectors:"

end

Function AddDetector(ctrlName) : ButtonControl
	String ctrlName

	wave/T detect_lb = $"root:Packages:stem_chop:detect_lb"
	wave detect_sw = $"root:Packages:stem_chop:detect_sw"
	
	variable n = DimSize(detect_lb, 0)
	InsertPoints n, 1, detect_lb, detect_sw
	detect_lb[n][0] = "<name>"
	detect_lb[n][1] = "0"
	detect_lb[n][2] = "0"
	detect_sw[n][] = 2
end

Function RemoveDetector(ctrlName) : ButtonControl
	String ctrlName

	wave/T detect_lb = $"root:Packages:stem_chop:detect_lb"
	wave detect_sw = $"root:Packages:stem_chop:detect_sw"
	
	wavestats/q detect_sw
	DeletePoints V_maxRowLoc, 1, detect_lb, detect_sw

End

Function DefaultDetectors(ctrlName) : ButtonControl
	String ctrlName

	string fol = GetDataFolder(1)
	SetDataFolder root:Packages:stem_chop:

	ControlInfo det_select_default

	switch(V_Value)
	
	case 1:	// Titan STEM defaults
		wave/T name = $"titan_stem_detect_name"
		wave param = $"titan_stem_detect_p"
		break

	case 2:  // Titan EFSTEM
		wave/T name = $"titan_efstem_detect_name"
		wave param = $"titan_efstem_detect_p"
		break

	case 3: // standard STEM
		wave/T name = $"standard_stem_detect_name"
		wave param = $"standard_stem_detect_p"	
		break

	default:
		printf "How the hell did you get here?  Error in DefaultDetectors.\r"
		return 0
	endswitch

	Make/O/T/N=(numpnts(name), 3) detect_lb
	Make/O/B/U/N=(numpnts(name), 3), detect_sw
	SetDimLabel 1, 1, 'inner angle (mr)', detect_lb
	SetDimLabel 1, 2, 'outer angle (mr)', detect_lb
	
	SetDimLabel 1, 0, 'detector name', detect_lb
	detect_lb[][0] = name[p]
	detect_lb[][1] = num2str(param[p][0])
	detect_lb[][2] = num2str(param[p][1])
	detect_sw = 2

	SetDataFolder fol

End

function ReadDetectors()

	wave/T detect_lb = $"root:Packages:stem_chop:detect_lb"
	
	string fol = GetDataFolder(1)
	SetDataFolder root:Packages:stem_chop:
	
	variable ndet = DimSize(detect_lb, 0)
	Make/O/T/N=(ndet) detect_name
	Make/O/N=(ndet, 2) detect_p
	
	detect_name = detect_lb[p][0]
	detect_p[][0] = str2num(detect_lb[p][1])
	detect_p[][1] = str2num(detect_lb[p][2])
	
	SetDataFolder fol
end

function thickness() : Panel

	string fol = GetDataFolder(1)
	SetDataFolder root:Packages:stem_chop:
	Make/O/T/N=(numpnts(thick_p)) thick_lb
	Make/O/B/U/N=(numpnts(thick_p)) thick_sw
	wave thick_p = $"thick_p"
	thick_lb = num2str(thick_p)
	thick_sw = 2
	
	SetDimLabel 0, 0, 'intermediate thicnkess', thick_lb

	// draw the panel
	PauseUpdate; Silent 1		// building window...
	//NewPanel /W=(605.25,68,1218.75,491.75) as "Detectors"
	ListBox thi_lb,pos={15,53},size={338,362}
	ListBox thi_lb,listWave=root:Packages:stem_chop:thick_lb
	ListBox thi_lb,selWave=root:Packages:stem_chop:thick_sw,mode= 7
	CheckBox thi_yn,pos={366,59},size={140,14},title="intermediate thicknesses?"
	CheckBox thi_yn,value= 0	
	Button thi_addthi,pos={366,83},size={85,20},proc=AddThickness,title="Add thickness"
	Button thi_removedet,pos={366,107},size={105,20},proc=RemoveThickness,title="Remove thickness"
	ValDisplay thi_slice,pos={366,152},size={110,14},title="Slice thickness"
	ValDisplay thi_slice,limits={0,0,0},barmisc={0,1000},value=#"root:Packages:stem_chop:sim_p[4]"
	SetVariable thi_interval,pos={366,174},size={100,16},title="Slice interval"
	SetVariable thi_interval,limits={1,inf,1},value=root:Packages:stem_chop:thick_interval
	SetVariable thi_end,pos={366,195},size={120,16},title="End thickness"
	SetVariable thi_end,limits={1,inf,1},value=root:Packages:stem_chop:thick_end
	Button thi_gen,pos={366,218},size={120,20},proc=GenerateThicknesses,title="Generate thicknesses"

	SetDataFolder fol

end

function ThicknessLabels()

	SetDrawLayer/K UserBack
	DrawText 15,49,"Intermediate Thicknesses"
	DrawText 483,138,"A"
	DrawText 495,182,"A"
	DrawText 366,267,"Note: The exit-face thickness is"
	DrawText 366,284,"always reported."

end

function AddThickness(crtlName) : ButtonControl
	string crtlName
	
	wave/T thick_lb = $"root:Packages:stem_chop:thick_lb"
	wave thick_sw = $"root:Packages:stem_chop:thick_sw"
	
	variable n = numpnts(thick_lb)
	InsertPoints n, 1, thick_lb, thick_sw
	thick_lb[n] = "NaN"
	thick_sw[n] = 2
	
end

function RemoveThickness(crtlName) : ButtonControl
	string crtlName
	
	wave/T thick_lb = $"root:Packages:stem_chop:thick_lb"
	wave thick_sw = $"root:Packages:stem_chop:thick_sw"
	
	wavestats/q thick_sw
	DeletePoints V_maxRowLoc, 1, thick_lb, thick_sw
	
end

function ReadThickness()

	wave/T thick_lb = $"root:Packages:stem_chop:thick_lb"
	
	string fol = GetDataFolder(1)
	SetDataFolder root:Packages:stem_chop:
	
	variable nthick = DimSize(thick_lb, 0)
	Make/O/N=(nthick) thick_p
	thick_p = str2num(thick_lb)
	
	SetDataFolder fol
end


Function GenerateThicknesses(crtlName) : ButtonControl
	string crtlName

	wave sim_p = $"root:Packages:stem_chop:sim_p"
	variable slice_t = sim_p[4]
	NVAR slice_interval = $"root:Packages:stem_chop:thick_interval"
	NVAR thick_end = $"root:Packages:stem_chop:thick_end"
	wave/T thick_lb = $"root:Packages:stem_chop:thick_lb"
	wave thick_sw = $"root:Packages:stem_chop:thick_sw"
	
	variable nthick = floor(thick_end / (slice_t*slice_interval))
	redimension/n=(nthick) thick_lb, thick_sw
	thick_sw = 2
	thick_lb = num2str(slice_t*slice_interval*(p+1))	

	return 0
End


function SetPhonons()

	wave sim_p = $"root:Packages:stem_chop:sim_p"
	ControlInfo sim_phononYN
	if(!V_value)	// no phonons
		sim_p[5] = -1
	endif
end

// Show or hide any kind of control in the top window
Function ShowControl( name, disable )
	String name
	Variable disable
	
	// What kind of control is it?
	ControlInfo $name
	Variable type = v_flag
	switch( abs(type) )
		case 1:		// button
			Button $name disable=disable
			break
		
		case 2:		// checkbox
			CheckBox $name disable=disable
			break
		
		case 3:		// popup menu
			PopupMenu $name disable=disable
			break
		
		case 4:
			ValDisplay $name disable=disable
			break
		
		case 5:
			SetVariable $name disable=disable
			break
		
		case 6:
			Chart $name disable=disable
			break
		
		case 7:
			Slider $name disable=disable
			break
		
		case 8:
			TabControl $name disable=disable
			break
		
		case 9:
			GroupBox $name disable=disable
			break
		
		case 10:
			TitleBox $name disable=disable
			break
		
		case 11:
			ListBox $name disable=disable
			break
	endswitch
End

function mic_param() : Panel

	//PauseUpdate; Silent 1		// building window...
	//NewPanel /W=(483.75,49.25,1097.25,473) as "Microscope Parameters Working"
	SetVariable mic_kv,pos={44,54},size={110,19},title="voltage",fSize=12
	SetVariable mic_kv,limits={0,Inf,1},value= root:Packages:stem_chop:stem_p[0]
	SetVariable mic_cir,pos={44,81},size={180,19},title="aperture inner radius"
	SetVariable mic_cir,fSize=12
	SetVariable mic_cir,limits={0,Inf,1},value= root:Packages:stem_chop:stem_p[1]
	SetVariable mic_cor,pos={44,108},size={190,19},title="aperture outer radius"
	SetVariable mic_cor,fSize=12
	SetVariable mic_cor,limits={0,Inf,1},value= root:Packages:stem_chop:stem_p[2]
	ListBox mic_aber,pos={309,54},size={255,213},frame=4
	ListBox mic_aber,listWave=root:Packages:stem_chop:aber_text
	ListBox mic_aber,selWave=root:Packages:stem_chop:aber_selwave,mode= 5
	Button mic_default_aber,pos={396,276},size={80,20},proc=TitanDefaultAberrations,title="Titan defaults"

End

function MicParamLabels()

	SetDrawLayer/K UserBack
	DrawText 229, 98,"mrad"
	DrawText 242, 125,"mrad"
	DrawText 162,69,"kV"
end

