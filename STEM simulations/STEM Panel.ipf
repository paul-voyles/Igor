#pragma rtGlobals=1		// Use modern global access method.

// first version completed 1/21/10
// bug in CTF1D calling the wrong button function and aberrations not updating properly fixed 1/25/10
// added source size support, 2/12/10
// modified for compatibility with Stem Chop Panel.ipf on 12/14/10 pmv

Menu "Macros"
	"Convolution STEM", stem_panel()
	help = {"Convolution STEM PSF, CTF, and image calculations"}
end


function STEMPanelSetup()
	
	string start_fol = GetDataFolder(1)
	
	NewDataFolder/O/S root:Packages
	NewDataFolder/O/S STEM
	
	string/G start_folder = start_fol
	string/G outname = "temp"
	variable/G keV, Cc, dE, ap, ds
	keV = 200
	Cc = 1.4
	dE = 0.73
	ap = 17.5
	ds = 0.5

	Make/O/N=(12, 2) aber, aber_input, aber_default
	Make/O/T/N=(12, 4) aber_text
	aber_default[][0] = {0, 0, 22.56, 22.08, 0.1198, 0.9018, 0.04964, 28.43, 11.84, 8.456, 0.622, 2.811}
	aber_default[][1] = {0, 0, -20.1, -7.5, 0, -170.1, 20.9, -120.6, 153.8, 76.1, 0, -125.5} 
	aber_input = aber_default
	aber = aber_default
	aber_text[][0] = {"C1", "A1", "A2", "B2", "C3", "A3", "S3", "A4", "D4", "B4", "C5", "A5"}
	aber_text[][1] = num2str(aber_default[p][0])
	aber_text[][2] = num2str(aber_default[p][1])
	aber_text[][3] = {"nm", "nm", "nm", "nm", "um", "um", "um", "um", "um", "um", "mm", "mm"}
	SetDimLabel 1, 0, Aberration, aber_text
	SetDimLabel 1, 1, Coefficient, aber_text
	SetDimLabel 1, 2, Angle, aber_text
	SetDimLabel 1, 3, Units, aber_text
	Make/U/B/O/N=(12,4) aber_selwave
	aber_selwave[][0] = 0
	aber_selwave[][1] = 2
	aber_selwave[][2] = 2
	aber_selwave[][3] = 0
	
	variable/G nk
	nk = 1024
	
	variable/G df_min, df_max, df_steps
	df_min = 0
	df_max = 200
	df_steps = 50
	
	SetDataFolder $start_fol
end

Function TitanDefaultAberrations2(ba) : ButtonControl
	STRUCT WMButtonAction &ba

	switch( ba.eventCode )
		case 2: // mouse up
			wave/t aber_text = $"root:Packages:STEM:aber_text"
			wave aber_default = $"root:Packages:STEM:aber_default"
			aber_text[][1,2] = num2str(aber_default[p][q-1])
			break
	endswitch

	return 0
End

function SetAberrations2()

	wave aber = $"root:Packages:STEM:aber"
	wave/t aber_text = $"root:Packages:STEM:aber_text"
	
	aber = str2num(aber_text[p][q+1])
	Note/K aber, "units=0"
end

function/S PrintAberrations(aber)
	wave aber
	
	string p
	sprintf p, "C1 = %g nm; A1 = %g nm, %g deg; A2 = %g nm, %g deg", aber[0][0], aber[1][0], aber[1][1], aber[2][0], aber[2][1]
	sprintf p, "%s; B2 = %g nm, %g deg; C3 = %g um; A3 = %g um, %g deg", p, aber[3][0], aber[3][1], aber[4][0], aber[5][0], aber[5][1]
	sprintf p, "%s; S3 = %g, um, %g deg; A4 = %g um; %g deg; D4 = %g um, %g deg", p, aber[6][0], aber[6][1], aber[7][0], aber[7][1], aber[8][0], aber[8][1]
	sprintf p, "%s; B4 = %g um, %g deg; C5 = %g mm; A5 = %g mm, %g deg", p, aber[9][0], aber[9][1], aber[10][0], aber[11][0], aber[11][1]
	
	return p
	
end

Function STEMPanelGo(ctrlName) : ButtonControl
	String ctrlName

	string fol = GetDataFolder(1)
	
	SetDataFolder $"root:Packages:STEM:"

	SVAR outn = $"outname"

	NVAR Cc = $"Cc"
	NVAR dE = $"dE"
	NVAR keV = $"keV"
	NVAR ap = $"ap"
	NVAR nk = $"nk"
	NVAR ds = $"ds"
	SetAberrations2()
	wave aber = $"aber"
	string out_total
	
	string p = PrintAberrations(aber)
	sprintf p, "%s; Cc = %g mm; dE = %g eV; ds = %g; ap = %g mrad; keV = %g keV", p, Cc, dE, ds, ap, keV

	if(!cmpstr(ctrlName, "psf1d_go"))
		if(Cc != 0 && dE != 0 || ds != 0)
			STEMPSF1DIncoh(aber, keV, Cc, dE, ds, ap, nk)
			wave out = $"probe1DIncoh"
		else
			STEMPSF1DCoh(aber, keV, ap, nk)
			wave out = $"probe1DCoh"
		endif
		out_total = fol+outn+"_psf1d"
	elseif(!cmpstr(ctrlName, "ctf1d_go"))
		if(Cc != 0 && dE != 0 || ds != 0)
			STEMCTF1DIncoh(aber, keV, Cc, dE, ds, ap, nk)
			wave out = $"CTF1DIncoh"
		else
			STEMCTF1DCoh(aber, keV, ap, nk)
			wave out = $"CTF1DCoh"
		endif
		out_total = fol+outn+"_ctf1d"
	elseif(!cmpstr(ctrlName, "psf2d_go"))
		if(Cc != 0 && dE != 0 || ds != 0)
			STEMPSF2DIncoh(aber, keV, Cc, dE, ds, ap, nk)
			wave out = $"probe2DIncoh"
		else
			STEMPSF2DCoh(aber, keV, ap, nk)
			wave out = $"probe2DCoh"
		endif
		out_total = fol+outn+"_psf2d"
	elseif(!cmpstr(ctrlName, "ctf2d_go"))
		if(Cc != 0 && dE != 0 || ds != 0)
			STEMCTF2DIncoh(aber, keV, Cc, dE, ds, ap, nk)
			wave out = $"CTF2DIncoh"
		else
			STEMCTF2DCoh(aber, keV, ap, nk)
			wave out = $"CTF2DCoh"
		endif
		out_total = fol+outn+"_ctf2d"
	else
		printf "Unknown option %s in STEMPanelGo.  Oops.\r", ctrlName
		return 0
	endif

	Note out, p
	Duplicate/O out $out_total
	Killwaves out
	SetDataFolder $fol

End


Function ZconGo(ba) : ButtonControl
	STRUCT WMButtonAction &ba

	switch( ba.eventCode )
		case 2: // mouse up

			string fol = GetDataFolder(1)
			ControlInfo xyz_wave_picker
			wave xyz = $(fol+S_value)
			if(!WaveExists(xyz))
				printf "Cannot find xyz wave %s", S_value
				return 0
			endif

			SetDataFolder "root:Packages:STEM:"			
			SVAR outn = $"outname"

			NVAR Cc = $"Cc"
			NVAR dE = $"dE"
			NVAR keV = $"keV"
			NVAR ap = $"ap"
			NVAR nk = $"nk"
			NVAR ds = $"ds"
			SetAberrations2()
			wave aber = $"aber"
			string out_total
	

			string p = PrintAberrations(aber)
			sprintf p, "%s; Cc = %g mm; dE = %g eV; ap = %g mrad; keV = %g keV", p, Cc, dE, ap, keV
			sprintf p, "xyz = %s; %s", NameOfWave(xyz), p

			zcon(xyz, aber, keV, Cc, dE, ds, ap)
			wave out = $"zcon_im"
			Note/K out, p
			
			out_total = fol+outn+"_zcon"
			
			Duplicate/O out, $out_total
			Killwaves out
			SetDataFolder $fol
			break
	endswitch

	return 0
End




Function ctfdf_go(ctrlName) : ButtonControl
	String ctrlName

	SVAR fol = $"start_folder"
	SVAR outn = $"outname"

	NVAR C5 = $"C5"
	NVAR Cs = $"Cs"
	NVAR Cc = $"Cc"
	NVAR dE = $"dE"
	NVAR keV = $"keV"
	NVAR ap = $"ap"
	NVAR nk = $"nk"
	NVAR kmax = $"kmax"
	NVAR df_min = $"df_min"
	NVAR df_max = $"df_max"
	NVAR df_steps = $"df_steps"
	
	Make/O/N=(nk, df_steps) ctfdf
	SetScale/I x 0, kmax, "", ctfdf
	SetScale/I y df_min, df_max, "", ctfdf
	//stemctf_df(ctfdf, C5, Cs, Cc, dE, keV, ap)
	string p
	sprintf p, "C5 = %g mm; Cs = %g mm; Cc = %g mm; dE = %g eV; ap = %g mrad; kV = %g kV", C5, Cs, Cc, dE, ap, keV
	Note ctfdf, p
	
	p = fol+outn+"_ctfdf"
	Duplicate/O ctfdf $p
	Killwaves ctfdf

End


Function STEM_Panel() : Panel

	wave test = $"root:Packages:STEM:aber"
	if(!WaveExists(test))
		STEMPanelSetup()
	endif

	PauseUpdate; Silent 1		// building window...
	NewPanel /W=(385,66,995,419) as "Convolution STEM"
	SetDrawLayer UserBack
	DrawText 110,43,"keV"
	DrawText 252,43,"mm"
	DrawText 252,65,"eV"
	DrawText 110,65,"mrad"
	SetDrawEnv fillpat= 0
	DrawRect 6,20,299,344
	DrawText 6,19,"Microscope Properties"
	DrawText 324,19,"Calculation"
	SetDrawEnv fillpat= 0
	DrawRect 323,20,599,344
	DrawText 209,91,"Å"
	SetVariable keV,pos={15,27},size={90,16},title="voltage"
	SetVariable keV,limits={0,inf,10},value= root:Packages:STEM:keV
	SetVariable Cc,pos={175,27},size={70,16},title="Cc"
	SetVariable Cc,limits={0,inf,0.1},value= root:Packages:STEM:Cc
	SetVariable dE,pos={175,49},size={70,18},title="\\F'Symbol'D\\F'MS Sans Serif'E"
	SetVariable dE,limits={0,inf,0.1},value= root:Packages:STEM:dE
	SetVariable ap,pos={15,49},size={90,16},title="aperture"
	SetVariable ap,help={"probe convergence half angle"}
	SetVariable ap,limits={0,inf,0.1},value= root:Packages:STEM:ap
	SetVariable outname,pos={336,299},size={250,16},title="output wave:"
	SetVariable outname,value= root:Packages:STEM:outname
	SetVariable psf1d_nk,pos={406,38},size={110,16},title="# of k points"
	SetVariable psf1d_nk,format="%d"
	SetVariable psf1d_nk,limits={50,inf,25},value= root:Packages:STEM:nk
	Button psf1d_go,pos={380,83},size={60,25},proc=STEMPanelGo,title="PSF 1D"
	Button psf1d_go,fSize=12
	ListBox aber,pos={25,98},size={255,213},frame=4
	ListBox aber,listWave=root:Packages:STEM:aber_text
	ListBox aber,selWave=root:Packages:STEM:aber_selwave,mode= 5
	Button default_aber,pos={112,316},size={80,20},proc=TitanDefaultAberrations2,title="Titan defaults"
	Button ctf1d_go,pos={380,125},size={60,25},proc=STEMPanelGo,title="CTF 1D"
	Button psf2d_go,pos={481,82},size={60,25},proc=STEMPanelGo,title="PSF 2D"
	Button ctf2d_go,pos={482,125},size={60,25},proc=STEMPanelGo,title="CTF 2D"
	Button zcon_go,pos={411,169},size={100,25},proc=ZconGo,title="Z Contrast Image"
	PopupMenu xyz_wave_picker,pos={332,220},size={258,21},bodyWidth=200,title="XYZ wave:"
	PopupMenu xyz_wave_picker,mode=1,popvalue="none",value= #"\"none;\"+Wavelist(\"*\",\";\",\"DIMS:2,MINCOLS:6,MAXCOLS:6\")"
	SetVariable ds,pos={84,75},size={120,16},title="source FWHM",format="%g"
	SetVariable ds,limits={0,inf,0.05},value= root:Packages:STEM:ds
EndMacro