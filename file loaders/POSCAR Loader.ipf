#pragma rtGlobals=3		// Use modern global access method and strict wave access.
#include "STEM"
//
// begun 04-17-15 pmv
//
// LoadPOSCAR: interactive loader with prompt to select a file and enter by hand supercell 
// parameters and Debye-Waller factors.  Call LoadPOSCARWork
//
// LoadPOSCARWork: non-interactive version which takes command line arguments or can
// be used programmatically.
//
// 04-19-15: major bug fix in coordinate loading.  minor improvement in dialog cancellation

// Add load POSCAR menu item
Menu "Load Waves"
	"Load POSCAR structure . . .", LoadPOSCAR()
	help = {"Load VASP POSCAR file into Igor / Kirkland XYZ wave."}
end

function LoadPOSCAR()

	variable fnum
	Open/R/D/T=".xyz" fnum

	if(!strlen(S_fileName))
		return -1
	endif
	
	variable ax, by, cz
	Prompt ax, "Enter ax in Angstroms:"
	DoPrompt "Supercell size", ax
	if(V_flag)
		return -1
	endif
	Prompt by, "Enter by in Angstroms:"
	DoPrompt "Supercell size", by
	if(V_flag)
		return -1
	endif
	Prompt cz, "Enter cz in Angstroms:"
	DoPrompt "Supercell size", cz
	if(V_flag)
		return -1
	endif
	
	make/o/n=(0, 2) dwf
	string s
	variable dw_z, dw_f
	Prompt s, "Enter Z and DWF separated by a comma, cancel when finished:"
	do
		DoPrompt "Debye-Waller factors", s
		if(V_flag)
			break
		endif
		InsertPoints 0, 1, dwf
		sscanf s, "%d, %g", dw_z, dw_f
		dwf[0][0] = dw_z
		dwf[0][1] = dw_f
		s = ""
	while(1)
		
	LoadPOSCARWork(S_fileName,ax, by, cz, dwf)

	variable n = ItemsInList(S_filename, ":")
	S_filename = StringFromList(n-1, S_filename, ":")
	S_filename = StringFromList(0, S_filename, ".")
	Rename poscar_xyz $S_filename

	killwaves dwf
end


function LoadPOSCARWork(file, ax, by, cz, dwf)
	string file
	variable ax, by, cz
	wave dwf
	
	LoadWave/D/J/L={0, 2,0,0,0}/V={" ","", 0, 0}/B="F=-2,N=Zs;N=fx;N=fy;N=fz;"/A/O file
	wave/T Zs=$"Zs"
	wave fx=$"fx"
	wave fy=$"fy"
	wave fz=$"fz"
	
	// strips trailing spaces from Zs
	variable i
	string s
	for(i=0; i<numpnts(Zs); i+=1)
		sscanf Zs[i], "%s", s
		Zs[i] = s
	endfor
	
	make/o/N=(numpnts(Zs), 6) poscar_xyz
	SetCell(poscar_xyz, ax, by, cz)
	
	poscar_xyz[][0] = SymboltoZ(Zs[p])
	poscar_xyz[][1] = fx[p]
	poscar_xyz[][2] = fy[p]
	poscar_xyz[][3] = fz[p]
	poscar_xyz[][4] = 1
	poscar_xyz[][5] = NaN
	
	for(i=0; i<DimSize(dwf, 0); i+=1)
		poscar_xyz[][5] = (poscar_xyz[p][0] == dwf[i][0] ? dwf[i][1] : poscar_xyz[p][5] )
	endfor

	wavestats/Q poscar_xyz
	if(V_numNaNs != 0)
		printf "Error.  %d unassigned Debye-Waller factors.\r", V_numNaNs
	endif
	
	Killwaves Zs, fx, fy, fz

end
