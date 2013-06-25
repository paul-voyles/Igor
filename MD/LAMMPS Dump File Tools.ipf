#pragma rtGlobals=1		// Use modern global access method.
#include "XYZ"

// functions to deal with LAMMPS dump files
//
// v0 begun 7-1-12, pmv 


// LAMMPS Menu
Menu "LAMMPS"
	"Load LAMMPS thermo data from an out file . . .", LMPLoadOut()
	help = {"Load thermo data from a LAMMPS output file . . ."}
	"Write an xyz wave to a LAMMPS structure file . . .", LMPWriteModelFile()
	help = {"Write an xyz structure wave to a LAMMPS .dat format text file.  User must provide atom type number for each atomic number."}
	"Read dump file timesteps . . .", LMPListTimeStep()
	help = {"List the timesteps in a LAMMPS dump file to the history area."}
	"Load one dump file configuration . . .", LMPReadConfig()
	help = {"Load one the model coordinates in one configuration from a LAMMPS dump file into an xyz wave."}
	"Write config for Voronoi analysis . . .", WriteVorFiles()
	help = {"Write the para.in and structure files needed for a Voronoi analysis of an xyz wave using the Ames Lab code."}
	"Write time series for Voronoi analysis . . ..", VorVTime()
	help = {"Write a series of para and structure files for Vornoi analysis of configurations from a single dump file."}
	
	
end


// lists timestamps in a dump file so the user can retrieve a particular configuration
function LMPListTimeStepWork(file)
	string file
	
	Grep/Q/INDX/E="ITEM: TIMESTEP" file 
	wave W_index = $"W_index"
	W_index += 1
	
	variable fnum, i=0, j=0
	Open/R fnum as file
	string l
	do
		FReadline fnum, l
		if(strlen(l) == 0)
			break
		endif
		if(i == W_index[j])
			printf "%s", l
			j+=1
		endif
		i+=1
	while(1)

	close fnum	
	
	Killwaves W_Index
end
	
// Menu item interface to LMPListTimeStepWork
function LMPListTimeStep()

	variable fnum
	Open/R/D/F="All Files:.*;" fnum
	
	if(!strlen(S_FileName))
		return -1
	endif
	
	LMPListTimeStepWork(S_Filename)
	
end

// utility function which returns the line number at the start of a particular configuration
// out of a dump file
function LMPFindTimeStep(file, timestep, spos)
	string file
	variable timestep, spos
	
	wave W_index = $"W_index"
	if(!WaveExists(W_index))
		Grep/Q/INDX/E="ITEM: TIMESTEP" file
		wave W_index = $"W_index"
		W_Index += 1
	endif
	
	variable fnum, i=0, j=0, n
	Open/R fnum as file
	FSetPos fnum, spos
	string l
	do
		FReadline fnum, l
		if(strlen(l) == 0)
			n=-1
			break
		endif
		if(i == W_index[j])
			if(str2num(l) == timestep)
				//n = i-1
				FStatus fnum
				n = V_filePos
				//printf "Found it on line %d\r", n
				break
			endif
			//printf "%s", l
			j+=1
		endif
		i+=1
	while(1)

	close fnum	
	
	return n

end

// Menu item interface to LMPReadConfigWork()
function LMPReadConfig()
	
	variable fnum
	Open/R/D/F="All Files:.*;" fnum
	
	if(!strlen(S_Filename))
		return -1
	endif
	
	variable timestep
	Prompt timestep, "Enter the timestep of the desired configuration: "
	DoPrompt "Timestep", timestep
	if(V_Flag)
		return -1
	endif
	
	LMPReadConfigWork(S_Filename, timestep, 1)
	printf "Read model from file %s\r", S_Filename
	
	Killwaves W_Index
end

// Loads a particular timestep from a LMP dump file.  Prompts the user to sets atomic 
// numbers if typeYN is non-zero
function LMPReadConfigWork(file, timestep, typeYN, [spos])
	string file
	variable timestep, typeYN, spos
		
	variable nstart, fnum, i, natoms, bm, bp, ax, by, cz
	variable id, type, xs, ys, zs
	string l
		
	nstart = LMPFindTimeStep(file, timestep, spos)
	if(nstart < 0)
		printf "Timestep %d not found in file %s.  Exiting.", timestep, file
		return -1
	endif

	Open/R fnum as file

	FSetPos fnum, nstart
	//for(i=0; i<nstart; i+=1)
	//	FReadline fnum, l
	//endfor
	
	//FReadline fnum, l  // should be text "ITEM: TIMESTEP"
	//FReadline fnum, l // should be number, selected time step
	FReadline fnum, l // should be text "ITEM: NUMBER OF ATOMS"
	FReadline fnum, l
	natoms = str2num(l)
	Make/O/N=(natoms, 6) lmp_xyz
	FReadline fnum, l // should be text "ITEM: BOX BOUNDS pp pp pp"
	FReadline fnum, l // x box, negative, positive
	sscanf l, "%g %g", bm, bp
	ax = bp - bm
	FReadline fnum, l // y box, negative, positive
	sscanf l, "%g %g", bm, bp
	by = bp - bm
	FReadline fnum, l // z box, negative, positive
	sscanf l, "%g %g", bm, bp
	cz = bp - bm
	SetCell(lmp_xyz, ax, by, cz)

	FReadline fnum, l // "ITEM: ATOMS id, type, xs, ys, zs"
	for(i=0; i<natoms; i+=1)
		FReadline fnum, l
		sscanf l, "%g %g %g %g %g", id, type, xs, ys, zs
		lmp_xyz[i][0] = type
		lmp_xyz[i][1] = xs*ax
		lmp_xyz[i][2] = ys*by
		lmp_xyz[i][3] = zs*cz
	endfor
	
	// periodic boundary conditions to get all the atoms back inside the box
	lmp_xyz[][1] = (lmp_xyz[p][1] < 0 ? lmp_xyz[p][1] + ax : lmp_xyz[p][1])
	lmp_xyz[][1] = (lmp_xyz[p][1] > ax ? lmp_xyz[p][1] - ax : lmp_xyz[p][1])
	lmp_xyz[][2] = (lmp_xyz[p][2] < 0 ? lmp_xyz[p][2] + by : lmp_xyz[p][2])
	lmp_xyz[][2] = (lmp_xyz[p][2] > by ? lmp_xyz[p][2] - by : lmp_xyz[p][2])
	lmp_xyz[][3] = (lmp_xyz[p][3] < 0 ? lmp_xyz[p][3] + cz : lmp_xyz[p][3])
	lmp_xyz[][3] = (lmp_xyz[p][3] > cz ? lmp_xyz[p][3] - cz : lmp_xyz[p][3])

	// set occupancy	
	lmp_xyz[][4] = 1		
	
	if(typeYN)
		// count the atom types
		variable ntypes=1, itype
		Make/o/n=(natoms) atypes = lmp_xyz[p][0]
		sort atypes, atypes
		itype = atypes[0]
		for(i=1; i<natoms; i+=1)
			if(atypes[i] != itype)
				ntypes += 1
				itype = atypes[i]
			endif
		endfor
		make/o/n=(ntypes) typeZ
		for(i=0; i<ntypes; i+=1)
			sprintf l, "Enter the atomic number for atom type %d", i+1
			Prompt itype, l
			DoPrompt "Z numbers", itype
			typeZ[i] = itype
		endfor
	
		for(i=0; i<ntypes; i+=1)
			lmp_xyz[][0] = (lmp_xyz[p][0] == i+1 ? typeZ[i] : lmp_xyz[p][0]) 
		endfor
		Killwaves atypes, typeZ
	endif
	
	FStatus fnum
	
	close fnum
	
	return V_filePos
end


//
function LMPWriteModelFile()

	variable fnum
	Open/D/F="All Files:.*;" fnum
	
	if(!strlen(S_Filename))
		return -1
	endif

	variable i
	Prompt i, "Select a model wave:", popup, WaveList("*", ";", "DIMS:2,MAXCOLS:6,MINCOLS:6")
	DoPrompt "Model selection", i
	wave xyz = $StringFromList(i-1, WaveList("*", ";", "DIMS:2,MAXCOLS:6,MINCOLS:6"))

	string n
	Prompt n, "Enter the model comment string:"
	DoPrompt "Model comment", n

	// count the atom types
	variable ntypes=1, itype
	variable natoms = DimSize(xyz, 0)
	Make/o/n=(natoms) atypes = xyz[p][0]
	Make/o/n=1 Zlist
	sort atypes, atypes
	itype = atypes[0]
	Zlist = itype
	for(i=1; i<natoms; i+=1)
		if(atypes[i] != itype)
			ntypes += 1
			itype = atypes[i]
			InsertPoints 0, 1, Zlist
			Zlist[0] = itype
		endif
	endfor
	
	string p	
	make/o/n=(ntypes) Ztypes
	Ztypes = p+1
	for(i=0; i<ntypes; i+=1)
		sprintf p, "Enter the atom type number for atomic number %d", ZList[i]
		itype = Ztypes[i]
		Prompt itype, p
		DoPrompt "Z numbers", itype
		Ztypes[i] = itype
	endfor
	
	Sort Ztypes, Zlist
	
	LMPWriteModelFileWork(S_Filename, n, xyz, ZList)
	
	Killwaves atypes, Zlist, Ztypes
end


// Ztypes is a list of the z-numbers in the model.  The entry in the list is atom type 1,
// the second is atom type 2, etc. in the output LAMMPS file.
function LMPWriteModelFileWork(file, name, xyz, ztypes)
	string file, name
	wave xyz, ztypes
	
	variable f
	Open/Z f as file
	if(V_flag)
		printf "Cannot open file %s.  Exiting.\r", file
		return 0
	endif
	
	variable natoms = dimsize(xyz, 0)
	Duplicate/O xyz xyz_temp
	Recenter(xyz, 0, 0, 0)
	Make/o/N=(natoms) num, Ztype_list, xx, yy, zz
	num = p+1
	xx = xyz[p][1]
	yy = xyz[p][2]
	zz = xyz[p][3]
	
	// Assign atom types
	Ztype_list = nan
	variable i
	for(i=0; i<numpnts(ztypes); i+=1)
		Ztype_list = (xyz[p][0] == ztypes[i] ? i+1 : Ztype_list[p] )
	endfor
	if(numtype(sum(Ztype_list)))
		printf "Some Z number not found in ztypes input.  Exiting.\r"
		close f
		return 0
	endif
	
	fprintf f, "%s\n", name
	fprintf f, "\n"
	fprintf f, "    %d atoms\n", natoms  
	fprintf f, "     %d atom types\n", numpnts(ztypes)
	fprintf f, "%g %g xlo xhi\n", -0.5*GetAX(xyz), 0.5*GetAX(xyz)
	fprintf f, "%g %g ylo yhi\n", -0.5*GetBY(xyz), 0.5*GetBY(xyz)
	fprintf f, "%g %g zlo zhi\n", -0.5*GetCZ(xyz), 0.5*GetCZ(xyz)
	fprintf f, "\n"
	fprintf f, "\n"
	fprintf f, "Atoms\n"
	fprintf f, "\n"
	
	wfprintf f, "%d %d %g %g %g \n", num, Ztype_list, xx, yy, zz
	fprintf f, "\n" 
	
	close f
	
	Killwaves xyz_temp, Ztype_list, xx, yy, zz, num
	
end

// Menu item interface to WriteVorPara and WriteVorStructureFileWork
function WriteVorFiles()

	string xyzn
	variable steps, cut
	
	Prompt xyzn, "Select an xyz wave", popup Wavelist("*", ";", "DIMS: 2,MAXCOLS:6,MINCOLS:6")
	DoPrompt "", xyzn
	if(V_flag)
		return -1
	endif
	wave xyz = $xyzn
	
	steps = 1  // for now hard-code to one configuration
	cut = 3.7
	
	Prompt cut, "Enter the nearest-neighbor cut-off distance: "
	DoPrompt "", cut
	if(V_flag)
		return -1
	endif
	
	variable f
	Open/D/F="All Files:.*;" f
	if(!strlen(S_filename))
		return -1
	endif
	
	string base = StringFromList(ItemsInList(S_filename, ":")-1, S_filename, ":")
	base = StringFromList(0, base, ".")
	string paras = RemoveFromList(StringFromList(ItemsInList(S_filename, ":") -1, S_filename, ":"), S_filename, ":")
	paras += (base+".vorp")
	
	WriteVorPara(paras, xyz, steps, cut)
	WriteVorStructureFileWork(S_filename, xyz)

end

// Writes the para.in file required by the Ames Lab Voronoi analysis program, as
// modified by Jinwoo Hwang
function WriteVorPara(file, xyz, steps, cut)
	string file
	wave xyz
	variable steps, cut
	
	variable fnum
	Open/F="All Files:.*;" fnum as file
	if(!strlen(S_filename))
		printf "Cannot open file %s for writing.  Exiting.\r", file
		return -1
	endif
	
	if( GetAX(xyz) != GetBY(xyz) || GetAX(xyz) != GetCZ(xyz) )
		printf "Voronoi code probably only works for cubic boxes.  Proceed at your own risk.\r"
	endif
	
	variable natoms = DimSize(xyz, 0)
	
	variable i, j, itype
	make/n=8/o type_count = 0
	make/o/n=(natoms) atype
	atype = xyz[p][0]
	sort atype atype
	itype = atype[0]
	j = 0
	for(i=0; i<natoms; i+=1)
		if(atype[i] == itype)
			type_count[j] += 1
		else
			itype = atype[i]
			j+=1
			type_count[j] += 1
		endif
	endfor
	
	fprintf fnum, "# steps, # total atoms, # atoms of type 1, # atoms of type 2\n"
	fprintf fnum, "%d \t %d \t %d \t %d \t\n", steps, natoms, type_count[0], type_count[1]
	fprintf fnum, "# box size, # neighbor cut-off\n"
	fprintf fnum, "%g \t %g\n", GetAX(xyz), cut
	
	close fnum
	
	Killwaves type_count, atype
end
	
// Writes a structure file in the format required by the Ames Lab Voronoi analysis program,
// as modified by Jinwoo Hwang	
function WriteVorStructureFileWork(file, xyz)
	wave xyz
	string file
	
	Duplicate/O xyz vor_xyz
	Recenter(vor_xyz, 0, 0, 0)

	//sort atoms by atomic number
	variable npts = DimSize(xyz, 0)
	Make/O/N=(npts) zn, xa, ya, za
	zn = vor_xyz[p][0]
	xa = vor_xyz[p][1]
	ya = vor_xyz[p][2]
	za = vor_xyz[p][3]

	Sort zn, zn, xa, ya, za

	xa = xa/(0.5*GetAX(xyz))
	ya = ya/(0.5*GetBY(xyz))
	za = za/(0.5*GetCZ(xyz))
	
	
	variable f
	Open/T=".vor" f as file
	if(!strlen(S_filename))
		printf "Error opening file.\r"
		return 0
	endif
	
	wfprintf f, "%g\t %g\t %g\n" xa, ya, za

	close f
	
	Killwaves vor_xyz, zn, xa, ya, za
	
end

function VorVTime()

	variable fnum
	Open/R/D/F="All Files:.*;" fnum
	
	if(!strlen(S_Filename))
		return -1
	endif

	string base
	variable start_t, stop_t, tstep, cut

	NewPath/O/M="Select the output folder" ptemp
	PathInfo ptemp
	string outp = S_Path
	
	Prompt base, "Enter the base name (no file extension) for the files: "
	DoPrompt "", base
	if(V_flag)
		return -1
	endif
	
	Prompt start_t, "Enter the first time step: "
	DoPrompt "", start_t
	if(V_flag)
		return -1
	endif
	
	Prompt stop_t, "Entering the last time step: "
	DoPrompt "", stop_t
	if(V_flag)
		return -1
	endif
	
	Prompt tstep, "Enter the # of time steps between calculations: "
	DoPrompt "", tstep
	if(V_flag)
		return -1
	endif

	Prompt cut, "Enter the nearest-neighbor cut-off distance: "
	DoPrompt "", cut
	if(V_flag)
		return -1
	endif
	

	VorVTimeWork(S_Filename, outp, base, start_t, stop_t, tstep, cut)
end


function VorVTimeWork(dfile, outp, base, start_t, stop_t, tstep, cut)
	string dfile, outp, base
	variable start_t, stop_t, tstep, cut
	
	string curfol = GetDataFolder(1)
	
	NewDataFolder/O/S root:Packages
	NewDataFolder/O/S LAMMPS
	NewDataFolder/O/S vor_temp	
	
	string w
	
	variable fnum
	sprintf w, "%s%s_vor.sh", outp, base
	Open fnum as w
	if(!strlen(S_Filename))
		return -1
	endif
	fprintf fnum, "#!/bin/bash\n"
	
	variable t, spos = 0
	string this_out
	for(t = start_t; t<=stop_t; t+=tstep)
	
		spos = LMPReadConfigWork(dfile, t, 0, spos = spos)
		if(spos==-1)
			return -1
		endif
		wave xyz = $"lmp_xyz"
	
		sprintf this_out, "%s_t%d", base, t

		sprintf w, "%s%s.vorm", outp, this_out
		WriteVorStructureFileWork(w, xyz)
	
		sprintf w, "%s%s.vorp", outp, this_out
		WriteVorPara(w, xyz, 1, cut)
	
		fprintf fnum, "vor %s.vorp %s.vorm %s_vor\n", this_out, this_out, this_out 
	
	endfor
	
	close fnum
	
	SetDataFolder $curfol
	KillDataFolder root:Packages:LAMMPS:vor_temp
	
end


function LMPLoadOut()

	variable f
	Open/R/D/F="output files (*.out):.out;" f
	if(!strlen(S_Filename))
		return 0
	endif
	
	string base
	Prompt base, "Base wave name for thermo data:"
	DoPrompt "", base
	if(V_Flag)
		return -1
	endif

	LMPLoadOutWork(S_Filename, base)
	
end

function LMPLoadOutWork(ofile, base)
	string ofile, base

	string curfol = GetDataFolder(1)
	NewDataFolder/O/S root:Packages
	NewDataFolder/O/S LAMMPS
	NewDataFolder/O/S ltemp

		
	LoadWave/G/A/B="N=step;N=temp;N=press;N=vol;N=pot;N=tot;" ofile

	string temp = WaveList("temp*", ";", "")
	concatenate/np/o WaveList("step*", ";", ""), $(base+"_step")
	concatenate/np/o WaveList("temp*", ";", ""), $(base+"_temp")
	concatenate/np/o WaveList("press*", ";", ""), $(base+"_press")
	concatenate/np/o WaveList("vol*", ";", ""), $(base+"_vol")
	concatenate/np/o WaveList("pot*", ";", ""), $(base+"_pot")
	concatenate/np/o WaveList("tot*", ";", ""), $(base+"_tot")

	Duplicate/O $(base+"_step") $(curfol+base+"_Step")
	Duplicate/O $(base+"_temp") $(curfol+base+"_Temp")
	Duplicate/O $(base+"_press") $(curfol+base+"_Press")
	Duplicate/O $(base+"_vol") $(curfol+base+"_Vol")
	Duplicate/O $(base+"_pot") $(curfol+base+"_PotEng")
	Duplicate/O $(base+"_tot") $(curfol+base+"_TotEng")
	
	SetDataFolder $curfol
	KillDataFolder root:Packages:LAMMPS:ltemp
end