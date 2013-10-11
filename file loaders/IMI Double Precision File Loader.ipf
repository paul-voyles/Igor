#pragma rtGlobals=1		// Use modern global access method.
#pragma rtGlobals=1		// Use modern global access method.
//
// IMI file loader
//
// added load folder, 6/21/13 pmv
// added save file, 7/8/13 pmv
// added save 3D image stack to separate files, 7/9/13 pmv
// added save menu items, 7/9/13 pmv
// fixed end-of-line \r vs \n bug 10/11/13 pmv

// Add load IMI menu item.
Menu "Load Waves"
	"Load IMI file . . .", ReadIMI()
	help = {"Load an IMI .dat file."}
	"Load folder of IMI files . . . ", LoadIMIFolder()
	help = {"Load all the IMI .dat files in a the selected folder."}
end

// Add save IMI menu items.
Menu "Save Waves"
	"Save IMI file . . .", WriteIMI()
	help = {"Save a 2D wave to an IMI .dat file."}
	"Save series to IMI files . . .", WriteIMIStack()
	help = {"Save a 3D wave to a series of IMI .dat files."}
end

// Presents an open file dialog for the user to select a .dat file
// that is then read into a wave of the appropriate type.
function ReadIMI()
	
	variable fnum
	Open/R/D/T=".dat" fnum

	if(!strlen(S_fileName))
		return -1
	endif
	
	printf "Reading data from file %s\r", S_fileName
	
	LoadIMI(S_fileName)

	variable n = ItemsInList(S_filename, ":")
	S_filename = StringFromList(n-1, S_filename, ":")
	S_filename = StringFromList(0, S_filename, ".")
	
	Duplicate/O imi_read $S_filename
	Killwaves imi_read

end

// Read an IMI .dat file into a wave of the appropriate
// size and data type.  file is a string with the path the file to open.
function LoadIMI(file)
	string file

	variable fnum
	
	// open the file, check that all is well
	Open/R fnum as file
	
	if(!strlen(S_fileName))
		return -1
	endif

	// Make a dummy wave for the input data.  Will be redimensioned to the right data type and size later.
 	Make/o imi_read
 	
 	string l
 	variable nx, ny
 	
 	// Read the file header
 	FReadLine fnum, l
 	FReadLine fnum, l
 	printf "comment: %s\r", l
 	FReadline fnum, l
 	sscanf l, "%d %d", nx, ny
 	FReadline fnum, l
 	//FBinRead/f=1 fnum, t
 	//FBinRead fnum, t
 	
 	Redimension/D/N=(nx, ny) imi_read
 	FBinRead/B=3 fnum, imi_read
 
 	close fnum

end


// The user selects a folder and the function reads all
// the IMI format image files in that folder into 2D waves.
function LoadIMIFolder()

	NewPath/O/M="Select a directory of .dat files" dirpath
	PathInfo dirpath
	string dirname = S_Path
	if(!strlen(dirname))
		return 0
	endif
	
	string gfx_name, gfx_path
	variable i=0
	do
		gfx_name = IndexedFile(dirpath, i, ".dat")
		if(!strlen(gfx_name))
			return 1
		endif 
		gfx_path = dirname  + gfx_name
		LoadIMI(gfx_path)
		gfx_name = StringFromList(0, gfx_name, ".")
		Rename imi_read $gfx_name
		i+=1
	while(1)
	
end

// Save a 2D data set, like an image or a spectrum profile to a
// Gatan fixed format .gfx file.
function WriteIMI()

	string w
	Prompt w, "Select a 2D Wave", popup, WaveList("*", ";", "DIMS:2")
	DoPrompt "Select a Wave", w
	
	if(V_flag)
		return -1
	endif
	
	wave dat = $w

	variable fnum
	Open/D/T=".dat" fnum
	
	if(!strlen(S_filename))
		return -1
	endif

	WriteIMIWork(dat, S_filename)
	
end

function WriteIMIWork(im, file)
	wave im
	string file
	
	variable fnum
	Open/T=".dat" fnum as file
	
	if(!strlen(S_fileName))
		return -1
	endif
	
	fprintf fnum, "P9\n"
	fprintf fnum, "#written from image %s by Igor Pro, %s\n", NameOfWave(im), date()
	fprintf fnum, "%d %d\n", DimSize(im, 0), DimSize(im, 1)
	
	wavestats/Q im
	fprintf fnum, "%d\n", round(V_max)
	
	if(wavetype(im) == 0x04)
		FBinWrite/f=5 fnum, im
	else
		Duplicate/O im, imi_write_temp
		Redimension/D imi_write_temp
		FBinWrite/f=5 fnum, imi_write_temp
		Killwaves imi_write_temp
	endif
	
	Close fnum
	
end

function WriteIMIStack()

	string w
	Prompt w, "Select a 3D Wave", popup, WaveList("*", ";", "DIMS:3")
	DoPrompt "Select a Wave", w
	
	if(V_flag)
		return -1
	endif
	
	wave dat = $w

	NewPath/O/M="Select a destination folder" dirpath
	PathInfo dirpath
	string dirname = S_Path
	if(!strlen(dirname))
		return 0
	endif

	Prompt w, "Enter the base output filename:"
	DoPrompt "Enter base name", w

	WriteIMIStackWork(dat, w, dirname)

end

function WriteIMIStackWork(im, base_file, fol)
	wave im
	string base_file, fol
	
	if(DimSize(im, 2) > 1000)
		printf "Cannot currently export more than 1000 images at a time because I'm too lazy to program it.  Sorry.\r"
		return -1
	endif
	
	variable i
	string file
	for(i=0; i<DimSize(im, 2); i+=1)
		sprintf file, "%s%s_%03d.dat", fol, base_file, i
		ImageTransform/P=(i) getPlane im
		wave one = $"M_ImagePlane"
		WriteIMIWork(one, file)
	endfor
	
	Killwaves M_ImagePlane
	
end
	