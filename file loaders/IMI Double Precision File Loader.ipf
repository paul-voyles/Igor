#pragma rtGlobals=1		// Use modern global access method.
#pragma rtGlobals=1		// Use modern global access method.
//
// IMI file loader
//
// added load folder, 6/21/13 pmv

// Add load SER menu item.
Menu "Load Waves"
		"Load IMI file . . .", ReadIMI()
		help = {"Load an IMI .dat file."}
		"Load folder of IMI files . . . ", LoadIMIFolder()
		help = {"Load all the IMI .dat files in a the selected folder."}
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
