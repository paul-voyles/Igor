#pragma rtGlobals=1		// Use modern global access method.

// Load FEI .bin type files.  File is just raw binary data with a 10-byte offset to the beginning
// of the data block.  Because there is no header, the user must specify the number of pixels
// and the data type.  Currently only single-precision (4-byte) floating point images are supported
//
// v1 pmv, 05/29/12


// Add load .bin menu item.
Menu "Load Waves"
		"Load TIA .bin file . . .", ReadTIABin()
		help = {"Load an FEI .bin file containing 4-byte floating point data."}
end

// Presents an open file dialog for the user to select a .bin file
// that is then read into a wave of the appropriate type.
// Currently only single-precision (4-byte) floating point images are supported
function ReadTIABin()
	
	variable fnum
	Open/R/D/T=".bin" fnum

	if(!strlen(S_fileName))
		return -1
	endif
	
	printf "Reading data from file %s\r", S_fileName
	
	variable px, py
	Prompt px,"x size:"
	Prompt py,"y size"
	DoPrompt "Enter image size in pixels",px, py
	if(V_Flag)
		Abort "The user pressed Cancel"
	endif	
	
	variable type = 1 // currently hard-coded for 4-byte float
	
	LoadTIABin(S_fileName, px, py, type)

	variable n = ItemsInList(S_filename, ":")
	S_filename = StringFromList(n-1, S_filename, ":")
	S_filename = StringFromList(0, S_filename, ".")
	
	Duplicate/O tia_read $S_filename
	Killwaves tia_read	

end

// Read a FEI .bin file into a wave of the appropriate
// size and data type.  file is a string with the path the file to open.
// px and py are the size of the image in pixels, and type is the data type
// currently only single precision (4-byte) floating point images are supported.
function LoadTIABin(file, px, py, type)
	string file
	variable px, py, type

	variable fnum
	// open the file, check that all is well
	Open/R fnum as file
	
	if(!strlen(S_fileName))
		return -1
	endif

	// Make a dummy wave for the input data.  Will be redimensioned to the right data type and size later.
 	Make/o/N=(px, py) tia_read
 	
 	FSetPos fnum, 10
 	FBinRead fnum, tia_read
 	
 	Close fnum
 	
end