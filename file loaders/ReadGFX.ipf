#pragma rtGlobals=1		// Use modern global access method.

//Functions for reading Gatan Fixed Format (.gfx) image files into Igor 2D waves.
//
// ReadGFX():		read a single GFX file with open dialog
// LoadGFXFolder():	automatically load an entire folder of GFX files
// LoadGFX(file):		utility function that actually reads the file

// Add load gfx menu item.
Menu "Load Waves"
		"Load GFX file . . .", ReadGFX()
		help = {"Load a DigitalMicrograph fixed format .gfx file."}
		"Load folder of GFX files . . . ", LoadGFXFolder()
		help = {"Load all the DigitalMicrograph fixed format .gfx files in a the selected folder."}
end

// The user selects a folder and the function reads all
// the GFX format image files in that folder into 2D waves.
function LoadGFXFolder()

	NewPath/O/M="Select a directory of gfx files" dirpath
	PathInfo dirpath
	string dirname = S_Path
	if(!strlen(dirname))
		return 0
	endif
	
	string gfx_name, gfx_path
	variable i=0
	do
		gfx_name = IndexedFile(dirpath, i, ".gfx")
		if(!strlen(gfx_name))
			return 1
		endif 
		gfx_path = dirname + gfx_name
		
		LoadGFX(gfx_path)
		//gfx_name = StringFromList(0, gfx_name, ".")
		gfx_name = "image_"+  num2str(i)
		Rename gfx_read $gfx_name
		i+=1
	while(1)
	
end


// LoadGFXFolder but with a path passed as a parameter
function LoadGFXFixedFolder(dirname)
	string dirname
	
	NewPath/O/Q dirpath, dirname

	string gfx_name, gfx_path
	variable i=0
	do
		gfx_name = IndexedFile(dirpath, i, ".gfx")
		if(!strlen(gfx_name))
			return 1
		endif 
		gfx_path = dirname + ":" + gfx_name
		

		LoadGFX(gfx_path)
		
		gfx_name = StringFromList(0, gfx_name, ".")
				print gfx_name
		
		Rename gfx_read $gfx_name
		
		i+=1
	while(1)
	
end



// Presents an open file dialog for the user to select a GFX file
// that is then read into a 2D wave of the appropriate type.
function ReadGFX()
	
	variable fnum
	Open/R/D/T=".gfx" fnum

	if(!strlen(S_fileName))
		return -1
	endif
	
	LoadGFX(S_fileName)

	variable n = ItemsInList(S_filename, ":")
	S_filename = StringFromList(n-1, S_filename, ":")
	S_filename = StringFromList(0, S_filename, ".")
	Rename gfx_read $S_filename

end

// Read a Gatan Fixed Format image file into a 2D wave of the appropriate
// size and data type.  file is a string with the path the file to open.
// procedure changed such that pixels from a failed job would return NaN instead of the number from previous job. cz 11-28-16
function LoadGFX(file)
	string file

	variable fnum
	Open/R fnum as file	
 	
	variable endian, size_x, size_y, size_z, byte_depth, enum_type
	
	FBinRead/F=3/U fnum, endian
	FBinRead/F=3 fnum, size_x
	FBinRead/F=3 fnum, size_y
	FBinRead/F=3 fnum, size_z
	FBinRead/F=3 fnum, byte_depth
	FbinRead/F=3 fnum, enum_type
		
//	if(size_z != 1)
//		printf "3D images are not supported.  No wave loaded.\r"
//		Close fnum
//		return -1
//	endif	
		
	Make/O gfx_read
	
	if(!strlen(S_fileName))
		Redimension/W/N=(size_x, size_y) gfx_read
		gfx_read=0
		return -1
	endif

	variable r = 0
	switch (enum_type)
	
		case 0:
			printf "Data type NULL_DATA is not supported.  No wave loaded.\r"
			break
		
		case 1: // signed 16-bit int
			if(size_z==1)
				printf "Loading a (%d, %d) signed integer 2-byte image from file %s\r", size_x, size_y, S_fileName
				Redimension/W/N=(size_x, size_y) gfx_read
				FBinRead fnum, gfx_read
				r=1
			else
				printf "Loading a (%d, %d, %d) signed integer 2-byte image from file %s\r", size_x, size_y, size_z, S_fileName
				Redimension/W/N=(size_x, size_y, size_z) gfx_read
				FBinRead fnum, gfx_read
				r=1				
			endif
			break
	
		case 2:	// real 4 byte
			if(size_z==1)
				printf "Loading a (%d, %d) real 4 byte image from file %s\r", size_x, size_y, S_fileName
				Redimension/R/N=(size_x, size_y) gfx_read
				FBinRead fnum, gfx_read
				r=1
			else
				printf "Loading a (%d, %d, %d) real 4 byte image from file %s\r", size_x, size_y, size_z, S_fileName
				Redimension/R/N=(size_x, size_y, size_z) gfx_read
				FBinRead fnum, gfx_read
				r=1			
			endif	
			break

		case 3:
			printf "Data type COMPLEX8_DATA is not supported.  No wave loaded.\r"		
			break
			
		case 4:
			printf "Data typeOBSELETE_DATA is not supported.  No wave loaded.\r"		
			break

		case 5:
			printf "Data type PACKED_DATA is not supported.  No wave loaded.\r"		
			break
		
		case 6: // unsigned 8-bit int
			printf "Loading a (%d, %d) unsigned 1-byte image from file %s\r", size_x, size_y, S_fileName
			Redimension/U/B/N=(size_x, size_y) gfx_read
			FBinRead fnum, gfx_read
			r=1
			break
	
		case 7:
			printf "Data type SIGNED_INT32_DATA is not supported.  No wave loaded.\r"		
			break		

		case 8:
			printf "Data type RGB_DATA is not supported.  No wave loaded.\r"		
			break		

		case 9:
				printf "Data type SIGNED_INT8_DATA is not supported.  No wave loaded.\r"		
			break		

		case 10:	// unsigned 16-bit int
			printf "Loading a (%d, %d) unsigned 2-byte image from file %s\r", size_x, size_y, S_fileName
			Redimension/W/U/N=(size_x, size_y) gfx_read
			FBinRead fnum, gfx_read
			r=1
			break	
		
		case 11:
			printf "Data type UNSIGNED_INT32_DATA is not supported.  No wave loaded.\r"		
			break		

		case 12:
				printf "Data type REAL8_DATA is not supported.  No wave loaded.\r"		
			break		

		case 13:
			printf "Data type COMPLEX16_DATA is not supported.  No wave loaded.\r"		
			break		

		case 14:
			printf "Data type BINARY_DATA is not supported.  No wave loaded.\r"		
			break		

		case 15:
			printf "Data type RGBA_FLOAT32_DATA is not supported.  No wave loaded.\r"		
			break		

		case 16:
			printf "Data type RGB_UINT16_DATA is not supported.  No wave loaded.\r"		
			break		

		case 17:
			printf "Data type RGB_FlOAT64_DATA is not supported.  No wave loaded.\r"		
			break		

		case 18:
			printf "Data type RGBA_FLOAT64_DATA is not supported.  No wave loaded.\r"		
			break		

		case 19:
			printf "Data type RGBA_UINT16_DATA is not supported.  No wave loaded.\r"		
			break		

		case 20:
			printf "Data type RGB_UINT8_DATA is not supported.  No wave loaded.\r"		
			break		

		case 21:
			printf "Data type RGBA_UINT8_DATA is not supported.  No wave loaded.\r"		
			break		

		case 22:
			printf "Data type LAST_DATA is not supported.  No wave loaded.\r"		
			break		

		case 23:
			printf "Data type OS_RGBA_UINT8_DATA is not supported.  No wave loaded.\r"		
			break		

		default:
			printf "Do not recognize data type %g.\r", enum_type
	
	endswitch
	
	// DM puts (0,0) in the top left.  Igor puts in the bottom left.  Reverse the y coordinates
	// of the wave so it looks right on the screen.

//	if(r)
//		Duplicate/O gfx_read gfx_tmp
//		gfx_read = gfx_tmp[p][(DimSize(gfx_read, 1) - q)]
//		KillWaves gfx_tmp	
//	endif
		
	Close fnum
	
end 


// gfx data types from http://www.gatan.com/~software/importexport.html
//	enum	type
//	--------	------
//	0		NULL_DATA
//	1		SIGNED_INT16_DATA
//	2		REAL4_DATA
//	3		COMPLEX8_DATA
//	4		OBSELETE_DATA
//	5		PACKED_DATA
//	6		UNSIGNED_INT8_DATA
//	7		SIGNED_INT32_DATA
//	8		RGB_DATA 
//	9		SIGNED_INT8_DATA
//	10		UNSIGNED_INT16_DATA
//	11		UNSIGNED_INT32_DATA
//	12		REAL8_DATA
//	13		COMPLEX16_DATA
//	14		BINARY_DATA
//	15		RGBA_FLOAT32_DATA
//	16		RGB_UINT16_DATA
//	17		RGB_FLOAT64_DATA
//	18		RGBA_FLOAT64_DATA
//	19		RGBA_UINT16_DATA
//	20		RGB_UINT8_DATA
//	21		RGBA_UINT8_DATA
//	22		LAST_DATA
//	23		OS_RGBA_UINT8_DATA = RGB_DATA
