#pragma rtGlobals=1		// Use modern global access method.

//Functions for reading and writing Gatan Fixed Format (.gfx) image files.
//
// ReadGFX():		read a single GFX file with open dialog
// LoadGFXFolder():	automatically load an entire folder of GFX files
// LoadGFX(file):		utility function that actually reads the file
// WriteGFX():		save a 2D wave containing an image or spectrum profile to a GFX file
// WriteGFXSI():		save a spectrum image imported from a TIA .ser to a GFX file
// WriteGFXWork():	utility routine that actually writes the GFX file format
//
// 09/08/10 pmv:
//	v2, at least.  Added save GFX file functions, renamed to file to GFX.ipf from ReadGFX.ipf 

// Add load gfx menu item.
Menu "Load Waves"
		"Load GFX file . . .", ReadGFX()
		help = {"Load a DigitalMicrograph fixed format .gfx file."}
		"Load folder of GFX files . . . ", LoadGFXFolder()
		help = {"Load all the DigitalMicrograph fixed format .gfx files in a the selected folder."}
end

// Add save GFX menu item.
Menu "Save Waves"
	"Save image or Spectrum Profile to GFX . . .", WriteGFX()
	help = {"Save a 2D data set (image or spectrum profile) to fixed format .gfx file."}
	"Save Spectrum Image to GFX . . .", WriteGFXSI()
	help = {"Save a 3D data set (spectrum image) to a fixed format .gfx file."}
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
		gfx_path = dirname  + gfx_name
		LoadGFX(gfx_path)
		gfx_name = StringFromList(0, gfx_name, ".")
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

// Save a 2D data set, like an image or a spectrum profile to a
// Gatan fixed format .gfx file.
function WriteGFX()

	string w
	Prompt w, "Select a 2D Wave", popup, WaveList("*", ";", "DIMS:2")
	DoPrompt "Select a Wave", w
	
	if(V_flag)
		return -1
	endif
	
	wave dat = $w

	variable fnum
	Open/D/T=".gfx" fnum
	
	if(!strlen(S_filename))
		return -1
	endif

	WriteGFXWork(dat, S_filename)
	
end

// Save a Spectrum Image imported from a .ser file to a Gatan fixed
// format file.
function WriteGFXSI()

	string w
	prompt w, "Select a 2D Wave", popup, WaveList("*", ";", "DIMS:2")
	DoPrompt "Select a Wave", w
	
	if(V_flag)
		return -1
	endif
	
	wave dat = $w
	wave xpos = $w+"_xpos"
	if(!WaveExists(xpos))
		printf "xpos wave %s not found.  Exiting.\r", w+"_xpos"
		return -1
	endif
	wave ypos = $w+"_ypos"
	if(!WaveExists(ypos))
		printf "ypos wave %s not found.  Exiting.\r", w+"_ypos"
		return -1
	endif
	
	if(WaveType(dat) != 2)
		printf "Only implemented for single precision floating point waves.  Exiting.\r"
		return -1
	endif
	
	variable ep, xp, yp
	ep = DimSize(dat, 0)
	FindLevel/Q/R=(1,inf) xpos, xpos[0]
	yp = V_LevelX
	xp = numpnts(xpos) / yp
	
	//printf "Writing (%d, %d, %d) spectrum image.\r", xp, yp, ep
	
	make/o/n=(xp, yp, ep) gfx_si_temp
	variable n, i, j
	n = 0
	for(i=0; i<xp; i+=1)
		for(j=0; j<yp; j+=1)
			gfx_si_temp[i][j][] = dat[r][n]
			n+=1
		endfor
	endfor
	
	make/o/n=(yp, xp, ep) gfx_si_new
	for(i=0; i<xp; i+=1)
		for(j=0; j<yp; j+=1)
			gfx_si_new[j][i][] = gfx_si_temp[i][j]
		endfor
	endfor
	printf "Writing (%d, %d, %d) spectrum image.\r", yp, xp, ep
	
	variable fnum
	Open/D/T=".gfx" fnum
	
	if(!strlen(S_filename))
		return -1
	endif

	WriteGFXWork(gfx_si_new, S_filename)
	killwaves gfx_si_temp
	killwaves gfx_si_new
	
end

// Utility routine for WriteGFX and WriteGFXSI.  Takes the wave in dat and writes it
// to a Gatan fixed format file "file".  Currently only implemented for single-precision
// floating point data.
function WriteGFXWork(dat, file)
	wave dat
	string file
	
	variable fnum
	Open/T=".gfx" fnum as file
	
	if(!strlen(S_fileName))
		return -1
	endif
	
	variable endian = 65535, size_x, size_y, size_z, byte_depth, enum_type
	
	size_x = DimSize(dat, 0)
	size_y = max(DimSize(dat, 1), 1)
	size_z = max(DimSize(dat, 2), 1)
	//byte_depth = 
	//enum_type = 

	variable wtype = WaveType(dat)
	variable npnts = size_x*size_y*size_z
	// printf "%d total points.\r", npnts
	
	FBinWrite/F=3/U fnum, endian
	FBinWrite/F=3 fnum, size_x
	FBinWrite/F=3 fnum, size_y
	FBinWrite/F=3 fnum, size_z

	
	make/o gfx_write
	
	switch(WaveType(dat))

		// real, 4-byte data
		case 2:
			byte_depth = 4
			enum_type = 2
			Redimension/R/N=(npnts) gfx_write
			gfx_write = dat
			//FBinWrite/F=3 fnum, byte_depth
			//FBinWrite/F=3 fnum, enum_type
			//FBinWrite/F=4 fnum, gfx_write		
		break

		default:
			printf "Sorry, only single-precision floating point currently supported.\r"
		break
		
	endswitch

	FBinWrite/F=3 fnum, byte_depth
	FBinWrite/F=3 fnum, enum_type
	FBinWrite/F=0 fnum, gfx_write

	close fnum
	
	Killwaves gfx_write
end

// Read a Gatan Fixed Format image file into a 2D wave of the appropriate
// size and data type.  file is a string with the path the file to open.
function LoadGFX(file)
	string file

	variable fnum	
	Open/R fnum as file
	
	if(!strlen(S_fileName))
		return -1
	endif
 	
	variable endian, size_x, size_y, size_z, byte_depth, enum_type
	
	FBinRead/F=3/U fnum, endian
	FBinRead/F=3 fnum, size_x
	FBinRead/F=3 fnum, size_y
	FBinRead/F=3 fnum, size_z
	FBinRead/F=3 fnum, byte_depth
	FbinRead/F=3 fnum, enum_type
		
	Make/O gfx_read

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
			printf "Load a (%d, %d) signed, 4-byte integer image from file %s\r", size_x, size_y, S_fileName
			Redimension/Y=0x20/N=(size_x, size_y) gfx_read
			FBinRead fnum, gfx_read
			r=1
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
