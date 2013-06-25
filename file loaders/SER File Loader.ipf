#pragma rtGlobals=1		// Use modern global access method.
//
// SER file loader
//
// based on file structure description in "Emispec Support - ES Vision Series File.mht"
// provided by FEI
// begun 07-22-09, pmv; based on ReadGFX.ipf
// v0.1, without description & unit strings, or spectrum image position tags, 07-23-09
// v1.0, fully functional 08-05-09

// Add load SER menu item.
Menu "Load Waves"
		"Load SER file . . .", ReadSER()
		help = {"Load an FEI .ser file."}
		"Load folder of SER files . . . ", LoadSERFolder()
		help = {"Load all the FEI .ser files in a the selected folder."}
end

// The user selects a folder and the function reads all
// the SER format image files in that folder into 2D waves.
function LoadSERFolder()

	NewPath/O/M="Select a directory of ser files" dirpath
	PathInfo dirpath
	string dirname = S_Path
	if(!strlen(dirname))
		return 0
	endif
	
	string ser_name, ser_path
	variable i=0
	do
		ser_name = IndexedFile(dirpath, i, ".ser")
		if(!strlen(ser_name))
			return 1
		endif 
		ser_path = dirname + ser_name
		LoadSER(ser_path)
		ser_name = StringFromList(0, ser_name, ".")
		Duplicate/O ser_read $ser_name
		Killwaves ser_read
	
		wave tmp = $"acquisition_time"
		if(WaveExists(acquisition_time))
			Duplicate/O tmp $(ser_name+"_ac_time")
			Killwaves tmp
		endif
		wave tmp = $"x_si_pos"
		if(WaveExists(x_si_pos))
			Duplicate/O tmp $(ser_name + "_xpos")
			Killwaves tmp
		endif
		wave tmp = $"y_si_pos"
		if(WaveExists(y_si_pos))
			Duplicate/O tmp $(ser_name+"_ypos")
			Killwaves tmp
		endif

		i+=1
	while(1)
	
end


// Presents an open file dialog for the user to select a SER file
// that is then read into a wave of the appropriate type.
function ReadSER()
	
	variable fnum
	Open/R/D/T=".ser" fnum

	if(!strlen(S_fileName))
		return -1
	endif
	
	printf "Reading data from file %s\r", S_fileName
	
	LoadSER(S_fileName)

	variable n = ItemsInList(S_filename, ":")
	S_filename = StringFromList(n-1, S_filename, ":")
	S_filename = StringFromList(0, S_filename, ".")
	
	Duplicate/O ser_read $S_filename
	Killwaves ser_read
	
	wave tmp = $"acquisition_time"
	if(WaveExists(acquisition_time))
		Duplicate/O tmp $(S_filename+"_ac_time")
		Killwaves tmp
	endif
	wave tmp = $"x_si_pos"
	if(WaveExists(x_si_pos))
		Duplicate/O tmp $(S_filename + "_xpos")
		Killwaves tmp
	endif
	wave tmp = $"y_si_pos"
	if(WaveExists(y_si_pos))
		Duplicate/O tmp $(S_filename+"_ypos")
		Killwaves tmp
	endif
	

end

// Read a FEI SER file into a wave of the appropriate
// size and data type.  file is a string with the path the file to open.
function LoadSER(file)
	string file

	// internal function variables, counters, etc.
 	variable fnum, iset, tmp, cpos, xpos, ypos, tag_offset_array_pos
 	string n
 	 
 	// variables describing the file being loaded
 	variable dimension		// 1 for 1D data like spectra, 2 for 2D data like images
 	variable tag_type			// 1 for time-series or time-stamp; 2 for a time & position series like a specturm image
 	variable nsets			// # of data sets in the series.  1 for just an image of spectrum, >1 for a spectrum image
 	variable offset_array_pos	// position in the file of the array that contains the starting positions for the data & tags
 	variable dtime			// time stamp data, seconds since Jan. 1, 1970
	string description	= ""		// whatever is stored in the DimensionArray description
	string unit_string	= ""		// whatever is stored in the DimensionArray units

	// open the file, check that all is well
	Open/R fnum as file
	
	if(!strlen(S_fileName))
		return -1
	endif

	// Make a dummy wave for the input data.  Will be redimensioned to the right data type and size later.
 	Make/o ser_read
 	
 	
 	// Read the file header
 	FBinRead/F=2 fnum, tmp	// ByteOrder, should always be litte-endian
 	FBinRead/F=2 fnum, tmp  // SeriesID, should always be 0x0197
 	FBinRead/F=2 fnum, tmp  // SeriesVersion, should always be 0x0210
 
 	FBinRead/F=3 fnum, tmp  // DataTypeID, determines the dimension of the fundamental data set
  	if(tmp == 16674)
 		dimension = 2
 		printf "2D data, "
 	elseif (tmp == 16672)
  		dimension = 1
 		printf "1D data, "
 	else
 		printf "Error.  Unknown data dimensionality.\r"
 		return 0
 	endif
 	
 	FBinRead/F=3 fnum, tmp  // TagTypeID, determines wheter the tag is just time, or position & time
 	if(tmp == 16722)
 		tag_type = 1
 		printf "time-stamped or time-series data only, "
 	elseif (tmp == 16706)
 		tag_type = 2
 		printf "time- and position-stamped spectrum image, "
 	else
 		printf "TagTypeID has unknown value.  Exiting.\r"
 		return 0
 	endif
 	
 	FBinRead/F=3 fnum, tmp		// TotalNumberElements, which is sometimes incorrect
 	FBinRead/F=3 fnum, nsets  	// ValidNumberElements, which seems more reliable
 	if (nsets == 1)
 		printf "single data set.\r"
 	else
 		printf "%d data sets.\r", nsets
 	endif

 	FBinRead/F=3 fnum, offset_array_pos
 	FBinRead/f=3 fnum, tmp		// NumberDimensions
 	

 	// Read the dimension array.  Don't need most of this - just the description and unit strings
 	FBinRead/f=3 fnum, tmp	// DimensionSize
 	FBinRead/f=5 fnum, tmp	// CalibrationOffset
	FBinRead/f=5 fnum, tmp	// CalibrationDelta
 	FBinRead/f=3 fnum, tmp	// CalibrationElement
 	FBinRead/f=3 fnum, tmp	// DescriptionLength
 	description = PadString(description, tmp, 0)
 	FBinRead fnum, description
 	FBinRead/f=3 fnum, tmp	// UnitsLength
 	unit_string = PadString(description, tmp, 0)
 	FBinRead fnum, unit_string


 	// Read the data 	
 	FSetPos fnum, offset_array_pos	// jump to beginning of the offset array
 	
 	// Single image
 	if (dimension == 2 && nsets == 1)
 		FBinRead/F=3 fnum, cpos		// read the image file position
		FStatus fnum
 		tag_offset_array_pos = V_FilePos
 		ReadOneImage(fnum, cpos)	// load it
 		Duplicate/O ser_oneimage ser_read
 		Killwaves ser_oneimage
 	
 	// Single spectrum
 	elseif (dimension == 1 && nsets == 1)
 		FBinRead/F=3 fnum, cpos			// read the spectrum file position
 		FStatus fnum
 		tag_offset_array_pos = V_FilePos
 		ReadOneSpectrum(fnum, cpos)	// load it
 		Duplicate/O ser_onespec ser_read
 		Killwaves ser_onespec
 	
 	// Spectrum image of images
 	elseif (dimension == 2 && nsets > 1)
 		for(iset = 0; iset<nsets; iset+=1)
 			FSetPos fnum, offset_array_pos + 4*iset
 			FbinRead/F=3 fnum, cpos
			if( iset == (nsets-1))
				FStatus fnum
				tag_offset_array_pos = V_FilePos
			endif
 			ReadOneImage(fnum, cpos)
 			
 			if(iset == 0)
 				Duplicate/O ser_oneimage ser_read
 				Redimension/N=(dimsize(ser_oneimage, 0), dimsize(ser_oneimage, 1), nsets) ser_read
 			endif
 			wave ser_oneimage = $"ser_oneimage"
 			ser_read[][][iset] = ser_oneimage[p][q]
 		endfor
 		killwaves ser_oneimage
 	
 	// Spectrum image of spectra
 	elseif (dimension == 1 && nsets > 1)
 		for(iset = 0; iset<nsets; iset+=1)
 			FSetPos fnum, offset_array_pos + 4*iset
 			FbinRead/F=3 fnum, cpos
			if( iset == (nsets-1))
				FStatus fnum
				tag_offset_array_pos = V_FilePos
			endif
 			ReadOneSpectrum(fnum, cpos)
 			
 			if(iset == 0)
 				Duplicate/O ser_onespec ser_read
 				Redimension/N=(numpnts(ser_onespec), nsets) ser_read
 			endif
 			wave ser_onespec = $"ser_onespec"
 			ser_read[][iset] = ser_onespec[p]
 		endfor
 		killwaves ser_onespec
 	
 	
 	else
 		printf "Something very strange has happened.\r"
 		printf "Trying to load a data set with dimension = %d and nsets = %d.\r", dimension, nsets
 		printf "Doesn't work.  Exiting.\r"
 		return 0
 	endif
 	
 	// Read the tags and apply an appropriate wave scaling
	if (nsets == 1)
	 	FSetPos fnum, offset_array_pos + 4*nsets
		FBinRead/f=3 fnum, cpos
		FSetPos fnum, cpos
		if( tag_type == 1)
			FBinRead/f=3 fnum, tmp		// TagTypeID again
			FBinRead/f=3 fnum, dtime  	// acquisition time in seconds since 1/1/1970
			dtime += date2secs(1970, 1, 1)
			n = note(ser_read)
			n = ReplaceStringByKey("acquisition_time", n, (secs2time(dtime, 0)+" "+secs2date(dtime, 0)), "=", ";")
			Note/K ser_read, n
		else // tag_type = 2
			printf "Single data set, tag_type 2, position & time.  This code is untested.\r"
			FBinRead/f=3 fnum, tmp		// TagTypeID again
			FBinRead/f=3 fnum, dtime
			FBinRead/f=5 fnum, xpos
			FBinRead/f=5 fnum, ypos
			n = note(ser_read)
			n = ReplaceNumberByKey("position_x", n, xpos, "=", ";")
			n = ReplaceNumberByKey("position_y", n, ypos, "=", ";")
			Note/K ser_read, n
		endif
	
	elseif (nsets > 1)
		Make/o/n=(nsets) acquisition_time
		SetScale/P x 0, 1, "dat", acquisition_time
		if(tag_type == 2)
			Make/O/N=(nsets) x_si_pos, y_si_pos
		endif
			
		for(iset=0; iset<nsets; iset+=1)
		 	FSetPos fnum, tag_offset_array_pos + 4*iset
 			FBinRead/F=3 fnum, cpos
 			FSetPos fnum, cpos
	 		FBinRead/f=3 fnum, tmp
 			FBinRead/f=3 fnum, dtime	
			dtime += date2secs(1970, 1, 1)
  			acquisition_time[iset] = dtime
			if(iset == 0)
				n = note(ser_read)
				n = ReplaceStringByKey("acquisition_time", n, (secs2time(dtime, 0)+" "+secs2date(dtime, 0)), "=", ";")
				Note/K ser_read, n
			endif

 			if(tag_type == 2)
 				FBinRead/f=5 fnum, xpos
 				FBinRead/f=5 fnum, ypos
 				x_si_pos[iset] = xpos
 				y_si_pos[iset] = ypos
			endif
  		endfor
	
	
	else
		printf "nsets = %d.  Shouldn't be able to get here.  Exiting.\r"
		return 0
	endif

	// Store descriptions and notes fields in wave note.
 	n = note(ser_read)
	n = ReplaceStringByKey("description", n, description, "=", ";")
	n = ReplaceStringByKey("units", n, unit_string, "=", ";")
	Note/K ser_read, n
	
 	close fnum

end


function ReadOneImage(fnum, dataarray_startpos)
	variable fnum, dataarray_startpos
	
	variable cal_startx, cal_deltax, cal_elx, cal_starty, cal_deltay, cal_ely
	variable data_type, size_x, size_y
	
	Make/O ser_oneimage

 	FSetPos fnum, dataarray_startpos
 	FBinRead/F=5 fnum, cal_startx
 	FBinRead/F=5 fnum, cal_deltax
 	FBinRead/F=3 fnum, cal_elx
 	FBinRead/F=5 fnum, cal_starty
 	FBinRead/F=5 fnum, cal_deltay
 	FBinRead/F=3 fnum, cal_ely
 	FBinread/F=2 fnum, data_type
 	FBinRead/F=3 fnum, size_x
 	FBinRead/F=3 fnum, size_y
 		
 	switch(data_type)
 		
 			case 1:  // unsigned, 1-btye integer
				//printf "Data type unsigned, 1-byte integer.\r"
				Redimension/U/B/N=(size_x, size_y) ser_oneimage
				FBinRead fnum, ser_oneimage
 				break
 				
 			case 2: // unsigned, 2-btye integer
				//printf "Data type unsigned, 2-byte integer.\r"
				Redimension/U/W/N=(size_x, size_y) ser_oneimage
				FBinRead fnum, ser_oneimage
 				break
 				
 			case 3: // unsigned, 4-byte integer
				//printf "Data type unsigned, 4-byte integer.\r"
				Redimension/U/I/N=(size_x, size_y) ser_oneimage
				FBinRead fnum, ser_oneimage
 				break
 				
 			case 4: // signed, 1-byte integer
				//printf "Data type signed, 1-byte integer.\r"
				Redimension/B/N=(size_x, size_y) ser_oneimage
				FBinRead fnum, ser_oneimage
 				break
 								
 			case 5: // signed 2-byte integer
				//printf "Data type signed, 2-byte integer.\r"
				Redimension/W/N=(size_x, size_y) ser_oneimage
				FBinRead fnum, ser_oneimage
 				break
 				 			
 			case 6: // signed 4-byte integer
				//printf "Data type signed, 4-byte integer.\r"
				Redimension/I/N=(size_x, size_y) ser_oneimage
				FBinRead fnum, ser_oneimage
 				break
 			
 			case 7: // 4-byte float
				//printf "Data type 4-byte float.\r"
				Redimension/S/N=(size_x, size_y) ser_oneimage
				FBinRead fnum, ser_oneimage
 				break
 			
 			case 8: // 8-byte float
				//printf "Data type 8-byte float.\r"
				Redimension/D/N=(size_x, size_y) ser_oneimage
				FBinRead fnum, ser_oneimage
 				break
 			
 			case 9:  // 8-byte complex
				//printf "Data type 8-byte complex.\r"
				Redimension/S/C/N=(size_x, size_y) ser_oneimage
				FBinRead fnum, ser_oneimage
 				break
 			
 			case 10: // 16-byte complex
				//printf "Data type 16-byte complex.\r"
				Redimension/C/D/N=(size_x, size_y) ser_oneimage
				FBinRead fnum, ser_oneimage
 				break
 			
 			default:
 				printf "Unknown data type in file.  Exiting.\r"
 				return 0
 		endswitch
 		
 		SetScale/P x (cal_startx - cal_deltax*cal_elx), cal_deltax, "", ser_oneimage
 		SetScale/P y (cal_starty - cal_deltay*cal_ely), cal_deltay, "", ser_oneimage		
end	


function ReadOneSpectrum(fnum, dataarray_startpos)
	variable fnum, dataarray_startpos

	variable cal_start, cal_delta, cal_el
	variable data_type, size
	
	Make/O ser_onespec

 	FSetPos fnum, dataarray_startpos
 	FBinRead/F=5 fnum, cal_start
 	FBinRead/F=5 fnum, cal_delta
 	FBinRead/F=3 fnum, cal_el
 	FBinread/F=2 fnum, data_type
 	FBinRead/F=3 fnum, size
 		
 	switch(data_type)
 		
 			case 1:  // unsigned, 1-btye integer
				//printf "Data type unsigned, 1-byte integer.\r"
				Redimension/U/B/N=(size) ser_onespec
				FBinRead fnum, ser_onespec
 				break
 				
 			case 2: // unsigned, 2-btye integer
				//printf "Data type unsigned, 2-byte integer.\r"
				Redimension/U/W/N=(size) ser_onespec
				FBinRead fnum, ser_onespec
 				break
 				
 			case 3: // unsigned, 4-byte integer
				//printf "Data type unsigned, 4-byte integer.\r"
				Redimension/U/I/N=(size) ser_onespec
				FBinRead fnum, ser_onespec
 				break
 				
 			case 4: // signed, 1-byte integer
				//printf "Data type signed, 1-byte integer.\r"
				Redimension/B/N=(size) ser_onespec
				FBinRead fnum, ser_onespec
 				break
 								
 			case 5: // signed 2-byte integer
				//printf "Data type signed, 2-byte integer.\r"
				Redimension/W/N=(size) ser_onespec
				FBinRead fnum, ser_onespec
 				break
 				 			
 			case 6: // signed 4-byte integer
				//printf "Data type signed, 4-byte integer.\r"
				Redimension/I/N=(size) ser_onespec
				FBinRead fnum, ser_onespec
 				break
 			
 			case 7: // 4-byte float
				//printf "Data type 4-byte float.\r"
				Redimension/S/N=(size) ser_onespec
				FBinRead fnum, ser_onespec
 				break
 			
 			case 8: // 8-byte float
				//printf "Data type 8-byte float.\r"
				Redimension/D/N=(size) ser_onespec
				FBinRead fnum, ser_onespec
 				break
 			
 			case 9:  // 8-byte complex
				//printf "Data type 8-byte complex.\r"
				Redimension/S/C/N=(size) ser_onespec
				FBinRead fnum, ser_onespec
 				break
 			
 			case 10: // 16-byte complex
				//printf "Data type 16-byte complex.\r"
				Redimension/C/D/N=(size) ser_onespec
				FBinRead fnum, ser_onespec
 				break
 			
 			default:
 				printf "Unknown data type in file.  Exiting.\r"
 				return 0
 		endswitch
 		
 		SetScale/P x (cal_start - cal_delta*cal_el), cal_delta, "", ser_onespec
end	
	
	