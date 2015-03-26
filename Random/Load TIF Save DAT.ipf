#pragma rtGlobals=1		// Use modern global access method.


// The user selects a folder of TIFs and a place from the DATs to be saved.
// The function reads all the TIF images and saves them as DATs.
function LoadTIFFolderSaveDAT()
	
	NewPath/O/Q/M="Select a directory of TIF files" TIFpath
	PathInfo TIFpath
	string TIFname = S_Path
	if(!strlen(TIFname))
		return 0
	endif
	
	NewPath/O/Q/M="Select a directory to save DAT files" DATpath
	PathInfo DATpath
	string DATname = S_Path
	if(!strlen(DATname))
		return 0
	endif
	
	string base_name
	Prompt base_name, "Enter the base filename of outputed DAT files:"
	DoPrompt "Enter base name", base_name
	
	// Get a semicolon-separated list of all TIF files in the folder
	String list_TIF = IndexedFile(TIFpath, -1, ".TIF")
	// Sort list using combined alpha and numeric sort
	list_TIF = SortList(list_TIF, ";", 16)
	// Process the list
	Variable numItems = ItemsInList(list_TIF)

	string TIF_name, DAT_name, DAT_path
	variable i=0
	for(i=0; i<numItems; i+=1)
		TIF_name = StringFromList(i, list_TIF)
		DAT_name = StringFromList(i, list_TIF)

		if(!strlen(TIF_name))
			return 1
		endif 
		if(!strlen(DAT_name))
			return 1
		endif 
		
		//string temp_wave
		//sprintf temp_wave, "%s_%03d.dat", base_name, i
		
		ImageLoad/Q/T=tiff/P=TIFpath TIF_name
		wave image = $TIF_name
		//Rename image $temp_wave
		
		sprintf DAT_path, "%s%s_%03d.dat", DATname, base_name, i
		
		WriteIMIWork(image, DAT_path)
		
		killwaves image
	endfor
	
end

