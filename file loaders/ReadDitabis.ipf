#pragma rtGlobals=1		// Use modern global access method.

// Add load IPC menu item.
Menu "Load Waves"
		"Load DITABIS image . . .", ReadIPM()
		help = {"Load a Ditabis IPC image file."}
end


// Presents an open file dialog for the user to select a IPC file
// that is then read into a 2D wave of the appropriate type.
function ReadIPM()
	
	variable fnum
	Open/R/D/T=".ipc.ipl.iph.ipr" fnum

	if(!strlen(S_fileName))
		return -1
	endif
	
	if(LoadIPM(S_fileName))
		variable n = ItemsInList(S_filename, ":")
		S_filename = StringFromList(n-1, S_filename, ":")
		S_filename = StringFromList(0, S_filename, ".")
		Duplicate/O ipc $S_filename
		Killwaves ipc
	endif
	
end



function LoadIPM(file)
	string file
	
	// open the file
	variable f
	Open/R f as file
	
	if(!strlen(S_filename))
		return -1
	endif
	
	string n
	variable ns = ItemsInList(S_filename, ":")
	S_filename = StringFromList(ns-1, S_filename, ":")
	sprintf n, "file=%s", S_filename
	
	// read and parse the text header
	variable header, xpixel, ypixel, data_depth, xres, yres, thumb_zoom
	variable mag, offset, offset_yn, gain, laser
	string line, t, val
	header = -1
	xpixel = -1
	ypixel = -1
	data_depth = -1
	
	do
		FReadline f, line
		if(!cmpstr(line, ""))
			printf "FReadline failed, which is very wierd.\r"
			close f
			return -1
		endif

		t = StringFromList(0, line, " ")
		if(!cmpstr(t, "COMMENT"))
			break
		endif

		t = StringFromList(0, line, "=")
		
		if(!cmpstr(t, "CREATED "))
			val = line[9,strlen(line)-4]
			sprintf n, "%s;date=%s", n, val
		endif
		if(!cmpstr(t, "HEADER "))
			sscanf line, "HEADER = %g", header
		endif
		if(!cmpstr(t, "XPIXEL "))
			sscanf line, "XPIXEL = %g", xpixel
		endif
		if(!cmpstr(t, "YPIXEL "))
			sscanf line, "YPIXEL = %g", ypixel
		endif
		if(!cmpstr(t, "BYTE PER PIXEL "))
			sscanf line, "BYTE PER PIXEL = %g", data_depth
		endif
		if(!cmpstr(t, "XRESOLUTION "))
			sscanf line, "XRESOLUTION = %g", xres
		endif
		if(!cmpstr(t, "YRESOLUTION "))
			sscanf line, "BYTE PER PIXEL = %g", yres
		endif
		if(!cmpstr(t, "THUMB-NAIL-ZOOM "))
			sscanf line, "THUMB-NAIL-ZOOM = %g", thumb_zoom
		endif
		if(!cmpstr(t, "OFFSET "))
			sscanf line, "OFFSET = %g", offset
			sprintf n, "%s;offset=%g", n, offset
		endif
		if(!cmpstr(t, "OFFSET CORRECTION "))
			sscanf line, "OFFSET CORRECTION = %g", offset_yn
		endif
		if(!cmpstr(t, "MAGNIFICATION "))
			sscanf line, "MAGNIFICATION = %g", mag
			sprintf n, "%s;mag=%g", n, mag 
		endif
		if(!cmpstr(t, "GAIN "))
			sscanf line, "GAIN = %g", gain
			sprintf n, "%s;gain=%g", n, gain
		endif
	
	
	while(1)
	
	if(header == -1)
		printf "Required header tag HEADER not found.  Cannot proceed.\r"
		close f
		return 0
	endif
	if(xpixel == -1)
		printf "Required header tag XPIXEL not found.  Cannot proceed.\r"
		close f
		return 0
	endif
	if(ypixel == -1)
		printf "Required header tag YPIXEL not found.  Cannot proceed.\r"
		close f
		return 0
	endif
	if(data_depth == -1)
		printf "Required header tag BYTE PER PIXEL not found.  Cannot proceed.\r"
		close f
		return 0
	endif
	
	printf "Loaded a (%d, %d) pixel, %d bit image.\r", xpixel, ypixel, data_depth*8
	
	// Read the image into a wave
	FSetPos f, header
	Make/O/N=(ypixel,xpixel) ipc

	if(data_depth == 2)
		Redimension/U/W ipc
		FBinRead/b=3/U/f=2 f, ipc 
	elseif(data_depth == 3)
		printf "3-byte data depth - I don't think the read will work properly.\r"
		Redimension/U/I ipc
		FBinRead/b=3/U/f=2 f, ipc
	elseif(data_depth == 4)
		Redimension/U/I ipc
		FBinRead/b=3/U/f=3 f, ipc
	endif
	
	Note ipc, n
	close f
	
	return 1
end