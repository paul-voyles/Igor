#pragma rtGlobals=1		// Use modern global access method.


function SeparateRAWseries()
	
	NewPath/O/Q/M="Select a directory of RAW files" dirpath
	PathInfo dirpath
	string dirname = S_Path
	if(!strlen(dirname))
		return 0
	endif
	
	NewPath/O/Q/M="Select a directory of RPL files" rplpath
	PathInfo rplpath
	string rplname = S_Path
	if(!strlen(rplname))
		return 0
	endif
	
	string raw_name, raw_path, channel_num, channel_path, file, rpl_name, rpl_path
	
	string base_name
	Prompt base_name, "Enter the base output filename:"
	DoPrompt "Enter base name", base_name
	
	// Get a semicolon-separated list of all RAW files in the folder
	String list_raw = IndexedFile(dirpath, -1, ".RAW")
	// Sort list using combined alpha and numeric sort
	list_raw = SortList(list_raw, ";", 16)
	// Process the list
	Variable numItems = ItemsInList(list_raw)
	
	// Get a semicolon-separated list of all RPL files in the folder
	String list_rpl = IndexedFile(rplpath, -1, ".RPL")
	// Sort list using combined alpha and numeric sort
	list_rpl = SortList(list_rpl, ";", 16)
	
	variable i=0
	variable j=0
	for(i=0; i<numItems; i+=1)
		//Load RAW 3D image
		raw_name = StringFromList(i, list_raw)
		rpl_name = StringFromList(i, list_rpl)
		
		if(!strlen(raw_name))
			return 1
		endif 
		if(!strlen(rpl_name))
			return 1
		endif 
		
		raw_path = dirname  + raw_name
		rpl_path = rplname  + rpl_name
		print "...Loading, separating, and saving image ", i+1, " of ", numItems
		LoadRAW(raw_path, rpl_path)
		raw_name = StringFromList(0, raw_name, ".")
		raw_name = "EDS_" + raw_name
		Rename raw_read $raw_name
		wave im = $raw_name
		
		//if i=0 open rpl and load size_z and make directories of that size
		if (i==0)
			variable size_z = DimSize(im, 2)
			MakeChannelDirectories(size_z)
			string/G ChannelDir
		endif
		
		//seperate and save each channel image in DAT format in corresponding channels folder
		j=0
		for(j=0; j<size_z; j+=1)
			
			imagetransform/p=(j) getplane im
			wave slice =$"M_ImagePlane"
			
			channel_num = "channel" + num2istr(j) + ":"
			channel_path =  ChannelDir + channel_num
			
			sprintf file, "%s%s_%03d.dat", channel_path, base_name, i
			
			WriteIMIWork(slice, file)
			
			killwaves M_ImagePlane
		endfor
		
		//kill image to make room for more (IGOR could not handle many EDS spectrum images
		killwaves $raw_name
	endfor
	killstrings ChannelDir
	killpath ChannelSeries
end



function MakeChannelDirectories(num)
	
	variable num
	
	NewPath/O/Q/M="Select a directory for the Channel Directories" channelpath
	PathInfo channelpath
	string channelname = S_Path
	if(!strlen(channelname))
		return 0
	endif
	
	string/G ChannelDir
	ChannelDir = channelname + "ChannelSeries:"
	NewPath/C/Q ChannelSeries, ChannelDir
	
	string raw_name, channel_num, channel_path
	
	variable i=0
	variable j=0
	
	//makes each channels directory in the ChannelSeries directory
	do
		raw_name = IndexedFile(dirpath, i, ".RAW")
		if(!strlen(raw_name))
			return 1
		endif 
		
		j=0
		for(j=0; j<num; j+=1)
			
			channel_num = "channel" + num2istr(j) + ":"
			channel_path =  ChannelDir + channel_num
			NewPath/C/Q channel_num, channel_path
			
			killpath channel_num
		endfor	
		i+=1
	while(1)
end



function ReassembleAndSum()
	
	//define spectrum image size
	variable size_x = 256
	variable size_y = 256
	variable size_z = 2048
	
	//identify what images in series failed to NRR in order to leave those images out of average
	variable numIgnores = 2		//# of image to be ignored. do not comment out if not needed, just make this = 0.
	variable ignoreimage1 = 2		//first image to ignore. comment out if not needed.
	variable ignoreimage2 = 39		//second image to ignore. comment out if not needed.
	
	//select channel directory
	NewPath/O/Q/M="Select a directory that contains the Channel Directories" channelpath
	PathInfo channelpath
	string channelname = S_Path
	if(!strlen(channelname))
		return 0
	endif
	string ChannelDir = channelname
	
	//make directory for reassembled EDS series
	NewPath/O/Q/M="Select a directory for the NRR EDS series" nrredspath
	PathInfo nrredspath
	string nrredsname = S_Path
	if(!strlen(nrredsname))
		return 0
	endif
	string nrredsDir
	nrredsDir = nrredsname + "NRR_EDS:"
	NewPath/C/Q NRR_EDS, nrredsDir
	
	string base_name
	Prompt base_name, "Enter the base filename of input files:"
	DoPrompt "Enter base name", base_name
	
	string channel_num
	string channel_path = ChannelDir +"channel0:"
	NewPath/O/Q InitialPath, channel_path
	
	// Get a semicolon-separated list of all RAW files in the folder
	String list = IndexedFile(InitialPath, -1, ".DAT")
	// Sort list using combined alpha and numeric sort
	list = SortList(list, ";", 16)
	// Process the list to get total number of images
	Variable numItems = ItemsInList(list)
	
	//make wave for the sum spectrum image
	make/o/n=(size_x, size_y, size_z) sum_wave
	sum_wave = 0
	
	string file, filesave
	
	variable i=0
	variable j=0
	//load each spectrum image chnnel by channel, average the spectrum images, and save as an IGOR binary file
	for(i=0; i<numItems; i+=1)
	//for(i=0; i<2; i+=1)
		if(i==ignoreimage1 || i==ignoreimage2)
			print "ignored image", i+1
		endif
		if(i!=ignoreimage1 && i!=ignoreimage2)
			//make wave for image i
			sprintf file, "%s_%03d.dat", base_name, i
			make/o/n=(size_x, size_y, size_z) $file
			wave im =$file
			sprintf filesave, "%s_%03d", "EDS", i
			
			for(j=0; j<size_z; j+=1)
				channel_num = "channel" + num2istr(j) + ":"
				channel_path =  ChannelDir + channel_num + file
			
				//load DAT file
				LoadIMI(channel_path)
	
				Duplicate/O imi_read slice
				Killwaves imi_read
				redimension/s slice
				
				//insert DAT slice into correct plane of of EDS SI
				imagetransform /INSI=slice /INSX=0 /INSY=0 /P=(j) insertImage im
			
				killwaves slice
			endfor
			
			//set inf in image series to zero, so summing works
			InfToZero(im)
			
			//save EDS spectrum image # i
			Save/O/C/P=NRR_EDS im as filesave
			//updates sum spectrum image
			sum_wave = sum_wave + im
		
			killwaves im
			print "completed image", i+1, " of ", numItems
		endif
	endfor
end



