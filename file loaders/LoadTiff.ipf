#pragma rtGlobals=3		// Use modern global access method and strict wave access.
// #include "gfx"
#include <All IP Procedures>
#include <Image Saver>

// functions to read a folder of tif files into a single 3D matrix.  Developed for use with ECM experiments.
// v1 from Li He / Pei Zhang 2015
//
// v2 pmv 8/13/15
//    Added user interface function LoadTiffFolder.  Previous utility function rename LoadTiffFolderWork
//    new UI function lets user select the folder, loads the first tif to get pixel dimensions, and count the
//    number of tiff files automatically.
//    
//   Dramatic speed up in loading large numbers of files.  The problem was in calling IndexedFile to get
//    the file name of every image.  That requires an OS call to get every file name, and the OS has to
//    count all the files in the folder from the first one to get to the nth file.  Thus, the process scales as
//    n^2.  Now, we use one IndexedFile call to get a list of all the files in the directory at once, stored in
//    a string, then pull names in order from the string.  It's all inside Igor, so avoids the OS calls, and it
//    scales at n.

function LoadTiffFolder()
	
	NewPath/O/Q/M="Select a directory of gfx files" dirpath
	PathInfo dirpath
	string dirname = S_Path
	if(!strlen(dirname))
		return 0
	endif
	
	string tiff_name, tiff_path
	tiff_name = IndexedFile(dirpath, 0, ".tiff")
	if(!strlen(tiff_name))
		printf "First file in the directory is not a tif.  Exiting.\r"
		return 1
	endif 
	tiff_path = dirname  + tiff_name

	imageload/q/t=tiff/N=name tiff_path
	wave test = $"name"
	variable nx, ny, nframes
	nx = DimSize(test, 0)
	ny = DimSize(test, 1)
	Killwaves test
	
	string all_tiffs = IndexedFile(dirpath, -1, ".tiff")
	nframes = ItemsInList(all_tiffs)
		
	printf "Loading %d (%d x %d) pixel images from directory %s.\r", nframes, nx, ny, dirname	
	LoadTiffFolderWork(dirname, nx, ny, nframes)
end	

function LoadTiffFolderWork(dirname, x_dimension, y_dimension, framenumber)
	string dirname
	variable x_dimension, y_dimension, framenumber
	
	print time()+"\r"
	NewPath/O/Q dirpath, dirname

	make/o/n=(x_dimension, y_dimension, framenumber) data
	if(!WaveExists(data))
		printf "Cannot allocate memory to data wave.  Exiting.\r"
		return 0
	endif
	
	string all_tiffs, tiff_name, tiff_path
	all_tiffs = IndexedFile(dirpath, -1, ".tiff")
	
	variable i=0
	printf "Loading image "
	for(i=0; i<framenumber; i+=1)
		if(!mod(i, 100))
			printf "%d . . .", i
		endif
		if(!mod(i, 1000))
			printf "\r"
		endif
		tiff_name = StringFromList(i, all_tiffs)
		tiff_path = dirname + tiff_name
		imageload/q/t=tiff/N=name tiff_path		
		wave temp=$"name"
		data[][][i]=temp[p][q]
		killwaves temp
	endfor
	
	print time()+"\r"
end
