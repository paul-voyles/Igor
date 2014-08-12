#pragma rtGlobals=1		// Use modern global access method.

// Utility functions for image data analysis
//
// PowerSpectrum2D		normalized power spectrum of an image in 2 dimensions
// PowerSpectrum1DInt	normalized power spectrum of an image with annular integration
// PowerSpectrum1DAv	normalized power spectrum of an image with annular averaging
// AnnularIntegrate		performs an annular integration about the center of an image
// AnnularAverage			performs an annular average in strips of a specified width about a specified point
// ImageDamp			applies a damping function to an image prior to calculating the power spectrum
// ImageDampDM			applies the Voyles/Treacy/Gibson exponential image damping function
// CountFiles				counts the files in a directory
// LoadTEM				load a TIFF image and massages the name to make Igor happy
// LEOMagCal			sets the image wave scaling for an image from the UW LEO at a given mag
// ClipImage				clips the input to the specified size.
//
// PMV and WGS, 6/25/03.
//
// renamed from "FEM Analysis Utilities.ipf", 8/12/14.  Not compatible with "STEM" group of procedures.

// Calculates the 2D power spectrum from an input image, with appropriate normalization
// for Parseval's theorem to work.  Output is in the wave im_ps2D.  Wave scaling is 
// preserved.
function PowerSpectrum2D(im)
	wave im
	
	//Copy image prior to FFT
	Duplicate/O im im_fft
	Duplicate/O im im_ps2d
	
	//Perform the FFT
	Redimension/s/c im_fft
	FFT im_fft
	
	// Calculate the normalized power spectrum.
	// Normalization to fit Parseval's theorem: variance = integral of power spectrum
	variable scale = (DimSize(im, 0)*DimSize(im, 1))^2
	im_ps2d=Real(im_fft[p][q]*conj(im_fft[p][q]))
	im_ps2d /= scale

	// set the units of the power spectrum appropriately.  (0, 0) is in the center of the image
	variable dm = 1 / ( DimDelta(im, 0) * DimSize(im, 0) )
	variable dm0 = -dm*DimSize(im, 0) / 2
	variable dn = 1 / (DimDelta(im, 1) * DimSize(im, 1) )
	variable dn0 = -dn*DimSize(im, 1)/2
	SetScale/P x dm0, dm, "", im_ps2D
	SetScale/P y dn0, dn, "", im_ps2D

	// Clean up the unneeded complex FFT wave
	Killwaves/Z im_fft

end

// 1D power spectrum of the input wave, computed as an annular average
// of the 2D power spectrum in strips of strip_width.  Value is returned in
// the wave im_ps1D, which is appropriately scaled.
function PowerSpectrum1DInt(im, strip_width)
	wave im	
	variable strip_width

	PowerSpectrum2D(im)
	wave im_ps2d = $"im_ps2D"
	if(!WaveExists(im_ps2D))
		printf "PowerSpectrum2D has failed somehow. Exiting.\r"
		return 0
	endif
	
	// PowerSpectrum2D scales im_ps2D so that (0, 0) in the right place,
	// so we can tell AnnularAverage to work about (0, 0)
	//AnnularAverage(im_ps2D, 0, 0, strip_width)
	AnnularIntegrate(im_ps2D, strip_width)
	wave annular_int = $"annular_int"
	if(!WaveExists(annular_int))
		printf "AnnularAverage has failed somehow.  Exiting.\r"
		return 0
	endif
	
	Duplicate/O annular_int im_ps1D
	
	Killwaves/Z annular_int, im_ps2D

end
	
// 1D power spectrum of the input wave, computed as an annular average
// of the 2D power spectrum in strips of strip_width.  Value is returned in
// the wave im_ps1D, which is appropriately scaled.
function PowerSpectrum1DAv(im, strip_width)
	wave im	
	variable strip_width

	PowerSpectrum2D(im)
	wave im_ps2d = $"im_ps2D"
	if(!WaveExists(im_ps2D))
		printf "PowerSpectrum2D has failed somehow. Exiting.\r"
		return 0
	endif
	
	// PowerSpectrum2D scales im_ps2D so that (0, 0) in the right place,
	// so we can tell AnnularAverage to work about (0, 0)
	//AnnularAverage(im_ps2D, 0, 0, strip_width)
	AnnularAverage(im_ps2D, 0, 0, strip_width)
	wave annular_int = $"annular_av"
	if(!WaveExists(annular_int))
		printf "AnnularAverage has failed somehow.  Exiting.\r"
		return 0
	endif
	
	Duplicate/O annular_int im_ps1D
	
	Killwaves/Z annular_int, im_ps2D

end
	

// Does a rotational average about the specified point, xx and yy with a bin width of strip_width.
// xx and yy are in real units, taking into account the wave scaling.  strip_width is in pixels.
// average is placed the wave annular_av, and the function returns the number of pixels in 
// the annular average.
function AnnularAverage(dat, xx, yy, strip_width)
	wave dat
	variable xx, yy, strip_width
	
	// translate center point from wave scaling units to pixels
	variable cp = round( ((xx - DimOffset(dat, 0))/DimDelta(dat, 0)) )
	variable cq = round( ((yy - DimOffSet(dat, 1))/DimDelta(dat, 1)) )
	
	// find the radius of the annular average accounting for strip_width
	variable sx = DimSize(dat, 0)
	variable sy = DimSize(dat, 1)
	variable rad = min( min(cp, abs(cp-sx)), min(cq, abs(cq-sy)) )
	rad = floor( rad / strip_width)
	
	//printf "The point closest to (%g, %g) is (%g, %g).\r", xx, yy, cp*DimDelta(dat, 0) + DimOffSet(dat, 0), cq*DimDelta(dat, 1) + DimOffSet(dat, 1)
	//printf "Annular average runs out to %g pixels\r", rad

	// make the output and pixel-counting waves
	Make/O/N=(rad) annular_av, npix
	SetScale/P x 0, DimDelta(dat, 0)*strip_width, "", annular_av
	annular_av = 0
	npix = 0
	
	// perform the annular average
	variable i, j, rpix
	for(i=0; i<sx; i+=1)
		for(j=0; j<sy; j+=1)
			rpix = round( sqrt( (i-cp)^2 + (j-cq)^2) / strip_width )
			if(rpix < rad)
				annular_av[rpix] += dat[i][j]
				npix[rpix] += 1
			endif
		endfor
	endfor
	
	// normalize by the number of pixels in each bin
	annular_av /= npix

	// clean up temporary waves
	Killwaves/Z npix

	// return the radius of the annular average
	return rad	
	
end

// annular integration only about the center of the image
function AnnularIntegrate(dat, strip_width)
	wave dat
	variable strip_width

	variable sx = DimSize(dat, 0)
	variable sy = DimSize(dat, 1)
	
	variable xc = floor(sx / 2)
	variable yc = floor(sy / 2)
	variable rad = floor( min(xc, yc) / strip_width)
	
	Make/O/N=(rad) annular_int
	SetScale/P x 0, DimDelta(dat, 0)*strip_width, "", annular_int
	annular_int = 0
	
	variable ii, jj, r
	for(ii=0; ii<sx; ii+=1)
		for(jj=0; jj<sy; jj+=1)
			r = floor( sqrt( (ii - xc)^2 + (jj-yc)^2 ) / strip_width)
			if(r < rad)
				annular_int[r] += dat[ii][jj]
			endif
		endfor
	endfor

end
	

// At this point just applies a Hamming window to the wave using the
// function ImageWindow.  Split out into a separate function mostly
// so we can change it later if we want.
// Overwrite the input wave with the damped version.
function ImageDamp(im)
	wave im
	
	ImageWindow/O hamming im
	return V_Value	// power spectrum correction value set by ImageWindow
end


function ImageDampDM(im)
	wave im
	
	variable rad, sigma, sx, sy, pc, qc
	sx = DimSize(im, 0)
	sy = DimSize(im, 1)
	pc = floor(sx/2)
	qc = floor(sy/2)
	
	rad = min(pc, qc) - 30
	sigma = 10
	
	variable correction = Pi*(rad+sigma)^2 / (DimSize(im, 0)*DimSize(im, 1))

	imagestats im
	im = (  ( (p-pc)^2+(q-qc)^2 > 1e4) ? ( im[p][q] * exp( -( (sqrt((p-pc)^2+(q-qc)^2) - 100)^2 / 200) ) + (1-exp( -(( sqrt((p-pc)^2+(q-qc)^2) - 100)^2 / 200)))*V_avg ) : im[p][q]) 
	//orig_clip = tert(iradius > 100 , orig_clip * exp(-((iradius-100)**2)/200) + (1-exp(-((iradius-100)**2)/200)) * mean, orig_clip)
end


// counts the number of files in a directory specified by pathname.
function CountFiles(pathname)
	string pathname
	
	string filename
	variable num_file = 0
	do
		filename=IndexedFile($pathname,num_file,"????")
		if (strlen(filename)==0)  //no more files in folder
			break
		endif		
		num_file+=1
	while (1)

	return num_file
end


// loads a TIFF format file and massages the wave name to make
// Igor happy: spaces and periods are replaced by "_" and ".tif" is stripped.
function/S loadTEM(filename, pathname)
	string filename, pathname
		
	//load the image
	Imageload/O/T=TIFF/P=$pathname filename 
	string wname = StringFromList(0, S_waveNames, ";")
	Redimension/S $wname
	
	//strip off ".tif" file extension
	string new_wname, working
	variable changed = 0
	working = wname
	variable pos = strsearch(working, ".tif", 0)
	if(pos > 0)	
		new_wname = working[0,pos-1]
		working = new_wname
		changed = 1
	endif
	
	new_wname = ""
	for(pos=0; pos<strlen(working); pos+=1)
		if(cmpstr(working[pos], ".") == 0)
			new_wname[pos] = "_"
			changed = 1
		else
			new_wname[pos] = working[pos]
		endif
	endfor
	
	if(changed)
		Duplicate/O $wname $new_wname
		Killwaves/Z $wname
		wname = new_wname
	endif

	return wname
end


function LoadTEMFolder()

	NewPath/O/M="Select a folder of .tif files" temp_path
	
	string filename
	variable index
	index=0
	do
		filename=IndexedFile(temp_path,index,"????")
		if (strlen(filename)==0)  //no more files in folder
			break
		endif
		
		wave im = $loadTEM(filename, "temp_path") //loads the TEMimage
		if(!WaveExists(im))
			printf "Image %s not loaded.  Exiting.\r", filename
			return 0
		endif
		index+=1
	while(1)
	
end

// Sets the wave scaling for a TEM image from the UW LEO 912
// angstrom / pixel calibration from Stratton, 2/19/03
function LeoMagCal(im, bin, mag)
	wave im
	variable bin, mag
	
	variable ang_per_pixel = (1.31e5 / mag)*bin
	
	SetScale/P x 0, ang_per_pixel, "", im
	SetScale/P y 0, ang_per_pixel, "", im
	
end


// clips the image im to the specified subimage
// the original image is overwritten.
function ClipImage(im, bot, left, top, right)
	wave im
	variable bot, left, top, right
	
	Make/O/N=( (top - bot), (right - left)) cliptmp
	SetScale/P x 0, DimDelta(im, 0), "", cliptmp
	SetScale/P y 0, DimDelta(im, 1), "", cliptmp
	
	cliptmp = im[p+left][q+bot]
	
	Duplicate/O cliptmp im
	Killwaves/Z cliptmp
end
