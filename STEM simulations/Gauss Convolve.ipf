#pragma rtGlobals=1		// Use modern global access method.

// This is code to include the source size of the incoherent electron emitter into STEM image simulations.
// I used Paul Voyles "Source Size Convlove" code as a starting point.
// Due to edge effects of the convolution process, the image needs to be embeded in a larger image before convolution.
// This code gives you the option to embed the original image in zero or in a non-zero value (= average of the outside ring of pixels). Non-zero is more accurate than zero.
// If your original image obeys periodic boundary conditions, you can also embed the image in a period array of the original image (this is the most accurate way).
// "EmbedZeroAndConvolveImage", "EmbedNonZeroAndConvolveImage", and "PeriodicAndConvolveImage" take a single image as the input.
// "EmbedandConvolveStack" takes a stack of images as the input.
// Before each procedure there are explanations of what that procedure does.
// Andrew B. Yankovich (6/24/14)
//
// Updated PeriodicContinue to preserve the wave scaling.  4/27/15, pmv



// takes an image im and crops it with boundary xi yi xf yf.  
// xi yi xf yf are all within the final cropped image.
function cropimage(im, xi, yi, xf, yf)
	wave im
	variable xi, yi, xf, yf

	DeletePoints xf+1,10000, im
	DeletePoints 0,xi, im
	DeletePoints/M=1 yf+1,10000, im
	DeletePoints/M=1 0,yi, im
end



// Takes an image (im) and convolves it with a Gaussian source function with FWHM (ds).  
// ds must be in the units of the pixel calibration.
// Uses the wave scaling of im.  
// Preserves the total intensity of the image (sum of all pixels) after scaling.
function SourceSizeConvolve(im, ds)
	wave im
	variable ds
	
	variable im_norm = sum(im)
	Duplicate/O im $(NameofWave(im)+"_ss")
	wave im_ss = $(NameofWave(im)+"_ss")
	variable fds = ds / (2*sqrt(2*ln(2)))	// real-space standard deviation for Gaussian with FWHM ds
	fds = 1/(2*Pi*fds)	// FT standard deviation
	//printf "fds = %g\r", fds
	Redimension/C im_ss
	wave/C res = $(NameofWave(im)+"_ss")
	FFT res
	res *=  cmplx(Gauss(x, 0.0, fds, y, 0.0, fds), 0.0)
	IFFT/C res
	Redimension/R res
	variable ss_norm = sum(res) 
	res *= im_norm / ss_norm
	SetScale/p x dimoffset(im, 0), dimdelta(im, 0), "", res
	SetScale/P y dimoffset(im, 1), dimdelta(im, 1), "", res
end



// This function embeds a single image (im) in the center of a larger image with magnification size (mag_x and mag_y) 
// with intensity value (value). Make the magnification an odd number.
function Embed(im, mag_x, mag_y, value)
	
	wave im
	variable mag_x, mag_y, value 

	variable xsize = DimSize(im,0)
	variable ysize = DimSize(im,1)
	variable xstart =((mag_x-1)/2)*xsize 
	variable ystart =((mag_y-1)/2)*ysize 
	Make/o/n=(mag_x*xsize, mag_y*ysize) $(NameofWave(im)+"_embed")
	
	wave/C new_image = $(NameofWave(im)+"_embed")
	
	new_image = value
	imagetransform /INSI=im /INSX=(xstart) /INSY=(ystart) insertImage new_image
end



// This function periodic continues a single image (im) with a magnification (mag_x and mag_y). Make mag an odd number.
function PeriodicContinue(im, mag_x, mag_y)
	
	wave im
	variable mag_x
	variable mag_y

	variable xsize = DimSize(im,0)
	variable ysize = DimSize(im,1)
	variable zsize = DimSize(im,2)
	Make/o/n=((mag_x*xsize), (mag_y*ysize), zsize) $(NameofWave(im)+"_periodic")
	
	wave/C new_image = $(NameofWave(im)+"_periodic")
	
	variable i,j
	for(i=0; i<mag_x; i+=1)
		for(j=0; j<mag_y; j+=1)
			imagetransform /INSI=im /INSX=(i*xsize) /INSY=(j*ysize) insertImage new_image
		endfor
	endfor
	
	SetScale/P x DimSize(im, 0), DimDelta(im, 0), new_image
	SetScale/P y DimSize(im, 1), DimDelta(im, 1), new_image
	
end



// This function embeds an image(im) in the center of a larger image with magnification size (mag_x and mag_y) 
// and with intensity value (value = 0). Make the magnification an odd number.
// It convolves that embeded image with a gaussian source size (ss) and preserves the intensity of the original image.
// ss must have the units of the pixel calibration of im.
function EmbedZeroAndConvolveImage(im, mag_x, mag_y, ss)
	
	wave im
	variable mag_x, mag_y, ss
	
	variable xscale = DimDelta(im, 0)
	variable yscale = DimDelta(im, 1)
	variable xsize = DimSize(im,0)
	variable ysize = DimSize(im,1)
	
	variable Im_Sum_old = sum(im)
	
	variable value = 0
	
	Embed(im, mag_x, mag_y, value)
	wave embed =$(NameofWave(im)+"_embed")
		
	setscale/p x, 0, xscale, embed
	setscale/p y, 0, yscale, embed

	SourceSizeConvolve(embed, ss)
	wave embed_ss =$(NameofWave(im)+"_embed_ss")
	
	Duplicate/o embed_ss im_ss
	
	cropimage(im_ss, ((mag_x-1)/2)*xsize, ((mag_y-1)/2)*ysize, (((mag_x-1)/2)*xsize)+xsize-1, (((mag_y-1)/2)*ysize)+ysize-1)
	
	variable Im_Sum_new = sum(im_ss)
	
	im_ss = im_ss * (Im_Sum_old/Im_Sum_new)
	
	setscale/p x, 0, xscale, im_ss
	setscale/p y, 0, yscale, im_ss
	
	Duplicate/o im_ss $(NameofWave(im)+"_ss")
	killwaves im_ss
end



// This function embeds an image(im) in the center of a larger image with magnification size (mag_x and mag_y) 
// and with intensity value = average of outer rim of im. Make the magnification an odd number.
// It convolves that embeded image with a gaussian source size (ss) and preserves the intensity of the original image.
// ss must have the units of the pixel calibration of im.
function EmbedNonZeroAndConvolveImage(im, mag_x, mag_y, ss)
	
	wave im
	variable mag_x, mag_y, ss

	variable xscale = DimDelta(im, 0)
	variable yscale = DimDelta(im, 1)
	variable xsize = DimSize(im,0)
	variable ysize = DimSize(im,1)
	
	variable Im_Sum_old = sum(im)
	
	Duplicate/o im LeftCol
	cropimage(LeftCol, 0, 0, xsize-1, 0)
	Duplicate/o im BotRow
	cropimage(BotRow, xsize-1, 1, xsize-1, ysize-2)
	Duplicate/o im TopRow
	cropimage(TopRow, 0, 1, 0, ysize-2)
	Duplicate/o im RightCol
	cropimage(RightCol, 0, ysize-1, xsize-1, ysize-1)	
	variable value = (sum(TopRow)+sum(BotRow)+sum(LeftCol)+sum(RightCol))/((2*xsize)+(2*ysize)-4)
	
	Embed(im, mag_x, mag_y, value)
	wave embed =$(NameofWave(im)+"_embed")
		
	setscale/p x, 0, xscale, embed
	setscale/p y, 0, yscale, embed

	SourceSizeConvolve(embed, ss)
	wave embed_ss =$(NameofWave(im)+"_embed_ss")
	
	Duplicate/o embed_ss im_ss
	
	cropimage(im_ss, ((mag_x-1)/2)*xsize, ((mag_y-1)/2)*ysize, (((mag_x-1)/2)*xsize)+xsize-1, (((mag_y-1)/2)*ysize)+ysize-1)
	
	variable Im_Sum_new = sum(im_ss)
		
	im_ss = im_ss * (Im_Sum_old/Im_Sum_new)
	
	setscale/p x, 0, xscale, im_ss
	setscale/p y, 0, yscale, im_ss
	
	Duplicate/o im_ss $(NameofWave(im)+"_ss")
	killwaves im_ss, TopRow, BotRow, LeftCol, RightCol
end



// This function periodic continues an image(im) in the center of a larger image with magnification size = mag_x and mag_y. 
// Make mag an odd number.
// It convolves that periodic image with a gaussian source size (ss) and preserves the intensity of the original image.
// ss must have the units of the pixel calibration of im.
function PeriodicAndConvolveImage(im, mag_x, mag_y, ss)
	
	wave im
	variable mag_x, mag_y, ss

	variable xscale = DimDelta(im, 0)
	variable yscale = DimDelta(im, 1)
	variable xsize = DimSize(im,0)
	variable ysize = DimSize(im,1)
	
	variable Im_Sum_old = sum(im)
	//print "old sum zero embed = ", Im_Sum_old
	
	PeriodicContinue(im, mag_x, mag_y)
	wave periodic =$(NameofWave(im)+"_periodic")
		
	setscale/p x, 0, xscale, periodic
	setscale/p y, 0, yscale, periodic

	SourceSizeConvolve(periodic, ss)
	wave periodic_ss =$(NameofWave(im)+"_periodic_ss")
	
	Duplicate/o periodic_ss im_ss
	
	cropimage(im_ss, (floor(mag_x/2)*(xsize)), (floor(mag_y/2)*(ysize)), (((floor(mag_x/2))+1)*xsize)-1, (((floor(mag_y/2))+1)*ysize)-1)
	
	variable Im_Sum_new = sum(im_ss)
	//print "new sum zero embed = ", Im_Sum_new
		
	im_ss = im_ss * (Im_Sum_old/Im_Sum_new)
	
	setscale/p x, 0, xscale, im_ss
	setscale/p y, 0, yscale, im_ss
	
	Duplicate/o im_ss $(NameofWave(im)+"_ss")
	killwaves im_ss
end



// This function embeds a stack of images(im_stack) in the center of a larger stack of images with magnification size (mag_x and mag_y).
// Make mag an odd number.
// If you want intensity value to be 0, make type = 0.
// If you want intensity value to be average of outside rim of image, make type = 1.
// If you want to periodic continue your image before convolving with the ss, make type = 2.
// It convolves each layer in the stack with a gaussian source size (ss), and preserves the intensity of the original image.
// ss must have the units of the pixel calibration of im.
function EmbedandConvolveStack(im_stack, mag_x, mag_y, type, ss)
	
	wave im_stack
	variable mag_x, mag_y, type, ss
	
	variable xscale = DimDelta(im_stack, 0)
	variable yscale = DimDelta(im_stack, 1)
	variable xsize = DimSize(im_stack,0)
	variable ysize = DimSize(im_stack,1)
	variable zsize = DimSize(im_stack,2)
	
	if (type != 0 && type != 1 && type != 2)
		print "type must be 0,1 or 2" 
	endif
	
	if (type == 0)
		print "Embeded in zero background" 
		Make/o/n=(mag_x*xsize, mag_y*ysize, zsize) $(NameofWave(im_stack)+"_embed")
		Make/o/n=(mag_x*xsize, mag_y*ysize, zsize) $(NameofWave(im_stack)+"_embed_ss")
		Make/o/n=(xsize, ysize, zsize) $(NameofWave(im_stack)+"_ss")
	endif
	
	if (type == 1)
		print "Embeded in non-zero background" 
		Make/o/n=(mag_x*xsize, mag_y*ysize, zsize) $(NameofWave(im_stack)+"_embed")
		Make/o/n=(mag_x*xsize, mag_y*ysize, zsize) $(NameofWave(im_stack)+"_embed_ss")
		Make/o/n=(xsize, ysize, zsize) $(NameofWave(im_stack)+"_ss")
	endif
	
	if (type == 2)
		print "Using a periodic continued image" 
		Make/o/n=(mag_x*xsize, mag_y*ysize, zsize) $(NameofWave(im_stack)+"_periodic")
		Make/o/n=(mag_x*xsize, mag_y*ysize, zsize) $(NameofWave(im_stack)+"_periodic_ss")
		Make/o/n=(xsize, ysize, zsize) $(NameofWave(im_stack)+"_ss")
	endif
	

	
	variable i
	for(i=0; i<zsize; i+=1)
		imagetransform/p=(i) getplane im_stack
		wave M_ImagePlane =$"M_ImagePlane"
		setscale/p x, 0, xscale, M_ImagePlane
		setscale/p y, 0, yscale, M_ImagePlane
		
		if (type == 0)
			EmbedZeroAndConvolveImage(M_ImagePlane, mag_x, mag_y, ss)
			
			wave im_embed =$"M_ImagePlane_embed"
			wave im_embed_ss =$"M_ImagePlane_embed_ss"
			wave im_ss =$"M_ImagePlane_ss"
		
			imagetransform /INSI=im_embed /INSX=0 /INSY=0 /P=(i) insertImage $(NameofWave(im_stack)+"_embed")
			imagetransform /INSI=im_embed_ss /INSX=0 /INSY=0 /P=(i) insertImage $(NameofWave(im_stack)+"_embed_ss")
			imagetransform /INSI=im_ss /INSX=0 /INSY=0 /P=(i) insertImage $(NameofWave(im_stack)+"_ss")
		
			killwaves M_ImagePlane, M_ImagePlane_embed, M_ImagePlane_embed_ss, M_ImagePlane_ss
		endif
		
		if (type == 1)
			EmbedNonZeroAndConvolveImage(M_ImagePlane, mag_x, mag_y, ss)
			
			wave im_embed =$"M_ImagePlane_embed"
			wave im_embed_ss =$"M_ImagePlane_embed_ss"
			wave im_ss =$"M_ImagePlane_ss"
		
			imagetransform /INSI=im_embed /INSX=0 /INSY=0 /P=(i) insertImage $(NameofWave(im_stack)+"_embed")
			imagetransform /INSI=im_embed_ss /INSX=0 /INSY=0 /P=(i) insertImage $(NameofWave(im_stack)+"_embed_ss")
			imagetransform /INSI=im_ss /INSX=0 /INSY=0 /P=(i) insertImage $(NameofWave(im_stack)+"_ss")
		
			killwaves M_ImagePlane, M_ImagePlane_embed, M_ImagePlane_embed_ss, M_ImagePlane_ss
		endif
		if (type == 2)
			PeriodicAndConvolveImage(M_ImagePlane, mag_x, mag_y, ss)
			
			wave im_periodic =$"M_ImagePlane_periodic"
			wave im_periodic_ss =$"M_ImagePlane_periodic_ss"
			wave im_ss =$"M_ImagePlane_ss"
		
			imagetransform /INSI=im_periodic /INSX=0 /INSY=0 /P=(i) insertImage $(NameofWave(im_stack)+"_periodic")
			imagetransform /INSI=im_periodic_ss /INSX=0 /INSY=0 /P=(i) insertImage $(NameofWave(im_stack)+"_periodic_ss")
			imagetransform /INSI=im_ss /INSX=0 /INSY=0 /P=(i) insertImage $(NameofWave(im_stack)+"_ss")
		
			killwaves M_ImagePlane, M_ImagePlane_periodic, M_ImagePlane_periodic_ss, M_ImagePlane_ss
		endif
		
	endfor
end