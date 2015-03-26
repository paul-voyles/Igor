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
function PeriodicContinueStack(im_stack, mag_x, mag_y)
	
	wave im_stack
	variable mag_x
	variable mag_y

	variable xscale = DimDelta(im_stack, 0)
	variable yscale = DimDelta(im_stack, 1)

	variable xsize = DimSize(im_stack,0)
	variable ysize = DimSize(im_stack,1)
	variable zsize = DimSize(im_stack,2)
	Make/o/n=((mag_x*xsize), (mag_y*ysize), zsize) $(NameofWave(im_stack)+"_periodic")
	
	wave/C new_image = $(NameofWave(im_stack)+"_periodic")
	
	variable k
	for(k=0; k<zsize; k+=1)
		imagetransform/p=(k) getplane im_stack
		wave M_ImagePlane =$"M_ImagePlane"
	
		variable i,j
		for(i=0; i<mag_x; i+=1)
			for(j=0; j<mag_y; j+=1)
				imagetransform /INSI=M_ImagePlane /INSX=(i*xsize) /INSY=(j*ysize)/P=(k) insertimage new_image
			endfor
		endfor
	endfor
	
	setscale/p x, 0, xscale, new_image
	setscale/p y, 0, yscale, new_image
	
	killwaves M_ImagePlane
end


// This function periodic continues a single image (im) with a magnification (mag_x and mag_y). Make mag an odd number.
function PeriodicContinue(im, mag_x, mag_y)
	
	wave im
	variable mag_x
	variable mag_y

	variable xsize = DimSize(im,0)
	variable ysize = DimSize(im,1)
	Make/o/n=((mag_x*xsize), (mag_y*ysize)) $(NameofWave(im)+"_periodic")
	
	wave/C new_image = $(NameofWave(im)+"_periodic")
	
	variable i,j
	for(i=0; i<mag_x; i+=1)
		for(j=0; j<mag_y; j+=1)
			imagetransform /INSI=im /INSX=(i*xsize) /INSY=(j*ysize) insertImage new_image
		endfor
	endfor
	
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




function VariableSourceSizePeriodic(im, mag_x, mag_y, ss_i, ss_f, ss_num)
	
	wave im
	variable mag_x, mag_y, ss_i, ss_f, ss_num
	variable ss
	
	make/o/n=(dimsize(im,0), dimsize(im,1), ss_num) SS_stack
	
	variable i
	for(i=0; i<ss_num; i+=1)
		ss = ss_i + (i*(ss_f-ss_i)/ss_num)
		PeriodicAndConvolveImage(im, mag_x, mag_y, ss)
		wave im_ss = $(NameofWave(im)+"_ss")
		imagetransform/p=(i)/INSW=im_ss insertimage SS_stack
	endfor
	
	killwaves $(NameofWave(im)+"_periodic"), $(NameofWave(im)+"_periodic_ss")
end



function SourceSizeDetermineWidth(sim_im, exp_im, sim_x, sim_y, sim_size, exp_x, exp_y, exp_size, mag_x, mag_y, ss_i, ss_f, ss_num)

	wave sim_im, exp_im
	variable sim_x, sim_y, sim_size, exp_x, exp_y, exp_size, mag_x, mag_y, ss_i, ss_f, ss_num
	
	VariableSourceSizePeriodic(sim_im, mag_x, mag_y, ss_i, ss_f, ss_num)
	wave sim_im_stack = $"SS_stack"
	setscale/p x, 0, dimdelta(sim_im, 0), sim_im_stack
	setscale/p y, 0, dimdelta(sim_im, 1), sim_im_stack
	
	variable xe_i, ye_i, xe_f, ye_f
	xe_i = exp_x - (exp_size/2)
	ye_i = exp_y - (exp_size/2)
	xe_f = exp_x + (exp_size/2)
	ye_f = exp_y + (exp_size/2)
	
	CurveFit/M=2/W=2/Q gauss2D, exp_im[xe_i,xe_f][ye_i,ye_f]
	
	wave W_coef =$"W_coef"
	Duplicate/o W_coef Exp_coef
	wave W_sigma =$"W_sigma"
	Duplicate/o W_sigma Exp_sigma
	variable exp_width = (W_coef[3] + W_coef[5]) /2
	//print exp_width

	variable xs_i, ys_i, xs_f, ys_f
	xs_i = sim_x - (sim_size/2)
	ys_i = sim_y - (sim_size/2)
	xs_f = sim_x + (sim_size/2)
	ys_f = sim_y + (sim_size/2)
	
	make/o/n=(dimsize(sim_im_stack,2)) width_diff
	
	variable sim_width
	variable diff, diff_min, ss_min
	diff_min = 10000
 	variable i=0
	for(i=0; i<dimsize(sim_im_stack,2); i+=1)
		imagetransform/p=(i) getplane sim_im_stack
		wave sim_tmp =$"M_ImagePlane"
		setscale/p x, 0, dimdelta(sim_im, 0), sim_tmp
		setscale/p y, 0, dimdelta(sim_im, 1), sim_tmp
		
		CurveFit/M=2/W=2/Q gauss2D, sim_tmp[xs_i,xs_f][ys_i,ys_f]
		wave W_coef =$"W_coef"
		sim_width = (W_coef[3] + W_coef[5]) /2
		//print sim_width
		diff = abs(sim_width - exp_width)
		width_diff[i] = diff
		if(diff<diff_min)
			diff_min = diff
			ss_min = ss_i + (i*(ss_f-ss_i)/ss_num)
		endif
	endfor

	print "The source size is ", ss_min
	killwaves W_sigma, W_coef, M_ImagePlane, M_Covar
end

// assume exp pacbed patterns were acquired with same expose time
function PACBEDChiCompare(im_e_s, im_e_ns, im_s)
	wave im_e_s, im_e_ns, im_s
	
	variable x_pos, y_pos, radius
	
	// crop exp nosample image
	ImageThreshold/q/m=1/I im_e_ns
	ImageAnalyzeParticles/q/m=3/a=50 stats M_ImageThresh
	wave W_xmin = $"W_xmin"
	wave W_xmax = $"W_xmax"
	wave W_ymin = $"W_ymin"
	wave W_ymax = $"W_ymax"
	x_pos = (W_xmax[0] + W_xmin[0])/2
	y_pos = (W_ymax[0] + W_ymin[0])/2
	radius = (((W_ymax[0] - W_ymin[0])/2) + ((W_xmax[0] - W_xmin[0])/2))/2
	duplicate/o im_e_ns im_e_ns_c
	cropimage(im_e_ns_c, x_pos - 2*radius, y_pos - 2*radius, x_pos + 2*radius, y_pos + 2*radius)
	
	// crop exp sample image
	//ImageThreshold/q/m=1/I im_exp_sample
	//ImageAnalyzeParticles/q/m=3/a=50 stats M_ImageThresh
	//wave W_xmin = $"W_xmin"
	//wave W_xmax = $"W_xmax"
	//wave W_ymin = $"W_ymin"
	//wave W_ymax = $"W_ymax"
	//x_pos = (W_xmax[0] + W_xmin[0])/2
	//y_pos = (W_ymax[0] + W_ymin[0])/2
	//radius = (((W_ymax[0] - W_ymin[0])/2) + ((W_xmax[0] - W_xmin[0])/2))/2
	duplicate/o im_e_s im_e_s_c
	cropimage(im_e_s_c, x_pos - 2*radius, y_pos - 2*radius, x_pos + 2*radius, y_pos + 2*radius)
	
	// crop simulated image stack
	imagetransform/p=0 getplane im_s
	wave M_ImagePlane = $"M_ImagePlane"
	ImageThreshold/q/m=1/I M_ImagePlane
	ImageAnalyzeParticles/q/m=3/a=50 stats M_ImageThresh
	wave W_xmin = $"W_xmin"
	wave W_xmax = $"W_xmax"
	wave W_ymin = $"W_ymin"
	wave W_ymax = $"W_ymax"
	x_pos = (W_xmax[0] + W_xmin[0])/2
	y_pos = (W_ymax[0] + W_ymin[0])/2
	radius = (((W_ymax[0] - W_ymin[0])/2) + ((W_xmax[0] - W_xmin[0])/2))/2
	duplicate/o im_s im_s_c
	cropimage(im_s_c, x_pos - 2*radius, y_pos - 2*radius, x_pos + 2*radius, y_pos + 2*radius)
	
	// resample experimental PACBED patterns to match simulated PACBED sampling
	duplicate/o im_e_ns_c im_e_ns_c_r
	Resample/UP=(dimsize( im_s_c,0)+1)/DOWN=(dimsize( im_e_ns_c,0)+1)/DIM=0  im_e_ns_c_r
	Resample/UP=(dimsize( im_s_c,0)+1)/DOWN=(dimsize( im_e_ns_c,0)+1)/DIM=1  im_e_ns_c_r
	duplicate/o im_e_s_c im_e_s_c_r
	Resample/UP=(dimsize( im_s_c,0)+1)/DOWN=(dimsize( im_e_ns_c,0)+1)/DIM=0  im_e_s_c_r
	Resample/UP=(dimsize( im_s_c,0)+1)/DOWN=(dimsize( im_e_ns_c,0)+1)/DIM=1  im_e_s_c_r
	
	//scale exp PACBED patterns
	imagetransform/p=0 getplane im_s_c
	wave M_ImagePlane = $"M_ImagePlane"
	Wavestats/q M_ImagePlane
	variable sim_sum = V_sum
	Wavestats/q im_e_ns_c_r
	variable exp_ns_sum = V_sum
	duplicate/o im_e_ns_c_r im_e_ns_n
	im_e_ns_n = im_e_ns_c_r *(sim_sum/exp_ns_sum)
	duplicate/o im_e_s_c_r im_e_s_n
	im_e_s_n = im_e_s_c_r *(sim_sum/exp_ns_sum)
	
	duplicate/o im_e_s_n im_chi
	duplicate/o im_s_c im_chi_stack
	im_chi_stack = 0
	variable chi_min = 10000000 
	variable chi_tmp, pos_min
	Make/o/n=(dimsize(im_s_c,2)) chi
	variable i
	for (i=0; i<dimsize(im_s_c,2); i+=1)
		imagetransform/p=(i) getplane im_s_c
		wave M_ImagePlane = $"M_ImagePlane"
		im_chi = ((im_e_s_n - M_ImagePlane)^2)/M_ImagePlane
		ImageTransform/o/p=(i) removeZplane im_chi_stack  
		ImageTransform/o/p=(i)/insw=im_chi insertZplane im_chi_stack  
		wavestats/q im_chi
		chi_tmp = V_sum
		chi[i] = chi_tmp
		if (chi_tmp < chi_min)
			chi_min = chi_tmp
			pos_min = i
		endif
	endfor
	
	killwaves M_ImagePlane, W_ImageObjArea, W_SpotX, W_SpotY, W_circularity, W_rectangularity, W_ImageObjPerimeter, W_xmin, W_xmax, W_ymin, W_ymax, M_RawMoments
end