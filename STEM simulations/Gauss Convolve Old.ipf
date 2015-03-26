#pragma rtGlobals=1		// Use modern global access method.

// Function to convolve an image with a two-dimensional Guassian
// v1 04-20-12 pmv

// takes an image im and convolves it with a Gaussian source
// function with FWHM ds.  Uses the wave scaling of im.  Preserves
// the total intensity of the image (sum of all pixels) after scaling.
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
	res *= ss_norm / im_norm
	SetScale/p x dimoffset(im, 0), dimdelta(im, 0), "", res
	SetScale/P y dimoffset(im, 1), dimdelta(im, 1), "", res

	// added this to retain total intensity before and after convolution (Andy Yankovich)
	variable sum_old, sum_new
	sum_old = sum(im)
	sum_new = sum(res)
	res = res * (sum_old / sum_new)
	

end


// takes an image im and crops it with boundary xi yi xf yf.  
//xi yi xf yf are all within the final cropped image.

function cropimage(im, xi, yi, xf, yf)

wave im
variable xi, yi, xf, yf

DeletePoints xf+1,10000, im
DeletePoints 0,xi, im
DeletePoints/M=1 yf+1,10000, im
DeletePoints/M=1 0,yi, im

end



// uses gaussian fit widths as a comparison
// takes simulated image and experimental image and calculate the microscope source size.
// inputs:
// sim_im is the simulated image embedded in a 0 background.
// exp_im is the experimental image to compare to, cropped, and averaged if desired.
// xi, yi, xf, yf are the cropping pixel boundaries to give same area as the experimental crop with the centroid nearest the middle pixel.
// ds_i is the source size start, ds_f is the source size finish, and ds_del is the delta ds.

function SourceSizeDetermineWidth(sim_im, exp_im, xi, yi, xf, yf, ds_i, ds_f, ds_del)


	wave sim_im, exp_im
	variable xi, yi, xf, yf, ds_i, ds_f, ds_del

	variable n_ds
	n_ds = 1 + (ds_f- ds_i) / ds_del
	print "number of source size calculations = ", n_ds
	variable ds
	make/o/n=(n_ds) xdiff, ydiff

	duplicate/o exp_im exp_fit, exp_res
	CurveFit/N/Q/NTHR=0 Gauss2D  exp_im /D=exp_fit /R=exp_res /A=0 

	wave W_coef =$"W_coef"
	Duplicate/o W_coef Exp_coef
	wave W_sigma =$"W_sigma"
	Duplicate/o W_sigma Exp_sigma

	killwaves W_coef, W_sigma


	variable source_size, xdiff_old, ydiff_old, xdiff_new, ydiff_new, av_old, av_new
 
	SourceSizeConvolve(sim_im, ds_i)
	wave im_ss = $(NameofWave(sim_im)+"_ss")
	cropimage(im_ss, xi, yi, xf, yf)
	duplicate/o im_ss sim_fit, sim_res, sim_fit_keep, sim_res_keep
	CurveFit/N/Q/NTHR=0 Gauss2D  im_ss /D=sim_fit /R=sim_res /A=0 
	wave W_coef =$"W_coef"
	Duplicate/o W_coef Sim_coef_keep
	wave W_sigma =$"W_sigma"
	Duplicate/o W_sigma Sim_sigma_keep

	xdiff_old = abs(W_coef[3] - Exp_coef[3])
	ydiff_old = abs(W_coef[5] - Exp_coef[5])
	av_old = (xdiff_old + ydiff_old)/2
	xdiff[0] = xdiff_old
	ydiff[0] = ydiff_old

	source_size =  ds_i
	sim_fit_keep = sim_fit
	sim_res_keep = sim_res

	killwaves W_coef, W_sigma, sim_fit, sim_res
	

	variable i=0
	for(i=1; i<n_ds; i+=1)
		ds = ds_i + (i * ds_del)
		//print ds
		SourceSizeConvolve(sim_im, ds)
		wave im_ss = $(NameofWave(sim_im)+"_ss")
		cropimage(im_ss, xi, yi, xf, yf)
		duplicate/o im_ss sim_fit, sim_res
		CurveFit/N/Q/NTHR=0 Gauss2D  im_ss /D=sim_fit /R=sim_res /A=0 
		wave W_coef =$"W_coef"
		wave W_sigma =$"W_sigma"
	
		xdiff_new = abs(W_coef[3] - Exp_coef[3])
		ydiff_new = abs(W_coef[5] - Exp_coef[5])
		av_new = (xdiff_new + ydiff_new)/2
	
		xdiff[i] = xdiff_new
		ydiff[i] = ydiff_new
	
		if(av_new < av_old)
			xdiff_old = xdiff_new
			ydiff_old = ydiff_new
			av_old = av_new
			Sim_coef_keep = W_coef
			Sim_sigma_keep = W_sigma
			source_size = ds
			sim_fit_keep = sim_fit
			sim_res_keep = sim_res
		endif
	
		killwaves W_coef, W_sigma, sim_fit, sim_res
	endfor

	print "The source size is ", source_size

end

// uses chi squared as a comparison
// takes simulated image and experimental image and calculate the microscope source size.
// inputs:
// sim_im is the simulated image embedded in a 0 background.
// exp_im is the experimental image to compare to, cropped, and averaged if desired. Also needs to be same size and pixel size as the simulated cropped area.
// xi, yi, xf, yf are the cropping pixel boundaries to give same area and pixel size as the experimental crop with the centroid nearest the middle pixel.
// ds_i is the source size start, ds_f is the source size finish, and ds_del is the delta ds.



function SourceSizeDetermineChi(sim_im, exp_im, xi, yi, xf, yf, ds_i, ds_f, ds_del)


	wave sim_im, exp_im
	variable xi, yi, xf, yf, ds_i, ds_f, ds_del

	variable n_ds
	n_ds = 1 + (ds_f- ds_i) / ds_del
	print "number of source size calculations = ", n_ds
	variable ds
	make/o/n=(n_ds) chi_all

	variable source_size, chi, chi_new
 
	SourceSizeConvolve(sim_im, ds_i)
	wave im_ss = $(NameofWave(sim_im)+"_ss")
	cropimage(im_ss, xi, yi, xf, yf)
	duplicate/o im_ss sq_diff, sim_keep
	sq_diff = (im_ss - exp_im)^2
	chi = sum(sq_diff)
	chi_all[0] = chi

	source_size =  ds_i

	variable i=0
	for(i=1; i<n_ds; i+=1)
		ds = ds_i + (i * ds_del)
		//print ds
		SourceSizeConvolve(sim_im, ds)
		wave im_ss = $(NameofWave(sim_im)+"_ss")
		cropimage(im_ss, xi, yi, xf, yf)
	
		sq_diff = (im_ss - exp_im)^2
		chi_new = sum(sq_diff)
		chi_all[i] = chi_new
	
		if(chi_new < chi)
			chi = chi_new
			sim_keep = im_ss
			source_size = ds
		endif
	
	endfor

	killwaves sq_diff, GaN_t34_embed_ss
	print "The source size is ", source_size

end


// This function embeds a single image(im) in the center of a larger image with magnification size (mag) 
//with intensity value (value). Make the magnification an odd number.
function Embed(im, mag, value)
	
	wave im
	variable mag, value 

	variable xsize = DimSize(im,0)
	variable ysize = DimSize(im,1)
	variable xstart =((mag-1)/2)*xsize 
	variable ystart =((mag-1)/2)*ysize 
	Make/o/n=(mag*xsize, mag*ysize) $(NameofWave(im)+"_embed")
	
	wave/C new_image = $(NameofWave(im)+"_embed")
	
	new_image = value
	imagetransform /INSI=im /INSX=(xstart) /INSY=(ystart) insertImage new_image
end


// This function embeds a stack of images(im) in the center of a larger stack of images
// with magnification size (mag) with intensity value (value). Make the magnification an odd number.
// It also convolves each layer in the stack with a gaussian source size (ss)
// px_calib is the pixel calibration of the simulated images in the units you want your source size to be in.
function EmbedandConvolveStack(im, mag, value, ss, px_calib)
	
	wave im
	variable mag, value, ss, px_calib

	variable xsize = DimSize(im,0)
	variable ysize = DimSize(im,1)
	variable zsize = DimSize(im,2)
	Make/o/n=(mag*xsize, mag*ysize, zsize) $(NameofWave(im)+"_embed_ss")
	
	wave/C new = $(NameofWave(im)+"_embed_ss")
	
	variable i
	for(i=0; i<zsize; i+=1)
		imagetransform/p=(i) getplane im
		wave M_ImagePlane =$"M_ImagePlane"
		
		Embed(M_ImagePlane, mag, value)
		wave M_ImagePlane_embed =$"M_ImagePlane_embed"
		
		setscale/p x, 0, px_calib, M_ImagePlane_embed
		setscale/p y, 0, px_calib, M_ImagePlane_embed
		
		SourceSizeConvolve(M_ImagePlane_embed, ss)
		wave im_ss =$"M_ImagePlane_embed_ss"
		
		imagetransform /INSI=im_ss /INSX=0 /INSY=0 /P=(i) insertImage new
		
		killwaves M_ImagePlane, M_ImagePlane_embed, M_ImagePlane_embed_ss
		
	endfor
	setscale/p x, 0, px_calib, new
	setscale/p y, 0, px_calib, new
end


// 
function StackAtomIntensity(im)
	
	wave im
	variable zsize = DimSize(im,2)
	
	Make/o/n=(zsize) IntensityVSslice
	
	variable i
	for(i=0; i<zsize; i+=1)
		
		Imagestats/P=(i) im
		IntensityVSslice[i] = V_avg
		
	endfor
end
