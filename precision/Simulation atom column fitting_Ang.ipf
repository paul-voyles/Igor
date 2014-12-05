#pragma rtGlobals=1		// Use modern global access method.


// takes an image im and crops it with boundary xi yi xf yf.  
// xi yi xf yf are all within the final cropped image.
function cropimage(im, xi, yi, xf, yf)

wave im
variable xi, yi, xf, yf

variable scalex = DimDelta(im, 0)
variable scaley = DimDelta(im, 1)

variable offsetx = DimOffset(im, 0)
variable offsety = DimOffset(im, 1)

variable xi_pixel = ((xi-offsetx)/scalex)
variable xf_pixel = ((xf-offsetx)/scalex) + 1
variable yi_pixel = ((yi-offsety)/scaley)
variable yf_pixel = ((yf-offsety)/scaley) + 1

DeletePoints xf_pixel,10000, im
DeletePoints 0,xi_pixel, im
DeletePoints/M=1 yf_pixel,10000, im
DeletePoints/M=1 0,yi_pixel, im

end



// This function fits the single peak positions that are in x_loc_s and y_loc_s. 
// x_loc_s and y_loc_s must be positions in Ang of the guessed atom position pixel
// This only works for 2D Gaussian fits. 
// Dumbbells need a different routine.
// Inputs: image, x_loc_s, y_loc_s, and gaussian fit size x size.
function GaussianFit(image, x_loc_s, y_loc_s, size_s)
	
	wave image, x_loc_s, y_loc_s
	variable size_s
	variable num_peaks_s = DimSize(x_loc_s,0)
	
	variable scalex = DimDelta(image, 0)
	variable scaley = DimDelta(image, 1)
	
	duplicate/O x_loc_s x_start
	duplicate/O x_loc_s x_finish
	duplicate/O y_loc_s y_start
	duplicate/O y_loc_s y_finish
	
	variable size_x = (round((size_s/2)/scalex))*scalex
	variable size_y = (round((size_s/2)/scaley))*scaley
	
	x_start = x_loc_s - size_x
	x_finish = x_loc_s + size_x - scalex
	y_start = y_loc_s - size_y
	y_finish = y_loc_s + size_y - scaley

	Make/D/O/N=(num_peaks_s) z0_s, A_s, x0_s, xW_s, y0_s, yW_s, cor_s
	Make/D/O/N=(num_peaks_s) sigma_z0_s, sigma_A_s, sigma_x0_s, sigma_xW_s, sigma_y0_s, sigma_yW_s, sigma_cor_s
	
	variable i=0
	for(i=0; i<num_peaks_s; i+=1)
		//crop out, save, and setscale of area to be fitted
		Duplicate/O image im
		cropimage(im, x_start[i], y_start[i], x_finish[i], y_finish[i])
		if( i == 0 )
			Make/O/D/N=(Dimsize(im,0), Dimsize(im,1), (num_peaks_s)) im_s
		endif
		imagetransform /INSI=im /INSX=0 /INSY=0 /P=(i) insertImage im_s
		setscale/P x, x_start[i], scalex, im
		setscale/P y, y_start[i], scaley, im
		
		Make/O/N=(DimSize(im,0), DimSize(im,1)) residual
		
		Duplicate/O/D im noise
		noise = sqrt(im)
		
		CurveFit/N/NTHR=0/W=2/Q gauss2D, im/R=residual/W=noise/D
		
		wave W_coef = $"W_coef"
		wave W_sigma = $"W_sigma"
		
		z0_s[i] = W_coef[0]
		A_s[i] = W_coef[1]
		x0_s[i] = W_coef[2]
		xW_s[i] = W_coef[3]
		y0_s[i] = W_coef[4]
		yW_s[i] = W_coef[5]
		cor_s[i] = W_coef[6]
		sigma_z0_s[i] = W_sigma[0]
		sigma_A_s[i] = W_sigma[1]
		sigma_x0_s[i] = W_sigma[2]
		sigma_xW_s[i] = W_sigma[3]
		sigma_y0_s[i] = W_sigma[4]
		sigma_yW_s[i] = W_sigma[5]
		sigma_cor_s[i] = W_sigma[6]	
		
		//save the fit images
		wave fit_im = $"fit_im"
		if( i == 0)
			Make/O/D/N=(Dimsize(fit_im,0), Dimsize(fit_im,1), (num_peaks_s)) fits_s
		endif	
		imagetransform /INSI=fit_im /INSX=0 /INSY=0 /P=(i) insertImage fits_s
		killwaves  fit_im
		
		//save the residual images
		if (i==0)
			Make/O/N=(DimSize(im,0), DimSize(im,1), num_peaks_s) residual_s
		endif	
		imagetransform /INSI=residual /INSX=0 /INSY=0 /P=(i) insertImage residual_s
		killwaves residual
		
		killwaves im, noise
		killwaves W_coef, W_sigma
	endfor	

killwaves x_start, x_finish, y_start, y_finish
end



//Creates the dumbbell 2 x 2DGaussian function.
Function Gaus2Dx2(w,x,y) : FitFunc
	Wave w
	Variable x
	Variable y

	//CurveFitDialog/ These comments were created by the Curve Fitting dialog. Altering them will
	//CurveFitDialog/ make the function less convenient to work with in the Curve Fitting dialog.
	//CurveFitDialog/ Equation:
	//CurveFitDialog/ f(x,y) = z0+(A*exp((-1/(2*(1-(cor1)^2)))*((((x-x1)/xw1)^2)+(((y-y1)/yw1)^2)-((2*cor1*(x-x1)*(y-y1))/(xw1*yw1)))))+ (B*exp((-1/(2*(1-(cor2)^2)))*((((x-x2)/xw2)^2)+(((y-y2)/yw2)^2)-((2*cor2*(x-x2)*(y-y2))/(xw2*yw2)))))
	//CurveFitDialog/ End of Equation
	//CurveFitDialog/ Independent Variables 2
	//CurveFitDialog/ x
	//CurveFitDialog/ y
	//CurveFitDialog/ Coefficients 13
	//CurveFitDialog/ w[0] = z0
	//CurveFitDialog/ w[1] = A
	//CurveFitDialog/ w[2] = B
	//CurveFitDialog/ w[3] = cor1
	//CurveFitDialog/ w[4] = cor2
	//CurveFitDialog/ w[5] = x1
	//CurveFitDialog/ w[6] = x2
	//CurveFitDialog/ w[7] = y1
	//CurveFitDialog/ w[8] = y2
	//CurveFitDialog/ w[9] = xw1
	//CurveFitDialog/ w[10] = xw2
	//CurveFitDialog/ w[11] = yw1
	//CurveFitDialog/ w[12] = yw2

	return w[0]+(w[1]*exp((-1/(2*(1-(w[3])^2)))*((((x-w[5])/w[9])^2)+(((y-w[7])/w[11])^2)-((2*w[3]*(x-w[5])*(y-w[7]))/(w[9]*w[11])))))+ (w[2]*exp((-1/(2*(1-(w[4])^2)))*((((x-w[6])/w[10])^2)+(((y-w[8])/w[12])^2)-((2*w[4]*(x-w[6])*(y-w[8]))/(w[10]*w[12])))))
End


// This function fits the dumbbell peak positions that are in x_loc_d and y_loc_d.
// This only works for 2x2D Gaussian fits. 
// Single atom columns need a different routine.
// Inputs: image, x_loc_d, y_loc_d, and gaussian fit size x size.
function DumbbellGaussianFit(image, x_loc_d, y_loc_d, size_d_x, size_d_y)
	
	wave image, x_loc_d, y_loc_d
	variable size_d_x, size_d_y
	variable num_peaks_d = DimSize(x_loc_d,0)
	
	variable scalex = DimDelta(image, 0)
	variable scaley = DimDelta(image, 1)
	
	variable offsetx = DimOffset(image, 0)
	variable offsety = DimOffset(image, 1)
	
	Make/D/N=(num_peaks_d/2)/O x_start
	Make/D/N=(num_peaks_d/2)/O x_finish
	Make/D/N=(num_peaks_d/2)/O y_start
	Make/D/N=(num_peaks_d/2)/O y_finish
	
	variable size_x = (round((size_d_x/2)/scalex))*scalex
	variable size_y = (round((size_d_y/2)/scaley))*scaley

	variable xpos1
	variable xpos2
	variable ypos1
	variable ypos2

	Make/D/N=(num_peaks_d/2)/O xc, yc
	variable i=0
	for(i=0; i<DimSize(x_loc_d,0); i+=2)
		xpos1 = x_loc_d[i] - offsetx
		xpos2 = x_loc_d[i+1] - offsetx
		ypos1 = y_loc_d[i] - offsety
		ypos2 = y_loc_d[i+1] - offsety
		
		xc[i/2] = ((round((xpos1 + xpos2)/(2*scalex)))*scalex)+offsetx
		yc[i/2] = ((round((ypos1 + ypos2)/(2*scaley)))*scaley)+offsety

		x_start[i/2] = xc[i/2] - size_x
		x_finish[i/2] = xc[i/2] + size_x - scalex
		y_start[i/2] = yc[i/2] - size_y
		y_finish[i/2] = yc[i/2] + size_y - scaley
	endfor	

	Make/D/O/N=(num_peaks_d) z0_d, A_d, x0_d, xW_d, y0_d, yW_d, cor_d
	Make/D/O/N=(num_peaks_d) sigma_z0_d, sigma_A_d, sigma_x0_d, sigma_xW_d, sigma_y0_d, sigma_yW_d, sigma_cor_d
	
	variable V_FitError
	
	// Make initial guesses for fit.
	variable cor1s = 0.01, cor2s = 0.01
	variable xw1s = 0.25, xw2s = 0.25, yw1s = 0.25, yw2s = 0.25
	
	i=0
	for(i=0; i<num_peaks_d; i+=2)
		//crop out, save, and setscale of area to be fitted
		Duplicate/O image im
		cropimage(im, x_start[i/2], y_start[i/2], x_finish[i/2], y_finish[i/2])
		if( i == 0 )
			Make/O/N=(Dimsize(im,0), Dimsize(im,1), (num_peaks_d/2)) im_d
		endif
		imagetransform /INSI=im /INSX=0 /INSY=0 /P=(i/2) insertImage im_d
		setscale/P x, x_start[i/2], scalex, im
		setscale/P y, y_start[i/2], scaley, im
		
		Make/O/N=(DimSize(im,0), DimSize(im,1)) residual
		
		Duplicate/O/D im noise
		noise = sqrt(im)
		
		// Make initial guesses for fit.
		Make/D/N=13/O W_coef
		Make/D/O/N=27 M_WaveStats
		Make/D/O/N=(7) W_sigma
		wavestats/Q/W im
		variable z0s = M_WaveStats[10]
		variable As = image(x_loc_d[i])(y_loc_d[i])
		variable Bs = image(x_loc_d[i+1])(y_loc_d[i+1])
		W_coef[0] = {z0s, As, Bs, cor1s, cor2s, x_loc_d[i], x_loc_d[i+1], y_loc_d[i], y_loc_d[i+1], xw1s, xw2s, yw1s, yw2s}
		
		V_FitError = 0
		FuncFitMD/N/NTHR=0/W=2/Q  Gaus2Dx2, W_coef, im/R=residual/W=noise/D
		
		if (V_FitError !=0)
			print "error in fit - slice #", i, "- V_FitError =", V_FitError
		endif
		
		z0_d[i] = W_coef(0)
		A_d[i] = W_coef(1)
		x0_d[i] = W_coef(5)
		xW_d[i] = W_coef(9)
		y0_d[i] = W_coef(7)
		yW_d[i] = W_coef(11)
		cor_d[i] = W_coef(3)
		z0_d[i+1] = W_coef(0)
		A_d[i+1] = W_coef(2)
		x0_d[i+1] = W_coef(6)
		xW_d[i+1] = W_coef(10)
		y0_d[i+1] = W_coef(8)
		yW_d[i+1] = W_coef(12)
		cor_d[i+1] = W_coef(4)
		sigma_z0_d[i] = W_sigma(0)
		sigma_A_d[i] = W_sigma(1)
		sigma_x0_d[i] = W_sigma(5)
		sigma_xW_d[i] = W_sigma(9)
		sigma_y0_d[i] = W_sigma(7)
		sigma_yW_d[i] = W_sigma(11)
		sigma_cor_d[i] = W_sigma(3)
		sigma_z0_d[i+1] = W_sigma(0)
		sigma_A_d[i+1] = W_sigma(2)
		sigma_x0_d[i+1] = W_sigma(6)
		sigma_xW_d[i+1] = W_sigma(10)
		sigma_y0_d[i+1] = W_sigma(8)
		sigma_yW_d[i+1] = W_sigma(12)
		sigma_cor_d[i+1] = W_sigma(4)
		
		//save the fit images
		wave fit_im = $"fit_im"
		if( i == 0)
			Make/O/D/N=(Dimsize(fit_im,0), Dimsize(fit_im,1), (num_peaks_d/2)) fits_d
		endif	
		imagetransform /INSI=fit_im /INSX=0 /INSY=0 /P=(i/2) insertImage fits_d
		killwaves  fit_im
		
		//save the residual images
		if (i==0)
			Make/O/N=(DimSize(im,0), DimSize(im,1), num_peaks_d/2) residual_d
		endif	
		imagetransform /INSI=residual /INSX=0 /INSY=0 /P=(i/2) insertImage residual_d
		killwaves residual
		
		killwaves im, noise
		killwaves W_coef, M_WaveStats, W_sigma
	endfor

	killwaves x_start, x_finish, y_start, y_finish, xc, yc
end



// This function fits all the atoms in a whole image 2D Gaussian functions.
// The atom position need to be supplied in x_loc_s, y_loc_s, x_loc_d, and y_loc_d
// x_loc_s and y_loc_s are the x and y positional of all the single atom columns.
// x_loc_d and y_loc_d are the x and y positional of all the dumbbell atom columns.
// size_s and size_d are the sides of a square to
function FitImage(image, x_loc_s, y_loc_s, x_loc_d, y_loc_d, size_s, size_d_x, size_d_y)
	wave image, x_loc_s, y_loc_s, x_loc_d, y_loc_d
	variable size_s, size_d_x, size_d_y

	variable num_peaks_s  = DimSize(x_loc_s,0)
	variable num_peaks_d  = DimSize(x_loc_d,0)

	if( num_peaks_s != 0)
		GaussianFit(image, x_loc_s, y_loc_s, size_s)
	endif
	
	if( num_peaks_d != 0)
		DumbbellGaussianFit(image, x_loc_d, y_loc_d, size_d_x, size_d_y)
	endif
	
end



// Does the FitImage function on a stack of images of the same structure. 
// startimage is the image number to start the analysis at. #s start at 0. Make 0 if you don't want to cut any of the beginning images.
function FitStack(image_stack, startimage, x_loc_s, y_loc_s, x_loc_d, y_loc_d, size_s, size_d_x, size_d_y)

	wave image_stack, x_loc_s, y_loc_s, x_loc_d, y_loc_d
	variable startimage, size_s, size_d_x, size_d_y
	
	//chops the image stack to start analysis at an image other than the first image.
	Duplicate/O image_stack Image_stack_new
	if( startimage != 0)
		DeletePoints/M=2 0, startimage, image_stack_new
	endif
	
	variable num_images = DimSize(image_stack_new, 2)
	variable num_peaks_s  = DimSize(x_loc_s,0)
	variable num_peaks_d  = DimSize(x_loc_d,0)
	variable num_peaks = num_peaks_s + num_peaks_d
	
	if( num_peaks_s == 0)
		print "No single columns"
	else
		Make/O/D/N=(num_peaks_s, num_images) z0_s_stack, A_s_stack, x0_s_stack, xW_s_stack, y0_s_stack, yW_s_stack, cor_s_stack 
		Make/O/D/N=(num_peaks_s, num_images) sig_z0_s_stack, sig_A_s_stack, sig_x0_s_stack, sig_xW_s_stack, sig_y0_s_stack, sig_yW_s_stack, sig_cor_s_stack 
	endif
	
	if( num_peaks_d == 0)
		print "No dumbell columns"
	else
		Make/O/D/N=(num_peaks_d, num_images) z0_d_stack, A_d_stack, x0_d_stack, xW_d_stack, y0_d_stack, yW_d_stack, cor_d_stack 
		Make/O/D/N=(num_peaks_d, num_images) sig_z0_d_stack, sig_A_d_stack, sig_x0_d_stack, sig_xW_d_stack, sig_y0_d_stack, sig_yW_d_stack, sig_cor_d_stack
	endif
	
	variable scalex = DimDelta(image_stack, 0)
	variable scaley = DimDelta(image_stack, 1)
	variable offsetx = DimOffset(image_stack, 0)
	variable offsety = DimOffset(image_stack, 1)
	
	variable i=0
	for(i=0; i<num_images; i+=1)
		Imagetransform/p=(i) getplane image_stack_new
		wave image = $"M_ImagePlane"
		setscale/P x, offsetx, scalex, image
		setscale/P y, offsety, scaley, image
		
		FitImage(image, x_loc_s, y_loc_s, x_loc_d, y_loc_d, size_s, size_d_x, size_d_y)
		killwaves M_ImagePlane
		
		//save fit parameters for each atom column in each frame of stack
		if( num_peaks_s != 0)
			wave z0_s = $"z0_s"
			wave A_s = $"A_s"
			wave x0_s = $"x0_s"
			wave xW_s = $"xW_s"
			wave y0_s = $"y0_s"
			wave yW_s = $"yW_s"
			wave cor_s = $"cor_s"
			wave sigma_z0_s = $"sigma_z0_s"
			wave sigma_A_s = $"sigma_A_s"
			wave sigma_x0_s = $"sigma_x0_s"
			wave sigma_xW_s = $"sigma_xW_s"
			wave sigma_y0_s = $"sigma_y0_s"
			wave sigma_yW_s = $"sigma_yW_s"
			wave sigma_cor_s = $"sigma_cor_s"
			wave residual_s = $"residual_s"
			wave fits_s = $"fits_s"
			wave im_s = $"im_s"
			
			z0_s_stack[][i] = z0_s[p] 
			A_s_stack[][i] = A_s[p] 
			x0_s_stack[][i] = x0_s[p]  
			xW_s_stack[][i] = xW_s[p]  
			y0_s_stack[][i] = y0_s[p]  
			yW_s_stack[][i] = yW_s[p]  
			cor_s_stack[][i] = cor_s[p] 
			sig_z0_s_stack[][i] = sigma_z0_s[p] 
			sig_A_s_stack[][i] = sigma_A_s[p] 
			sig_x0_s_stack[][i] = sigma_x0_s[p]  
			sig_xW_s_stack[][i] = sigma_xW_s[p]  
			sig_y0_s_stack[][i] = sigma_y0_s[p]  
			sig_yW_s_stack[][i] = sigma_yW_s[p]  
			sig_cor_s_stack[][i] = sigma_cor_s[p]  
		endif
		if( num_peaks_d != 0)
			wave z0_d = $"z0_d"
			wave A_d = $"A_d"
			wave x0_d = $"x0_d"
			wave xW_d = $"xW_d"
			wave y0_d = $"y0_d"
			wave yW_d = $"yW_d"
			wave cor_d = $"cor_d"
			wave sigma_z0_d = $"sigma_z0_d"
			wave sigma_A_d = $"sigma_A_d"
			wave sigma_x0_d = $"sigma_x0_d"
			wave sigma_xW_d = $"sigma_xW_d"
			wave sigma_y0_d = $"sigma_y0_d"
			wave sigma_yW_d = $"sigma_yW_d"
			wave sigma_cor_d = $"sigma_cor_d"
			wave residual_d = $"residual_d"
			wave fits_d = $"fits_d"
			wave im_d = $"im_d"
			
			z0_d_stack[][i] = z0_d[p] 
			A_d_stack[][i] = A_d[p] 
			x0_d_stack[][i] = x0_d[p]  
			xW_d_stack[][i] = xW_d[p]  
			y0_d_stack[][i] = y0_d[p]  
			yW_d_stack[][i] = yW_d[p]  
			cor_d_stack[][i] = cor_d[p] 
			sig_z0_d_stack[][i] = sigma_z0_d[p] 
			sig_A_d_stack[][i] = sigma_A_d[p] 
			sig_x0_d_stack[][i] = sigma_x0_d[p]  
			sig_xW_d_stack[][i] = sigma_xW_d[p]  
			sig_y0_d_stack[][i] = sigma_y0_d[p]  
			sig_yW_d_stack[][i] = sigma_yW_d[p]  
			sig_cor_d_stack[][i] = sigma_cor_d[p]  
		endif
		
		//save the fits, residuals, and original image fit areas.
		if( num_peaks_s != 0)
			if (i == 0)
				variable m=0
				for(m=0; m<num_peaks_s; m+=1)
					string imS
					sprintf imS, "image_s_%d", m
					Make/O/N=(DimSize(im_s, 0), DimSize(im_s, 1), num_images) $imS
				endfor
			endif
			m=0
			for(m=0; m<num_peaks_s; m+=1)
				sprintf imS, "image_s_%d", m
				Imagetransform/p=(m) getplane im_s
				wave im = $"M_ImagePlane"
				wave ImageS = $imS
				imagetransform /INSI=im /INSX=0 /INSY=0 /P=(i) insertImage ImageS
				killwaves M_ImagePlane
			endfor
			
			if (i == 0)
				m=0
				for(m=0; m<num_peaks_s; m+=1)
					string fitS
					sprintf fitS, "fit_s_%d", m
					Make/O/D/N=(DimSize(fits_s, 0), DimSize(fits_s, 1), num_images) $fitS
				endfor
			endif
			m=0
			for(m=0; m<num_peaks_s; m+=1)
				sprintf fitS, "fit_s_%d", m
				Imagetransform/p=(m) getplane fits_s
				wave fit = $"M_ImagePlane"
				wave FitsS = $fitS
				imagetransform /INSI=fit /INSX=0 /INSY=0 /P=(i) insertImage FitsS
				killwaves M_ImagePlane
			endfor
		
			if (i == 0)
				m=0
				for(m=0; m<num_peaks_s; m+=1)
					string residS
					sprintf residS, "residual_s_%d", m
					Make/O/N=(DimSize(residual_s, 0), DimSize(residual_s, 1), num_images) $residS
				endfor
			endif
			m=0
			for(m=0; m<num_peaks_s; m+=1)
				sprintf residS, "residual_s_%d", m
				Imagetransform/p=(m) getplane residual_s
				wave resid = $"M_ImagePlane"
				wave ResidualS = $residS
				imagetransform /INSI=resid /INSX=0 /INSY=0 /P=(i) insertImage ResidualS
				killwaves M_ImagePlane
			endfor
		endif
		if( num_peaks_d != 0)
			if (i == 0)
				m=0
				for(m=0; m<num_peaks_d/2; m+=1)
					string imD
					sprintf imD, "image_d_%d", m
					Make/O/N=(DimSize(im_d, 0), DimSize(im_d, 1), num_images) $imD
				endfor
			endif
			m=0
			for(m=0; m<num_peaks_d/2; m+=1)
				sprintf imD, "image_d_%d", m
				Imagetransform/p=(m) getplane im_d
				wave im = $"M_ImagePlane"
				wave ImageD = $imD
				imagetransform /INSI=im /INSX=0 /INSY=0 /P=(i) insertImage ImageD
				killwaves M_ImagePlane
			endfor
			
			if (i == 0)
				m=0
				for(m=0; m<num_peaks_d/2; m+=1)
					string fitD
					sprintf fitD, "fit_d_%d", m
					Make/O/D/N=(DimSize(fits_d, 0), DimSize(fits_d, 1), num_images) $fitD
				endfor
			endif
			m=0
			for(m=0; m<num_peaks_d/2; m+=1)
				sprintf fitD, "fit_d_%d", m
				Imagetransform/p=(m) getplane fits_d
				wave fit = $"M_ImagePlane"
				wave FitsD = $fitD
				imagetransform /INSI=fit /INSX=0 /INSY=0 /P=(i) insertImage FitsD
				killwaves M_ImagePlane
			endfor
		
			if (i == 0)
				m=0
				for(m=0; m<num_peaks_d/2; m+=1)
					string residD
					sprintf residD, "residual_d_%d", m
					Make/O/N=(DimSize(residual_d, 0), DimSize(residual_d, 1), num_images) $residD
				endfor
			endif
			m=0
			for(m=0; m<num_peaks_d/2; m+=1)
				sprintf residD, "residual_d_%d", m
				Imagetransform/p=(m) getplane residual_d
				wave resid = $"M_ImagePlane"
				wave ResidualD = $residD
				imagetransform /INSI=resid /INSX=0 /INSY=0 /P=(i) insertImage ResidualD
				killwaves M_ImagePlane
			endfor
		endif

		killwaves residual_s, fits_s, im_s
		killwaves residual_d, fits_d, im_d
		killwaves z0_s, A_s, x0_s, xW_s, y0_s, yW_s, cor_s, sigma_z0_s, sigma_A_s, sigma_x0_s, sigma_xW_s, sigma_y0_s, sigma_yW_s, sigma_cor_s
		killwaves z0_d, A_d, x0_d, xW_d, y0_d, yW_d, cor_d, sigma_z0_d, sigma_A_d, sigma_x0_d, sigma_xW_d, sigma_y0_d, sigma_yW_d, sigma_cor_d
	endfor
end



