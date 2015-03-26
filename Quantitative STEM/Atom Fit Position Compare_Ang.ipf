#pragma rtGlobals=1		// Use modern global access method.



function CompareDumbbellFitPositions(image_s, image_d, x_loc_d, y_loc_d, size_d_x, size_d_y)

	wave image_s, image_d, x_loc_d, y_loc_d
	variable size_d_x, size_d_y

	DumbbellGaussianFit(image_s, x_loc_d, y_loc_d, size_d_x, size_d_y)
	wave x0_d = $"x0_d"
	wave y0_d = $"y0_d"
	duplicate/o x0_d x0_orig_n
	duplicate/o y0_d y0_orig_n
	killwaves z0_d, A_d, x0_d, xW_d, y0_d, yW_d, cor_d
	killwaves sigma_z0_d, sigma_A_d, sigma_x0_d, sigma_xW_d, sigma_y0_d, sigma_yW_d, sigma_cor_d

	DumbbellGaussianFitDouble(image_d, x_loc_d, y_loc_d, size_d_x, size_d_y)
	wave x0_d = $"x0_d"
	wave y0_d = $"y0_d"
	duplicate/o x0_d x0_avg_n
	duplicate/o y0_d y0_avg_n
	killwaves z0_d, A_d, x0_d, xW_d, y0_d, yW_d, cor_d
	killwaves sigma_z0_d, sigma_A_d, sigma_x0_d, sigma_xW_d, sigma_y0_d, sigma_yW_d, sigma_cor_d
	killwaves im_d, fits_d, residual_d

	Duplicate/o x0_orig_n xdiff, ydiff
	xdiff = x0_orig_n - x0_avg_n
	ydiff = y0_orig_n - y0_avg_n

	Wavestats/q/w xdiff 
	wave M_WaveStats = $"M_WaveStats"
	print "AverageX =", V_avg
	print "StdvX =", V_sdev
	print "RMSX =", V_rms

	Wavestats/q/w ydiff 
	wave M_WaveStats = $"M_WaveStats"
	print "AverageY =", V_avg
	print "StdvY =", V_sdev
	print "RMSY =", V_rms

	killwaves M_WaveStats

end



function RMS(x0_orig_n, x0_avg_n, y0_orig_n, y0_avg_n)

	wave x0_orig_n, x0_avg_n, y0_orig_n, y0_avg_n
	
	Duplicate/o x0_orig_n xdiff_scaled, ydiff_scaled
	xdiff_scaled = x0_orig_n - x0_avg_n
	ydiff_scaled = y0_orig_n - y0_avg_n
	
	Wavestats/q/w xdiff_scaled
	wave M_WaveStats = $"M_WaveStats"
	print "AverageX =", V_avg
	print "StdvX =", V_sdev
	print "RMSX =", V_rms

	Wavestats/q/w ydiff_scaled
	wave M_WaveStats = $"M_WaveStats"
	print "AverageY =", V_avg
	print "StdvY =", V_sdev
	print "RMSY =", V_rms

	killwaves M_WaveStats

end




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



// This function fits the dumbbell peak positions that are in x_loc_d and y_loc_d.
// This only works for 2x2D Gaussian fits. 
// Single atom columns need a different routine.
// Inputs: image, x_loc_d, y_loc_d, and gaussian fit size x size.
function DumbbellGaussianFitDouble(image, x_loc_d, y_loc_d, size_d_x, size_d_y)
	
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
			Make/O/D/N=(Dimsize(im,0), Dimsize(im,1), (num_peaks_d/2)) im_d
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


