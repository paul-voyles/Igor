#pragma rtGlobals=1		// Use modern global access method.



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


// This function finds the peak positions using the analyze particle IGOR function and reports the locations in x_loc and y_loc.
// For Si dumbells, it reports the center of the dumbells. To change the mimimun area of the particles change the A flag in the ImageAnalyzeParticle function.
// Inputs: image
function PeakPositions(image)
	
	wave image
	
	//NewImage image
	
	ImageThreshold/I/M=(1)/Q image
	ImageAnalyzeParticles /E/W/Q/F/M=3/A=50/EBPC stats, root:M_ImageThresh
	
	duplicate/O W_xmin x_loc	
	duplicate/O W_ymin y_loc
	duplicate/O W_xmin xmin	
	duplicate/O W_ymin ymin
	duplicate/O W_xmax xmax	
	duplicate/O W_ymax ymax
	
	x_loc = (xmax + xmin) / 2
	y_loc = (ymax + ymin) / 2
	
	//appendtograph/t y_loc vs x_loc
	//ModifyGraph mode=2
	
	killwaves W_ImageObjArea, W_SpotX, W_SpotY, W_circularity, W_rectangularity, W_ImageObjPerimeter, M_Moments, M_RawMoments
	killwaves W_BoundaryX, W_BoundaryY, W_BoundaryIndex, W_xmin, W_xmax, W_ymin, W_ymax, xmin, xmax, ymin, ymax
end


// This function fits the single peak positions that are in x_loc_s and y_loc_s.
// This only works for 2D Gaussian fits. 
// Dumbbells need a different routine.
// Inputs: image, x_loc_s, y_loc_s, and gaussian fit size x size.
function GaussianFit(image, x_loc, y_loc, size, wiggle)
	
	wave image, x_loc, y_loc
	variable size, wiggle
	variable num_peaks = DimSize(x_loc,0)
	
	Make/o/n=(DimSize(x_loc,0)) x_start
	Make/o/n=(DimSize(x_loc,0)) x_finish
	Make/o/n=(DimSize(x_loc,0)) y_start
	Make/o/n=(DimSize(x_loc,0)) y_finish
	
	x_start = x_loc - round(size/2)
	x_finish = x_loc + round(size/2)
	y_start = y_loc - round(size/2)
	y_finish = y_loc + round(size/2)

	Make/D/O/N=(num_peaks) z0, A, x0, xW, y0, yW, cor
	Make/D/O/N=(num_peaks) sigma_z0, sigma_A, sigma_x0, sigma_xW, sigma_y0, sigma_yW, sigma_cor
	
	variable V_FitError
	
	variable i=0
	for(i=0; i<num_peaks; i+=1)
		//crop out, save, and setscale of area to be fitted
		Duplicate/O image im
		cropimage(im, round(x_start[i]), round(y_start[i]), round(x_finish[i]-1), round(y_finish[i]-1))
		setscale/P x, x_start[i], 1, im
		setscale/P y, y_start[i], 1, im
		
		Make/O/N=(DimSize(im,0), DimSize(im,1)) residual
		
		Duplicate/O/D im noise
		noise = sqrt(im)
		
		V_FitError = 0
		if(wiggle == 0)
			CurveFit/N/NTHR=0/W=2/Q gauss2D, im/R=residual/W=noise
		else	
			OneGaussFitWiggle(image, x_loc[i], y_loc[i], size, wiggle)
		endif
		
		wave W_coef = $"W_coef"
		wave W_sigma = $"W_sigma"
		
		if (V_FitError !=0)
			print "error in fit - atom #", i, "- V_FitError =", V_FitError
			z0[i] = NaN
			A[i] = NaN
			x0[i] = NaN
			xW[i] = NaN
			y0[i] = NaN
			yW[i] = NaN
			cor[i] = NaN
			sigma_z0[i] = NaN
			sigma_A[i] = NaN
			sigma_x0[i] = NaN
			sigma_xW[i] = NaN
			sigma_y0[i] = NaN
			sigma_yW[i] = NaN
			sigma_cor[i] = NaN	
		else
			z0[i] = W_coef[0]
			A[i] = W_coef[1]
			x0[i] = W_coef[2]
			xW[i] = W_coef[3]
			y0[i] = W_coef[4]
			yW[i] = W_coef[5]
			cor[i] = W_coef[6]
			sigma_z0[i] = W_sigma[0]
			sigma_A[i] = W_sigma[1]
			sigma_x0[i] = W_sigma[2]
			sigma_xW[i] = W_sigma[3]
			sigma_y0[i] = W_sigma[4]
			sigma_yW[i] = W_sigma[5]
			sigma_cor[i] = W_sigma[6]	
		endif
		
		//save the residual images
		wave residual = $"residual"
		if (i==0)
			Make/O/N=(DimSize(residual,0), DimSize(residual,1), num_peaks) residual_stack
		endif	
		imagetransform /INSI=residual /INSX=0 /INSY=0 /P=(i) insertImage residual_stack
		killwaves residual
		
	endfor	
	killwaves W_coef, W_sigma
	killwaves x_start, x_finish, y_start, y_finish
	killwaves im, noise
end


// Fits the region size x size about the center position (x_loc, y_loc) to a 2D Guassian,
// but repeats the fit shifting the region around from -wiggle to +wiggle in x and y.  
// The fit with the minimum chisq is considered the best, and the coefficients are maintained
// in waves W_Coef and W_Sigma. All other fitting coefficients are discarded.
function OneGaussFitWiggle(image, x_loc, y_loc, size, wiggle)
	wave image
	variable x_loc, y_loc, size, wiggle

	variable x_s, x_f, y_s, y_f
	
	variable V_chisq_min = Inf
	
	variable V_FitError
	Make/O/N=7 keep_coef, keep_sigma
	Make/O/N=(size,size)  keep_res
	variable i, j, m
	m=0
	for(i=-1.0*wiggle; i<=wiggle; i+=1)
		for(j=-1.0*wiggle; j<=wiggle; j+=1)
			
			x_s = round(x_loc) - round(size/2)+i
			x_f = round(x_loc) + round(size/2)+i
			y_s = round(y_loc) - round(size/2)+j
			y_f = round(y_loc) + round(size/2)+j

			//crop out, save, and setscale of area to be fitted
			Duplicate/O image im
			cropimage(im, x_s, y_s, x_f-1, y_f-1)
			setscale/P x, x_s, 1, im
			setscale/P y, y_s, 1, im
		
			Make/O/N=(DimSize(im,0), DimSize(im,1)) residual
		
			Duplicate/O/D im noise
			noise = sqrt(im)

			V_FitError = 0
			
			CurveFit/N/NTHR=0/W=2/Q gauss2D, im/R=residual/W=noise
			
			wave W_coef = $"W_coef"
			wave W_sigma = $"W_sigma"
			
			if(V_FitError == 0)	  // if an error occured during fitting, ignore the results
				if(V_chisq < V_chisq_min)	// keep solutions with lower chisq
					if( (W_coef[2] > x_s) && (W_coef[2] < x_f) )	// keep solutions with a center inside the fit box
						if( (W_coef[4] > y_s) && (W_coef[4] < y_f) )
							keep_coef = W_coef 
							keep_sigma = W_sigma 
							keep_res = residual
							V_chisq_min = V_chisq
							m=1
						endif
					endif
				endif
			endif
		endfor
	endfor
	
	Duplicate/o keep_coef W_coef
	Duplicate/o keep_sigma W_sigma
	Duplicate/o keep_res residual
	
	if(m == 0)
		W_coef = NaN
		W_sigma = NaN
		residual = NaN
	endif
	
	Killwaves keep_coef, keep_sigma, keep_res	
end


// This one is being used in the Precision and Stack functions. It should work for all images if the inputs are correct.
// Inputs: x_space & y_space=spacing you expect, space_delta = +/- error in spacings (2), x_angle & y_angle 
// = expected angle in deg of vector in image(0 deg points right), angle_delta = +/- err in angle (5)
function Separation(x_space, y_space, space_delta, x_angle, y_angle, angle_delta)

	variable x_space, y_space, space_delta, x_angle, y_angle, angle_delta
	
	wave tmp_x0 = $"x0"
	wave tmp_y0 = $"y0"
		if(!waveexists(tmp_x0))
			print "You are an idoit, x0 doesn't exist.\r"
			return 0
		endif
		if(!waveexists(tmp_y0))
			print "You are an idoit, y0 doesn't exist.\r"
			return 0
		endif
	
	variable num_peaks = DimSize(tmp_x0,0)
	
	Make/O/N=0 x_distances, y_distances

	variable i, j, xm, ym, d, angle
		for(i=0; i<num_peaks-1; i+=1)
			for(j=i+1; j<num_peaks; j+=1)
				xm = abs(tmp_x0(j) - tmp_x0(i))
				ym = abs(tmp_y0(j) - tmp_y0(i))
				d = sqrt((xm^2) + (ym^2))
				angle = atan(ym/xm) * 180 / Pi
				
				if(d > (x_space - space_delta) && d < (x_space + space_delta))
					if(angle > (x_angle - angle_delta) && angle < (x_angle + angle_delta))
						InsertPoints 0, 1, x_distances
						x_distances[x2pnt(x_distances,0)] = d
					endif
				endif
			endfor
		endfor

	//Make/O/N=27 M_WaveStats
	wavestats/Q/W x_distances
	wave M_WaveStats = $"M_WaveStats"
	
	variable xaverage = M_WaveStats(3)
	variable xstdev = M_WaveStats(4)
	
//	print "x average =" 
//	print xaverage
//	print "x stdev =" 
//	print xstdev
	
		for(i=0; i<num_peaks-1; i+=1)
			for(j=i+1; j<num_peaks; j+=1)
				xm = abs(tmp_x0(j) - tmp_x0(i))
				ym = abs(tmp_y0(j) - tmp_y0(i))
				d = sqrt((xm^2) + (ym^2))
				angle = atan(ym/xm) * 180 / Pi
				
				if(d > (y_space - space_delta) && d < (y_space + space_delta))
					if(angle > (y_angle - angle_delta) && angle < (y_angle + angle_delta))
						InsertPoints 0, 1, y_distances
						y_distances[x2pnt(y_distances,0)] = d
					endif
				endif
			endfor
		endfor

	wavestats/Q/W y_distances
	
	variable yaverage = M_WaveStats(3)
	variable ystdev = M_WaveStats(4)
	
//	print "y average =" 
//	print yaverage
//	print "y stdev =" 
//	print ystdev

	make/o/n=4 prec
	prec[0] = xaverage
	prec[1] = xstdev
	prec[2] = yaverage
	prec[3] = ystdev
	
	//killwaves M_WaveStats
end



// Does the peak find, 2dGaus fit, and x and y precision calculations on one image of single atom columns.
// Outputs sep_out with x and y avergae and stdev.
function Precision(image, size, x_space, y_space, space_delta, x_angle, y_angle, angle_delta)

	wave image
	variable size, x_space, y_space, space_delta, x_angle, y_angle, angle_delta
	
	PeakPositions(image)
	
	wave x_loc =$"x_loc"
	wave y_loc =$"y_loc"
	if(!WaveExists(x_loc) || !WaveExists(y_loc))
		print "PeakPositions failed to create x_loc and/or y_loc.  Exiting.\r"
		return 0
	endif
	
	GaussianFit(image, x_loc, y_loc, size, 1)	// hard coded for 1 fitting box wiggle
	
	Separation(x_space, y_space, space_delta, x_angle, y_angle, angle_delta)	
end



// Does the peak find, and 2dGaus fit on one image of single atom columns.
// Outputs Gaussian Fit parameters for each atom column position.
function FitPosition(unblurred, blurred, size)

	wave unblurred, blurred
	variable size
	
	variable xsize = DimSize(unblurred,0)
	variable ysize = DimSize(unblurred,1)
	variable scalex = DimDelta(unblurred, 0)
	variable scaley = DimDelta(unblurred, 1)
	variable offsetx = DimOffset(unblurred, 0)
	variable offsety = DimOffset(unblurred, 1)
	
	duplicate/o unblurred unblurred_n
	duplicate/o blurred blurred_n
	
	setscale/p x, 0, 1, blurred_n
	setscale/p y, 0, 1, blurred_n
	setscale/p x, 0, 1, unblurred_n
	setscale/p y, 0, 1, unblurred_n
	
	//Blur image_blurred_n
	//Smooth/E=2/F=0 40, image_blurred_n
	
	PeakPositions(blurred_n)
	
	wave x_loc =$"x_loc"
	wave y_loc =$"y_loc"
	if(!WaveExists(x_loc) || !WaveExists(y_loc))
		print "PeakPositions failed to create x_loc and/or y_loc.  Exiting.\r"
		return 0
	endif
	
	GaussianFit(unblurred_n, x_loc, y_loc, size, 1)	// hard coded for zero fitting box wiggle
	
	duplicate/o x0 x0_rescaled
	duplicate/o y0 y0_rescaled
	
	//rescale x_0, y_0
	
	NewImage image_unblurred
	appendtograph/t y0_rescaled vs x0_rescaled
	ModifyGraph mode=3,marker=19,msize=1
	
	killwaves cor, yW, xW, A, z0, sigma_cor, sigma_yW, sigma_xW, sigma_A, sigma_z0, sigma_x0, sigma_y0, unblurred_n, blurred_n
	
end

