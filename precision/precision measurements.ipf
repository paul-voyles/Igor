#pragma rtGlobals=1		// Use modern global access method.
//
// Functions for processing STEM image stacks. Peak finding, peak fitting, and precision 
// measurements.
//
// List of Functions:
//
// CountstoElectrons: converts images from counts to number of electrons in each pixel. inputs are "average" which is the image in counts, 
// "numSamples" which is the image of number of frames used in pixel for the averaging, and "gain" which is the gain of the HAADF detector.
//
// PeakPositions: Input image wave name. Uses the IGOR ImageAnalyzeParticle function to find the 
// positions of all the atoms, not counting the atoms on the perimiter of the image.
//  
// DumbbellPeakPositions: This function finds the peak positions using the analyze particle IGOR function and reports the locations in x_loc
// and y_loc. For Si dumbells, it reports the center of each atom in a dumbbell. It has a higher sensitivity to distiguish atoms in a dumbbell.
// Inputs: image wavename
//
// GaussianFit: Input wave name, x_loc, y_loc, and width of fit area in pixels. This functions 
// conducts a 2D Gaussian fit for each peak found in the PeakPositions function which are inputed in x_loc and y_loc.
//
// GaussianFitInitialGuess: same as GaussianFit, including inputs, but this time no weighting is done to the fit and the code automatically 
// determines intial guesses for all the fitting parameters.
//
// Gaus2Dx2(w,x,y) : FitFunc: This function creates the fitting function of the sum of 2 2D Gaussians.
//
// HAADFDumbbellGaussianFit: Same inputs as GaussianFit, but this time tries to fit the 2 times 2D Gaussian fit function to an HAADF dumbbell.
// No weighting is included here because this was mainly used for simulated images.
//
// ABFDumbbellGaussianFit: Same inputs as GaussianFit, but this time tries to fit the 2 times 2D Gaussian fit function to an ABF dumbbell.
// No weighting is included here because this was mainly used for simulated images.
//
// DumbbellGaussianFit: input wave name, x_loc, y_loc, and separation of the dumbell in pixels. This function fits the dumbbell peak positions that
// are in x_loc and y_loc that were found in DumbbellPeakPositions(image). This only works for dumbbell 2D Gaussian fits that are aligned vertically.
//
// Separation: This calculates all the specific separations in the x and y directions. This one is being used in the Precision and Stack functions. 
// It should work for all images if the inputs are correct. Outputs all the spacings from an image, as well as the average and stdev of both x and y spacings.
// Inputs: x_space & y_space=spacing you expect, space_delta = +/- error in spacings (2), x_angle & y_angle 
// = expected angle in deg of vector in image(0 deg points right), angle_delta = +/- err in angle (5)
// 
// Precision: Does the peak find, 2dGaus fit, and x and y precision calculations on one image of single atom columns.
// Outputs sep_out with x and y avergae and stdev.
// inputs: image, size, x_space, y_space, space_delta, x_angle, y_angle, angle_delta
//
// DumbbellPrecision: Does the peak find, 2 x 2dGaus fit, and x and y precision calculations on one image of dumbbell atom columns.
// Outputs sep_out with x and y avergae and stdev.
// inputs: image, dum_sep, x_space, y_space, space_delta, x_angle, y_angle, angle_delta
//
// StackPrecision: Does the peak find, 2dGaus fit, and x and y precision calculation on a stack of images.
// Outputs xav, yav, xstd, ystd, and fit parameters for all column position in each image in stack.
// inputs: image_stack, fit size, x_space, y_space, space_delta, x_angle, y_angle, angle_delta
//
// DumbbellStackPrecision: Does the dumbbell peak find, dumbbell Gaus fit, and x and y precision calculation on a stack of images.
// Outputs xav, yav, xstd, ystd, and fit parameters for all column position in each image in stack.
// inputs: image_stack, dum_sep, x_space, y_space, space_delta, x_angle, y_angle, angle_delta.
//
// SimulationStackPosition, HAADFSimStackPosition, HAADFDumSimStackPosition, ABFDumSimStackPosition: variation of the above 
// two function for simulated images.
//
// StackContrast: Claculates contrast as a function of stack image using image stats, not the fit parameters.
//
// FloatingWindow: Creates an average image from a floating window within a stack, "image_stack" with width "float_width".
//
// OneGaussFitWiggle: Fits a 2D Gaussian to a sub-window of the image, but "wiggles" the image around a bit to find the best fit.
//
// Started 4/30/12, ABY
//Updates: 5/2/12: added weighting in fit to be sqrt of counts. Created Stack functions to calculate precision vs. stack frame.
//5/9/12: Made the Stack function save all the fitting parameter from each stack frame. Prior version replaced each fit with each frame in stack.   
//5/15/12: Made the SimulationPosition, SimulationStackPosition, HAADFSimStackPosition, ABFDumbbellSimStackPosition functions
//5/18/12 Added the function FloatingWindow.
//7/10/12: Consolidated functions, erased unused ones, and added comments to make easier to use in the future.
// 10/11/12:  Added OneGaussFitWiggle, modified GaussianFit to use it, and modified other calls to GaussianFit for consistency, pmv.
// 6/22/13: Modified GaussianFit and OneGaussFitWiggle to discard fits that fail to converge or have a center outside the fit box.  pmv
// 12/5/14: Added GaussianFitStack: find Gaussian atom fit parameters for every image in a series.  pmv
// 12/5/14: Added option to fit to fixed atom center position to GaussianFit for use in GaussianFitStack.  pmv
// 12/9/14: modify PeakPositions to use the current data folder, not hard-coded for root data folder.  pmv
// 12/11/14: added StackMean.  Very similar to SumIntensity.  pmv
// 12/11/14: added "no noise" optional parameter for GaussianFit and GaussianFitStack.  Useful for images with negative pixels, which
// otherwise are NaN in the sqrt(N) noise image and generate errors from the fit code.  pmv.

// This function converts images in counts to electrons
// inputs: original image in HAADF counts, numSamples, and cm = average count in HAADF probe image,
// c0 = average count HAADF dark image, p = probe current in electrons/sec, t = pixel dwell time in seconds.
function CountstoElectrons(average, cm, c0, p, t)
	
	wave average 
	variable cm, c0 , p, t


	Duplicate/O average average_AbsInt
	average_AbsInt = (average - c0)/(cm - c0)
	Duplicate/O average average_electrons
	average_electrons = average_AbsInt * p * t
end


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


// takes an stack of residual images and outputs two 1D waves (res_av and res_stdv) which are the avgs and stdevs from each ind residual.
function ResidualAverage(res)

	wave res
	
	variable size_x = DimSize(res,0)
	variable size_y = DimSize(res,1)
	variable size_z = DimSize(res,2)
	
	make/o/n=(size_z) res_av, res_stdv
	
	variable i=0
	for(i=0; i<size_z; i+=1)
		Imagetransform/p=(i) getplane res
		wave im = $"M_ImagePlane"
		Wavestats/Q/W im 
		wave M_WaveStats = $"M_WaveStats"
		res_av[i] = M_WaveStats[3]
		res_stdv[i] = M_WaveStats[4]
	endfor	
	
	killwaves M_ImagePlane, M_WaveStats
	
end




// takes an image im and crops it with boundary xi yi xf yf defined by the maximum in the numSamples image.  
function NumSamplesCrop(av, numSamples)

	wave av, numSamples
	variable xi=256, yi=256, xf=0, yf=0
	variable xi_temp, yi_temp, xf_temp, yf_temp

	Wavestats/Q/W numSamples 
	wave M_WaveStats = $"M_WaveStats"
	variable maxNum = M_WaveStats[12]
	variable xsize = DimSize(av,0)
	variable ysize = DimSize(av,1)
	
	variable i=0
	variable j=0
	for(i=0; i<xsize; i+=1)
		for(j=0; j<ysize; j+=1)
			if(numSamples[i][j] == maxNum)
				if(i < xi)
					xi = i
				endif
				if(i > xf)
					xf = i
				endif
				if(j < yi)
					yi = j
				endif	
				if(j > yf)
					yf = j
				endif
			endif
		endfor
	endfor	
	
	duplicate/o numSamples numSamples_crop
	duplicate/o av av_crop
	cropimage(numSamples_crop, xi, yi, xf, yf)
	cropimage(av_crop, xi, yi, xf, yf)
	
	Wavestats/Q/W numSamples_crop
	wave M_WaveStats = $"M_WaveStats"
	variable minimum = M_WaveStats[10]
	
	do
		make/o/n=(1,DimSize(av_crop,1)) row1
		make/o/n=(1,DimSize(av_crop,1)) row2
		make/o/n=(DimSize(av_crop,0)) col1
		make/o/n=(DimSize(av_crop,0)) col2
		row1 = numSamples_crop[0][q]
		row2 = numSamples_crop[DimSize(av_crop,0)-1][q]
		col1 = numSamples_crop[p][0]
		col2 = numSamples_crop[p][DimSize(av_crop,1)-1]
		variable mean_row1 = mean(row1)
		variable mean_row2 = mean(row2)
		variable mean_col1 = mean(col1)
		variable mean_col2 = mean(col2)
		variable minimumRow_mean = min(mean_row1, mean_row2)
		variable minimumCol_mean = min(mean_col1, mean_col2)
		variable minimum_mean = min(minimumRow_mean, minimumCol_mean)
		if(mean_row1 == minimum_mean)
			DeletePoints 0,1, numSamples_crop
			DeletePoints 0,1, av_crop
			xi = xi +1
		endif
		if(mean_row2 == minimum_mean)
			DeletePoints (DimSize(av_crop,0)-1),1, numSamples_crop
			DeletePoints (DimSize(av_crop,0)-1),1, av_crop
			xf = xf -1
		endif
		if(mean_col1 == minimum_mean)
			DeletePoints/M=1 0,1, numSamples_crop
			DeletePoints/M=1 0,1, av_crop
			yi = yi +1
		endif
		if(mean_col2 == minimum_mean)
			DeletePoints/M=1 (DimSize(av_crop,1)-1),1, numSamples_crop		
			DeletePoints/M=1 (DimSize(av_crop,1)-1),1, av_crop
			yf = yf -1
		endif
		
		killwaves row1, row2, col1, col2
		
		Wavestats/Q/W numSamples_crop
		wave M_WaveStats = $"M_WaveStats"
		minimum = M_WaveStats[10]
	while(minimum < maxNum)
	print "xi=", xi, "xf=", xf, "yi=", yi, "yf=", yf
	killwaves M_WaveStats
end


// This function converts images to an absolute intensity scale that can be compared to simulations.
// inputs: im is the image, Cm and Co are the average pixel intensity of the HAADF probe image and dark image respectively.
function AbsoluteIntensity(im, Cm, Co)
	
	wave im
	variable Cm, Co

	Duplicate/O im $(NameofWave(im)+"_AbsInt")
	wave/C AbsInt = $(NameofWave(im)+"_AbsInt")
	AbsInt = (im - Co) / (Cm - Co)
	
end



// This function finds the peak positions using the analyze particle IGOR function and reports the locations in x_loc and y_loc.
// For Si dumbells, it reports the center of the dumbells. To change the mimimun area of the particles change the A flag in the ImageAnalyzeParticle function.
// Inputs: image
function PeakPositions(image)
	
	wave image
	
	NewImage image
	
	//Manual Threshold
	//ImageThreshold/I/M=(0)/Q/T=22500 image
	//ImageAnalyzeParticles /E/W/Q/F/M=3/A=1/EBPC stats, root:M_ImageThresh
	
	// Add ROI to ignore pixels in input image that are NaN
	duplicate/O image roi_tmp
	Redimension/b/u roi_tmp
	roi_tmp = numtype(image)
	
	//Auto Threshold
	ImageThreshold/I/M=(1)/Q/R={roi_tmp, 2} image
	ImageAnalyzeParticles /E/W/Q/F/M=3/A=15/EBPC stats, M_ImageThresh
	
	duplicate/O W_xmin x_loc	
	duplicate/O W_ymin y_loc
	duplicate/O W_xmin xmin	
	duplicate/O W_ymin ymin
	duplicate/O W_xmax xmax	
	duplicate/O W_ymax ymax
	
	x_loc = (xmax + xmin) / 2
	y_loc = (ymax + ymin) / 2
	
	appendtograph/t y_loc vs x_loc
	ModifyGraph mode=2
	
	killwaves W_ImageObjArea, W_SpotX, W_SpotY, W_circularity, W_rectangularity, W_ImageObjPerimeter, M_Moments, M_RawMoments
	killwaves W_BoundaryX, W_BoundaryY, W_BoundaryIndex, W_xmin, W_xmax, W_ymin, W_ymax, xmin, xmax, ymin, ymax
	Killwaves roi_tmp
	Killwaves M_ImageTresh, M_Particle
end



// This function finds the peak positions using the analyze particle IGOR function and reports the locations in x_loc and y_loc.
// For Si dumbells, it reports the center of each atom in a dumbbell. It has a higher sensitivity to distiguish atoms in a dumbbell.
// Inputs: image
function DumbbellPeakPositions(image)
	
	wave image
	
	//NewImage image
	
	ImageThreshold/I/M=(1)/Q image
	ImageAnalyzeParticles /E/W/Q/F/M=3/A=100/EBPC stats, root:M_ImageThresh
	
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



function ABFPeakPositions(image)
	
	wave image
	
	//NewImage image
	
	ImageThreshold/M=(1)/Q image
	ImageAnalyzeParticles /E/W/Q/M=3/A=15 stats, root:M_ImageThresh
	
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



// This function fits the peak positions that are in x_loc and y_loc that were found in the previous routine.
// This only works for 2D Gaussian fits. Si dumbbells need a different fitting routine.
// Inputs: image, x_loc, y_loc, and gaussian fit size x size.
function GaussianFit(image, x_loc, y_loc, size, wiggle, [fix_xy, no_noise])
	
	wave image, x_loc, y_loc
	variable size, wiggle, fix_xy, no_noise
	variable half_size = size / 2
	variable num_peaks = DimSize(x_loc,0)
	variable success, V_FitError

	Duplicate/O image noise
	noise = sqrt(image)

	duplicate/O x_loc x_start
	duplicate/O x_loc x_finish
	duplicate/O y_loc y_start
	duplicate/O y_loc y_finish
	
	x_start = x_loc - half_size
	x_finish = x_loc + half_size
	y_start = y_loc - half_size
	y_finish = y_loc + half_size

	Make/O/N=(num_peaks) z0, A, x0, xW, y0, yW, cor
	Make/O/N=(num_peaks) sigma_z0, sigma_A, sigma_x0, sigma_xW, sigma_y0, sigma_yW, sigma_cor
	Make/O/N=6 W_Coef, W_Sigma

	variable i
	for(i=0; i<num_peaks; i+=1)
		
		success = 1
		if(fix_xy == 1)  // fit with held x and y center positions
			if(wiggle == 0)
				if(no_noise == 1)  // leave out Poisson noise weighting from the fit
					V_FitError = 0
					CurveFit/M=2/W=2/Q/O gauss2D, image[x_start[i],x_finish[i]][y_start[i],y_finish[i]]  /I=1 /D  // generate initial guesses
					W_Coef[2] = x_loc[i]  // replace guessed (x,y) position with fixed position
					W_Coef[4] = y_loc[i]
					CurveFit/M=2/W=2/Q/H="0010100" gauss2D, kwCWave = W_Coef, image[x_start[i],x_finish[i]][y_start[i],y_finish[i]] /I=1 /D  // do the fit
					if(V_FitError)
						success = 0
						printf "V_FitError = %d.\t", V_FitError
					endif
				else
					V_FitError = 0
					CurveFit/M=2/W=2/Q/O gauss2D, image[x_start[i],x_finish[i]][y_start[i],y_finish[i]] /W=noise /I=1 /D  // generate initial guesses
					W_Coef[2] = x_loc[i]  // replace guessed (x,y) position with fixed position
					W_Coef[4] = y_loc[i]
					CurveFit/M=2/W=2/Q/H="0010100" gauss2D, kwCWave = W_Coef, image[x_start[i],x_finish[i]][y_start[i],y_finish[i]] /W=noise /I=1 /D  // do the fit
					if(V_FitError)
						success = 0
						printf "V_FitError = %d.\t", V_FitError
					endif
				endif
			else
				printf "Fixed (x,y) atom center position & wiggle fit not yet implemented.\r"
				success = 0
			endif				
		
		else
			if(no_noise == 1)
				printf "Fitting without noise and without fixed (x,y) position is not implemented yet.\r"
				success = 0
			else		
				if(wiggle == 0)
					V_FitError =0
					CurveFit/M=2/W=2/Q gauss2D, image[x_start[i],x_finish[i]][y_start[i],y_finish[i]] /W=noise /I=1 /D
					if(V_FitError)
						success = 0
						printf "V_FitError = %d.\t", V_FitError
					endif
				else	
					success = OneGaussFitWiggle(image, x_loc[i], y_loc[i], size, wiggle)
				endif
			endif
		endif
		
		if(W_coef[2] <= x_start[i] || W_coef[2] >= x_finish[i])
			printf "x0 out of range.\t"
			success = 0
		endif
		if(W_coef[4] <= y_start[i] || W_coef[4] >= y_finish[i])
			printf "y0 out of range.\t"
			success = 0
		endif
		
		if(success)
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
		else
			printf "Fit to point %d at (%g, %g) failed.\r", i, x_loc[i], y_loc[i]
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
		endif

	endfor	

	killwaves W_sigma, W_coef, M_covar, x_finish, x_start, y_finish, y_start, noise
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
			
			CurveFit/N/NTHR=0/W=2/Q gauss2D, im/R=residual//W=noise
			
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


// This function fits the peak positions that are in x_loc and y_loc that were found in the previous routine.
// This only works for 2D Gaussian fits. Si dumbbells need a different fitting routine.
// Inputs: image, x_loc, y_loc, and gaussian fit size x size.
function OldGaussianFit(image, x_loc, y_loc, size, wiggle)
	
	wave image, x_loc, y_loc
	variable size, wiggle
	variable half_size = size / 2
	variable num_peaks = DimSize(x_loc,0)

	Duplicate/O image noise
	noise = sqrt(image)

	duplicate/O x_loc x_start
	duplicate/O x_loc x_finish
	duplicate/O y_loc y_start
	duplicate/O y_loc y_finish
	
	x_start = x_loc - half_size
	x_finish = x_loc + half_size
	y_start = y_loc - half_size
	y_finish = y_loc + half_size

	Make/O/N=(num_peaks) z0, A, x0, xW, y0, yW, cor
	Make/O/N=(num_peaks) sigma_z0, sigma_A, sigma_x0, sigma_xW, sigma_y0, sigma_yW, sigma_cor

	variable numx=DimSize(image,0)
	variable numy=DimSize(image,1)  
	variable numz=DimSize(x_loc,0)
	
	//Make/O/N=(numx, numy,  numz) Residual_stack

	variable i
	for(i=0; i<num_peaks; i+=1)

		if(wiggle == 0)
			CurveFit/M=2/W=2/Q gauss2D, image[x_start[i],x_finish[i]][y_start[i],y_finish[i]] /W=noise /I=1 /D/R
		else	
			OneGaussFitWiggle(image, x_loc[i], y_loc[i], size, wiggle)
		endif
		
		wave W_coef = $"W_coef"
		wave W_sigma = $"W_sigma"
		//wave residual = $"residual"
		
		//Residual_stack[][][i] = residual[p][q]
		
		//killwaves residual
		
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
	endfor	

	//Make/O/N=(numx, numy) Residual_Sum
	//imagetransform sumplanes Residual_stack
	//wave M_SumPlanes = $"M_SumPlanes"
	//Residual_Sum = M_SumPlanes
	//killwaves M_SumPlanes
	
	killwaves W_sigma, W_coef, M_covar
end

// Fits the region size x size about the center position (x_loc, y_loc) to a 2D Guassian,
// but repeats the fit shifting the region around from -wiggle to +wiggle in x and y.  
// The fit with the minimum chisq is considered the best, and the coefficients are maintained
// in waves W_Coef and W_Sigma. All other fitting coefficients are discarded.
function OldOneGaussFitWiggle(image, x_loc, y_loc, size, wiggle)
	wave image
	variable x_loc, y_loc, size, wiggle
	
	Duplicate/O image noise
	Make/O/N=(DimSize(image,0),DimSize(image,1))  residual, single_residual
	
	noise = sqrt(image)

	variable x_start, x_finish, y_start, y_finish, half_size
	half_size = size/2
	
	variable V_chisq_min = Inf
	
	variable i, j, V_FitError
	for(i=-1.0*wiggle; i<=wiggle; i+=1)
		for(j=-1.0*wiggle; j<=wiggle; j+=1)

			x_start = (x_loc - half_size) + i
			x_finish = (x_loc + half_size) + i
			y_start = (y_loc - half_size) + j
			y_finish = (y_loc + half_size) + j

			V_FitError = 0
			CurveFit/N/M=2/W=2/Q gauss2D,image[x_start,x_finish][y_start,y_finish] /W=noise /I=1 /D/R=single_residual
		
			wave W_Coef = $"W_Coef"
			wave W_Sigma = $"W_Sigma"
			//printf "(x, y) = (%g, %g);  chisq = %g\r", W_coef[2], W_coef[4], V_chisq
		
			if(!V_FitError)	  // if an error occured during fitting, ignore the results
				if(V_chisq < V_chisq_min)	// keep solutions with lower chisq
					if( (W_coef[2] > x_start) && (W_coef[2] < x_finish) )	// keep solutions with a center inside the fit box
						if( (W_coef[4] > y_start) && (W_coef[4] < y_finish) )
							Duplicate/O W_Coef keep_coef
							Duplicate/O W_Sigma keep_sigma
							V_chisq_min = V_chisq
							residual = single_residual
						endif
					endif
				endif
			endif
		endfor
	endfor
	
	Duplicate/O keep_coef W_coef
	Duplicate/O keep_sigma W_sigma
	
	Killwaves keep_coef, keep_sigma, noise, single_residual
	
end

//No noise included, but initial guess of all parameters determined.
function GaussianFitInitialGuess(image, x_loc, y_loc, size)
	
	wave image, x_loc, y_loc
	variable size
	variable half_size = size / 2
	variable num_peaks = DimSize(x_loc,0)
	variable xW_i = half_size
	variable yW_i = half_size
	variable cor_i = 0.01 
	
	//Duplicate/O image noise
	//noise = sqrt(image)

	duplicate/O x_loc x_start
	duplicate/O x_loc x_finish
	duplicate/O y_loc y_start
	duplicate/O y_loc y_finish
	
	x_start = x_loc - half_size
	x_finish = x_loc + half_size
	y_start = y_loc - half_size
	y_finish = y_loc + half_size

	Make/O/N=(7) W_sigma
	Make/O/N=(num_peaks) z0, A, x0, xW, y0, yW, cor
	Make/O/N=(num_peaks) sigma_z0, sigma_A, sigma_x0, sigma_xW, sigma_y0, sigma_yW, sigma_cor

	Make/D/O/N=7 coef

	variable i
	for(i=0; i<num_peaks; i+=1)
		Make/O/N=27 M_WaveStats
		Wavestats/Q/W image
		variable z0_i = M_WaveStats(10), A_i = M_WaveStats(12)
		
		coef[0] = {z0_i, A_i, x_loc(i), xW_i, y_loc(i), yW_i, cor_i}
		CurveFit/G/NTHR=0/W=2/Q Gauss2D kwCWave=coef,  image[x_start(i),x_finish(i)][y_start(i),y_finish(i)] /D/R 
		
		z0[x2pnt(z0,i)] = coef(0)
		A[x2pnt(A,i)] = coef(1)
		x0[x2pnt(x0,i)] = coef(2)
		xW[x2pnt(xW,i)] = coef(3)
		y0[x2pnt(y0,i)] = coef(4)
		yW[x2pnt(yW,i)] = coef(5)
		cor[x2pnt(cor,i)] = coef(6)
		sigma_z0[x2pnt(sigma_z0,i)] = W_sigma(0)
		sigma_A[x2pnt(sigma_A,i)] = W_sigma(1)
		sigma_x0[x2pnt(sigma_x0,i)] = W_sigma(2)
		sigma_xW[x2pnt(sigma_xW,i)] = W_sigma(3)
		sigma_y0[x2pnt(sigma_y0,i)] = W_sigma(4)
		sigma_yW[x2pnt(sigma_yW,i)] = W_sigma(5)
		sigma_cor[x2pnt(sigma_cor,i)] = W_sigma(6)	
	endfor	

//appendtograph/t y0 vs x0
//ModifyGraph mode=2
//RemoveFromGraph/Z y_loc
killwaves W_sigma, coef, M_WaveStats
end



//Created the dumbbell 2 x 2DGaussian function for the Si fits
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



//No noise included
function HAADFDumbbellGaussianFit(image, x_loc, y_loc, size)
	
	wave image, x_loc, y_loc
	variable size
	variable num_peaks = DimSize(x_loc,0)
	

	Make/O/N=(7) W_sigma
	Make/O/N=(num_peaks) z0, A, x0, xW, y0, yW, cor
	Make/O/N=(num_peaks) sigma_z0, sigma_A, sigma_x0, sigma_xW, sigma_y0, sigma_yW, sigma_cor

	Make/D/N=13/O W_coef

	Make/O/N=27 M_WaveStats
	wavestats/Q/W image
	variable z0s = M_WaveStats(10), As = M_WaveStats(12), Bs = M_WaveStats(12), cor1s = 0.01, cor2s = 0.01
	variable xw1s = size / 3, xw2s = size / 3, yw1s = size / 3, yw2s = size / 3

	variable i
	for(i=0; i<num_peaks; i+=1)
		
		W_coef[0] = {z0s,As,Bs,cor1s,cor2s,x_loc(0),x_loc(1),y_loc(0),y_loc(1),xw1s,xw2s,yw1s,yw2s}
		FuncFitMD/N/NTHR=0/W=2/Q  Gaus2Dx2 W_coef  image /D/R
		
		z0[x2pnt(z0,0)] = W_coef(0)
		A[x2pnt(A,0)] = W_coef(1)
		x0[x2pnt(x0,0)] = W_coef(5)
		xW[x2pnt(xW,0)] = W_coef(9)
		y0[x2pnt(y0,0)] = W_coef(7)
		yW[x2pnt(yW,0)] = W_coef(11)
		cor[x2pnt(cor,0)] = W_coef(3)
		z0[x2pnt(z0,1)] = W_coef(0)
		A[x2pnt(A,1)] = W_coef(2)
		x0[x2pnt(x0,1)] = W_coef(6)
		xW[x2pnt(xW,1)] = W_coef(10)
		y0[x2pnt(y0,1)] = W_coef(8)
		yW[x2pnt(yW,1)] = W_coef(12)
		cor[x2pnt(cor,1)] = W_coef(4)
		sigma_z0[x2pnt(sigma_z0,2*i)] = W_sigma(0)
		sigma_A[x2pnt(sigma_A,2*i)] = W_sigma(1)
		sigma_x0[x2pnt(sigma_x0,2*i)] = W_sigma(5)
		sigma_xW[x2pnt(sigma_xW,2*i)] = W_sigma(9)
		sigma_y0[x2pnt(sigma_y0,2*i)] = W_sigma(7)
		sigma_yW[x2pnt(sigma_yW,2*i)] = W_sigma(11)
		sigma_cor[x2pnt(sigma_cor,2*i)] = W_sigma(3)
		sigma_z0[x2pnt(sigma_z0,((2*i)+1))] = W_sigma(0)
		sigma_A[x2pnt(sigma_A,((2*i)+1))] = W_sigma(2)
		sigma_x0[x2pnt(sigma_x0,((2*i)+1))] = W_sigma(6)
		sigma_xW[x2pnt(sigma_xW,((2*i)+1))] = W_sigma(10)
		sigma_y0[x2pnt(sigma_y0,((2*i)+1))] = W_sigma(8)
		sigma_yW[x2pnt(sigma_yW,((2*i)+1))] = W_sigma(12)
		sigma_cor[x2pnt(sigma_cor,((2*i)+1))] = W_sigma(4)
	endfor	

killwaves W_sigma, W_coef, M_WaveStats
end



//No noise included
function ABFDumbbellGaussianFit(image, x_loc, y_loc, size)
	
	wave image, x_loc, y_loc
	variable size
	variable num_peaks = DimSize(x_loc,0)
	

	Make/O/N=(7) W_sigma
	Make/O/N=(num_peaks) z0, A, x0, xW, y0, yW, cor
	Make/O/N=(num_peaks) sigma_z0, sigma_A, sigma_x0, sigma_xW, sigma_y0, sigma_yW, sigma_cor

	Make/D/N=13/O W_coef

	Make/O/N=27 M_WaveStats
	wavestats/Q/W image
	variable z0s = M_WaveStats(12), As = M_WaveStats(10), Bs = M_WaveStats(10), cor1s = 0.01, cor2s = 0.01
	variable xw1s = size / 3, xw2s = size / 3, yw1s = size / 3, yw2s = size / 3

	variable i
	for(i=0; i<num_peaks; i+=1)
		
		W_coef[0] = {z0s,As,Bs,cor1s,cor2s,x_loc(0),x_loc(1),y_loc(0),y_loc(1),xw1s,xw2s,yw1s,yw2s}
		FuncFitMD/N/NTHR=0/W=2/Q  Gaus2Dx2 W_coef  image /D/R
		
		z0[x2pnt(z0,0)] = W_coef(0)
		A[x2pnt(A,0)] = W_coef(1)
		x0[x2pnt(x0,0)] = W_coef(5)
		xW[x2pnt(xW,0)] = W_coef(9)
		y0[x2pnt(y0,0)] = W_coef(7)
		yW[x2pnt(yW,0)] = W_coef(11)
		cor[x2pnt(cor,0)] = W_coef(3)
		z0[x2pnt(z0,1)] = W_coef(0)
		A[x2pnt(A,1)] = W_coef(2)
		x0[x2pnt(x0,1)] = W_coef(6)
		xW[x2pnt(xW,1)] = W_coef(10)
		y0[x2pnt(y0,1)] = W_coef(8)
		yW[x2pnt(yW,1)] = W_coef(12)
		cor[x2pnt(cor,1)] = W_coef(4)
		sigma_z0[x2pnt(sigma_z0,2*i)] = W_sigma(0)
		sigma_A[x2pnt(sigma_A,2*i)] = W_sigma(1)
		sigma_x0[x2pnt(sigma_x0,2*i)] = W_sigma(5)
		sigma_xW[x2pnt(sigma_xW,2*i)] = W_sigma(9)
		sigma_y0[x2pnt(sigma_y0,2*i)] = W_sigma(7)
		sigma_yW[x2pnt(sigma_yW,2*i)] = W_sigma(11)
		sigma_cor[x2pnt(sigma_cor,2*i)] = W_sigma(3)
		sigma_z0[x2pnt(sigma_z0,((2*i)+1))] = W_sigma(0)
		sigma_A[x2pnt(sigma_A,((2*i)+1))] = W_sigma(2)
		sigma_x0[x2pnt(sigma_x0,((2*i)+1))] = W_sigma(6)
		sigma_xW[x2pnt(sigma_xW,((2*i)+1))] = W_sigma(10)
		sigma_y0[x2pnt(sigma_y0,((2*i)+1))] = W_sigma(8)
		sigma_yW[x2pnt(sigma_yW,((2*i)+1))] = W_sigma(12)
		sigma_cor[x2pnt(sigma_cor,((2*i)+1))] = W_sigma(4)
	endfor	

killwaves W_sigma, W_coef, M_WaveStats
end


//No noise included
function ABFDumbbellGaussianFitAng(image, x_loc, y_loc, size)
	
	wave image, x_loc, y_loc
	variable size
	variable num_peaks = DimSize(x_loc,0)
	

	Make/O/N=(7) W_sigma
	Make/O/N=(num_peaks) z0, A, x0, xW, y0, yW, cor
	Make/O/N=(num_peaks) sigma_z0, sigma_A, sigma_x0, sigma_xW, sigma_y0, sigma_yW, sigma_cor

	Make/D/N=13/O W_coef
	Make/O/N=27 M_WaveStats
	wavestats/Q/W image
	variable z0s = M_WaveStats(12), As = M_WaveStats(10), Bs = M_WaveStats(10), cor1s = 0.01, cor2s = 0.01
	variable xw1s = size / 3, xw2s = size / 3, yw1s = size / 3, yw2s = size / 3

	variable i
	for(i=0; i<num_peaks; i+=1)
		
		W_coef[0] = {z0s,As,Bs,cor1s,cor2s,x_loc(0),x_loc(1),y_loc(0),y_loc(1),xw1s,xw2s,yw1s,yw2s}
		FuncFitMD/N/NTHR=0/W=2/Q  Gaus2Dx2 W_coef  image
		
		z0[x2pnt(z0,0)] = W_coef(0)
		A[x2pnt(A,0)] = W_coef(1)
		x0[x2pnt(x0,0)] = W_coef(5)
		xW[x2pnt(xW,0)] = W_coef(9)
		y0[x2pnt(y0,0)] = W_coef(7)
		yW[x2pnt(yW,0)] = W_coef(11)
		cor[x2pnt(cor,0)] = W_coef(3)
		z0[x2pnt(z0,1)] = W_coef(0)
		A[x2pnt(A,1)] = W_coef(2)
		x0[x2pnt(x0,1)] = W_coef(6)
		xW[x2pnt(xW,1)] = W_coef(10)
		y0[x2pnt(y0,1)] = W_coef(8)
		yW[x2pnt(yW,1)] = W_coef(12)
		cor[x2pnt(cor,1)] = W_coef(4)
		sigma_z0[x2pnt(sigma_z0,2*i)] = W_sigma(0)
		sigma_A[x2pnt(sigma_A,2*i)] = W_sigma(1)
		sigma_x0[x2pnt(sigma_x0,2*i)] = W_sigma(5)
		sigma_xW[x2pnt(sigma_xW,2*i)] = W_sigma(9)
		sigma_y0[x2pnt(sigma_y0,2*i)] = W_sigma(7)
		sigma_yW[x2pnt(sigma_yW,2*i)] = W_sigma(11)
		sigma_cor[x2pnt(sigma_cor,2*i)] = W_sigma(3)
		sigma_z0[x2pnt(sigma_z0,((2*i)+1))] = W_sigma(0)
		sigma_A[x2pnt(sigma_A,((2*i)+1))] = W_sigma(2)
		sigma_x0[x2pnt(sigma_x0,((2*i)+1))] = W_sigma(6)
		sigma_xW[x2pnt(sigma_xW,((2*i)+1))] = W_sigma(10)
		sigma_y0[x2pnt(sigma_y0,((2*i)+1))] = W_sigma(8)
		sigma_yW[x2pnt(sigma_yW,((2*i)+1))] = W_sigma(12)
		sigma_cor[x2pnt(sigma_cor,((2*i)+1))] = W_sigma(4)
	endfor	

killwaves W_sigma, W_coef, M_WaveStats
end



//This function fits the dumbbell peak positions that are in x_loc and y_loc that were found in DumbbellPeakPositions(image).
//This only works for dumbbell 2D Gaussian fits that are aligned vertically.
//Inputs: image, x_loc, y_loc, and dumbbell separation.
function DumbbellGaussianFit(image, x_loc, y_loc, dum_sep)
	
	wave image, x_loc, y_loc
	variable dum_sep
	variable half_size = dum_sep / 2
	variable num_peaks = 2 * DimSize(x_loc,0)

	Duplicate/O image noise
	noise = sqrt(image)

	duplicate/O x_loc x_start
	duplicate/O x_loc x_finish
	duplicate/O y_loc y_start
	duplicate/O y_loc y_finish
	
	x_start = x_loc - half_size - 1
	x_finish = x_loc + half_size + 1
	y_start = y_loc - (dum_sep)
	y_finish = y_loc + (dum_sep)

	Make/O/N=(7) W_sigma
	Make/O/N=(num_peaks) z0, A, x0, xW, y0, yW, cor
	Make/O/N=(num_peaks) sigma_z0, sigma_A, sigma_x0, sigma_xW, sigma_y0, sigma_yW, sigma_cor

	//Define starting positions for fit.
	Make/O/N=27 M_WaveStats
	wavestats/Q/W image
	variable z0s = M_WaveStats(10), As = M_WaveStats(12), Bs = M_WaveStats(12), cor1s = 0.1, cor2s = 0.1
	variable xw1s = dum_sep / 3, xw2s = dum_sep / 3, yw1s = dum_sep / 3, yw2s = dum_sep / 3

	variable i
	for(i=0; i<num_peaks; i+=1)
		
		Make/D/N=13/O W_coef
		W_coef[0] = {z0s,As,Bs,cor1s,cor2s,x_loc(i),x_loc(i),(y_loc(i)-half_size),(y_loc(i)+half_size),xw1s,xw2s,yw1s,yw2s}
		FuncFitMD/N/NTHR=0/Q Gaus2Dx2 W_coef  image[x_start(i),x_finish(i)][y_start(i),y_finish(i)] /W=noise /I=1 /D/R
		
		z0[x2pnt(z0,2*i)] = W_coef(0)
		A[x2pnt(A,2*i)] = W_coef(1)
		x0[x2pnt(x0,2*i)] = W_coef(5)
		xW[x2pnt(xW,2*i)] = W_coef(9)
		y0[x2pnt(y0,2*i)] = W_coef(7)
		yW[x2pnt(yW,2*i)] = W_coef(11)
		cor[x2pnt(cor,2*i)] = W_coef(3)
		z0[x2pnt(z0,((2*i)+1))] = W_coef(0)
		A[x2pnt(A,((2*i)+1))] = W_coef(2)
		x0[x2pnt(x0,((2*i)+1))] = W_coef(6)
		xW[x2pnt(xW,((2*i)+1))] = W_coef(10)
		y0[x2pnt(y0,((2*i)+1))] = W_coef(8)
		yW[x2pnt(yW,((2*i)+1))] = W_coef(12)
		cor[x2pnt(cor,((2*i)+1))] = W_coef(4)
		sigma_z0[x2pnt(sigma_z0,2*i)] = W_sigma(0)
		sigma_A[x2pnt(sigma_A,2*i)] = W_sigma(1)
		sigma_x0[x2pnt(sigma_x0,2*i)] = W_sigma(5)
		sigma_xW[x2pnt(sigma_xW,2*i)] = W_sigma(9)
		sigma_y0[x2pnt(sigma_y0,2*i)] = W_sigma(7)
		sigma_yW[x2pnt(sigma_yW,2*i)] = W_sigma(11)
		sigma_cor[x2pnt(sigma_cor,2*i)] = W_sigma(3)
		sigma_z0[x2pnt(sigma_z0,((2*i)+1))] = W_sigma(0)
		sigma_A[x2pnt(sigma_A,((2*i)+1))] = W_sigma(2)
		sigma_x0[x2pnt(sigma_x0,((2*i)+1))] = W_sigma(6)
		sigma_xW[x2pnt(sigma_xW,((2*i)+1))] = W_sigma(10)
		sigma_y0[x2pnt(sigma_y0,((2*i)+1))] = W_sigma(8)
		sigma_yW[x2pnt(sigma_yW,((2*i)+1))] = W_sigma(12)
		sigma_cor[x2pnt(sigma_cor,((2*i)+1))] = W_sigma(4)
	endfor	

//appendtograph/t y0 vs x0
//ModifyGraph mode=2
//RemoveFromGraph/Z y_loc
killwaves W_sigma, W_coef
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

// Gaussian fit on every image in a stack, starting from the same initial guess positions without peak finding
function GaussianFitStack(st, x_loc, y_loc, size, wiggle [fix_xy, no_noise])
	wave st, x_loc, y_loc
	variable size, wiggle, fix_xy, no_noise
	
	string cur_fol = GetDataFolder(1)
	NewDataFolder/O/S root:Packages
	NewDataFolder/O/S GaussianFitStack
	
	variable natoms = numpnts(x_loc)
	variable nimages = DimSize(st, 2)
	
	Make/O/N=(natoms, nimages) z0_st, A_st, x0_st, xW_st, y0_st, yW_st, cor_st
	Make/O/N=(natoms, nimages) sig_z0_st, sig_A_st, sig_x0_st, sig_xW_st, sig_y0_st, sig_yW_st, sig_cor_st
	
	variable i
	for(i=0; i<nimages; i+=1)
		ImageTransform/P=(i) getplane st
		wave im = $"M_ImagePlane"
		
		if(fix_xy == 1)
			if(no_noise == 1)
				GaussianFit(im, x_loc, y_loc, size, wiggle, fix_xy = fix_xy, no_noise = no_noise)
			else
				GaussianFit(im, x_loc, y_loc, size, wiggle, fix_xy = fix_xy)
			endif
		else
			GaussianFit(im, x_loc, y_loc, size, wiggle)
		endif
		
		wave z0 = $"z0"
		wave A = $"A"
		wave x0 = $"x0"
		wave xW = $"xW"
		wave y0 = $"y0"
		wave yW = $"yW"
		wave cor = $"cor"
		wave sigma_z0 = $"sigma_z0"
		wave sigma_A = $"sigma_A"
		wave sigma_x0 = $"sigma_x0"
		wave sigma_xW = $"sigma_xW"
		wave sigma_y0 = $"sigma_y0"
		wave sigma_yW = $"sigma_yW"
		wave sigma_cor = $"sigma_cor"

		z0_st[][i] = z0[p] 
		A_st[][i] = A[p] 
		x0_st[][i] = x0[p]  
		xW_st[][i] = xW[p]  
		y0_st[][i] = y0[p]  
		yW_st[][i] = yW[p]  
		cor_st[][i] = cor[p] 
		sig_z0_st[][i] = sigma_z0[p] 
		sig_A_st[][i] = sigma_A[p] 
		sig_x0_st[][i] = sigma_x0[p]  
		sig_xW_st[][i] = sigma_xW[p]  
		sig_y0_st[][i] = sigma_y0[p]  
		sig_yW_st[][i] = sigma_yW[p]  
		sig_cor_st[][i] = sigma_cor[p]  
	
	endfor
	
	killwaves z0, A, x0, xW, y0, yW, cor, sigma_z0, sigma_A, sigma_x0, sigma_xW, sigma_y0, sigma_yW, sigma_cor	
	
	Duplicate/O z0_st $(cur_fol+"z0_st")
	Duplicate/O A_st $(cur_fol+"A_st")
	Duplicate/O x0_st $(cur_fol+"x0_st")
	Duplicate/O xW_st $(cur_fol+"xW_st")
	Duplicate/O y0_st $(cur_fol+"y0_st")
	Duplicate/O yW_st $(cur_fol+"yW_st")
	Duplicate/O cor_st $(cur_fol+"cor_st")
	Duplicate/O sig_z0_st $(cur_fol+"sig_z0_st")
	Duplicate/O sig_A_st $(cur_fol+"sig_A_st")
	Duplicate/O sig_x0_st $(cur_fol+"sig_x0_st")
	Duplicate/O sig_xW_st $(cur_fol+"sig_xW_st")
	Duplicate/O sig_y0_st $(cur_fol+"sig_y0_st")
	Duplicate/O sig_yW_st $(cur_fol+"sig_yW_st")
	Duplicate/O sig_cor_st $(cur_fol+"sig_cor_st")
	
	SetDataFolder $cur_fol
	
end


// Does the peak find, 2dGaus fit, and x and y precision calculations on 1 image of dumbbells.
// Outputs sep_out with x and y avergae and stdev.
function DumbbellPrecision(image, dum_sep, x_space, y_space, space_delta, x_angle, y_angle, angle_delta)

	wave image
	variable dum_sep, x_space, y_space, space_delta, x_angle, y_angle, angle_delta
	
	DumbbellPeakPositions(image)
	
	wave x_loc =$"x_loc"
	wave y_loc =$"y_loc"
	if(!WaveExists(x_loc) || !WaveExists(y_loc))
		print "PeakPositions failed to create x_loc and/or y_loc.  Exiting.\r"
		return 0
	endif
	
	DumbbellGaussianFit(image, x_loc, y_loc, dum_sep)
	
	Separation(x_space, y_space, space_delta, x_angle, y_angle, angle_delta)	
end



// Does the peak find, 2dGaus fit, and x and y precision calculation on a stack of images.
// Outputs xav, yav, xstd, ystd, and fit parameters for all column position in each image in stack.
function StackPrecision(image_stack, size, x_space, y_space, space_delta, x_angle, y_angle, angle_delta)

	wave image_stack
	variable size, x_space, y_space, space_delta, x_angle, y_angle, angle_delta
	
	Make/O/N=(DimSize(image_stack, 2)) xav, yav, xstd, ystd

	Display  xav
	ModifyGraph rgb(xav)=(16384,48896,65280);DelayUpdate
	ErrorBars/T=0 xav Y,wave=(xstd,xstd)
	AppendToGraph xav
	ModifyGraph rgb(xav#1)=(0,0,0);DelayUpdate
	Label left "Average X Seperation (pixels)";DelayUpdate
	Label bottom "Frame"

	Display  xstd
	ModifyGraph rgb(xstd)=(0,0,0);DelayUpdate
	Label left "X Stdev (pixels)";DelayUpdate
	Label bottom "Frame"

	Display  yav
	ModifyGraph rgb(yav)=(16384,48896,65280);DelayUpdate
	ErrorBars/T=0 yav Y,wave=(ystd,ystd)
	AppendToGraph yav
	ModifyGraph rgb(yav#1)=(0,0,0)
	Label left "Average Y Seperation (pixels)";DelayUpdate
	Label bottom "Frame"

	Display  ystd
	ModifyGraph rgb(ystd)=(0,0,0);DelayUpdate
	Label left "Y Stdev (pixels)";DelayUpdate
	Label bottom "Frame"
	
	variable i
	for(i=0; i<DimSize(image_stack, 2); i+=1)
	//for(i=0; i<4; i+=1)
		Imagetransform/p=(i) getplane image_stack
		wave im = $"M_ImagePlane"
		Precision(im, size, x_space, y_space, space_delta, x_angle, y_angle, angle_delta)
	
		//save averages and stdevs
		wave prec = $"prec"
		xav[i] = prec[0]
		xstd[i] = prec[1]
		yav[i] = prec[2]
		ystd[i] = prec[3]
		
		//save fit parameters for each atom column in each frame of stack
		wave z0 = $"z0"
		wave A = $"A"
		wave x0 = $"x0"
		wave xW = $"xW"
		wave y0 = $"y0"
		wave yW = $"yW"
		wave cor = $"cor"
		wave sigma_z0 = $"sigma_z0"
		wave sigma_A = $"sigma_A"
		wave sigma_x0 = $"sigma_x0"
		wave sigma_xW = $"sigma_xW"
		wave sigma_y0 = $"sigma_y0"
		wave sigma_yW = $"sigma_yW"
		wave sigma_cor = $"sigma_cor"
		variable num_peaks = DimSize(z0,0)
		if (i == 0)
			Make/O/N=(num_peaks, DimSize(image_stack, 2)) z0_stack, A_stack, x0_stack, xW_stack, y0_stack, yW_stack, cor_stack 
			Make/O/N=(num_peaks, DimSize(image_stack, 2)) sig_z0_stack, sig_A_stack, sig_x0_stack, sig_xW_stack, sig_y0_stack, sig_yW_stack, sig_cor_stack 
		endif
		z0_stack[][i] = z0[p] 
		A_stack[][i] = A[p] 
		x0_stack[][i] = x0[p]  
		xW_stack[][i] = xW[p]  
		y0_stack[][i] = y0[p]  
		yW_stack[][i] = yW[p]  
		cor_stack[][i] = cor[p] 
		sig_z0_stack[][i] = sigma_z0[p] 
		sig_A_stack[][i] = sigma_A[p] 
		sig_x0_stack[][i] = sigma_x0[p]  
		sig_xW_stack[][i] = sigma_xW[p]  
		sig_y0_stack[][i] = sigma_y0[p]  
		sig_yW_stack[][i] = sigma_yW[p]  
		sig_cor_stack[][i] = sigma_cor[p]  
	
		killwaves z0, A, x0, xW, y0, yW, cor, sigma_z0, sigma_A, sigma_x0, sigma_xW, sigma_y0, sigma_yW, sigma_cor	
	endfor	


end



// Does the dumbbell peak find, dumbbell Gaus fit, and x and y precision calculation on a stack of images.
// Outputs xav, yav, xstd, ystd, and fit parameters for all column position in each image in stack.
function DumbbellStackPrecision(image_stack, dum_sep, x_space, y_space, space_delta, x_angle, y_angle, angle_delta)

	wave image_stack
	variable dum_sep, x_space, y_space, space_delta, x_angle, y_angle, angle_delta
	
	Make/O/N=(DimSize(image_stack, 2)) xav, yav, xstd, ystd
	
	Display  xav
	ModifyGraph rgb(xav)=(16384,48896,65280);DelayUpdate
	ErrorBars/T=0 xav Y,wave=(xstd,xstd)
	AppendToGraph xav
	ModifyGraph rgb(xav#1)=(0,0,0);DelayUpdate
	Label left "Average X Seperation (pixels)";DelayUpdate
	Label bottom "Frame"

	Display  xstd
	ModifyGraph rgb(xstd)=(0,0,0);DelayUpdate
	Label left "X Stdev (pixels)";DelayUpdate
	Label bottom "Frame"

	Display  yav
	ModifyGraph rgb(yav)=(16384,48896,65280);DelayUpdate
	ErrorBars/T=0 yav Y,wave=(ystd,ystd)
	AppendToGraph yav
	ModifyGraph rgb(yav#1)=(0,0,0)
	Label left "Average Y Seperation (pixels)";DelayUpdate
	Label bottom "Frame"

	Display  ystd
	ModifyGraph rgb(ystd)=(0,0,0);DelayUpdate
	Label left "Y Stdev (pixels)";DelayUpdate
	Label bottom "Frame"
	
	variable i
	for(i=0; i<DimSize(image_stack, 2); i+=1)
	//for(i=0; i<4; i+=1)
		Imagetransform/p=(i) getplane image_stack
		wave im = $"M_ImagePlane"
		
		DumbbellPrecision(im, dum_sep, x_space, y_space, space_delta, x_angle, y_angle, angle_delta)
	
		//save averages and stdevs
		wave prec = $"prec"
		xav[i] = prec[0]
		xstd[i] = prec[1]
		yav[i] = prec[2]
		ystd[i] = prec[3]
		
		//save fit parameters for each atom column in each frame of stack
		wave z0 = $"z0"
		wave A = $"A"
		wave x0 = $"x0"
		wave xW = $"xW"
		wave y0 = $"y0"
		wave yW = $"yW"
		wave cor = $"cor"
		wave sigma_z0 = $"sigma_z0"
		wave sigma_A = $"sigma_A"
		wave sigma_x0 = $"sigma_x0"
		wave sigma_xW = $"sigma_xW"
		wave sigma_y0 = $"sigma_y0"
		wave sigma_yW = $"sigma_yW"
		wave sigma_cor = $"sigma_cor"
		variable num_peaks = DimSize(z0,0)
		if (i == 0)
			Make/O/N=(num_peaks, DimSize(image_stack, 2)) z0_stack, A_stack, x0_stack, xW_stack, y0_stack, yW_stack, cor_stack 
			Make/O/N=(num_peaks, DimSize(image_stack, 2)) sig_z0_stack, sig_A_stack, sig_x0_stack, sig_xW_stack, sig_y0_stack, sig_yW_stack, sig_cor_stack 
		endif
		z0_stack[][i] = z0[p] 
		A_stack[][i] = A[p] 
		x0_stack[][i] = x0[p]  
		xW_stack[][i] = xW[p]  
		y0_stack[][i] = y0[p]  
		yW_stack[][i] = yW[p]  
		cor_stack[][i] = cor[p] 
		sig_z0_stack[][i] = sigma_z0[p] 
		sig_A_stack[][i] = sigma_A[p] 
		sig_x0_stack[][i] = sigma_x0[p]  
		sig_xW_stack[][i] = sigma_xW[p]  
		sig_y0_stack[][i] = sigma_y0[p]  
		sig_yW_stack[][i] = sigma_yW[p]  
		sig_cor_stack[][i] = sigma_cor[p]  
	
		killwaves z0, A, x0, xW, y0, yW, cor, sigma_z0, sigma_A, sigma_x0, sigma_xW, sigma_y0, sigma_yW, sigma_cor	
	endfor
end



// Does the peak find, and 2dGaus fit on one image of single atom columns.
// Outputs Gaussian Fit parameters for each atom column position.
function SimulationPosition(image, size)

	wave image
	variable size
	
	PeakPositions(image)
	
	wave x_loc =$"x_loc"
	wave y_loc =$"y_loc"
	if(!WaveExists(x_loc) || !WaveExists(y_loc))
		print "PeakPositions failed to create x_loc and/or y_loc.  Exiting.\r"
		return 0
	endif
	
	GaussianFit(image, x_loc, y_loc, size, 0)	// hard coded for zero fitting box wiggle
end



// Does the peak find and 2dGaus fit on a stack of images.
// Outputs Gaus fit parameters for each atom position in each image in stack.
function SimulationStackPosition(image_stack, size)

	wave image_stack
	variable size
	
	variable i
	for(i=1; i<DimSize(image_stack, 2); i+=1)
	//for(i=0; i<4; i+=1)
		Imagetransform/p=(i) getplane image_stack
		wave im = $"M_ImagePlane"
		
		SimulationPosition(im, size)
	
		//save fit parameters for each atom column in each frame of stack
		wave z0 = $"z0"
		wave A = $"A"
		wave x0 = $"x0"
		wave xW = $"xW"
		wave y0 = $"y0"
		wave yW = $"yW"
		wave cor = $"cor"
		wave sigma_z0 = $"sigma_z0"
		wave sigma_A = $"sigma_A"
		wave sigma_x0 = $"sigma_x0"
		wave sigma_xW = $"sigma_xW"
		wave sigma_y0 = $"sigma_y0"
		wave sigma_yW = $"sigma_yW"
		wave sigma_cor = $"sigma_cor"
		variable num_peaks = DimSize(z0,0)
		if (i == 1)
			Make/O/N=(num_peaks, DimSize(image_stack, 2)-1) z0_stack, A_stack, x0_stack, xW_stack, y0_stack, yW_stack, cor_stack 
			Make/O/N=(num_peaks, DimSize(image_stack, 2)-1) sig_z0_stack, sig_A_stack, sig_x0_stack, sig_xW_stack, sig_y0_stack, sig_yW_stack, sig_cor_stack 
		endif
		z0_stack[][i-1] = z0[p] 
		A_stack[][i-1] = A[p] 
		x0_stack[][i-1] = x0[p]  
		xW_stack[][i-1] = xW[p]  
		y0_stack[][i-1] = y0[p]  
		yW_stack[][i-1] = yW[p]  
		cor_stack[][i-1] = cor[p] 
		sig_z0_stack[][i-1] = sigma_z0[p] 
		sig_A_stack[][i-1] = sigma_A[p] 
		sig_x0_stack[][i-1] = sigma_x0[p]  
		sig_xW_stack[][i-1] = sigma_xW[p]  
		sig_y0_stack[][i-1] = sigma_y0[p]  
		sig_yW_stack[][i-1] = sigma_yW[p]  
		sig_cor_stack[][i-1] = sigma_cor[p]  
	
		killwaves z0, A, x0, xW, y0, yW, cor, sigma_z0, sigma_A, sigma_x0, sigma_xW, sigma_y0, sigma_yW, sigma_cor	
	endfor	
	
end



// Does the peak find and 2dGaus fit on a stack of images.
// Outputs Gaus fit parameters for each atom position in each image in stack.
function HAADFSimStackPosition(image_stack, size)

	wave image_stack
	variable size
	
	//Find initial x_loc and y_loc to create a fitting box window from
	Imagetransform/p=(DimSize(image_stack, 2)-1) getplane image_stack
	wave im_0 = $"M_ImagePlane"
	PeakPositions(im_0)
	wave x_loc =$"x_loc"
	wave y_loc =$"y_loc"
	if(!WaveExists(x_loc) || !WaveExists(y_loc))
		print "PeakPositions failed to create x_loc and/or y_loc.  Exiting.\r"
		return 0
	endif
	
	//Make/O/N=(DimSize(image_stack, 2)) xav, yav, xstd, ystd
	
	variable i
	for(i=1; i<DimSize(image_stack, 2); i+=1)
		Imagetransform/p=(i) getplane image_stack
		wave im = $"M_ImagePlane"
		
		//This is to fit using gaussian fit where it inputs all the fit parameters initial guess with the initial positions all the same.
		GaussianFitInitialGuess(im, x_loc, y_loc, size)
		
		//This is to fit using gaussian fot where it inputs a box and not all the fit parameters initial guess
		//GaussianFit(im, x_loc, y_loc, size)
		
		//save fit parameters for each atom column in each frame of stack
		wave Res_M_ImagePlane = $"Res_M_ImagePlane"
		wave fit_M_ImagePlane = $"fit_M_ImagePlane"
		wave z0 = $"z0"
		wave A = $"A"
		wave x0 = $"x0"
		wave xW = $"xW"
		wave y0 = $"y0"
		wave yW = $"yW"
		wave cor = $"cor"
		wave sigma_z0 = $"sigma_z0"
		wave sigma_A = $"sigma_A"
		wave sigma_x0 = $"sigma_x0"
		wave sigma_xW = $"sigma_xW"
		wave sigma_y0 = $"sigma_y0"
		wave sigma_yW = $"sigma_yW"
		wave sigma_cor = $"sigma_cor"
		variable num_peaks = DimSize(z0,0)
		if (i == 1)
			Make/O/N=(num_peaks, DimSize(image_stack, 2)-1) z0_stack, A_stack, x0_stack, xW_stack, y0_stack, yW_stack, cor_stack 
			Make/O/N=(num_peaks, DimSize(image_stack, 2)-1) sig_z0_stack, sig_A_stack, sig_x0_stack, sig_xW_stack, sig_y0_stack, sig_yW_stack, sig_cor_stack 
			Duplicate/O image_stack res_stack	, fit_stack
		endif
		z0_stack[][i-1] = z0[p] 
		A_stack[][i-1] = A[p] 
		x0_stack[][i-1] = x0[p]  
		xW_stack[][i-1] = xW[p]  
		y0_stack[][i-1] = y0[p]  
		yW_stack[][i-1] = yW[p]  
		cor_stack[][i-1] = cor[p] 
		sig_z0_stack[][i-1] = sigma_z0[p] 
		sig_A_stack[][i-1] = sigma_A[p] 
		sig_x0_stack[][i-1] = sigma_x0[p]  
		sig_xW_stack[][i-1] = sigma_xW[p]  
		sig_y0_stack[][i-1] = sigma_y0[p]  
		sig_yW_stack[][i-1] = sigma_yW[p]  
		sig_cor_stack[][i-1] = sigma_cor[p]  
		
		res_stack[][][i-1] = Res_M_ImagePlane[p][q] 
		fit_stack[][][i-1] = fit_M_ImagePlane[p][q] 
		
		killwaves z0, A, x0, xW, y0, yW, cor, sigma_z0, sigma_A, sigma_x0, sigma_xW, sigma_y0, sigma_yW, sigma_cor	
	endfor	
killwaves M_ImagePlane, Res_M_ImagePlane, fit_M_ImagePlane
end



// Does the peak find and 2dGaus fit on a stack of images.
// Outputs Gaus fit parameters for each atom position in each image in stack.
function HAADFDumSimStackPosition(image_stack, size)

	wave image_stack
	variable size
	
	//Find initial x_loc and y_loc to create a fitting box window from
	Imagetransform/p=(DimSize(image_stack, 2)-1) getplane image_stack
	wave im_0 = $"M_ImagePlane"
	PeakPositions(im_0)
	wave x_loc =$"x_loc"
	wave y_loc =$"y_loc"
	if(!WaveExists(x_loc) || !WaveExists(y_loc))
		print "PeakPositions failed to create x_loc and/or y_loc.  Exiting.\r"
		return 0
	endif
	
	//Make/O/N=(DimSize(image_stack, 2)) xav, yav, xstd, ystd
	
	variable i
	for(i=1; i<DimSize(image_stack, 2); i+=1)
	//for(i=1; i<2; i+=1)
		Imagetransform/p=(i) getplane image_stack
		wave im = $"M_ImagePlane"
		
		HAADFDumbbellGaussianFit(im, x_loc, y_loc, size)
		
		//save fit parameters for each atom column in each frame of stack
		wave Res_M_ImagePlane = $"Res_M_ImagePlane"
		wave fit_M_ImagePlane = $"fit_M_ImagePlane"
		wave z0 = $"z0"
		wave A = $"A"
		wave x0 = $"x0"
		wave xW = $"xW"
		wave y0 = $"y0"
		wave yW = $"yW"
		wave cor = $"cor"
		wave sigma_z0 = $"sigma_z0"
		wave sigma_A = $"sigma_A"
		wave sigma_x0 = $"sigma_x0"
		wave sigma_xW = $"sigma_xW"
		wave sigma_y0 = $"sigma_y0"
		wave sigma_yW = $"sigma_yW"
		wave sigma_cor = $"sigma_cor"
		variable num_peaks = DimSize(z0,0)
		if (i == 1)
			Make/O/N=(num_peaks, DimSize(image_stack, 2)-1) z0_stack, A_stack, x0_stack, xW_stack, y0_stack, yW_stack, cor_stack 
			Make/O/N=(num_peaks, DimSize(image_stack, 2)-1) sig_z0_stack, sig_A_stack, sig_x0_stack, sig_xW_stack, sig_y0_stack, sig_yW_stack, sig_cor_stack 
			Duplicate/O image_stack res_stack	, fit_stack
		endif
		z0_stack[][i-1] = z0[p] 
		A_stack[][i-1] = A[p] 
		x0_stack[][i-1] = x0[p]  
		xW_stack[][i-1] = xW[p]  
		y0_stack[][i-1] = y0[p]  
		yW_stack[][i-1] = yW[p]  
		cor_stack[][i-1] = cor[p] 
		sig_z0_stack[][i-1] = sigma_z0[p] 
		sig_A_stack[][i-1] = sigma_A[p] 
		sig_x0_stack[][i-1] = sigma_x0[p]  
		sig_xW_stack[][i-1] = sigma_xW[p]  
		sig_y0_stack[][i-1] = sigma_y0[p]  
		sig_yW_stack[][i-1] = sigma_yW[p]  
		sig_cor_stack[][i-1] = sigma_cor[p]  
		
		res_stack[][][i-1] = Res_M_ImagePlane[p][q] 
		fit_stack[][][i-1] = fit_M_ImagePlane[p][q] 
		
		killwaves z0, A, x0, xW, y0, yW, cor, sigma_z0, sigma_A, sigma_x0, sigma_xW, sigma_y0, sigma_yW, sigma_cor	
	endfor	
killwaves M_ImagePlane, Res_M_ImagePlane, fit_M_ImagePlane
end



function ABFDumSimStackPositionPixel(image_stack, size)

	wave image_stack
	variable size
	
	//Find initial x_loc and y_loc to create a fitting box window from
	Imagetransform/p=(DimSize(image_stack, 2)-1) getplane image_stack
	wave im_0 = $"M_ImagePlane"
	ABFPeakPositions(im_0)
	wave x_loc =$"x_loc"
	wave y_loc =$"y_loc"
	if(!WaveExists(x_loc) || !WaveExists(y_loc))
		print "PeakPositions failed to create x_loc and/or y_loc.  Exiting.\r"
		return 0
	endif
	
	variable i
	for(i=1; i<DimSize(image_stack, 2); i+=1)
	//for(i=1; i<2; i+=1)
		Imagetransform/p=(i) getplane image_stack
		wave im = $"M_ImagePlane"
		
		ABFDumbbellGaussianFit(im, x_loc, y_loc, size)
		
		//save fit parameters for each atom column in each frame of stack
		wave Res_M_ImagePlane = $"Res_M_ImagePlane"
		wave fit_M_ImagePlane = $"fit_M_ImagePlane"
		wave z0 = $"z0"
		wave A = $"A"
		wave x0 = $"x0"
		wave xW = $"xW"
		wave y0 = $"y0"
		wave yW = $"yW"
		wave cor = $"cor"
		wave sigma_z0 = $"sigma_z0"
		wave sigma_A = $"sigma_A"
		wave sigma_x0 = $"sigma_x0"
		wave sigma_xW = $"sigma_xW"
		wave sigma_y0 = $"sigma_y0"
		wave sigma_yW = $"sigma_yW"
		wave sigma_cor = $"sigma_cor"
		variable num_peaks = DimSize(z0,0)
		if (i == 1)
			Make/O/N=(num_peaks, DimSize(image_stack, 2)-1) z0_stack, A_stack, x0_stack, xW_stack, y0_stack, yW_stack, cor_stack 
			Make/O/N=(num_peaks, DimSize(image_stack, 2)-1) sig_z0_stack, sig_A_stack, sig_x0_stack, sig_xW_stack, sig_y0_stack, sig_yW_stack, sig_cor_stack 
			Duplicate/O image_stack res_stack	, fit_stack
		endif
		z0_stack[][i-1] = z0[p] 
		A_stack[][i-1] = A[p] 
		x0_stack[][i-1] = x0[p]  
		xW_stack[][i-1] = xW[p]  
		y0_stack[][i-1] = y0[p]  
		yW_stack[][i-1] = yW[p]  
		cor_stack[][i-1] = cor[p] 
		sig_z0_stack[][i-1] = sigma_z0[p] 
		sig_A_stack[][i-1] = sigma_A[p] 
		sig_x0_stack[][i-1] = sigma_x0[p]  
		sig_xW_stack[][i-1] = sigma_xW[p]  
		sig_y0_stack[][i-1] = sigma_y0[p]  
		sig_yW_stack[][i-1] = sigma_yW[p]  
		sig_cor_stack[][i-1] = sigma_cor[p]  
		
		res_stack[][][i-1] = Res_M_ImagePlane[p][q] 
		fit_stack[][][i-1] = fit_M_ImagePlane[p][q] 
		
		killwaves z0, A, x0, xW, y0, yW, cor, sigma_z0, sigma_A, sigma_x0, sigma_xW, sigma_y0, sigma_yW, sigma_cor	
	endfor	
killwaves M_ImagePlane, Res_M_ImagePlane, fit_M_ImagePlane
end



function ABFDumSimStackPositionAng(image_stack, x_loc, y_loc, size)

	wave image_stack, x_loc, y_loc
	variable size
	
	variable i
	for(i=1; i<DimSize(image_stack, 2); i+=1)
		Imagetransform/p=(i) getplane image_stack
		wave im = $"M_ImagePlane"
		
		ABFDumbbellGaussianFitAng(im, x_loc, y_loc, size)
		
		//save fit parameters for each atom column in each frame of stack
		wave z0 = $"z0"
		wave A = $"A"
		wave x0 = $"x0"
		wave xW = $"xW"
		wave y0 = $"y0"
		wave yW = $"yW"
		wave cor = $"cor"
		wave sigma_z0 = $"sigma_z0"
		wave sigma_A = $"sigma_A"
		wave sigma_x0 = $"sigma_x0"
		wave sigma_xW = $"sigma_xW"
		wave sigma_y0 = $"sigma_y0"
		wave sigma_yW = $"sigma_yW"
		wave sigma_cor = $"sigma_cor"
		variable num_peaks = DimSize(z0,0)
		if (i == 1)
			Make/O/N=(num_peaks, DimSize(image_stack, 2)-1) z0_stack, A_stack, x0_stack, xW_stack, y0_stack, yW_stack, cor_stack 
			Make/O/N=(num_peaks, DimSize(image_stack, 2)-1) sig_z0_stack, sig_A_stack, sig_x0_stack, sig_xW_stack, sig_y0_stack, sig_yW_stack, sig_cor_stack 
		endif
		z0_stack[][i-1] = z0[p] 
		A_stack[][i-1] = A[p] 
		x0_stack[][i-1] = x0[p]  
		xW_stack[][i-1] = xW[p]  
		y0_stack[][i-1] = y0[p]  
		yW_stack[][i-1] = yW[p]  
		cor_stack[][i-1] = cor[p] 
		sig_z0_stack[][i-1] = sigma_z0[p] 
		sig_A_stack[][i-1] = sigma_A[p] 
		sig_x0_stack[][i-1] = sigma_x0[p]  
		sig_xW_stack[][i-1] = sigma_xW[p]  
		sig_y0_stack[][i-1] = sigma_y0[p]  
		sig_yW_stack[][i-1] = sigma_yW[p]  
		sig_cor_stack[][i-1] = sigma_cor[p]  
		
		killwaves z0, A, x0, xW, y0, yW, cor, sigma_z0, sigma_A, sigma_x0, sigma_xW, sigma_y0, sigma_yW, sigma_cor	
	endfor	
killwaves M_ImagePlane
end



// Claculates contrast as a function of stack image
function StackContrast(image_stack)

	wave image_stack
	
	Make/O/N=(DimSize(image_stack, 2)) minimum, maximum, contrast
	Make/O/N=27 M_WaveStats
	
	variable i
	for(i=0; i<DimSize(image_stack, 2); i+=1)
	//for(i=0; i<4; i+=1)
		Imagetransform/p=(i) getplane image_stack
		wave im = $"M_ImagePlane"
		
		wavestats/Q/W  im
		
		variable minim = M_WaveStats(10)
		variable maxim = M_WaveStats(12)
		
		
		minimum[i] = minim
		maximum[i] = maxim
		contrast[i] = (maxim - minim) / minim
		
	endfor
	
	Display  contrast
	Label left "Contrast";DelayUpdate
	Label bottom "Frame"
	Display  maximum
	Label left "Maximum";DelayUpdate
	Label bottom "Frame"
	Display  minimum
	Label left "Minimum";DelayUpdate
	Label bottom "Frame"

killwaves M_WaveStats
end



//Creates an average image from a floating window within a stack, "image_stack" with width "float_width".
function FloatingWindow(image_stack, float_width)

	wave image_stack
	variable float_width
	
	Make/O/N=(DimSize(image_stack, 0), DimSize(image_stack, 1), DimSize(image_stack, 2)-float_width) float_stack
	Make/O/N=(DimSize(image_stack, 0), DimSize(image_stack, 1), float_width) float_window
	
	variable i, j
	for(i=0; i<DimSize(image_stack, 2)-float_width; i+=1)
		
		for(j=0; j<float_width; j+=1)
			Imagetransform/p=(i+j) getplane image_stack
			wave PlaneA = $"M_ImagePlane"
		
			float_window[][][j] = PlaneA[p][q]
		endfor
		
		imagetransform sumPlanes float_window
		wave PlaneB = $"M_SumPlanes"	
		float_stack[][][i] = PlaneB[p][q]/float_width
	
	endfor

killwaves M_ImagePlane, M_SumPlanes, float_window
end


//Calculates the summed intensity in each frame of a series and outputs it in Total_Intensity
function SumIntensity(image_stack)

	wave image_stack
	variable numFrames=DimSize(image_stack, 2)
	
	Make/O/N=(numFrames) Sum_Intensity
	
	Make/O/N=27 M_WaveStats
	
	variable i, j
	for(i=0; i<DimSize(image_stack, 2); i+=1)
		
		Imagetransform/p=(i) getplane image_stack
		wave PlaneA = $"M_ImagePlane"
		
		wavestats/Q/W PlaneA
		
		Sum_Intensity[i] = M_WaveStats(23)
	
	endfor

	killwaves M_ImagePlane, M_WaveStats
end

// takes an image series (stack) and calculates the mean of each image.
// Places the results in stack_mean.
function StackMean(st)
	wave st
	
	make/o/n=(DimSize(st, 2)) stack_mean
	
	variable i
	for(i=0; i<DimSize(st, 2); i+=1)
		Imagetransform/P=(i) getplane st
		wave pl = $"M_ImagePlane"
		wavestats/q/m=1 pl
		stack_mean[i] = V_avg
	endfor
	
	Killwaves pl
	
end