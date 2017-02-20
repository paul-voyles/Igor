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


// This function converts images in counts to electrons
// inputs: original image in counts, numSamples, and HAADF gain in electrons/counts

function stackJDDCPeakPositions(imstack,startpoint,threshold)

	variable threshold
	wave imstack, startpoint
	variable numframes = DimSize(imstack,2)
	variable i
	Make/O/N=(numframes,2) PrecList
	
	for (i=0; i<numframes; i++)
		Imagetransform/p=(i) getplane imstack
		wave im = $"M_ImagePlane"
		JDDCPeakPositions(im,startpoint[i][0],startpoint[i][1],threshold)
		wave prec = $"prec"
		PrecList[i][0]=prec[1]*21.16
		PrecList[i][1]=prec[3]*21.19
	endfor
	
end



function JDDCPeakPositions(average,startX,startY,threshold)

	variable startX, startY, threshold
	wave average
	duplicate/O average average_new //backup the original image
	
	//crop image according to start point and image size
	variable endX=startX + 149
	variable endY=startY + 149
	cropimage(average_new, startX, startY, endX, endY)
	
	//detect peaks with thresholdpeakpositions then filter peaks too close to border
	ThresholdPeakPositions(average_new, threshold)
	wave x_loc = $"x_loc"
	wave y_loc = $"y_loc"
	// Consider how to use SortColumn operate to make it work better
	variable numpeaks = DimSize(x_loc,0)
	variable i
	for (i=0; i<numpeaks ; i++)
		if(x_loc(i)<17 || x_loc(i)>135 || y_loc(i)<20 || y_loc(i)>135)
			deletepoints i, 1, x_loc
			deletepoints i, 1, y_loc
			i = i -1;
			numpeaks = numpeaks - 1;
		endif
	endfor
	
	//Typical gaussian fit and separation anlysis
	Gaussianfit(average_new,x_loc,y_loc,14,1)
	wave x0 = $"x0"
	wave y0 = $"y0"
	separation(19,19,3,0,-90,10)
	wave prec = $"prec"
	variable x_prec = prec[1]
	variable y_prec = prec[3]
	x_prec = x_prec * 21.16
	y_prec = y_prec * 21.19
	printf "x precision is %g pm\n", x_prec
	printf "y precision is %g pm\n", y_prec 
	
	//Display results	
	newimage average_new
	appendtograph/t y0 vs x0
	modifygraph mode=2
	modifygraph lsize=3
end

function CountstoElectrons(average, numSamples, gain)
	
	variable gain
	
	wave average 
	wave numSamples 

	Duplicate/O average average_withnum
	average_withnum = average * numSamples
	Duplicate/O average average_electrons
	average_electrons = average_withnum * gain
end

// This function converts a whole image series from counts to number of electrons
// gain should be a number in the units of electrons/counts, which can be cauculated with current, expo time and dark/bright level
function ConvertElectrons(img_stack,gain)

	variable gain
	wave img_stack
	
	Duplicate/O img_stack img_stack_e
	img_stack_e = img_stack * gain
	
end


function TiPeakPositions(image,threshold_1, threshold_2)
	
	variable threshold_1, threshold_2
	wave image
	
	//NewImage image
	
	ImageThreshold/I/T=(threshold_1)/M=0/Q image
	ImageAnalyzeParticles /E/W/Q/F/M=3/A=2/EBPC stats, root:M_ImageThresh
	
	duplicate/O W_xmin x_loc_1	
	duplicate/O W_ymin y_loc_1
	duplicate/O W_xmin xmin	
	duplicate/O W_ymin ymin
	duplicate/O W_xmax xmax	
	duplicate/O W_ymax ymax
	
	x_loc_1 = (xmax + xmin) / 2
	y_loc_1 = (ymax + ymin) / 2
	
	killwaves W_ImageObjArea, W_SpotX, W_SpotY, W_circularity, W_rectangularity, W_ImageObjPerimeter, M_Moments, M_RawMoments
	killwaves W_BoundaryX, W_BoundaryY, W_BoundaryIndex, W_xmin, W_xmax, W_ymin, W_ymax, xmin, xmax, ymin, ymax
	
	ImageThreshold/I/T=(threshold_2)/M=0/Q image
	ImageAnalyzeParticles /E/W/Q/F/M=3/A=2/EBPC stats, root:M_ImageThresh
	
	duplicate/O W_xmin x_loc_2	
	duplicate/O W_ymin y_loc_2
	duplicate/O W_xmin xmin	
	duplicate/O W_ymin ymin
	duplicate/O W_xmax xmax	
	duplicate/O W_ymax ymax
	
	x_loc_2 = (xmax + xmin) / 2
	y_loc_2 = (ymax + ymin) / 2
	
	killwaves W_ImageObjArea, W_SpotX, W_SpotY, W_circularity, W_rectangularity, W_ImageObjPerimeter, M_Moments, M_RawMoments
	killwaves W_BoundaryX, W_BoundaryY, W_BoundaryIndex, W_xmin, W_xmax, W_ymin, W_ymax, xmin, xmax, ymin, ymax
	
	variable num_peaks_1, num_peaks_2
	num_peaks_1 = DimSize(x_loc_1, 0)
	num_peaks_2 = DimSize(x_loc_2, 0)
	//printf "%g %g\r" num_peaks_1, num_peaks_2
	
	variable j,k
	variable i=1
	Make/O/N=0 x_loc, y_loc
	for(j=0; j<num_peaks_2; j+=1)
		for(k=0; k<num_peaks_1; k+=1)
			//printf "%g %g\r" abs(x_loc_2(j) - x_loc_1(k)), abs(y_loc_2(j) - y_loc_1(k)) 
			if(abs(x_loc_2(j) - x_loc_1(k)) < 2 && abs(y_loc_2(j) - y_loc_1(k)) < 2)
				i=0
			endif
		endfor
		//printf "%g\r" i
		if(i==1)
			InsertPoints 0,1,x_loc
			x_loc[0]=x_loc_2(j)
			InsertPoints 0,1,y_loc
			y_loc[0]=y_loc_2(j)
		endif
		i=1
	endfor
	appendtograph/t y_loc vs x_loc
	modifygraph mode=2
	modifygraph lsize=3
	killwaves x_loc_1, y_loc_1, x_loc_2, y_loc_2,M_ImageThresh, M_Particle
	
end


function ThresholdPeakPositions(image, threshold)

	wave image
	variable threshold
	
	//NewImage image
	
	ImageThreshold/I/T=(threshold)/M=0/Q image
	ImageAnalyzeParticles /E/W/Q/F/M=3/A=2/EBPC stats, root:M_ImageThresh
	
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
	ModifyGraph lsize=3
	
	killwaves W_ImageObjArea, W_SpotX, W_SpotY, W_circularity, W_rectangularity, W_ImageObjPerimeter, M_Moments, M_RawMoments
	killwaves W_BoundaryX, W_BoundaryY, W_BoundaryIndex, W_xmin, W_xmax, W_ymin, W_ymax, xmin, xmax, ymin, ymax
end

// This function finds the peak positions using the analyze particle IGOR function and reports the locations in x_loc and y_loc.
// For Si dumbells, it reports the center of the dumbells. To change the mimimun area of the particles change the A flag in the ImageAnalyzeParticle function.
// Inputs: image
function PeakPositions(image)
	
	wave image
	
	//NewImage image
	
	ImageThreshold/I/T=4600/M=0/Q image //method=2 image histogrtam is a simple bimodal distribution, can only pick Sr atoms from STO image
	ImageAnalyzeParticles /E/W/Q/F/M=3/A=2/EBPC stats, root:M_ImageThresh
	
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
	ModifyGraph lsize=2
	
	killwaves W_ImageObjArea, W_SpotX, W_SpotY, W_circularity, W_rectangularity, W_ImageObjPerimeter, M_Moments, M_RawMoments
	killwaves W_BoundaryX, W_BoundaryY, W_BoundaryIndex, W_xmin, W_xmax, W_ymin, W_ymax, xmin, xmax, ymin, ymax
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
function GaussianFit(image, x_loc, y_loc, size, wiggle)
	
	wave image, x_loc, y_loc
	variable size, wiggle
	variable half_size = size / 2
	variable num_peaks = DimSize(x_loc,0)
	variable success, V_FitError
	variable chisq_keep
	//newimage image

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

	Make/O/N=(num_peaks) z0, A, x0, xW, y0, yW, cor, chisq
	Make/O/N=(num_peaks) sigma_z0, sigma_A, sigma_x0, sigma_xW, sigma_y0, sigma_yW, sigma_cor

	variable i
	for(i=0; i<num_peaks; i+=1)
		
		if(wiggle == 0)
			V_FitError =0
			CurveFit/M=2/W=2/Q gauss2D, image[x_start[i],x_finish[i]][y_start[i],y_finish[i]] /W=noise /I=1 /D
			if(V_FitError)
				success = 0
			endif
		else	
			success = OneGaussFitWiggle(image, x_loc[i], y_loc[i], size, wiggle)
		endif
		
		wave W_coef = $"W_coef"
		wave W_sigma = $"W_sigma"
		
		if(W_coef[2] <= x_start[i] || W_coef[2] >= x_finish[i])
			success = 0
		endif
		if(W_coef[4] <= y_start[i] || W_coef[4] >= y_finish[i])
			success = 0
		endif
		
		if(success)
			chisq[i] = chisq_keep //how to report chisq in wave???
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
	//appendtograph/t y0 vs x0
	//modifygraph mode=2
	//modifygraph lsize=3

	killwaves W_sigma, W_coef, M_covar
end

//This function is used to fit each image in a stack with the same initial guess in x_loc, y_loc
//12/06/2016 cz
function StackGaussianFit(im_stack, x_loc, y_loc, size)

	wave im_stack, x_loc, y_loc
	variable size
	variable num_frames = DimSize(im_stack,2)
	variable num_peaks = Dimsize(x_loc,0)
	variable total_peaks = num_frames*num_peaks
	variable half_size = size/2
	
	duplicate/O x_loc x_start
	duplicate/O x_loc x_finish
	duplicate/O y_loc y_start
	duplicate/O y_loc y_finish
	
	x_start = x_loc - half_size
	x_finish = x_loc + half_size
	y_start = y_loc - half_size
	y_finish = y_loc + half_size

	Make/O/N=(num_peaks,num_frames) z0_stack, A_stack, x0_stack, xW_stack, y0_stack, yW_stack, cor_stack
	Make/O/N=(num_peaks,num_frames) sigma_z0_stack, sigma_A_stack, sigma_x0_stack, sigma_xW_stack, sigma_y0_stack, sigma_yW_stack, sigma_cor_stack
	
	variable i
	for(i=0; i<num_frames; i+=1)
		Imagetransform/p=(i) getplane im_stack
		wave im = $"M_ImagePlane"
		GaussianFit(im, x_loc, y_loc, size, 1) //1 is used as default wiggle size as it doesn't matter much
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
		z0_stack[][i] = z0[p]
		A_stack[][i] = A[p]
		x0_stack[][i] = x0[p]
		xW_stack[][i] = xW[p]
		y0_stack[][i] = y0[p]
		yW_stack[][i] = yW[p]
		cor_stack[][i] = cor[p]
		sigma_z0_stack[][i] = sigma_z0[p]
		sigma_A_stack[][i] = sigma_A[p]
		sigma_x0_stack[][i] = sigma_x0[p]
		sigma_xW_stack[][i] = sigma_xW[p]
		sigma_y0_stack[][i] = sigma_y0[p]
		sigma_yW_stack[][i] = sigma_yW[p]
	endfor
	killwaves z0, A, x0, xW, y0, yW, cor, sigma_z0, sigma_A, sigma_x0, sigma_xW, sigma_y0, sigma_yW
end

// Fits the region size x size about the center position (x_loc, y_loc) to a 2D Guassian,
// but repeats the fit shifting the region around from -wiggle to +wiggle in x and y.  
// The fit with the minimum chisq is considered the best, and the coefficients are maintained
// in waves W_Coef and W_Sigma. All other fitting coefficients are discarded.
function OneGaussFitWiggle(image, x_loc, y_loc, size, wiggle)
	wave image
	variable x_loc, y_loc, size, wiggle
	
	Duplicate/O image noise
	noise = sqrt(image)

	variable x_start, x_finish, y_start, y_finish, half_size
	half_size = size/2
	
	
	variable V_chisq_min = Inf
	
	variable i, j, V_FitError, success = 0
	for(i=-1.0*wiggle; i<=wiggle; i+=1)
		for(j=-1.0*wiggle; j<=wiggle; j+=1)

			x_start = (x_loc - half_size) + i
			x_finish = (x_loc + half_size) + i
			y_start = (y_loc - half_size) + j
			y_finish = (y_loc + half_size) + j

			V_FitError = 0
			CurveFit/N/M=2/W=2/Q gauss2D,image[x_start,x_finish][y_start,y_finish] /W=noise /I=1 /D
		
			wave W_Coef = $"W_Coef"
			//printf "(x, y) = (%g, %g);  chisq = %g\r", W_coef[2], W_coef[4], V_chisq
		
			if(!V_FitError)	  // if an error occured during fitting, ignore the results
				if(V_chisq < V_chisq_min)	// keep solutions with lower chisq
					if( (W_coef[2] > x_start) && (W_coef[2] < x_finish) )	// keep solutions with a center inside the fit box
						if( (W_coef[4] > y_start) && (W_coef[4] < y_finish) )
							Duplicate/O W_Coef keep_coef
							Duplicate/O W_Sigma keep_sigma
							V_chisq_min = V_chisq
							success = 1
						endif
					endif
				endif
			endif
		endfor
	endfor
	
	if(success)
		printf "(x, y) = (%g, %g);  chisq = %g\r", W_coef[2], W_coef[4], V_chisq_min
		Duplicate/O keep_coef W_coef
		Duplicate/O keep_sigma W_sigma
		variable chisq_keep
		chisq_keep = V_chisq_min
		Killwaves keep_coef, keep_sigma, noise
	endif
	
	return success
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
		FuncFitMD/N/NTHR=0/Q Gaus2Dx2 W_coef  image[x_start(i),x_finish(i)][y_start(i),y_finish(i)] /W=noise /I=1 /D
		
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
				if(x_angle<0||y_angle<0)
					xm = tmp_x0(j) - tmp_x0(i) 
					ym = tmp_y0(j) - tmp_y0(i)
				else
					xm = abs(tmp_x0(j) - tmp_x0(i)) 
					ym = abs(tmp_y0(j) - tmp_y0(i))
				endif
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
				if(x_angle<0||y_angle<0)
					xm = tmp_x0(j) - tmp_x0(i) 
					ym = tmp_y0(j) - tmp_y0(i)
				else
					xm = abs(tmp_x0(j) - tmp_x0(i)) 
					ym = abs(tmp_y0(j) - tmp_y0(i))
				endif
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
	
	ThresholdPeakPositions(image, 12000)
	
	wave x_loc =$"x_loc"
	wave y_loc =$"y_loc"
	if(!WaveExists(x_loc) || !WaveExists(y_loc))
		print "PeakPositions failed to create x_loc and/or y_loc.  Exiting.\r"
		return 0
	endif
	
	GaussianFit(image, x_loc, y_loc, size, 1)	// hard coded for zero fitting box wiggle
	
	Separation(x_space, y_space, space_delta, x_angle, y_angle, angle_delta)	
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
	Label left "Average X Seperation (A)";DelayUpdate
	Label bottom "Frame"

	Display  xstd
	ModifyGraph rgb(xstd)=(0,0,0);DelayUpdate
	Label left "X Stdev (A)";DelayUpdate
	Label bottom "Frame"

	Display  yav
	ModifyGraph rgb(yav)=(16384,48896,65280);DelayUpdate
	ErrorBars/T=0 yav Y,wave=(ystd,ystd)
	AppendToGraph yav
	ModifyGraph rgb(yav#1)=(0,0,0)
	Label left "Average Y Seperation (A)";DelayUpdate
	Label bottom "Frame"

	Display  ystd
	ModifyGraph rgb(ystd)=(0,0,0);DelayUpdate
	Label left "Y Stdev (A)";DelayUpdate
	Label bottom "Frame"
	
	variable i
	for(i=0; i<DimSize(image_stack, 2); i+=1)
	//for(i=0; i<4; i+=1)
		Imagetransform/p=(i) getplane image_stack
		wave im = $"M_ImagePlane"
		//Precision(im, size, x_space, y_space, space_delta, x_angle, y_angle, angle_delta)
		//Used fixed x_loc and y_loc in this scheme to save some time
		GaussianFit(im, x_loc, y_loc, size, 1)		
		Separation(x_space, y_space, space_delta, x_angle, y_angle, angle_delta)	
	
		//save averages and stdevs
		wave prec = $"prec"
		//xav[i] = prec[0]*0.109
		//xstd[i] = prec[1]*0.109
		//yav[i] = prec[2]*0.109
		//ystd[i] = prec[3]*0.109
		//use output in pixels for both precision and separations
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



function ABFDumSimStackPosition(image_stack, size)

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

function cropimage(im, xi, yi, xf, yf)

wave im
variable xi, yi, xf, yf

variable scalex = DimDelta(im, 0)
variable scaley = DimDelta(im, 1)

variable xi_pixel = (xi/scalex)
variable xf_pixel = (xf/scalex) + 1
variable yi_pixel = (yi/scaley)
variable yf_pixel = (yf/scaley) + 1

DeletePoints xf_pixel,10000, im
DeletePoints 0,xi_pixel, im
DeletePoints/M=1 yf_pixel,10000, im
DeletePoints/M=1 0,yi_pixel, im

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

// Integrated Intensity calculation of an all atomic columns in x0 & y0 that have been fit with Gaussian functions and outputed x0, y0, xW, yW, cor, z0 and A waves.
// If "pixel_ size" wave is present in the root directory, the code will automatically find it and do the integrated intensities with these pixel sizes included.
// If "pixel_ size" wave is not present in the root directory, the integrated intensities will be calculated with a pixel size of 1.
// "pixel_size" wave is a 2 row / 1 column wave where the first cell is the x _pixel size and the second cell is the y_pixel size.
// size is the area that the integrated intensity is calculted in.
// If visibility_map = 0, no visibility map will be displayed. If it is = 1 then a visibility map will be shown.
// If intensity_map = 0, no intensity map will be shown. If it is = 1 then a visibility map with be shown.
function IntegratedGaussianIntensity(image, size, intensity_map, visibility_map)
	
	wave image
	variable size, intensity_map, visibility_map
	wave cor = $"cor"
	wave z0 = $"z0"
	wave A = $"A"
	wave xW = $"xW"
	wave yW = $"yW"
	wave x0 = $"x0"
	wave y0 = $"y0"
	if(!WaveExists(cor) || !WaveExists(z0) || !WaveExists(A) || !WaveExists(xW) || !WaveExists(yW))
		print "Failed to find GaussianFitting parameters while calculate intensity.\r"
		return 0
	endif
	
	variable pixel = 0
	wave pixel_size = $"pixel_size"
	if(!WaveExists(pixel_size))
		pixel = 0
		print "Failed to find the pixel_size wave. Calculting integrated intensities with pixel size 1."
	else
		pixel = 1
		print "Found the wave pixel_size and used it for the integrated intensity calculation."
		variable pixel_x = pixel_size[0]
		variable pixel_y = pixel_size[1]
	endif
	
	variable peak_num = DimSize(cor,0)
	Make/N = (peak_num)/D/O Intensity
	Make/N = (peak_num)/D/O Visibility
	
	// calculate integrated intensity for each atomic column
	variable i
	for (i = 0; i < peak_num; i +=1)
		if (pixel == 0)
			Intensity[i] = (z0[i]*(size^2))+(2*A[i]*sqrt(1-(cor[i]^2))*pi*abs(xW[i])*abs(yW[i]))
		endif
		if (pixel ==1)
			Intensity[i] = (pixel_x*pixel_y*z0[i]*(size^2))+(2*A[i]*sqrt(1-(cor[i]^2))*pi*abs(xW[i]*pixel_x)*abs(yW[i]*pixel_y))
		endif
	endfor
	
	// calculate visibility of each atomic column compared to its nearest neighbors.
	variable intensity_tmp, intensity_tmp2, x0_tmp, y0_tmp, x0_tmp2, y0_tmp2, d
	variable n, o, p
	for (n = 0; n < peak_num; n +=1)
		intensity_tmp = Intensity[n]
		x0_tmp = x0[n]
		y0_tmp = y0[n]
		make/o/n=0 int_tmp2log
		for (o = 0; o < peak_num; o +=1)
			if( n != o)
				intensity_tmp2 = Intensity[o]
				x0_tmp2 = x0[o]
				y0_tmp2 = y0[o]
				d = sqrt(((x0_tmp2 - x0_tmp)^2) + ((y0_tmp2 - y0_tmp)^2))
				if ( d <= 4 * (mean(xW)+mean(yW)))
					InsertPoints 0, 1, int_tmp2log
					int_tmp2log[x2pnt(int_tmp2log,0)] = intensity_tmp2
				endif
			endif
		endfor
		Visibility[n] = (1- (intensity_tmp / mean(int_tmp2log))) * 100
		killwaves int_tmp2log
	endfor
	
	variable minimum, maximum
	ColorTab2Wave YellowHot
	wave colors = $"M_colors"
	variable color_size = dimsize(colors,0)
	
	// show intensity map if = 1.
	if ( intensity_map == 1)
		WaveStats/q Intensity
		minimum = V_min
		maximum = V_max
		NewImage image
		variable  x_position, y_position
		variable j,k
		for (j = 0; j < peak_num; j +=1)
			intensity_tmp = Intensity[j]
			x_position = x0[j]/(dimsize(image,0))
			y_position = y0[j]/(dimsize(image,1))
			variable color_1, color_2, color_3
			for (k = 0; k < color_size; k +=1)
				if (intensity_tmp >= minimum + (k*(maximum - minimum)/color_size))
					if (intensity_tmp <= minimum + ((k+1) * (maximum - minimum)/color_size))
						color_1= colors[k][0]
						color_2= colors[k][1]
						color_3= colors[k][2]
						SetDrawEnv fillbgc= (color_1,color_2,color_3),fillfgc= (color_1,color_2,color_3),linethick= 0.00;DelayUpdate
					endif
				endif
			endfor
			DrawOval x_position-0.01, y_position-0.01, x_position+0.01, y_position+0.01
			ColorScale/C/N=text0/F=0/A=MC/X=0.00/Y=0.00  ctab={minimum,maximum,YellowHot,0},lblMargin=0;DelayUpdate
			ColorScale/C/N=text0 "Intensity"
		endfor
	endif
	
	// show visibility map if = 1.
	if ( visibility_map == 1)
		WaveStats/q Visibility
		minimum = V_min
		maximum = V_max
		NewImage image
		variable visibility_tmp, x_position2, y_position2
		variable l, m
		for (l = 0; l < peak_num; l +=1)
			visibility_tmp = visibility[l]
			x_position2 = x0[l]/(dimsize(image,0))
			y_position2 = y0[l]/(dimsize(image,1))
			variable color_4, color_5, color_6
			for (m = 0; m < color_size; m +=1)
				if (visibility_tmp >= minimum + (m*(maximum - minimum)/color_size))
					if (visibility_tmp <= minimum + ((m+1) * (maximum - minimum)/color_size))
						color_4= colors[m][0]
						color_5= colors[m][1]
						color_6= colors[m][2]
						SetDrawEnv fillbgc= (color_4,color_5,color_6),fillfgc= (color_4,color_5,color_6),linethick= 0.00;DelayUpdate
					endif
				endif
			endfor
			DrawOval x_position2-0.01, y_position2-0.01, x_position2+0.01, y_position2+0.01
			ColorScale/C/N=text0/F=0/A=MC/X=0.00/Y=0.00  ctab={minimum,maximum,YellowHot,0},lblMargin=0;DelayUpdate
			ColorScale/C/N=text0 "Visibility"
		endfor	
	endif
	killwaves M_colors
end

