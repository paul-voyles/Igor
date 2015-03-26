#pragma rtGlobals=1		// Use modern global access method.
//
// Functions for processing STEM image stacks. Peak finding, peak fitting, and precision 
// measurements.
//
// List of Functions:
//
// PeakPositions: Input image wave name. Uses the IGOR ImageAnalyzeParticle function to find the 
// positions of all the atoms, not counting the atoms on the perimiter of the image.
//  
// GaussianFit: Input wave name, x_loc, y_loc, and width of fit area in pixels. This functions 
// conducts a 2D Gaussian fit for each peak found in the PeakPositions function which are inputed in x_loc and y_loc.
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
// Started 12/2/14, JF


// This function converts images in counts to electrons
// inputs: original image in HAADF counts, numSamples, and cm = average count in HAADF probe image,
// c0 = average count HAADF dark image, p = probe current, t = pixel dwell time.
function CountstoElectrons(average, numSamples, cm, c0, p, t)
	
	wave average, numSamples 
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
	
	//Auto Threshold
	ImageThreshold/I/M=(1)/Q image
	ImageAnalyzeParticles /E/W/Q/F/M=3/A=15/EBPC stats, root:M_ImageThresh
	
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
end

// This function revised algorithm of x_loc and y_loc in PeakPosition above.
// 12/2/14 by JFeng
function ActPeakPositions(image)
	
	wave image
	
	NewImage image
	
	ImageThreshold/I/M=(1)/Q image
	ImageAnalyzeParticles /E/W/Q/F/M=3/A=10/EBPC stats, root:M_ImageThresh
	
	duplicate/O W_xmin x_loc	
	duplicate/O W_ymin y_loc
	wave W_xmin = $"W_xmin"
	wave W_ymin = $"W_ymin"
	wave W_xmax = $"W_xmax"
	wave W_ymax = $"W_ymax"
	
	variable n = 0 
	for (n = 0; n < V_NumParticles; n+=1)
		variable xMax = W_xmin[n]
		variable yMax = W_ymin[n]
		
		variable tmpMax = image[W_xmin[n]][W_ymin[n]]
		
		variable i = W_xmin[n]
		do
			variable j = W_ymin[n]
			do
				if (tmpMax < image[i][j])
					tmpMax = image[i][j]
					xMax = i
					yMax = j
				endif
				j += 1
			while (j <= W_ymax[n])
			i += 1
		while (i <= W_xmax[n])	
		
		x_loc[n] = xMax
		y_loc[n] = yMax
		
	endfor
	
	appendtograph/t y_loc vs x_loc
	ModifyGraph mode=2
	
	killwaves W_ImageObjArea, W_SpotX, W_SpotY, W_circularity, W_rectangularity, W_ImageObjPerimeter, M_Moments, M_RawMoments
	killwaves W_BoundaryX, W_BoundaryY, W_BoundaryIndex, W_xmin, W_xmax, W_ymin, W_ymax
end


// Does the peak find, 2dGaus fit, and x and y precision calculations on one image of single atom columns.
// Outputs sep_out with x and y avergae and stdev.
//x_sep_1 = s_0,  x_sep_2 = s_1,  y_sep = s_2
function Precision(image, size, row_num, x_sep_1, x_sep_2, y_sep)

	wave image
	variable size, row_num, x_sep_1, x_sep_2, y_sep	
	
	//ActPeakPositions(image)
	
	wave x_loc =$"x_loc"
	wave y_loc =$"y_loc"
	if(!WaveExists(x_loc) || !WaveExists(y_loc))
		print "PeakPositions failed to create x_loc and/or y_loc.  Exiting.\r"
		return 0
	endif
	
	GaussianFit(image, x_loc, y_loc, size, 1)	// hard coded for 1 fitting box wiggle
	PixelToAng(row_num, x_sep_1, x_sep_2, y_sep)	//calculate pixels sizes of experiment based on its thickness 
														//and the cooresponding simulation results
	Separation(row_num, 1)
	GauInteg(row_num, size)	
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
			CurveFit/N/NTHR=0/W=2/Q gauss2D, im/R=residual//W=noise
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

// transfer from pixel to Angstrom
// row_num is the number of rows, x_sep_1, x_sep_2 and y_sep are the 3 kinds of separations with unit of Angstorm
function PixelToAng(row_num, x_sep_1, x_sep_2, y_sep)
	
	variable row_num, x_sep_1, x_sep_2, y_sep
	wave x0 = $"x0"
	wave y0 = $"y0"
	variable peaks_num = DimSize(x0, 0)
	variable column_num = peaks_num / row_num
	
	variable x_start, y_start, x_end, y_end, PixelSize_X, PixelSize_Y, i
	variable x_start_ave, x_end_ave, x_lentgh
	x_start = 0
	x_end = 0
	for (i = 0; i < peaks_num; i += column_num)
		x_start += x0(i)
	endfor
	for (i = column_num - 1; i < peaks_num; i += column_num)
		x_end += x0(i)
	endfor
	x_start_ave = x_start / row_num
	x_end_ave = x_end / row_num
	x_lentgh = (x_sep_1 + x_sep_2) * (column_num - 1) / 2
	PixelSize_X = x_lentgh / (x_end_ave - x_start_ave)
	
	y_start = 0
	y_end = 0
	for (i = 0; i < column_num; i += 1)
		y_start += y0(i)
	endfor
	for (i = peaks_num - column_num; i < peaks_num; i +=1)
		y_end += y0(i)
	endfor
	PixelSize_Y = y_sep * (row_num - 1) / (y_end / column_num - y_start / column_num)
	Make /N = 2 /D/O PixelSize
	PixelSize[0] = PixelSize_X
	PixelSize[1] = PixelSize_Y
	print "PixelSize_X = ", PixelSize[0]
	print "PixelSize_Y = ", PixelSize[1]
end


// This one is being used in the Precision and Stack functions. It should work for all images if the inputs are correct.
// Inputs: row_num is the num of rows in the image. Require each row has the same number of columns 
function Separation(row_num, pixel_size)

	variable row_num, pixel_size
	
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
	
	variable column_num = num_peaks / row_num
	print "column_num =", column_num
	print "row_num =", row_num
	Make/O/N=0 x_distance_1, x_distance_2
	Make/O/D/N = ((row_num - 1) * column_num) y_distance
	
	//To calculate x_distance_1 and x_distance_2
	variable i, j, xm, ym, d
	variable n = 0
	for (i = 0; i < row_num; i += 1)
		for (j = 0; j < column_num - 1; j += 1)
			if (mod(i, 2) == 0)	// for even index row
				if (mod(j, 2) == 0)
					xm = abs(tmp_x0(n) - tmp_x0(n+1))		
					ym = abs(tmp_y0(n) - tmp_y0(n+1))		
					d = sqrt((xm^2) + (ym^2))
					InsertPoints 0, 1, x_distance_1
					x_distance_1[0] = d
					n += 1
				else
					xm = abs(tmp_x0(n) - tmp_x0(n+1))	
					ym = abs(tmp_y0(n) - tmp_y0(n+1))
					d = sqrt((xm^2) + (ym^2))
					InsertPoints 0, 1, x_distance_2
					x_distance_2[0] = d
					n += 1
				endif
			else		// for odd index row
				if (mod(j, 2) == 0)
					xm = abs(tmp_x0(n) - tmp_x0(n+1))
					ym = abs(tmp_y0(n) - tmp_y0(n+1))
					d = sqrt((xm^2) + (ym^2))
					InsertPoints 0, 1, x_distance_2
					x_distance_2[0] = d
					n += 1
				else
					xm = abs(tmp_x0(n) - tmp_x0(n+1))
					ym = abs(tmp_y0(n) - tmp_y0(n+1))
					d = sqrt((xm^2) + (ym^2))
					InsertPoints 0, 1, x_distance_1
					x_distance_1[0] = d
					n += 1
				endif 
			endif
		endfor
		n += 1
	endfor
	
	//convert x_distance_1 and x_distance_2
	variable x1_num = DimSize(x_distance_1, 0)
	Duplicate x_distance_1, x_distance_1_copy
	Duplicate x_distance_2, x_distance_2_copy
	for (i = 0; i < x1_num; i += 1)
		x_distance_1[i] = x_distance_1_copy[x1_num - 1 - i]
		x_distance_2[i] = x_distance_2_copy[x1_num - 1 - i]
	endfor
	
	killwaves x_distance_1_copy, x_distance_2_copy
	
	x_distance_1 = pixel_size * x_distance_1
	wavestats/Q/W x_distance_1
	wave M_WaveStats = $"M_WaveStats"
	variable x_dis_1_ave = M_WaveStats(3)
	variable x_dis_1_std = M_WaveStats(4)
	
	x_distance_2 = pixel_size * x_distance_2
	wavestats/Q/W x_distance_2
	variable x_dis_2_ave = M_WaveStats(3)
	variable x_dis_2_std = M_WaveStats(4)
	
	//To calculate y_distance
	n = 0
	for (i = 0; i < row_num - 1; i += 1)
		for (j = 0; j < column_num; j += 1)
			xm = abs(tmp_x0(column_num * i + j) - tmp_x0(column_num * (i + 1) + j))
			ym = abs(tmp_y0(column_num * i + j) - tmp_y0(column_num * (i + 1) + j))
			y_distance[n] = sqrt((xm^2) + (ym^2))	
			n += 1
		endfor
	endfor
	
	y_distance = pixel_size * y_distance
	wavestats/Q/W y_distance
	variable y_dis_ave = M_WaveStats(3)
	variable y_dis_std = M_WaveStats(4)
	
	//Write Prec
	Wave PixelSize = $"PixelSize"
	variable PixelSize_X, PixelSize_Y
	PixelSize_X = PixelSize[0]
	PixelSize_Y = PixelSize[1]
	Make/O/N = (6, 2) Prec
	Prec[0][0] = x_dis_1_ave
	Prec[1][0] = x_dis_1_std
	Prec[2][0] = x_dis_2_ave
	Prec[3][0] = x_dis_2_std
	Prec[4][0] = y_dis_ave
	Prec[5][0] = y_dis_std
	Prec[0][1] = x_dis_1_ave * PixelSize_X
	Prec[1][1] = x_dis_1_std * PixelSize_X
	Prec[2][1] = x_dis_2_ave * PixelSize_X
	Prec[3][1] = x_dis_2_std * PixelSize_X
	Prec[4][1] = y_dis_ave * PixelSize_Y
	Prec[5][1] = y_dis_std * PixelSize_Y
	
	// Draw Histogram for x_distance_1, x_distance_2 and y_distance
	variable bin_num_dx1 = round(DimSize(x_distance_1,0) / 4)
	variable bin_num_dx2 = round(DimSize(x_distance_2,0) / 4)
	variable bin_num_dy = round(DimSize(y_distance,0) / 4)
	print bin_num_dx1
	Make/N=(3, 4)/D GauFit_coef
	// X_distance_1
	Make/N=(bin_num_dx1)/O x_distance_1_Hist
	Histogram/B=1 x_distance_1,x_distance_1_Hist
	//CurveFit/M=2/W=0 gauss, x_distance_1_Hist/D
	//wave W_coef = $"W_coef"
	//GauFit_coef[0][0] = W_coef[0]
	//GauFit_coef[0][1] = W_coef[1]
	//GauFit_coef[0][2] = W_coef[2]
	//GauFit_coef[0][3] = W_coef[3]
	
	// X_distance_2
	Make/N=(bin_num_dx2)/O x_distance_2_Hist
	Histogram/B=1 x_distance_2,x_distance_2_Hist
	//CurveFit/M=2/W=0 gauss, x_distance_2_Hist/D
	//wave W_coef = $"W_coef"
	//GauFit_coef[1][0] = W_coef[0]
	//GauFit_coef[1][1] = W_coef[1]
	//GauFit_coef[1][2] = W_coef[2]
	//GauFit_coef[1][3] = W_coef[3]
	
	// Y_distance
	Make/N=(bin_num_dy)/O y_distance_Hist
	Histogram/B=1 y_distance,y_distance_Hist
	//CurveFit/M=2/W=0 gauss, y_distance_Hist/D
	//wave W_coef = $"W_coef"
	//GauFit_coef[2][0] = W_coef[0]
	//GauFit_coef[2][1] = W_coef[1]
	//GauFit_coef[2][2] = W_coef[2]
	//GauFit_coef[2][3] = W_coef[3]
	
	//killwaves W_coef, M_Covar, W_sigma
	
	// Another way to calculate x_distance, don't distinguish x_distance_1 and x_distance_2.
	// Use this to check whether the method distinguishing x_distance_1 and x_distance_2 is right.
	// The results shows all Seps calculated by two methods are the same. So should be right.
	//make /O/D/N=(row_num*(column_num-1)) x_distance
	//n = 0
	//for (i = 0; i < row_num; i+=1)
	//	for (j = 0; j < column_num - 1; j+=1)
	//		x_distance[n] = pixel_size * sqrt((tmp_x0[column_num*i+j] - tmp_x0[column_num*i+j+1])^2 + (tmp_y0[column_num*i+j] - tmp_y0[column_num*i+j+1])^2)
	//		n += 1
	//	endfor
	//endfor
end


// Intensity calculation
function GauInteg(row_num, size)
	
	variable row_num, size
	wave cor = $"cor"
	wave z0 = $"z0"
	wave A = $"A"
	wave xW = $"xW"
	wave yW = $"yW"
	Wave PixelSize = $"PixelSize"
	if(!WaveExists(cor) || !WaveExists(z0) || !WaveExists(A) || !WaveExists(xW) || !WaveExists(yW) || !WaveExists(PixelSize))
		print "Failed to find GaussianFitting parameters while calculate intensity.\r"
		return 0
	endif
	variable peaks_num = DimSize(cor,0)
	variable column_num = peaks_num / row_num
	Make/N = (row_num, column_num)/D/O Intensity
	variable i, j
	variable k = 0
	variable PixelSize_X = PixelSize[0]
	variable PixelSize_Y = PixelSize[1]
	for (i = 0; i < row_num; i +=1)
		for (j = 0; j < column_num; j +=1)
			Intensity[i][j] = z0[k]*size^2*PixelSize_X*PixelSize_Y+2*A[k]*sqrt(1-cor[k]^2)*pi*abs(xW[k]*PixelSize_X)*abs(yW[k]*PixelSize_Y)
			k += 1
		endfor
	endfor
	variable MeanIntensity = Mean(Intensity)
	Duplicate Intensity IntensityVisibility
	for (i = 0; i < row_num; i +=1)
		for (j = 0; j < column_num; j +=1)
			IntensityVisibility[i][j] = (1 - Intensity[i][j] / MeanIntensity) * 100		//%
		endfor
	endfor
end

//More Precise measure after Precision
//Crop the target column which gives by (x, y)
// x is row_num, y is column_num 0-index
function MPrecision(row_num, x, y)

	variable row_num, x, y
	wave x0 = $"x0"
	wave Intensity = $"Intensity"
	wave x_distance_1 = $"x_distance_1"
	wave x_distance_2 = $"x_distance_2"
	wave y_distance = $"y_distance"
	if(!WaveExists(x0) || !WaveExists(Intensity) || !WaveExists(x_distance_1) || !WaveExists(x_distance_2) || !WaveExists(y_distance))
		print "Failed to find wave in function MPrecision(x, y).\r"
		return 0
	endif
	
	variable peaks_num = DimSize(x0, 0)
	variable column_num = peaks_num / row_num
	variable MeanIntensity = (Sum(Intensity) - Intensity[x][y]) / (peaks_num - 1)
	Duplicate Intensity IntensityVisibility_crop
	variable i, j
	for (i = 0; i < row_num; i +=1)
		for (j = 0; j < column_num; j +=1)
			IntensityVisibility_crop[i][j] = (1 - Intensity[i][j] / MeanIntensity) * 100		//%
		endfor
	endfor
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
