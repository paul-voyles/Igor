#pragma rtGlobals=1		// Use modern global access method.



// This function converts images to an absolute intensity scale that can be compared to simulations.
// inputs: im is the image, Cm and Co are the average pixel intensity of the HAADF probe image and dark image respectively.
function AbsoluteIntensity(im, Cm, Co)
	
	wave im
	variable Cm, Co

	Duplicate/O im $(NameofWave(im)+"_AbsInt")
	wave/C AbsInt = $(NameofWave(im)+"_AbsInt")
	AbsInt = (im - Co) / (Cm - Co)
	
end


function CSHalfAtomInt(IntVsSlice)
	
	wave IntVsSlice
	variable numSlices = DimSize(IntVsSlice,0)
	variable numCS = (2*(numSlices-1)) +1

	Interpolate2/T=2/N=(numCS)/E=2/Y=$(NameofWave(IntVsSlice)+"_CS") $(NameofWave(IntVsSlice))
	wave/C CS = $(NameofWave(IntVsSlice)+"_CS")
	
	Make/o/n=(numSlices-2) CS_Half_Int
	
	variable i
	for(i=0; i<(numSlices-2); i+=1)
		
		CS_Half_Int[i] = CS[(2*i)+3]
		
	endfor	
	
end


function AtomCount(AtomInt, CS_Half_Int)
	
	wave AtomInt, CS_Half_Int
	
	variable numAtomCol = DimSize(AtomInt, 0)
	variable numSlices = DimSize(CS_Half_Int, 0)
	Make/o/n=(numAtomCol) Atom_Count
	
	variable int
	
	variable i, j = 0
	for(i=0; i<numAtomCol; i+=1)
	
		int= AtomInt[i]
		
		if( int <= CS_Half_Int[0])
			Atom_Count[i] = 1
		endif
		
		if( int > CS_Half_Int[numSlices-1])
			print "Intensity out of range:", i
		endif
		
		j=0
		for(j=0; j<numSlices; j+=1)
			if( int > CS_Half_Int[j])
				if( int <= CS_Half_Int[j+1])
					Atom_Count[i] = j+2
				endif
			endif
		endfor
	endfor
end
	


// This function finds the peak positions using the analyze particle IGOR function and reports the locations in x_loc and y_loc.
// For Si dumbells, it reports the center of the dumbells. To change the mimimun area of the particles change the A flag in the ImageAnalyzeParticle function.
// Inputs: image
function PeakPositions(image)
	
	wave image
	
	NewImage image
	
	ImageThreshold/I/M=(1)/Q image
	ImageAnalyzeParticles /E/W/Q/F/M=3/A=0/EBPC stats, root:M_ImageThresh
	
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



// This function fits the peak positions that are in x_loc and y_loc that were found in the previous routine.
// This only works for 2D Gaussian fits. Si dumbbells need a different fitting routine.
// Inputs: image, x_loc, y_loc, and gaussian fit size x size.
function GaussianFit(image, x_loc, y_loc, size, wiggle)
	
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
			CurveFit/M=2/W=2/Q gauss2D, image[x_start[i],x_finish[i]][y_start[i],y_finish[i]] /W=noise /I=1 /D//R=residual
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
function OneGaussFitWiggle(image, x_loc, y_loc, size, wiggle)
	wave image
	variable x_loc, y_loc, size, wiggle
	
	Duplicate/O image noise
	//Make/O/N=(DimSize(image,0),DimSize(image,1))  residual, single_residual
	
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
			CurveFit/N/M=2/W=2/Q gauss2D,image[x_start,x_finish][y_start,y_finish] /W=noise /I=1 /D//R=single_residual
		
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
							//residual = single_residual
						endif
					endif
				endif
			endif
		endfor
	endfor
	
	Duplicate/O keep_coef W_coef
	Duplicate/O keep_sigma W_sigma
	
	Killwaves keep_coef, keep_sigma, noise//, single_residual
	
end

// This function rounds the atom positions found from fitting to the nearest pixel position.
// input x0 and y0, atom fit positions.
// outputs x0_rnd and y0_rnd which are the rounded positions. 
function PositionRound(x0,y0)
	
	wave x0, y0
	
	duplicate/o x0 x0_rnd
	duplicate/o y0 y0_rnd
	
	variable num, x0_i, x0_f, y0_i, y0_f
	num = DimSize(x0,0)
	
	variable i
	for(i=0; i<num; i+=1)
		x0_i = x0[i]
		y0_i = y0[i]
		x0_f = round(x0_i)
		y0_f = round(y0_i)
		x0_rnd[i] = x0_f
		y0_rnd[i] = y0_f
	endfor
end


function PositionRoundDown(x0,y0)
	
	wave x0, y0
	
	duplicate/o x0 x0_floor
	duplicate/o y0 y0_floor
	
	variable num, x0_i, x0_f, y0_i, y0_f
	num = DimSize(x0,0)
	
	variable i
	for(i=0; i<num; i+=1)
		x0_i = x0[i]
		y0_i = y0[i]
		x0_f = floor(x0_i)
		y0_f = floor(y0_i)
		x0_floor[i] = x0_f
		y0_floor[i] = y0_f
	endfor
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


// This function crops the image (im) with a window size (size) that must be an odd number around each atom position inputed from x0 and y0.
// outputs a stack of cropped images (crop_stack) with the order being the order of x0, y0. 
function CropAndStack(im,size,x0,y0)
	
	wave im, x0, y0
	variable size
	
	variable num 
	num = DimSize(x0,0)
	
	Make/o/n=(size, size, num) crop_stack
	Make/o/n=(num) xi_all, yi_all, xf_all, yf_all
	variable xi, yi, xf, yf
	
	variable i
	for(i=0; i<num; i+=1)
		
		Duplicate/o im im_crop
		xi = x0[i] - ((size-1)/2)
		yi = y0[i] - ((size-1)/2)
		xf = x0[i] + ((size-1)/2)
		yf = y0[i] + ((size-1)/2)
		
		xi_all[i] = xi
		yi_all[i] = yi
		xf_all[i] = xf
		yf_all[i] = yf
		
		cropimage(im_crop, xi, yi, xf, yf)
		
		crop_stack[][][i] = im_crop[p][q]
		
		killwaves im_crop
	endfor
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



function AddAtomNum(im, x0, y0, Atom_Count)
	
	wave im, x0, y0, Atom_Count
	NewImage im
	SetDrawEnv fsize=6, xcoord=prel, ycoord=prel,  textrgb= (65280,0,0), save
	
	variable size = DimSize(Atom_Count, 0)
	
	wavestats/Q/W Atom_Count
	wave M_WaveStats = $"M_WaveStats"
	variable numatoms = M_WaveStats(12)
	
	variable x, y
	
	variable i, j 
	for(i=0; i<size; i+=1)
		j=0
		for(j=0; j<numatoms; j+=1)
			if( Atom_Count[i] == j+1 )
				x = (x0[i] / 256) - 0.0005
				y = (y0[i] / 256) + 0.008
				string n 
				sprintf n, "%d", j+1
				
				//TextBox/C/N=text0/A=MC n

                		SetDrawLayer UserFront
				DrawText x,y, n
			endif
		endfor
	endfor
	
	killwaves M_WaveStats
	//Legend/C/N=text0/J/A=MC "\\K(39168,0,0)\\W5199 Atoms\r\\K(65280,21760,0)\\W5198 Atoms\r\\K(65280,43520,0)\\W5197 Atoms\r\\K(65280,65280,0)\\W5196 Atoms";DelayUpdate
	
end




function AddAtomInt(im, x0, y0, Atom_Count)
	
	wave im, x0, y0, Atom_Count
	NewImage im
	
	variable size = DimSize(Atom_Count, 0)
	
	variable l, t, r, b
	
	wavestats/Q/W Atom_Count
	wave M_WaveStats = $"M_WaveStats"
	variable numatoms = M_WaveStats(12)
	
	ColorTab2Wave Geo32
	wave M_colors = $"M_colors"
	Variable N= DimSize(M_colors,0)
	
	Make/O/N=(numatoms,3) colors
	
	variable delta = N/(numatoms-1)
	
	variable m=0
	for(m=0; m<numatoms; m+=1)
		colors[m][] = M_colors[floor(delta*m)][q]
	endfor
	
	variable i=0
	for(i=0; i<size; i+=1)
		variable w=0
		for(w=0; w<numatoms; w+=1)
			if( Atom_Count[i] == w+1 )
				variable d = colors[w][0]
				variable e = colors[w][1]
				variable f = colors[w][2]
				SetDrawEnv fillfgc= (d,e,f), xcoord=prel, ycoord=prel,  linethick= 0.00, save
				l = (x0[i] / 256) - 0.006
				t = (y0[i] / 256) - 0.006
				r = (x0[i] / 256) + 0.006
				b = (y0[i] / 256) + 0.006
				DrawOval l,t,r,b
			endif
		endfor
	endfor
	
	//draw atom number legend
	w=0
	for(w=0; w<numatoms; w+=1)
		d = colors[w][0]
		e = colors[w][1]
		f = colors[w][2]
		
		string h
		sprintf h, "\\K(%d,%d,%d)\\W519\\K(0,0,0) %d atom", d, e, f, w+1
		
		if (w == 0)
			TextBox/C/N=text0/A=MC h
		else
			AppendText/N=text0 h
		endif
	endfor
	
	
	//killwaves M_color, M_WaveStats, colors
	
end


function AddAtomIntSpecific(im, x0, y0, Atom_Count, colors)
	
	wave im, x0, y0, Atom_Count, colors
	NewImage im
	
	variable size = DimSize(Atom_Count, 0)
	
	variable l, t, r, b
	
	wavestats/Q/W Atom_Count
	wave M_WaveStats = $"M_WaveStats"
	variable numatoms = M_WaveStats(12)
	
	variable i=0
	for(i=0; i<size; i+=1)
		variable w=0
		for(w=0; w<numatoms; w+=1)
			if( Atom_Count[i] == w+2 )
				variable d = colors[w][0]
				variable e = colors[w][1]
				variable f = colors[w][2]
				SetDrawEnv fillfgc= (d,e,f), xcoord=prel, ycoord=prel,  linethick= 0.00, save
				l = (x0[i] / 256) - 0.01
				t = (y0[i] / 256) - 0.01
				r = (x0[i] / 256) + 0.01
				b = (y0[i] / 256) + 0.01
				DrawOval l,t,r,b
			endif
		endfor
	endfor
	
	//draw atom number legend
	w=0
	for(w=0; w<numatoms-1; w+=1)
		d = colors[w][0]
		e = colors[w][1]
		f = colors[w][2]
		
		string h
		sprintf h, "\\K(%d,%d,%d)\\W519\\K(%d,%d,%d) %d atom", d, e, f, d, e, f, w+2
		
		if (w == 0)
			TextBox/C/N=text0/A=MC h
		else
			AppendText/N=text0 h
		endif
	endfor
	
end

				