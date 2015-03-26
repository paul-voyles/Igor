#pragma rtGlobals=1 // Use modern global access method.

// Created by Andrew B. Yankovich during Feb 2015.
// Use with EDS Spectrum Image Tools procedure.

// Included are functions to create phantom atomic resolution EDS spectrum images.
// These were written to test EDS denoising code to see if they intorduce any artifacts and work properly. 
// Current state of code does not work well for creating peaks with the inputted peak heights if the peaks are not well separated. This problem has to do with overlap of peaks messing with the predicted peak heights.



// Creates all parts of the Phantom at once. 
// Need to input x, y, and z sizes.
// "offset" is the energy value of first energy channel, and "dispersion" is the energy/channel. 
// "offset" and "dispersion" should both be in the same units of eV of keV.
// "maximum" is the max height of the continuum background intensity.
// "params" give the position, size, and amplitude of the characteristic peaks for 1 element. The wave is a 10xn wave where n is the number of peaks for that element.
// params: row 1 = z position of peak, row 2 = size x, row 3 = size y, row 4 = size z, row 5 = cor xy (0), row 6 = cor xz (0), row 7 = cor yz (0), row 8 = desired aplitude of peak, row 9 = automatically filled amplitude (leave 0), row 10 = fraction of characteristic peak that is the delocalized component (ex: 0.7)
// "positions" gives the spatial positions of one element in the x and y directions. The wave format is a n x 3 wave where n is the number of atom positions, col 1 = x positions, col 2 = y positions, col 3 = occupanct between 0 and 1.
// This function will automatically find the positions and params waves in the root directory. It searchs for waves with "_params" and "_positions" with any prefix before these names. Make sure the matching params and positions waves have the same prefix so they can be automatically grouped by this function. 
function MakePhantomSpectrumImage(size_x, size_y, size_z, offset, dispersion, maximum, a, b, c)
	variable size_x, size_y, size_z, offset, dispersion, maximum, a, b, c
	
	// Create the spectrum image template
	MakeSpectrumImageTemplate(size_x, size_y, size_z, offset, dispersion)
	wave im = $"im"
	
	// Add continuum background and then creates it integrated spectrum
	print "computing continuum background"
	MakeContinuumBackground(im, maximum, a, b, c)
	IntegrateSpectrum(im)
	wave IntegratedSpectrum = $"IntegratedSpectrum"
	duplicate/o IntegratedSpectrum IS_back
	killwaves IntegratedSpectrum
	
	// creates wavelists of param and position waves, sorts them alphabetically
	string params_list = wavelist("*_params",";","")
	string positions_list = wavelist("*_positions",";","")
	params_list = SortList(params_list, ";", 4)
	positions_list = SortList(positions_list, ";", 4)
	if ( ItemsInList(params_list , ";") != ItemsInList(positions_list , ";"))
		print " number of parameters waves does not equal the number of posiitons waves. Exiting..."
		return 0
	endif
	if ( ItemsInList(params_list , ";") == ItemsInList(positions_list , ";"))
		variable num_elements = ItemsInList(params_list , ";")
		print "number of inputed elements = ", num_elements
	endif
	
	string params, positions
	
	// count the number of peaks in all param waves
	variable count = 0
	variable d, e, i, j, k, m, n
	for(d=0; d<num_elements; d+=1)
		params = stringfromlist(d, params_list, ";")
		count = count + dimsize($params,1)
	endfor
	
	// documents all peak positions and delocalized to total peak height ratio 
	variable count_2 = 0
	Make/o/n=(count) all_peak_positions
	Make/o/n=(count) all_peak_ratios
	for(d=0; d<num_elements; d+=1)
		params = stringfromlist(d, params_list, ";")
		wave params_tmp = $params
		for(m=0; m<DimSize(params_tmp,1); m+=1)
			all_peak_positions[count_2] = params_tmp[0][m]
			all_peak_ratios[count_2] = params_tmp[9][m]
			count_2 = count_2 + 1
		endfor
	endfor
	
	// calculates each peak intensity including all the other peaks tails, and then adjusts each peaks amplitude to achieve the desired peak height
	count_2 = 0
	variable occupancy_tmp
	Make/o/n=(DimSize(im,0), DimSize(im,1)) im_tmp = 0
	Make/o/n=(10) coef
	Make/o/n=(count) peak_amp_initial = 0
	for(d=0; d<num_elements; d+=1)
		params = stringfromlist(d, params_list, ";")
		positions = stringfromlist(d, positions_list, ";")
		wave params_tmp = $params
		wave positions_tmp = $positions
		for(m=0; m<DimSize($params,1); m+=1)
			count_2 = count_2 + 1
			coef[0] = params_tmp[1][m]		// sx
			coef[1] = params_tmp[2][m]		// sy
			coef[2] = params_tmp[3][m]		// sz
			coef[3] = params_tmp[4][m]		// cxy
			coef[4] = params_tmp[5][m]		// cxz
			coef[5] = params_tmp[6][m]		// cyz
			coef[8] = params_tmp[0][m]		// z0
			coef[9] = (params_tmp[7][m]-IS_back[params_tmp[0][m]]) * (1- params_tmp[9][m])		//amplitude
			for(e=0; e<DimSize(all_peak_positions,0); e+=1)
				k=all_peak_positions[e]
				for(n=0; n<DimSize(positions_tmp,0); n+=1)
					coef[6] = positions_tmp[n][0]			// x0
					coef[7] = positions_tmp[n][1]			// y0
					occupancy_tmp = positions_tmp[n][2]	// occupancy
					for(i=(positions_tmp[n][0] - 5*params_tmp[1][m]); i<(positions_tmp[n][0] + 5*params_tmp[1][m]); i+=1)
						for(j=(positions_tmp[n][1] - 5*params_tmp[2][m]); j<(positions_tmp[n][1] + 5*params_tmp[2][m]); j+=1)
							if( i >= 0 && j>=0)
								if( i <= 255 && j<=255)
									im_tmp[i][j] = im_tmp[i][j] + (occupancy_tmp * gauss3d(i,j,k, coef) / (DimSize(im, 0)*DimSize(im,1)))
								endif
							endif
						endfor
					endfor
				endfor
				Wavestats/q im_tmp
				peak_amp_initial[e] = peak_amp_initial[e] + V_sum
				im_tmp=0
			endfor
		endfor
	endfor
	killwaves im_tmp
	
	// update aplitudes in param waves with new estimates 
	count_2 = 0
	for(d=0; d<num_elements; d+=1)
		params = stringfromlist(d, params_list, ";")
		wave params_tmp = $params
		for(m=0; m<DimSize(params_tmp,1); m+=1)
			params_tmp[8][m] = (((params_tmp[7][m]-IS_back[params_tmp[0][m]]) * (1- params_tmp[9][m]))^2)/peak_amp_initial[count_2]
			count_2 = count_2 + 1
		endfor
	endfor
	
	IntegrateSpectrum(im)
	wave IntegratedSpectrum = $"IntegratedSpectrum"
	duplicate/o IntegratedSpectrum IS_pre
	killwaves IntegratedSpectrum
	
	// calculate and fill the localized characteistic signal for each peak.
	count_2 = 0
	for(d=0; d<num_elements; d+=1)
		params = stringfromlist(d, params_list, ";")
		positions = stringfromlist(d, positions_list, ";")
		wave params_tmp = $params
		wave positions_tmp = $positions
		for(m=0; m<DimSize(params_tmp,1); m+=1)
			coef[0] = params_tmp[1][m]		// sx
			coef[1] = params_tmp[2][m]		// sy
			coef[2] = params_tmp[3][m]		// sz
			coef[3] = params_tmp[4][m]		// cxy
			coef[4] = params_tmp[5][m]		// cxz
			coef[5] = params_tmp[6][m]		// cyz
			coef[8] = params_tmp[0][m]		// z0
			coef[9] = params_tmp[8][m]		// new amplitude
		
			for(n=0; n<DimSize(positions_tmp,0); n+=1)
				coef[6] = positions_tmp[n][0]			// x0
				coef[7] = positions_tmp[n][1]			// y0
				occupancy_tmp = positions_tmp[n][2]	// occupancy
		
				for(i=positions_tmp[n][0] - 5*params_tmp[1][m]; i<positions_tmp[n][0] + 5*params_tmp[1][m]; i+=1)
					for(j=positions_tmp[n][1] - 5*params_tmp[2][m]; j<positions_tmp[n][1] + 5*params_tmp[2][m]; j+=1)
						for(k=params_tmp[0][m] - 5*params_tmp[3][m]; k<params_tmp[0][m] + 5*params_tmp[3][m]; k+=1)
							if( i >= 0 && j>=0 && k>=0)
								if( i <= 255 && j<=255 && k<=2047)
									im[i][j][k] =  im[i][j][k] + (occupancy_tmp * gauss3d(i,j,k, coef) / (DimSize(im, 0)*DimSize(im,1)))
								endif
							endif
						endfor
					endfor
				endfor
			endfor
			
			string name = "IS_" + num2str(count_2)
			
			IntegrateSpectrum(im)
			wave IntegratedSpectrum = $"IntegratedSpectrum"
			duplicate/o IntegratedSpectrum IS_peak
			IS_peak = IS_peak-IS_pre
			duplicate/o IS_peak $name
			IS_pre = IntegratedSpectrum 
			killwaves IntegratedSpectrum, IS_peak
			
			count_2 = count_2 +1
		endfor
	endfor
	IntegrateSpectrum(im)
	wave IntegratedSpectrum = $"IntegratedSpectrum"
	duplicate/o IntegratedSpectrum IS_back_peak
	killwaves IntegratedSpectrum
	
	// calculate and fill the delocalized characteistic signal for each peak.
	string IS_tmp
	variable int_tmp
	for(d=0; d<count; d+=1)
		IS_tmp = "IS_" + num2str(d)
		wave wv_tmp = $IS_tmp
		for(k=0; k<dimsize(im,2); k+=1)
			if (wv_tmp[k] > 0)
				int_tmp = (all_peak_ratios[d] * (wv_tmp[k]/(1-all_peak_ratios[d])))/(dimsize(im,0)*dimsize(im,1))
				for(i=0; i<dimsize(im,0); i+=1)
					for(j=0; j< dimsize(im,1); j+=1)
						im[i][j][k] =  im[i][j][k] + int_tmp
					endfor
				endfor
			endif
		endfor
	endfor	
	
	IntegrateSpectrum(im)
	wave IntegratedSpectrum = $"IntegratedSpectrum"
	duplicate/o IntegratedSpectrum IS_final
	killwaves IntegratedSpectrum, IS_pre
	
	// displays the integrated intensity components
	Display/K=0  root:IS_final
	AppendToGraph IS_back_peak
	AppendToGraph IS_back
	ModifyGraph mode=7,hbFill=2,rgb(IS_final)=(30583,30583,30583);DelayUpdate
	ModifyGraph rgb(IS_back_peak)=(52428,52428,52428),rgb(IS_back)=(0,0,0)
	
	killwaves coef, peak_amp_initial//, all_peak_positions,  all_peak_ratios
end



// Creates a spectrum image "im" with the given x, y, and z sizes.
// "offset" is the energy value of first energy channel, and "dispersion" is the energy/channel. 
// "offset" and "dispersion" should both be in the same units of eV.
function MakeSpectrumImageTemplate(size_x, size_y, size_z, offset, dispersion)
	variable size_x, size_y, size_z, offset, dispersion
	
	Make/o/n=(size_x, size_y, size_z) im
	im = 0
	setscale/p z, offset, dispersion, im
end



// Adss a continuum background to the spectrum image "im".
// "maximum" is the max height of the continuum background intensity.
// a, b and c are paramters from the fit of the exp background to the function used below for variable "sum_int"
function MakeContinuumBackground(im, maximum, a, b, c)
	wave im
	variable maximum, a, b, c
	
	variable offset = DimOffset(im, 2) 
	variable dispersion = DimDelta(im,2)
	variable x, sum_int
	variable i, j, k
	
	for(i=0; i<DimSize(im,2); i+=1)
		x=i-(abs(offset)/dispersion)
		if( i > (abs(offset)/dispersion))
			sum_int = a*(exp(-b/x))*((x-c)/x)		// background fit function
			for(j=0; j<DimSize(im,0); j+=1)
				for(k=0; k<DimSize(im,1); k+=1)
					im[j][k][i] = (sum_int/(DimSize(im,0)*DimSize(im,1)))
				endfor
			endfor
		endif
	endfor
	// rescale "im" to have a maximum intensity of the given "maximum".
	wavestats/q im
	im = im*maximum/(V_max*DimSize(im,0)*DimSize(im,1))
end



// returns the value of a 1D Gaussian at posiiton x.
// A is the amplitude. x0 is the peak position. sx is the peak width.
function Gaussian1D(A, x, x0, sx)
	variable A, x, x0, sx
	
	return A*exp(-((x-x0)/sx)^2)
end



// returns the value of a 3D gaussian at position x,y,z.
// "coef" wave includes the 3D Gaussian paramters. The size and order of "coef" is indicated below.
function gauss3d(x,y,z,coef)
	variable x,y,z
	wave coef
	
	variable A
	variable x0, y0, z0
	variable sx, sy, sz
	variable cxy, cxz, cyz
	A = coef[9]		// Aplitude
	x0 = coef[6]	// x position of peak
	y0 = coef[7]	// y position of peak
	z0 = coef[8]	// z position of peak
	sx = coef[0]	// size of peak in x direction
	sy = coef[1]	// size of peak in y direction
	sz = coef[2]		// size of peak in z direction
	cxy = coef[3]	// cor xy
	cxz = coef[4]	// cor xz
	cyz = coef[5]	// cor yz
	return A*exp( -1/(2* (-1+cxy^2+cxz^2+cyz^2-2*cxy*cxz*cyz)) * ( (cyz^2-1)*(x-x0)^2/sx^2 + (cxz^2-1)*(y-y0)^2/sy^2 + (cxy^2-1)*(z-z0)^2/sz^2 + (2*cxy*(x-x0)*(y-y0)-2*cxz*cyz*(x-x0)*(y-y0))/(sx*sy) + (2*cxz*(x-x0)*(z-z0)-2*cxy*cyz*(x-x0)*(z-z0))/(sx*sz) + (2*cyz*(y-y0)*(z-z0)-2*cxy*cxz*(y-y0)*(z-z0))/(sy*sz) ) )
end



// This function creates a position wave for multiple nanoparticles.
// "template" in the wave that has all the possible atom positions for the desired element.
// "positions" is a 3 x n wave, where n is the number of nanoparticles in the field of view.
// row 1 of "positions" is the x positions of the NP center in pixels.
// row 2 of "positions" is the y positions of the NP center  in pixels.
// row 3 of "positions" is the radius of the NP in pixels.
// outputed is "NP_positions" where column 1 is the x positions of all atoms, column 2 is the y positions of all atoms, and column 3 is the occupancy of all atoms.
// occupancy is between 0.1 and 1, where 1 is near the NP center and 0.1 is near the NP surface.
function NanoparticleGenerator(positions, template)
	wave positions, template
	
	Make/o/n=(0,3) NP_positions
	
	variable x_tmp, y_tmp, radius_pos, x_pos, y_pos, dist_sq
	variable i, j
	for(i=0; i<DimSize(positions,1); i+=1)
		x_pos = positions[0][i]
		y_pos = positions[1][i]
		radius_pos = positions[2][i]
		for(j=0; j<DimSize(template,0); j+=1)
			x_tmp = template[j][0]
			y_tmp = template[j][1]
			dist_sq = ((x_tmp - x_pos)^2) + ((y_tmp - y_pos)^2)
			if( dist_sq <= (radius_pos^2))
				InsertPoints/M=0 0, 1, NP_positions
				NP_positions[0][0] = x_tmp
				NP_positions[0][1] = y_tmp
				
				if(0  <= dist_sq/(radius_pos^2) && dist_sq/(radius_pos^2) < 0.1)
					NP_positions[0][2] = 1
				endif
				if(0.1 <= dist_sq/(radius_pos^2) && dist_sq/(radius_pos^2) < 0.2)
					NP_positions[0][2] = 0.9
				endif
				if(0.2  <= dist_sq/(radius_pos^2) && dist_sq/(radius_pos^2) <  0.3)
					NP_positions[0][2] = 0.8
				endif
				if(0.3  <= dist_sq/(radius_pos^2) && dist_sq/(radius_pos^2) <  0.4)
					NP_positions[0][2] = 0.7
				endif
				if(0.4 <= dist_sq/(radius_pos^2) && dist_sq/(radius_pos^2) < 0.5)
					NP_positions[0][2] = 0.6
				endif
				if(0.5  <= dist_sq/(radius_pos^2) && dist_sq/(radius_pos^2) <  0.6)
					NP_positions[0][2] = 0.5
				endif
				if(0.6  <= dist_sq/(radius_pos^2) && dist_sq/(radius_pos^2) <  0.7)
					NP_positions[0][2] = 0.4
				endif
				if(0.7  <= dist_sq/(radius_pos^2) && dist_sq/(radius_pos^2) < 0.8)
					NP_positions[0][2] = 0.3
				endif
				if(0.8  <= dist_sq/(radius_pos^2) && dist_sq/(radius_pos^2) <  0.9)
					NP_positions[0][2] = 0.2
				endif
				if(0.9  <= dist_sq/(radius_pos^2) && dist_sq/(radius_pos^2) <= 1)
					NP_positions[0][2] = 0.1
				endif
			endif
		endfor
	endfor
	
	Display NP_positions[][1] vs NP_positions[][0]
	SetAxis left 0,256;DelayUpdate
	SetAxis bottom 0,256
	ModifyGraph mode=3,marker=19
end


// creates a new positions wave that is a collection of randomly picked positions from "host_positions.
// percent is the percent of the host positions you want to include. percent is between 0 and 100.
// output is "random_positions"
function RandomOccupancyeGenerator(host_positions, percent)
	wave host_positions
	variable percent
	
	make/o/n=(round((dimsize(host_positions,0))*(percent/100)),3) random_positions
	duplicate/o host_positions host_pos_tmp
	
	variable random_num
	variable i
	for (i=0; i<DimSize(random_positions,0); i+=1)
		
		random_num = floor(abs(enoise(dimsize(host_pos_tmp,0))))
		random_positions[i][0] = host_pos_tmp[random_num][0]
		random_positions[i][1] = host_pos_tmp[random_num][1]
		random_positions[i][2] = 1
		
		deletepoints random_num, 1, host_pos_tmp
	endfor

	Display random_positions[][1] vs random_positions[][0]
	SetAxis left 0,256;DelayUpdate
	SetAxis bottom 0,256
	ModifyGraph mode=3,marker=19
	
	killwaves host_pos_tmp
end


// This function decreasdes the host_position occupancies based on the occupancies of the impurity_positions.
// decrease variable is the maximum decrease in the host occupancies, between 0 and 1.
function DecreaseHostOccupancy(host_positions, impurity_positions, decrease)
	wave host_positions, impurity_positions
	variable decrease
	
	variable host_x, host_y, imp_x, imp_y, imp_occ
	
	variable i, j
	for (i=0; i<DimSize(host_positions,0); i+=1)
		host_x = host_positions[i][0]
		host_y = host_positions[i][1]
		for (j=0; j<DimSize(impurity_positions,0); j+=1)
			imp_x = impurity_positions[j][0]
			imp_y = impurity_positions[j][1]
			
			if(imp_x == host_x && imp_y == host_y)
				host_positions[i][2]  = host_positions[i][2] - (decrease*impurity_positions[j][2])
			endif
			
		endfor
	endfor
	
end

// this function resets the occupancies in the "positions" wave to whatever value you input as "num".
function ResetOccupancy(positions, num)
	wave positions
	variable num
	
	variable i
	for (i=0; i<DimSize(positions,0); i+=1)
		positions[i][3] = num
	endfor
end
