#pragma rtGlobals=1		// Use modern global access method.
// Created/updated on 12/8/14.  ABY
//
// These functions are for calculating average lattice paramters, creating perfect 2D lattices of atomic column positions from those 
// parameters, and fitting the perfect lattices to experimental lattices by minimizing the rms between all the points.
//
// 12-09-14 Added checking for a few more assumed pre-existing waves; modified SortGrid to not return fitted atom positions if
// matching grid atoms are not found.  pmv
// 12-09-14 Modified BasisParameters and CreateGridWithBasis to allow any number of basis atoms, not just one.  pmv
// 12-10-14 bug fix for number of atoms in CreateGridWithBasis. pmv
// 12-10-14 change CreateGridWithBasis to shift the created grid to the center of the atom positions, then shift to the nearest matching
// grid / atom pair.  pmv


// This functions finds the average lattice parameters(a,b) of a two-dimensional grid of atomic columns positions with only 1 atom basis
// listed in x0 and y0. 
// The order of the atoms positions in the x0, y0 list of atom positions can be in any order, but x0[i] must correspond to y0[i].
// a_space and b_space are guesses of what the seperations are for the lattice paramters. space_delta is the +/- separation error window
// for the selection process. 
// a_angle and b_angle are guesses for the angles of the lattice parameters. angle_delta is the +/- angle error window for the selection process. 
// horizantal is 0 degrees. straight down is -90 degrees.
// +/- 10% seems to work for the space_delta and angle_delta.
// outputed is ax, ay, bx, and by, which are the x and y components of all the found interatomic separations.
// outputed is a_param and b_param. row 0 = average ax or bx. row 1 = average ay or by, row 2 = average angle.
// the average values and number of separations found are printed in the history. 
// lattice parameters and angles must be in positive x and positive y; may mean angle is negative
function GridParameters(x0, y0, a_space, b_space, space_delta, a_angle, b_angle, angle_delta)

	wave x0, y0
	variable a_space, b_space, space_delta, a_angle, b_angle, angle_delta
	
	variable num_peaks = DimSize(x0,0)
	Make/O/N=0 by, bx, ay, ax
	
	// find and save the a spacings
	variable i, j, xm, ym, d, angle
	i=0
	j=0
	for(i=0; i<num_peaks; i+=1)
		for(j=0; j<num_peaks; j+=1)
			xm = (x0(j) - x0(i))
			ym = (y0(j) - y0(i))
			d = sqrt((xm^2) + (ym^2))
			angle = atan(ym/xm) * 180 / Pi
			
			if(d > (a_space - space_delta) && d < (a_space + space_delta))
				if(angle > (a_angle - angle_delta) && angle < (a_angle + angle_delta))
					if(xm > 0)
						InsertPoints 0, 1, ax
						InsertPoints 0, 1, ay
						ax[x2pnt(ax,0)] = xm
						ay[x2pnt(ay,0)] = ym
					endif
				endif
			endif
		endfor
	endfor

	// calculate and save the a parameters
	wavestats/Q/W ax
	wave M_WaveStats = $"M_WaveStats"
	variable ax_average = M_WaveStats(3)
	variable ax_stdev = M_WaveStats(4)
	wavestats/Q/W ay
	wave M_WaveStats = $"M_WaveStats"
	variable ay_average = M_WaveStats(3)
	variable ay_stdev = M_WaveStats(4)
	make/o/n=3 a_param
	a_param[0] = ax_average
	a_param[1] = ay_average
	a_param[2] = atan(ay_average/ax_average) * 180 / Pi
	print "a lattice parameter in ( x , y , angle ) = (", a_param[0], ",", a_param[1], ",", a_param[2], ") and found", dimsize(ax,0), "seperations."
	
	//find and save the b spacings
	i=0
	j=0
	for(i=0; i<num_peaks; i+=1)
		for(j=0; j<num_peaks; j+=1)
			xm = (x0(j) - x0(i))
			ym = (y0(j) - y0(i))
			d = sqrt((xm^2) + (ym^2))
			angle = atan(ym/xm) * 180 / Pi
			
			if(d > (b_space - space_delta) && d < (b_space + space_delta))
				if(angle > (b_angle - angle_delta) && angle < (b_angle + angle_delta))
					if(ym > 0)
						InsertPoints 0, 1, bx
						InsertPoints 0, 1, by
						bx[x2pnt(bx,0)] = xm
						by[x2pnt(by,0)] = ym
					endif
				endif
			endif
		endfor
	endfor

	// calculate and save the b parameters
	wavestats/Q/W bx
	wave M_WaveStats = $"M_WaveStats"
	variable bx_average = M_WaveStats(3)
	variable bx_stdev = M_WaveStats(4)
	wavestats/Q/W by
	wave M_WaveStats = $"M_WaveStats"
	variable by_average = M_WaveStats(3)
	variable by_stdev = M_WaveStats(4)
	make/o/n=3 b_param
	b_param[0] = bx_average
	b_param[1] = by_average
	b_param[2] = atan(by_average/bx_average) * 180 / Pi
	print "b lattice parameter in ( x , y , angle ) = (", b_param[0], ",", b_param[1], ",", b_param[2], ") and found", dimsize(bx,0), "seperations."
	
	killwaves M_WaveStats
end


// This functions finds a single average lattice separation (c) of a two-dimensional grid of atomic columns positions listed in x0 and y0, which can act as a second atom for a 2 atom basis 2D unit cell. 
// The order of the atoms positions in the x0, y0 list of atom positions can be in any order, but x0[i] must correspond to y0[i].
// c_space is the guess for the lattice paramter seperation. space_delta is the +/- separation error window for the selection process. 
// c_angle is the guess for the angle of the lattice parameters. angle_delta is the +/- angle error window for the selection process. 
// horizantal is 0 degrees. straight down is -90 degrees.
// +/- 10% seems to work for the space_delta and angle_delta.
// outputed is cx and cy, which are the x and y components of all the found interatomic separations.
// outputed is c_param. row 0 = average cx. row 1 = average cy, row 2 = average angle.
// the average values and number of separations found are printed in the history. 
function BasisParameters(x0, y0, c_space, space_delta, c_angle, angle_delta)
	
	wave x0, y0
	variable c_space, space_delta, c_angle, angle_delta
	
	variable num_peaks = DimSize(x0,0)
	Make/O/N=0 cy, cx

	//find and save the c spacings
	variable i, j, xm, ym, d, angle
	i=0
	j=0
	for(i=0; i<num_peaks; i+=1)
		for(j=0; j<num_peaks; j+=1)
			xm = (x0(j) - x0(i))
			ym = (y0(j) - y0(i))
			d = sqrt((xm^2) + (ym^2))
			angle = atan(ym/xm) * 180 / Pi
			
			if(d > (c_space - space_delta) && d < (c_space + space_delta))
				if(angle > (c_angle - angle_delta) && angle < (c_angle + angle_delta))
					if(ym > 0)
						InsertPoints 0, 1, cx
						InsertPoints 0, 1, cy
						cx[x2pnt(cx,0)] = xm
						cy[x2pnt(cy,0)] = ym
					endif
				endif
			endif
		endfor
	endfor

	if(numpnts(cx) > 0)
	
		// calculate and save the c parameters
		wavestats/Q/W cx
		wave M_WaveStats = $"M_WaveStats"
		variable cx_average = M_WaveStats(3)
		variable cx_stdev = M_WaveStats(4)
		wavestats/Q/W cy
		wave M_WaveStats = $"M_WaveStats"
		variable cy_average = M_WaveStats(3)
		variable cy_stdev = M_WaveStats(4)

		wave c_param = $"c_param"
		if(!WaveExists(c_param))
			make/o/n=(3, 1) c_param
		else
			InsertPoints/M=1 0, 1, c_param
		endif

		c_param[0][0] = cx_average
		c_param[1][0] = cy_average
		c_param[2][0] = atan(cy_average/cx_average) * 180 / Pi
		print "c lattice basis in ( x , y , angle ) = (", c_param[0], ",", c_param[1], ",", c_param[2], ") and found", dimsize(cx,0), "seperations."
		killwaves M_WaveStats
	
	else
		printf "No matching interatomic distances found within tolerances.\r"
	endif
	
end


// print the vector in (length, angle) parameter space from point (x1, x2) to (y1, y2).  Use with hcsr() and vcsr() functions to define
// vectors on an image - EstimateVector(hcsr(A), hcsr(B), vcsr(A), vcsr(B))
function EstimateVector(x1, x2, y1, y2)
	variable x1, x2, y1, y2

	variable xm, ym, d, angle
	xm = (x1 - x2)
	ym = (y1 - y2)
	d = sqrt((xm^2) + (ym^2))
	angle = atan(ym/xm) * 180 / Pi

	printf "Vector from point (%g, %g) to (%g, %g) is d = %g, angle = %g\r", x1, x2, y1, y2, d, angle

end

// This function creates a perfect 2D lattice of positions from what is in a_param and b_param. 
// This function only creates a 2D lattice with a 1 atom basis.
// na (or nb) is the number of positions the lattice is extended in the a direction (or b direction).
// This function creates waves called x_lat and y_lat which are the positions of the perfect lattice.
// This function automatically initializes(shifts) the position of the gird so that the atom position closest to (0,0) in (x_lat, y_lat) is equal to the atom position closest to (0,0) in (x0, y0).
// This function also automatically sorts x_lat and y_lat so that (x_lat[i], y_lat[i]) corresponds to the same atom position (i) in x0, y0. These soted outputs are found in x_lat_sort, y_lat_sort.
function CreateGrid(na, nb)

	variable na, nb
	variable numatoms = na * nb

	//open parameter waves
	wave tmp_a_param = $"a_param"
	wave tmp_b_param = $"b_param"
		if(!waveexists(tmp_a_param))
			print "You are an idoit, a_param doesn't exist.\r"
			return 0
		endif
		if(!waveexists(tmp_b_param))
			print "You are an idoit, b_param doesn't exist.\r"
			return 0
		endif
	
	Make/o/n=(numatoms) y_lat, x_lat
	
	//populate lattice with a and b parameters
	variable i, j, k
	k=0
	for(i=0; i<na; i+=1)
		for(j=0; j<nb; j+=1)
			x_lat[k] = (i * tmp_a_param[0]) + (j * tmp_b_param[0])
			y_lat[k] = (i * tmp_a_param[1]) + (j * tmp_b_param[1])
			k = k+1
		endfor
	endfor
	
	wave x0 = $"x0"
	wave y0 = $"y0"
	
	//calculate the position in x0, y0 that is closest to the origin to which x_lat and y_lat will be initialized
	variable x_initial, y_initial =10000000
	variable initial_distance = sqrt((x_initial^2) + (y_initial^2))
	variable distance
	i=0
	for(i=0; i<DimSize(x0,0); i+=1)
		distance = sqrt((x0[i]^2) + (y0[i]^2))
		if(distance < initial_distance)
			x_initial = x0[i]
			y_initial = y0[i]
			initial_distance = distance
		endif
	endfor
	
	//Initialize x_lat and y_lat
	InitializeGrid(x_lat, y_lat, x_initial, y_initial)
	//Sort x_lat and y_lat waves to be in same aroder as x0 and y0. outputs are x_lat_sort and y_lat_sort
	SortGrid(x_lat, y_lat, x0, y0, 2)
end


// This function creates a perfect 2D lattice of positions from what is in a_param, b_param, and c_param. 
// This function creates a 2D lattice with a 2 atom basis. The first atom basis is determined by a_param and b_param. The second atom in the basis is determined by c_param.
// na (or nb) is the number of positions the lattice is extended in the a direction (or b direction).
// This function creates waves called x_lat and y_lat which are the positions of the perfect lattice.
// This function automatically initializes(shifts) the position of the gird so that the atom position closest to (0,0) in (x_lat, y_lat) is equal to the atom position closest to (0,0) in (x0, y0).
// This function also automatically sorts x_lat and y_lat so that (x_lat[i], y_lat[i]) corresponds to the same atom position (i) in x0, y0. These soted outputs are found in x_lat_sort, y_lat_sort.
function CreateGridWithBasis(na, nb)
	variable na, nb

	//open parameter waves
	wave tmp_a_param = $"a_param"
	wave tmp_b_param = $"b_param"
	wave tmp_c_param = $"c_param"
	wave x0 = $"x0"
	wave y0 = $"y0"

	if(!waveexists(tmp_a_param))
		print "You are an idoit, a_param doesn't exist.\r"
		return 0
	endif
	if(!waveexists(tmp_b_param))
		print "You are an idoit, b_param doesn't exist.\r"
		return 0
	endif
	if(!waveexists(tmp_c_param))
		print "You are an idoit, c_param doesn't exist.\r"
		return 0
	endif
	if(!waveexists(x0))
		print "Cannot find wave x0.  Exiting.\r"
		return 0
	endif
	if(!waveexists(y0))
		print "Cannot find wave y0.  Exiting.\r"
		return 0
	endif
	
	
	variable numbasis = DimSize(tmp_c_param, 1)
	variable numatoms = na * nb * (numbasis+1)

	Make/o/n=(numatoms) y_lat, x_lat
	x_lat = 0
	y_lat = 0
	
	//populate lattice with a and b parameters with the c basis
	variable i, j, k, l
	l = 0
	for(i=0; i<na; i+=1)
		for(j=0; j<nb; j+=1)
			x_lat[l] = (i * tmp_a_param[0]) + (j * tmp_b_param[0])
			y_lat[l] = (i * tmp_a_param[1]) + (j * tmp_b_param[1])
			l += 1
			for(k=0; k<numbasis; k+=1)
				x_lat[l] = x_lat[l-(k+1)] + tmp_c_param[0][k]
				y_lat[l] = y_lat[l-(k+1)] + tmp_c_param[1][k]
				l+=1
			endfor
		endfor
	endfor
	
	
	// shift the calculated grid to cover the same range of coordinates as the fitted atom positions
	// then to exactly overlap with one particular atom
	variable x0_mid, y0_mid
	variable gridx_mid, gridy_mid
	
	wavestats/q/M=0 x0
	x0_mid = (V_max + V_min)/2
	wavestats/q/M=0 y0
	y0_mid = (V_max + V_min)/2
	wavestats/q/M=0 x_lat
	gridx_mid = (V_max + V_min)/2
	wavestats/q/M=0 y_lat
	gridy_mid = (V_max + V_min)/2
	
	// rough center the grid on the atom positions
	x_lat = x_lat + (x0_mid - gridx_mid)
	y_lat = y_lat + (y0_mid - gridy_mid)

	// now find the closest grid / atom pair
	Make/O/N=(numpnts(x_lat), numpnts(x0)) dist_t
	dist_t = (x_lat[p] - x0[q])^2 + (y_lat[p] - y0[q])^2
	wavestats/q/M=0 dist_t
	variable shift_x = x0[V_minColLoc] - x_lat[V_minRowLoc]
	variable shift_y = y0[V_minColLoc] - y_lat[V_minRowLoc]
	x_lat = x_lat + shift_x
	y_lat = y_lat + shift_y
		
	// Sort x_lat and y_lat waves to be in same order as x0 and y0. outputs are x_lat_sort and y_lat_sort
	SortGrid(x_lat, y_lat, x0, y0, 2)
	
	Killwaves dist_t
	
end


// grid with starting point.  (sx ,sy) is the origin for the lattice.  However, it runs from sa to sa+na unit cells
// along the a parameter and sb to sb+nb along the b parameter.  Uses global waves a_param, b_param, and 
// c_param.
function CreateGridFromPoint(sx, sy, sa, sb, na, nb)
	variable sx, sy, sa, sb, na, nb
	
	//open parameter waves
	wave tmp_a_param = $"a_param"
	wave tmp_b_param = $"b_param"
	wave tmp_c_param = $"c_param"
	wave x0 = $"x0"
	wave y0 = $"y0"

	if(!waveexists(tmp_a_param))
		print "You are an idoit, a_param doesn't exist.\r"
		return 0
	endif
	if(!waveexists(tmp_b_param))
		print "You are an idoit, b_param doesn't exist.\r"
		return 0
	endif
	if(!waveexists(tmp_c_param))
		print "You are an idoit, c_param doesn't exist.\r"
		return 0
	endif
	if(!waveexists(x0))
		print "Cannot find wave x0.  Exiting.\r"
		return 0
	endif
	if(!waveexists(y0))
		print "Cannot find wave y0.  Exiting.\r"
		return 0
	endif
	
	
	variable numbasis = DimSize(tmp_c_param, 1)
	variable numatoms = na * nb * (numbasis+1)

	Make/o/n=(numatoms) y_lat, x_lat
	x_lat = 0
	y_lat = 0
	
	// adjust starting points for non-zero sa, sb
	sx += sa*tmp_a_param[0] + sb*tmp_b_param[0]
	sy += sa*tmp_a_param[1] + sb*tmp_b_param[1]
	
	//populate lattice with a and b parameters with the c basis
	variable i, j, k, l
	l = 0
	for(i=0; i<na; i+=1)
		for(j=0; j<nb; j+=1)
			x_lat[l] = sx + (i * tmp_a_param[0]) + (j * tmp_b_param[0])
			y_lat[l] = sy + (i * tmp_a_param[1]) + (j * tmp_b_param[1])
			l += 1
			for(k=0; k<numbasis; k+=1)
				x_lat[l] = x_lat[l-(k+1)] + tmp_c_param[0][k]
				y_lat[l] = y_lat[l-(k+1)] + tmp_c_param[1][k]
				l+=1
			endfor
		endfor
	endfor
	
	// Sort x_lat and y_lat waves to be in same order as x0 and y0. outputs are x_lat_sort and y_lat_sort
	// allowable error is set to 20% of the shorter lattice parameter.
	variable err = 0.2*min( sqrt(tmp_a_param[0]^2 + tmp_a_param[1]^2), sqrt(tmp_b_param[0]^2 + tmp_b_param[1]^2) )
	SortGrid(x_lat, y_lat, x0, y0, err)	
	
end

// This function shifts the entire prefect lattice in x_lat, y_lat by shifting (x_lat[0], y_lat[0]) to the position (x_initial, y_initial)
function InitializeGrid(x_lat, y_lat, x_initial, y_initial)
	
	wave x_lat, y_lat
	variable x_initial, y_initial
	
	variable xdiff, ydiff
	xdiff = x_initial - x_lat[0]
	ydiff = y_initial - y_lat[0]
	x_lat = x_lat + xdiff
	y_lat = y_lat + ydiff
end


// This function sorts x_lat and y_lat so that (x_lat[i], y_lat[i]) corresponds to the same atom position (i) in x0, y0. These soted outputs are found in x_lat_sort, y_lat_sort.
// error is the +/- range around atom poisitions that the code identify 2 positions as the same position.
function SortGrid(x_lat, y_lat, x0, y0, error)

	wave x_lat, y_lat, x0, y0
	variable error
	
	variable npts = numpnts(x0)
	variable nlat = numpnts(x_lat)
	make/o/n=(npts) x_lat_sort, y_lat_sort
	x_lat_sort = NaN
	y_lat_sort = NaN
	
	variable count = 0
	variable i, j = 0
	for(i=0; i<npts; i+=1)
		for(j=0; j<nlat; j+=1)
			if(x_lat[j] > (x0[i]-error) && x_lat[j] < (x0[i]+error))
				if(y_lat[j] > (y0[i]-error) && y_lat[j] < (y0[i]+error))
					x_lat_sort[i] = x_lat[j]
					y_lat_sort[i] = y_lat[j]
					count = count + 1
				endif
			endif
		endfor
	endfor
	print "sorting process correctly found ", count, "of ", npts, "atoms columns"
end

function TruncateGrid(xs, xe, ys, ye, xl, yl)
	variable xe, xs, ye, ys
	wave xl, yl
	
	Duplicate/O xl x_lat_trunc
	Duplicate/O yl y_lat_trunc
	
	Sort x_lat_trunc, x_lat_trunc, y_lat_trunc

	FindLevel/Q/P x_lat_trunc, xs
	if(!V_Flag)
		DeletePoints 0, V_LevelX-1, x_lat_trunc, y_lat_trunc
	endif
	
	FindLevel/Q/P x_lat_trunc, xe
	if(!V_Flag)
		DeletePoints V_LevelX, numpnts(x_lat_trunc), x_lat_trunc, y_lat_trunc
	endif
	
	Sort y_lat_trunc, x_lat_trunc, y_lat_trunc

	FindLevel/Q/P y_lat_trunc, ys
	if(!V_Flag)
		DeletePoints 0, V_LevelX-1, x_lat_trunc, y_lat_trunc
	endif
	
	FindLevel/Q/P y_lat_trunc, ye
	if(!V_Flag)
		DeletePoints V_LevelX, numpnts(x_lat_trunc), x_lat_trunc, y_lat_trunc
	endif
	
	
end

// This function finds the position of x_lat_sort and y_lat_sort that minimizes the rms between the positions in (x_lat_sort, y_lat_sort) and (x0, y0).
// range is the number of pixels in x and y you want to scan over (range x range) = area scanned.
// nmpnts is the number of points you want split that range up into. (nmpnts x nmpnts) = number of points in area.
// Outputs wave rms_im which shows the rms value for each position in the scan. 
// The atom column position for the first atom in (x_lat,y_lat) is saved for each point in rms_im and can be found in x_initial_im and y_initial_im.
// The x_lat_sort and y_lat_sort waves are shifted (or initialized) to the position that minimizes the rms between it and (x0, y0) and saved in x_lat_final and y_lat_final. 
// The minimum rms value and the initialized (x,y) position of the first atom position in the list are printed in the history.
function RMS(x_lat_sort, y_lat_sort, x0, y0, range, nmpnts)
	
	wave x_lat_sort, y_lat_sort, x0, y0
	variable range, nmpnts
	
	duplicate/o x_lat_sort, x_lat_tmp
	duplicate/o y_lat_sort, y_lat_tmp
	duplicate/o x_lat_sort, diff
	
	variable step = range/nmpnts
	
	make/o/n=(nmpnts,nmpnts) rms_im
	make/o/n=(nmpnts,nmpnts) x_initial_im
	make/o/n=(nmpnts,nmpnts) y_initial_im
	
	//shift perfect lattice everywhere within range x range and calculate rms at each position. Saves rms and position information.
	variable i=0
	variable j=0
	for(i=0; i<nmpnts; i+=1)
		for(j=0; j<nmpnts; j+=1)
			x_lat_tmp = x_lat_sort - ((nmpnts/2)*step)+(i*step)
			y_lat_tmp = y_lat_sort - ((nmpnts/2)*step) + (j*step)
			diff = sqrt(((x_lat_tmp- x0)^2) + ((y_lat_tmp-y0)^2))
			Wavestats/q diff 
			rms_im[i][j] = V_rms
			x_initial_im[i][j] = x_lat_tmp[0]
			y_initial_im[i][j] = y_lat_tmp[0]
		endfor
	endfor
	
	//Find minimum rms from the saved values.
	variable rms_min = 1000000
	variable x_pos_min, y_pos_min
	i=0
	j=0
	for(i=0; i<nmpnts; i+=1)
		for(j=0; j<nmpnts; j+=1)
			if(rms_im[i][j] < rms_min)
				rms_min = rms_im[i][j]
				x_pos_min = x_initial_im[i][j]
				y_pos_min = y_initial_im[i][j]
				
			endif
		endfor
	endfor
	print "minimum rms = ", rms_min, "occurs when first atom column is located at ( x , y ) = (", x_pos_min, ",", y_pos_min, ")"
	print "x_lat_final and y_lat_final have been initialized to this position"
	duplicate/o x_lat_sort, x_lat_final
	duplicate/o y_lat_sort, y_lat_final
	//Initialize x_lat_final, y_lat_final to the found minimum rms value.
	InitializeGrid(x_lat_final, y_lat_final, x_pos_min, y_pos_min)
	
	killwaves diff, x_lat_tmp, y_lat_tmp
end




// calculates the displacement vector for every points in (x_lat,y_lat) to the nearest
// point in (x0,y0).  (x_lat,y_lat) and (x0,y0) do not have to be the same size.
function GridDisplacement(x_lat, y_lat, x0, y0)
	wave x_lat, y_lat, x0, y0
	
	duplicate/o x0 dx, dy
	
	if(numpnts(x_lat) != numpnts(y_lat))
		printf "Must have the same number of x and y coordinates for grid in x_lat and y_lat.\r"
		return 0
	endif
	if(numpnts(x0) != numpnts(y0))
		printf "Must have the same number of x and y coordinates for atom location in x0 and y0.\r"
		return 0
	endif
	
	variable npts1 = numpnts(x_lat)
	variable npts2 = numpnts(x0)
	
	make/n=(npts1, npts2)/o pdist
	
	 pdist = (x_lat[p] - x0[q])^2 + (y_lat[p] - y0[q])^2

	variable i
	for(i=0; i<npts2; i+=1)
		matrixop/o one_dist = col(pdist, i)
		 wavestats/m=1/q one_dist
		 dx[i] = x0[i] - x_lat[V_minloc]
		 dy[i] = y0[i] - y_lat[V_minloc]
	
	endfor
	
	Killwaves pdist, one_dist
end


// takes (dx, dy) waves from GridDisplacement and makes a wave for an Igor ArrowPlot.
// sc is the scale - arrow plot uses pixels, not graph units, so this is a scale factor between
// the physical distance in (dx, dy) and the displaced distance on the graph.  To make an arrow
// plot, first make a graph of the points, e.g. display gy vs gx.  Then do
//
// Modifygraph arrowmarker(gy) = {arrow_dat, 1, 3, 3, 0}
//
// See Igor help in the PDF manual (search index for vector plot) for more info.			
function ArrowPlotData(dx, dy, sc)
	wave dx, dy
	variable sc
	
	make/o/n=(numpnts(dx), 2) arrow_dat
	arrow_dat[][0] = sc*sqrt(dx[p]^2 + dy[p]^2)
	arrow_dat[][1] = atan2(dy[p], dx[p])+Pi
	
end


function AngleAvgDisplErr(dx, dy, sigma_x0, sigma_y0)
	
	wave dx, dy, sigma_x0, sigma_y0
	
	variable num
	num = numpnts(dx)
	
	Duplicate/o dx weighted_sigma
	Duplicate/o dx abs_dx
	Duplicate/o dy abs_dy
	
	abs_dx = abs(dx)
	abs_dy = abs(dy)
	
	variable i
	for(i=0; i<num; i+=1)
		weighted_sigma[i] = ((abs_dx[i] * sigma_x0[i]) + (abs_dy[i] * sigma_y0[i])) / (abs_dx[i] + abs_dy[i])
	endfor
	
	killwaves abs_dx, abs_dy
end