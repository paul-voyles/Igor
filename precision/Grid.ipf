#pragma rtGlobals=1		// Use modern global access method.


function GridParameters(a_space, b_space, space_delta, a_angle, b_angle, angle_delta)

	variable a_space, b_space, space_delta, a_angle, b_angle, angle_delta
	
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
	
	Make/O/N=0 ax, ay, bx, by

	variable i, j, xm, ym, d, angle
		for(i=0; i<num_peaks-1; i+=1)
			for(j=i+1; j<num_peaks; j+=1)
				xm = abs(tmp_x0(j) - tmp_x0(i))
				ym = abs(tmp_y0(j) - tmp_y0(i))
				d = sqrt((xm^2) + (ym^2))
				angle = atan(ym/xm) * 180 / Pi
				
				if(d > (a_space - space_delta) && d < (a_space + space_delta))
					if(angle > (a_angle - angle_delta) && angle < (a_angle + angle_delta))
						InsertPoints 0, 1, ax
						InsertPoints 0, 1, ay
						ax[x2pnt(ax,0)] = abs(tmp_x0(j) - tmp_x0(i))
						ay[x2pnt(ay,0)] = abs(tmp_y0(j) - tmp_y0(i))
					endif
				endif
			endfor
		endfor

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
	
	
	
		for(i=0; i<num_peaks-1; i+=1)
			for(j=i+1; j<num_peaks; j+=1)
				xm = abs(tmp_x0(j) - tmp_x0(i))
				ym = abs(tmp_y0(j) - tmp_y0(i))
				d = sqrt((xm^2) + (ym^2))
				angle = atan(ym/xm) * 180 / Pi
				
				if(d > (b_space - space_delta) && d < (b_space + space_delta))
					if(angle > (b_angle - angle_delta) && angle < (b_angle + angle_delta))
						InsertPoints 0, 1, bx
						InsertPoints 0, 1, by
						bx[x2pnt(bx,0)] = abs(tmp_x0(j) - tmp_x0(i))
						by[x2pnt(by,0)] = abs(tmp_y0(j) - tmp_y0(i))
					endif
				endif
			endfor
		endfor

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
	
	
	//killwaves M_WaveStats
end



function CreateGrid(na, nb)

	variable na, nb
	variable numatoms = na * nb

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
	
	Make/o/n=(numatoms) x_lat, y_lat
	
	variable i, j, k
	k=0
	for(i=0; i<na; i+=1)
			for(j=0; j<nb; j+=1)
				x_lat[k] = (i * tmp_a_param(0)) + (j * tmp_b_param(0))
				y_lat[k] = (i * tmp_a_param(1)) + (j * tmp_b_param(1))
				k = k+1
			endfor
	endfor
	
end


function InitializeGrid(x_initial, y_initial, x_lat, y_lat)

	variable x_initial, y_initial
	wave x_lat, y_lat
	
	x_lat = x_lat + x_initial
	y_lat = y_lat + y_initial
	
end


function RMSX(range, nmpnts, x_lat, x0)

	variable range, nmpnts
	wave x_lat, x0
	
	variable npts = numpnts(x_lat)
	make/n=(npts)/o xdiff 
	
	make/n=(nmpnts-1)/o xrms, x_lat0
	
	variable shift
	shift = (x0[0] - (range/2)) - x_lat[0]
	x_lat = x_lat + shift
	
	variable i
	for(i=0; i<nmpnts-1; i+=1)
		xdiff = abs(x_lat - x0)
		
		Wavestats/q/w xdiff 
		wave M_WaveStats = $"M_WaveStats"
		xrms[i] = M_WaveStats(5)
		x_lat0[i] = x_lat[0]
		
		x_lat = x_lat + range/nmpnts
	endfor
	
	
end


function RMSY(range, nmpnts, y_lat, y0)

	variable range, nmpnts
	wave y_lat, y0
	
	variable npts = numpnts(y_lat)
	make/n=(npts)/o ydiff 
	
	make/n=(nmpnts-1)/o yrms, y_lat0
	
	variable shift
	shift = (y0[0] - (range/2)) - y_lat[0]
	y_lat = y_lat + shift
	
	variable i
	for(i=0; i<nmpnts-1; i+=1)
		ydiff = abs(y_lat - y0)
		
		Wavestats/q/w ydiff 
		wave M_WaveStats = $"M_WaveStats"
		yrms[i] = M_WaveStats(5)
		y_lat0[i] = y_lat[0]
		
		y_lat = y_lat + range/nmpnts
	endfor
	
	
end

// calculates the displacement vector for every points in g2 to the nearest
// point in g1.  g2 and g1 do not have to be the same size.
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