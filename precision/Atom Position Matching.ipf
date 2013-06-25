#pragma rtGlobals=1		// Use modern global access method.

// calculates the displacement vector for every points in g2 to the nearest
// point in g1.  g2 and g1 do not have to be the same size.
function GridDisplacement(g1x, g1y, g2x, g2y)
	wave g1x, g1y, g2x, g2y
	
	duplicate/o g2x dx, dy
	
	if(numpnts(g1x) != numpnts(g1y))
		printf "Must have the same number of x and y coordinates for grid 1.\r"
		return 0
	endif
	if(numpnts(g2x) != numpnts(g2y))
		printf "Must have the same number of x and y coordinates for grid 2.\r"
		return 0
	endif
	
	variable npts1 = numpnts(g1x)
	variable npts2 = numpnts(g2x)
	
	make/n=(npts1, npts2)/o pdist
	
	 pdist = (g1x[p] - g2x[q])^2 + (g1y[p] - g2y[q])^2

	variable i
	for(i=0; i<npts2; i+=1)
		matrixop/o one_dist = col(pdist, i)
		 wavestats/m=1/q one_dist
		 dx[i] = g2x[i] - g1x[V_minloc]
		 dy[i] = g2y[i] - g1y[V_minloc]
	
	endfor
	
	Killwaves pdist, one_dist
end

// generates a 2D lattice with basis vectors a =(ax, ay) and b = (bx, by), starting from
// the origin point (x0, y0).  Generates points from (namin to namax) along a and 
// (nbmin to nbmax) along b.  namin and nbmin can be negative.  Output results are in
// waves gx and gy.
function MakeGrid(ax, ay, bx, by, x0, y0, namin, namax, nbmin, nbmax)
	variable ax, ay, bx, by, x0, y0, namin, namax, nbmin, nbmax
	
	variable npts = ((namax - namin) + 1) * ((nbmax - nbmin)+1)
	make/o/n=(npts) gx, gy
	
	variable i, j, k=0
	for(i=namin; i<=namax; i+=1)
		for(j=nbmin; j<=nbmax; j+=1)
			gx[k] = x0 + ax*i + bx*j
			gy[k] = y0 + ay*i + by*j
			k+=1
		endfor
	endfor
	
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
	