#pragma rtGlobals=1		// Use modern global access method.

// Functions for annular / azimuthal processing of 2D images.
//
//  AnnularAverage
//  AnnularIntegrate
//  AnnularAverageCorners
//  Angle
//  ExcludedSectorAnnularAverage
//
// started 10-09-03 pmv
// AnnularAverageCorners added, AnnularAverage set to ignore pixels that are nan 08-06-04 pmv

// Does a rotational average about the specified point, xx and yy with a bin width of strip_width.
// xx and yy are in real units, taking into account the wave scaling.  strip_width is in pixels.
// average is placed the wave annular_av, and the function returns the number of pixels in 
// the annular average.
function AnnularAverage(dat, xx, yy, strip_width)
	wave dat
	variable xx, yy, strip_width
	
	// translate center point from wave scaling units to pixels
	variable cp = round( ((xx - DimOffset(dat, 0))/DimDelta(dat, 0)) )
	variable cq = round( ((yy - DimOffSet(dat, 1))/DimDelta(dat, 1)) )
	
	// find the radius of the annular average accounting for strip_width
	variable sx = DimSize(dat, 0)
	variable sy = DimSize(dat, 1)
	variable rad = max( max(cp, abs(cp-sx)), max(cq, abs(cq-sy)) )
	rad = floor( rad / strip_width)
	
	// make the output and pixel-counting waves
	Make/O/N=(rad) annular_av, npix
	SetScale/P x 0, DimDelta(dat, 0)*strip_width, "", annular_av
	annular_av = 0
	npix = 0
	
	// perform the annular average
	variable i, j, rpix
	for(i=0; i<sx; i+=1)
		for(j=0; j<sy; j+=1)
			if(!NumType(dat[i][j]))
				rpix = round( sqrt( (i-cp)^2 + (j-cq)^2) / strip_width )
				if(rpix < rad)
					annular_av[rpix] += dat[i][j]
					npix[rpix] += 1
				endif
			endif
		endfor
	endfor
	
	// normalize by the number of pixels in each bin
	annular_av /= npix

	// clean up temporary waves
	//Killwaves/Z npix

	// return the radius of the annular average
	return rad	
end


// Does a rotational integral about the specified point, xx and yy with a bin width of strip_width.
// xx and yy are in real units, taking into account the wave scaling.  strip_width is in pixels.
// average is placed the wave annular_av, and the function returns the number of pixels in 
// the annular average.
function AnnularIntegral(dat, xx, yy, strip_width)
	wave dat
	variable xx, yy, strip_width
	
	// translate center point from wave scaling units to pixels
	variable cp = round( ((xx - DimOffset(dat, 0))/DimDelta(dat, 0)) )
	variable cq = round( ((yy - DimOffSet(dat, 1))/DimDelta(dat, 1)) )
	
	// find the radius of the annular average accounting for strip_width
	variable sx = DimSize(dat, 0)
	variable sy = DimSize(dat, 1)
	variable rad = min( min(cp, abs(cp-sx)), min(cq, abs(cq-sy)) )
	rad = floor( rad / strip_width)
	
	// make the output and pixel-counting waves
	Make/O/N=(rad) annular_int
	SetScale/P x 0, DimDelta(dat, 0)*strip_width, "", annular_int
	annular_int= 0
	
	// perform the annular average
	variable i, j, rpix
	for(i=0; i<sx; i+=1)
		for(j=0; j<sy; j+=1)
			if(!NumType(dat[i][j]))
				rpix = round( sqrt( (i-cp)^2 + (j-cq)^2) / strip_width )
				if(rpix < rad)
					annular_int[rpix] += dat[i][j]
				endif
			endif
		endfor
	endfor
	
	// return the radius of the annular average
	return rad	
	
end


function Angle(xx, yy, xc, yc)
	variable xx, yy, xc, yc
	
	variable xd = (xx - xc)
	variable yd = (yy - yc)
	
	variable angle = atan(yd / xd) 
	
	if(yd > 0)
		if(xd >= 0)
			// do nothing
		else
			angle += Pi
		endif
	else
		if(xd >= 0)
			angle += 2*Pi
		else
			angle += Pi
		endif
	endif
	
	return angle
end


function ExcludedSectorAnnularAverage(dat, xc, yc, tmin, tmax)
	wave dat
	variable xc, yc, tmin, tmax
	
	// translate center point from wave scaling units to pixels
	variable cp = round( ((xc - DimOffset(dat, 0))/DimDelta(dat, 0)) )
	variable cq = round( ((yc - DimOffSet(dat, 1))/DimDelta(dat, 1)) )
	
	// find the radius of the annular average accounting for strip_width
	variable sx = DimSize(dat, 0)
	variable sy = DimSize(dat, 1)
	variable rad = min( min(cp, abs(cp-sx)), min(cq, abs(cq-sy)) )
	//rad = floor( rad / strip_width)
	
	// make the output and pixel-counting waves
	Make/O/N=(rad) annular_av, npix
	SetScale/P x 0, DimDelta(dat, 0), "", annular_av
	annular_av= 0
	npix = 0
	
	// perform the annular average
	variable i, j, rpix, t
	for(i=0; i<sx; i+=1)
		for(j=0; j<sy; j+=1)
			rpix = round( sqrt( (i-cp)^2 + (j-cq)^2) )
			if(rpix < rad)
				t = angle( i, j, cp, cq )
				if(t < tmax & t > tmin)
					annular_av[rpix] += dat[i][j]
					npix[rpix] += 1
				endif
			endif
		endfor
	endfor
	
	annular_av /= npix
	Killwaves npix
	
	// return the radius of the annular average
	return rad	

end	


function AnnularAverageCorners(dat, xx, yy, strip_width)
	wave dat
	variable xx, yy, strip_width
	
	// translate center point from wave scaling units to pixels
	variable cp = round( ((xx - DimOffset(dat, 0))/DimDelta(dat, 0)) )
	variable cq = round( ((yy - DimOffSet(dat, 1))/DimDelta(dat, 1)) )
	
	// find the radius of the annular average accounting for strip_width
	variable sx = DimSize(dat, 0)
	variable sy = DimSize(dat, 1)
	variable c1, c2, c3, c4
	c1 = cp^2+cq^2
	c2 = (cp-sx)^2 + cq^2
	c3 = cp^2 + (cq-sy)^2
	c4 = (cp-sx)^2 + (cq-sy)^2
	
	variable rad = sqrt(max( max( c1, c2 ), max(c3, c4 )))
	rad = floor( rad / strip_width)
	
	//printf "The point closest to (%g, %g) is (%g, %g).\r", xx, yy, cp*DimDelta(dat, 0) + DimOffSet(dat, 0), cq*DimDelta(dat, 1) + DimOffSet(dat, 1)
	//printf "Annular average runs out to %g pixels\r", rad

	// make the output and pixel-counting waves
	Make/O/N=(rad) annular_av, npix
	SetScale/P x 0, DimDelta(dat, 0)*strip_width, "", annular_av
	annular_av = 0
	npix = 0
	
	// perform the annular average
	variable i, j, rpix
	for(i=0; i<sx; i+=1)
		for(j=0; j<sy; j+=1)
			if(!NumType(dat[i][j]))  // dat[i][j] is not nan
				rpix = round( sqrt( (i-cp)^2 + (j-cq)^2) / strip_width )
				if(rpix < rad)
					annular_av[rpix] += dat[i][j]
					npix[rpix] += 1
				endif
			endif
		endfor
	endfor
	
	// normalize by the number of pixels in each bin
	annular_av /= npix

	// clean up temporary waves
	Killwaves/Z npix

	// return the radius of the annular average
	return rad	
	
end
