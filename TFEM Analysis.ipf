#pragma rtGlobals=1		// Use modern global access method.

// functions for time FEM analysis
// started 05/29/13 pmv
// following "relaxation v2.ipf" by Li He

// Calculates the time auto-correlation function g2(t) from a single, 1-D time series
// of data.  Takes time scaling from the input wave scaling.  Same algorithm as Li He's
// version, but somewhat faster implementation using wave operations instead of loops.
function OneG2t(tseries)
	wave tseries
	
	// variable ptimer = startMSTimer
	
	variable tp = numpnts(tseries)
	
	make/o/n=(tp) t1, g2t
	make/o/n=(2*tp) t2
	t2[0, tp-1] = tseries[p]
	t2[tp, ] = 0

	variable i
	for(i=0; i<tp; i+=1)
		Multithread t1 = tseries[p]*t2[p+i]
		g2t[i] = sum(t1)
		g2t[i] /= ( (tp-i)*mean(t2, 0, tp-i-1)^2)
	endfor
	
	 //variable av = mean(tseries)
	 //g2t /= ( (tp-p)*av^2)
		
	SetScale/P x 0, deltax(tseries), "", g2t
	Killwaves t1, t2
	
	// variable ptime = StopMSTimer(ptimer)
	// printf "Execution time is %g ms.\r", ptime/1000
end



// Calculates g2t for every pixel (kx, ky) in a stack of diffraction patterns.  dat
// is a 3D wave containing DPs in x and y and time in z.  Takes time and k scaling
// from the input wave.  Points at k radius greater than max_kr are not calculated
// to save time.  Warning: can be very slow for large data sets.
function VecKG2t(dat, max_kr)
	wave dat
	variable max_kr
	
	max_kr = max_kr^2
	
	variable nkx, nky, nkz, kx, ky, i, j
	nkx = DimSize(dat, 0)
	nky = DimSize(dat, 1)
	nkz = DimSize(dat, 2)
	
	Make/O/n=(nkx, nky, nkz) Veck_g2t
	Veck_g2t = 0
	
	for(i=0; i<nkx; i+=1)
		kx = DimOffset(dat, 0) + DimDelta(dat, 0)*i
		if(!mod(i, 10))
			printf "row = %d\r", i
		endif
		for(j=0; j<nky; j+=1)
			ky = DimOffSet(dat, 1) + DimDelta(dat, 1)*j
			if( (kx^2 + ky^2) <= max_kr )
				ImageTransform/O/beam={(i), (j)} getbeam dat
				wave tseries = $"W_Beam"
				OneG2t(tseries)
				wave g2t = $"g2t"
				ImageTransform/beam={(i), (j)}/D=g2t setbeam Veck_g2t
			endif
		endfor
	endfor
	
	SetScale/P x DimOffset(dat, 0), DimDelta(dat, 0), "", Veck_g2t
	SetScale/P y DimOffset(dat, 1), DimDelta(dat, 1), "",Veck_g2t
	SetScale/P z 0, DimDelta(dat, 2), "", Veck_g2t
	
	Killwaves tseries, g2t
end


// Calculates g2t for every scalar k by an annular average at each time.  vecg2t is the
// output from VecKG2t above.  xc and yc are pattern center positions (in scaled units),
// strip_width is annular integration width in pixels, xbmin, xbmax, ybmin, ybmax are
// beamstop block positions, also in pixels.
function ScalarKG2t(vecg2t, xc, yc, strip_width, xbmin, xbmax, ybmin, ybmax)
	wave vecg2t
	variable xc, yc, strip_width, xbmin, xbmax, ybmin, ybmax
	
	ImageTransform/PTYP=0/P=0 getplane vecg2t
	wave onet = $"M_ImagePlane"
	
	variable rad, nt, i
	rad = AnnularAverage(onet, xc, yc, strip_width)
	nt = DimSize(vecg2t, 2)

	make/O/N=(nt, rad) Scalark_g2t
	
	for(i=0; i<nt; i+=1)
		ImageTransform/PTYP=0/P=(i) getplane vecg2t
		wave onet = $"M_ImagePlane"
		onet[xbmin, xbmax][ybmin,ybmax] = nan
		AnnularAverage(onet, xc, yc, strip_width)
		wave aav = $"annular_av"
		Scalark_g2t[i][] = aav[q]
	endfor
	
	SetScale/P x 0, DimDelta(vecg2t, 2), "", Scalark_g2t
	SetScale/P y 0, deltax(aav), "", Scalark_g2t
	
	Killwaves M_ImagePlane, annular_av, npix
end		


// Fitting function for a single stretched exponential	 
Function relaxationfunc(w,t) : FitFunc
	Wave w
	Variable t

	//CurveFitDialog/ These comments were created by the Curve Fitting dialog. Altering them will
	//CurveFitDialog/ make the function less convenient to work with in the Curve Fitting dialog.
	//CurveFitDialog/ Equation:
	//CurveFitDialog/ f(t) =A+B*exp(-2*(t/tau)^beta)
	//CurveFitDialog/ End of Equation
	//CurveFitDialog/ Independent Variables 1
	//CurveFitDialog/ t
	//CurveFitDialog/ Coefficients 4
	//CurveFitDialog/ w[0] = A
	//CurveFitDialog/ w[1] = B
	//CurveFitDialog/ w[2] = tau
	//CurveFitDialog/ w[3] = beta
	
	return w[0]+w[1]*exp(-2*(t/w[2])^w[3])
End
	
	
// fitting function for the sum of two stretched exponentials with different
// relaxation times and exponents.	 
Function relax_two_tau(w,t) : FitFunc
	Wave w
	Variable t

	//CurveFitDialog/ These comments were created by the Curve Fitting dialog. Altering them will
	//CurveFitDialog/ make the function less convenient to work with in the Curve Fitting dialog.
	//CurveFitDialog/ Equation:
	//CurveFitDialog/ f(t) = A+B1*exp(-2*(t/tau1)^beta1) + B2*exp(-2*(t/tau2)^beta2)
	//CurveFitDialog/ End of Equation
	//CurveFitDialog/ Independent Variables 1
	//CurveFitDialog/ t
	//CurveFitDialog/ Coefficients 7
	//CurveFitDialog/ w[0] = A
	//CurveFitDialog/ w[1] = B1
	//CurveFitDialog/ w[2] = tau1
	//CurveFitDialog/ w[3] = beta1
	//CurveFitDialog/ w[4] = B2
	//CurveFitDialog/ w[5] = tau2
	//CurveFitDialog/ w[6] = beta2

	return w[0]+w[1]*exp(-2*(t/w[2])^w[3]) + w[4]*exp(-2*(t/w[5])^w[6])
End
	 