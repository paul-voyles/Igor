#pragma rtGlobals=1		// Use modern global access method.

function wavelen(kV)
	variable kV
	
	return 1.2398 / sqrt(kV * (2*511 + kV))
end

// returns Scherzer aperture half-angle in mrad.
// takes accelerating voltage in kilovolts and
// spherical aberration coefficient in mm.
function ScherzerAperture(kV, Cs)
	variable kV, Cs

	Cs *= 10^6	// into nm
	variable wl  = wavelen(kV) 	// in nm

	return 1000*(6*wl / Cs)^0.25

end

// returns Scherzer defocus in nm.
// takes accelerating voltage in kilovolts and
// spherical aberration coefficient in mm.
function ScherzerDefocus(kV, Cs)
	variable kV, Cs

	Cs *= 10^6	// into nm
	variable wl  = wavelen(kV) 	// in nm

	return sqrt( 1.5*Cs*wl )

end

// kV in kilovolts
// df in nm
// Cs in mm
// ap in mrad
// rmin, rmax in nm
function ProbeProfile(kV, df, Cs, ap, rmin, rmax, npts)
	variable kV, df, Cs, ap, rmin, rmax, npts
	
	variable aber_pnts = 100

	variable wl  = wavelen(kV)				 	// in nm
	variable kmax = 0.001*ap / wl				// in nm^-1
	Cs *= 10^6	// in nm
	
	Make/O/N=(npts) profile
	SetScale/i x rmin, rmax, "", profile
	
	Make/O/C/N=(aber_pnts) aber, int
	SetScale/i x 0, kmax, "", aber, int
	
	aber = Pi*wl*x^2*(0.5*Cs*wl^2*x^2 - df)
	aber = x*exp( -cmplx(0, 1)*aber[p])
	
	variable i, r
	variable/C one_p
	for(i=0; i<npts; i+=1)
		r = pnt2x(profile, i)
		int = cmplx(bessj(0, 2*Pi*x*r, 0), 0)
		int *= aber[p]
		one_p = area(int, 0, kmax)
		profile[i] = magsqr(one_p)
	endfor

	profile[1, npts-1] /= profile[0]
	profile[0] = 1
	
	KillWaves aber, int
end
	
	
function ProbeFWHM(kV, df, Cs, ap)
	variable kV, df, Cs, ap
	
	variable rmax = 0.5	// nm
	variable rpnts = 100
	
	ProbeProfile(kV, df, Cs, ap, 0, rmax, rpnts)
	wave profile = $"profile"
	
	Wavestats/Q profile
	profile /= V_max
	FindLevels/Q profile, 0.5
	
	if(V_flag == 0)
		printf "The entire profile is equal to 0.5!  Cannot continue.\r"
		return nan
	elseif (V_flag == 2)
		printf "The probe profile does not cross 0.5.  Cannot continue.\r"
		return nan
	endif
	
	if(V_LevelsFound != 1)
		printf "The probe profile crosses 0.5 more than once.  FWHM is not a good measure.\r"
		return nan
	endif
	
	wave fl = $"W_FindLevels"
	variable FWHM = 2*fl[0]
	KillWaves fl
	
	return FWHM
	
	
end
	