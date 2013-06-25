#pragma rtGlobals=1		// Use modern global access method.

// Procedures for aligning low-loss spectra using the ZLP and measuring the plasmon
// energy for a series of spectra.
//
// begun 11/14/12, pmv


function PeakMaxFind(spec, se_min, se_max,peak_frac)
	wave spec
	variable se_min, se_max, peak_frac
	
	variable e_min, e_max
	
	
	Wavestats/R=(se_min, se_max)/Q spec
	FindLevel/Q/R=(se_min, V_maxloc) spec, peak_frac*V_max
	if(V_flag)
		printf "Cannot find left edge of ZLP in wave %s.  Exiting.\r", NameofWave(spec)
		return 0
	endif
	e_min =V_LevelX
	
	FindLevel/Q/R=(V_maxloc,se_max) spec, peak_frac*V_max
	if(V_flag)
		printf "Cannot find right edge of ZLP in wave %s.  Exiting.\r", NameofWave(spec)
		return 0
	endif
	e_max =V_LevelX
	
	// printf "Fitting from %g to %g eV.\r", e_min, e_max
	
	CurveFit/NTHR=0/Q/N=1/W=2 gauss spec(e_min, e_max)
	
	wave wc = $"W_Coef"
	variable mid = wc[2]
	Killwaves W_coef, W_sigma
	
	return mid
end

function ZLPPeakFind(spec)
	wave spec
	
	variable mid = PeakMaxFind(spec, leftx(spec), rightx(spec), 0.8)
	return mid
end
	
function ShiftZLP(spec)
	wave spec
	
	variable mid = ZLPPeakFind(spec)
	
	SetScale/P x (leftx(spec)-mid), deltax(spec), spec
	
end

function PlasmonPeakFind(spec, e_min, e_max)
	wave spec
	variable e_min, e_max
	
	variable ep
	ep = PeakMaxFind(spec, e_min, e_max, 0.9)
	
	return ep
	
end


function FolderShiftZLP()

	variable i
	string wl = WaveList("*", ";", "")
	
	for(i=0; i<ItemsInList(wl, ";"); i+=1)
		wave spec = $(StringFromList(i, wl, ";"))
		ShiftZLP(spec)
	endfor

end

function FolderPlasmonPeakFind(e_min, e_max)
	variable e_min, e_max
	
	string wl = WaveList("*", ";", "")
	variable nspec = ItemsInList(wl, ";")
	
	make/o/t/n=(nspec) spec_names
	make/O/N=(nspec) ep
	
	variable i
	for(i=0; i<nspec; i+=1)
		wave spec = $(StringFromList(i, wl, ";"))
		spec_names[i] = StringFromList(i, wl, ";'")
		ep[i] = PlasmonPeakFind(spec, e_min, e_max)
	endfor
	
end
		

function CCDBoundaryFix(spec, win)
	wave spec
	variable win
	
	CurveFit/Q/W=2/NTHR=0 line  spec[(1023-win),1023]
	wave W_Coef =$"W_Coef"
	variable low_e_cts = W_Coef[0] + W_coef[1]*pnt2x(spec, 1023.5)
	
	CurveFit/Q/W=2/NTHR=0 line  spec[1024,1024+win]
	variable high_e_cts = W_Coef[0] + W_coef[1]*pnt2x(spec,1023.5)
	
	spec[1024,] -= (high_e_cts - low_e_cts)
	Killwaves W_Coef, W_Sigma
	
end
	
function CCDBoundaryFixFolder(win)
	variable win

	variable i
	string wl = WaveList("*", ";", "")
	
	for(i=0; i<ItemsInList(wl, ";"); i+=1)
		wave spec = $(StringFromList(i, wl, ";"))
		CCDBoundaryFix(spec, win)
	endfor

end

function FolderDoubleGaussFit(coef, e_min, e_max)
	wave coef
	variable e_min, e_max
	
	string wl = WaveList("*", ";", "")
	variable nspec = ItemsInList(wl, ";")-1
	
	make/o/t/n=(nspec) spec_names
	make/O/N=(nspec) ep_g1, ep_g2
	make/O/N=(nspec, 8) all_coef
	
	variable i
	for(i=0; i<nspec; i+=1)
		duplicate/O coef W_Coef
		if(cmpstr(StringFromList(i,wl, ";"), NameofWave(coef)))
			wave spec = $(StringFromList(i, wl, ";"))
			FuncFit/NTHR=0/Q/W=2 DoubleGauss W_coef  spec(e_min, e_max)	
			spec_names[i] = StringFromList(i, wl, ";'")
			all_coef[i][] = W_coef[q]
			ep_g1[i] = W_coef[2]
			ep_g2[i] = W_coef[5]
		endif
	endfor
	
	Killwaves W_coef, W_sigma
	
end

function DoubleGaussContributions(w, e_min, e_max)
	wave w
	variable e_min, e_max
	
	make/o/n=300 dg_line, dg_g1, dg_g2, dg_tot
	setscale/I x e_min, e_max, dg_line, dg_g1, dg_g2
	
	dg_line = w[0] + w[7]*x
	dg_g1 = dg_line + w[1]*exp( -(x-w[2])^2 / w[3]^2 )
	dg_g2 = dg_line + w[4]*exp( -(x-w[5])^2 / w[6]^2 )
	dg_tot = DoubleGauss(w, x)
	
end
	
Function DoubleGauss(w,x) : FitFunc
	Wave w
	Variable x

	//CurveFitDialog/ These comments were created by the Curve Fitting dialog. Altering them will
	//CurveFitDialog/ make the function less convenient to work with in the Curve Fitting dialog.
	//CurveFitDialog/ Equation:
	//CurveFitDialog/ f(x) = y0 + b*x + A1*exp( -(x-x01)^2 / w1^2 ) + a2*exp( -(x-x02)^2 / w2^2 )
	//CurveFitDialog/ End of Equation
	//CurveFitDialog/ Independent Variables 1
	//CurveFitDialog/ x
	//CurveFitDialog/ Coefficients 8
	//CurveFitDialog/ w[0] = y0
	//CurveFitDialog/ w[1] = A1
	//CurveFitDialog/ w[2] = x01
	//CurveFitDialog/ w[3] = w1
	//CurveFitDialog/ w[4] = a2
	//CurveFitDialog/ w[5] = x02
	//CurveFitDialog/ w[6] = w2
	//CurveFitDialog/ w[7] = b

	return w[0] + w[7]*x + w[1]*exp( -(x-w[2])^2 / w[3]^2 ) + w[4]*exp( -(x-w[5])^2 / w[6]^2 )
End
