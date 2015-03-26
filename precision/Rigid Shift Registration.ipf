#pragma rtGlobals=1		// Use modern global access method.

// Functions for rigid shift registration by cross correlation
//
// v1 pmv 08-02-12
//
// Use PairwiseRigidShift
//
// 01-28-14 removed Imagetransform swap in TwoImageFindShift. no longer needed due to change in behavior of MatrixOp
//               added Hanning window before MatrixOp correlation to reduce aliasing artifacts

// f1 is reference, f2 is shift.  Reported shift works in imagetransform/trns
// fit_win is NxN box to fit to cross correlation image to find shift with
// sub-pixel accuracy.  max_shift is maximum allowed shift without an 
// error - must be less than one unit cell.  (sx, sy) return the x and y shift
// respectively, all in pixels.
function TwoImageFindShift(f1, f2, fit_win, max_shift, sx, sy, windowYN)
	wave f1, f2
	variable fit_win, max_shift
	variable &sx, &sy
	variable windowYN
	
	if(windowYN)
		ImageWindow/O Bartlett f1
		ImageWindow/O Bartlett f2
	endif
	
	Matrixop/O cr = correlate(f1, f2, 0)
	// Imagetransform swap cr

	wavestats/q/m=1 cr
	
	variable dist
	dist = sqrt( (V_maxrowloc - DimSize(f1, 0)/2)^2 +  (V_maxcolloc - DimSize(f1, 1)/2)^2 )
	if(dist > max_shift)
		printf "Rough shift estimate of %g is larger than maximum shift of %g.  Exiting.\r", dist, max_shift
		sx = 0
		sy = 0
		return 0
	endif	
	
	CurveFit/Q/NTHR=0/N/W=2 Gauss2D  cr[ (V_maxrowloc - fit_win),(V_maxrowloc+fit_win)][(V_maxcolloc - fit_win),(V_maxcolloc+fit_win)] /D 
	wave W_Coef = $"W_Coef"
	
	sx = -(W_Coef[2] - DimSize(f1,0)/2)
	sy = -(W_Coef[4] - DimSize(f1,1)/2)
	
	
	killwaves cr, W_coef, W_sigma, fit_cr
	
end
	
function ShiftTest(f1, f2, fit_win, max_shift)
	wave f1, f2
	variable fit_win, max_shift

	variable sx, sy
	
	TwoImageFindShift(f1, f2, fit_win, max_shift, sx, sy, 0)
	
	printf "Shift is (%g, %g) pixels.\r", sx, sy
	
end

function GlobalRigidShift(st, ref, fit_win, max_shift)
	wave st
	variable ref, fit_win, max_shift
		
	duplicate/o st st_algn
	variable nim = DimSize(st, 2)
	variable i, sx, sy
	
	Imagetransform/P=(ref) getplane st
	wave tmp = $"M_ImagePlane"
	Duplicate/O tmp ref_im
	
	for(i=0; i<nim; i+=1)	
		if(i==ref)
			st_algn[][][i] = ref_im[p][q]
		else
		
			Imagetransform/P=(i) getplane st
			wave tmp = $"M_ImagePlane"
			duplicate/O tmp shift_im
			TwoImageFindShift(ref_im, shift_im, fit_win, max_shift, sx, sy, 0)
			ImageInterpolate/TRNS={scaleshift, sx, 1, sy, 1} resample shift_im
			wave int_im = $"M_InterpolatedImage"
			st_algn[][][i] = int_im[p][q]
		
		endif
	endfor
	
	Killwaves ref_im, shift_im, M_ImagePlane, M_InterpolatedImage
	
end

// st is the 3D image stack.  Ref is the reference image out of the stack that all the
// other images are aligned to.  fit_win and max_shift are defined in TwoImageFindShift
// currently there's a bug - this only works for ref = 0
function PairwiseRigidShift(st, ref, fit_win, max_shift, windowYN)
	wave st
	variable ref, fit_win, max_shift, windowYN
	
	duplicate/O st, st_algn
	
	variable nim = DimSize(st, 2)
	make/O/N=(nim) sx, sy
	variable sxt, syt, i
	
	for(i=1; i<nim; i+=1)
		ImageTransform/p=(i-1) getplane st
		wave tmp = $"M_ImagePlane"
		duplicate/o tmp f1
		ImageTransform/p=(i) getplane st 
		wave tmp = $"M_ImagePlane"
		duplicate/O tmp f2
		TwoImageFindShift(f1, f2, fit_win, max_shift, sxt, syt, windowYN)
		sx[i] = sxt
		sy[i] = syt
	endfor
	
	for(i=0; i<nim; i+=1)
		if(i == ref)
			st_algn[][][i] = st[p][q][i]
		else
			Imagetransform/p=(i) getplane st
			wave f3 = $"M_ImagePlane"
			
			sxt = sum(sx, ref, i)
			syt = sum(sy, ref, i)
			
			ImageInterpolate/TRNS={scaleshift, sxt, 1, syt, 1} resample f3
			wave shift_im = $"M_InterpolatedImage"
			st_algn[][][i] = shift_im[p][q]
		endif
	endfor
	
	Killwaves f1, f2, M_ImagePlane, M_InterpolatedImage	
	
end