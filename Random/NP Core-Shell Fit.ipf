#pragma rtGlobals=1		// Use modern global access method.
#include <All IP Procedures>
#include <Image Saver>


// creates the fit function for the NP sphere function on a support.
Function NPsphere(w,x,y) : FitFunc
	Wave w
	Variable x
	Variable y

	//CurveFitDialog/ These comments were created by the Curve Fitting dialog. Altering them will
	//CurveFitDialog/ make the function less convenient to work with in the Curve Fitting dialog.
	//CurveFitDialog/ Equation:
	//CurveFitDialog/ f(x,y) = (P1*2*sqrt((r^2)-((x-xc)^2)-((y-yc)^2)))+P2
	//CurveFitDialog/ End of Equation
	//CurveFitDialog/ Independent Variables 2
	//CurveFitDialog/ x
	//CurveFitDialog/ y
	//CurveFitDialog/ Coefficients 5
	//CurveFitDialog/ w[0] = P1
	//CurveFitDialog/ w[1] = P2
	//CurveFitDialog/ w[2] = r
	//CurveFitDialog/ w[3] = xc
	//CurveFitDialog/ w[4] = yc

	return (w[0]*2*sqrt((w[2]^2)-((x-w[3])^2)-((y-w[4])^2)))+w[1]
End



// creates the mask for the NP core
Function npCoremask(image, t)
	Wave image
	variable t
	
	duplicate/o image mask
	
	variable numx=DimSize(image,0)
	variable numy=DimSize(image,1)  
	
	variable p
	
	variable i, j
	i=0
	for(i=0; i<numx; i+=1)
		j=0
		for(j=0; j<numy; j+=1)
			p= image[i][j]
			if( p >= t )
				mask[i][j] = 0
			else
				mask[i][j] = 1
			endif
		endfor
	endfor
End


function npCoreShellFit(image, t, P1_i, P2_i, r, xc, yc)
	
	wave image
	variable t, P1_i, P2_i, r, xc, yc
	
	npCoremask(image, t)
	
	wave mask = $"mask"
	
	variable V_FitError
		
	// Make initial guesses for fit.
	Make/D/N=5/O W_coef
	W_coef[0] = {P1_i, P2_i, r, xc, yc}
		
	V_FitError = 0
	FuncFitMD/H="00111"/NTHR=0 NPsphere W_coef  image /M=mask /D /R
	//FuncFitMD/NTHR=0 NPsphere W_coef  image /M=mask /D /R
	
	if (V_FitError !=0)
		print "error in fit - V_FitError =", V_FitError
	endif
		
end
