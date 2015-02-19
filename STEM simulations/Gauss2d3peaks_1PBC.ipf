#pragma rtGlobals=1		// Use modern global access method.

// takes an image im and crops it with boundary xi yi xf yf.  
// xi yi xf yf are all within the final cropped image.
function cropimage(im, xi, yi, xf, yf)

wave im
variable xi, yi, xf, yf

variable scalex = DimDelta(im, 0)
variable scaley = DimDelta(im, 1)

variable xi_pixel = (xi/scalex)
variable xf_pixel = (xf/scalex) + 1
variable yi_pixel = (yi/scaley)
variable yf_pixel = (yf/scaley) + 1

DeletePoints xf_pixel,10000, im
DeletePoints 0,xi_pixel, im
DeletePoints/M=1 yf_pixel,10000, im
DeletePoints/M=1 0,yi_pixel, im

end



//Creates the Gauss2d*3 function.
Function Gauss2d3peaks(w, x, y) : FitFunc
	Wave w
	Variable x
	Variable y
	
	//Revised form Gauss2D*2
	//CurveFitDialog/ These comments were created by the Curve Fitting dialog. Altering them will
	//CurveFitDialog/ make the function less convenient to work with in the Curve Fitting dialog.
	//CurveFitDialog/ Equation:
	//CurveFitDialog/ f(x,y) = z0+(A*exp((-1/(2*(1-(cor1)^2)))*((((x-x1)/xw1)^2)+(((y-y1)/yw1)^2)-((2*cor1*(x-x1)*(y-y1))/(xw1*yw1)))))+ (B*exp((-1/(2*(1-(cor2)^2)))*((((x-x2)/xw2)^2)+(((y-y2)/yw2)^2)-((2*cor2*(x-x2)*(y-y2))/(xw2*yw2)))))
	//CurveFitDialog/ End of Equation
	//CurveFitDialog/ Independent Variables 2
	//CurveFitDialog/ x
	//CurveFitDialog/ y
	//CurveFitDialog/ Coefficients 19
	//CurveFitDialog/ w[0] = z0
	//CurveFitDialog/ w[1] = A
	//CurveFitDialog/ w[2] = B
	//CurveFitDialog/ w[3] = cor1
	//CurveFitDialog/ w[4] = cor2
	//CurveFitDialog/ w[5] = x1
	//CurveFitDialog/ w[6] = x2
	//CurveFitDialog/ w[7] = y1
	//CurveFitDialog/ w[8] = y2
	//CurveFitDialog/ w[9] = xw1
	//CurveFitDialog/ w[10] = xw2
	//CurveFitDialog/ w[11] = yw1
	//CurveFitDialog/ w[12] = yw2
	//CurveFitDialog/ w[13] = C
	//CurveFitDialog/ w[14] = cor3
	//CurveFitDialog/ w[15] = x3
	//CurveFitDialog/ w[16] = y3
	//CurveFitDialog/ w[17] = xw3
	//CurveFitDialog/ w[18] = yw3

	return w[0]+(w[1]*exp((-1/(2*(1-(w[3])^2)))*((((x-w[5])/w[9])^2)+(((y-w[7])/w[11])^2)-((2*w[3]*(x-w[5])*(y-w[7]))/(w[9]*w[11])))))+ (w[2]*exp((-1/(2*(1-(w[4])^2)))*((((x-w[6])/w[10])^2)+(((y-w[8])/w[12])^2)-((2*w[4]*(x-w[6])*(y-w[8]))/(w[10]*w[12])))))+(w[13]*exp((-1/(2*(1-(w[14])^2)))*((((x-w[15])/w[17])^2)+(((y-w[16])/w[18])^2)-((2*w[14]*(x-w[15])*(y-w[16]))/(w[17]*w[18])))))
	
End


function Gauss2D3PeakFit(image, x_loc, y_loc, size_x, size_y)

	wave image, x_loc, y_loc
	variable size_x, size_y
	variable num_peaks = DimSize(x_loc, 0)
	
	variable scalex = DimDelta(image, 0)
	variable scaley = DimDelta(image, 1)
	
	Make/D/N = (num_peaks/3)/O x_start
	Make/D/N = (num_peaks/3)/O x_finish
	Make/D/N = (num_peaks/3)/O y_start
	Make/D/N = (num_peaks/3)/O y_finish
	
	variable x_size = (round((size_x/2)/scalex))*scalex
	variable y_size = (round((size_y/2)/scaley))*scaley
	
	Make/D/N=(num_peaks/3)/O xc, yc
	
	variable i = 0
	
	for(i=0; i<num_peaks; i+=3)
		xc[i/3] = (round((x_loc[i] + x_loc[i+1]+x_loc[i+2])/(3*scalex)))*scalex
		yc[i/3] = (round((y_loc[i] + y_loc[i+1]+y_loc[i+2])/(3*scaley)))*scaley

		x_start[i/3] = xc[i/3] - x_size
		x_finish[i/3] = xc[i/3] + x_size - scalex
		y_start[i/3] = yc[i/3] - y_size
		y_finish[i/3] = yc[i/3] + y_size - scaley
	endfor	
	
	Make/D/O/N=(num_peaks) z0, A, x0, xW, y0, yW, cor
	Make/D/O/N=(num_peaks) sigma_z0, sigma_A, sigma_x0, sigma_xW, sigma_y0, sigma_yW, sigma_cor
	
	variable V_FitError
	
	// Make initial guesses for fit.
	variable cor1s = 0.01, cor2s = 0.01, cor3s = 0.01
	variable xw1s = 0.25, xw2s = 0.25, xw3s = 0.25, yw1s = 0.25, yw2s = 0.25, yw3s = 0.25
	
	i=0
	for(i=0; i<num_peaks; i+=3)
		//crop out, save, and setscale of area to be fitted
		Duplicate/O image im
		cropimage(im, x_start[i/3], y_start[i/3], x_finish[i/3], y_finish[i/3])
		if( i == 0 )
			Make/O/N=(Dimsize(im,0), Dimsize(im,1), (num_peaks/3)) im_t
		endif
		imagetransform /INSI=im /INSX=0 /INSY=0 /P=(i/3) insertImage im_t
		setscale/P x, x_start[i/3], scalex, im
		setscale/P y, y_start[i/3], scaley, im
		
		Make/O/N=(DimSize(im,0), DimSize(im,1)) residual
		
		Duplicate/O/D im noise
		noise = sqrt(im)
		
		//make initial guesses for fitting.
		Make/D/N=19/O W_coef
		Make/D/O/N=27 M_WaveStats
		Make/D/O/N=(7) W_sigma
		wavestats/Q/W im
		variable z0s = M_WaveStats[10]
		variable As = image(x_loc[i])(y_loc[i])
		variable Bs = image(x_loc[i])(y_loc[i+1])
		variable Cs = image(x_loc[i+1])(y_loc[i])
		W_coef[0] = {z0s, As, Bs, cor1s, cor2s, x_loc[i], x_loc[i+1], y_loc[i], y_loc[i+1], xw1s, xw2s, yw1s, yw2s, Cs, cor3s, x_loc[i+2], y_loc[i+2], xw3s, yw3s}
		
		V_FitError = 0
		FuncFitMD/N/NTHR=0/W=2/Q  Gauss2d3peaks, W_coef, im/R=residual/W=noise/D
		
		if (V_FitError !=0)
			print "error in fit - slice #", i, "- V_FitError =", V_FitError
		endif
		
		z0[i] = W_coef(0)
		A[i] = W_coef(1)
		x0[i] = W_coef(5)
		xW[i] = W_coef(9)
		y0[i] = W_coef(7)
		yW[i] = W_coef(11)
		cor[i] = W_coef(3)
		z0[i+1] = W_coef(0)
		A[i+1] = W_coef(2)
		x0[i+1] = W_coef(6)
		xW[i+1] = W_coef(10)
		y0[i+1] = W_coef(8)
		yW[i+1] = W_coef(12)
		cor[i+1] = W_coef(4)
		z0[i+2] = W_coef(0)
		A[i+2] = W_coef(13)
		x0[i+2] = W_coef(15)
		xW[i+2] = W_coef(17)
		y0[i+2] = W_coef(16)
		yW[i+2] = W_coef(18)
		cor[i+2] = W_coef(14)
		sigma_z0[i] = W_sigma(0)
		sigma_A[i] = W_sigma(1)
		sigma_x0[i] = W_sigma(5)
		sigma_xW[i] = W_sigma(9)
		sigma_y0[i] = W_sigma(7)
		sigma_yW[i] = W_sigma(11)
		sigma_cor[i] = W_sigma(3)
		sigma_z0[i+1] = W_sigma(0)
		sigma_A[i+1] = W_sigma(2)
		sigma_x0[i+1] = W_sigma(6)
		sigma_xW[i+1] = W_sigma(10)
		sigma_y0[i+1] = W_sigma(8)
		sigma_yW[i+1] = W_sigma(12)
		sigma_cor[i+1] = W_sigma(4)
		sigma_z0[i+2] = W_sigma(0)
		sigma_A[i+2] = W_sigma(13)
		sigma_x0[i+2] = W_sigma(15)
		sigma_xW[i+2] = W_sigma(17)
		sigma_y0[i+2] = W_sigma(16)
		sigma_yW[i+2] = W_sigma(18)
		sigma_cor[i+2] = W_sigma(14)
		
		//save the fit images
		wave fit_im = $"fit_im"
		if( i == 0)
			Make/O/D/N=(Dimsize(fit_im,0), Dimsize(fit_im,1), (num_peaks/3)) fits_t
		endif	
		imagetransform /INSI=fit_im /INSX=0 /INSY=0 /P=(i/3) insertImage fits_t
		killwaves  fit_im
		
		//save the residual images
		if (i==0)
			Make/O/N=(DimSize(im,0), DimSize(im,1), (num_peaks/3)) residual_t
		endif	
		imagetransform /INSI=residual /INSX=0 /INSY=0 /P=(i/3) insertImage residual_t
		killwaves residual
		
		killwaves im, noise
		killwaves W_coef, M_WaveStats, W_sigma
	endfor

	killwaves x_start, x_finish, y_start, y_finish, xc, yc
end 

function FitStack(image_stack, startimage, x_loc, y_loc, size_x, size_y, cl)

	wave image_stack, x_loc, y_loc
	variable startimage, size_x, size_y, cl
	
	setscale/P x, 0, 0.1268, image_stack
	setscale/P y, 0, 0.1210, image_stack

	//chops the image stack to start analysis at an image other than the first image.
	Duplicate/O image_stack Image_stack_new
	if( startimage != 0)
		DeletePoints/M=2 0, startimage, image_stack_new
	endif
	
	variable num_images = DimSize(image_stack_new, 2)
	variable num_peaks  = DimSize(x_loc,0)
	
	Make/O/D/N=(num_peaks, num_images) z0_stack, A_stack, x0_stack, xW_stack, y0_stack, yW_stack, cor_stack 
	Make/O/D/N=(num_peaks, num_images) sig_z0_stack, sig_A_stack, sig_x0_stack, sig_xW_stack, sig_y0_stack, sig_yW_stack, sig_cor_stack

	variable scalex = DimDelta(image_stack, 0)
	variable scaley = DimDelta(image_stack, 1)
	
	variable i=0
	for(i=0; i<num_images; i+=1)
		Imagetransform/p=(i) getplane image_stack_new
		wave image = $"M_ImagePlane"
		setscale/P x, 0, scalex, image
		setscale/P y, 0, scaley, image
		
		Gauss2D3PeakFit(image, x_loc, y_loc, size_x, size_y)
		killwaves M_ImagePlane
		
		//save fit parameters for each atom column in each frame of stack
		
		wave z0 = $"z0"
		wave A = $"A"
		wave x0 = $"x0"
		wave xW = $"xW"
		wave y0 = $"y0"
		wave yW = $"yW"
		wave cor = $"cor"
		wave sigma_z0 = $"sigma_z0"
		wave sigma_A = $"sigma_A"
		wave sigma_x0 = $"sigma_x0"
		wave sigma_xW = $"sigma_xW"
		wave sigma_y0 = $"sigma_y0"
		wave sigma_yW = $"sigma_yW"
		wave sigma_cor = $"sigma_cor"
		wave residual = $"residual_t"
		wave fits = $"fits_t"
		wave im = $"im_t"
			
		z0_stack[][i] = z0[p] 
		A_stack[][i] = A[p] 
		x0_stack[][i] = x0[p]  
		xW_stack[][i] = xW[p]  
		y0_stack[][i] = y0[p]  
		yW_stack[][i] = yW[p]  
		cor_stack[][i] = cor[p] 
		sig_z0_stack[][i] = sigma_z0[p] 
		sig_A_stack[][i] = sigma_A[p] 
		sig_x0_stack[][i] = sigma_x0[p]  
		sig_xW_stack[][i] = sigma_xW[p]  
		sig_y0_stack[][i] = sigma_y0[p]  
		sig_yW_stack[][i] = sigma_yW[p]  
		sig_cor_stack[][i] = sigma_cor[p]  

		//save the fits, residuals, and original image fit areas.
		if (i == 0)
			variable m=0
			for(m=0; m<num_peaks/3; m+=1)
				string imT
				sprintf imT, "image_t_%d", m
				Make/O/N=(DimSize(im_t, 0), DimSize(im_t, 1), num_images) $imT
			endfor
		endif
		m=0
		for(m=0; m<num_peaks/3; m+=1)
			sprintf imT, "image_t_%d", m
			Imagetransform/p=(m) getplane im_t
			wave im = $"M_ImagePlane"
			wave ImageT = $imT
			imagetransform /INSI=im /INSX=0 /INSY=0 /P=(i) insertImage ImageT
			killwaves M_ImagePlane
		endfor
		
		if (i == 0)
			m=0
			for(m=0; m<num_peaks/3; m+=1)
				string fitT
				sprintf fitT, "fit_t_%d", m
				Make/O/D/N=(DimSize(fits_t, 0), DimSize(fits_t, 1), num_images) $fitT
			endfor
		endif
		m=0
		for(m=0; m<num_peaks/3; m+=1)
			sprintf fitT, "fit_t_%d", m
			Imagetransform/p=(m) getplane fits_t
			wave fit = $"M_ImagePlane"
			wave FitsT = $fitT
			imagetransform /INSI=fit /INSX=0 /INSY=0 /P=(i) insertImage FitsT
			killwaves M_ImagePlane
		endfor
		
		if (i == 0)
			m=0
			for(m=0; m<num_peaks/3; m+=1)
				string residT
				sprintf residT, "residual_t_%d", m
				Make/O/N=(DimSize(residual_t, 0), DimSize(residual_t, 1), num_images) $residT
			endfor
		endif
		m=0
		for(m=0; m<num_peaks/3; m+=1)
			sprintf residT, "residual_t_%d", m
			Imagetransform/p=(m) getplane residual_t
			wave resid = $"M_ImagePlane"
			wave ResidualT = $residT
			imagetransform /INSI=resid /INSX=0 /INSY=0 /P=(i) insertImage ResidualT
			killwaves M_ImagePlane
		endfor
		
		killwaves residual_t, fits_t, im_t
		killwaves z0, A, x0, xW, y0, yW, cor, sigma_z0, sigma_A, sigma_x0, sigma_xW, sigma_y0, sigma_yW, sigma_cor
	endfor
	
	GauInteg_Mn(cor_stack, yW_stack, xW_stack, A_stack, z0_stack, size_x, size_y)
	Save/G/P=home Mn_Inte as "Mn_Inte_"+num2str(cl)+".txt"
	
	Column_dis_La(x0_stack, y0_stack)			//calculate separation directly 08/25/2014, JF
	Save/G/P=home Colu_dis_0 as "Colu_dis_"+num2str(cl)+"_LaLa_0.txt"	
	Save/G/P=home Colu_dis_1 as "Colu_dis_"+num2str(cl)+"_LaLa_1.txt"	
	Save/G/P=home Colu_dis_2 as "Colu_dis_"+num2str(cl)+"_LaLa_2.txt"	
	Save/G/P=home Colu_dis_3 as "Colu_dis_"+num2str(cl)+"_LaLa_3.txt"	
	//Save/G/P=home Colu_dis_4 as "Colu_dis_"+num2str(cl)+"_LaLa_4.txt"	
	//Save/G/P=home Colu_dis_5 as "Colu_dis_"+num2str(cl)+"_LaLa_5.txt"	
	//Save/G/P=home Colu_dis_6 as "Colu_dis_"+num2str(cl)+"_LaLa_6.txt"	
end


function GauInteg_Mn(cor_stack, yW_stack, xW_stack, A_stack, z0_stack, size_x, size_y)		
	wave cor_stack, yW_stack, xW_stack, A_stack, z0_stack		
	variable size_x, size_y	
	variable num_layers = DimSize(xW_stack,1)
	make /O/D/N=(1,num_layers) Mn_Inte
	variable j=0
	for (j=0; j<num_layers; j+=1)
		Mn_Inte[1][j]=z0_stack[2][j]*size_x*size_y+2*A_stack[2][j]*sqrt(1-cor_stack[2][j]*cor_stack[2][j])*pi*abs(xW_stack[2][j])*abs(yW_stack[2][j])
	endfor
end

function Column_dis_La(x0_stack, y0_stack)
	wave x0_stack, y0_stack
	variable num_layers = DimSize(x0_stack,1)
	make /O/D/N=(1,num_layers) Colu_dis_0, Colu_dis_1, Colu_dis_2, Colu_dis_3 //Colu_dis_4, Colu_dis_5, Colu_dis_6
	variable n=0
	for (n=0; n<num_layers; n+=1)
		Colu_dis_0[1][n]=sqrt((x0_stack[0][n]-x0_stack[3][n])^2+(y0_stack[0][n]-y0_stack[3][n])^2)			//adding for LaLa_0 on 08/25/2014 
	endfor
	for (n=0; n<num_layers; n+=1)
		Colu_dis_1[1][n]=sqrt((x0_stack[2][n]-x0_stack[5][n])^2+(y0_stack[2][n]-y0_stack[5][n])^2)			//adding for LaLa_1 on 08/25/2014 
	endfor
	for (n=0; n<num_layers; n+=1)
		Colu_dis_2[1][n]=sqrt((x0_stack[0][n]-x0_stack[2][n])^2+(y0_stack[0][n]-y0_stack[2][n])^2)			//adding for LaLa_2 on 08/25/2014 
	endfor
	for (n=0; n<num_layers; n+=1)
		Colu_dis_3[1][n]=sqrt((x0_stack[3][n]-x0_stack[5][n])^2+(y0_stack[3][n]-y0_stack[5][n])^2)			//adding for LaLa_3 on 08/25/2014 
	endfor
	//for (n=0; n<num_layers; n+=1)
	//	Colu_dis_4[1][n]=sqrt((x0_stack[0][n]-x0_stack[2][n])^2+(y0_stack[0][n]-y0_stack[2][n])^2)			//adding for LaLa_4 on 08/25/2014 
	//endfor
	//for (n=0; n<num_layers; n+=1)
	//	Colu_dis_5[1][n]=sqrt((x0_stack[3][n]-x0_stack[5][n])^2+(y0_stack[3][n]-y0_stack[5][n])^2)			//adding for LaLa_5 on 08/25/2014 
	//endfor
	//for (n=0; n<num_layers; n+=1)
	//	Colu_dis_6[1][n]=sqrt((x0_stack[6][n]-x0_stack[8][n])^2+(y0_stack[6][n]-y0_stack[8][n])^2)			//adding for LaLa_6 on 08/25/2014 
	//endfor
end