#pragma rtGlobals=1		// Use modern global access method.
// Requires to be opened with some Igor procedures as below.
// Open Annular Average2.ipf, DP Annular Average UI v1.3.ipf and 
// Auto Find DP Center2.ipf
//**************************************************************************
//  First load .ser file by using SER file loader
//  For STEM FEM data, there are many DP images in single .SER file
//  It calls AnnularAverageCorners to calculate the variance of diffraction data
//  Also the center of DP should be found before this function is called
//  which is AutoCenterFind(im, xmin, ymin, xmax, ymax, size, radius, angle)
//  Need to know the beam stop range (x_min, y_min) to (x_max, y_max)
//  This is done through calling BlockBeamStop, but this function does not
//  output the cursor value, therefore, we may need to add a output statement
//  to give (x_min, y_min) to (x_max, y_max)
//  If we need to know the number of pixels of certain ring, we have to comment 
//  Killwaves/Z npix in function AnnularAverage in Annular Average2.ipf procedure
//**********************************************************
//Contains function
//function Calc_STEM_Var(dat, x_c, y_c, strip_width, x_rect_min, x_rect_max, y_rect_min, y_rect_max)
//this function calculates the variance in input data file, which is 3-D wave, the layer number indicates
//the image index in *.SER file
//(x_c,y_c) is the center of DP pattern in each image
//strip_width is integration stip width, normally 1 or 2 in pixel
//x_rect_min, x_rect_max the rectangle beam-stopper minimum and maximum x coordinates
//y_rect_min, y_rect_max the rectangle beam-stopper minimum and maximum y coordinates
//it gives var_dp wave, which stores the variance results

//change log
// This procedure was written by Feng Yi on 12/14/2009
// This procedure has been rewritten by Jinwoo Hwang to do a smart thickness filtering  - as of 05/31/2010
//03/28/11 function Calc_STEM_Var_bias_corr was added. It uses a correct bias image to fix the errors in the bias images acquired before Sep. 2010 in Titan. -JWH 
//03/28/11 function bias_subtraction was added. subtract a bias image from a stack of images. -JWH
//03/28/11 function thickness_count was found. It looks like it is just a copy of calc_stem_var. -JWH
// 04/18/13 bias correction and other currently unused functions removed, Calc_STEM_Var cleaned up a bit.  PMV


//Arguments
// dat = 3D stack of nanodiffraction patterns
// name_var  = output basename string
// guide_img = HAADF simultaneous acquired image
// x_c , y_c = center of diffraction pattern in pixels
// strip_width = width of the bins in the annular average, typically 2
// x_rect_min, x_rect_max, y_rect_min, y_rect_max = box for the beam stop
// thickness_start = minimum thickness in the guide image in nm
// thickness_end = maximum thickness in the guide image in nm
// slope = proptionality between HAADF counts and thickness
// offset = black level in HAADF image
// allow = thickness range in nm to be allowed into variance calculation
function Calc_STEM_Var(dat, name_var, guide_img, x_c, y_c, strip_width, x_rect_min, x_rect_max, y_rect_min, y_rect_max, thickness_start, thickness_end, slope, offset, allow)
	wave dat
	string name_var
	wave guide_img
	variable x_c, y_c, strip_width
	variable x_rect_min, x_rect_max, y_rect_min, y_rect_max, thickness_start, thickness_end, slope, offset, allow
	
	
	// required calibration or global variables
	variable G_Gain = 7.46 // CCD camera gain (4.72 cts/primary electron for UW Titan, 7.46 cts/elec for ER-C Titan-S)
	variable kcal_base = 0.03516  // width of one pixel in the DP in inverse units (1/Angstroms or 1/nm)
								   // 0.002727985 1/Ang for UW Titan at binning 4, CL=512 mm, EFSTEM mode
								   // 0.01758 1/nm for Julich Titan-S at binning 1, CL= 	510 mm, EFSTEM mode
	
	variable kcal = kcal_base*strip_width
	variable num_images
	variable num_dps
	variable index_image, i, j, k, ii
	variable x_dim, y_dim  //x and y dimension of wave dat
	variable radius
	string name_ave  //annular average wave name
	string name_var0
	
	name_var0 = name_var
	
	
	make/o annular_av, annular_int, npix //dummy waves
	
	num_images=DimSize(dat, 2)  //how many DP images in wave dat
	x_dim=DimSize(dat, 0)
	y_dim=DimSize(dat, 1)
	
	
	num_dps = DimSize(guide_img, 0) * DimSize(guide_img, 1)   //total pixel numbers in the guide image
	print num_images, num_dps
	make/o/n=(num_dps) guide_pix_int
	
	k=0
	for (i=0; i<DimSize(guide_img, 1); i+=1)
		for (j=0; j<DimSize(guide_img, 0); j+=1)	
			guide_pix_int[k] = guide_img[j][i]
			k+=1
		endfor
	endfor
	
	variable int_min, int_max, thickness
	
	for (ii=thickness_start; ii<thickness_end; ii+=2*allow)
	
		thickness = ii
		print "thickness=",  thickness
			
		int_min = (thickness - allow)*slope + offset
		int_max = (thickness + allow)*slope +offset
		
		print "int_min=", int_min
		print "int_max=", int_max
			
		make/o/n=(x_dim,y_dim) temp_image

		//first, calculate the annular average and integrated intensity for the 1st image
		i=0
		temp_image=dat[p][q][i]
		temp_image[x_rect_min,x_rect_max][y_rect_min,y_rect_max] = nan
		radius = AnnularAverage(temp_image, x_c, y_c, strip_width)
		radius = AnnularIntegral(temp_image, x_c, y_c, strip_width)
		
		
		duplicate/o annular_av, i_square, i_ave, var_dp, i_total 
		//i_squre, i_ave are i^2 and i per pixel over  many diffraction patterns
		//i_total is the integrated i at certain radius over many diffraction pattern
		i_square=0
		i_ave=0
		var_dp=0
		i_total=0
		
		variable img_count	
		
		img_count = 0
		for (i=0; i<num_images; i+=1)
			
			if(guide_pix_int[i] >= int_min)
				if(guide_pix_int[i] < int_max)
				//if(i!=18)
		 			name_ave = NameofWave(dat)
		 			name_ave=name_ave+"_ave_"+num2str(i)
					temp_image=dat[p][q][i]
					//setup the beam stopper area to NaN
					temp_image[x_rect_min,x_rect_max][y_rect_min,y_rect_max] = nan
			
					//calculate the annular average excluding beam stopper area
					//and stores the annular avarage data in annular_av
					radius = AnnularAverage(temp_image, x_c, y_c, strip_width)
					radius = AnnularIntegral(temp_image, x_c, y_c, strip_width)
					

					
					i_square = i_square + annular_av^2
					i_ave = i_ave + annular_av
					i_total = i_total + annular_int
					duplicate/o annular_av, $name_ave
					SetScale/P x 0,kcal,"", $name_ave
					img_count +=1 
				//endif
				endif
			endif
					
		endfor
	
		print "number of dps chosen=", img_count
		i_square = i_square / img_count
		i_ave = i_ave /img_count
		i_total = i_total/img_count
	
		
		// this is the correct way to calculate the shot noise
		var_dp =i_square/i_ave^2-1-G_Gain/i_ave/npix
	
		SetScale/P x 0,kcal,"", var_dp
		name_var = name_var0 + "_" + num2str(thickness) + "nm_"+ num2str(img_count) + "ea"
		duplicate/o var_dp, $name_var
		


	endfor
	
	//kill unnessary waves
	killwaves/z  annular_av, i_square, i_ave, i_total, var_dp
	killwaves/z annular_int, npix
	killwaves/z temp_image

end
