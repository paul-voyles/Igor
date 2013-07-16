#pragma rtGlobals=1		// Use modern global access method.
//
// Code for computing angular correlations in nanodiffraction patterns.
// v1.0 06-03-10 pmv
// v1.1 12-05-10 pmv  added elliptical diffraction distortion correction to
//                             ImageC2P, and supporting functions ellipX, ellipY
//                             and EllipticalSampling
//
//        12-10-10 pmv  fixed bug so that StackPowerSpectrumVariance actually
//                             calculates the variance, not the mean
//
// ImageC2P: resample an image from Cartesian (x,y) coordinates to 
// polar (r, theta) coordinates using interp2D
//
// AngularCorrelationAllK: compute the (theta, theta+delta) angular autocorrelation
// function, normalized as in Wochner et al. PNAS 106, 11511 (2009)


// take an image with sampling of square pixels (x,y) and
// resample it to (r, theta).  Preserves pixels set to NaN to block the beam stop.
// theta range is always 0 to 2 Pi, in radians.

//change log
//03/28/11 Function StackImageC2P_thickness was added. It does ImageC2P with thickness filtering -JWH
//03/28/11 Function remove_kt_nz was added. It removes blank images from a stack of kt images. -JWH
//03/28/11 Function remove_ps_nz was added. It removes blank images from a stack of ps images. -JWH
//03/28/11 Function stack_images was added. It manually stacks images -JWH
//03/28/11 Function stack_images2 was added It manually stacks images. -JWH
//03/28/11 Function extract_image was added It extracts an image from a stack of images. -JWH
//03/28/11 Function random_disk was added It simulates randomly distributed speckle pairs. Used for brute-force simulation of speckles. -JWH
//03/28/11 Function random_disk_gaussian was added. Same as random_disk function but speckles are gaussians. -JWH
//03/28/11 Function AngularCorrelationAllK_debug was added. Used for debugging. No bugs found -JWH
//03/28/11 Functions counthighfold and counthigh2fold were added. They count the number of ACs that have values high signals in a certain angle range. -JWH
//03/28/11 Function function cosine was added. It makes cosine waves. -JWH
//03/28/11 Function beamstop was added. It masks the beam stop in an image stack. -JWH
//03/28/11 Function  aligned_v2d was added. It rotates and aligns DPs at the brightest speckles of each image before calculating V2D of an image stack -JWH.





function ImageC2P(im, rmin, rmax, rstep, tstep, ellip, t_off)
	wave im
	variable rmin, rmax, rstep, tstep, ellip, t_off
	
	variable nr = floor( (rmax-rmin)/rstep )
	variable nt = floor( 2*Pi  / tstep )

	// ensure that the number of points in k and theta are both even for later FT
	if (mod(nr, 2))
		nr += 1
	endif
	if(mod(nt, 2))
		nt += 1
	endif
	
	make/d/o/n=(nr, nt) im_polar
	setscale/I x rmin, rmax, im_polar
	setscale/I y 0, 2*Pi, im_polar

	//im_polar = interp2d(im, ellip*x*cos(y), x*sin(y))	
	im_polar = interp2d(im, ellipX(x, y, ellip, t_off),  ellipY(x, y, ellip, t_off))
end

// ellipX and ellipY translate (r, theta) to (x, y), respectively, with elliptical distortion. 
// ellip is the ratio of the major to minor axis lengths.  t_off is the rotation angle of the 
// ellipse axes with respect to the coordinate axes.
function ellipX(r, theta, ellip, t_off)
	variable r, theta, ellip, t_off
	
	return cos(t_off)*(2*r*ellip / (1+ellip))*cos(theta) - sin(t_off)*(2*r/(1+ellip))*sin(theta)
	//return cos(t_off)*ellip*r*cos(theta) - sin(t_off)*r*sin(theta)
	
end

function ellipY(r, theta, ellip ,t_off)
	variable, r, theta, ellip, t_off

	return sin(t_off)*(2*r*ellip / (1+ellip))*cos(theta) + cos(t_off)*(2*r/(1+ellip))*sin(theta)
	//return sin(t_off)*ellip*r*cos(theta)+cos(t_off)*r*sin(theta)
	
end

// generates the (x,y) grid of points with elliptical distortion for a particular (r, theta), 
// ellipticity, and rotation angle t_off.  Use it to overlay the grid of points on the data to
// work out the ellipticity and rotation angle.
function EllipseSampling(rmin, rmax, rstep, tstep, ellip, t_off)
	variable rmin, rmax, rstep, tstep, ellip, t_off
	
	variable nr = floor( (rmax-rmin)/rstep )
	variable nt = floor( 2*Pi  / tstep )

	// ensure that the number of points in k and theta are both even for later FT
	if (mod(nr, 2))
		nr += 1
	endif
	if(mod(nt, 2))
		nt += 1
	endif

	Make/o/n=(nr*nt) ellipse_x, ellipse_y
	
	ellipse_x = ellipX(mod(p, nr)*rstep - rmin, floor(p/nr)*(2*Pi/nt), ellip, t_off)
	ellipse_y = ellipY(mod(p, nr)*rstep - rmin, floor(p/nr)*(2*Pi/nt), ellip, t_off)

end	

// Computes the angular autocorrelation function C(k, delta), which delta
// is the angle between theta and theta prime.  Ignores pixels set to NaN
// to account for the beam stop.	
function AngularCorrelationAllK(imp)
	wave imp
	
	duplicate/O imp imp_rot
	duplicate/O imp, ang_c
	Redimension/D ang_c
	
	// angular correlation: rotate the image by single angle steps
	variable i, j
	for(i=0; i<DimSize(imp, 1); i+=1)
		imp_rot = imp
		Imagetransform/G=(i) rotateCols imp_rot
		Matrixop/o overlap = imp_rot*imp
		
		for(j=0; j<DimSize(imp, 0); j+=1)
			Imagetransform/g=(j) getrow overlap
			wave one_k = $"W_ExtractedRow"
			wavestats/m=1/q one_k
			if(V_sum > 0 && V_npnts > 0)		
				ang_c[j][i] = V_sum / V_npnts
			else
				//printf "Not enough data.\r"
				ang_c[j][i] = 0
			endif
		endfor
	
	endfor
	
	// normalization by unrotated image
	variable n
	//Matrixop/o overlap = imp*imp
	for(j=0; j<DimSize(imp, 0); j+=1)
		Imagetransform/g=(j) getrow imp
		wave one_k = $"W_ExtractedRow"
		wavestats/m=1/q one_k
		if(V_sum > 0 && V_npnts > 0)		
			n = V_sum / V_npnts
			ang_c[j][] = (ang_c[j][q] - n^2)/n^2
		else
			//printf "Not enough data.\r"
			ang_c[j][] = 0
		endif
	endfor
	
	Killwaves imp_rot, overlap, W_ExtractedRow
end


function AngularPowerSpectrum(dat)
	wave dat
	
	FFT/rows/dest=ft_tmp dat
	Duplicate/O ft_tmp ang_ps
	Redimension/R ang_ps
	ang_ps = magsqr(ft_tmp)
	setscale/p y 0, 1, ang_ps
	killwaves ft_tmp
end


function StackImageC2P(imst,  rmin, rmax, rstep, tstep, ellip, t_off)
	wave imst
	variable rmin, rmax, rstep, tstep, ellip, t_off
	
	ImageTransform/P=0 getPlane imst
	wave imt = $"M_ImagePlane"
	
	ImageC2P(imt, rmin, rmax, rstep, tstep, ellip, t_off)
	wave imp = $"im_polar"
	make/O/D/N=(DimSize(imp, 0), DimSize(imp, 1), DimSize(imst, 2)) stack_polar
	SetScale/P x DimOffSet(imp, 0), DimDelta(imp, 0), "", stack_polar
	SetScale/P y DimOffset(imp, 1), DimDelta(imp, 1), "", stack_polar
	stack_polar[][][0] = imp[p][q]
	
	variable i
	for(i=1; i<DimSize(imst,2); i+=1)
		ImageTransform/P=(i) getPlane imst
		ImageC2P(imt, rmin, rmax, rstep, tstep, ellip, t_off)
		stack_polar[][][i] = imp[p][q]
	endfor
	
	string name_wave
	name_wave = nameofwave(imst)
	name_wave = name_wave+"_kt"
	duplicate/o stack_polar $name_wave
	
	Killwaves imt, imp, stack_polar
	
end


function StackImageC2P_thickness(imst, guide_img,  rmin, rmax, rstep, tstep, ellip, t_off, thickness_start, thickness_end, slope, offset, allow)
	wave imst
	wave guide_img
	variable rmin, rmax, rstep, tstep, ellip, t_off
	variable  thickness_start, thickness_end, slope, offset, allow
	
	ImageTransform/P=0 getPlane imst
	wave imt = $"M_ImagePlane"
	
	ImageC2P(imt, rmin, rmax, rstep, tstep, ellip, t_off)
	wave imp = $"im_polar"
	make/O/D/N=(DimSize(imp, 0), DimSize(imp, 1), DimSize(imst, 2)) stack_polar
	SetScale/P x DimOffSet(imp, 0), DimDelta(imp, 0), "", stack_polar
	SetScale/P y DimOffset(imp, 1), DimDelta(imp, 1), "", stack_polar
	
	variable i, j, k
	make/o/n=(DimSize(imst, 2)) guide_pix_int
	k=0
	for (i=0; i<DimSize(guide_img, 1); i+=1)
		for (j=0; j<DimSize(guide_img, 0); j+=1)	
			guide_pix_int[k] = guide_img[j][i]
			k+=1
		endfor
	endfor
	
	if( guide_pix_int[0] >=	(thickness_start - allow)*slope + offset)  
		if( (thickness_end + allow)*slope + offset > guide_pix_int[0])  
			//print (thickness_start - allow)*slope + offset,  guide_pix_int[0], (thickness_end + allow)*slope + offset
			stack_polar[][][0] = imp[p][q]
		endif
	endif	
	if( guide_pix_int[0] <	(thickness_start - allow)*slope + offset)  
		//print "else", (thickness_start - allow)*slope + offset,  guide_pix_int[0], (thickness_end + allow)*slope + offset
		stack_polar[][][0] = nan	
	endif	
	if( (thickness_end + allow)*slope + offset <= guide_pix_int[0])  
		//print "else", (thickness_start - allow)*slope + offset,  guide_pix_int[0], (thickness_end + allow)*slope + offset
		stack_polar[][][0] = nan
	endif
		
	for(i=1; i<DimSize(imst,2); i+=1)
	if( guide_pix_int[i] >=	(thickness_start - allow)*slope + offset)  
		if( (thickness_end + allow)*slope + offset > guide_pix_int[i]) 
			//print (thickness_start - allow)*slope + offset,  guide_pix_int[i], (thickness_end + allow)*slope + offset
			ImageTransform/P=(i) getPlane imst
			ImageC2P(imt, rmin, rmax, rstep, tstep, ellip, t_off)
			stack_polar[][][i] = imp[p][q]
		endif
	endif	
	if( guide_pix_int[i] <	(thickness_start - allow)*slope + offset)  
		//print "else", (thickness_start - allow)*slope + offset,  guide_pix_int[i], (thickness_end + allow)*slope + offset
		stack_polar[][][i] = nan	
	endif
	if( (thickness_end + allow)*slope + offset <= guide_pix_int[i])  
		//print "else", (thickness_start - allow)*slope + offset,  guide_pix_int[i], (thickness_end + allow)*slope + offset
		stack_polar[][][i] = nan	
	endif		
	endfor
	
	string name_wave
	name_wave = nameofwave(imst)
	name_wave = name_wave+"_kt"
	duplicate stack_polar $name_wave
	
	Killwaves imt, imp
	
end




		
function StackAngularCorrelationAllK(imst)
	wave imst
	
	Duplicate/O imst stack_ac
	
	variable i 
	for(i=0; i<DimSize(imst, 2); i+=1)
		ImageTransform/P=(i) getPlane imst
		wave imp = $"M_ImagePlane"
		AngularCorrelationAllK(imp)
		wave ang_c = $"ang_c"
		stack_ac[][][i] = ang_c[p][q]
	endfor
	
	string name_wave
	name_wave = nameofwave(imst)
	name_wave = name_wave+"_ac"
	duplicate stack_ac $name_wave
	
	Killwaves imp, ang_c, stack_ac
	
end	

function StackAngularPowerSpectrum(imst)
	wave imst
	
	ImageTransform/p=0 getPlane imst
	wave im = $"M_ImagePlane"
	AngularPowerSpectrum(im)
	wave ang_ps = $"ang_ps"
	
	Make/O/D/N=(DimSize(ang_ps, 0), DimSize(ang_ps, 1), DimSize(imst, 2)) stack_ps
	SetScale/P x DimOffset(ang_ps, 0), DimDelta(ang_ps, 0), stack_ps
	SetScale/P y DimOffset(ang_ps, 1), DimDelta(ang_ps, 1), stack_ps
	stack_ps[][][0] = ang_ps[p][q]
	
	variable i
	for(i=1; i<DimSize(imst, 2); i+=1)
		ImageTransform/p=(i) getPlane imst
		AngularPowerSpectrum(im)
		stack_ps[][][i] = ang_ps[p][q]
	endfor
	
	string name_wave
	name_wave = nameofwave(imst)
	name_wave = name_wave+"_ps"
	duplicate stack_ps $name_wave
	
	Killwaves im, ang_ps
	
end

function remove_kt_nz(dat)
	wave dat
	
	variable i, j, k
	variable img_sum
	variable img_count
	img_sum = 0
	img_count = 0
	for(k=0; k<DimSize(dat, 2); k+=1)	
		//for(i=0; i<DimSize(dat, 0); i+=1)   
		for(i=80; i<82; i+=1) 
			//for(j=0; j<DimSize(dat, 1); j+=1)   
			for(j=300; j<302; j+=1) 
				img_sum +=dat[i][j][k]
			endfor
		endfor
		//print img_sum
		if(img_sum == img_sum) 
			img_count +=1		
		endif		
		img_sum = 0
	endfor	
	print img_count	
	make/o/d/n=(DimSize(dat, 0), DimSize(dat, 1), img_count) dat2
	
	img_sum = 0
	img_count = 0
	for(k=0; k<DimSize(dat, 2); k+=1)	
		//for(i=0; i<DimSize(dat, 0); i+=1)   
		for(i=80; i<82; i+=1) 
			//for(j=0; j<DimSize(dat, 1); j+=1)   
			for(j=300; j<302; j+=1) 
				img_sum +=dat[i][j][k]
			endfor
		endfor
		if(img_sum == img_sum)    
			dat2[][][img_count] = dat[p][q][k]	
			img_count +=1
		endif
		img_sum = 0
	endfor	
	
	string name_wave
	name_wave = nameofwave(dat)
	name_wave = name_wave+"_nz"
	duplicate dat2 $name_wave	 

end	



function remove_ps_nz(dat)
	wave dat
	
	variable i, j, k
	variable img_sum
	variable img_count
	img_sum = 0
	img_count = 0
	for(k=0; k<DimSize(dat, 2); k+=1)	
		for(i=60; i<80; i+=1) 
			for(j=350; j<380; j+=1) 
				img_sum +=dat[i][j][k]
			endfor
		endfor
		//if(img_sum !=0) 
		if(img_sum==img_sum)
			img_count +=1		
		endif		
		img_sum = 0
	endfor	
		
	make/o/d/n=(DimSize(dat, 0), DimSize(dat, 1), img_count) dat2
	
	img_sum = 0
	img_count = 0
	for(k=0; k<DimSize(dat, 2); k+=1)	
		for(i=60; i<80; i+=1) 
			for(j=350; j<380; j+=1) 
				img_sum +=dat[i][j][k]
			endfor
		endfor
		//if(img_sum !=0)   
		if(img_sum==img_sum)
			dat2[][][img_count] = dat[p][q][k]	
			img_count +=1
		endif
		img_sum = 0
	endfor	
	
	string name_wave
	name_wave = nameofwave(dat)
	name_wave = name_wave+"_nz"
	duplicate dat2 $name_wave	 

end	



function StackPowerSpectrumVariance(stack_ps)
	wave stack_ps
	
	make/o/d/n=(DimSize(stack_ps, 0), DimSize(stack_ps, 1)) stack_ps_v
	SetScale/P x DimOffset(stack_ps, 0), DimDelta(stack_ps, 0), stack_ps_v
	SetScale/P y DimOffSet(stack_ps, 1), DimOffSet(stack_ps, 1), stack_ps_v
	
	make/o/d/n=(DimSize(stack_ps, 0), DimSize(stack_ps, 1)) stack_ps_avg
	SetScale/P x DimOffset(stack_ps, 0), DimDelta(stack_ps, 0), stack_ps_avg
	SetScale/P y DimOffSet(stack_ps, 1), DimOffSet(stack_ps, 1), stack_ps_avg
	
	make/o/d/n=(DimSize(stack_ps, 0), DimSize(stack_ps, 1)) stack_ps_sdev
	SetScale/P x DimOffset(stack_ps, 0), DimDelta(stack_ps, 0), stack_ps_sdev
	SetScale/P y DimOffSet(stack_ps, 1), DimOffSet(stack_ps, 1), stack_ps_sdev
	
	variable i, j
	for(i=0; i<DimSize(stack_ps, 0); i+=1)
		for(j=0; j<DimSize(stack_ps, 1); j+=1)
			ImageTransform/beam={(i), (j)} getBeam stack_ps
			wave b = $"W_Beam"
			wavestats/q b
			stack_ps_v[i][j] = V_sdev^2	
			stack_ps_avg[i][j] = V_avg
			stack_ps_sdev[i][j] = V_sdev
		endfor
	endfor
	
	string name_wave1
	name_wave1 = nameofwave(stack_ps)
	name_wave1 = name_wave1+"_v"
	duplicate stack_ps_v $name_wave1
	
	string name_wave2
	name_wave2 = nameofwave(stack_ps)
	name_wave2 = name_wave2+"_avg"
	duplicate stack_ps_avg $name_wave2
	
	string name_wave3
	name_wave3 = nameofwave(stack_ps)
	name_wave3 = name_wave1+"_sdev"
	duplicate stack_ps_sdev $name_wave3
end

// stacks stacks of images, but doesn't work when a stack contains only one image
function stack_images(dat1, dat2)   
	wave dat1, dat2
	variable num_layers, i, j
	num_layers = dimsize(dat1, 2) + dimsize(dat2, 2)
	make/o/d/n=(DimSize(dat1, 0), DimSize(dat1, 1), num_layers) temp
	
	for(i=0; i<DimSize(dat1, 2); i+=1)
		temp[][][i] = dat1[p][q][i]
		//print i
	endfor
	//print dimsize(dat1,2), dimsize(dat2,2)
	for(i=DimSize(dat1, 2); i<num_layers; i+=1)
		j = i-dimsize(dat1, 2)
		//print i,j
		temp[][][i] = dat2[p][q][j]
	endfor
	duplicate/o temp temp1
end	


//revising to handle one-image-stack - error when stacking 2+1
function stack_images2(dat1, dat2)
	wave dat1, dat2
	variable num_layers, i, j, k
	variable num_layer1, num_layer2

	if(dimsize(dat1, 2) == 0)
		num_layer1 = 1
	else
		num_layer1 = dimsize(dat1, 2)
	endif
	if(dimsize(dat2, 2) == 0)
		num_layer2 = 1
	else
		num_layer2 = dimsize(dat2, 2)
	endif
	//num_layers = dimsize(dat1, 2) + dimsize(dat2, 2)
	num_layers = num_layer1 + num_layer2

	print num_layer1 , num_layer2

	make/o/d/n=(DimSize(dat1, 0), DimSize(dat1, 1), num_layers) temp   //somehow "stack of 3 gives error"

	//for(i=0; i<DimSize(dat1, 2); i+=1)
	for(i=0; i<num_layer1; i+=1)
		temp[][][i] = dat1[p][q][i]
		//print i
	endfor

	//for(i=DimSize(dat1, 2); i<num_layers; i+=1)
	for(i=num_layer1; i<num_layers; i+=1)
		j = i-num_layer1
		//print i,num_layer1, j
		temp[][][i] = dat2[p][q][j]
	endfor
	duplicate/o temp temp1

	


end	



function extract_image(dat, img_num)
	wave dat
	variable img_num
	variable num_layers, i, j
	num_layers = dimsize(dat, 2) 
	make/o/d/n=(DimSize(dat, 0), DimSize(dat, 1)) temp
	
	for(i=0; i<DimSize(dat, 0); i+=1)
		for(j=0; j<DimSize(dat, 1); j+=1)
			temp[i][j] = dat[i][j][img_num]
		endfor
	endfor


end	

function random_disk(dat)
	wave dat
	variable center_x, center_y, pix_x, pix_y
	variable i, j, k
	variable dist, radius
	variable num_disk, kr, kr_pix
	variable dp_int
	variable rand0, rand_num, kx, ky
	variable ii
	variable num_dp
	string name_wave
	
	rand0 = abs(floor(enoise(10000)))
	//print rand0
	SetRandomSeed rand0

//-------------------------------------------------------------------------------
	num_disk = 12	 // number of disk pairs
	kr = 0.4  //inv ang
	num_dp = 20
//--------------------------------------------------------------------------------

	center_x = 255
	center_y = 255
	pix_x = 512
	pix_y = 512
	radius = 5
	kr_pix = floor(kr/0.00545597)
	dp_int = 1000
	

for(ii=0; ii<num_dp; ii+=1)	
	
	dat = 0
	dat[center_x][center_y] = dp_int
	
	//for(i=0; i<pix_x; i+=1)
	//	for(j=0; j<pix_y; j+=1)
	//		dist = sqrt( (i-center_x)^2 + (j-center_y)^2 )
	//		if(dist<=radius)
	//			dat[i][j] = dp_int	
	//		endif		
	//	endfor
	//endfor
	

		
	make/o/d/n=(num_disk) rand_wave
	for(i=0; i<num_disk; i+=1)	
		rand_num = enoise(pi/2.0) +pi/2.0    // enoise for even distribution, gnoise for gaussain 
		//rand_wave[i] = rand_num
		kx = cos(rand_num)*kr_pix
		ky = abs(sqrt(kr_pix^2-kx^2))
		//print kx, ky
		dat[kx+center_x][ky+center_y] = dp_int		
		dat[-kx+center_x][-ky+center_y] = dp_int	
		
		for(j=0; j<pix_x; j+=1)
			for(k=0; k<pix_y; k+=1)
				dist = sqrt( (j- (kx+center_x) )^2 + (k- (ky+center_y) )^2 )
				if(dist<=radius)
					dat[j][k] = dp_int	
				endif	
				dist = sqrt( (j- (-kx+center_x) )^2 + (k- (-ky+center_y) )^2 )
				if(dist<=radius)
					dat[j][k] = dp_int	
				endif		
			endfor
		endfor	
		
	endfor	

	name_wave = "dp_"+num2str(ii)
	duplicate/o dat $name_wave
	
	
endfor			
	
end	


function random_disk_gaussian(dat)
	wave dat
	variable center_x, center_y, pix_x, pix_y
	variable i, j, k
	variable dist, fwhm, sigma
	variable num_disk, kr, kr_pix
	variable dp_int
	variable rand0, rand_num, kx, ky
	variable ii
	variable num_dp
	string name_wave
	
	rand0 = abs(floor(enoise(10000)))
	//print rand0
	SetRandomSeed rand0

//-------------------------------------------------------------------------------
	num_disk = 12	 // number of disk pairs
	kr = 0.4  //inv ang
	num_dp = 500
//--------------------------------------------------------------------------------

	center_x = 255
	center_y = 255
	pix_x = 512
	pix_y = 512
	fwhm = 9    //in pixel
	sigma = fwhm/2.35482
	kr_pix = floor(kr/0.00545597)
	dp_int = 500
	

for(ii=0; ii<num_dp; ii+=1)	
	
	dat = 200
	dat[center_x][center_y] = dp_int
	
	
	make/o/d/n=(num_disk) rand_wave
	
	variable rand_num_2

	for(i=0; i<num_disk; i+=1)	
		rand_num = enoise(pi/2.0) +pi/2.0    // enoise for even distribution, gnoise for gaussain
		rand_num_2 =  abs(enoise(1))
		rand_wave[i] = rand_num_2
		kx = cos(rand_num)*kr_pix
		ky = abs(sqrt(kr_pix^2-kx^2))
		//print kx, ky
		//dat[kx+center_x][ky+center_y] = dp_int		
		//dat[-kx+center_x][-ky+center_y] = dp_int	
		
		for(j=0; j<pix_x; j+=1)
			for(k=0; k<pix_y; k+=1)
				dist = sqrt( (j- (kx+center_x) )^2 + (k- (ky+center_y) )^2 )
				if(dist<= 3.0*fwhm/2)
					dat[j][k] = dat[j][k] + rand_num_2 *dp_int * exp( -(  ((j- (kx+center_x) )^2)/(2*sigma^2) + ((k- (ky+center_y) )^2)/(2*sigma^2)  )  )	
				endif	
				dist = sqrt( (j- (-kx+center_x) )^2 + (k- (-ky+center_y) )^2 )
				if(dist<= 3.0*fwhm/2)
					dat[j][k] = dat[j][k] + (1-rand_num_2) * dp_int * exp( -(  ((j- (-kx+center_x) )^2)/(2*sigma^2) + ((k- (-ky+center_y) )^2)/(2*sigma^2)  )  )
				endif		
			endfor
		endfor	
		
	endfor	

	name_wave = "dp_"+num2str(ii)
	duplicate/o dat $name_wave
	
	
endfor			
	
end	


function AngularCorrelationAllK_debug(imp)
	wave imp
	
	duplicate/O imp imp_rot
	duplicate/O imp, ang_c
	Redimension/D ang_c
	
	ang_c= 0
	
	// angular correlation: rotate the image by single angle steps
	variable i, j
//	for(i=0; i<DimSize(imp, 1); i+=1)
	i=100
		imp_rot = imp
		Imagetransform/G=(i) rotateCols imp_rot   
		Matrixop/o overlap = imp_rot*imp   
		
//		for(j=0; j<DimSize(imp, 0); j+=1)
		j=43
			Imagetransform/g=(j) getrow overlap
			//wave one_k = $"W_ExtractedRow"   //debug
			duplicate/o $"W_ExtractedRow" one_k    //debug
			
			//debug
			if(i==100)
				duplicate/o one_k one_k_100
			endif
			if(i==528)
				duplicate/o one_k one_k_528
			endif	
			
			
			wavestats/m=1/q one_k
			if(V_sum > 0 && V_npnts > 0)		
				ang_c[j][i] = V_sum    //debug

				//ang_c[j][i] = V_sum / V_npnts  //debug
			else
				//printf "Not enough data.\r"
				ang_c[j][i] = 0
			endif
//		endfor
//	endfor
	
	//print ang_c[43][100]     //debug
	//print [43][528]
	
//	// normalization by unrotated image
//	variable n
//	//Matrixop/o overlap = imp*imp
//	for(j=0; j<DimSize(imp, 0); j+=1)
//		Imagetransform/g=(j) getrow imp
//		wave one_k = $"W_ExtractedRow"
//		wavestats/m=1/q one_k
//		if(V_sum > 0 && V_npnts > 0)		
//			n = V_sum / V_npnts
//			ang_c[j][] = (ang_c[j][q] - n^2)/n^2
//		else
//			//printf "Not enough data.\r"
//			ang_c[j][] = 0
//		endif
//	endfor
	
//	Killwaves imp_rot, overlap, W_ExtractedRow
end



function counthighfold(dat_ps, dat_ac, pix_x, pix_y, thrash)
	wave dat_ps
	wave dat_ac
	variable pix_x, pix_y, thrash
	variable i, j, k
	variable n, int, count
	string ps_name, ac_name

	ps_name = "ps_"
	ac_name = "ac_"
	n=dimsize(dat_ps, 2)
	count = 0
	for(i=0; i<n; i+=1)
		int = dat_ps[pix_x][pix_y][i]
		if(int > thrash)
			count = count+1
			extract_image(dat_ps, i)
			duplicate/o temp $ps_name+num2str(count)
			extract_image(dat_ac, i)
			duplicate/o temp $ac_name+num2str(count)
			print "nth image chosen=", i
		endif
	endfor
	imagetransform stackimages  $ps_name+"1"
	duplicate/o M_stack  $ps_name+"stack_filtered"
	killwaves M_stack
	imagetransform stackimages  $ac_name+"1"
	duplicate/o M_stack  $ac_name+"stack_filtered"
	killwaves M_stack
	
end


function counthigh2fold(dat, k_point) 
	wave dat
	variable k_point
	variable i, j, k
	variable count
	
	make/o/d/n=(DimSize(dat, 2)) int_sum, int_avg, int_avg_sq
	
	for(k=0; k<DimSize(dat, 2); k+=1)	
		int_sum[k] = 0
		count = 0
		for(j=0; j<DimSize(dat, 1); j+=1) 
			if(dat[k_point][j][k] == dat[k_point][j][k])
				int_sum[k] += dat[k_point][j][k]
				count +=1
			endif
		endfor
		int_avg[k] = int_sum[k] / count
		int_avg_sq[k] = int_avg[k]^2
	endfor
					
	make/o/d/n=(DimSize(dat, 1), DimSize(dat, 2)) corr
	variable jj
	
	make/o/d/n=(DimSize(dat, 2)) ncount
	ncount = 0
	
	for(k=0; k<DimSize(dat, 2); k+=1)	
		for(j=0; j<DimSize(dat, 1); j+=1) 
			if(dat[k_point][j][k] == dat[k_point][j][k])
			
				jj = j+DimSize(dat, 1)/2
				if(jj > DimSize(dat, 1)-1)
					jj -= DimSize(dat, 1)
				endif
				if(dat[k_point][jj][k] == dat[k_point][jj][k])	
					corr[j][k] = dat[k_point][j][k] * dat[k_point][jj][k]		
					if(corr[j-1][k] < int_avg_sq[k])
						if(corr[j][k] >= int_avg_sq[k])
							ncount[k] +=1	
							endif
					endif								
				endif
			endif
		endfor
	endfor
	
	ncount /=2
	
	variable ncount_avg
	for(k=0; k<DimSize(dat, 2); k+=1)	
		ncount_avg += ncount[k]
	endfor
	ncount_avg /= DimSize(dat, 2)
	print "ncount_avg=", ncount_avg
end


function cosine(freq)
	variable freq
	variable i
	variable num_points
	
	num_points = 628
	
	make/o/d/n=(num_points) cosine_wave
	for(i=0; i<num_points; i+=1)	
		cosine_wave[i] = cos(freq*i/100)
	endfor
		
	setscale/p x 0, 0.01, cosine_wave
	FFT/dest=ft_cosine_wave cosine_wave
	setscale/p x 0, 1, ft_cosine_wave

end



function beamstop(dat)
	wave dat
	variable i, j, k
	variable x_start, x_end, x_width
	variable y_start, y_end
	
	x_start = 252
	x_width = 314-252
	x_end = x_start + x_width
	
	y_start = 233
	y_end = 511
	
	for(k=0; k<DimSize(dat, 2); k+=1)	
		for(i=0; i<DimSize(dat, 0); i+=1) 
			if(i >= x_start)
				if(i <= x_end)
					for(j=0; j<DimSize(dat, 1); j+=1) 
						if(j >= y_start)
							if(j <= y_end)
								dat[i][j][k] = nan
							endif
						endif			
					endfor
				endif
			endif
		endfor
	endfor
end

function aligned_v2d(imp, align_k_pix)
	wave imp
	variable align_k_pix
	
	variable i, j, k
	variable pix_int, pix_int_max, rot_pix
	

	
	duplicate/o imp imp_rot
	
	for(j=0; j<DimSize(imp, 2); j+=1)   // angular loop
		pix_int_max = 0	
		for(i=0; i<DimSize(imp, 1); i+=1)   // angular loop
			pix_int = imp[align_k_pix][i][j] 
			if(pix_int == pix_int)
				if(pix_int > pix_int_max)
					pix_int_max = pix_int
					rot_pix = i
				endif
			endif
		endfor
		extract_image(imp, j)

		duplicate/o temp imp_rot_layer
		Imagetransform/G=(rot_pix) rotateCols imp_rot_layer
		imp_rot[][][j] = imp_rot_layer[p][q]
			
	endfor
	print j, pix_int_max, rot_pix

	
end	