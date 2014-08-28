#pragma rtGlobals=1		// Use modern global access method.
// Functions for analysis of 3D FFTs of model structures
//
// begun 03-05-09 pmv

#include <GizmoSlicer>

// To generate paramfile from FT (512 pixels) for input into IFT program:

//     Import your 3D FT wave you created with Fortran's 3dft.f90 and python's lines_to_3d_wave.py
//         It is encouraged that you rename your FT wave to "ft" because that's way easy to type and you can then just copy and paste the below functions for the most part.
//     Create a new fitting function in Curve Fitting with these parameters:
//         Name:  gauss3df
//         IN ORDER: Fit Coefficients:  sx, sy, sz, cxy, cxz, cyz, x0, y0, z0
//                           Independent Variables:  x, y, z
//         Use the function:  f(x,y,z) = exp(  -1/(2* (-1+cxy^2+cxz^2+cyz^2-2*cxy*cxz*cyz)) * ( (cyz^2-1)*(x-x0)^2/sx^2 + (cxz^2-1)*(y-y0)^2/sy^2 + (cxy^2-1)*(z-z0)^2/sz^2 + (2*cxy*(x-x0)*(y-y0)-2*cxz*cyz*(x-x0)*(y-y0))/(sx*sy) + (2*cxz*(x-x0)*(z-z0)-2*cxy*cyz*(x-x0)*(z-z0))/(sx*sz) + (2*cyz*(y-y0)*(z-z0)-2*cxy*cxz*(y-y0)*(z-z0))/(sy*sz) )  )

//     Run FindSpots(ft, 15)   --- this generates "particles"

//     Run gauss_fit_particles(ft, particles)  --- this fits all the spots found in FindSpots and writes the parameters to particle_coefs
//         You should go through and make sure all the particles look good. Take a look at res_sum for each of the fits (output from gauss_fit_particles),
//         and fix any particles that contain multiple g-vectors (if you think you can).
//         The best way to "fix" them so far (imo) is to run ShowXYZLineProfile(ft, particles, r, 0) where "r" is the spot you want to fix. This creates three graphs, one for each x,y,z direction.
//         You can take the gaussian sigmas (widths) and manually type them into the W_coef wave. Then regenerate the particle fit with these sigmas (modify them if you don't like the
//         x,y, or z fit widths by the way) by running regenerateSpotParameters(ft, r, particles, particle_coefs, W_coef).
//         This "fixing" is attempted automatically in a function called: try_automagic_spot_fix(ft, particles, particle_coefs, r, cut)  (set cut == 0 first). However, be careful using this, I don't know if or how well it will even work.
//         To undo the "fixing" run:  gauss_fit_particle(ft, particles, r, particle_coefs, 0, 0)   and   create_output_particle(ft, particles, particle_coefs, r, out, out_fortran)
//	     You can also try to 
//     Hopefully now you have good fits for all your particles.

//     Run create_output_wave(ft, particles, particle_coefs)  --- this generates "out" which is a modified version of "particles" for looking at
//         It also creates "out_fortran" which is the out-wave you should pass to SaveParamfile for generating the 256 pixel IFT input parameters
//         If you just want to regenerate the "out" wave for particle "r", then use create_output_particle(ft, particles, particle_coefs, r)
//             This regenerates "out_fortran" as well.

//     Change the spot name in the function SaveParamfile
//     Run SaveParamfile(out_fortran)  --- this is the paramfile you will use for the IFT and post-IFT batch analysis
//     Change the modelfile in the paramfile you saved to the correct one

//     Copy to ACI, run IFT by passing paramfile to slurm.sh (you may need to modify slurm.sh first, so take a look)

//     Another great function is: AppendGvecsToVk(ft, out, kexp, vkexp) after you have done the above. In order to use this function you have
//         to import the V(k) experimental data. I call the k wave "kexp" and the vk wave "vkexp". This plots the g-vectors you found above on
//         top of V(k) along with error bars indicating the size of the spot. Now you can visually see where all your g-vecors are!
//         You may want to change the function to only plot the top e.g. 20 spots (two lines, a make and a for loop).

// FindSpots(ft,15)
// gauss_fit_particles(ft, particles)
// There is currently some bug so run this next:  gauss_fit_particle(ft, particles, 0, particle_coefs, 0, 0)
// create_output_wave(ft, particles, particle_coefs)
// search_for_multiple_spots(ft, particles, particle_coefs, out, out_fortran) // searches for potential window'd images that may have multiple g-vectors in them. Be careful, this erases any hand made modifications you have made!
// create_output_wave(ft, particles, particle_coefs)

// window_edge_checking(out)
// cut_ft_smaller_for_isosurface(ft)
// create_windowed_spots(ft, out, rr)
// create_unwindowed_spots(ft, out, rr)
// AppendGvecsToVk(ft, out, kexp, vkexp)

function IsoAverage3D(dat, strip_width)
	wave dat
	variable strip_width
	
	variable x0, y0, z0, dx, dy, dz
	x0 = DimOffSet(dat, 0)
	y0 = DimOffSet(dat, 1)
	z0 = DimOffset(dat, 2)
	dx = DimDelta(dat, 0)
	dy = DimDelta(dat, 1)
	dz = DimDelta(dat, 2)
	
	// assume dat is centered at (0,0,0) and extends from (-x0 to +x0), etc.
	Make/O/N=(min( min(DimSize(dat, 0), DimSize(dat, 1)), Dimsize(dat, 2)) / (2*strip_width) +1) iso_av, npix
	SetScale/I x, 0, min( min(abs(x0), abs(y0)), abs(z0)), iso_av, npix
	iso_av = 0
	npix = 0
	
	variable i, j, k, r
	for(i=0; i<DimSize(dat, 0); i+=1)
		for(j=0; j<DimSize(dat, 1); j+=1)
			for(k=0; k<DimSize(dat, 2); k+=1)
				r = sqrt( (x0+dx*i)^2 + (y0+dy*j)^2 + (z0+dz*k)^2 )
				if(r < rightx(iso_av))
					npix[x2pnt(npix, r)] += 1
					iso_av[x2pnt(iso_av, r)] += dat[i][j][k]
				endif
			endfor
		endfor
	endfor	
		
	iso_av /= npix	
	Killwaves npix
end

function RadiusChop(dat, rmin, rmax)
	wave dat
	variable rmin, rmax
	
	Duplicate/O dat radius_chop
	radius_chop = 0

	variable x0, y0, z0, dx, dy, dz
	x0 = DimOffSet(dat, 0)
	y0 = DimOffSet(dat, 1)
	z0 = DimOffset(dat, 2)
	dx = DimDelta(dat, 0)
	dy = DimDelta(dat, 1)
	dz = DimDelta(dat, 2)
	
	variable i, j, k, r
	for(i=0; i<DimSize(dat, 0); i+=1)
		for(j=0; j<DimSize(dat, 1); j+=1)
			for(k=0; k<DimSize(dat, 2); k+=1)
				r = sqrt( (x0+dx*i)^2 + (y0+dy*j)^2 + (z0+dz*k)^2 )
				if( (r <= rmax) && (r >= rmin) )
					radius_chop[i][j][k] = dat[i][j][k]
				endif
			endfor
		endfor
	endfor	
end

Function makeSphere(pointsx,pointsy)
	Variable pointsx,pointsy
	
	Variable i,j,rad
	Make/O/n=(pointsx,pointsy,3) sphereData
	Variable anglePhi,angleTheta
	Variable dPhi,dTheta
	
	
	dPhi=2*pi/(pointsx-1)
	dTheta=pi/(pointsy-1)
	Variable xx,yy,zz
	Variable sig
	
	for(j=0;j<pointsy;j+=1)
		angleTheta=j*dTheta
		zz=sin(angleTheta)
		if(angleTheta>pi/2)
			sig=-1
		else
			sig=1
		endif
		for(i=0;i<pointsx;i+=1)
			anglePhi=i*dPhi
			xx=zz*cos(anglePhi)
			yy=zz*sin(anglePhi)
			sphereData[i][j][0]=xx
			sphereData[i][j][1]=yy
			sphereData[i][j][2]=sig*sqrt(1-xx*xx-yy*yy)
		endfor
	endfor
End

Function makeSphereColorWave(pointsx,pointsy,dat,radius,thick)
	Variable pointsx,pointsy
	wave dat
	variable radius, thick
	
	Variable i,j,rad
	Make/O/N=(pointsx, pointsy) sphereGray
	SetScale/I x 0, 2*Pi, "phi", sphereGray
	SetScale/I y 0, Pi, "theta", sphereGray
	Make/O/n=(pointsx,pointsy,4) sphereColor
	Variable anglePhi,angleTheta
	Variable dPhi,dTheta
	string n
	
	dPhi=2*pi/(pointsx-1)
	dTheta=pi/(pointsy-1)
	Variable xx,yy,zz
	Variable sig
	
	for(j=0;j<pointsy;j+=1)
		angleTheta=j*dTheta
		zz=sin(angleTheta)
		if(angleTheta>pi/2)
			sig=-1
		else
			sig=1
		endif
		for(i=0;i<pointsx;i+=1)
			anglePhi=i*dPhi
			
			xx = radius*sin(angleTheta)*cos(anglePhi)
			yy = radius*sin(angleTheta)*sin(anglePhi)
			zz = radius*cos(angleTheta)
			sphereGray[i][j] = VolumeSum(dat, thick, xx, yy, zz)
		endfor
	endfor

	n = ""
	n = ReplaceNumberByKey("radius", n, radius, "=")
	n = ReplaceNumberByKey("thick",  n, thick, "=")
	Note/k sphereGray, n

	wavestats/Q sphereGray
	sphereColor[][][0,2] = (sphereGray[p][q] - V_min) / (V_max - V_min)
	sphereColor[][][3] = 1

End

function VolumeStdDev(dat,radius, xc, yc, zc, temp)
	wave dat, temp
	variable radius, xc, yc, zc
	variable s, ii, jj, kk //,rp, xp, yp, zp, dist,
	//Make/O/N=44395 temp
	//Make/O/N=91125 temp
	
	//xp = round( (xc - DimOffset(dat, 0)) / DimDelta(dat, 0))
	//yp =  round( (yc - DimOffset(dat, 1)) / DimDelta(dat, 1))
	//zp =  round( (zc - DimOffset(dat, 2)) / DimDelta(dat, 2))
	//rp = ceil(radius / DimDelta(dat, 0))
	//rp = radius/2
	//print xp,xc,yp,yc,zp,zc
	//Make/N=((2*radius+1)^3) temp
	//if(rp<=1)
	//	printf "Radius too small.  oops.\r"
	//endif
	
	s = 0
	for(ii=xc-radius; ii<=xc+radius; ii+=1)
		for(jj=yc-radius; jj<=yc+radius; jj+=1)
			for(kk=zc-radius; kk<=zc+radius; kk+=1)
				//dist = (ii-xp)^2 + (jj-yp)^2 + (kk-zp)^2
				//if(dist < rp^2)
					//insertPoints 0,1,temp
					temp[s]= dat[ii][jj][kk]
					s += 1
				//endif
			endfor
		endfor
	endfor
	
	wavestats/Q temp
	//print numpnts(temp)
	//Killwaves temp
	//print V_sdev
	return V_sdev
end

function StdDev3d(dat, radius)
	wave dat
	variable radius
	variable ii, jj, kk, xdim, ydim, zdim
	Duplicate/O dat, stddevdat
	
	Make/O/N=((2*radius+1)^3) temp
	
	xdim = DimSize(dat,0)
	ydim = DimSize(dat,1)
	zdim = DimSize(dat,2)
	printf "%g,%g,%g\r",xdim,ydim,zdim
	for(ii=0; ii<xdim; ii+=1)
		for(jj=0; jj<ydim; jj+=1)
			for(kk=0; kk<zdim; kk+=1)
				stddevdat[ii][jj][kk] = VolumeStdDev(dat,radius,ii,jj,kk,temp)
				//printf "%g,%g,%g\r", ii,jj,kk
			endfor
			printf "%g,%g,%g\r", ii,jj,kk
		endfor
	endfor
end
	
function VolumeSum(dat, radius, xc, yc, zc)
	wave dat
	variable radius, xc, yc, zc
	
	variable rp, xp, yp, zp, dist, s, ii, jj, kk
	
	xp = round( (xc - DimOffset(dat, 0)) / DimDelta(dat, 0))
	yp =  round( (yc - DimOffset(dat, 1)) / DimDelta(dat, 1))
	zp =  round( (zc - DimOffset(dat, 2)) / DimDelta(dat, 2))
	rp = ceil(radius / DimDelta(dat, 0))
	
	if(rp<=1)
		printf "Radius too small.  oops.\r"
	endif
	
	//printf "Averaging %d pixels around position (%d, %d, %d).\r", rp, xp, yp, zp
	
	s=0
	for(ii=xp-rp; ii<=xp+rp; ii+=1)
		for(jj=yp-rp; jj<=yp+rp; jj+=1)
			for(kk=zp-rp; kk<=zp+rp; kk+=1)
				dist = (ii-xp)^2 + (jj-yp)^2 + (kk-zp)^2
				if(dist < rp^2)
					s += dat[ii][jj][kk]
				endif
			endfor
		endfor
	endfor
	
	return s
	
end

function Sphere2XYZ(phi, theta, sphereGray)
	variable phi, theta
	wave sphereGray
	
	variable radius
	radius = NumberByKey("radius", note(sphereGray),  "=")
	if(numtype(radius) == 2)
		printf "Cannot find radius of sphere-surface wave.  Exiting.\r"
		return 0
	endif
	printf "(x, y, z) = (%g, %g, %g)\r", radius*sin(theta)*cos(phi), radius*sin(theta)*sin(phi), radius*cos(theta)
	
end

function FindSpots(dat, numstdevs)
// This finds spots (or g-vectors) in the FT wave "dat"
// numstdevs is the number of stdevs above the mean
// of the FT *after* the center has been cut out.
// 15 seems to be a good number. This is used in
// ImageAnalyzeParticles as a threshold to look
// for spots.
	wave dat
	variable numstdevs
	variable rmin
	variable above
	
	printf "Processing IsoAverage3D...\r"
	setscale x,-1.5,1.5,dat
	setscale y,-1.5,1.5,dat
	setscale z,-1.5,1.5,dat
	IsoAverage3D(dat,2)
	Display iso_av
	rmin = FindFirstMin(iso_av)
	
	duplicate/o dat dat2

	variable x0, y0, z0, dx, dy, dz
	x0 = DimOffSet(dat, 0)
	y0 = DimOffSet(dat, 1)
	z0 = DimOffset(dat, 2)
	dx = DimDelta(dat, 0)
	dy = DimDelta(dat, 1)
	dz = DimDelta(dat, 2)
	
	printf "Doing radius chop on input data to cut out middle using cutoff %g...\r", rmin
	variable i, j, k, rr
	for(i=0; i<DimSize(dat, 0); i+=1)
		for(j=0; j<DimSize(dat, 1); j+=1)
			for(k=0; k<DimSize(dat, 2); k+=1)
				rr = sqrt( (x0+dx*i)^2 + (y0+dy*j)^2 + (z0+dz*k)^2 )
				if( ( rr < rmin) )
					dat2[i][j][k] = 0
				endif
			endfor
		endfor
	endfor
	WaveStats/Q dat2
	above = V_avg + numstdevs*V_sdev
	
	printf "Setting up spots for ImageAnalyzeParticles using cutoff=%g...\r", above
	ImageThreshold/I/T=(above)/M=0 dat2
	duplicate/o M_ImageThresh, spots
	Killwaves M_ImageThresh
	
	printf "Processing ImageAnalyzeParticles on spots...\r"
	ImageAnalyzeParticles/A=4 stats spots
	wave M_3DParticleInfo
	// Save the particles in "deleteme".
	// This wave will be deleted after the duplicate spots have been
	// removed, which is done after sorting based on max intensity.
	// If two spots have very very similar intensty (within the deviation
	// in symmetry in the FT) then things could go wrong here (i.e. 
	// you might miss a spot and/or get two of the same one).
	Make/O/N=(dimsize(M_3DParticleInfo,0),dimsize(M_3DParticleInfo,1)+1) deleteme
	for( i=0; i<dimsize(M_3DParticleInfo,0); i+=1)
		for( j=0; j<dimsize(M_3DParticleInfo,1); j+=1)
			deleteme[i][j] = M_3DParticleInfo[i][j]
		endfor
	endfor
	variable minRow, maxRow, minCol, maxCol, minLayer, maxLayer
	variable maxint, sumint
	printf "Finding intensity of spots...\r"
	for(rr=0; rr<dimsize(deleteme,0); rr+=1)
		deleteme[rr][6] = 0.0
		maxint = 0
		sumint = 0.0
		for( i=deleteme[rr][0]; i<=deleteme[rr][1]; i+=1)
			for( j=deleteme[rr][2]; j<=deleteme[rr][3]; j+=1)
				for( k=deleteme[rr][4]; k<=deleteme[rr][5]; k+=1)
					sumint += dat[i][j][k]
					if(maxint < dat[i][j][k])
						maxint = dat[i][j][k]
						deleteme[rr][6] = maxint
						// Save the center as well
						deleteme[rr][7] = i
						deleteme[rr][8] = j
						deleteme[rr][9] = k
					endif
				endfor
			endfor
		endfor
		deleteme[rr][11] = sumint
	endfor
	// I tried to normalize below by the volume of the spot, 
	// but if we do that it doesn't help me to understand what's going on any better.
	// I basically just get a stdev around the background it seems like.
	//deleteme[][11] /= ( (deleteme[p][1]-deleteme[p][0]+1) * (deleteme[p][3]-deleteme[p][2]+1) * (deleteme[p][5]-deleteme[p][4]+1) )
	printf "Sorting spots by max intensity\r"
	MDsort(deleteme,6,1)
	// Set gvecs
	for(rr=0; rr<dimsize(deleteme,0); rr+=1)
		deleteme[rr][10] = sqrt( (dimsize(dat,0)/2-deleteme[rr][7])^2 + (dimsize(dat,1)/2-deleteme[rr][8])^2 + (dimsize(dat,2)/2-deleteme[rr][9])^2 )*dimdelta(dat,0)
	endfor
	
	printf "Calculating background...\r"
	Duplicate/O dat2, background
	for( rr=0; rr<dimsize(deleteme,0); rr+=1)
		for( i=deleteme[rr][0]; i<=deleteme[rr][1]; i+=1)
			for( j=deleteme[rr][2]; j<=deleteme[rr][3]; j+=1)
				for( k=deleteme[rr][4]; k<=deleteme[rr][5]; k+=1)
					background[i][j][k] = 0.0
				endfor
			endfor
		endfor
	endfor
	printf "Mean of background: %g\r", mean(background)
	
	printf "Copying into particles...\r"
	// Now that "deleteme" is done being changed, copy the particles
	// we want (i.e. remove duplicates) into "particles"
	Make/D/O/N=(dimsize(deleteme,0)/2,dimsize(deleteme,1)) particles
	for( i=0; i<dimsize(deleteme,0); i+=2)
		for( j=0; j<dimsize(deleteme,1); j+=1)
			particles[i/2][j] = deleteme[i][j]
		endfor
	endfor
	
	make/o/n=(dimsize(particles,0)) spot_intensities
	spot_intensities = particles[p][6]
	make/o/n=(dimsize(particles,0)) gvecs_unsorted
	gvecs_unsorted = particles[p][10]
	Display  spot_intensities vs gvecs_unsorted
	ModifyGraph mode=3,marker=2,rgb=(0,0,0)
	Display  spot_intensities
	ModifyGraph mode=3,marker=2,rgb=(0,0,0)
	
	killwaves deleteme
	killwaves dat2, background, spots // Comment this line if you want to know more about what's going on.
	return 0
end

function gauss_fit_particle(dat, particles, r, particle_coefs, hold_sigmas, hold_centers)
	// set hold_sigmas to 0 or 1
	// hold_sigmas == 0 means to not hold the sigmas in the fit
	// hold_sigmas == 1 means to hold the sigmas in the fit
	// use hold_sigmas == 1 when you are trying to fix a spot's fitted sigmas manually
	
	// set hold_centers to 0 or 1
	// hold_centers == 0 means to not hold the center positions in the fit
	// hold_centers == 1 means to hold the center positions in the fit
	// use hold_centers == 1 when you are trying to fix a spot's fitted center positions manually
	wave dat, particles, particle_coefs
	variable r, hold_sigmas, hold_centers
	variable i, j, k
	variable xwn, ywn, zwn
	wave W_coef,image,image_fit
	variable res_sum
	variable best_i, best_res_sum
	best_res_sum = 99999999.0
	for( i=1; i<=3; i+=1)
		//printf "Running 3D gauss fit on particle %g with extra=%g\r",r,i
		res_sum = fitSpotTo3DGauss(dat, particles, r, i, hold_sigmas, hold_centers) //0 means do not hold sigmas
		xwn = round(sqrt(2*ln(10))*W_coef[0])
		ywn = round(sqrt(2*ln(10))*W_coef[1])
		zwn = round(sqrt(2*ln(10))*W_coef[2])
		if(abs(res_sum) < abs(best_res_sum) && xwn > 2 && ywn > 2 && zwn >2)
			best_res_sum = res_sum
			best_i = i
		endif
	endfor
	if( best_i != 3)
		res_sum = fitSpotTo3DGauss(dat, particles, r, best_i, hold_sigmas, hold_centers)
	endif
	xwn = round(sqrt(2*ln(10))*W_coef[0])
	ywn = round(sqrt(2*ln(10))*W_coef[1])
	zwn = round(sqrt(2*ln(10))*W_coef[2])
	variable V_chisq
	printf "Particle %g's best extra fit = %g with a  res_sum = %g and chisq = %g\r", r, best_i, best_res_sum, V_chisq
	//printf "Widths: %g %g %g\r", xwn, ywn, zwn
	//print W_coef
	make/o/n=(xwn*2,ywn*2) image_layer
	image_layer = image[p][q][round(W_coef[8])]
	//NewImage  root:image_layer
	make/o/n=(xwn*2,ywn*2) image_layer_fit
	image_layer_fit = image_fit[p][q][round(W_coef[8])]
	//AppendMatrixContour image_layer_fit
	
	// Save the correct fit parameters for each particle
	particle_coefs[r][0] = W_coef[0]
	particle_coefs[r][1] = W_coef[1]
	particle_coefs[r][2] = W_coef[2]
	particle_coefs[r][3] = W_coef[3]
	particle_coefs[r][4] = W_coef[4]
	particle_coefs[r][5] = W_coef[5]
	particle_coefs[r][6] = W_coef[6]
	particle_coefs[r][7] = W_coef[7]
	particle_coefs[r][8] = W_coef[8]
	// Also save the position of the max intensity in this spot so we know
	// where it lies in relation to the rest of the FT
	variable found = 0
	for( i=0; i<=dimsize(image,0); i+=1)
		for( j=0; j<=dimsize(image,1); j+=1)
			for( k=0; k<=dimsize(image,2); k+=1)
				if( 1 == image[i][j][k] )
					particle_coefs[r][9]   = i   // x center
					particle_coefs[r][10] = j   // y center
					particle_coefs[r][11] = k  // z center
					found = 1
				endif
			endfor
		endfor
	endfor
	if( !found)
		printf "ERROR! COULD NOT FIND CORRECT INTENSITY IN PARTICLE %g!\r", r
	endif
end

function gauss_fit_particles(dat, particles)
	wave dat, particles
	variable r, i, j, k
	variable xwn, ywn, zwn
	wave W_coef,image,image_fit
	variable res_sum
	variable best_i, best_res_sum = 99999999.0
	make/o/n=(dimsize(particles,0), 12) particle_coefs
	
	for( r=0; r<dimsize(particles,0); r+=1)
		gauss_fit_particle(dat, particles, r, particle_coefs, 0, 0)
	endfor
end

function regenerateSpotParameters(dat, r, particles, particle_coefs, W_coef)
	// W_coef is used in gauss_fit_particle in this case
	variable r
	wave dat, particle_coefs, W_coef, particles
	gauss_fit_particle(dat, particles, r, particle_coefs, 1, 0) // set hold == 1 so that the sigmas are held constant, but don't hold centers constant
end

function create_output_wave(dat, particles, particle_coefs)
	wave dat, particles, particle_coefs

	wave image
	wave W_coef = $"W_coef"
	variable r, i, j, k
	make/o/n=(dimsize(particles,0),dimsize(particles,1)+9) out
	duplicate/o out, out_fortran
	variable x_start_pix, y_start_pix, z_start_pix, found, x_start_pix_part, y_start_pix_part, z_start_pix_part
	variable offset_x, offset_y, offset_z
	for( r=0; r<dimsize(particles,0); r+=1)
		create_output_particle(dat, particles, particle_coefs, r, out, out_fortran)
	endfor
end

function create_output_particle(dat, particles, particle_coefs, r, out, out_fortran)
	wave dat, particles, particle_coefs, out, out_fortran
	variable r
	wave W_coef = $"W_coef"
	variable i, j, k
	//make/o/n=(dimsize(particles,0),dimsize(particles,1)+9) out
	variable x_start_pix, y_start_pix, z_start_pix, found, x_start_pix_part, y_start_pix_part, z_start_pix_part
	variable offset_x, offset_y, offset_z
	// Locate max intensity in particle so I can match it to the position in dat
	found = 0
	for( i=particles[r][0]; i<=particles[r][1]; i+=1)
		for( j=particles[r][2]; j<=particles[r][3]; j+=1)
			for( k=particles[r][4]; k<=particles[r][5]; k+=1)
				if( particles[r][6] == dat[i][j][k] )
					x_start_pix = i
					y_start_pix = j
					z_start_pix = k
					found = 1
				endif
			endfor
		endfor
	endfor
	if( !found)
		printf "ERROR! COULD NOT FIND CORRECT INTENSITY IN DATA %g!\r", r
	endif
	//print x_start_pix, y_start_pix, z_start_pix
	x_start_pix_part = particle_coefs[r][9]
	y_start_pix_part = particle_coefs[r][10]
	z_start_pix_part = particle_coefs[r][11]
	//print x_start_pix_part, y_start_pix_part, z_start_pix_part
	
	// Set the necessary parameters in an output wave which I will later save into a parameter file for Fortran
	// to read in.
	// First, increase the sigmas of the fit by 1.5 or 2 or something so that the fit becomes
	// a windowing function instead. What you multiply by determines how big the final spot is.
	W_coef[] = particle_coefs[r][p]
	//print W_coef
	W_coef[0] *= 1.3
	W_coef[1] *= 1.3
	W_coef[2] *= 1.3
	/// Now we know where the maximum intensity of the spot is in both the data and in the particle (i.e. "image")
	// So now we need to correct for where the gaussian fit found the center
	// which is stored in W_coef [6] to [8]
	offset_x = x_start_pix_part - round(W_coef[6])
	offset_y = y_start_pix_part - round(W_coef[7])
	offset_z = z_start_pix_part - round(W_coef[8])
	out[r][0] = offset_x + x_start_pix - x_start_pix_part + round(W_coef[6] - sqrt(2*ln(10))*W_coef[0])
	out[r][1] = offset_x + x_start_pix - x_start_pix_part + round(W_coef[6] + sqrt(2*ln(10))*W_coef[0])
	out[r][2] = offset_y + y_start_pix - y_start_pix_part + round(W_coef[7] - sqrt(2*ln(10))*W_coef[1])
	out[r][3] = offset_y + y_start_pix - y_start_pix_part + round(W_coef[7] + sqrt(2*ln(10))*W_coef[1])
	out[r][4] = offset_z + z_start_pix - z_start_pix_part + round(W_coef[8] - sqrt(2*ln(10))*W_coef[2])
	out[r][5] = offset_z + z_start_pix - z_start_pix_part + round(W_coef[8] + sqrt(2*ln(10))*W_coef[2])
	out[r][6] = particles[r][6]
	out[r][7] = x_start_pix - x_start_pix_part + round(W_coef[6])
	out[r][8] = y_start_pix - y_start_pix_part + round(W_coef[7])
	out[r][9] = z_start_pix - z_start_pix_part + round(W_coef[8])
	out[r][10] = particles[r][10]
	
	// Set the gaussian fit parameters for the fitting window so Fortran knows what they are.
	out[r][11] = W_coef[0]
	out[r][12] = W_coef[1]
	out[r][13] = W_coef[2]
	out[r][14] = W_coef[3]
	out[r][15] = W_coef[4]
	out[r][16] = W_coef[5]
	out[r][17] = x_start_pix - x_start_pix_part + W_coef[6]
	out[r][18] = y_start_pix - y_start_pix_part + W_coef[7]
	out[r][19] = z_start_pix - z_start_pix_part + W_coef[8]
	
	if(0 == 1)
	// Create the final particle to look at.
	// It is necessary that only 1 g-vector is included in this image!
	// If not, you will probably have to manually correct. Good luck.
	// This is the final image with the Gaussian Window applied.
	// It takes a while to generate, be a bit patient.
	// You really should look at this for every particle you do an IFT on.
	// To do this for individual particles, you have to modify the commented
	// for loop above to only go through a single r.
	print "Creating final window'd particle...\r"
	duplicate/o dat, dat_dup
	dat_dup *= gauss3d(p,q,r,W_coef)
	duplicate/o/R=[out[r][0],out[r][1]][out[r][2],out[r][3]][out[r][4],out[r][5]] dat_dup,part
	killwaves dat_dup
	// It's easiest to just look at the center layer. The issue is that
	// out[r][9] is the layer in the full FT, so you have to correct for
	// the z starting position of the particle (which is out[r][4]).
	make/o/n=(out[r][1]-out[r][0]+1,out[r][3]-out[r][2]+1) part_layer
	part_layer = part[p][q][out[r][9]-out[r][4]]
	endif
	
	if( 1 == 1)
	for( i=0; i<dimsize(out,1); i+=1)
		out_fortran[r][i] = out[r][i]
	endfor
	// Modify for generating paramfile if you are doing the final run.
	// The issue is that we run the IFT with 256 pixels and not 512,
	// so we have to 1/2 the x,y,z coordinates and the sigmas.
	out_fortran[r][0] = round((offset_x + x_start_pix - x_start_pix_part + W_coef[6] - sqrt(2*ln(10))*W_coef[0])/2)+1
	out_fortran[r][1] = round((offset_x + x_start_pix-x_start_pix_part+W_coef[6] + sqrt(2*ln(10))*W_coef[0])/2)+1
	out_fortran[r][2] = round((offset_y + y_start_pix-y_start_pix_part+W_coef[7] - sqrt(2*ln(10))*W_coef[1])/2)+1
	out_fortran[r][3] = round((offset_y + y_start_pix-y_start_pix_part+W_coef[7] + sqrt(2*ln(10))*W_coef[1])/2)+1
	out_fortran[r][4] = round((offset_z + z_start_pix-z_start_pix_part+W_coef[8] - sqrt(2*ln(10))*W_coef[2])/2)+1
	out_fortran[r][5] = round((offset_z + z_start_pix-z_start_pix_part+W_coef[8] + sqrt(2*ln(10))*W_coef[2])/2)+1
	out_fortran[r][7] = round((x_start_pix - x_start_pix_part + W_coef[6]))//2)
	out_fortran[r][8] = round((y_start_pix - y_start_pix_part + W_coef[7]))//2)
	out_fortran[r][9] = round((z_start_pix - z_start_pix_part + W_coef[8]))//2)
	out_fortran[r][11] = W_coef[0]/2 //sx
	out_fortran[r][12] = W_coef[1]/2
	out_fortran[r][13] = W_coef[2]/2
	out_fortran[r][14] = W_coef[3] //cxy
	out_fortran[r][15] = W_coef[4]
	out_fortran[r][16] = W_coef[5]
	out_fortran[r][17] = (x_start_pix-x_start_pix_part+W_coef[6])/2+1 //x0
	out_fortran[r][18] = (y_start_pix-y_start_pix_part+W_coef[7])/2+1
	out_fortran[r][19] = (z_start_pix-z_start_pix_part+W_coef[8])/2+1
	else
	// this changes W_coef to be the windowing function instead of the fit
	W_coef[6] = out[r][17]
	W_coef[7] = out[r][18]
	W_coef[8] = out[r][19]
	endif
end

function fitSpotTo3DGauss(dat, particles, r, extra, hold_sigmas, hold_centers)
	// If you use either holds, the wave W_coef must be set for sigmas or centers outside this function
	// Does a fit on particle r
	// extra specifies the extra amount you multiply the initial fit window by to get a better window (it multiplies)
	// Updated variables include:  temp, temp_fit, residuals, and gauss_fit
	// image_layer and image_layer_fit are NOT updated here (they are updated/made in the iterating gauss fit function over all particles)
	wave dat, particles
	variable r, extra, hold_sigmas, hold_centers
	variable sx,sy,sz
	
	// Set up first fit based on the inputs
	variable xw, yw, zw
	variable print_true = 0
	variable i,j,k
	if(hold_sigmas == 0 && hold_centers == 0) // no holding
		Make/D/N=9/O W_coef
	else
		wave W_coef = $"W_coef"
	endif
	variable x0,y0,z0
	//if(hold_centers != 0)  // holding centers
	//	//wave W_coef = $"W_coef"
	//	x0 = W_coef[6]
	//	y0 = W_coef[7]
	//	z0 = W_coef[8]
	//endif
	xw = particles[r][1]-particles[r][0] // right edge - left edge = width
	yw = particles[r][3]-particles[r][2]
	zw = particles[r][5]-particles[r][4]
	x0 = particles[r][7]-particles[r][0] - (xw/2 - xw*extra/2) // center - left edge + extra addon width = new center
	y0 = particles[r][8]-particles[r][2] - (yw/2 - yw*extra/2)
	z0 = particles[r][9]-particles[r][4] - (zw/2 - zw*extra/2)
	variable xcl, xcr, ycl, ycr, zcl, zcr
	xcl = x0 - (particles[r][7] - particles[r][0])/2  // the center should stay within the original spot region found by ImageAnalyzeParticles
	xcr = x0 + (particles[r][7] - particles[r][0])/2
	ycl = y0 - (particles[r][8] - particles[r][2])/2
	ycr = y0 + (particles[r][8] - particles[r][2])/2
	zcl = z0 - (particles[r][9] - particles[r][4])/2  // the center should stay within the original spot region found by ImageAnalyzeParticles
	zcr = z0 + (particles[r][9] - particles[r][4])/2
	//printf( " %g %g %g\r" ) (xw/2 - xw*extra/2), (yw/2 - yw*extra/2), (zw/2 - zw*extra/2)
	//printf( " %g %g %g\r" ) x0, y0, z0
	//printf( " %g %g %g %g %g %g\r\r" ) xcl, xcr, ycl, ycr, zcl, zcr
	if(hold_sigmas == 0) // no holding
		//Make/D/N=9/O W_coef
		W_coef[0] = {xw/2, yw/2, zw/2, 0.1, 0.1, 0.1, xw*extra/2, yw*extra/2, zw*extra/2}
		Make/O/T/N=18 T_Constraints
		T_Constraints = {"K0 > 0","K0 < 32","K1 > 0","K1 < 30","K2 > 0","K2 < 24","K3 > -.5","K3 < .5","K4 > -.5","K4 < .5","K5 > -.5","K5 < .5","K6 > 0","K6 < 32","K7 > 0","K7 < 30","K8 > 0","K8 < 24"}
		T_Constraints[1]   = "K0 < " + num2str(xw*extra)
		T_Constraints[3]   = "K1 < " + num2str(yw*extra)
		T_Constraints[5]   = "K2 < " + num2str(zw*extra)
		//T_Constraints[13] = "K6 < " + num2str(xw*extra)
		//T_Constraints[15] = "K7 < " + num2str(yw*extra)
		//T_Constraints[17] = "K8 < " + num2str(zw*extra)
		T_Constraints[12] = "K6 > " + num2str( xcl ) // centers, which should not deviate from the original spot region found in ImageAnalyzeParticles
		T_Constraints[13] = "K6 < " + num2str( xcr )
		T_Constraints[14] = "K7 > " + num2str( ycl )
		T_Constraints[15] = "K7 < " + num2str( ycr )
		T_Constraints[16] = "K8 > " + num2str( zcl )
		T_Constraints[17] = "K8 < " + num2str( zcr )
	else // holding sigmas
		//wave W_coef = $"W_coef"
		sx = W_coef[0]
		sy = W_coef[1]
		sz = W_coef[2]
		W_coef[0] = {xw/2, yw/2, zw/2, 0.1, 0.1, 0.1, xw*extra/2, yw*extra/2, zw*extra/2}
		W_coef[0] = sx
		W_coef[1] = sy
		W_coef[2] = sz
		printf "Using sigmas for fit: %g %g %g\r", W_coef[0], W_coef[1], W_coef[2]
		Make/O/T/N=12 T_Constraints
		T_Constraints = {"K3 > -.5","K3 < .5","K4 > -.5","K4 < .5","K5 > -.5","K5 < .5","K6 > 0","K6 < 32","K7 > 0","K7 < 30","K8 > 0","K8 < 24"}
		//T_Constraints[7] = "K6 < " + num2str(xw*extra)
		//T_Constraints[9] = "K7 < " + num2str(yw*extra)
		//T_Constraints[11] = "K8 < " + num2str(zw*extra)
		T_Constraints[6] = "K6 > " + num2str( xcl ) // centers, which should not deviate from the max intensity point by more than 1/4 of the window size
		T_Constraints[7] = "K6 < " + num2str( xcr )
		T_Constraints[8] = "K7 > " + num2str( ycl )
		T_Constraints[9] = "K7 < " + num2str( ycr )
		T_Constraints[10] = "K8 > " + num2str( zcl )
		T_Constraints[11] = "K8 < " + num2str( zcr )
	endif
	
	if(hold_centers != 0)  // holding centers
		//wave W_coef = $"W_coef"
		W_coef[6] = x0
		W_coef[7] = y0
		W_coef[8] = z0
		printf "Using center positions from particles for fit: %g %g %g\r", W_coef[6], W_coef[7], W_coef[8]
		Make/O/T/N=12 T_Constraints
		T_Constraints = {"K0 > 0","K0 < 32","K1 > 0","K1 < 30","K2 > 0","K2 < 24","K3 > -.5","K3 < .5","K4 > -.5","K4 < .5","K5 > -.5","K5 < .5"}
		T_Constraints[1]   = "K0 < " + num2str(xw*extra)
		T_Constraints[3]   = "K1 < " + num2str(yw*extra)
		T_Constraints[5]   = "K2 < " + num2str(zw*extra)
	endif

	// Create the image to fit
	make/o/n=(xw*extra,yw*extra,zw*extra) image
	for( i=particles[r][7]-xw*extra/2; i<particles[r][7]+xw*extra/2; i +=1)
		for( j=particles[r][8]-yw*extra/2; j<particles[r][8]+yw*extra/2; j+=1)
			for( k=particles[r][9]-zw*extra/2; k<particles[r][9]+zw*extra/2; k+=1)
				image[i-particles[r][7]+xw*extra/2][j-particles[r][8]+yw*extra/2][k-particles[r][9]+zw*extra/2] = dat[i][j][k]
			endfor
		endfor
	endfor
	variable wm = particles[r][6]
	image /= wm // Need to normalize to 1 because there is no intensity in the gauss3df function.
	// Do the fit
	duplicate/o image,image_fit
	image_fit = NaN
	if(hold_sigmas == 0) // 0 means false, 1 means true
		FuncFitMD/Q/H="000000000"/NTHR=0 gauss3df W_coef  image /D=image_fit /C=T_Constraints
	elseif(hold_centers != 0)
		FuncFitMD/Q/H="000000111"/NTHR=0 gauss3df W_coef  image /D=image_fit /C=T_Constraints
	else
		FuncFitMD/Q/H="111000000"/NTHR=0 gauss3df W_coef  image /D=image_fit /C=T_Constraints
	endif
	
	// Resize based on outputs of sigmas from fit and then find the correct center position
	// for the new image by using the same fitting coefficients but allowing only the center
	// position of the gaussian fit to change.
	// NO! This works poorly because if there is another g-vector in the image that results 
	// in a better fit (e.g. if the other spot has a better shape) then the center position gets
	// moved. It is much better to calculate how much the resizing changed the center position
	// by and fix it manually.
	variable xwn, ywn, zwn
	xwn = round(sqrt(2*ln(10))*W_coef[0])
	ywn = round(sqrt(2*ln(10))*W_coef[1])
	zwn = round(sqrt(2*ln(10))*W_coef[2])
	make/o/n=(xwn*2,ywn*2,zwn*2) image
	image = 0
	for( i=particles[r][7]-xwn; i<particles[r][7]+xwn; i +=1)
		for( j=particles[r][8]-ywn; j<particles[r][8]+ywn; j+=1)
			for( k=particles[r][9]-zwn; k<particles[r][9]+zwn; k+=1)
				image[i-particles[r][7]+xwn][j-particles[r][8]+ywn][k-particles[r][9]+zwn] = dat[i][j][k]
			endfor
		endfor
	endfor
	image /= wm
	duplicate/o image,image_fit
	//image_fit = NaN
	//T_Constraints = {"K6 > 0","K6 < 32","K7 > 0","K7 < 30","K8 > 0","K8 < 24"}
	//T_Constraints[1] = "K6 < " + num2str(xw*extra)
	//T_Constraints[3] = "K7 < " + num2str(yw*extra)
	//T_Constraints[5] = "K8 < " + num2str(zw*extra)
	//FuncFitMD/Q/H="111111000"/NTHR=0 gauss3df W_coef  image /D=image_fit /C=T_Constraints
	//printf( "Reduced V_chisq = %g\r" ) V_chisq/ (xwn*ywn*zwn - 6)
	// Fix manually:
	//print W_coef
	//printf( "%g %g %g\r" ) particles[r][7]-xw*extra/2, particles[r][7]-xwn, xwn - round(xw*extra/2)
	//printf( "%g %g %g\r" ) particles[r][8]-yw*extra/2, particles[r][8]-ywn, ywn - round(yw*extra/2)
	//printf( "%g %g %g\r" ) particles[r][9]-zw*extra/2, particles[r][9]-zwn, zwn - round(zw*extra/2)
	W_coef[6] += xwn - round(xw*extra/2)
	W_coef[7] += ywn - round(yw*extra/2)
	W_coef[8] += zwn - round(zw*extra/2)
	//print W_coef
	//printf("\r")
	
	// Create the residual image and an image of just the gaussian fit, and calculate the total residual
	variable res_sum = 0.0
	duplicate/o image, residuals
	//residuals = 0
	for( i=0; i<xwn*2; i += 1)
		for( j=0; j<ywn*2; j+=1)
			for( k=0; k<zwn*2; k+=1)
				residuals[i][j][k] = image[i][j][k] - gauss3d(i,j,k,W_coef)
			endfor
		endfor
	endfor
	for( i=xwn/4; i<xwn*3/4; i += 1)
		for( j=ywn/4; j<ywn*3/4; j+=1)
			for( k=zwn/4; k<zwn*3/4; k+=1)
				res_sum += abs(residuals[i][j][k])
			endfor
		endfor
	endfor
	
	duplicate/o image, gauss_fit
	gauss_fit = 0
	for( i=0; i<xwn*2; i += 1)
		for( j=0; j<ywn*2; j+=1)
			for( k=1; k<zwn*2; k+=1)
				gauss_fit[i][j][k] = gauss3d(i,j,k,W_coef)
			endfor
		endfor
	endfor
	
	// Extract the center layer and show it on screen if print_true > 0
	// Show some other stuff too (optional based on if you comment it)
	if(print_true == 1)
		printf "Sum of residual: %g\r", res_sum//sum(residuals)
		printf "V_chisq=%g\r", V_chisq
		print W_coef
		make/o/n=(xwn*2,ywn*2) image_layer
		image_layer = image[p][q][round(W_coef[8])]
		make/o/n=(xwn*2,ywn*2) image_layer_fit
		image_layer_fit = image_fit[p][q][round(W_coef[8])]
		//NewImage  root:image_layer
		//AppendMatrixContour image_layer_fit
		//NewImage  root:residuals
		//WMAppend3DImageSlider();
		//NewImage  root:gauss_fit
		//WMAppend3DImageSlider();
	endif
	return res_sum
	//return V_chisq/ (xwn*ywn*zwn - 6)
end

function SaveParamfile(out_fortran)
	wave out_fortran
	// Particles wave should be output from FindSpots
	// columns 0-5 are the xmin, xmax, ymin, ymax, zmin, zmax
	// column 6 is the intensity of the brightest pixel, whose coordinates
	// are stored in columns 7-9. Column 10 is irrelevant.
	Make/O/N=(DimSize(out_fortran, 0)) xmin, ymin, zmin, xmax, ymax, zmax, xc, yc, zc, gvecs, x0, y0, z0, sx, sy, sz, cxy, cxz, cyz
	xmin = out_fortran[p][0]
	xmax = out_fortran[p][1]
	ymin = out_fortran[p][2]
	ymax = out_fortran[p][3]
	zmin = out_fortran[p][4]
	zmax = out_fortran[p][5]
	xc = out_fortran[p][7]
	yc = out_fortran[p][8]
	zc = out_fortran[p][9]
	gvecs = out_fortran[p][10]
	sx = out_fortran[p][11]
	sy = out_fortran[p][12]
	sz = out_fortran[p][13]
	cxy = out_fortran[p][14]
	cxz = out_fortran[p][15]
	cyz = out_fortran[p][16]
	x0 = out_fortran[p][17]
	y0 = out_fortran[p][18]
	z0 = out_fortran[p][19]
	
	variable i, f
	Open/T=".txt" f
	if(!strlen(S_filename))
		printf "Error opening file.\r"
		return -1
	endif
	
	fprintf f, "%s\n", "model-filename goes here"
	fprintf f, "%g\n", min(20,DimSize(out_fortran,0))
	for( i=0; i<DimSize(out_fortran,0); i+=1)   // Save all the particles, but only set Fortran to calculate the first 20
		fprintf f, "t3_gauss_spot%g_\n", i
		fprintf f, "%g %g %g\n", x0[i], y0[i], z0[i]
		fprintf f, "%g %g %g\n", sx[i], sy[i], sz[i]
		fprintf f, "%g %g %g\n", cxy[i], cxz[i], cyz[i]
		fprintf f, "%g %g %g %g\n", xmin[i]+1, xmax[i]+1, xc[i]+1, gvecs[i]
		fprintf f, "%g %g %g %g\n", ymin[i]+1, ymax[i]+1, yc[i]+1, gvecs[i]
		if(i != DimSize(particles,0)-1)
			fprintf f, "%g %g %g %g\n", zmin[i]+1, zmax[i]+1, zc[i]+1, gvecs[i]
		else
			fprintf f, "%g %g %g %g", zmin[i]+1, zmax[i]+1, zc[i]+1, gvecs[i]
		endif
	endfor
	close f
	killwaves xmin, xmax, ymin, ymax, zmin, zmax, xc, yc, zc
end

function create_windowed_spots(ft, out, rr)
	wave out, ft
	variable rr
	variable npix = dimsize(ft,0)
	Duplicate/O/R=[out[rr][0],out[rr][1]] [out[rr][2],out[rr][3]] [out[rr][4],out[rr][5]] ft, wind_im
	Duplicate/O/R=[npix-out[rr][1],npix-out[rr][0]] [npix-out[rr][3],npix-out[rr][2]] [npix-out[rr][5],npix-out[rr][4]] ft, wind_im_opp
	Make/O/N=9 coefs
	coefs[0] = out[rr][11] //sigmas
	coefs[1] = out[rr][12]
	coefs[2] = out[rr][13]
	coefs[3] = out[rr][14] //cxy's
	coefs[4] = out[rr][15]
	coefs[5] = out[rr][16]
	coefs[6] = out[rr][17] //x0's
	coefs[7] = out[rr][18]
	coefs[8] = out[rr][19]
	variable i, j, k
	for( i=0; i<dimsize(wind_im,0); i+=1)
		for( j=0; j<dimsize(wind_im,1); j+=1)
			for( k=0; k<dimsize(wind_im,2); k+=1)
				wind_im[i][j][k] *= gauss3d(i+out[rr][0],j+out[rr][2],k+out[rr][4],coefs)
			endfor
		endfor
	endfor
	coefs[6] = npix-coefs[6] // fix x0's
	coefs[7] = npix-coefs[7]
	coefs[8] = npix-coefs[8]
	coefs[3] *= -1 // fix cxy's
	coefs[4] *= -1
	coefs[5] *= -1
	for( i=0; i<dimsize(wind_im_opp,0); i+=1)
		for( j=0; j<dimsize(wind_im_opp,1); j+=1)
			for( k=0; k<dimsize(wind_im_opp,2); k+=1)
				wind_im_opp[i][j][k] *= gauss3d(i+(npix-out[rr][1]),j+(npix-out[rr][3]),k+(npix-out[rr][5]),coefs)
			endfor
		endfor
	endfor
	killwaves coefs
end

function create_unwindowed_spots(ft, out, rr)
	wave out, ft
	variable rr
	variable npix = dimsize(ft,0)
	Duplicate/O/R=[out[rr][0],out[rr][1]] [out[rr][2],out[rr][3]] [out[rr][4],out[rr][5]] ft, unwind_im
	Duplicate/O/R=[npix-out[rr][1],npix-out[rr][0]] [npix-out[rr][3],npix-out[rr][2]] [npix-out[rr][5],npix-out[rr][4]] ft, unwind_im_opp
end

function MDsort(w,keycol, reversed)
	Wave w
	variable keycol, reversed
	variable type
	type = Wavetype(w)
	make/Y=(type)/free/n=(dimsize(w,0)) key
	make/free/n=(dimsize(w,0)) valindex
	if(type == 0)
		Wave/t indirectSource = w
		Wave/t output = key
		output[] = indirectSource[p][keycol]
	else
		Wave indirectSource2 = w
		multithread key[] = indirectSource2[p][keycol]
 	endif
	valindex=p
 	if(reversed)
 		sort/a/r key,key,valindex
 	else
		sort/a key,key,valindex
 	endif
	if(type == 0)
		duplicate/free indirectSource, M_newtoInsert
		Wave/t output = M_newtoInsert
	 	output[][] = indirectSource[valindex[p]][q]
	 	indirectSource = output
	else
		duplicate/free indirectSource2, M_newtoInsert
	 	multithread M_newtoinsert[][] = indirectSource2[valindex[p]][q]
		multithread indirectSource2 = M_newtoinsert
 	endif 
end

function FindFirstMin(wave1d)
	wave wave1d
	variable minint, xx
	minint = wavemax(wave1d)
	variable i, j
	FindAPeak 0,1,5,wave1d
	if( V_Flag == 1)
		return 0.0
	endif
	
	for( j=0; j<=dimsize(wave1d,0); j+=1)
		if( wave1d[j] == wave1d[x2pnt(wave1d,V_peakX)] )
			i = j
			break
		endif
	endfor
	//print j
	duplicate/o/R=[0,j] wave1d, wave1d_temp
	for( j=0; j<=dimsize(wave1d_temp,0); j+=1)
		if( wave1d[j] == wavemin(wave1d_temp) )
			xx = j
			break
		endif
	endfor
	killwaves wave1d_temp
	
	//for( i=j; i<dimsize(wave1d,0); i+=1)
	//	if(wave1d[i] <= minint)
	//		minint = wave1d[i]
	//		xx = i
	//	else
	//		break
	//	endif
	//endfor
	return xx*dimdelta(wave1d,0)
end

function gauss3d(x,y,z,coef)
	variable x,y,z
	wave coef
	variable x0, y0, z0
	variable sx, sy, sz
	variable cxy, cxz, cyz
	x0 = coef[6]
	y0 = coef[7]
	z0 = coef[8]
	sx = coef[0]
	sy = coef[1]
	sz = coef[2]
	cxy = coef[3]
	cxz = coef[4]
	cyz = coef[5]
	//return exp( -1/(2*(1-cxy^2)*(1-cxz^2)*(1-cyz^2)) * (  ((x-x0)/sx)^2 + ((y-y0)/sy)^2 + ((z-z0)/sz)^2 - (2*cxy*(x-x0)*(y-y0))/(sx*sy) - (2*cxz*(x-x0)*(z-z0))/(sx*sz) - (2*cyz*(y-y0)*(z-z0))/(sy*sz) ))
	return exp(  -1/(2* (-1+cxy^2+cxz^2+cyz^2-2*cxy*cxz*cyz)) * ( (cyz^2-1)*(x-x0)^2/sx^2 + (cxz^2-1)*(y-y0)^2/sy^2 + (cxy^2-1)*(z-z0)^2/sz^2 + (2*cxy*(x-x0)*(y-y0)-2*cxz*cyz*(x-x0)*(y-y0))/(sx*sy) + (2*cxz*(x-x0)*(z-z0)-2*cxy*cyz*(x-x0)*(z-z0))/(sx*sz) + (2*cyz*(y-y0)*(z-z0)-2*cxy*cxz*(y-y0)*(z-z0))/(sy*sz) )  )
end

function hanning3d(w,x,y,z)
	// Input x,y,z as pixels from 0 to npix (i.e. as p,q,r)
	wave w
	variable x,y,z
	variable xn, yn, zn, tx, ty, tz, xrot, yrot, zrot
	variable dx, dy, dz
	variable xe,ye,ze
	variable cxy, cxz, cyz
	xe = 0.25
	ye = 1
	ze = 0.5
	cxy = 0
	
	dx = dimdelta(w,0)
	dy = dimdelta(w,1)
	dz = dimdelta(w,2)
	xn = (x - dimsize(w,0)/2)*dx
	yn = (y - dimsize(w,1)/2)*dy
	zn = (z - dimsize(w,2)/2)*dz
	
	xrot = 0//pi/3
	yrot = 0//pi/8
	zrot = pi/3
	tx = pi*xn/xe
	ty = pi*yn/ye
	tz = pi*zn/ze
	wave M_product
	make/o matRotx= { {1,0,0}, {0, cos(xrot), sin(xrot)}, {0, -sin(xrot), cos(xrot)} }
	make/o matRoty= { {cos(yrot), 0, -sin(yrot)}, {0,1,0}, {sin(yrot), 0, cos(yrot)} }
	make/o matRotz= { {cos(zrot), sin(zrot), 0}, {-sin(zrot), cos(zrot), 0}, {0,0,1} }
	//print matRotx
	//print matRoty
	//print matRotz
	make/o matK = {tx, ty, tz}
	MatrixMultiply matRotz/T, matRoty/T, matRotx/T, matK
	tx = M_product[0]
	ty = M_product[1]
	tz = M_product[2]
	
	//make/o matRot={{cos(xrot), -sin(xrot)},{sin(xrot), cos(xrot)}}
	//make/o matK = {xn,yn}
	//MatrixMultiply matRot/T, matK
	//xn = M_product[0]
	//yn = M_product[1]
	
	//make/o matRot={{cos(yrot), -sin(yrot)},{sin(xrot), cos(yrot)}}
	//make/o matK = {yn,zn}
	//MatrixMultiply matRot/T, matK
	//yn = M_product[0]
	//zn = M_product[1]
	
	//if( xn > xe || xn < -xe)
	//	return 0
	//elseif( yn > ye || yn < -ye)
	//	return 0
	//elseif( zn > ze || zn < -ze)
	//	return 0
	//endif

	return real( 0.25*(exp(cmplx(0,tx)) + exp(-cmplx(0,tx))) +0.5) * real( 0.25*(exp(cmplx(0,ty)) + exp(-cmplx(0,ty))) +0.5) * real( 0.25*(exp(cmplx(0,tz)) + exp(-cmplx(0,tz))) +0.5)
	//return 0.5*(cos(tx)+1) * 0.5*(cos(ty)+1) * 0.5*(cos(tz)+1)
	// - (2*cxy*(x)*(y))/(xe*ye)
end

function my_gauss1d(x,x0,sx)
	variable x,x0,sx
	return exp( -(x-x0)^2/sx^2)
end

function my_gauss2d(x,y,x0,y0,sx,sy,cxy)
	variable x,y,x0,y0,sx,sy,cxy
	return exp( -1/(2*(1-cxy^2)) * ( (x-x0)^2/sx^2 + (y-y0)^2/sy^2 - (2*cxy*(x-x0)*(y-y0))/(sx*sy) ) )
end

function tukey(kx, ky)
	variable kx,ky
	variable w, kc, k, kxn, kyn, a, eps
	//a = 0
	//make/o matRot={{cos(a), -sin(a)},{sin(a), cos(a)}}
	//make/o matK = {kx,ky}
	//MatrixMultiply matRot/T, matK
	//wave M_product
	
	//kxn = M_product[0] - 127
	//kyn = M_product[1] - 127
	eps = 6
	kxn = kx - 127
	kyn = ky - 127
	k = sqrt(kxn^2 + kyn^2)
	kc = 108* sqrt( (kxn^2+kyn^2)/(kxn^2+kyn^2*eps) )
	w = 20
	if(abs(k) < kc)
		return 1
	endif
	if(kc + w <= abs(k))
		return 0
	endif
	return cos( (pi*(abs(k) - kc))/(2*w) )^2
end


function ShowXYZLineProfiles(dat, particles, cut)
	wave dat, particles
	variable cut
	variable r
	Make/O/N=4 coefx, coefy, coefz
	variable delta_x, delta_y, delta_z
	Duplicate/O/R=(0,19) particles, output_particles
	for( r=0; r<20; r+=1)
	delta_x = particles[r][1] - particles[r][0]
	delta_y = particles[r][3] - particles[r][2]
	delta_z = particles[r][5] - particles[r][4]
	ImageLineProfile3D(dat,{0,511},{particles[r][8],particles[r][8]},{particles[r][9],particles[r][9]},0,0)
	if(cut)
		Duplicate/O/R=(particles[r][7]-2*delta_x,particles[r][7]+2*delta_x) LineIntensity, LineIntensityx
	else
		Duplicate/O LineIntensity, LineIntensityx
	endif
	//Display LineIntensityx //
	CurveFit/Q/O/M=2/W=0 gauss, kwCWave=coefx, LineIntensityx/D
	coefx[3] = particles[r][7]
	CurveFit/Q/M=2/W=0/H="0010"/TBOX=(0x300) gauss, kwCWave=coefx, LineIntensityx/D
	ImageLineProfile3D(dat,{particles[r][7],particles[r][7]},{0,511},{particles[r][9],particles[r][9]},0,0)
	if(cut)
		Duplicate/O/R=(particles[r][8]-2*delta_y,particles[r][8]+2*delta_y) LineIntensity, LineIntensityy
	else
		Duplicate/O LineIntensity, LineIntensityy
	endif
	//Display LineIntensityy //
	CurveFit/Q/O/M=2/W=0 gauss, kwCWave=coefy, LineIntensityy/D
	coefy[3] = particles[r][8]
	CurveFit/Q/M=2/W=0/H="0010"/TBOX=(0x300) gauss, kwCWave=coefy, LineIntensityy/D
	ImageLineProfile3D(dat,{particles[r][7],particles[r][7]},{particles[r][8],particles[r][8]},{0,511},0,0)
	if(cut)
		Duplicate/O/R=(particles[r][9]-2*delta_z,particles[r][9]+2*delta_z) LineIntensity, LineIntensityz
	else
		Duplicate/O LineIntensity, LineIntensityz
	endif
	//Display LineIntensityz //
	CurveFit/Q/O/M=2/W=0 gauss, kwCWave=coefz, LineIntensityz/D
	coefz[3] = particles[r][9]
	CurveFit/Q/M=2/W=0/H="0010"/TBOX=(0x300) gauss, kwCWave=coefz, LineIntensityz/D
	make/O/N=3 mean_width
	variable mean_total_width
	mean_width[0] = coefx[3]
	mean_width[1] = coefy[3]
	mean_width[2] = coefz[3]
	mean_total_width = ceil(wavemax(mean_width)*1.2)
	print r, mean_total_width
	
	output_particles[r][0] = floor((output_particles[r][7] - mean_total_width))///2)
	output_particles[r][1] = floor((output_particles[r][7] + mean_total_width))///2)
	output_particles[r][2] = floor((output_particles[r][8] - mean_total_width))///2)
	output_particles[r][3] = floor((output_particles[r][8] + mean_total_width))///2)
	output_particles[r][4] = floor((output_particles[r][9] - mean_total_width))///2)
	output_particles[r][5] = floor((output_particles[r][9] + mean_total_width))///2)
	//output_particles[r][7] = floor(output_particles[r][7]/2)
	//output_particles[r][8] = floor(output_particles[r][8]/2)
	//output_particles[r][9] = floor(output_particles[r][9]/2)
	endfor
end

function ShowXYZLineProfile(dat, particles, r, cut)
	wave dat, particles
	variable cut
	variable r
	Make/O/N=4 coefx, coefy, coefz
	variable delta_x, delta_y, delta_z
	Duplicate/O/R=(0,19) particles, output_particles
	delta_x = particles[r][1] - particles[r][0]
	delta_y = particles[r][3] - particles[r][2]
	delta_z = particles[r][5] - particles[r][4]
	ImageLineProfile3D(dat,{0,511},{particles[r][8],particles[r][8]},{particles[r][9],particles[r][9]},0,0)
	if(cut)
		Duplicate/O/R=(particles[r][7]-2*delta_x,particles[r][7]+2*delta_x) LineIntensity, LineIntensityx
	else
		Duplicate/O LineIntensity, LineIntensityx
	endif
	Display LineIntensityx //
	CurveFit/Q/O/M=2/W=0 gauss, kwCWave=coefx, LineIntensityx/D
	coefx[3] = particles[r][7]
	CurveFit/Q/M=2/W=0/H="0010"/TBOX=(0x300) gauss, kwCWave=coefx, LineIntensityx/D
	ImageLineProfile3D(dat,{particles[r][7],particles[r][7]},{0,511},{particles[r][9],particles[r][9]},0,0)
	if(cut)
		Duplicate/O/R=(particles[r][8]-2*delta_y,particles[r][8]+2*delta_y) LineIntensity, LineIntensityy
	else
		Duplicate/O LineIntensity, LineIntensityy
	endif
	Display LineIntensityy //
	CurveFit/Q/O/M=2/W=0 gauss, kwCWave=coefy, LineIntensityy/D
	coefy[3] = particles[r][8]
	CurveFit/Q/M=2/W=0/H="0010"/TBOX=(0x300) gauss, kwCWave=coefy, LineIntensityy/D
	ImageLineProfile3D(dat,{particles[r][7],particles[r][7]},{particles[r][8],particles[r][8]},{0,511},0,0)
	if(cut)
		Duplicate/O/R=(particles[r][9]-2*delta_z,particles[r][9]+2*delta_z) LineIntensity, LineIntensityz
	else
		Duplicate/O LineIntensity, LineIntensityz
	endif
	Display LineIntensityz //
	CurveFit/Q/O/M=2/W=0 gauss, kwCWave=coefz, LineIntensityz/D
	coefz[3] = particles[r][9]
	CurveFit/Q/M=2/W=0/H="0010"/TBOX=(0x300) gauss, kwCWave=coefz, LineIntensityz/D
	make/O/N=3 mean_width
	variable mean_total_width
	mean_width[0] = coefx[3]
	mean_width[1] = coefy[3]
	mean_width[2] = coefz[3]
	mean_total_width = ceil(wavemax(mean_width)*1.2)
	print r, mean_total_width
end

function try_automagic_spot_fix(dat, particles, particle_coefs, r, cut)
	// This function uses ShowXYZLineProfile functions to generate new sigmas for the fit.
	// It then allows the other parameters (center positions and tilting) to change to optimize the fit.
	wave dat, particles, particle_coefs
	variable cut
	variable r
	wave out = $"out"
	wave out_fortran = $"out_fortran"
	ShowXYZLineProfile(dat, particles, r, cut)
	wave W_coef, coefx, coefy, coefz
	W_coef[0] = coefx[3]
	W_coef[1] = coefy[3]
	W_coef[2] = coefz[3]
	print W_coef
	regenerateSpotParameters(dat, r, particles, particle_coefs, W_coef)
	create_output_particle(dat, particles, particle_coefs, r, out, out_fortran)
end

// This function creates 3 waves in the x y and z directions that
// are line profiles through the center spot.
// I don't remember the specifics anymore, but it works decently.
// xwave, ywave, and zwave are 2-tuples.
// cut is T/F
Function ImageLineProfile3D(dat,xwave,ywave,zwave,cut,center)
	wave dat,xwave,ywave,zwave
	variable cut, center
	variable i,j,k,xx,yy,zz,dt,t0,tf,aa,bb,cc,npix
	npix = dimsize(dat,0)
	//npix = max(max( abs(xwave[1]-xwave[0]), abs(ywave[1]-ywave[0]) ), abs(zwave[1]-zwave[0]) )
	
	// assume delta_t = 1 from x0 to xf to start with
	aa = xwave[1]-xwave[0]
	bb = ywave[1]-ywave[0]
	cc = zwave[1]-zwave[0]
	//print aa,bb,cc, npix

	i = 0
	do
		xx = round(aa*i/npix + xwave[0])
		yy = round(bb*i/npix + ywave[0])
		zz = round(cc*i/npix + zwave[0])
		//print xx,yy,zz
		if( xx > 0 && yy > 0 && zz > 0 && xx < npix-1 && yy < npix-1 && zz < npix-1)
			i -= 1
		else
			break
		endif
	while(1)
	t0 = i
	//print round(aa*t0/npix + xwave[0]), round(bb*t0/npix + ywave[0]), round(cc*t0/npix + zwave[0])
	i = 0
	do
		xx = round(aa*i/npix + xwave[0])
		yy = round(bb*i/npix + ywave[0])
		zz = round(cc*i/npix + zwave[0])
		//print xx,yy,zz
		if( xx > 0 && yy > 0 && zz > 0 && xx < npix-1 && yy < npix-1 && zz < npix-1)
			i += 1
		else
			break
		endif
	while(1)
	tf = i
	if(tf == t0)
		tf = npix
	endif
	//print round(aa*tf/npix + xwave[0]), round(bb*tf/npix + ywave[0]), round(cc*tf/npix + zwave[0])
	dt = (tf-t0)/npix
	//print t0,tf,dt
	
	//Make/O/N=(ceil(tf-t0+1)) LineIntensity
	Make/O/N=(npix) LineIntensity
	j = 0
	for( i=t0; i<tf; i+=dt)
		xx = round(aa*i/npix + xwave[0])
		yy = round(bb*i/npix + ywave[0])
		zz = round(cc*i/npix + zwave[0])
		//print xx,yy,zz
		LineIntensity[j] = dat[xx][yy][zz]
		j += 1
	endfor
	if(cut)
		//Make/O/N=51 temp
		//for( i=center-25; i<=center+25; i+=1)
		//	temp[i-center+25] = LineIntensity[i]
		//endfor
		Duplicate/O/R=(center-25,center+25) LineIntensity, temp
		Duplicate/O temp, LineIntensity
		killwaves temp
	endif
	//Display LineIntensity
end

function MaskFromXYZMinMax(dim1,dim2,dim3,xyzminmax)
	variable dim1,dim2,dim3
	wave xyzminmax
	
	Make/O/B/U/N=(dim1,dim2,dim3) particle_mask
	particle_mask = 0
	variable i,j,k,n
	for( n=0; n<dimsize(xyzminmax,0); n+=1)
		printf "Calculating particle %g...\r", n
		print xyzminmax[n][0],  xyzminmax[n][1]
		print xyzminmax[n][2],  xyzminmax[n][3]
		print xyzminmax[n][4],  xyzminmax[n][5]
		
		for(i=xyzminmax[n][0]; i<=xyzminmax[n][1]; i+=1)
			for(j=xyzminmax[n][2]; j<=xyzminmax[n][3]; j+=1)
				for(k=xyzminmax[n][4]; k<=xyzminmax[n][5]; k+=1)
					particle_mask[i][j][k] = 255
				endfor
			endfor
		endfor
	endfor
end

function AppendGvecsToVk(ft, out, kexp, vkexp)
	wave ft, out, kexp, vkexp
	// You must run FindSpots(ft,15) first
	// Run the full fits on all the spots too to generate "out"
	
	setscale/P x,kexp[0],kexp[1]-kexp[0],vkexp
	//make/o/n=(dimsize(out,0)) gvecs, gvecs_id, gvecs_size
	make/o/n=(20) gvecs, gvecs_id, gvecs_size
	gvecs_id = p
	// if you generated "out" use the better below line, otherwise approximate with the one after
	gvecs_size = sqrt( out[p][11]^2 + out[p][12]^2 + out[p][13]^2 )/2  // use sigmas to get gvec size
	//gvecs_size = (out[p][1]-out[p][0])/2
	gvecs_size*=dimdelta(ft,0)
	variable rr, i
	// if you generated "out" use the below loop, otherwise just comment it
	//for(rr=0; rr<dimsize(out,0); rr+=1)
	//	out[rr][10] = sqrt( (dimsize(ft,0)/2-out[rr][17])^2 + (dimsize(ft,1)/2-out[rr][18])^2 + (dimsize(ft,2)/2-out[rr][19])^2 )*dimdelta(ft,0)
	//endfor
	gvecs = out[p][10]
	duplicate/o gvecs, gvecs_y
	setscale/P x,kexp[0],kexp[1]-kexp[0],vkexp
	variable x1, x2
	//for(rr=0; rr<dimsize(out,0); rr+=1)
	for(rr=0; rr<20; rr+=1)
		for(i=0; i<dimsize(kexp,0); i+=1)
			if( kexp[i] > gvecs[rr] )
				x1 = i - 1
				x2 = i
				//print rr, x1,x2, kexp[x1],kexp[x2], vkexp[x1],vkexp[x2], gvecs[rr]
				break
			endif
		endfor
		// if you generated "out" use the better below line, otherwise approximate with the one after
		gvecs_y[rr] = linearInterp(kexp[x1],vkexp[x1],kexp[x2],vkexp[x2],gvecs[rr])
		//gvecs_y[rr] = vkexp[x2pnt(vkexp,gvecs[rr])]
		//print gvecs_y[rr]
	endfor
	//Sort gvecs gvecs_y,gvecs,gvecs_id,gvecs_size
	Display  vkexp
	AppendToGraph gvecs_y vs gvecs
	ModifyGraph mode(gvecs_y)=3,marker(gvecs_y)=19,msize(gvecs_y)=1.5;DelayUpdate
	ModifyGraph rgb(gvecs_y)=(0,0,0)
	Legend/C/N=text0/J/F=0/A=MC "\\s(vkexp) V(k) experimental\r\\s(gvecs_y) These represent |g|\rfor various g-vectors. The\ry-axis has been scaled to";DelayUpdate
	AppendText "match V(k), the x-positions\rare |g|."
	Legend/C/N=text0/J "\\s(vkexp) V(k) experimental\r\\s(gvecs_y) These represent |g|\rfor various g-vectors found in\rthe FT. The y-axis has been ";DelayUpdate
	AppendText/N=text0 "scaled to match V(k), the x-\rpositions are |g|."
	ErrorBars gvecs_y X,wave=(gvecs_size,gvecs_size)
end

function linearInterp(x1,y1,x2,y2,x0)
	variable x1,x2,y1,y2,x0
	//print (y2-y1)/(x2-x1)*(x0-x1)+y1
	//print (y2-y1)/(x2-x1)*(x0-x2)+y2
	return (y2-y1)/(x2-x1)*(x0-x1)+y1
end

function cut_ft_smaller_for_isosurface(ft)
	wave ft
	make/o/n=(256,256,256) ft_iso
	variable i,j,k
	
	for(i=128; i<256+128; i+=1)
		printf "%g\r", i-128
		for(j=128; j<256+128; j+=1)
			//printf "  %g\r", j-128
			for(k=128; k<256+128; k+=1)
				//printf "    %g\r", k
				ft_iso[i-128][j-128][k-128] = ft[i][j][k]
			endfor
		endfor
	endfor
end

function window_edge_checking(out)
	wave out
	Make/O/N=9 W_coef
	variable i
	variable allowance = 0.5 // in percent
	for( i=0; i<dimsize(out,0); i+=1 )
		W_coef[0] = out[i][11]/2
		W_coef[1] = out[i][12]/2
		W_coef[2] = out[i][13]/2
		W_coef[3] = out[i][14]
		W_coef[4] = out[i][15]
		W_coef[5] = out[i][16]
		W_coef[6] = out[i][17]
		W_coef[7] = out[i][18]
		W_coef[8] = out[i][19]
		
		if( gauss3d(out[i][0],out[i][8],out[i][9],W_coef)*100 > allowance || gauss3d(out[i][7],out[i][2],out[i][9],W_coef)*100 > allowance || gauss3d(out[i][7],out[i][8],out[i][4],W_coef)*100 > allowance || gauss3d(out[i][1],out[i][8],out[i][9],W_coef)*100 > allowance || gauss3d(out[i][7],out[i][3],out[i][9],W_coef)*100 > allowance || gauss3d(out[i][7],out[i][8],out[i][5],W_coef)*100 > allowance )
		
		printf( "%g  " ) i //, gauss3d(out[i][0],out[i][8],out[i][9],W_coef)*100, gauss3d(out[i][7],out[i][2],out[i][9],W_coef)*100, gauss3d(out[i][7],out[i][8],out[i][4],W_coef)*100, gauss3d(out[i][1],out[i][8],out[i][9],W_coef)*100, gauss3d(out[i][7],out[i][3],out[i][9],W_coef)*100, gauss3d(out[i][7],out[i][8],out[i][5],W_coef)*100
		
		endif
		
		if( gauss3d(out[i][0],out[i][8],out[i][9],W_coef)*100 > allowance )
			printf( "xi gauss = %g  ") gauss3d(out[i][0],out[i][8],out[i][9],W_coef)*100
		endif
		if( gauss3d(out[i][7],out[i][2],out[i][9],W_coef)*100 > allowance )
			printf( "yi gauss = %g  ") gauss3d(out[i][7],out[i][2],out[i][9],W_coef)*100
		endif
		if( gauss3d(out[i][7],out[i][8],out[i][4],W_coef)*100 > allowance )
			printf( "zi gauss = %g  ") gauss3d(out[i][7],out[i][8],out[i][4],W_coef)*100
		endif
		
		if( gauss3d(out[i][1],out[i][8],out[i][9],W_coef)*100 > allowance )
			printf( "xf gauss = %g  ") gauss3d(out[i][1],out[i][8],out[i][9],W_coef)*100
		endif
		if( gauss3d(out[i][7],out[i][3],out[i][9],W_coef)*100 > allowance )
			printf( "yf gauss = %g  ") gauss3d(out[i][7],out[i][3],out[i][9],W_coef)*100
		endif
		if( gauss3d(out[i][7],out[i][8],out[i][5],W_coef)*100 > allowance )
			printf( "zf gauss = %g  ") gauss3d(out[i][7],out[i][8],out[i][5],W_coef)*100
		endif
		
		if( gauss3d(out[i][0],out[i][8],out[i][9],W_coef)*100 > allowance || gauss3d(out[i][7],out[i][2],out[i][9],W_coef)*100 > allowance || gauss3d(out[i][7],out[i][8],out[i][4],W_coef)*100 > allowance || gauss3d(out[i][1],out[i][8],out[i][9],W_coef)*100 > allowance || gauss3d(out[i][7],out[i][3],out[i][9],W_coef)*100 > allowance || gauss3d(out[i][7],out[i][8],out[i][5],W_coef)*100 > allowance )
		
		printf( "\r" ) //, gauss3d(out[i][0],out[i][8],out[i][9],W_coef)*100, gauss3d(out[i][7],out[i][2],out[i][9],W_coef)*100, gauss3d(out[i][7],out[i][8],out[i][4],W_coef)*100, gauss3d(out[i][1],out[i][8],out[i][9],W_coef)*100, gauss3d(out[i][7],out[i][3],out[i][9],W_coef)*100, gauss3d(out[i][7],out[i][8],out[i][5],W_coef)*100
		
		endif
		
	endfor
end

function compare_centers(out,particles)
	wave out, particles
	// Calculates the distance between the two center positions. Units are number of pixels.
	variable rr, x1,x2,y1,y2,z1,z2, dist
	for( rr=0; rr<dimsize(out,0); rr+=1)
		compare_center(out,particles,rr)
	endfor
end

function compare_center(out,particles,rr)
	wave out, particles
	variable rr
	variable x1,x2,y1,y2,z1,z2, dist
	x1 = out[rr][7]
	y1 = out[rr][8]
	z1 = out[rr][9]
	x2 = particles[rr][7]
	y2 = particles[rr][8]
	z2 = particles[rr][9]
	dist = sqrt( (x1-x2)^2 + (y1-y2)^2 + (z1-z2)^2 )
	printf( "Dist for spot %g = %g\r") rr, dist
end

function try_automagic_fix_holdcenter(dat, particles, particle_coefs, r, W_coef)
	// This function uses the center particle position from "particles" and allows the other parameters to vary in the fitting.
	// This function also doesn't work.
	wave dat, particles, particle_coefs, W_coef
	variable r
	wave out = $"out"
	wave out_fortran = $"out_fortran"
	W_coef[6] = particles[r][7]
	W_coef[7] = particles[r][8]
	W_coef[8] = particles[r][9]
	gauss_fit_particle(dat, particles, r, particle_coefs, 0, 1) // set hold == 1 so that the centers are held constant, but don't hold sigmas constant
	create_output_particle(dat, particles, particle_coefs, r, out, out_fortran)
end

function regenerate_and_show_spot(ft, particles, particle_coefs, out, out_fortran, r)
	wave ft, particles, particle_coefs, out, out_fortran
	variable r
	gauss_fit_particle(ft, particles, r, particle_coefs, 0, 0)  // change sigma hold here to 1 if modifying W_coef manually
	create_output_particle(ft, particles, particle_coefs, r, out, out_fortran)
	create_windowed_spots(ft, out, r)
	//NewImage  root:wind_im
	//WMAppend3DImageSlider();
end

function search_for_multiple_spots(ft, particles, particle_coefs, out, out_fortran)//, r)
	wave ft, particles, particle_coefs, out, out_fortran
	variable r
	for( r=0; r<dimsize(particles,0); r+=1)

	regenerate_and_show_spot(ft, particles, particle_coefs, out, out_fortran, r)
	wave wind_im
	ImageThreshold /Q/I/M=1 wind_im
	variable i, thresh
	thresh = V_threshold
	for( i=1.0; i<=3.0; i+=0.1)
		ImageThreshold /Q/I/T=(thresh*i) wind_im
		ImageAnalyzeParticles/A=0 stats M_ImageThresh
		if( dimsize(M_3DParticleInfo,0) > 1)
			printf( "    Spot %g may have more than one g-vector in it, you should check. Use Threshold = %g.\r" ) r, thresh*i
		endif
	endfor
	
	endfor
end