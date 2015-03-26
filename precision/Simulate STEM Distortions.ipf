#pragma rtGlobals=3		// Use modern global access method and strict wave access.
// functions for generating simulated STEM distortions, including sample drift,
// scan noise, and acoustic noise
//
// begun 01-20-2014, pmv
//
// Change ac_time wave to double precision, 02-17-14 pmv
// Fix standard deviation vs variance bug in FramePNoise 02-27-14 pmv
// Add output of distortions without Poisson noise as well as noise-added version 04-15-14 pmv

// input is the image to distort in im and the number of distorted frames to generate in
// steps.  Input image calibration is used: distances and velocities are assumed to be in
// the input image calibration units.  Input image is assumed on an absolute intensity scale
// 0 to 1.  Details of the applied distortion are set by variables inside the funciton, not input
// parameters.
function STEMDistortionStack(im, steps)
	wave im
	variable steps
	
	// parameters describing the image series
	variable nx = 128, ny = 128, npix = nx*ny
	variable x0=10.112, y0=10.112			// starting point in the full frame for the resampled sub-frame
	variable dwell = 16e-6 				// pixel dwell time in seconds
	variable flyback = 60e-6 			// line flyback time in seconds
	variable interframe = 1e-3 			// frame-to-frame delay time in seconds
	variable dose_per_pixel = 832		// incident electrons per pixel in each frame

	// variables describing the distortions
	variable slow_rand_f = 0.08  	// frequency (Hz) of random walk excuted by sample.
							  	// should have a period longer than the frame time
	variable slow_rand_mag = 0.3 	// maximum magntidue of the random walk, in the distance units
								// of the (calibrated) input image
	variable vel_x = 0.007, vel_y = 0.007		// x- and y- components of constant velocity linear drift.  Velocity in units
									// of distance / second, for calibrated distance units of the input image

	// waves for adjustable sets of high-frequency, zero-mean distortions
	// line frequency 60 Hz, plus harmonics at 30 Hz, 120 Hz, and 180 Hz.  Magntidue for harmonics
	// should be smaller than the main line at 60 Hz.
	// There is also a small bump in the instability spectrum around 400 Hz.  The bandwidth of the image
	// single in principle goes up to frequency = 1/(2*dwell) ~ 30 kHz.
	Make/o/N=(6) freq, mag
	freq = {30, 120, 160, 430, 1000, 2725}
	mag = {0.015, 0.01, 0.007, 0.004, 0.004, 0.004}  // distortions based on UW-Titan measurements
	
	// mag = {0.05, 0.01, 0.007, 0.01, 0.004, 0.004}  // much larger distortions

	
	// other useful variables
	variable dx = DimDelta(im, 0), dy = DimDelta(im, 1)
	variable i

	// set up the perfect grid sampling
	PerfectGrid(x0, y0, dx, dy, nx, ny)
	wave grid_x = $"grid_x"
	wave grid_y = $"grid_y"

	// set up the acquisition timing wave
	AcTime(nx, ny, steps, dwell, flyback, interframe)
	wave ac_time = $"ac_time"
	wavestats/q ac_time
	variable total_time = V_max

	// set up the wave in which to accumulate all the shifts.
	// Initialized with the perfect grid in the first frame
	make/o/n=(npix, steps) cum_sx, cum_sy
	cum_sx[][0] = grid_x[p]
	cum_sy[][0] = grid_y[p]	

	// generate the time series shifts for the slow random walk and apply them, accumulating
	// shifts from frame to frame, seeded by the perfect grid in the first frame
	TimeSeriesShift(total_time, slow_rand_f, slow_rand_mag)
	wave tseries_x = $"tseries_x"
	wave tseries_y = $"tseries_y"
	CumulativeShiftPerFrame(ac_time, tseries_x, tseries_y, cum_sx, cum_sy)
	
	// add a constant velocity linear drift
	FixedVelocityShift(ac_time,cum_sx, cum_sy, vel_x, vel_y)	

	// add zero-mean shifts from electronic noise and jitter
	for(i=0; i<numpnts(freq); i+=1)
		TimeSeriesShift(total_time, freq[i], mag[i])
		ShiftPerFrame(ac_time, tseries_x, tseries_y, cum_sx, cum_sy)
	endfor
	
	// apply all the accumulated shfits
	ResampleFrames(im, nx, ny, steps, cum_sx, cum_sy)
	wave resample_st = $"resample_st"
	
	duplicate/o resample_st resample_st_pn
	resample_st *= dose_per_pixel
	
	// add Poisson noise
	for(i=0; i<steps; i+=1)
		ImageTransform/P=(i)/PTYP=0 getPlane resample_st_pn
		wave M_ImagePlane = $"M_ImagePlane"
		FramePNoise(M_ImagePlane, dose_per_pixel)
		ImageTransform/P=(i)/PTYP=0/D=M_ImagePlane setPlane resample_st_pn
	endfor
	  
	Killwaves M_ImagePlane, ac_time, grid_x, grid_y, cum_sx, cum_sy
	Killwaves tseries_x, tseries_y, freq, mag
end


function ResampleImage(im, nx, ny, new_x, new_y)
	wave im
	variable nx, ny
	wave new_x, new_y
	
	if( (numpnts(new_x) != nx*ny) || (numpnts(new_y) != nx*ny) )
		printf "Error: non-conforming output image size and pixel position list.\r"
		return 0
	endif
	
	make/O/N=(nx, ny) resample_im
	
	resample_im = interp2D(im, new_x[Index2dTo1d(nx, ny, p, q)], new_y[Index2DTo1d(nx, ny, p, q)])
	
	SetScale/P x DimOffset(im, 0), DimDelta(im, 0), resample_im
	SetScale/P y DimOffset(im, 1), DimDelta(im, 1), resample_im
	
end

function Index2dTo1d(nx, ny, pp, qq)
	variable nx, ny, pp, qq
	
	return pp+qq*nx
	
end

function PerfectGrid(x0, y0, dx, dy, nx, ny)
	variable x0, y0, dx, dy, nx, ny
	
	make/o/n=(nx*ny) grid_x, grid_y
	
	grid_x = x0 + dx*mod(p, nx)
	grid_y = y0 + dy*floor(p / ny)
	
end

function Shift2D(dx, dy, new_x, new_y)
	variable dx, dy
	wave new_x, new_y
	
	new_x += dx
	new_y += dy
	
end

function ResampleFrames(im, nx, ny, steps, cum_x, cum_y)
	wave im
	variable nx, ny, steps
	wave cum_x, cum_y
	
	make/o/n=(nx*ny) new_x, new_y
	make/o/n=(nx, ny, steps) resample_st
	SetScale/P x 0, DimDelta(im, 0), resample_st
	SetScale/P y 0, DimDelta(im, 1), resample_st
	
	variable i
	for(i=0; i<steps; i+=1)
		new_x = cum_x[p][i]
		new_y = cum_y[p][i]
		ResampleImage(im, nx, ny, new_x, new_y)
		wave resample_im = $"resample_im"
		resample_st[][][i] = resample_im[p][q]
	endfor
	
	Killwaves resample_im, new_x, new_y
end

function ResampleRandomWalkPerFrame(im, x0, y0, nx, ny, mag, steps)
	wave im
	variable x0, y0, nx, ny, mag, steps
	
	PerfectGrid(x0, y0, DimDelta(im, 0), DimDelta(im,1), nx, ny)
	wave grid_x = $"grid_x"
	wave grid_y = $"grid_y"
	duplicate/O grid_x new_x
	duplicate/O grid_y new_y
	
	variable tot_vx, tot_vy, diff_vx, diff_vy
	tot_vx = 0
	tot_vy = 0

	Make/O/N=(nx, ny, steps) rand_walk_st
	Make/O/N=2 r_vec
	
	variable i, diff_mag, r_mag
	for(i=0; i<steps; i+=1)
		RandVec2D(r_vec, mag)		
		tot_vx += r_vec[0]
		tot_vy += r_vec[1]
		new_x = grid_x + tot_vx
		new_y = grid_y + tot_vy
		ResampleImage(im, nx, ny, new_x, new_y)
		wave res = $"resample_im"
		rand_walk_st[][][i] = res[p][q]
	endfor
	
end

function RandVec2D(vec, mag)
	wave vec
	variable mag	
	
	variable vx, vy, diff_mag, r_mag
	vx = enoise(1)
	vy = enoise(1)
	r_mag = mag*(enoise(0.5)+0.5)
	diff_mag = sqrt(vx^2 + vy^2)
	vec[0] = sqrt(r_mag / diff_mag)*vx
	vec[1] = sqrt(r_mag / diff_mag)*vy
	
end
		
function AcTime(nx, ny, steps, dwell, flyback, interframe)
	variable nx, ny, steps, dwell, flyback, interframe
	
	make/d/o/n=(nx*ny*steps) ac_time
	
	ac_time = 0
	variable i=1
	do
		ac_time[i] = ac_time[i-1] + dwell
		if(!mod(i, nx))
			ac_time[i] += flyback
		endif
		if(!mod(i, nx*ny))
			ac_time[i] += interframe
		endif
		i+=1
	while(i<numpnts(ac_time))
		
end
		
function TimeSeriesShift(total_time, freq, mag)
	variable total_time, freq, mag
	
	variable pts_per = 40
	make/o/N=(total_time*pts_per*freq) tseries_x, tseries_y
	SetScale/I x 0, total_time, tseries_x, tseries_y
	make/o/n=(pts_per) x_per, y_per
	
	make/o/n=2 r_vec
	variable i, tend
	// loop over period
	for(i=0; i<numpnts(tseries_x); i+=pts_per)
		RandVec2D(r_vec, mag)
		x_per = r_vec[0]*sin(2*Pi*p/pts_per)
		y_per = r_vec[1]*sin(2*Pi*p/pts_per)
		
		tend = i+(pts_per-1)
		if(tend >= numpnts(tseries_x))
			tend = numpnts(tseries_x) -1
		endif
		
		tseries_x[i,tend] = x_per[p-i]
		tseries_y[i,tend] = y_per[p-i]
	endfor
	
	killwaves r_vec, x_per, y_per
end

// tracks sample motions slower than the frame time
Function CumulativeShiftPerFrame(ac_time, tseries_x, tseries_y, cum_sx, cum_sy)
	wave ac_time, tseries_x, tseries_y, cum_sx, cum_sy
	
	wavestats/q ac_time
	if(V_max > rightx(tseries_x))
		printf "Total acquisition time is %g.  Total shfit time series time is %g.  Don't have enough shift time series data.\r", V_max, rightx(tseries_x)
		return 0
	endif  

	variable npix = DimSize(cum_sx, 0)
	variable steps = DimSize(cum_sx, 1)

	// set up the first frame
	cum_sx[][0] = cum_sx[p][0] + tseries_x(ac_time[p])
	cum_sy[][0] = cum_sy[p][0] + tseries_y(ac_time[p])
	

	// loop over frames
	variable i
	for(i=1; i<steps; i+=1)
		cum_sx[][i] = cum_sx[p][i-1] + tseries_x(ac_time[npix*i+p])
		cum_sy[][i] = cum_sy[p][i-1] + tseries_y(ac_time[npix*i+p])
	endfor

end	

// tracks net, fixed velocity drift.
// vel_x, vel_y are velocities, length/sec, in length units of the image scaling
function FixedVelocityShift(ac_time,cum_sx, cum_sy, vel_x, vel_y)
	wave ac_time, cum_sx, cum_sy
	variable vel_x, vel_y
	
	variable npix = Dimsize(cum_sx, 0)
	variable steps = DimSize(cum_sx, 1)
	
	variable i
	for(i=0; i<steps; i+=1)
		cum_sx[][i] += vel_x*ac_time[npix*i + p]
		cum_sy[][i] += vel_y*ac_time[npix*i + p]
	endfor
	
end
	

// tracks instabilities faster than the frame time
Function ShiftPerFrame(ac_time, tseries_x, tseries_y, cum_sx, cum_sy)
	wave ac_time, tseries_x, tseries_y, cum_sx, cum_sy
	
	wavestats/q ac_time
	if(V_max > rightx(tseries_x))
		printf "Total acquisition time is %g.  Total shfit time series time is %g.  Don't have enough shift time series data.\r", V_max, rightx(tseries_x)
		return 0
	endif  

	variable npix = DimSize(cum_sx, 0)
	variable steps = DimSize(cum_sx, 1)

	// loop over frames
	variable i
	for(i=0; i<steps; i+=1)
		cum_sx[][i] = cum_sx[p][i] + tseries_x(ac_time[npix*i+p])
		cum_sy[][i] = cum_sy[p][i] + tseries_y(ac_time[npix*i+p])
	endfor

end	



// Adds Poisson noise, assuming input image is on absolute scale 0 to 1
// i.e. it's a multislice simulation.  dose_per_pixel is electrons incident on
// the sample in one pixel in one frame
function FramePNoise(im, dose_per_pixel)
	wave im
	variable dose_per_pixel
	
	im = poissonNoise(dose_per_pixel*im[p][q])
	
end