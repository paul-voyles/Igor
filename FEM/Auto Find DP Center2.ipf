#pragma rtGlobals=1		// Use modern global access method.

// Functions for automated center-finding in ring diffraction patterns, based
// on the algorithm in T. C. Peterson et al. Ultramicoscopy 103, 275 (2005).
// This algorithm assumes that the pattern is circularlly symmetric about it's
// true center; therefore, a pattern rotated about the true center should be identical
// to the original pattern.  A pattern rotated about some other point will be difference.
// The mean square deviation is adopted as a measure of that difference, and it
// is systematically minimized as a function of position in the image to find the
// true center.
//
// v1.0 pmv 06-06-06
//
//  pmv 06-28-06: changed MeanSquareDiff to calculate mean sq. diff per pixel in common 
//        to the rotated and original image.  Otherwise, a very wrong center could give
//        a small difference because only a few pixels overlap.  This is important for images 
//        where the center is outside the image.


// Automatically find the center of the diffraction pattern.  Arguments are:
// 	the image wave containing the DP
// 	the rectangle in which to search for the center defined by (xmin, ymin) and (xmax, ymax)
//	the size in pixels of the centers search wave 
//	the radius out from the center to consider ()
//	the rotation angle 
// size = 7 seems about optimum.  The radius should be either the diameter of the zero beam 
// or the radius of the first ring.  A rotation angle of 30 deg seems to work if the entire ring is
// present in the pattern.  If only a fraction of the total arc is availabe, a smaller angle might be better.
function/S AutoCenterFind(im, xmin, ymin, xmax, ymax, size, radius, angle)
	wave im
	variable xmin, ymin, xmax, ymax
	variable size, radius, angle
	
	string centers_w
	variable done = 0, iteration = 1
	
	do
		sprintf centers_w, "centers_%d", iteration
		Make/O/N=(size, size) $centers_w
		wave cent = $centers_w
		done = CentersWaveScale(cent, xmin, xmax, ymin, ymax)
		printf "Working on iteration %d with x step %d and ystep %d.\r", iteration, DimDelta(cent, 0), DimDelta(cent, 1)
		cent = RecenterTest(im, x, y, radius, angle)
		WaveStats/Q cent
		printf "\tProvisional center is (%d, %d)\r", V_minRowLoc, V_minColLoc
		xmin = V_minRowLoc - DimDelta(cent, 0)
		xmax = V_minRowLoc + DimDelta(cent, 0)
		ymin = V_minColLoc - DimDelta(cent, 0)
		ymax = V_minColLoc + DimDelta(cent, 0)
		iteration+=1
	while(!done)
	
	printf "Final center is (%d, %d)\r", V_minRowLoc, V_minColLoc
	
	sprintf centers_w, "%d,%d", V_minRowLoc, V_minColLoc	
	return centers_w
end
	
// Calculate the mean-square deviation assuming the image center is 
// at (xc, yc)	
function RecenterTest(dat, xc, yc, radius, angle)
	wave dat
	variable xc, yc, radius, angle
	
	//Recenter(dat, xc, yc)   //to avoid bug. not needed right now  - JWH 
	wave centered = $"centered"

	variable crop_size = ceil(2.1*radius)
	Make/O/N=(crop_size, crop_size) crop_centered
	variable x_offset, y_offset
	x_offset = floor(DimSize(centered, 0)/2)-crop_size/2
	y_offset = floor(DimSize(centered, 1)/2)-crop_size/2
	crop_centered = centered[p+x_offset][q+y_offset]
	ImageRotate/A=(angle) crop_centered
	wave crop_centered_rotated = $"M_RotatedImage"
	
	variable msd = MeanSquareDiff(crop_centered, crop_centered_rotated, radius)
	Killwaves centered, crop_centered, M_RotatedImage, square_diff
	
	return msd
end

// Calculate the mean square difference between two images, out to
// a certain radius, and normalized by the number of pixels in common
// to the original and rotated images.
function MeanSquareDiff(im1, im2, radius)
	wave im1, im2
	variable radius
	
	make/O/N=(radius*2, radius*2) square_diff
	variable xc1 = floor(DimSize(im1, 0) / 2) - radius
	variable yc1 = floor(DimSize(im1, 1) / 2) - radius
	variable xc2 = floor(DimSize(im2, 0) / 2) - radius
	variable yc2 = floor(DimSize(im2, 1) / 2) - radius
	
	square_diff = ( im1[xc1 + p][yc1 + q] - im2[xc2 + p][yc2 + q])
	FastOp square_diff = square_diff * square_diff
	square_diff = (( (p-radius)^2 + (q-radius)^2) < radius^2 ? square_diff[p][q] : NaN)
	WaveStats/Q/M=1 square_diff
	square_diff = ( (NumType(square_diff[p][q]) == 2) ? 0 : square_diff[p][q])
	
	return sum(square_diff, -inf, inf) / V_npnts

end

// Recenter the DP centers at the center of the n x n image - JWH
function Recenter(data, center_x, center_y, n)
	wave data
	variable center_x, center_y  //center of the DP
	variable n  // number of pixels in 1D in the output image.
	
	variable size_x, size_y, size_z, shift_x, shift_y
	size_x = DimSize(data, 0)
	size_y = DimSize(data, 1)
	size_z = DimSize(data, 2)
	
	//shift_x = floor(size_x/2) - center_x
	//shift_y = floor(size_y/2) - center_y
	shift_x = floor(n/2) - center_x
	shift_y = floor(n/2) - center_y
	
	//Make/O/N=(size_x + 2.0*abs(shift_x), size_y + 2.0*abs(shift_y)) centered
	Make/O/N=(n, n, size_z) centered
	//wave centered = $"centered"
	centered = nan

	variable new_cx, new_cy
	new_cx = DimSize(centered, 0)
	new_cy = DimSize(centered, 1)
	new_cx = floor(new_cx/2)
	new_cy = floor(new_cy/2)

	variable top,left,bottom, right
	//bottom = (new_cy - floor(size_y/2)) + shift_y
	bottom = shift_y
	top = bottom + size_y
	//left = (new_cx - floor(size_x/2)) + shift_x
	left = shift_x
	right = left + size_x
	
	variable i
	i=0
	for (i=0; i<size_z; i+=1)
		centered[left,right][bottom, top][i] = data[p-left][q-bottom][i]   
	endfor
	
	string name_wave
	name_wave = NameofWave(data)
	name_wave = name_wave+"_c"
	duplicate centered, $name_wave
	killwaves centered
end





function CentersWaveScale(centers, xmin, xmax, ymin, ymax)
	wave centers
	variable xmin, xmax, ymin, ymax
	
	variable nx, ny, xstep, ystep, one_pixel_x, one_pixel_y
	one_pixel_x = 0
	one_pixel_y = 0
	nx = DimSize(centers, 0) - 1
	ny = DimSize(centers, 1) - 1
	
	xstep = floor( (xmax - xmin) / nx )
	ystep = floor( (ymax - ymin) / ny )
	
	if(xstep <= 1)
		printf "x step size 1 pixel.\r"
		xstep = 1
		one_pixel_x = 1
		Redimension/N=( (xmax - xmin + 1), -1) centers
	endif
	if(ystep <= 1)
		printf "y step size 1 pixel.\r"
		ystep = 1
		one_pixel_y = 1
		Redimension/N=( -1, (ymax - ymin + 1) ) centers
	endif
	
	SetScale/P x xmin, xstep, "", centers
	SetScale/P y ymin, ystep, "", centers
	
	return (one_pixel_x & one_pixel_y)
end

