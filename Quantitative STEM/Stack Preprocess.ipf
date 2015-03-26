#pragma rtGlobals=1		// Use modern global access method.
#include "Image Analysis Utilities"
//
// Functions for pre-processing STEM image stacks
//
// requires FEM Analysis Utilities.ipf
//
// FrameDifference: count the number of identical pixels in consecutive frame in stack st
// DeleteDupes: delete frames from stack st which contain a fraction of identical pixels 
//                      greater than thresh.  0.1 seems to work for thresh.  Makes a new stack
//                      dupe_del
//
// OneFramePS: extract one frame from the stack and calculate the 2D power spectrum
// StackWindowPS: calculate the integrated power spectrum in the box xmin->xmax, ymin->ymax
//                            for every image in the stack
// DeleteLowPS: remove images from the stack whose integrated power spectrum in the box
//                       is lower than thresh
//  
//  started 3/3/12, pmv


function FrameDifference(st)
	wave st
	
	make/o/n=(DimSize(st, 0), DimSize(st, 1)) dt
	make/o/n=(DimSize(st, 2)) diff
	
	variable j
	for(j=1; j<DimSize(st, 2); j+=1)
		dt  = st[p][q][j] - st[p][q][j-1]
		dt = 1/dt
		wavestats/q dt
		diff[j-1] = V_numINFs
	endfor
	
	diff /= (DimSize(st, 0)*DimSize(st,1))
	Killwaves dt
end

function DeleteDupes(st, thresh)
	wave st
	variable thresh
	
	Duplicate/O st, dupe_del
	make/o/n=(DimSize(st, 0), DimSize(st, 1)) dt
	variable diff
	
	variable i=1, j=0
	do
		// printf "Working on image %d\r", i
		dt = dupe_del[p][q][i] - dupe_del[p][q][i-1]
		dt = 1/dt
		wavestats/q dt
		diff = V_numINFs / (DimSize(st, 0)*DimSize(st, 1))
		
		if(diff >= thresh)
//			Imagetransform/O/P=(i) removeZplane dupe_del
			Deletepoints/M=2 i, 1, dupe_del
			j+=1
			printf "Deleted image %d\r" i
		else
			i+=1
		endif	
	while(i<DimSize(dupe_del, 2))
	
	printf "Removed %d duplicate frames.\r", j

	Killwaves dt
end

function OneFramePS(st, n)
	wave st
	variable n
	
	imagetransform/p=(n) getplane st
	wave im = $"M_ImagePlane"
	redimension/s im
	PowerSpectrum2D(im)
	
	killwaves M_ImagePlane
	
end

function StackWindowPS(st, xmin, xmax, ymin, ymax)
	wave st
	variable xmin, xmax, ymin, ymax
	
	Make/o/N=(DimSize(st, 2)) stack_ps
	
	Make/O/N=(DimSize(st, 0), DimSize(st, 1)) mask
	mask = 0
	mask[xmin, xmax][ymin,ymax] = 1
	
	variable i
	for(i=0; i<DimSize(st, 2); i+=1)
		imagetransform/p=(i) getplane st
		wave im = $"M_ImagePlane"
		Redimension/S im
		PowerSpectrum2D(im)
		wave im_ps2d = $"im_ps2d"
		im_ps2d *= mask
		stack_ps[i] = sum(im_ps2d)
	
	endfor
	
	Killwaves im_ps2d, M_ImagePlane, mask

end	 

 function DeleteLowPS(st, xmin, xmax, ymin, ymax, thresh)
 	wave st
 	variable xmin, xmax, ymin, ymax, thresh
 	
 	Duplicate/O st, lowps_del
 	
	Make/O/N=(DimSize(st, 0), DimSize(st, 1)) mask
	mask = 0
	mask[xmin, xmax][ymin,ymax] = 1
	
	variable i=0, j=0, diff
	do
		imagetransform/p=(i) getplane lowps_del
		wave im = $"M_ImagePlane"
		Redimension/S im
		PowerSpectrum2D(im)
		wave im_ps2d = $"im_ps2d"
		im_ps2d *= mask
		diff = sum(im_ps2d)

		if(diff <= thresh)
//			Imagetransform/O/P=(i) removeZplane dupe_del
			Deletepoints/M=2 i, 1, lowps_del
			j+=1
			printf "Deleted image %d\r" i
		else
			i+=1
		endif	
	
	while(i<DimSize(lowps_del, 2))
	
	printf "Deleted %d  images from the stack.\r", j
	
	Killwaves im_ps2d, M_ImagePlane, mask

end