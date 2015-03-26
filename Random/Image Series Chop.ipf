#pragma rtGlobals=1		// Use modern global access method.


function EvenChop(im)
	
	wave im
	
	variable size_x = DimSize(im,0)
	variable size_y = DimSize(im,1)
	variable size_z = DimSize(im,2)

	Make/o/u/w/n=(size_x, size_y, size_z/2) im_even_chop
	
	variable i
	for(i=0; i<size_z/2; i+=1)
		imagetransform/p=((i*2)+1) getplane im
		wave M_ImagePlane = $"M_ImagePlane"
		imagetransform /INSI=M_ImagePlane /INSX=0 /INSY=0 /P=(i) insertImage im_even_chop
	endfor
	killwaves M_ImagePLane
end

function OddChop(im)
	
	wave im
	
	variable size_x = DimSize(im,0)
	variable size_y = DimSize(im,1)
	variable size_z = DimSize(im,2)

	Make/o/u/w/n=(size_x, size_y, size_z/2) im_odd_chop
	
	variable i
	for(i=0; i<size_z/2; i+=1)
		imagetransform/p=(i*2) getplane im
		wave M_ImagePlane = $"M_ImagePlane"
		imagetransform /INSI=M_ImagePlane /INSX=0 /INSY=0 /P=(i) insertImage im_odd_chop
	endfor
	killwaves M_ImagePLane
end