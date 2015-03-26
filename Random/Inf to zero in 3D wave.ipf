#pragma rtGlobals=1		// Use modern global access method.


function InfToZero(im_series)
	wave im_series
	variable x_size = DimSize(im_series, 0)
	variable y_size = DimSize(im_series, 1)
	variable z_size = DimSize(im_series, 2)
	
	variable value
	variable counter = 0
	
	variable i=0
	variable j=0
	variable k=0
	for(i=0; i<x_size; i+=1)
		j=0
		for(j=0; j<y_size; j+=1)
			k=0
			for(k=0; k<z_size; k+=1)

				value = im_series[i][j][k]
				if(value == inf)
					im_series[i][j][k] = 0
					counter = counter +1	
				endif
				
			endfor
		endfor
	endfor
	print "number of inf changes to zero = ", counter
end

function EraseFaultyImages(im_series)
	wave im_series
	
	variable z_size = DimSize(im_series, 2)
	
	variable value
	variable counter = 0
	
	variable i=0
	for(i=0; i<z_size; i+=1)

		imagetransform/p=(i-counter) getPlane im_series
		WaveStats/Q M_ImagePlane
		
		if(V_max > 1000)
			imagetransform/o/p=(i-counter) removeZplane im_series
			counter = counter +1	
		endif
		
		killwaves M_ImagePlane
	endfor
	print "number of slices removed = ", counter
end