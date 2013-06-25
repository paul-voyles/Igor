#pragma rtGlobals=1		// Use modern global access method.

// Function to convolve an image with a two-dimensional Guassian
// v1 04-20-12 pmv

// takes an image im and convolves it with a Gaussian source
// function with FWHM ds.  Uses the wave scaling of im.  Preserves
// the total intensity of the image (sum of all pixels) after scaling.
function SourceSizeConvolve(im, ds)
	wave im
	variable ds
	
	variable im_norm = sum(im)
	Duplicate/O im $(NameofWave(im)+"_ss")
	wave im_ss = $(NameofWave(im)+"_ss")
	variable fds = ds / (2*sqrt(2*ln(2)))	// real-space standard deviation for Gaussian with FWHM ds
	fds = 1/(2*Pi*fds)	// FT standard deviation
	//printf "fds = %g\r", fds
	Redimension/C im_ss
	wave/C res = $(NameofWave(im)+"_ss")
	FFT res
	res *=  cmplx(Gauss(x, 0.0, fds, y, 0.0, fds), 0.0)
	IFFT/C res
	Redimension/R res
	variable ss_norm = sum(res)
	res *= ss_norm / im_norm
	SetScale/p x dimoffset(im, 0), dimdelta(im, 0), "", res
	SetScale/P y dimoffset(im, 1), dimdelta(im, 1), "", res

end