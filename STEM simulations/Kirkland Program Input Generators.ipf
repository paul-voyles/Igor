#pragma rtGlobals=1		// Use modern global access method.

// added C5 to stem_p for probe forming and image forming 7/5/04.  
// updated to 64-bit, Nov 2008 version of autostem, 01-05-10.
// autoconfocal removed 01/05/10
// random number seed input added 02/24/10
// 01-03-11 added ability to generate output at multiple thicknesses using new version of autostem
// 01-09-11 added support for autopacbed program
// 07-07-11 fixed bug in OneAutoStemImageOut that printed thickness levels as integers instead of fp

// stem_p wave 5 points: kV, Cs, df, apeture inner angle, aperture outer angle
// image_p wave 6 points:  xi, xf, yi, yf, nx, ny
// sim_p wave  8 points: tf nx, tf ny, wf nx, wf ny, slice thickness, temperature (< 0 no phonons), # of phonon configurations, source size
// detect_angles, 2D wave, inner and outer angles for each detector
// detect_name, output image basename for each detector.
// outnum is the number in the total sequence, starting at 0.
// Initializing Cc and dE to 0.0 and 0.0, respectively for autostem calculation on Condor. This will not influence the running on Odie.  07/10/2014 JFeng
// 06-30-17 adapted to choose version between c and cpp cz

function ProtoOutput(directory, basename, modelname, stem_p, aber, sim_p, image_p, detect_p, detect_name, thick_p, outnum, version, detector_p, sensmap)
	string directory, basename, modelname, sensmap
	wave stem_p, aber, sim_p, image_p, detect_p, detector_p
	wave/t detect_name
	wave thick_p
	variable outnum
	variable version

	printf "Prototype output function.  Should never be called.\r"

end	


function OneAutostemImageOut(directory, basename, modelname, stem_p, aber, sim_p, image_p, detect_p, detect_name, thick_p, outnum, version, detector_p, sensmap)
	string directory, basename, modelname, sensmap
	wave stem_p, aber, sim_p, image_p, detect_p, detector_p
	wave/t detect_name
	wave thick_p
	variable outnum
	variable version
	
	if(image_p[5] == 1.0 || image_p[6] == 1.0)
		printf "One-dimensional images must be generated with a line-scan function.  Exiting.\r"
		return 0
	endif
	
	string filename
	sprintf filename, "%s%s_im%d.input", directory, basename, outnum
	
	variable f
	open/Z f as filename
	if(V_flag)
		printf "There was an error opening the file %s.\r", filename
		return 0
	endif
	
	variable ndetect = numpnts(detect_name)
	
	if (detector_p[3])
		fprintf f, "%s\n", sensmap //to use detector sensitivity map
   		fprintf f, "%d %d\n", detector_p[0], detector_p[1] // detector sensitivity map center position in unit of pixels
		fprintf f, "0\n"   // rotation to the CBED pattern. Positive value corresponds to clockwise rotation
	endif
	fprintf f, "%s\n", modelname
	fprintf f, "%d %d %d\n", sim_p[11], sim_p[12],sim_p[13]	// replicate unit cell
	fprintf f, "%f   %f   %f \n", stem_p[0], stem_p[1], stem_p[2]
	fprintf f, "%f   %f   %f \n", aber[0][0], aber[1][0], aber[1][1]  // C1, A1
	fprintf f, "%f   %f   %f   %f \n", aber[2][0], aber[2][1], aber[3][0], aber[3][1]  // A2, B2
	fprintf f, "%f   %f   %f   %f   %f\n", aber[4][0], aber[5][0], aber[5][1], aber[6][0], aber[6][1]  // C3, A3, S3
	fprintf f, "%f   %f   %f   %f   %f   %f\n" aber[7][0], aber[7][1], aber[8][0], aber[8][1], aber[9][0], aber[9][1] // A4, D4, B4
	fprintf f, "%f   %f   %f\n", aber[10][0], aber[11][0], aber[11][1]  // C5, A5
	fprintf f, "%d   %d\n", sim_p[0], sim_p[1]
	fprintf f, "%d   %d\n", sim_p[2], sim_p[3]
	fprintf f, "%f   %f\n", sim_p[9], sim_p[10]	// crystal tilt
	fprintf f, "n\n"	// 2-D image, not line scan

	if(numpnts(thick_p))  // use thickness levels
		fprintf f, "%d\n", (numpnts(thick_p)+1)
		variable nt
		for(nt = 0; nt<numpnts(thick_p); nt+=1)
			fprintf f, "%g\n", thick_p[nt]
		endfor
	else
		fprintf f, "1\n"  
	endif

	fprintf f, "%s_im%d_\n", basename, outnum
	fprintf f, "%d\n", ndetect
	variable i
	for(i=0; i<ndetect; i+=1)
		if (version == 1)
			fprintf f, "%f   %f\n", detect_p[i][0], detect_p[i][1]
			// same distance is used for all settings, OK to use same detector for different CL, cannot use multiple detector (like HAADF and ABF)
			if (detector_p[3])	
				fprintf f, "%d\n", detector_p[2] //distance (in px) between detector center and its inner edge
				fprintf f, "sensmap_adapted%d.tif\n", i+1 // fftshift detector to adapt to non-fftshifted detector map
			endif
		else
			fprintf f, "%f   %f m\n", detect_p[i][0], detect_p[i][1]
		endif
	endfor
	
	fprintf f, "%f   %f   %f   %f   %d   %d\n", image_p[0], image_p[1], image_p[2], image_p[3], image_p[4], image_p[5]
	fprintf f, "%f\n", sim_p[4]
	
	// phonons
	if(sim_p[5] < 0.0)
		fprintf f, "n\n"
	else
		fprintf f, "y\n"
		fprintf f, "%f\n", sim_p[5]
		fprintf f, "%d\n", sim_p[6]
		if(version == 1)
			fprintf f, "%ld\n", (2^31-1)*(enoise(0.5)+0.5)	// random seed and source size, only for C version
			fprintf f, "%f\n", sim_p[7]
		endif
	endif
	
	if(version == 1)
		fprintf f, "0.0   0.0\n", 	// Cc and dE 07/10/2014 by Jie Feng
	else
		fprintf f, "0.0\n",		// zero source size for c++ version
	endif 
	
	//fprintf f, "version code is %d\n", version
	
	close f
end


// stem_p wave 5 points: kV, Cs, df, apeture inner angle, aperture outer angle
// image_p wave 5 points:  xi, xf, yi, yf, np
// sim_p wave  8 points: tf nx, tf ny, wf nx, wf ny, slice thickness, temperature (< 0 no phonons), # of phonon configurations, source size
// detect_p, 2D wave, inner and outer angles for each detector
// outname - name of the text file to hold the image output
// outnum is the number in the total sequence, starting at 0.

function OneAutostemLineOut(directory, basename, modelname, stem_p, aber, sim_p, image_p, detect_p, thick_p, outname, outnum)
	string directory, basename, modelname
	wave stem_p, aber, sim_p, image_p, detect_p, thick_p
	string outname
	variable outnum
	
	if(image_p[5] == 1.0 || image_p[6] == 1.0)
		printf "One-dimensional images must be generated with a line-scan function.  Exiting.\r"
		return 0
	endif
	
	string filename
	sprintf filename, "%s\%s_im%d.input", directory, basename, outnum
	
	variable f
	open/Z f as filename
	if(V_flag)
		printf "There was an error opening the file %s.\r", filename
		return 0
	endif
	
	variable ndetect = numpnts(detect_name)
	
	fprintf f, "%s\n", modelname
	fprintf f, "1 1 1\n"	// replicate unit cell
	fprintf f, "%f   %f   %f \n", stem_p[0], stem_p[1], stem_p[2]
	fprintf f, "%f   %f   %f \n", aber[0][0], aber[1][0], aber[1][1]  // C1, A1
	fprintf f, "%f   %f   %f   %f \n", aber[2][0], aber[2][1], aber[3][0], aber[3][1]  // A2, B2
	fprintf f, "%f   %f   %f   %f   %f\n", aber[4][0], aber[5][0], aber[5][1], aber[6][0], aber[6][1]  // C3, A3, S3
	fprintf f, "%f   %f   %f   %f   %f   %f\n" aber[7][0], aber[7][1], aber[8][0], aber[8][1], aber[9][0], aber[9][1] // A4, D4, B4
	fprintf f, "%f   %f   %f\n", aber[10][0], aber[11][0], aber[11][1]  // C5, A5
	fprintf f, "%d   %d\n", sim_p[0], sim_p[1]
	fprintf f, "%d   %d\n", sim_p[2], sim_p[3]
	fprintf f, "0.0   0.0\n", 	// crystal tilt
	fprintf f, "y\n"	// line scan, not 2-D image
	fprintf f, "%s_im%d.txt\n", outname, outnum
	
	fprintf f, "%d\n", ndetect
	variable i
	for(i=0; i<ndetect; i+=1)
		fprintf f, "%f   %f\n", detect_p[i][0], detect_p[i][1]
	endfor
	
	fprintf f, "%f   %f   %f   %f   %d   %d\n", image_p[0], image_p[1], image_p[2], image_p[3], image_p[4]
	fprintf f, "%f\n", sim_p[4]
	
	// phonons
	if(sim_p[5] < 0.0)
		fprintf f, "n\n"
	else
		fprintf f, "y\n"
		fprintf f, "%f\n", sim_p[5]
		fprintf f, "%d\n", sim_p[6]
		fprintf f, "%ld\n", (2^31-1)*(enoise(0.5)+0.5)
		fprintf f, "%f\n", sim_p[7]
	endif
	
	fprintf f, "0.0   0.0\n", 	// Cc and dE 07/10/2014 by Jie Feng
	
	close f
end

