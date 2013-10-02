#pragma rtGlobals=1		// Use modern global access method.
#include "Kirkland Program Input Generators"

// functions to generate a set of condor inputs for autostem and autoconfocal:
// take a full MXN output image and break it up into parts to be 
// parcelled out a few lines at a time.
// 03-07-10 reassemble bug fix pmv
// 03-09-10 added parameters for email address and exectuables.
// 01-03-11 added configuration for output of STEM images at multiple thicknesses
// 01-09-11 added support for autopacbed program
// fixed nw=0 bug in reassemble macro writing in StemChopInputsAndReassemble
// 12-31-11 fixed another thickness-related bug in StemChopInputsAndReassemble
// 03-29-12 changed reassemble macro to read the number of pixels i nthe pacbed pattern
//               from the data, not assume it.  Now autopacbed can set that flexibly and 
//               reassemble will still work
// 10-02-13 updated condor cmd file for current (2013) version of Condor.  Changed Image_size
//               to Requested_memory, and upped it by 50% to prevent jobs from being kicked 
//               off nodes for requesting too little memory.

function MakeControlWaves()
	
	string curFol = GetDataFolder(1)
	NewDataFolder/O root:Packages
	NewDataFolder/O/S root:Packages:stem_chop

	// Make wave containing control information
	Make/O/N=3 stem_p
	Make/T/O/N=3 stem_p_labels = {"kV", "condenser inner angle", "condenser outer angle"}
	stem_p = {200.0, 0.0, 17.5}
	Make/O/N=(12, 2) aber, aber_input, aber_default
	Make/O/T/N=(12, 4) aber_text
	aber_default[][0] = {0, 0, 22.56, 22.08, 0.1198, 0.9018, 0.04964, 28.43, 11.84, 8.456, 0.622, 2.811}
	aber_default[][1] = {0, 0, -20.1, -7.5, 0, -170.1, 20.9, -120.6, 153.8, 76.1, 0, -125.5} 
	aber_input = aber_default
	aber = aber_default
	aber_text[][0] = {"C1", "A1", "A2", "B2", "C3", "A3", "S3", "A4", "D4", "B4", "C5", "A5"}
	aber_text[][1] = num2str(aber_default[p][0])
	aber_text[][2] = num2str(aber_default[p][1])
	aber_text[][3] = {"nm", "nm", "nm", "nm", "um", "um", "um", "um", "um", "um", "mm", "mm"}
	SetDimLabel 1, 0, Aberration, aber_text
	SetDimLabel 1, 1, Coefficient, aber_text
	SetDimLabel 1, 2, Angle, aber_text
	SetDimLabel 1, 3, Units, aber_text
	Make/U/B/O/N=(12,4) aber_selwave
	aber_selwave[][0] = 0
	aber_selwave[][1] = 2
	aber_selwave[][2] = 2
	aber_selwave[][3] = 0
	Make/O/N=9 sim_p
	Make/O/T/N=9 sim_p_labels = {"transmission nx", "transmission ny", "probe nx", "probe ny", "slice thickness", "temperature", "# of phonons", "source size", "memory"}
	Make/O/T/N=1 detect_name = {"detector names"}
	make/O/N=(1,2) detect_p
	Make/O/N=1 thick_p
	SetDimLabel 0, 0, Thickness, thick_p
	make/O/N=9 imageout_p
	Make/O/T/N=9 imageout_p_labels = {"xi", "xf", "yi", "yf", "nx", "ny", "# of chunks in x", "# of chunks in y", "thickness level y/n"}
	make/O/T/N=4 names
	make/O/T/N=4 filename_labels = {"model filename", "output basename", "stem cmd comment", "model file path"}
	make/O/T/N=4 sim_paths = { "/home/voyles/bin/autostem", "voyles@engr.wisc.edu", "/filespace/people/v/voyles/bin/autostem", "/filespace/people/v/voyles/simulations/working/"}
	make/O/T/N=4 sim_paths_labels = {"cluster executable", "email", "condor executable", "condor working directory"}
	
	// Make default detector waves.
	Make/O/T/N=(12) titan_stem_detect_name = { "cl38mm", "cl48mm", "cl60mm", "cl77mm", "cl102mm", "cl128mm", "cl160mm", "cl196mm", "cl245mm", "cl301mm", "cl378mm", "cl478mm" }
	Make/O/N=(12, 2) titan_stem_detect_p = { {226, 179, 143, 112, 84.42, 67.3, 53.9, 44.0, 35.3, 28.8, 23.0, 18.2 }, {1130, 896, 717, 550, 422, 337, 270, 220, 176, 144, 115, 91.0 } }
	Make/O/T/N=(14) titan_efstem_detect_name = { "cl48mm", "cl68mm", "cl91mm", "cl115mm", "cl150mm", "cl194mm", "cl248mm", "cl313mm", "cl405mm", "cl512mm", "cl651mm", "cl837mm", "cl1037mm", "cl1317mm" }
	Make/O/N=(14, 2) titan_efstem_detect_p = { {534, 377, 282, 223, 171, 132, 103, 81.9, 63.3, 50.1, 39.4, 30.6, 24.7, 19.5 }, { 1180, 830, 620, 491, 376, 291, 228, 180, 139, 110, 86.7, 67.4, 54.4, 42.8 } }
	Make/O/T/N=(13) standard_stem_detect_name = {"BF", "d15_30mr", "d30_50mr", "d30_100mr", "d50_250mr", "d75_250mr", "d100_250mr", "d50_500mr", "d75_500mr", "d100_500mr", "d125_500mr", "d150_500mr", "d175_500mr"}
	Make/O/N=(13, 2) standard_stem_detect_p = { {0, 15, 30, 30, 50, 75, 100, 50, 75, 100, 125, 150, 175}, {1.5, 30, 50, 100, 250, 250, 250, 500, 500, 500, 500, 500, 500} }
	
	// Make memory estimate global variable
	variable/G mem_estimate, xmainpix, ymainpix, xpixleft, ypixleft, xspatial, yspatial, njobs, thick_interval, thick_end
	SetFormula xmainpix, "floor(imageout_p[4] / imageout_p[6])"
	SetFormula ymainpix, "floor(imageout_p[5] / imageout_p[7])"
	SetFormula xpixleft, "mod(imageout_p[4], imageout_p[6])"
	SetFormula ypixleft, "mod(imageout_p[5], imageout_p[7])"
	SetFormula mem_estimate, "(8*sim_p[0]*sim_p[1] + 8*sim_p[2]*sim_p[3]*ymainpix)/(1024^2) + 8"
	
	SetFormula xspatial, "(imageout_p[1] - imageout_p[0]) / (imageout_p[4] - 1)"
	SetFormula yspatial, "(imageout_p[3] - imageout_p[2]) / (imageout_p[5] - 1)"
	
	SetFormula njobs, "(imageout_p[6]*imageout_p[7]) + (xpixleft > 0 ? imageout_p[7] : 0) + (ypixleft > 0 ? imageout_p[6] : 0) + ( (xpixleft > 0) && (ypixleft > 0) ? 1 : 0)"
			
	SetDataFolder curFol
	
//	Edit stem_p_labels, stem_p
//	Edit imageout_p_labels, imageout_p
//	Edit sim_p_labels, sim_p
//	Edit detect_name, detect_p
//	Edit filename_labels, names
	
end

function StemChop(target, program)
	variable target  // 1 = cluster, 2 = condor
	string program

	NewPath/Q/O/M="Chose an output directory" this_path
	Pathinfo this_path
	
	variable outnum
	
	wave/T names = $"root:Packages:stem_chop:names"
	wave stem_p = $"root:Packages:stem_chop:stem_p"
	wave sim_p = $"root:Packages:stem_chop:sim_p"
	wave detect_p = $"root:Packages:stem_chop:detect_p"
	wave/T detect_name = $"root:Packages:stem_chop:detect_name"
	wave thick_p = $"root:Packages:stem_chop:thick_p"
	wave imageout_p = $"root:Packages:stem_chop:imageout_p"
	wave/T sim_paths = $"root:Packages:stem_chop:sim_paths"
	outnum = StemChopInputsandReassemble(program, S_path, names[1],  names[0], stem_p, aber, sim_p, detect_p, detect_name, imageout_p, thick_p)

	if(!outnum)
		printf "Error writing the input or reassembly file.  Exiting.\r"
		return 0
	endif
	StemChopCmd(target, sim_paths, S_path, names[1], names[2], sim_p[8], outnum)
	StemChopDescription(program, S_path, names, stem_p, aber, sim_p, detect_p, detect_name, imageout_p, thick_p, outnum)
	SaveControlWaves(S_path)

end


function SaveControlWaves(directory)
	string directory
	
	wave/T names = $"root:Packages:stem_chop:names"
	wave stem_p = $"root:Packages:stem_chop:stem_p"
	wave aber = $"root:Packages:stem_chop:aber"
	wave sim_p = $"root:Packages:stem_chop:sim_p"
	wave detect_p = $"root:Packages:stem_chop:detect_p"
	wave/T detect_name = $"root:Packages:stem_chop:detect_name"
	wave thick_p = $"root:Packages:stem_chop:thick_p"
	wave imageout_p = $"root:Packages:stem_chop:imageout_p"

	directory += (names[1] + "_param.awav")

	Save/T names, stem_p, aber, sim_p, detect_p, detect_name, thick_p, imageout_p as directory
	
end

function StemChopDescription(program, directory, names, stem_p, aber, sim_p, detect_p, detect_name, imageout_p, thick_p, outnum)
	string program, directory
	wave/T names
	wave stem_p, aber, sim_p, detect_p
	wave/T detect_name
	wave imageout_p
	wave thick_p
	variable outnum
	
	variable f
	string desc_name
	sprintf desc_name, "%s%s.txt", directory, names[1]
	open/Z f as desc_name
	if(V_flag)
		printf "Cannot open description file.  Exiting.\r"
		return 0
	endif
	
	fprintf f, "Stem chop description generated %s %s.\r\n\r\n", time(), date() 
	
	fprintf f, "Program %s operating on model %s.\r\n", program, names[0]
	fprintf f, "%s\r\n", names[2]
	fprintf f, "Accelerating voltage = %g kV.\r\n\r\n", stem_p[0]
	
	fprintf f, "Probe lens parameters:\r\n"
	fprintf f, "inner angle = %g mrad.\r\n", stem_p[1]
	fprintf f, "Outer angle = %g mrad.\r\n\r\n", stem_p[2]
	fprintf f, "C1 = %g nm; A1 = %g nm, %g deg; A2 = %g nm, %g deg\r\n", aber[0][0], aber[1][0], aber[1][1], aber[2][0], aber[2][1]
	fprintf f, "B2 = %g nm, %g deg; C3 = %g um; A3 = %g um, %g deg\r\n", aber[3][0], aber[3][1], aber[4][0], aber[5][0], aber[5][1]
	fprintf f, "S3 = %g, um, %g deg; A4 = %g um; %g deg; D4 = %g um, %g deg\r\n", aber[6][0], aber[6][1], aber[7][0], aber[7][1], aber[8][0], aber[8][1]
	fprintf f, "B4 = %g um, %g deg; C5 = %g mm; A5 = %g mm, %g deg\r\n\r\n", aber[9][0], aber[9][1], aber[10][0], aber[11][0], aber[11][1]
		
	fprintf f, "Simulation parameters:\r\n"
	fprintf f, "Specimen transmission function: %d x %d pixels.\r\n", sim_p[0], sim_p[1]
	fprintf f, "Probe wave functions: %d x %d pixels.\r\n", sim_p[2], sim_p[3]
	fprintf f, "Slice thicknes = %g A.\r\n", sim_p[4]
	if(sim_p[5] >= 0)
		fprintf f, "Averaged over %d phonon configurations at a temperature of %g K.\r\n", sim_p[6], sim_p[5]
	else
		fprintf f, "No phonon disorder added.\r\n"
	endif
	fprintf f, "Source size = %g A.\r\n", sim_p[7]
	fprintf f, "Expected memory use = %g MB.\r\n\r\n", sim_p[8]
	
	fprintf f, "Output image parameters:\r\n"
	fprintf f, "Output image runs from x = %g to %g and y = %g to %g A.\r\n", imageout_p[0], imageout_p[1], imageout_p[2], imageout_p[3]
	fprintf f, "Image consists of %d x %d pixels.\r\n", imageout_p[4], imageout_p[5]
	if(imageout_p[8])
		fprintf f, "Output saved at multiple thickness:\r\n"
		variable nthick = numpnts(thick_p)
		variable nt
		for(nt=0; nt<nthick; nt+=1)
			fprintf f, "%g A\r\n", thick_p[nt]
		endfor
		fprintf f, "and the final model thickness.\r\n\r\n"
	else
		fprintf f, "Output saved only at final model thickness.\r\n"
	endif
	fprintf f, "Image is generated in %d x chunks and %d y chunks.\r\n", imageout_p[6], imageout_p[7]
	fprintf f, "This results in %d total subimages.\r\n\r\n", outnum
	
	fprintf f, "Detector parameters:\r\n"
	variable ndetect = numpnts(detect_name)
	variable i
	for(i=0; i<ndetect; i+=1)
		fprintf f, "Detector %s extends from %g to %g mrad.\r\n ", detect_name[i], detect_p[i][0], detect_p[i][1]
//		if(!cmpstr(program, "autostem"))
//			fprintf f, "mrad.\r\n"
//		else
//			fprintf f, "A.\r\n"
//		endif
	endfor
	

	close f
end


// image_p wave, N points: xi, xf, yi, yf, nx, ny, ncx, ncy

function StemChopInputsAndReassemble(program, directory, basename, modelname, stem_p, aber, sim_p, detect_p, detect_name, imageout_p, thick_p)
	string program
	string directory, basename, modelname
	wave stem_p, aber, sim_p, detect_p
	wave/t detect_name
	wave imageout_p
	wave thick_p
	
	variable xi = imageout_p[0]
	variable xf = imageout_p[1]
	variable yi = imageout_p[2]
	variable yf = imageout_p[3]
	variable nx = imageout_p[4]
	variable ny = imageout_p[5]
	variable nchunksx = imageout_p[6]
	variable nchunksy = imageout_p[7]
	variable thick_yn = imageout_p[8]
	
	//variable pnx = sim_p[2]
	//variable pny = sim_p[3]
	string imname

	FUNCREF ProtoOutput outfunc = OneAutostemImageOut

	if(nchunksx > nx || nchunksy > ny)
		printf "Don't be stupid!  Have at least one pixel per simulation.\r"
		return 0
	endif
	
	variable pacbed = !cmpstr(program, "autopacbed")
	
	variable f
	string reassemble
	sprintf reassemble, "%s%s_reassemble.ipf", directory, basename
	Open/Z f as reassemble
	if(V_Flag)
		printf "Cannot open reassembly ipf file.  Exiting.\r"
		return 0
	endif
	fprintf f, "Macro reassemble()\r"

	fprintf f, "NewPath/O/M=\"Select the directory of gfx files\" gfxp\r"
	fprintf f, "PathInfo gfxp\r"
	fprintf f, "if (!V_flag)\r"
	fprintf f, "return 0\r"
	fprintf f, "endif\r"
	
	// read # of pacbed pixels from first stored pattern
	if(pacbed)
		fprintf f, "variable npx, npy\r"
		sprintf imname, "%s_im0_pacbed_0" basename
		fprintf f, "LoadGFX(S_Path+\"%s.gfx\")\r", imname
		fprintf f, "npx = DimSize(gfx_read, 0)\r"
		fprintf f, "npy = DimSize(gfx_read, 1)\r"
	endif		

	variable ndetect = numpnts(detect_name)
	variable nthick = numpnts(thick_p)+1
	if(!thick_yn)
		nthick = 0
		make/o/n=0 thick_out_p
	else
		duplicate/o thick_p thick_out_p
	endif
	variable nd, nt, nw = 0
	if(thick_yn)
		Make/T/O/N=(ndetect*(nthick)) wname
	else
		Make/T/O/N=(ndetect) wname
	endif
	Make/T/O/N=(nthick) pname
	string onen
	for(nd=0; nd<ndetect; nd+=1)
		if(thick_yn)
			for(nt = 0; nt<nthick; nt += 1)
				if(pacbed && nd==0)
					sprintf onen, "%s_pacbed_t%d", basename, nt+1
					pname[nt] = onen
					fprintf f, "Make/O/N=(npx, npy) %s\r",  onen
				endif
				sprintf onen, "%s_%s_t%d", basename, detect_name[nd], nt+1
				wname[nw] = onen
				fprintf f, "Make/O/N=(%d, %d) %s\r", nx, ny, wname[nw]
				fprintf f, "SetScale/I x %f, %f, \"\", %s\r", xi, xf, wname[nw]
				fprintf f, "SetScale/I y %f, %f, \"\", %s\r", yi, yf, wname[nw]
				nw+=1				
			endfor
		else
			if(pacbed)
				sprintf onen, "%s_pacbed", basename
				pname[0] = onen
				fprintf f, "Make/O/N=(npx, npy) %s\r",  onen
			endif
			sprintf onen, "%s_%s", basename, detect_name[nd]
			wname[nw] = onen
			fprintf f, "Make/O/N=(%d, %d) %s\r", nx, ny, wname[nw]
			fprintf f, "SetScale/I x %f, %f, \"\", %s\r", xi, xf, wname[nw]
			fprintf f, "SetScale/I y %f, %f, \"\", %s\r", yi, yf, wname[nw]
			nw+=1
		endif
	endfor

	
	variable xic, xfc, yic, yfc, nxc, nyc
	variable xstep, ystep
	xstep = (xf - xi) / (nx -1)
	ystep = (yf - yi) / (ny -1)
	printf "Pixel size is %f x %f Angstroms.\r", xstep, ystep

	// can only do an integral number of pixels in each chunk.  find out how many pixels are
	// left over once we've divided the nx, ny pixels into nchunksx, nchunksy chunks.
	variable npx_leftover = mod(nx, nchunksx)
	variable npy_leftover = mod(ny, nchunksy)
	
	// number of pixels per chunk
	variable npx = floor(nx / nchunksx)
	variable npy = floor(ny / nchunksy)

	printf "%d by %d chunks of %d by %d pixels each, with %d pixels left over in x and %d pixels left over in y.\r", nchunksx, nchunksy, npx, npy, npx_leftover, npy_leftover
	

	Make/O/N=6 chop_image_p
	variable outnum = 0


	//Generate files for the part of the image that fits in the chunks
	chop_image_p[4] = npx
	chop_image_p[5] = npy

	variable ii, jj
	for(ii=0; ii<nchunksx; ii+=1 )
		for(jj=0; jj<nchunksy; jj+=1)
			chop_image_p[0] = xi + xstep*ii*npx
			chop_image_p[1] = chop_image_p[0] + xstep*(npx-1)
			chop_image_p[2] = yi + ystep*jj*npy
			chop_image_p[3] = chop_image_p[2] + ystep*(npy-1)
			outfunc(directory, basename, modelname, stem_p, aber, sim_p, chop_image_p, detect_p, detect_name, thick_out_p, outnum)
			nw = 0
			for(nd=0; nd<ndetect; nd+=1)
				if(thick_yn)
					for(nt=0; nt<nthick; nt+=1)
						if(pacbed && nd==0)
							sprintf imname, "%s_im%d_pacbed_%d" basename, outnum, nt
							fprintf f, "LoadGFX(S_Path+\"%s.gfx\")\r", imname
							fprintf f, "%s += gfx_read\r", pname[nt]
						endif
						sprintf imname, "%s_im%d_%.3d%.3d", basename, outnum, nd, nt
						fprintf f, "LoadGFX(S_Path+\"%s.gfx\")\r", imname
						fprintf f, "%s[%d, %d][%d, %d] =gfx_read[p-%d][%d - q]\r", wname[nw], ii*npx, (ii+1)*npx - 1, jj*npy, (jj+1)*npy -1, ii*npx, (jj+1)*npy -1	
						nw+=1			
					endfor
				else
					if(pacbed)
						sprintf imname "%s_im%d_pacbed", basename, outnum
						fprintf f, "LoadGFX(S_Path+\"%s.gfx\")\r", imname
						fprintf f, "%s += gfx_read\r", pname[0]
					endif
					sprintf imname, "%s_im%d_%.3d000", basename, outnum, nd
					fprintf f, "LoadGFX(S_Path+\"%s.gfx\")\r", imname
					fprintf f, "%s[%d, %d][%d, %d] =gfx_read[p-%d][%d - q]\r", wname[nw], ii*npx, (ii+1)*npx - 1, jj*npy, (jj+1)*npy -1, ii*npx, (jj+1)*npy -1
					nw+=1
				endif
			endfor
			outnum += 1
		endfor
	endfor
	
	//Generate files for the extra x chunk along all the existing y chunks
	if(npx_leftover)
		chop_image_p[4] = npx_leftover
		chop_image_p[5] = npy
		chop_image_p[0] = xi + xstep*nchunksx*npx
		chop_image_p[1] = chop_image_p[0] + xstep*(npx_leftover-1)
		for(jj=0; jj<nchunksy; jj+=1)
			chop_image_p[2] = yi + ystep*jj*npy
			chop_image_p[3] = chop_image_p[2] + ystep*(npy-1)
			outfunc(directory, basename, modelname, stem_p, aber, sim_p, chop_image_p, detect_p, detect_name, thick_out_p, outnum)
			nw = 0 
			for(nd=0; nd<ndetect; nd+=1)
				if(thick_yn)
					for(nt=0; nt<nthick; nt+=1)
						if(pacbed && nd==0)
							sprintf imname, "%s_im%d_pacbed_%d" basename, outnum, nt
							fprintf f, "LoadGFX(S_Path+\"%s.gfx\")\r", imname
							fprintf f, "%s += gfx_read\r", pname[nt]
						endif
						sprintf imname, "%s_im%d_%.3d%.3d", basename, outnum, nd, nt
						fprintf f, "LoadGFX(S_Path+\"%s.gfx\")\r", imname
						fprintf f, "%s[%d, %d][%d, %d] = gfx_read[p-%d][%d -q]\r", wname[nw], nchunksx*npx, nchunksx*npx +npx_leftover- 1, jj*npy, (jj+1)*npy-1, nchunksx*npx, (jj+1)*npy-1
						nw+=1				
					endfor
				else
					if(pacbed)
						sprintf imname "%s_im%d_pacbed", basename, outnum
						fprintf f, "LoadGFX(S_Path+\"%s.gfx\")\r", imname
						fprintf f, "%s += gfx_read\r", pname[0]
					endif
					sprintf imname, "%s_im%d_%.3d000", basename, outnum, nd
					fprintf f, "LoadGFX(S_Path+\"%s.gfx\")\r", imname
					fprintf f, "%s[%d, %d][%d, %d] = gfx_read[p-%d][%d -q]\r", wname[nw], nchunksx*npx, nchunksx*npx +npx_leftover- 1, jj*npy, (jj+1)*npy-1, nchunksx*npx, (jj+1)*npy-1
					nw+=1
				endif
			endfor
			outnum += 1
		endfor
	endif	
	
	//Generate files for the extra y chunk along all the existing x chunks
	if(npy_leftover)
		chop_image_p[4] = npx
		chop_image_p[5] = npy_leftover
		chop_image_p[2] = yi + ystep*nchunksy*npy
		chop_image_p[3] = chop_image_p[2] + ystep*(npy_leftover-1)
		for(ii=0; ii<nchunksx; ii+=1)
			chop_image_p[0] = xi + xstep*ii*npx
			chop_image_p[1] = chop_image_p[0] + xstep*(npx-1)
			outfunc(directory, basename, modelname, stem_p, aber, sim_p, chop_image_p, detect_p, detect_name, thick_out_p, outnum)
			nw=0
			for(nd=0; nd<ndetect; nd+=1)
				if(thick_yn)
					for(nt=0; nt<nthick; nt+=1)
						if(pacbed && nd==0)
							sprintf imname, "%s_im%d_pacbed_%d" basename, outnum, nt
							fprintf f, "LoadGFX(S_Path+\"%s.gfx\")\r", imname
							fprintf f, "%s += gfx_read\r", pname[nt]
						endif
						sprintf imname, "%s_im%d_%.3d%.3d", basename, outnum, nd, nt
						fprintf f, "LoadGFX(S_Path+\"%s.gfx\")\r", imname
						fprintf f, "%s[%d, %d][%d, %d] = gfx_read[p-%d][%d-q]\r", wname[nw], ii*npx, (ii+1)*npx - 1, nchunksy*npy, nchunksy*npy +npy_leftover-1, ii*npx, nchunksy*npy +npy_leftover-1
						nw+=1
					endfor					
				else
					if(pacbed)
						sprintf imname "%s_im%d_pacbed", basename, outnum
						fprintf f, "LoadGFX(S_Path+\"%s.gfx\")\r", imname
						fprintf f, "%s += gfx_read\r", pname[0]
					endif
					sprintf imname, "%s_im%d_%.3d000", basename, outnum, nd
					fprintf f, "LoadGFX(S_Path+\"%s.gfx\")\r", imname
					fprintf f, "%s[%d, %d][%d, %d] = gfx_read[p-%d][%d-q]\r", wname[nw], ii*npx, (ii+1)*npx - 1, nchunksy*npy, nchunksy*npy +npy_leftover-1, ii*npx, nchunksy*npy +npy_leftover-1
					nw+=1
				endif
			endfor
			outnum += 1
		endfor
	endif	

	//Generate the corner extra chunk: extra in x and y
	if(npx_leftover && npy_leftover)
		chop_image_p[4] = npx_leftover
		chop_image_p[5] = npy_leftover
		chop_image_p[0] = xi + xstep*nchunksx*npx
		chop_image_p[1] = chop_image_p[0] + xstep*(npx_leftover-1)
		chop_image_p[2] = yi + ystep*nchunksy*npy
		chop_image_p[3] = chop_image_p[2] + ystep*(npy_leftover-1)
		outfunc(directory, basename, modelname, stem_p, aber, sim_p, chop_image_p, detect_p, detect_name, thick_out_p, outnum)
		nw=0
		for(nd=0; nd<ndetect; nd+=1)
			if(thick_yn)
				for(nt=0; nt<nthick; nt+=1)
					if(pacbed && nd==0)
						sprintf imname, "%s_im%d_pacbed_%d" basename, outnum, nt
						fprintf f, "LoadGFX(S_Path+\"%s.gfx\")\r", imname
						fprintf f, "%s += gfx_read\r", pname[nt]
					endif
					sprintf imname, "%s_im%d_%.3d%.3d", basename, outnum, nd, nt
					fprintf f, "LoadGFX(S_Path+\"%s.gfx\")\r", imname
					fprintf f, "%s[%d, %d][%d, %d] = gfx_read[p-%d][%d - q]\r", wname[nw], nchunksx*npx, nchunksx*npx +npx_leftover- 1, nchunksy*npy, nchunksy*npy +npy_leftover-1, nchunksx*npx, nchunksy*npy +npy_leftover-1
					nw+=1
				endfor				
			else
				if(pacbed)
					sprintf imname "%s_im%d_pacbed", basename, outnum
					fprintf f, "LoadGFX(S_Path+\"%s.gfx\")\r", imname
					fprintf f, "%s += gfx_read\r", pname[0]
				endif
				sprintf imname, "%s_im%d_%.3d000", basename, outnum, nd
				fprintf f, "LoadGFX(S_Path+\"%s.gfx\")\r", imname
				fprintf f, "%s[%d, %d][%d, %d] = gfx_read[p-%d][%d - q]\r", wname[nw], nchunksx*npx, nchunksx*npx +npx_leftover- 1, nchunksy*npy, nchunksy*npy +npy_leftover-1, nchunksx*npx, nchunksy*npy +npy_leftover-1
				nw+=1
			endif
		endfor
		outnum += 1
	endif

	if(pacbed)
		if(thick_yn)
			for(nt=0; nt<nthick; nt+=1)
				fprintf f, "%s /= %d\r", pname[nt], outnum
			endfor
		else
			fprintf f, "%s /= %d\r", pname[0], outnum
		endif
	endif

	fprintf f, "Killwaves gfx_read\r"
	fprintf f, "end\r"
	close f

	printf "Wrote %d %s input files to directory %s.\r", outnum, program, directory
	Killwaves chop_image_p, wname, pname, thick_out_p
	return outnum
	
end


function StemChopCmd(target, sim_paths, directory, basename, comment, size, outnum)
	variable target  // 1 for cluster, 2 for condor
	wave/t sim_paths 
	string directory, basename, comment
	variable outnum, size
	
	variable f, g, nim
	string chopname

	// cluster chop
	if(target == 1)
		sprintf chopname, "%s%s.sh", directory, basename
		Open/Z f as chopname
		if(V_flag)
			printf "Cannot open qsub .sh file.  Exiting.\r"
			return 0
		endif

		fprintf f, "#!/bin/bash\n"
		
		for(nim=0; nim < outnum; nim+=1)
			sprintf chopname, "%s%s_im%d.sh", directory, basename, nim
			Open/Z g as chopname
			if(V_flag)
				printf "Cannot open image %d .sh file.  Exiting.\r", nim
				return 0
			endif
			
			fprintf f, "qsub %s_im%d.sh\n", basename, nim
			
			fprintf g, "#!/bin/bash\n"
			fprintf g, "#\n"
			fprintf g, "#$ -cwd\n"
			fprintf g, "#$ -m be\n"
			fprintf g, "#$ -M %s\n", sim_paths[1]
			fprintf g, "#$ -S /bin/bash\n"
			fprintf g, "#$ -e %s_im%d.err\n", basename, nim
			fprintf g, "#$ -i %s_im%d.input\n", basename, nim
			fprintf g, "#$ -o %s_im%d.out\n", basename, nim
			fprintf g, "#\n"
			fprintf g, "%s\n", sim_paths[0]
			close g
		endfor

	// Condor chop
	elseif (target == 2)
	
		sprintf chopname, "%s%s.cmd", directory, basename
		Open/Z f as chopname
		if(V_flag)
			printf "Cannot open condor .cmd file.  Exiting.\r"
			return 0
		endif
	
		string s
	
		fprintf f, "#%s\n\n", comment
	
		fprintf f, "#General Condor options\n"
		fprintf f, "Universe = standard\n"
		fprintf f, "\n"
	
		fprintf f, "#Job description\n"
		fprintf f, "Executable = %s\n", sim_paths[2]
		fprintf f, "Request_memory = %d \n", 1.5*size
		fprintf f, "\n"
		
		fprintf f, "#File locations\n"
		variable i = ItemsInList(directory, ":")
		directory = StringFromList((i-1), directory, ":")
		sprintf s, "%s%s", sim_paths[3], directory
		fprintf f, "Initialdir = %s\n", s
		sprintf s, "%s_im$(Process).input", basename
		fprintf f, "Input = %s\n", s
		sprintf s, "%s_im$(Process).out", basename
		fprintf f, "Output = %s\n", s
		sprintf s, "%s_im$(Process).err", basename
		fprintf f, "Error = %s\n", s
		sprintf s, "%s_im$(Process).log", basename
		fprintf f, "Log = %s\n", s
		fprintf f, "\n"
		
		fprintf f, "#Number of jobs from this cmd file:\n"
		fprintf f, "Queue %d\n", outnum
	endif
		
	close f

	
end