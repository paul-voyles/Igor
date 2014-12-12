#pragma rtGlobals=3		// Use modern global access method and strict wave access.

// Functions to coarse-grain a lattice image by breaking into cells, each of which contains one
// column, then assigning each cell one intensity.  Originally developed to track black spot defects
// in LAADF images of SiC.
//
// Basic usage: (1) define a set lattice positions by e.g. fitting or setting a grid.  (2) Define VoronoiCells for
// those positions.  (3) MakeAllMasks from the Voronoi cells.  (4) Coarse-grain the image (image
// stack) using MakeCoarseGrain (MakeCoarseGrainStack) and the appropriate MaskProc.
//
// VoronoiCells: from a set of identified columns positions in the image, divide the image into cells
// using a Voronoi tesselation.  Identify the vertices that are associated with each cell.
//
// VoronoiVertices: take the output of the built-in Igor routine to calculate the edges of a Voronoi tessellation
// of a set of points, and extract the vertices of the Voronoi polygons that fit entirely inside the a defined area
// 
// MakeAllMasks: take the identified vertices and generate a stack of mask waves.  Each wave in the stack
// is the same size as the original image and is 1 inside a particular Voronoi polygon and one outside it.  Ideally,
// the union of all the masks would be the area of the original image, although it doesn't quite work that way.
//
// PolyEdgesfromVoronoi: not used.
//
// AllVertNearPoints / AllVertCircuits / VertCircuit / GenerateMasks: utility routines for MakeAllMasks
//
// MaskCoarseGrain: take an image and a stack of masks.  Perform some operation defined in the function reference
// MaskProc to calculate the intensity inside each mask, and collect all those intensities into a coarse-grain representation
// of the input image.  MaskAverage is an example MaskProc.
//
// MaskCoarseGrainStack: same process as MaskCoarseGrain, but on every image in a stack of input images.
//
// ProtoMaskProc: prototype of a MaskProc.  Take an input image and a mask image.  Every MaskProc must have the
// same arrguments as the prototype.
//
// MaskAverage: assign the coarse-grain region to the average value of the input image inside the each mask / Voronoi polygon
//
// MaskZ0: assign the coarse-grain region to the z0 background value from a fit of Gaussian peaks to the original image.
// The Gaussian peak fit parameters enter through hard-coded global variables and global waves.
//
// MaskGauss: assign the coarse-grain region to the average of the Gaussian peak fitted to the original image.
//
// MaskGaussSub: assign the coarse-grain region to the average of the original image, minus just the peak part of the Gaussian
// fit (not the constant z0 background part).
//
// v1, 12/11/14 pmv.


function VoronoiCells(im, x0, y0)
	wave im, x0, y0
	
	Make/O/N=(numpnts(x0), 3) vor_tmp
	vor_tmp[][0] = x0[p]
	vor_tmp[][1] = y0[p]
	vor_tmp[][2] = 0
	
	Imagetransform voronoi vor_tmp
	wave cir = $"M_Circles"
	wave tri = $"M_DelaunayTriangles"
	wave vor_edges = $"M_VoronoiEdges"
	
	Killwaves vor_tmp, cir, tri
	
	VoronoiVertices(vor_edges, DimOffset(im, 0), DimOffset(im, 1), DimDelta(im, 0)*DimSize(im, 0) + DimDelta(im, 0), DimDelta(im, 1)*DimSize(im, 1) + DimDelta(im, 1))
	 
end

function MakeAllMasks(im, x0, y0, Voronoi_Vertices, vertices_per_point)
	wave im, x0, y0, Voronoi_Vertices
	variable vertices_per_point
	
	variable tol = 0.05
	
	AllVertNearPoints(x0, y0, Voronoi_Vertices, vertices_per_point, tol)
	wave vert_near_pts = $"vert_near_pts"
	AllVertCircuits(Voronoi_Vertices, vert_near_pts)
	wave vc_x =  $"all_vert_circuits_x"
	wave vc_y = $"all_vert_circuits_y"
	GenerateMasks(x0, y0, vc_x, vc_y, DimSize(im, 0), DimSize(im, 1))
	
end


// Maybe this isn't needed?
function PolyEdgesFromVoronoi(v_edges, start_i)
	wave v_edges
	variable start_i
	
	variable start_x, start_y, end_x, end_y, beg_x, beg_y
	variable done = 0, track_i
	start_x = v_edges[start_i][0]
	start_y = v_edges[start_i][1]
		
	make/o/n=(DimSize(v_edges, 0)) x_tmp, y_tmp, match
	make/o/n=(3, 2) poly_edges
	poly_edges[0][] = v_edges[start_i][q]
	poly_edges[1][] = v_edges[start_i+1][q]
	poly_edges[2][] = NaN
	
	end_x = v_edges[start_i+1][0]
	end_y = v_edges[start_i+1][1]
	variable j = start_i
	do

		x_tmp = ((v_edges[p][0] == end_x) ? 1 : 0)
		y_tmp = ((v_edges[p][1] == end_y) ? 1 : 0)
		match = x_tmp[p] && y_tmp[p]
		FindValue/V=1/S=0 match
		if(V_value == j)
			FindValue/V=1/S=(V_value) match
		endif
		
		if( (v_edges[V_value][0] == start_x) && (v_edges[V_value][1] == start_y) )
			done = 1
		else
			j = V_value
			InsertPoints/M=0 0, 3, poly_edges
			poly_edges[0][0] = end_x
			poly_edges[0][1] = end_y
			
			if(!numtype(v_edges[V_value+1][0]) )
				poly_edges[1][0] = v_edges[V_value+1][0]
				poly_edges[1][1] = v_edges[V_value+1][1]
			else
				poly_edges[1][0] = v_edges[V_value-1][0]
				poly_edges[1][1] = v_edges[V_value-1][1]
			endif
			end_x = poly_edges[1][0]
			end_y = poly_edges[1][1]
			poly_edges[2][0] = NaN
			poly_edges[2][1] = NaN
			
		endif
	
	while(!done)

end

function VoronoiVertices(v_edges, sx, sy, ex, ey)
	wave v_edges
	variable sx, sy, ex, ey
	
	make/o/n=(DimSize(v_edges, 0)) x_tmp, y_tmp
	x_tmp = v_edges[p][0]
	y_tmp = v_edges[p][1]
	
	sort x_tmp y_tmp, x_tmp

	
	variable i = 0, j =0, done1=0, done2 = 0

	do
		if(numtype(x_tmp[i]))
			break
		endif
		i+=1
	while(1)
	DeletePoints/M=0 i, numpnts(x_tmp), x_tmp, y_tmp

	FindLevel/P/Q x_tmp, sx
	if(!V_Flag)
		DeletePoints/M=0 0, V_LevelX+1, x_tmp, y_tmp
	endif
	FindLevel/P/Q x_tmp, ex
	if(!V_flag)
		DeletePoints/M=0 V_LevelX+1, numpnts(x_tmp), x_tmp, y_tmp
	endif

	sort y_tmp, x_tmp, y_tmp
	FindLevel/P/Q y_tmp, sy
	if(!V_Flag)
		DeletePoints/M=0 0, V_LevelX+1, x_tmp, y_tmp
	endif
	FindLevel/P/Q y_tmp, ey
	if(!V_flag)
		DeletePoints/M=0 V_LevelX+1, numpnts(x_tmp), x_tmp, y_tmp
	endif
	

	make/o/n=(0, 2) vor_vert
	i = 0
	do
		InsertPoints/M=0 0, 1, vor_vert
		vor_vert[0][0] = x_tmp[i]
		vor_vert[0][1] = y_tmp[i]
		do
			i+=1
			if(i >= numpnts(x_tmp))
				done2 = 1
			elseif( (x_tmp[i] != vor_vert[0][0]) && (y_tmp[i] != vor_vert[0][1]) )
				done2=1
			endif
		
		while(!done2)
	
		if(i >= numpnts(x_tmp))
			done1 = 1
		else
			done2 = 0
		endif
	
	while(!done1)
	
	Duplicate/O vor_vert Voronoi_Vertices

	Killwaves x_tmp, y_tmp, vor_vert
	
end


function AllVertNearPoints(px, py, vert_list, nvert, tol)
	wave px, py, vert_list
	variable nvert, tol
	
	Make/O/N=(nvert, numpnts(px)) vert_near_pts
	
	variable i
	for(i=0; i<numpnts(px); i+=1)
		VerticesNearPoint(px[i], py[i], vert_list, nvert, tol)
		wave vert_near = $"vert_near"
		vert_near_pts[][i] = vert_near[p]
	endfor
	
	Killwaves vert_near
end

function VerticesNearPoint(px, py, vert, nvert, tol)
	variable px, py
	wave vert
	variable nvert, tol
	
	make/o/n=(DimSize(vert, 0)) vert_d, ind
	vert_d = (px - vert[p][0])^2 + (py - vert[p][1])^2
	ind = p
	sort vert_d, ind, vert_d
	
	Make/O/N=(nvert) vert_near
	vert_near = ind[p]
	
	wavestats/q vert_d
	
	if( (V_max - V_min) / V_min > tol)
		printf "Tolerance exceeded\r"
	endif
	
	Killwaves vert_d, ind	
end
	
function AllVertCircuits(vert_pos, vert_near_pts)
	wave vert_pos, vert_near_pts
	
	Make/O/N=( (DimSize(vert_near_pts, 0)+1), DimSize(vert_near_pts, 1)) all_vert_circuits_x, all_vert_circuits_y
	Make/O/N=(DimSize(vert_near_pts, 0)) one_near
	
	variable i
	for(i=0; i<DimSize(vert_near_pts, 1); i+=1)
		one_near = vert_near_pts[p][i]
		VertCircuit(vert_pos, one_near)
		wave vert_circuit = $"vert_circuit"
		all_vert_circuits_x[][i] = vert_circuit[p][0]
		all_vert_circuits_y[][i] = vert_circuit[p][1]
	endfor
	
	Killwaves one_near, vert_circuit
end
	
function VertCircuit(vert_pos, near_vert)
	wave vert_pos, near_vert

	variable nvert = numpnts(near_vert)

	Make/o/n=(nvert+1, 2) vert_circuit
	
	Make/O/N=(nvert) near_vert_circuit, vert_d, ind
	near_vert_circuit = 0
	ind = 0
	vert_d = 0
	
	near_vert_circuit[0] = near_vert[0]
	
	variable i, j
	for(i=0; i<nvert; i+=1)
	
		vert_d = (vert_pos[near_vert_circuit[i]][0] - vert_pos[near_vert[p]][0])^2 + (vert_pos[near_vert_circuit[i]][1] - vert_pos[near_vert[p]][1])^2
		ind = near_vert[p]
		sort vert_d ind, vert_d
		
		for(j=1; j<nvert; j+=1)
			FindValue/S=0/V=(ind[j]) near_vert_circuit
			if(V_Value == -1)
				near_vert_circuit[i+1] = ind[j]
				break
			endif
		endfor
	endfor

	vert_circuit[0,nvert-1][0] = vert_pos[near_vert_circuit[p]][0]
	vert_circuit[0, nvert-1][1] = vert_pos[near_vert_circuit[p]][1]
	vert_circuit[nvert][0] = vert_pos[near_vert_circuit[0]][0]
	vert_circuit[nvert][1] = vert_pos[near_vert_circuit[0]][1]
	
	Killwaves vert_d, ind, near_vert_circuit
end
		
function GenerateMasks(px, py, all_vert_circuits_x, all_vert_circuits_y, sx, sy)
	wave px, py, all_vert_circuits_x, all_vert_circuits_y
	variable sx, sy
	
	make/O/N=(sx, sy, DimSize(all_vert_circuits_x, 1)) all_masks
	make/O/n=(DimSize(all_vert_circuits_x, 0)) vx, vy
	
	variable i
	for(i=0; i<DimSize(all_vert_circuits_x, 1); i+=1)
		vx = all_vert_circuits_x[p][i]
		vy = all_vert_circuits_y[p][i]
		ImageBoundaryToMask width =(sx), height = (sy), xwave = vx, ywave = vy, seedX = (px[i]), seedY = (py[i])
		wave roi = $"M_ROIMask"
		all_masks[][][i] = roi[p][q]
	endfor
	
	i=0
	duplicate/O px mask_x0
	duplicate/O py mask_y0
	do
		Imagetransform/P=(i) getplane all_masks
		wave M_ImagePlane = $"M_ImagePlane"
		wavestats/q/M=1 M_ImagePlane
		if(V_avg > 0.5) // mask includes more than half the image, so it's a bad mask
			DeletePoints/M=2 i, 1, all_masks
			DeletePoints/M=0 i, 1, mask_x0, mask_y0
		else
			i+=1
		endif
	
	while(i < DimSize(all_masks, 2))
	
	Killwaves vx, vy, roi, M_ImagePlane
end

function MaskCoarseGrain(im, all_masks, im_proc)
	wave im, all_masks
	FUNCREF ProtoMaskProc im_proc
	
	variable i, nmasks = Dimsize(all_masks, 2)
	Duplicate/O im coarse_im
	Duplicate/O im im_tmp
	coarse_im = 0
	variable/G mask_i
	
	for(i=0; i<nmasks; i+=1)
		ImageTransform/P=(i) GetPlane all_masks
		wave one_mask = $"M_ImagePlane"
		im_tmp = im
		mask_i = i
		im_proc(im_tmp, one_mask)
		MatrixOp/O	coarse_im = coarse_im + im_tmp
	endfor
	
	Killwaves one_mask, im_tmp
end

function MaskCoarseGrainStack(im_stack, all_masks, im_proc)
	wave im_stack, all_masks
	FUNCREF ProtoMaskProc im_proc
	
	variable i, nstack = DimSize(im_stack, 2)
	Duplicate/O im_stack coarse_im_stack
	variable/G stack_i
	
	for(i=0; i<nstack; i+=1)
		if(!mod(i, 25))
			printf "Working on image %d\r", i
		endif
		Imagetransform/P=(i) getplane im_stack
		wave M_ImagePlane = $"M_ImagePlane"
		Duplicate/O M_ImagePlane stack_tmp
		Killwaves M_ImagePlane
		stack_i = i
		MaskCoarseGrain(stack_tmp, all_masks, im_proc)
		wave coarse_im = $"coarse_im"
		Imagetransform/P=(i)/D=coarse_im setplane coarse_im_stack
		Killwaves coarse_im
	endfor
	
	Killwaves stack_tmp
	
end

// Mask procedures take two input images, im and mask, that are the same size.
// im contains the full frame image, mask contains the full frame mask, 1 inside
// the ROI and 0 outside.  The procedure replaces im with the result.
function ProtoMaskProc(im, mask)
	wave im, mask
	
	printf "Ended up in ProtoMaskProc.  Should never get here.\r"
end

function MaskAverage(im, mask)
	wave im, mask

	Duplicate/o mask mask_av_tmp
	mask_av_tmp = (mask[p][q] == 0 ? 0 : im[p][q])
		
	WaveStats/M=1/Q mask_av_tmp
	variable im_av = V_avg
	if(V_numNaNs)
		im = (mask[p][q] == 0 ? 0 : NaN)
	else
		WaveStats/M=1/Q mask
		variable mask_av = V_avg
		im = (mask[p][q] == 0 ? 0 : (im_av / mask_av) )
	endif
	
	Killwaves mask_av_tmp
	
end

function MaskZ0(im, mask)
	wave im, mask
	
	NVAR mask_i = $"mask_i"
	NVAR stack_i = $"stack_i"
	wave z0_st = $"z0_st"

	variable z0 = z0_st[mask_i][stack_i]
	
	if(!numtype(z0))
		MatrixOp/O im = z0 * mask
	else
		im = (mask[p][q] == 1 ? NaN : 0)
	endif

end

function MaskSubGauss(im, mask)
	wave im, mask

	NVAR mask_i = $"mask_i"
	NVAR stack_i = $"stack_i"

	wave z0_st = $"z0_st"
	wave A_st = $"A_st"
	wave cor_st = $"cor_st"
	wave x0_st = $"x0_st"
	wave y0_st = $"y0_st"
	wave xW_st = $"xW_st"
	wave yW_st = $"yW_st"

	if(numtype(z0_st[mask_i][stack_i]))
		im = (mask[p][q] == 1 ? NaN : 0)
	else
		Make/O/N=7 W_Coef
		W_Coef[0] = z0_st[mask_i][stack_i]
		W_Coef[1] = A_st[mask_i][stack_i]
		W_Coef[2] = x0_st[mask_i][stack_i]
		W_Coef[3] = xW_st[mask_i][stack_i]
		W_Coef[4] = y0_st[mask_i][stack_i]
		W_Coef[5] = yW_st[mask_i][stack_i]
		W_Coef[6] = cor_st[mask_i][stack_i]

		Duplicate/o mask mask_av_tmp
		mask_av_tmp = (mask[p][q] == 1 ? (im[p][q] - Gauss2D(W_coef, x, y) + z0_st[mask_i][stack_i]) : 0)
		WaveStats/M=1/Q mask_av_tmp
		variable im_av = V_avg

		WaveStats/M=1/Q mask
		variable mask_av = V_avg
		im = (mask[p][q] == 0 ? 0 : (im_av / mask_av) )
	
		Killwaves mask_av_tmp

	endif	

end

function MaskGauss(im, mask)
	wave im, mask

	NVAR mask_i = $"mask_i"
	NVAR stack_i = $"stack_i"

	wave z0_st = $"z0_st"
	wave A_st = $"A_st"
	wave cor_st = $"cor_st"
	wave x0_st = $"x0_st"
	wave y0_st = $"y0_st"
	wave xW_st = $"xW_st"
	wave yW_st = $"yW_st"

	if(numtype(z0_st[mask_i][stack_i]))
		im = (mask[p][q] == 1 ? NaN : 0)
	else
		Make/O/N=7 W_Coef
		W_Coef[0] = z0_st[mask_i][stack_i]
		W_Coef[1] = A_st[mask_i][stack_i]
		W_Coef[2] = x0_st[mask_i][stack_i]
		W_Coef[3] = xW_st[mask_i][stack_i]
		W_Coef[4] = y0_st[mask_i][stack_i]
		W_Coef[5] = yW_st[mask_i][stack_i]
		W_Coef[6] = cor_st[mask_i][stack_i]

		Duplicate/o mask mask_av_tmp
		mask_av_tmp = (mask[p][q] == 1 ? Gauss2D(W_coef, x, y)  : 0)
		WaveStats/M=1/Q mask_av_tmp
		variable im_av = V_avg

		WaveStats/M=1/Q mask
		variable mask_av = V_avg
		im = (mask[p][q] == 0 ? 0 : (im_av / mask_av) )
	
		Killwaves mask_av_tmp

	endif	

end
