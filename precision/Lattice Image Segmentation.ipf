#pragma rtGlobals=3		// Use modern global access method and strict wave access.

function StackGaussFit(image_stack, x_loc, y_loc, size, wiggle, fixposYN)
	wave image_stack, x_loc, y_loc
	variable size, wiggle, fixposYN
	
	if(fixposYN)
		printf "Fitting with fixed atom center positions not yet supported.  Exiting.\r"
		return 0
	endif
	
	Make/O/N=(DimSize(image_stack, 2)) xav, yav, xstd, ystd
	
	variable i
	for(i=0; i<DimSize(image_stack, 2); i+=1)
	//for(i=0; i<4; i+=1)
		Imagetransform/p=(i) getplane image_stack
		wave im = $"M_ImagePlane"
		GaussianFit(im, x_loc, y_loc, size, wiggle)
			
		//save fit parameters for each atom column in each frame of stack
		wave z0 = $"z0"
		wave A = $"A"
		wave x0 = $"x0"
		wave xW = $"xW"
		wave y0 = $"y0"
		wave yW = $"yW"
		wave cor = $"cor"
		wave sigma_z0 = $"sigma_z0"
		wave sigma_A = $"sigma_A"
		wave sigma_x0 = $"sigma_x0"
		wave sigma_xW = $"sigma_xW"
		wave sigma_y0 = $"sigma_y0"
		wave sigma_yW = $"sigma_yW"
		wave sigma_cor = $"sigma_cor"
		variable num_peaks = DimSize(z0,0)

		if (i == 0)
			Make/O/N=(num_peaks, DimSize(image_stack, 2)) z0_stack, A_stack, x0_stack, xW_stack, y0_stack, yW_stack, cor_stack 
			Make/O/N=(num_peaks, DimSize(image_stack, 2)) sig_z0_stack, sig_A_stack, sig_x0_stack, sig_xW_stack, sig_y0_stack, sig_yW_stack, sig_cor_stack 
		endif

		z0_stack[][i] = z0[p] 
		A_stack[][i] = A[p] 
		x0_stack[][i] = x0[p]  
		xW_stack[][i] = xW[p]  
		y0_stack[][i] = y0[p]  
		yW_stack[][i] = yW[p]  
		cor_stack[][i] = cor[p] 
		sig_z0_stack[][i] = sigma_z0[p] 
		sig_A_stack[][i] = sigma_A[p] 
		sig_x0_stack[][i] = sigma_x0[p]  
		sig_xW_stack[][i] = sigma_xW[p]  
		sig_y0_stack[][i] = sigma_y0[p]  
		sig_yW_stack[][i] = sigma_yW[p]  
		sig_cor_stack[][i] = sigma_cor[p]  
	
		killwaves z0, A, x0, xW, y0, yW, cor, sigma_z0, sigma_A, sigma_x0, sigma_xW, sigma_y0, sigma_yW, sigma_cor	
	endfor	
	
	
end

function StackFitUnwrap(dat, im_x, im_y)
	wave dat
	variable im_x, im_y
	
	variable nimages = DimSize(dat, 1)
	make/o/n=(im_x, im_y, nimages) stack_unwrap
	Make/O/n=(im_x, im_y) im_tmp
	make/o/n=(DimSize(dat, 0)) unwrap_tmp
	
	variable i
	for(i=0; i<nimages; i+=1)
		unwrap_tmp = dat[p][i]
		ImageTransform/D=unwrap_tmp fillImage im_tmp
		ImageTransform/P=(i)/D=stack_unwrap setPlane im_tmp
	endfor
	
	Killwaves im_tmp, unwrap_tmp
	
end

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

function MakeAllMasks(im, x0, y0, M_VoronoiEdges, Voronoi_Vertices, vertices_per_point)
	wave im, x0, y0, M_VoronoiEdges, Voronoi_Vertices
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

// this doesn't work
function VerticesBondGraph(vor_vert, cut)
	wave vor_vert
	variable cut
	
	variable max_con = 6  // maximum connections allowed for each vertex
	
	variable nvert = DimSize(vor_vert, 0)
	
	Make/O/N=(nvert, nvert) pair_dist
	
	pair_dist = (vor_vert[p][0] - vor_vert[q][0])^2 + (vor_vert[p][1] - vor_vert[q][1])^2
	pair_dist = sqrt(pair_dist)
	
	Make/O/N=100 pair_dist_Hist
	Histogram/B={0,(2*cut/100),100} pair_dist,pair_dist_Hist
	
	pair_dist = (pair_dist <= cut ? 1 : 0)
	
	Make/O/N=(nvert, 6) vert_bonds
	vert_bonds = NaN
	Make/O/N=(nvert) list_tmp
	
	variable i, j, k
	for(i =0; i<nvert; i+=1)
		list_tmp = pair_dist[p][i]
		j = 0
		k = 0
		do
			FindValue/S=(j)/V=1 list_tmp
			if(V_Value != -1)
				j = V_value + 1
				if(V_value != i)
					vert_bonds[i][k] = V_value
					k += 1
				endif
			endif
		while(V_Value != -1 && j < nvert) 
	endfor
	
	// Killwaves list_tmp, pair_dist
	
end

// not used
function AllVerticesBondCircuits(vert_bonds, max_d)
	wave vert_bonds
	variable max_d
	
	Make/O/N=(DimSize(vert_bonds, 0), max_d, max_d) bond_circuits
	bond_circuits = NaN

	variable i
	for(i=0; i<DimSize(vert_bonds, 0); i+=1)
		OneVertexBondCircuits(vert_bonds, max_d, i)
	endfor
	
end

// doesn't work
function OneVertexBondCircuitsRecur(vert_bonds, max_d, start_vert, circuit_i, member_i, vert_i)
	wave vert_bonds
	variable max_d, start_vert, circuit_i, member_i, vert_i
	
	wave bond_circuits = $"bond_circuits"
	
	if(member_i >= max_d)
		return -1
	endif
	if(circuit_i >= max_d)
		return -1
	endif
	
	variable j=0, new_vert, cc
	do
		printf "Starting on %g: circuit %g; current vertex %g; neighbor %g = %g\r", start_vert, circuit_i, vert_i, j, vert_bonds[vert_i][j]
		if(numtype(vert_bonds[vert_i][j]))
			return circuit_i + 1
		endif
		bond_circuits[start_vert][member_i][circuit_i] = vert_bonds[vert_i][j]
		if(vert_bonds[vert_i][j] != start_vert)
			cc = OneVertexBondCircuitsRecur(vert_bonds, max_d, start_vert, circuit_i, member_i+1, vert_bonds[vert_i][j])
		endif
		j+=1
	while(1)
	return cc
	
end

function OneVertexBondCircuits(vert_bonds, max_d, start_vert)
	wave vert_bonds
	variable max_d, start_vert
	
	
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
	
	for(i=0; i<nmasks; i+=1)
		ImageTransform/P=(i) GetPlane all_masks
		wave one_mask = $"M_ImagePlane"
		im_tmp = im
		im_proc(im_tmp, one_mask)
		coarse_im += im_tmp
	endfor
	
	Killwaves one_mask, im_tmp
end

function MaskCoarseGrainStack(im_stack, all_masks, im_proc)
	wave im_stack, all_masks
	FUNCREF ProtoMaskProc im_proc
	
	variable i, nstack = DimSize(im_stack, 2)
	Duplicate/O im_stack coarse_im_stack
	
	for(i=0; i<nstack; i+=1)
		if(!mod(i, 25))
			printf "Working on image %d\r", i
		endif
		Imagetransform/P=(i) getplane im_stack
		wave M_ImagePlane = $"M_ImagePlane"
		Duplicate/O M_ImagePlane stack_tmp
		Killwaves M_ImagePlane
		MaskCoarseGrain(stack_tmp, all_masks, im_proc)
		wave coarse_im = $"coarse_im"
		Imagetransform/P=(i)/D=coarse_im setplane coarse_im_stack
		Killwaves coarse_im
	endfor
	
	Killwaves stack_tmp
	
end

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