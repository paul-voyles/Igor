#pragma rtGlobals=1		// Use modern global access method.
#include "XYZ"

//
// Procedure for common neighbor analysis of a XYZ model structure
//
// Uses 3-index notation from Clarke and Jonsson, PRE 47, 3975 (1993)
// This is similar to Honeycutt-Anderson, but a little simpler.  In this notation
// fcc is entirely 421, hcp is half 421 and half 422, and icosahedron is 555.  This
// code reproduces all of results.
//
// begun 8/15/11 pmv
// v1.0, 8/18/11 pmv


// cutw is a NxN matrix, where N is the number of different components in the model
// It should have all the A-B cutoff distances, sort by atomic number, determined from
// the partial RDFs.
function NeighborList(xyz, cutw)
	wave xyz, cutw
	
	string curfol = GetDataFolder(1)
	NewDataFolder/O/S root:Packages
	NewDataFolder/O/S :cna
	
	variable natoms = DimSize(xyz, 0)
	variable ax = GetAX(xyz)
	variable by = GetBY(xyz)
	variable cz = GetCZ(xyz)
	
	//Set up
	variable nnmax = 20
	Make/O/N=(natoms, nnmax) $(NameofWave(xyz)+"_nn")
	wave nnlist = $(NameofWave(xyz)+"_nn")
	nnlist = NaN
	Duplicate/O cutw cutw2
	cutw2 = cutw^2
	wavestats/q/m=1 cutw2
	variable cutm = V_max
	
	//Set up interatomic distance, tracking waves, and zlist
	Make/o/N=(natoms) vx, vy, vz, nnd, nni, nnz
	variable nnn, i, j, iz, jz
	nnz = xyz[p][0]
	sort nnz, nnz
	make/o/n=1 zlist
	zlist[0] = nnz[0]
	j = 0
	for(i=1; i<natoms; i+=1)
		if(zlist[j] != nnz[i])
			InsertPoints j+1, 1, zlist
			j+=1
			zlist[j] = nnz[i]
		endif
	endfor
	
	for(i=0; i<natoms; i+=1)
		if(numpnts(zlist) == 1)
			iz = 0
		else
			FindLevel/P/Q zlist, xyz[i][0]
			iz = V_LevelX
		endif
		
		// calculate pair vectors
		vx = xyz[p][1] - xyz[i][1]
		vy = xyz[p][2] - xyz[i][2]
		vz = xyz[p][3] - xyz[i][3]
		
		//apply periodic boundary conditions
		vx = (vx[p] > 0.5*ax ? vx[p]-ax : vx[p])
		vx = (vx[p] < -0.5*ax ? vx[p]+ax : vx[p])
		vy = (vy[p] > 0.5*by ? vy[p]-by : vy[p])
		vy = (vy[p] < -0.5*by ? vy[p]+by : vy[p])
		vz = (vz[p] > 0.5*cz ? vz[p]-cz : vz[p])
		vz = (vz[p] < -0.5*cz ? vz[p]+cz : vz[p])
		
		// calculate square of neighbor distances
		nnd = ( vx[p]^2 +  vy[p]^2 +  vz[p]^2 )
		nni = p
		nnz = xyz[p][0]
		Sort nnd, nnd, nni, nnz
		FindLevel/P/Q nnd, cutm
		if(V_flag)
			if(nnd[0] > cutm)
				nnn = 0
			else
				nnn = numpnts(nnd)-1
			endif
		else
			nnn = floor(V_LevelX)
		endif
		// printf "Number of nearest neighbors = %d\r", nnn
		if(nnn > nnmax)
			printf "Too many candidate neighbors.  %d > %d max.  Exiting.\r" nnn, nnmax
			return 0
		endif
		
		make/o/n=(nnn+1) tnni
		for(j=1; j<nnn+1; j+=1)
			FindLevel/P/Q zlist, nnz[j]
			jz = V_LevelX
		
			if(nnd[j] <= cutw2[iz][jz])
				tnni[j] = nni[j]
			else
				tnni[j] = nan
			endif
		endfor
		sort tnni, tnni
		nnlist[i][0,nnn-1] = tnni[q+1]
		
	endfor
	
	wavestats/q/m=1 nnlist
	printf "Average number of nearest neighbors = %g\r" (V_npnts/natoms)
	
	Killwaves cutw2, nnd, nni, nnz, tnni, vx, vy, vz, zlist
	SetDataFolder $curfol
end

function AllCNA(xyz)
	wave xyz
	
	string curfol = GetDataFolder(1)
	NewDataFolder/O/S root:Packages
	NewDataFolder/O/S :cna

	wave nnlist = $(NameofWave(xyz)+"_nn")
	if(!WaveExists(nnlist))
		printf "Neighbor list must be run first.  Exiting.\r"
		return 0
	endif
		
	make/o/n=0 cna_w
	make/o/n=0 cna_r1, cna_r2
	
	variable i, j, m
	for(i=0; i<DimSize(nnlist, 0); i+=1)
		j = 0
		do
			if(numtype(nnlist[i][j]) == 2)
				break
			endif
			InsertPoints 0, 1, cna_r1, cna_r2
			cna_r1[0] = i
			cna_r2[0] = nnlist[i][j]
			InsertPoints 0, 1, cna_w
			cna_w[0] = OneCNA(i, j, nnlist)
			j+=1
		while(1)
	endfor

	Sort cna_w, cna_w, cna_r1, cna_r2

	Duplicate/O cna_w $(curfol+nameofwave(xyz)+"_cna")
	Duplicate/O cna_r1 $(curfol+nameofwave(xyz)+"_cnr1")
	Duplicate/O cna_r2 $(curfol+nameofwave(xyz)+"_cnr2")
	
	Killwaves cna_w, cna_r1, cna_r2
	
	SetDatafolder $curfol

end

// r1 is the first atom list, aa is the position in it's neighbor list (starting at zero) of the second atom
// in the root pair
// first digit is the number of common neighbors.
// second digit is the number of bonds within the set of common neighbors
// third digit is the longest circuit of bonds with the set of common neighbors
function OneCNA(r1, aa, nnlist)
	variable r1, aa
	wave nnlist
	
	variable nnmax = 20
	make/o/n=(nnmax) r1list, r2list, rtmp
	
	// r1 and r2 are the indices of the root pair
	variable cna, r2, r1n, r2n, i
	r2 = nnlist[r1][aa]
	
	// Find their neighbors
	r1list = nnlist[r1][p]
	r2list = nnlist[r2][p]
	rtmp = numtype(r1list)
	FindLevel/q/p rtmp, 2
	r1n = V_LevelX
	rtmp = numtype(r2list)
	FindLevel/q/p rtmp, 2
	r2n = V_LevelX

	// printf "Base pair is %d and %d.  ", r1, r2
	//printf "r1 has %d neighbors.\r", r1n
	//printf "r2 has %d neighbors.\r", r2n
	
	// Find their _common_ neighbors
	make/o/n=0 common
	variable ii, jj, kk, cn = 0
	for(ii=0; ii<r1n; ii+=1)
		for(jj=0; jj<r2n; jj+=1)
			if(r1list[ii] == r2list[jj])
				InsertPoints 0, 1, common
				common[0] = r2list[jj]
				jj=r2n
			endif
		endfor
	endfor
	cn = numpnts(common)
	cna = cn*100
	if(cn == 0)
		return cna
	endif
	
	// construct the neighbor list just for the common neighbor atoms, excluding the
	// original pair.
	make/o/n=(cn, nnmax) common_nn
	for(ii=0; ii<cn; ii+=1)
		common_nn[ii][] = nnlist[common[ii]][q]
	endfor
	variable onek
	for(ii=0; ii<cn; ii+=1)
		for(jj=0; jj<nnmax; jj+=1)
			onek = nan
			for(kk=0; kk<cn; kk+=1)
				if(numtype(common_nn[ii][jj]) == 2)
					kk=cn
					jj=nnmax
				elseif(common_nn[ii][jj] == common[kk])
					onek = kk
					kk=cn
				endif
			endfor
			common_nn[ii][jj] = onek
		endfor
	endfor
	
	for(ii=0; ii<cn; ii+=1)
		rtmp = common_nn[ii][p]
		sort rtmp, rtmp
		common_nn[ii][] = rtmp[q]
	endfor
	
	// Second CNA digit is the number of bonds in the common neighbor network.  Each bond is 
	// counted twice, so that's just half the number of non-NaN entries in common_nn.
	wavestats/q/m=1 common_nn
	cna += 5*V_npnts
	
	// Last digit in the longest bond path through the common neighbor bond network.  Need a recursive
	// search through all possible paths to find the longest one.
	cna += LongestPath(common_nn)
	
	// printf "CNA = %d\r", cna
	Killwaves common_nn, rtmp, r1list, r2list
	return cna
end


// find the longest circuit through the bond graph in nnlist.  (Turns out CNA only cares about chains,
// not circuits, but I was calculating both anyway.
function LongestPath(nnlist)
	wave nnlist
	
	variable lp = 0, i, j
	string p = ""
	string q
	variable nit
	
	for(i=0; i<DimSize(nnlist, 0); i+=1)
		make/O/T/N=1 one_atom_paths
		sprintf p, "%d;", i
		one_atom_paths[0] = p
		OneAtomPaths(one_atom_paths, nnlist, i, 0)
		for(j=0; j<numpnts(one_atom_paths); j+=1)
			q = one_atom_paths[j]
			//printf "%s\r", q
			nit = ItemsInList(q, ";")
			lp = max(lp, nit-2)
		endfor
		//printf "lp = %d\r", lp
	endfor
	
	Killwaves one_atom_paths
	return lp
	
end


function OneAtomPaths(paths,  nnlist, ni,listn)
	wave/t paths
	wave nnlist
	variable ni, listn
	
	string p = paths[listn]
	
	string next_i
	string last_i
	variable pos
	
	if(ItemsInList(p, ";") > 1)
		// get the last list item
		last_i = StringFromList(ItemsinList(p, ";")-1, p, ";")	

		// does the last item occur elsewhere in the list?
		pos = WhichListItem(last_i, p, ";")
		if(pos == 0) // last item is the same as first item, a good circuit, and we're done
			if(ItemsInList(p, ";") == 3) // list is only one bond long, like 0, 3, 0, and it's bad
				p=RemoveListItem(ItemsinList(p, ";")-1, p, ";")
			endif
			p = AddListItem("done", p, ";", ItemsInList(p, ";"))
			paths[listn] = p
			return 0
		elseif(pos > 0 && pos != (ItemsInList(p, ";")-1)) // last item is present, but is not the first item, and the circuit is bad
			p = RemoveListItem(ItemsInList(p, ";")-1, p, ";")  // remove the wrong bond from the list
			p = AddListItem("done", p, ";", ItemsInList(p, ";'"))
			paths[listn] = p
			return 0
		endif
	endif
	
	// this last item on the list is a new item.  Keep going by adding all of its neighbors
	// to new lists
	variable i=0	
	string q
	do
		if(numtype(nnlist[ni][i]) == 2)
			break
		endif
		sprintf next_i, "%d", nnlist[ni][i]
		q = AddListItem(next_i, p, ";", ItemsInList(p, ";"))
		if(i==0)
			paths[listn] = q
			OneAtomPaths(paths, nnlist, nnlist[ni][i], listn)
		else
			InsertPoints numpnts(paths), 1, paths
			listn+=1
			paths[listn] = q
			OneAtomPaths(paths, nnlist, nnlist[ni][i], listn)
		endif
		i+=1	
	while(1)

end
			 
			 
function distance(n1, n2, xyz)
	variable n1, n2
	wave xyz
	
	variable ax = GetAX(xyz)
	variable by = GetBY(xyz)
	variable cz = GetCZ(xyz)
	
	variable vx, vy, vz
	vx = xyz[n1][1] - xyz[n2][1]
	vy = xyz[n1][2] - xyz[n2][2]
	vz = xyz[n1][3] - xyz[n2][3]

	if(vx > 0.5*ax)
		vx -= ax
	endif
	if(vx < -0.5*ax)
		vx += ax
	endif

	if(vy > 0.5*by)
		vy -= by
	endif
	if(vy < -0.5*by)
		vy += by
	endif
	
	if(vz > 0.5*cz)
		vz -= cz
	endif
	if(vz < -0.5*cz)
		vz += cz
	endif
	
	variable d
	d = sqrt( vx^2 + vy^2 + vz^2)
	//d = sqrt( (xyz[n1][1] - xyz[n2][1])^2 + (xyz[n1][2] - xyz[n2][2])^2 + (xyz[n1][3] - xyz[n2][3])^2 )
	
	//printf "distance between %d and %d is %g\r", n1, n2, d
	return d
	
end

function CheckNN(nnlist, xyz)
	wave nnlist
	wave xyz
	
	variable d, dmin, dmax
	dmin = 2.88
	dmax = 2.951
	
	variable i, j, count = 0
	for(i=0; i<DimSize(nnlist,0); i+=1)
		if(!mod(count, 100))
			printf "step %d\r", count
		endif
		j=0
		do
			if(numtype(nnlist[i][j]) == 2)
				break
			endif
			d = distance(i, nnlist[i][j], xyz)
			if(d < dmin || d > dmax)
				printf "d of %g is outside bounds for atoms %d, %d.\r", d, i, nnlist[i][j]
			endif
			j+=1
			count+=1
		while(1)
	endfor
	
end

// make an xyz structure wave for the atoms ii and jj in list xyz.  Fails if those atoms are not
// neighbors.
function PullOneCNA(ii, jj, xyz)
	variable ii, jj
	wave xyz
	
	string curfol = GetDataFolder(1)
	SetDataFolder root:Packages:cna
	
	wave nnlist = $(NameOfWave(xyz) + "_nn")
	variable ni = -1, i
	for(i=0; i<DimSize(nnlist, 1); i+=1)
		if(nnlist[ii][i] == jj)
			ni = i
		endif
	endfor
	if(ni == -1)
		printf "Cannot find atom %d in the neighbor list of atom %d.  Exiting.\r"
		return 0
	endif
	
	
	variable cna
	cna	= OneCNA(ii, ni, nnlist)
	printf "CNA for atom %d and its neighbor %d is %d\r", ii, nnlist[ii][ni], cna
	wave common = $"common"
	
	make/o/n=(  (numpnts(common)+2), 6) cna_xyz
	
	for(i=0; i<numpnts(common); i+=1)
		cna_xyz[i][] = xyz[common[i]][q]
	endfor
	
	cna_xyz[i][] = xyz[ii][q]
	cna_xyz[i+1][] = xyz[nnlist[ii][ni]][q]
	cna_xyz[i][0] += 2
	cna_xyz[i+1][0] += 2
	
	string n = note(xyz)
	Note cna_xyz, n
	
	Duplicate/O cna_xyz   $(curfol+"cna_xyz")
	
	killwaves cna_xyz
	
	SetDataFolder $curfol
	
end
	
function CNAHist(cna)
	wave cna
	
	make/o/n=1 cna_type, cna_count
	cna_type[0] = cna[0]
	cna_count[0] = 1
	
	variable i, j = 0
	for(i=1; i<numpnts(cna); i+=1)
		if(cna[i] == cna_type[j])
			cna_count[j] += 1
		else
			InsertPoints numpnts(cna_type), 1, cna_type, cna_count
			j+=1
			cna_type[j] = cna[i]
			cna_count[j] += 1
		endif
	
	endfor

	sort/r cna_count, cna_count, cna_type
	wavestats/q cna_count
	Duplicate/O cna_count cna_norm
	cna_norm /= V_sum
	Duplicate/O cna_norm cna_cum
	cna_cum = sum(cna_norm, 0, p)
	
end
	