#pragma rtGlobals=1		// Use modern global access method.
//#include "XYZ"



function STEMpot(xyz, ax, by, nx, ny)
	wave xyz
	variable ax, by, nx, ny
	
	variable zpow = 1.7  // 2 for Rutherford, 1.7 for screening
	
	variable dx = ax / nx
	variable dy = by / ny
	
	Make/D/O/N=(nx, ny) pot
	SetScale/I x 0, ax, "", pot
	SetScale/i y 0, by, "", pot
	
	variable natoms = Dimsize(xyz, 0)
	Make/O/N=(natoms) iax, ibx, fax, iay, iby, fay
	make/O/D/N=(natoms) V1, V2, V3, V4
	
	iax = floor(xyz[p][1] / dx)
	ibx = mod(iax, nx) + 1
	fax = 1 - mod( (xyz[p][1] / dx), 1 )
	
	iay = floor( xyz[p][2] /dy )
	iby = mod(iay, ny) + 1
	fay = 1 - mod( (xyz[p][2]/dy) , 1)
	
	V1 = fax[p]*fay[p]*xyz[p][0]^zpow	
	V2 = (1 - fax[p])*fay[p]*xyz[p][0]^zpow
	V3 = fax[p]*(1-fay[p])*xyz[p][0]^zpow
	V4 = (1-fax[p])*(1-fay[p])*xyz[p][0]^zpow
	
	pot = 0
	variable i
	for(i=0; i<natoms; i+=1)
		pot[iax[i]][iay[i]] += V1[i]
		pot[ibx[i]][iay[i]] += V2[i]
		pot[iax[i]][iby[i]] += V3[i]
		pot[ibx[i]][iby[i]] += V4[i]
	endfor
	
	KillWaves fay, iby, iay, fax, ibx, iax, v1, v2, v3, v4 
end
	
	
function STEMpot3D(xyz, nx, ny, nz)
	wave xyz
	variable nx, ny, nz
	
	variable ax = GetAX(xyz)
	variable by = GetBY(xyz)
	variable cz = GetCZ(xyz)
	
	Make/o/N=(nx, ny, nz) pot3D
	SetScale/I x 0, ax, "", pot3D
	SetScale/I y 0, by, "", pot3D
	SetScale/I z 0, cz, "", pot3D
	
	variable zpow = 1.7  // 2 for Rutherford, 1.7 for screening
	
	variable dx = ax / nx
	variable dy = by / ny
	variable dz = cz / nz
		
	variable natoms = Dimsize(xyz, 0)
	Make/O/N=(natoms) iax, ibx, fax, iay, iby, fay, iaz, ibz, faz
	make/O/D/N=(natoms) V1, V2, V3, V4, V5, V6, V7, V8
	
	iax = floor( xyz[p][1] / dx)
	ibx = mod(iax, nx) + 1
	fax = 1 - mod( (xyz[p][1] / dx), 1 )
	
	iay = floor( xyz[p][2] /dy )
	iby = mod(iay, ny) + 1
	fay = 1 - mod( (xyz[p][2]/dy) , 1)
	
	iaz = floor( xyz[p][3] / dz )
	ibz = mod(iaz, nz ) + 1
	faz = 1 - mod( (xyz[p][3]/dz), 1)
	
	V1 =    fax[p]   *    fay[p]*      faz[p]   * xyz[p][0]^zpow	
	V2 =    fax[p]   *    fay[p]   * (1-faz[p]) * xyz[p][0]^zpow
	V3 = (1-fax[p]) *    fay[p]   *    faz[p]   * xyz[p][0]^zpow
	V4 = (1-fax[p]) *    fay[p]   * (1-faz[p]) * xyz[p][0]^zpow
	V5 =    fax[p]   * (1-fay[p]) *    faz[p]   * xyz[p][0]^zpow
	V6 =    fax[p]   * (1-fay[p]) * (1-faz[p]) * xyz[p][0]^zpow
	V7 = (1-fax[p]) * (1-fay[p]) *    faz[p]   * xyz[p][0]^zpow
	V8 = (1-fax[p]) * (1-fay[p]) * (1-faz[p]) * xyz[p][0]^zpow
	
	pot3D = 0
	variable i
	for(i=0; i<natoms; i+=1)
		pot3D[iax[i]][iay[i]][iaz[i]] += V1[i]
		pot3D[iax[i]][iay[i]][ibz[i]] += V2[i]
		pot3D[ibx[i]][iay[i]][iaz[i]] += V3[i]
		pot3D[ibx[i]][iay[i]][ibz[i]] += V4[i]
		pot3D[iax[i]][iby[i]][iaz[i]] += V5[i]
		pot3D[iax[i]][iby[i]][ibz[i]] += V6[i]
		pot3D[ibx[i]][iby[i]][iaz[i]] += V7[i]
		pot3D[ibx[i]][iby[i]][ibz[i]] += V8[i]
	endfor
	
	KillWaves fay, iby, iay, fax, ibx, iax, faz, iaz, ibz, v1, v2, v3, v4, v5, v6, v7, v8 
	
 	
end