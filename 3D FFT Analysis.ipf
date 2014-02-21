#pragma rtGlobals=1		// Use modern global access method.
// Functions for analysis of 3D FFTs of model structures
//
// begun 03-05-09 pmv

#include <GizmoSlicer>

function IsoAverage3D(dat, strip_width)
	wave dat
	variable strip_width
	
	variable x0, y0, z0, dx, dy, dz
	x0 = DimOffSet(dat, 0)
	y0 = DimOffSet(dat, 1)
	z0 = DimOffset(dat, 2)
	dx = DimDelta(dat, 0)
	dy = DimDelta(dat, 1)
	dz = DimDelta(dat, 2)
	
	// assume dat is centered at (0,0,0) and extends from (-x0 to +x0), etc.
	Make/O/N=(min( min(DimSize(dat, 0), DimSize(dat, 1)), Dimsize(dat, 2)) / (2*strip_width) +1) iso_av, npix
	SetScale/I x, 0, min( min(abs(x0), abs(y0)), abs(z0)), iso_av, npix
	iso_av = 0
	npix = 0
	
	variable i, j, k, r
	for(i=0; i<DimSize(dat, 0); i+=1)
		for(j=0; j<DimSize(dat, 1); j+=1)
			for(k=0; k<DimSize(dat, 2); k+=1)
				r = sqrt( (x0+dx*i)^2 + (y0+dy*j)^2 + (z0+dz*k)^2 )
				if(r < rightx(iso_av))
					npix[x2pnt(npix, r)] += 1
					iso_av[x2pnt(iso_av, r)] += dat[i][j][k]
				endif
			endfor
		endfor
	endfor	
		
	iso_av /= npix	
	Killwaves npix
end

function RadiusChop(dat, rmin, rmax)
	wave dat
	variable rmin, rmax
	
	Duplicate/O dat radius_chop
	radius_chop = 0

	variable x0, y0, z0, dx, dy, dz
	x0 = DimOffSet(dat, 0)
	y0 = DimOffSet(dat, 1)
	z0 = DimOffset(dat, 2)
	dx = DimDelta(dat, 0)
	dy = DimDelta(dat, 1)
	dz = DimDelta(dat, 2)
	
	variable i, j, k, r
	for(i=0; i<DimSize(dat, 0); i+=1)
		for(j=0; j<DimSize(dat, 1); j+=1)
			for(k=0; k<DimSize(dat, 2); k+=1)
				r = sqrt( (x0+dx*i)^2 + (y0+dy*j)^2 + (z0+dz*k)^2 )
				if( (r <= rmax) && (r >= rmin) )
					radius_chop[i][j][k] = dat[i][j][k]
				endif
			endfor
		endfor
	endfor	
	
end

Function makeSphere(pointsx,pointsy)
	Variable pointsx,pointsy
	
	Variable i,j,rad
	Make/O/n=(pointsx,pointsy,3) sphereData
	Variable anglePhi,angleTheta
	Variable dPhi,dTheta
	
	
	dPhi=2*pi/(pointsx-1)
	dTheta=pi/(pointsy-1)
	Variable xx,yy,zz
	Variable sig
	
	for(j=0;j<pointsy;j+=1)
		angleTheta=j*dTheta
		zz=sin(angleTheta)
		if(angleTheta>pi/2)
			sig=-1
		else
			sig=1
		endif
		for(i=0;i<pointsx;i+=1)
			anglePhi=i*dPhi
			xx=zz*cos(anglePhi)
			yy=zz*sin(anglePhi)
			sphereData[i][j][0]=xx
			sphereData[i][j][1]=yy
			sphereData[i][j][2]=sig*sqrt(1-xx*xx-yy*yy)
		endfor
	endfor
End

Function makeSphereColorWave(pointsx,pointsy,dat,radius,thick)
	Variable pointsx,pointsy
	wave dat
	variable radius, thick
	
	Variable i,j,rad
	Make/O/N=(pointsx, pointsy) sphereGray
	SetScale/I x 0, 2*Pi, "phi", sphereGray
	SetScale/I y 0, Pi, "theta", sphereGray
	Make/O/n=(pointsx,pointsy,4) sphereColor
	Variable anglePhi,angleTheta
	Variable dPhi,dTheta
	string n
	
	dPhi=2*pi/(pointsx-1)
	dTheta=pi/(pointsy-1)
	Variable xx,yy,zz
	Variable sig
	
	for(j=0;j<pointsy;j+=1)
		angleTheta=j*dTheta
		zz=sin(angleTheta)
		if(angleTheta>pi/2)
			sig=-1
		else
			sig=1
		endif
		for(i=0;i<pointsx;i+=1)
			anglePhi=i*dPhi
			
			xx = radius*sin(angleTheta)*cos(anglePhi)
			yy = radius*sin(angleTheta)*sin(anglePhi)
			zz = radius*cos(angleTheta)
			sphereGray[i][j] = VolumeSum(dat, thick, xx, yy, zz)
		endfor
	endfor

	n = ""
	n = ReplaceNumberByKey("radius", n, radius, "=")
	n = ReplaceNumberByKey("thick",  n, thick, "=")
	Note/k sphereGray, n

	wavestats/Q sphereGray
	sphereColor[][][0,2] = (sphereGray[p][q] - V_min) / (V_max - V_min)
	sphereColor[][][3] = 1

End

function VolumeSum(dat, radius, xc, yc, zc)
	wave dat
	variable radius, xc, yc, zc
	
	variable rp, xp, yp, zp, dist, s, ii, jj, kk
	
	xp = round( (xc - DimOffset(dat, 0)) / DimDelta(dat, 0))
	yp =  round( (yc - DimOffset(dat, 1)) / DimDelta(dat, 1))
	zp =  round( (zc - DimOffset(dat, 2)) / DimDelta(dat, 2))
	rp = ceil(radius / DimDelta(dat, 0))
	
	if(rp<=1)
		printf "Radius too small.  oops.\r"
	endif
	
	//printf "Averaging %d pixels around position (%d, %d, %d).\r", rp, xp, yp, zp
	
	s=0
	for(ii=xp-rp; ii<=xp+rp; ii+=1)
		for(jj=yp-rp; jj<=yp+rp; jj+=1)
			for(kk=zp-rp; kk<=zp+rp; kk+=1)
				dist = (ii-xp)^2 + (jj-yp)^2 + (kk-zp)^2
				if(dist < rp^2)
					s += dat[ii][jj][kk]
				endif
			endfor
		endfor
	endfor
	
	return s
	
end

function Sphere2XYZ(phi, theta, sphereGray)
	variable phi, theta
	wave sphereGray
	
	variable radius
	radius = NumberByKey("radius", note(sphereGray),  "=")
	if(numtype(radius) == 2)
		printf "Cannot find radius of sphere-surface wave.  Exiting.\r"
		return 0
	endif
	
	printf "(x, y, z) = (%g, %g, %g)\r", radius*sin(theta)*cos(phi), radius*sin(theta)*sin(phi), radius*cos(theta)
	
end
