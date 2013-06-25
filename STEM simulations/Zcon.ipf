#pragma rtGlobals=1		// Use modern global access method.
#include "XYZ"
#include "STEMpot"
#include "STEMpsf"

// 07-03-07 fixed bug in the scaling of psf2d that doubled the resolution of all simulated images pmv
// 01-20-10 updated to use new STEMPSF routines incorporating higher-order aberrations and chromatic
// aberration
// 02-12-10 updated to use new STEMPSF routines incorporating source size
//
// Aberration parameters are:
//	C1, defocus, in nm
//	A1, 2-fold astigmatism, in nm
//	A2, 3-fold astigmatism, in nm
//	B2, axial coma, in nm
//	C3, primary spherical aberration, in um
//	A3, 4-fold astigmasitm, in um
//	S3, star aberration, in um
//	A4, 5-fold astigmatism, in um
//	D4, 3-lobe aberration, in um
//	B4, axial coma, in um
//	C5, 5th order spherical aberration, in mm
//	A5, 5th order spherical aberration, in mm
//
// Aberrations must be specified in a 12x2 wave.  The first column is the aberration coefficient, the
// second is the rotation angle.  Rotation angles for axially symmetric aberrations (C1, C3, C5) are
// ignored.  Aberrations can be set to zero.
//
// Other required parameters are:
// Cc = chromatic aberration coefficient in mm
// dE = beam energy spread (FWHM) in eV
// keV = beam energy in kV
// ds = source function FWHM.  Source is assumed to be Gaussian
// ap = condenser aperture in mrad


function zcon(xyz, aber, keV, Cc, dE, ds, ap)
	wave xyz, aber
	variable keV, Cc, dE, ds, ap
	
	variable ax = GetAX(xyz)
	variable by = GetBY(xyz)
	variable rmax = max(ax, by)
	variable nk
	
	// STEMpsf must go from -rmax/2 to rmax/2 to get scaling correct
	// need to set nk to make that happen
	nk = floor(40*(0.001*ap / wavlen(keV))*rmax)
	if( mod(nk, 2) )
		nk += 1
	endif
	
	printf "Calculating PSF with %d pixels.\r", nk

	if(Cc == 0.0 || dE == 0.0)
		STEMPSF2DCoh(aber, keV, ap, nk)
		Redimension/D/C probe2DCoh
		wave/c psf2D = $"probe2DCoh"
	else
		STEMPSF2DIncoh(aber, keV, Cc, dE, ds, ap, nk)
		Redimension/D/C probe2DIncoh
		wave/c psf2D = $"probe2DIncoh"
	endif

	if(!waveexists(psf2d))
		printf "stempsf2d has failed for some reason.  cannot continue.\r"
		return 0
	endif
	
	printf "Calculating potential.\r"
	stempot(xyz, rmax, rmax, nk, nk)
	Redimension/C pot
	wave/c pot = $"pot"
	if(!waveexists(pot))
		printf "stempot has failed for some reason.  cannot continue.\r"
		return 0
	endif
	
	printf "Calculating image.\r"
	fft psf2d
	fft pot
	FastOp/c pot = pot*psf2d
	ifft/c pot
	
	duplicate/o pot zcon_im
	Redimension/R zcon_im
	zcon_im = real(pot)
	imagetransform swap zcon_im
	SetScale/P x DimOffSet(pot, 0), DimDelta(pot, 0), "", zcon_im
	SetScale/P y DimOffSet(pot, 1), DimDelta(pot, 1), "", zcon_im
	
	Killwaves psf2d, pot
end


function zcon3D(xyz, C5, Cs, df, keV, ap, nx, nz)
	wave xyz
	variable C5, Cs, df, keV, ap, nx, nz
	
	string curfol = GetDataFolder(1)
	NewDataFolder/O/S Packages
	NewDataFolder/O/S zcon
	
	variable ax = GetAX(xyz)
	variable by = GetBY(xyz)
	variable cz = GetCZ(xyz)
	variable rmax = max(ax, by)
	
	printf "Calculating 3D PSF . . .\r "
	Make/O/N=(nx, nx, nz) psf3D
	Setscale/I x -rmax/2, rmax/2, "", psf3D
	SetScale/I y -rmax/2, rmax/2, "", psf3D
	Setscale/I z 0, cz, "", psf3D
	stempsf3Dsym(psf3D, C5, Cs, df, 0, 0, keV, ap)
	wave/c psf3Dc = $"psf3D"
	
	printf "Calculating 3D scattering power . . .\r "
	stempot3D(xyz, nx, nx, nz)
	wave/c pot3D = $"pot3D"
	if(!waveexists(pot3D))
		printf "stempot has failed for some reason.  cannot continue.\r"
		return 0
	endif
	
	printf "FFTs . . .\r"
	fft psf3d
	fft pot3D

	printf "Multiplication of PSF by scattering . . .\r"
	pot3D *= psf3d

	printf "IFFT . . .\r"
	ifft pot3D
	imagetransform swap3D pot3D
	
	duplicate/o pot3D $curfol+"zcon_im"
	SetScale/P x DimOffSet(pot3D, 0), DimDelta(pot3D, 0), "", $curfol+"zcon_im"
	SetScale/P y DimOffSet(pot3D, 1), DimDelta(pot3D, 1), "", $curfol+"zcon_im"
	
	SetDataFolder $curfol
	Killwaves psf3d, pot3D
	printf "zcon3D done.\r"
end
