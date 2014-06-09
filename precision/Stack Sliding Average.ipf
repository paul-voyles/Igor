#pragma rtGlobals=3		// Use modern global access method and strict wave access.

// function to do full-image sliding-window average through a 3D data stack
// Intended for de-noising by averaging of aligned HRSTEM images.
//
// 06-09-14 begun pmv

// st is input 3D data stack.  win is the width of the sliding average window.
// It must be an odd number >= 3. The output st_av will contain st, averaged
// over a sliding window from -win/2 to +win/2.  
function StackSlidingAverage(st, win)
	wave st
	variable win
	
	if(win < 3)
		printf "Stack average only work for window of at least 3.\r"
		win = 3
	endif
	
	if(!mod(win, 2))
		printf "Stack average is well-defined only for odd window size.\r"
		win -= 1
		printf "Resetting window size to %d.\r", win
	endif
	variable half_win = (win-1)/2
	
	duplicate/O st st_av
	make/O/N=( Dimsize(st, 0), Dimsize(st, 1), win) win_av
	
	variable imin, imax, i, j
	imin = half_win
	imax = DimSize(st, 2) - half_win
	
	for(i=imin; i<imax; i+=1)
		win_av = st[p][q][r+i-half_win]
		Imagetransform averageImage win_av
		wave one_av = $"M_AveImage"
		Imagetransform/P=(i)/D=one_av setplane st_av
	endfor
	
	make/O/N=( Dimsize(st, 0), Dimsize(st, 1), half_win+1) win_av
	for(i=0; i<imin; i+=1)
		win_av = st[p][q][r+i]
		MatrixOp/O one_av =sumbeams(win_av)
		one_av /= (half_win+1)
		Imagetransform/P=(i)/D=one_av setplane st_av
	endfor
	
	for(i=imax; i<DimSize(st, 2); i+=1)
		win_av = st[p][q][r+i-half_win]
		MatrixOp/O one_av =sumbeams(win_av)
		one_av /= (half_win+1)
		Imagetransform/P=(i)/D=one_av setplane st_av
	endfor
	
	Killwaves $"M_AveImage", $"M_StdvImage", win_av
end
	