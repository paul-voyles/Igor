#pragma rtGlobals=1		// Use modern global access method.

function LoadStandardSTEMProfileWork(fname, n)
	string fname, n
	
	string col
	sprintf col, "N=%s_xp;N=%s_yp;N=%s_bf;N=%s_mr15_30;N=%s_mr30_50;N=%s_mr30_100;", n,n,n,n,n,n
	sprintf col, "%sN=%s_mr50_250;N=%s_mr75_250;N=%s_mr100_250;N=%s_mr50_500;N=%s_mr75_500;", col, n,n,n,n,n
	sprintf col, "%sN=%s_mr100_500;N=%s_mr125_500;N=%s_mr150_500;N=%s_mr175_500;", col, n,n,n,n
	
	LoadWave/A/O/G/B=col fname

end

function LoadStandardSTEMProfile(n)
	string n
	
	variable f
	Open/D/R/F="Data Files:.dat;" f
	
	if(!strlen(S_filename))
		return 0
	endif
	
	LoadStandardSTEMProfileWork(S_filename, n)
	
end	

function LoadStandardSTEMProfileFolder()

	NewPath/O/M="Select a directory of STEM profile files" dirpath
	PathInfo dirpath
	string dirname = S_Path
	if(!strlen(dirname))
		return 0
	endif
	
	string fname,fpath
	variable i=0
	do
		fname = IndexedFile(dirpath, i, ".dat")
		if(!strlen(fname))
			return 1
		endif 
		fpath = dirname  + fname
		fname = StringFromList(0, fname, ".")
		LoadStandardSTEMProfileWork(fpath, fname)
		i+=1
	while(1)

end