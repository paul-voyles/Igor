#pragma rtGlobals=1		// Use modern global access method.


function LoadVorStatDir()

	NewPath/O/M="Select a folder of VP output _stat.out files." p
	PathInfo p
	if(!strlen(S_Path))
		return -1
	endif
	
	LoadVorStatDirWork(S_Path)
	
end


function LoadVorStatDirWork(infol)
	string infol
	
	NewPath/O p, infol
	PathInfo p

	string f, n
	
	variable i=0
	do
		f = IndexedFile(p, i, ".out")
		if(!strlen(f))
			break
		endif
		
		n = ParseFilePath(3, f, ":", 0, 0)
		if(stringmatch(n, "*_stat*"))
			LoadVorStatWork((S_Path+f))
			wave vor_stat = $"vor_stat"
			SortVor(vor_stat)
	
			n = ReplaceString("_vor_stat", n, "")
			rename vor_stat $(n+"_vstat")
			rename vor_stat_total $(n+"_vtot")
			rename vor_count $(n+"_vct")
		endif
		
		i+=1
	while(1)
		
end

function LoadVorStat()
	
	variable f
	Open/D/R/F="All Files:.*;" f
	if(!strlen(S_filename))
		return 0
	endif

	LoadVorStatWork(S_Filename)
	
end

function LoadVorStatWork(file)
	string file
	
	LoadWave/Q/G/A/B="N='_skip_'; N=total_n; N=center1_n; N=center2_n; N=center3_n; N=n3;N=n4; N=n5; N=n6; N=n7; N=n8; N=n9; N=n10; N=nbr_1; N=nbr_2; N=nbr_3;" file
	
	Make/o/n=(numpnts(total_n), 8) vor_stat
	Make/o/n=(numpnts(total_n), 15) vor_stat_total
	
	wave total_n = $"total_n"
	wave center1_n = $"center1_n"
	wave center2_n = $"center2_n"
	wave center3_n  =$"center3_n"
	wave n3 = $"n3"
	wave n4 = $"n4"
	wave n5 = $"n5"
	wave n6 = $"n6"
	wave n7 = $"n7"
	wave n8 = $"n8"
	wave n9 = $"n9"
	wave n10 = $"n10"
	wave nbr_1 = $"nbr_1"
	wave nbr_2 = $"nbr_2"
	wave nbr_3 = $"nbr_3"	
	
	Sort/R total_n, total_n, center1_n, center2_n, center3_n, n3, n4, n5, n6, n7, n8, n9, n10, nbr_1, nbr_2, nbr_3
	
	vor_stat[][0] = n3[p]
	vor_stat[][1] = n4[p]
	vor_stat[][2] = n5[p]
	vor_stat[][3] = n6[p]
	vor_stat[][4] = total_n[p]
	vor_stat[][5] = center1_n[p]
	vor_stat[][6] = center2_n[p]
	vor_stat[][7] = center3_n[p]

	vor_stat_total[][0] = n3[p]
	vor_stat_total[][1] = n4[p]
	vor_stat_total[][2] = n5[p]
	vor_stat_total[][3] = n6[p]
	vor_stat_total[][4] = total_n[p]
	vor_stat_total[][5] = center1_n[p]
	vor_stat_total[][6] = center2_n[p]
	vor_stat_total[][7] = center3_n[p]

	Killwaves total_n, center1_n, center2_n, center3_n, n3, n4, n5, n6, n7, n8, n9, n10, nbr_1, nbr_2, nbr_3
	
end

function SortVor(vor_stat)
	wave vor_stat
	
	make/o/n=4/t vor_types = {"icosahedral", "crystal-like", "mixed", "other"}
	make/o/n=(4, 4) vor_count
	vor_count = 0

	variable i, used, n3, n4, n5,n6, ico, cryst, mix
	for(i=0; i<DimSize(vor_stat, 0); i+=1)
		ico = 0
		cryst = 0
		mix = 0
		n3 = vor_stat[i][0]
		n4 = vor_stat[i][1]
		n5 = vor_stat[i][2]
		n6 = vor_stat[i][3]
		
		//perfect and distorted ico
		if(CVorNN6(n3, n4, n5, 0, 0, 12))
			ico = 1
		endif
		
		if(CVorNN6(n3, n4, n5, 0, 1, 10))
			ico = 1
		endif
		if(CVorNN6(n3, n4, n5, 0, 2, 8))
			ico = 1
		endif
		
		// crystal-like
		if(CVorNN6( n3, n4, n5, 0, 4, 4))
			cryst = 1
		endif
		if(CVOrNN6(n3, n4, n5, 0, 4, 5))
			cryst = 1
		endif
		if(CVorNN6(n3, n4, n5, 0, 5, 2))
			cryst = 1
		endif

		// mixed
		if(CVorNN6(n3, n4, n5, 0, 3, 6))
			mix = 1
		endif
		if(CVorNN6(n3, n4, n5, 0, 3, 7))
			mix = 1
		endif
		
		if(ico == 1)
			vor_count[0][0] += vor_stat[i][4]
			vor_count[0][1] += vor_stat[i][5]
			vor_count[0][2] += vor_stat[i][6]
			vor_count[0][3] += vor_stat[i][7]
		elseif(cryst == 1)
			vor_count[1][0] += vor_stat[i][4]
			vor_count[1][1] += vor_stat[i][5]
			vor_count[1][2] += vor_stat[i][6]
			vor_count[1][3] += vor_stat[i][7]
		elseif(mix == 1)
			vor_count[2][0] += vor_stat[i][4]
			vor_count[2][1] += vor_stat[i][5]
			vor_count[2][2] += vor_stat[i][6]
			vor_count[2][3] += vor_stat[i][7]
		else
			//printf "other: VP = <%d %d %d %d>, %d times\r", n3, n4,n5, n6, vor_stat[i][4]
			vor_count[3][0] += vor_stat[i][4]
			vor_count[3][1] += vor_stat[i][5]
			vor_count[3][2] += vor_stat[i][6]
			vor_count[3][3] += vor_stat[i][7]
		endif
			
	endfor	
	
	
	
end

function CVor(n3, n4, n5, n6, ref3, ref4, ref5, ref6)
	variable n3, n4, n5, n6, ref3, ref4, ref5, ref6
	
	if( n3 == ref3 && n4 == ref4 && n5 == ref5 && n6 == ref6 )
		return 1
	else 
		return 0
	endif
		
end

function CVorNN6(n3, n4, n5, ref3, ref4, ref5)
	variable n3, n4, n5, ref3, ref4, ref5
	
	if(n3 == ref3 && n4 == ref4 && n5 == ref5)
		return 1
	else
		return 0
	endif
	
end