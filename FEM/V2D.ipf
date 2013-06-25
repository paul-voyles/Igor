#pragma rtGlobals=1		// Use modern global access method.

//change log
//03/28/11 function V2D_manual was added. It calculates V2D without using sumbeams.  
//             Sumbeams include NaN values in the calculation.  V2D_manual does not include NaNs in the calculation. -JWH


function V2D(s)
	wave s
	
	variable np = 1.0/DimSize(s, 2)

	Duplicate/O s v2d_t
	fastop v2d_t = v2d_t*v2d_t
	
	matrixop/o i2_2d = sumbeams(v2d_t)
	i2_2d = i2_2d*np
	matrixop/o i_2d = sumbeams(s)
	i_2d = np*i_2d
	Make/O/N=(DimSize(s, 0), Dimsize(s,1)) v_2d
	v_2d = i2_2d / (i_2d*i_2d) - 1.0
	
	duplicate/O i2_2d $(NameofWave(s)+"_i2_2d")
	duplicate/O i_2d $(NameofWave(s) + "_i_2d")
	duplicate/O v_2d $(NameofWave(s) + "_v_2d")
	
	Killwaves v2d_t, i2_2d, i_2d, v_2d
	
	
	
end	


function V2D_manual(dat)
	wave dat
	
	variable i, j, k
	variable i_sum, i_count,  i2_sum, i2_count
	
	make/o/d/n=(DimSize(dat, 0), DimSize(dat, 1)) i_avg, i2_avg, var
	
	var = 0
	
	for(i=0; i<DimSize(dat, 0); i+=1)
		for(j=0; j<DimSize(dat, 1); j+=1)	
			i_sum = 0
			i_count = 0
			i2_sum = 0
			i2_count = 0
			for(k=0; k<DimSize(dat, 2); k+=1)
				if(dat[i][j][k] == dat[i][j][k])
					i_sum += dat[i][j][k]
					i_count +=1
					i2_sum += dat[i][j][k]*dat[i][j][k]
					i2_count +=1
				endif	
			endfor
			i_avg[i][j] = i_sum / i_count
			i2_avg[i][j] = i2_sum / i2_count
			var[i][j] = i2_avg[i][j] / (i_avg[i][j] * i_avg[i][j]) -1
		endfor
	endfor
	
	
	
	print i_sum, i_count
	
	   

	
	
	
end
