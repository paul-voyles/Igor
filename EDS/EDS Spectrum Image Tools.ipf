#pragma rtGlobals=1 // Use modern global access method.


// Integrates a 3D spectrum image (SI) and outputs "IntegratedSpectrum" that is a 1D wave containg Energy vs Intensity.
function IntegrateSpectrum(SI)
	wave SI
	make/o/n=(dimsize(SI,2)) IntegratedSpectrum
	variable i=0
	for(i=0; i<dimsize(SI,2); i+=1)
		imagetransform/p=(i) getplane SI
		wave im =M_ImagePlane
		wavestats/Q im
		IntegratedSpectrum[i] = V_sum
		killwaves im
	endfor
	setscale/p x, DimOffset(SI, 2), DimDelta(SI,2), IntegratedSpectrum
end

// scans through a 3D wave and sets any pixel that is "inf" to zero.
// After NRR EDS spectrum images, some pixels were set to inf instead of zero. This will fix that issue so the spectrum image can be viewed.
function InfToZero(im_series)
	wave im_series
	variable x_size = DimSize(im_series, 0)
	variable y_size = DimSize(im_series, 1)
	variable z_size = DimSize(im_series, 2)
	
	variable value
	variable counter = 0
	
	variable i=0
	variable j=0
	variable k=0
	for(i=0; i<x_size; i+=1)
		j=0
		for(j=0; j<y_size; j+=1)
			k=0
			for(k=0; k<z_size; k+=1)

				value = im_series[i][j][k]
				if(value == inf)
					im_series[i][j][k] = 0
					counter = counter +1	
				endif
			endfor
		endfor
	endfor
end


// thresholds an integrated spectrum image by the given threshold value. 
// outputs a wave called spectrum_th that contains 0 where the integrated spectrum was less than threshold, and a 1 where the integrated spectrum is greater than the threshold.
function ThresholdSpectrum(IntegratedSpectrum, threshold)
	wave IntegratedSpectrum
	
	variable threshold
	
	duplicate/o IntegratedSpectrum spectrum_th
	
	variable int
	
	variable i=0
	for(i=0; i<2048; i+=1)
		int = IntegratedSpectrum[i]
		if(int < threshold)
			spectrum_th[i] = 0
		endif
		if(int >= threshold)
			spectrum_th[i] = 1
		endif
	endfor
	
	setscale/p x, -0.48, 0.01, spectrum_th
	
	Display/R spectrum_th; AppendToGraph IntegratedSpectrum
	ModifyGraph mode(spectrum_th)=5,rgb(spectrum_th)=(16385,49025,65535);DelayUpdate
	ModifyGraph rgb(IntegratedSpectrum)=(0,0,0)
	ModifyGraph nticks(bottom)=100
	ModifyGraph mode=7,hbFill=2
	ModifyGraph nticks(right)=0,axRGB(right)=(65535,65535,65535)
	
	Label bottom "Energy (keV)"
	Label left "# of Counts"
end


// truncates a 3D spectrum image. add more
function TruncateSpectrumImage(SI, peak_locations)
	wave SI, peak_locations
	
	variable start, finish
	variable x_size = DimSize(SI, 0)
	variable y_size = DimSize(SI, 1)
	
	make/o/n=(x_size,y_size, 0) SI_trunc
	
	variable i,j,k = 0
	for(i=0; i<dimsize(peak_locations,0); i+=1)
		j= abs(i+1-dimsize(peak_locations,0))
		start = peak_locations[j][0]
		finish = peak_locations[j][1]

		WindowSpectrumImage(SI, start, finish)
		wave im = window_image
		imagetransform/p=0/insw=im/o insertZplane SI_trunc
		killwaves im
	endfor
end



function WindowSpectrumImage(SI, start, finish)
	wave SI
	variable start, finish
	
	variable x_size = DimSize(SI, 0)
	variable y_size = DimSize(SI, 1)
	
	make/o/n=(x_size,y_size) window_image
	window_image = 0
	
	variable i=0
	for(i=start; i<finish+1; i+=1)
		imagetransform/p=(i) getplane SI
		wave im =M_ImagePlane
		
		window_image = window_image + im
		
		killwaves im
	endfor
end


function MultiWindowSpectrumImage(SI, spectrum_th, cutoff, color)
	wave SI, spectrum_th
	variable cutoff, color
	
	variable start, finish
	variable p1, p2
	
	string name, name2
	String name_list="name_list"
	
	variable count=0
	variable i=0
	variable j=0
	variable p,q=0
	for(i=0; i<2048; i+=1)
		p1 = spectrum_th[i]
		p2 = spectrum_th[i+1]
		
		if(p1==0 && p2==1)
			start = i+1
		endif
		
		if(p1==1 && p2==0)
			finish = i
			WindowSpectrumImage(SI, start, finish)
			wave im = window_image
			imagestats/Q im
			if(V_max >= cutoff)
				sprintf name, "SI_%g_%g", (start*0.01)-0.48, (finish*0.01)-0.48
				sprintf name2, "Window_%g_%g", (start*0.01)-0.48, (finish*0.01)-0.48
				
				Rename im $name
				NewImage/K=0  $name
				ModifyGraph nticks=0,axRGB=(65535,65535,65535)
				
				//add color schemes 
				if(color == 1)
					if(mod(count,8)==0)
						//ModifyImage $name,ctab= {*,*,Blue,1}
						Make/O/N=(1000, 3) colorwave1
						p=0
						q=0
						for(p=0; p<1000; p+=1)
							for(q=0; q<3; q+=1)
								if(q==0)
									colorwave1[p][q] = 0
								endif
								if(q==1)
									colorwave1[p][q] = 0
								endif
								if(q==2)
									colorwave1[p][q] = round((p/1000)*65535)
								endif
							endfor
						endfor
						setscale/I x, 0, V_max, colorwave1
						ModifyIMage $name, cindex=colorwave1
					endif	
					if(mod(count,8)==1)
						//ModifyImage $name,ctab= {*,*,Green,1}
						Make/O/N=(1000, 3) colorwave2
						p=0
						q=0
						for(p=0; p<1000; p+=1)
							for(q=0; q<3; q+=1)
								if(q==0)
									colorwave2[p][q] = 0
								endif
								if(q==1)
									colorwave2[p][q] = round((p/1000)*65535)
								endif
								if(q==2)
									colorwave2[p][q] = 0
								endif
							endfor
						endfor
						setscale/I x, 0, V_max, colorwave2
						ModifyIMage $name, cindex=colorwave2
					endif	
					if(mod(count,8)==2)
						//ModifyImage $name,ctab= {*,*,Red,1}
						Make/O/N=(1000, 3) colorwave3
						p=0
						q=0
						for(p=0; p<1000; p+=1)
							for(q=0; q<3; q+=1)
								if(q==0)
									colorwave3[p][q] = round((p/1000)*65535)
								endif
								if(q==1)
									colorwave3[p][q] = 0
								endif
								if(q==2)
									colorwave3[p][q] = 0
								endif
							endfor
						endfor
						setscale/I x, 0, V_max, colorwave3
						ModifyIMage $name, cindex=colorwave3
					endif	
					if(mod(count,8)==3)
						//ModifyImage $name,ctab= {*,*,Magenta,1}
						Make/O/N=(1000, 3) colorwave4
						p=0
						q=0
						for(p=0; p<1000; p+=1)
							for(q=0; q<3; q+=1)
								if(q==0)
									colorwave4[p][q] = round((p/1000)*65535)
								endif
								if(q==1)
									colorwave4[p][q] = 0
								endif
								if(q==2)
									colorwave4[p][q] = round((p/1000)*52428)
								endif
							endfor
						endfor
						setscale/I x, 0, V_max, colorwave4
						ModifyIMage $name, cindex=colorwave4
					endif	
					if(mod(count,8)==4)
						//ModifyImage $name,ctab= {*,*,Yellow,1}
						Make/O/N=(1000, 3) colorwave5
						p=0
						q=0
						for(p=0; p<1000; p+=1)
							for(q=0; q<3; q+=1)
								if(q==0)
									colorwave5[p][q] = round((p/1000)*65535)
								endif
								if(q==1)
									colorwave5[p][q] = round((p/1000)*65535)
								endif
								if(q==2)
									colorwave5[p][q] = 0
								endif
							endfor
						endfor
						setscale/I x, 0, V_max, colorwave5
						ModifyIMage $name, cindex=colorwave5
					endif	
					if(mod(count,8)==5)
						//ModifyImage $name,ctab= {*,*,Cyan,1}
						Make/O/N=(1000, 3) colorwave6
						p=0
						q=0
						for(p=0; p<1000; p+=1)
							for(q=0; q<3; q+=1)
								if(q==0)
									colorwave6[p][q] = 0
								endif
								if(q==1)
									colorwave6[p][q] = round((p/1000)*65535)
								endif
								if(q==2)
									colorwave6[p][q] = round((p/1000)*65535)
								endif
							endfor
						endfor
						setscale/I x, 0, V_max, colorwave6
						ModifyIMage $name, cindex=colorwave6
					endif	
					if(mod(count,8)==6)
						//ModifyImage $name,ctab= {*,*,Grays256,1}
						Make/O/N=(1000, 3) colorwave7
						p=0
						q=0
						for(p=0; p<1000; p+=1)
							for(q=0; q<3; q+=1)
								if(q==0)
									colorwave7[p][q] = round((p/1000)*43690)
								endif
								if(q==1)
									colorwave7[p][q] = round((p/1000)*43690)
								endif
								if(q==2)
									colorwave7[p][q] = round((p/1000)*43690)
								endif
							endfor
						endfor
						setscale/I x, 0, V_max, colorwave7
						ModifyIMage $name, cindex=colorwave7
					endif	
					if(mod(count,8)==7)
						//ModifyImage $name,ctab= {*,*,Gold,1}
						Make/O/N=(1000, 3) colorwave8
						p=0
						q=0
						for(p=0; p<1000; p+=1)
							for(q=0; q<3; q+=1)
								if(q==0)
									colorwave8[p][q] = round((p/1000)*65535)
								endif
								if(q==1)
									colorwave8[p][q] = round((p/1000)*21845)
								endif
								if(q==2)
									colorwave8[p][q] = 0
								endif
							endfor
						endfor
						setscale/I x, 0, V_max, colorwave8
						ModifyIMage $name, cindex=colorwave8
					endif	
				endif
				
				make/o/n=(2048) $name2=0
				wave temp = $name2
				temp=0
				for(j=start; j<finish+1; j+=1)
					temp[j] = 1
				endfor
				setscale/p x, 0, 10, temp
				
				sprintf name_list, "%s;%s", name_list, name2
				
				print name, "max=", V_max, "min=", V_min, "avg=", V_avg, "stdev=", V_sdev 
				count=count+1
			endif
		endif
	endfor
	
	killwaves window_image
	
	setscale/P x, -0.48, 0.01, IntegratedSpectrum
	Display IntegratedSpectrum
	ModifyGraph rgb(IntegratedSpectrum)=(0,0,0)
	ModifyGraph nticks(bottom)=100
	
	string name3
	variable k=0
	for(k=0; k<count; k+=1)
		name3 = StringFromList(k+1, name_list)	
		setscale/P x, -0.48, 0.01, $name3
		AppendToGraph/R $name3
		ModifyGraph mode($name3)=7,hbFill($name3)=4
		
		if(color==0)
			ModifyGraph rgb($name3)=(65535,0,0)
		endif
		
		if(color==1)
			if(mod(k,8)==0)
				ModifyGraph rgb($name3)=(0,0,65535)
			endif	
			if(mod(k,8)==1)
				ModifyGraph rgb($name3)=(0,65535,0)
			endif	
			if(mod(k,8)==2)
				ModifyGraph rgb($name3)=(65535,0,0)
			endif	
			if(mod(k,8)==3)
				ModifyGraph rgb($name3)=(65535,0,52428)
			endif	
			if(mod(k,8)==4)
				ModifyGraph rgb($name3)=(65535,65535,0)
			endif	
			if(mod(k,8)==5)
				ModifyGraph rgb($name3)=(0,65535,65535)
			endif	
			if(mod(k,8)==6)
				ModifyGraph rgb($name3)=(43690,43690,43690)
			endif	
			if(mod(k,8)==7)
				ModifyGraph rgb($name3)=(65535,21845,0)
			endif	
		endif
		
	endfor
	RemoveFromGraph IntegratedSpectrum
	AppendToGraph IntegratedSpectrum
	ModifyGraph rgb(IntegratedSpectrum)=(0,0,0)
	ModifyGraph nticks(right)=0,axRGB(right)=(65535,65535,65535)
	ModifyGraph mode(IntegratedSpectrum)=7,hbFill(IntegratedSpectrum)=2
	SetAxis bottom 0,10
	Label bottom "Energy (keV)"
	Label left "# of Counts"
	
end


//subtract calculated continuum background from experimental EDS spectrum image
function SubContinuumBackground(im, a, b, c)
	wave im
	variable a, b, c
	
	variable offset = DimOffset(im, 2) 
	variable dispersion = DimDelta(im,2)
	variable x, sum_int
	variable i, j, k
	
	for(i=0; i<DimSize(im,2); i+=1)
		x=i-(abs(offset)/dispersion)
		//x = offset + (i*dispersion)
		if( i > (abs(offset)/dispersion))
			sum_int = a*(exp(-b/x))*((x-c)/x)
			for(j=0; j<DimSize(im,0); j+=1)
				for(k=0; k<DimSize(im,1); k+=1)
					im[j][k][i] = im[j][k][i] - (sum_int/(DimSize(im,0)*DimSize(im,1)))
				endfor
			endfor
		endif
	endfor
	
end
