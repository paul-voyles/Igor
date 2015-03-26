#pragma rtGlobals=1		// Use modern global access method.

function WaveNameChange(image)
	
	wave image
	
	string S_name
	S_name = NameOfWave(image)
	
	String expr="([[:alpha:]]+)_([[:alpha:]]+)_([[:digit:]]+)_([[:alpha:]]+)"
	String base1, base2, number, ending
	SplitString/E=(expr) S_name, base1, base2, number, ending
	//print base1
	//print base2
	//print number
	//print ending

	String S_name_new 
  	S_name_new = base1 + "_" + base2 + "_" + ending + "_" + number
	//print S_name_new
	
	Rename image $S_name_new
	
end


function AllWaveNameChange(startimage, total)
	
	wave startimage
	variable total
	
	string S_name
	S_name = NameOfWave(startimage)
	
	String expr="([[:alpha:]]+)_([[:alpha:]]+)_([[:digit:]]+)_([[:alpha:]]+)"
	String base1, base2, number, ending
	SplitString/E=(expr) S_name, base1, base2, number, ending
	
	variable num
	num = str2num(number)
	//print num
	
	variable i=num
	for(i=num; i<total+1; i+=1)
	 
		String S_num
		S_num = num2str(i)
		//print S_num
	
		String S_name_call
		if (i>=0 && i<=9)
			S_name_call = base1 + "_" + base2 + "_00" + S_num + "_" + ending
		elseif (i>=10 && i<=99)
			S_name_call = base1 + "_" + base2 + "_0" + S_num + "_" + ending
		elseif (i>=100 && i<=999)
			S_name_call = base1 + "_" + base2 + "_" + S_num + "_" + ending
		endif
		//print S_name_call
	 
		wave image = $S_name_call
		WaveNameChange(image)
		
	endfor
		
end



function SideBySideStacks(Stack1,Stack2)
	
	wave Stack1
	wave Stack2
	
	variable sizex1 = DimSize(Stack1, 0)
	variable sizex2 = DimSize(Stack2, 0)
	variable sizey1 = DimSize(Stack1, 1)
	variable sizey2 = DimSize(Stack2, 1)
	variable sizez1 = DimSize(Stack1, 2)
	variable sizez2 = DimSize(Stack2, 2)
	
	If (sizey1 != sizey2)
		print "y dimensions don't match"
	endif
	
	If (sizez1 != sizez2)
		print "z dimensions don't match"
	endif
	
	Make/o/d/n=(sizex1+sizex2,sizey1,sizez1) Stack_SbS
	
	variable i=0
	for(i=0; i<sizez1; i+=1)
		Imagetransform/p=(i) getplane Stack1
		wave image1 = $"M_ImagePlane"
		Imagetransform/p=(i) getplane Stack2
		wave image2 = $"M_ImagePlane"
		
		imagetransform /INSI=image1 /INSX=0 /INSY=0 /P=(i) insertImage Stack_SbS
		imagetransform /INSI=image2 /INSX=(sizex1) /INSY=0 /P=(i) insertImage Stack_SbS
	endfor
	killwaves M_ImagePlane
end


function TopAndBottomStacks(Stack1,Stack2)
	
	wave Stack1
	wave Stack2
	
	variable sizex1 = DimSize(Stack1, 0)
	variable sizex2 = DimSize(Stack2, 0)
	variable sizey1 = DimSize(Stack1, 1)
	variable sizey2 = DimSize(Stack2, 1)
	variable sizez1 = DimSize(Stack1, 2)
	variable sizez2 = DimSize(Stack2, 2)
	
	If (sizex1 != sizex2)
		print "x dimensions don't match"
	endif
	
	If (sizez1 != sizez2)
		print "z dimensions don't match"
	endif
	
	Make/o/d/n=(sizex1,sizey1+sizey2,sizez1) Stack_TaB
	
	variable i=0
	for(i=0; i<sizez1; i+=1)
		Imagetransform/p=(i) getplane Stack1
		wave image1 = $"M_ImagePlane"
		Imagetransform/p=(i) getplane Stack2
		wave image2 = $"M_ImagePlane"
		
		imagetransform /INSI=image1 /INSX=0 /INSY=0 /P=(i) insertImage Stack_TaB
		imagetransform /INSI=image2 /INSX=0 /INSY=(sizey1) /P=(i) insertImage Stack_TaB
	endfor
	killwaves M_ImagePlane
end