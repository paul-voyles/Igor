#pragma rtGlobals=1		// Use modern global access method.


Menu "Analysis"
	"Waves average . . .", average_sdom()
	help = {"Compute the average and standard deviation of the mean of a set of waves."}
end


function average_sdom() : Panel

	string fol = GetDataFolder(1)
	string wl = Wavelist("*", ";", "DIMS:1")
	variable nwaves = ItemsInList(wl)
	NewDataFolder/O/S root:packages
	NewDataFolder/O/S av_sdom
	Make/O/T/N=(nwaves) listwave
	Make/O/B/U/N=(nwaves) selwave
	variable i
	for(i=0; i<nwaves; i+=1)
		listwave[i] = StringFromList(i, wl, ";")
	endfor
	SetDataFolder fol

	// build the panel
	PauseUpdate; Silent 1		// building window...
	NewPanel/K=1/W=(440.25,44,776.25,410.75) as "Average and SDOM"
	SetDrawLayer UserBack
	DrawText 6,20,"Select data waves:"
	ListBox pickwaves,pos={5,23},size={265,341},frame=4,listWave=root:packages:av_sdom:listwave
	ListBox pickwaves,selWave=root:packages:av_sdom:selwave,mode= 4
	Button average,pos={280,23},size={50,20},proc=PanelAverageSDOM,title="Average"
	Button SDOM,pos={280,57},size={50,20},proc=PanelAverageSDOM,title="SDOM"
	Button both,pos={280,91},size={50,20},proc=PanelAverageSDOM,title="Both"


End


Function PanelAverageSDOM(ctrlName) : ButtonControl
	String ctrlName

	// locate the data from the listbox
	wave/t listwave = $"root:packages:av_sdom:listwave"
	wave selwave = $"root:packages:av_sdom:selwave"
	
	variable nlist = numpnts(listwave)
	variable nwave = sum(selwave, -inf, inf)
	
	// check to see if all the selected waves are the same length
	variable i, npts=-1
	for(i=0; i<nlist; i+=1)
		if(selwave[i])
			wave w = $listwave[i]
			if(npts == -1)
				npts = numpnts(w)
			else
				if(numpnts(w) != npts)
					DoAlert 0, "All the waves you select must be the same length."
					return 0
				endif
			endif
		endif
	endfor

	// build temporary and output waves
	Make/O/N=(nwave) temp1d
	Make/O/N=(npts, nwave) temp2d
	Duplicate/O w average
	Duplicate/O w sdom
	
	variable j=0
	for(i=0; i<nlist; i+=1)
		if(selwave[i])
			wave w = $listwave[i]
			temp2d[][j] = w[p]
			j+=1
		endif
	endfor
	
	for(i=0; i<npts; i+=1)
			temp1d = temp2d[i][p]
			wavestats/q temp1d
			average[i] = V_avg
			sdom[i] = V_sdev
	endfor
	sdom /= sqrt(nwave)
	
	// figure out which button was pushed
	strswitch (ctrlName)
		case "Average":
			Killwaves sdom
			break
		case "SDOM":
			Killwaves average
			break
		case "Both":
			break
		default:
			printf "How do you get here?  There's an error in PanelAverageSDOM . . .\r"
			return 0
	endswitch
	
	//Killwaves temp1d, temp2d

End
