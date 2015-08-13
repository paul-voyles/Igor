#pragma rtGlobals=1		// Use modern global access m--ethod.
#include "SER File Loader"
#include "Annular Average2"
#include "DP Annular Average UI v1.3"
#include "Auto Find DP Center2"
#include "Average and SDOM"
#include "GFX"
//#include "Stem_Var_filter_bin8"

// The function "correlatefunc_k" is to calculate time autocorrelatlion g2(t) at a reciprocal space vector k from diffraction patttern time series.
            // data: diffraction pattern time series.
            // allowance: small change in k. Find all points from diffraction patterns in the range (k-allowance, k+allowance).
            // (xcenter, ycenter): diffraction pattern center position.
            // xblockstart, xblockend, yblockstart, yblockend: the location of beam stopper in the diffraction pattern.
            
// The function "relaxationfunc" is to fit g2(t) to 1+B*exp(-2*(t/tau)^beta).
// by Li He, updated on 01-22-15.
//annular_average(data, xcenter, ycenter, xblockstart, xblockend, yblockstart, yblockend, strip_width). // Calculate annular average. Titan CCD.
//correlatefunc_k_map(data,k,bin,xcenter,ycenter,xblockstart,xblockend,yblockstart,yblockend) // Calculate g2, correlation function.
//resampleg2(g2, name) // Re-sample g2 to be evenly distributed in log time scale.

Function correlatefunc_allPixels(im_stack, x1, x2, y1, y2, name)    // Calculate g2 from STEM image series to detect, e.g., BSD relaxation time.
																	// x1, y1, x2, y2 are selected regions.
	       wave im_stack
	       variable x1, y1, x2, y2	
	       string name      
	       
	       variable tp = DimSize(im_stack, 2), samplesize=(x2-x1)*(y2-y1)
	
              make/o/n=(tp) add_td=0
              make/o/n=(tp, samplesize) g2t_tpt_map=0
              make/o/n=(samplesize) I_ave, xmap, ymap
               
	       variable x, y, i, j, ave, pcount=0
	       string g2t_tpt_mapname="g2tmap_tpt_"+name
	       string g2t_tpt_name="g2t_tpt_"+name
	       string I_avename="I_ave_"+name
	       string xmapname="x_"+name
	       string ymapname="y_"+name
	                                                                                                                                                     

	       for(x=x1; x<x2; x+=1)
	       	for(y=y1; y<y2; y+=1)
	       		Imagetransform/beam={(x), (y)} getbeam im_stack
				wave W_beam = $"W_beam"
				ave=mean(W_beam)
				add_td=0
				for (j=0;j<tp;j+=1)
					for (i=0;i<(tp-j);i+=1)
						add_td[j]=add_td[j]+W_beam[i]*W_beam[i+j]
			      		endfor
				endfor

				g2t_tpt_map[][pcount]=(tp-p)*add_td[p]/sum(W_beam,0,tp-p-1)/sum(W_beam,p,tp-1)
			 	I_ave[pcount]=ave
			 	xmap[pcount]=x
			 	ymap[pcount]=y
			 	pcount=pcount+1
			 endfor
             endfor
             
             SetScale/P x 0, 0.948,"", g2t_tpt_map // STEM image series, 0.84 s per frame, total 0.948 s per frame. 
              
             redimension/n=(-1, pcount) g2t_tpt_map
             redimension/n=(pcount) I_ave, xmap, ymap
             imagetransform sumallrows g2t_tpt_map
             duplicate/o W_sumRows g2t
             g2t=g2t/pcount
             SetScale/P x 0,0.948,"", g2t
              
             duplicate/o g2t_tpt_map, $g2t_tpt_mapname
             duplicate/o g2t, $g2t_tpt_name

             duplicate/o I_ave, $I_avename
             duplicate/o xmap, $xmapname
             duplicate/o ymap, $ymapname
                
             killwaves g2t_tpt_map
             killwaves g2t, I_ave, add_td, xmap, ymap
             killwaves W_beam, W_sumrows
             print pcount
end

Function multi_tauBNL(data, k, bin, xcenter, ycenter, xa, ya, yb, xc, yc, xd, xe, ye)    // Multi-tau, the first segment is 1/16 total data length.
// BNL K2 data
	       wave data
	       variable k, bin, xcenter, ycenter, xa, ya, yb, xc, yc, xd, xe, ye
	      
	       //variable r=k/0.008052248 // Bin 8. Camera length 380mm.
	       variable r=k/0.01006531 // Bin 10. Camera length 380mm.
	       variable tp = DimSize(data, 2), neighbour=round((bin-1)/2)
	       variable samplesize=round(2*pi*r)
	       variable kbc=(yc-yb)/(xc-xa), kde=(ye-yc)/(xe-xd)
	
               make/o/n=(tp) add_td=0
               make/o/n=(tp,samplesize) g2t_tpt_map=0
               make/o/n=(samplesize) x_map, y_map, I_ave
               make/o/n=(samplesize) theta, xp, yp
               theta[]=p*2*pi/samplesize
               xp=round(xcenter+r*cos(theta))
               yp=round(ycenter+r*sin(theta))
               
	       variable i, j, sumi, sumip
	       variable pcount=0
	       string name="_"+nameofwave(data)+"_"+num2str(k)+"_bin"+num2str(bin)
	       string g2t_tpt_mapname="g2tmap_tpt"+name
	       string g2t_tpt_name="g2t_tpt"+name
	       string x_mapname="x_map"+name
	       string y_mapname="y_map"+name
	       string I_avename="I_ave"+name
	                                                                                                                                                               // Find the circular points on diffraction patterns with the reciprocal space vector k.

	       for(k=0; k<samplesize; k+=1)
		        // Remove beam stop shadow
		        if (xp[k]>xc && xp[k]<xd && yp[k]>ya && yp[k]<yc)
		              pcount=pcount
		       elseif (xp[k]>xa && xp[k]<xc && yp[k]>ya && yp[k]<yb)
		       	      pcount=pcount
		       elseif (xp[k]>xa && xp[k]<xc && yp[k]>yb && yp[k]<kbc*(xp[k]-xa)+yb)
		       	      pcount=pcount
		       elseif (xp[k]>xd && xp[k]<xe && yp[k]>ya && yp[k]<ye)
		              pcount=pcount
		        elseif (xp[k]>xd && xp[k]<xe && yp[k]>ye && yp[k]<kde*(xp[k]-xd)+yc)
		              pcount=pcount
		       elseif(xp[k]>xa && xp[k]<xe && yp[k]>ya && yp[k]<yb)
		        		pcount=pcount
		        else
		              make/o/n=(tp) bin_average=0
		              for(i=0-neighbour; i<neighbour+1; i+=1)
		                     for(j=0-neighbour; j<neighbour+1; j+=1)
			                        Imagetransform/beam={(xp[k]+i), (yp[k]+j)} getbeam data
			                        wave W_beam = $"W_beam"
			                        bin_average=bin_average+W_beam
			            endfor
			      endfor
			      W_beam=bin_average/bin^2
			      
			      add_td=0
			      for (j=0;j<tp/16;j+=1)
			            for (i=0;i<(tp/16-j);i+=1)
			                      add_td[j]=add_td[j]+W_beam[i]*W_beam[i+j]
			            endfor
			      endfor

			      g2t_tpt_map[0, tp/16-1][pcount]=(tp/16-p)*add_td[p]/sum(W_beam,0,tp/16-p-1)/sum(W_beam,p,tp/16-1)
			      
			      add_td=0
			      sumi=0
			      sumip=0
			      for (j=tp/16; j<tp/8; j+=1)
			      		for (i=0; i<(tp/8-j); i+=2)
			      			add_td[j]=add_td[j]+W_beam[i]*W_beam[i+j]
			      			sumi=sumi+W_beam[i]
			      			sumip=sumip+W_beam[i+j]
			      		endfor
			      	endfor
			      	
			      	g2t_tpt_map[tp/16, tp/8-1][pcount]=(tp/8-p)/2*add_td[p]/sumi/sumip
			      	
			      	add_td=0
			      sumi=0
			      sumip=0
			      for (j=tp/8; j<tp/4; j+=1)
			      		for (i=0; i<(tp/4-j); i+=4)
			      			add_td[j]=add_td[j]+W_beam[i]*W_beam[i+j]
			      			sumi=sumi+W_beam[i]
			      			sumip=sumip+W_beam[i+j]
			      		endfor
			      	endfor
			      
			      g2t_tpt_map[tp/8, tp/4-1][pcount]=(tp/4-p)/4*add_td[p]/sumi/sumip
			      	
			      	add_td=0
			      sumi=0
			      sumip=0
			      for (j=tp/4; j<tp/2; j+=1)
			      		for (i=0; i<(tp/2-j); i+=8)
			      			add_td[j]=add_td[j]+W_beam[i]*W_beam[i+j]
			      			sumi=sumi+W_beam[i]
			      			sumip=sumip+W_beam[i+j]
			      		endfor
			      	endfor
			      
			      g2t_tpt_map[tp/4, tp/2-1][pcount]=(tp/2-p)/8*add_td[p]/sumi/sumip
			      
			      add_td=0
			      sumi=0
			      sumip=0
			      for (j=tp/2; j<tp; j+=1)
			      		for (i=0; i<(tp-j); i+=16)
			      			add_td[j]=add_td[j]+W_beam[i]*W_beam[i+j]
			      			sumi=sumi+W_beam[i]
			      			sumip=sumip+W_beam[i+j]
			      		endfor
			      	endfor
			      
			      g2t_tpt_map[tp/2, tp-1][pcount]=(tp-p)/16*add_td[p]/sumi/sumip
			      
			      x_map[pcount]=xp[k]
			      y_map[pcount]=yp[k]
			      I_ave[pcount]=mean(W_beam)
			      
                               pcount=pcount+1
                        endif
                 endfor
                
                redimension/n=(-1, pcount) g2t_tpt_map
                redimension/n=(pcount) x_map, y_map, I_ave
                SetScale/P x 0, 0.0025,"", g2t_tpt_map        //K2, 400 frames/s
                 
                imagetransform sumallrows g2t_tpt_map
                duplicate/o W_sumRows g2t
                g2t=g2t/pcount
                SetScale/P x 0,0.0025,"", g2t
                
                duplicate/o g2t_tpt_map, $g2t_tpt_mapname
                duplicate/o g2t, $g2t_tpt_name
                
                duplicate/o x_map, $x_mapname
                duplicate/o y_map, $y_mapname
                duplicate/o I_ave, $I_avename
                
                killwaves g2t_tpt_map, g2t, x_map, y_map, I_ave, add_td, bin_average, W_beam, W_sumrows, xp, yp, theta
                print pcount
end

Function correlatefunc_imageSeries(im_stack, xgrain, ygrain, x1, x2, y1, y2, name)    // Calculate g2 from STEM image series to detect, e.g., BSD relaxation time. xgrain, ygrain is selected pixels array.
																	// x1, y1, x2, y2 are selected regions.
	       wave im_stack, xgrain, ygrain
	       variable x1, y1, x2, y2	
	       string name      
	       
	       variable tp = DimSize(im_stack, 2)
	       variable samplesize=dimsize(xgrain, 0)
	
              make/o/n=(tp) add_td=0
              make/o/n=(tp,samplesize) g2t_tpt_map=0
              make/o/n=(samplesize) I_ave, xmap, ymap
               
	       variable i, j, k, ave, pcount=0
	       string g2t_tpt_mapname="g2tmap_tpt_"+name
	       string g2t_tpt_name="g2t_tpt_"+name
	       string I_avename="I_ave_"+name
	       string xmapname="x_"+name
	       string ymapname="y_"+name
	                                                                                                                                                     

	       for(k=0; k<samplesize; k+=1)
	       	if(xgrain[k]>x1&&xgrain[k]<x2&&ygrain[k]>y1&&ygrain[k]<y2)
	       		Imagetransform/beam={(xgrain[k]), (ygrain[k])} getbeam im_stack
				wave W_beam = $"W_beam"
				ave=mean(W_beam)
				add_td=0
				for (j=0;j<tp;j+=1)
					for (i=0;i<(tp-j);i+=1)
						add_td[j]=add_td[j]+W_beam[i]*W_beam[i+j]
			      		endfor
				endfor

				g2t_tpt_map[][pcount]=(tp-p)*add_td[p]/sum(W_beam,0,tp-p-1)/sum(W_beam,p,tp-1)
			 	I_ave[pcount]=ave
			 	xmap[pcount]=xgrain[k]
			 	ymap[pcount]=ygrain[k]
			 	pcount=pcount+1
			 endif
             endfor
             
             SetScale/P x 0, 0.948,"", g2t_tpt_map // STEM image series, 0.84 s per frame, total 0.948 s per frame. 
              
             redimension/n=(-1, pcount) g2t_tpt_map
             redimension/n=(pcount) I_ave, xmap, ymap
             imagetransform sumallrows g2t_tpt_map
             duplicate/o W_sumRows g2t
             g2t=g2t/pcount
             SetScale/P x 0,0.948,"", g2t
              
             duplicate/o g2t_tpt_map, $g2t_tpt_mapname
             duplicate/o g2t, $g2t_tpt_name

             duplicate/o I_ave, $I_avename
             duplicate/o xmap, $xmapname
             duplicate/o ymap, $ymapname
                
             killwaves g2t_tpt_map
             killwaves g2t, I_ave, add_td, xmap, ymap
             killwaves W_beam, W_sumrows
             print pcount
end

Function intensity_dose(av_series, indexes, name)
	wave av_series, indexes
	string name
	
	variable N=numpnts(indexes), i
	make/o/n=(N) 'sum_0.4_0.49', 'norm_0.4_0.49', 'sum_0.59_0.67', 'norm_0.59_0.67', 'sum_0.69_0.77', 'norm_0.69_0.77'
	make/o/n=(N) time_s, total_dose
	
	for(i=0; i<N; i+=1)
		duplicate/o av_series av_tmp
		redimension/n=(-1, indexes[i]) av_tmp
		
		imagetransform sumallrows av_tmp
		wave W_sumRows=$"W_sumRows"
		SetScale/P x 0,0.006967,"", W_sumRows
		duplicate/o W_sumRows av_index
		av_index=W_sumRows/indexes[i]
		
		'sum_0.4_0.49'[i]=sum(av_index, 0.40, 0.49)
		'sum_0.59_0.67'[i]=sum(av_index, 0.59, 0.67)
		'sum_0.69_0.77'[i]=sum(av_index, 0.69, 0.77)
		
		string av_index_name=name+"_"+num2str(indexes[i])
		duplicate/o av_index $av_index_name
	endfor
	
	sort indexes, indexes, 'sum_0.4_0.49', 'sum_0.59_0.67', 'sum_0.69_0.77'
	'norm_0.4_0.49'='sum_0.4_0.49'/'sum_0.4_0.49'[0]
	'norm_0.59_0.67'='sum_0.59_0.67'/'sum_0.59_0.67'[0]
	'norm_0.69_0.77'='sum_0.69_0.77'/'sum_0.69_0.77'[0]
	
	time_s[]=0.2363*indexes[p]
	total_dose[]=6.6e-12/1.6e-19*time_s[p]/pi*4/1.6^2
	
	killwaves av_tmp, av_index, W_sumRows
end
		
Function taufromMap(g2t_map, name)
		wave g2t_map
		string name
		
		variable N=dimsize(g2t_map, 1)
		variable V_flag, V_LevelX, count=0, total=0, i //V_maxloc
		string tauname="tau_"+name
		
		make/o/n=(N) tau
		duplicate/o g2t_map g2t
		redimension/n=(-1,0) g2t
		
		for(i=0; i<N; i+=1)
			V_LevelX=nan
			g2t[]=g2t_map[p][i]
			duplicate/o g2t g2t_smth
			Smooth 50, g2t_smth
			//wavestats/q g2t
			findlevel/edge=2/q g2t_smth, 1+(g2t[0]-1)*exp(-2)
			if(V_flag==0 && V_LevelX>0.2363*3 && V_LevelX<200)
				count=count+1
				//if(V_LevelX<1)
					//string fitname=name+"_"+num2istr(i)
					//string parametername="parameter_"+name+"_"+num2istr(i)
					//fitg2_noGraph(g2t, 15, fitname)
					//wave parameter=$parametername
					//tau[i]=parameter[0]
				//else
				tau[i]=V_LevelX
				//endif
				total=total+V_LevelX
				else
					tau[i]=nan
			endif
		endfor
		
		
		print count
		print total/count
		duplicate/o tau $tauname
		killwaves tau, g2t, g2t_smth//, parameter
end
			
Function g2fromMap(g2t_map, name)
		wave g2t_map
		string name
		
		variable n=dimsize(g2t_map, 1)
		string g2tname="g2t_"+name
		
		imagetransform sumallrows g2t_map
		duplicate/o W_sumRows g2t
		g2t[]=g2t[p]/n
		SetScale/P x 0,0.2363,"", g2t
		
		duplicate/o g2t $g2tname
		
		killwaves g2t
end

Function I_points(data, x, y, name)
               wave data
               variable x, y
               string name
               
               variable zp = DimSize(data, 2)
               string Iname="I_"+name+"_"+num2istr(x)+"_"+num2istr(y)
                       
               redimension/s data
               Imagetransform/beam={(x), (y)} getbeam data
               wave W_beam = $"W_beam"
               SetScale/P x 0, 0.2276,"", W_beam
               duplicate/o W_beam $Iname
               
               killwaves W_beam
end

function resampleg2(g2, name, samplepts)
		wave g2
		string name
		variable samplepts
		
		string g2_rsname="g2t_rs_"+name
		
		make/o/n=(samplepts) g2_rs, t_rs		
		//make/o/n=(100) g2_rs, t_rs
		
		variable datalength=dimsize(g2, 0)
		t_rs[]=alog(log(0.948)+p/(samplepts-1)*(log(0.948*datalength)-log(0.948))) // BSD image series, 0.84 s per frame, totally 0.948 s per frame.
		//t_rs=alog(log(1.528)+p/(samplepts-1)*(log(1.528*50)-log(1.528))) // Cluster size map, interpixel distance 1.528 nm.
		//t_rs=alog(log(0.248587083)+p/19*(log(0.248587083*1000)-log(0.248587083)))
		//t_rs=alog(log(0.248587083)+p/99*(log(0.248587083*1000)-log(0.248587083))) // 0.18s/frame, 1000 frames, bin8, 200^2 pixels.
		//t_rs=alog(log(0.0025)+p/(samplepts-1)*(log(0.0025*48000)-log(0.0025))) // K2 camera, 0.0025fps, 2 min.
		//t_rs=alog(log(0.267647)+p/99*(log(0.267647*500)-log(0.267647))) // 0.2s/frame, 1000 frames, bin8, 200^2 pixels. g2 half t length.
		//t_rs=alog(log(0.26547)+p/99*(log(0.26547*1000)-log(0.26547))) // 0.2s, 350^2 pixels, search mode.
		//t_rs=alog(log(2)+p/99*(log(2*50)-log(2))) // Cluster size map, interpixel distance 2 nm.
		
		g2_rs[]=g2[t_rs[p]/0.948]
		//g2_rs[]=g2[t_rs[p]/1.528]
		//g2_rs[]=g2[t_rs[p]/0.248587083]
		//g2_rs[]=g2[t_rs[p]/0.0025]
		//g2_rs[]=g2[t_rs[p]/0.26547]
             //g2_rs[]=g2[t_rs[p]/2]
             
		duplicate/o g2_rs, $g2_rsname
		killwaves g2_rs
end

Function filterlongtau(g2t_map, name)
		wave g2t_map
		string name
		
		variable N=dimsize(g2t_map, 1), count=0, i, V_min
		string g2tname="g2t_"+name+"_cross 1"
		string g2t_rsname="g2t_rs_"+name
		
		duplicate/o g2t_map g2t
		redimension/n=(-1,0) g2t
		duplicate/o g2t, g2t_temp
		g2t=0
		
		for (i=0;i<N;i+=1)
			g2t_temp[]=g2t_map[p][i]
			wavestats/Q/Z/R=[0, 1000] g2t_temp // ERC.
			//wavestats/Q/Z/R=[0, 900] g2t_temp // 0.2s per frames, 1000 frames.
			//wavestats/Q/Z/R=[0, 1800] g2t_temp // K2, 0.005s per frame, 10s data
			//wavestats/Q/Z/R=[0, 2500] g2t_temp // 0.25min data
			//wavestats/Q/Z/R=[0, 5000] g2t_temp // 0.5min data
			//wavestats/Q/Z/R=[0, 10000] g2t_temp // 1min data
			//wavestats/Q/Z/R=[0, 20000] g2t_temp // 2min data
			//wavestats/Q/Z/R=[0, 39000] g2t_temp // 4min data
			if (V_min<1)
				g2t[]=g2t[p]+g2t_temp[p]
				count=count+1
			endif
		endfor
		
		g2t=g2t/count
		//SetScale/P x 0,3.07,"", g2t
		//SetScale/P x 0,2.1414,"", g2t
		
		make/o/n=(100) g2t_rs, t_rs
		
		//t_rs=alog(log(0.267647)+p/99*(log(0.267647*1000)-log(0.267647))) // 0.2s/frame, 1000 frames, bin8, 200^2 pixels.
		//t_rs=alog(log(0.26547)+p/99*(log(0.26547*1000)-log(0.26547))) // 0.2s/frame, 1000 frames, bin4, 350^2 pixels.
		//t_rs=alog(log(0.005)+p/99*(log(0.005*12000)-log(0.005))) // BNL K2 camera, 400 f/s. Combine two frames together. 0.005 s.
		t_rs=alog(log(0.2363)+p/99*(log(0.2363*1000)-log(0.2363))) // ERC, 0.1s per frame, 350^2 pixels.
		
		//g2t_rs[]=g2t[t_rs[p]/0.267647]
		//g2t_rs[]=g2t[t_rs[p]/0.26547]
             //g2t_rs[]=g2t[t_rs[p]/0.005]
             g2t_rs[]=g2t[t_rs[p]/0.2363]
             	//t_rs=alog(log(0.0584)+p/99*(log(0.0584*2100)-log(0.0584))) // BNL Orius camera, 0.05 s.
             	//g2t_rs[]=g2t[t_rs[p]/0.0584]
             
		duplicate/o g2t, $g2tname
		duplicate/o g2t_rs, $g2t_rsname
		killwaves g2t, g2t_rs, g2t_temp
		print count
end

Function tau_k_map(g2t_map, name, samplepts) //Fitting to extract tau from every pixels. "samplepts" refers to number of points in g2_rs.
		wave g2t_map
		string name
		variable samplepts
		
		variable i//, j
		variable M=dimsize(g2t_map,0), N=dimsize(g2t_map, 1)
		//string fitname, g2FitName, fitXname//, parameterName, g2NormName, , fitNormname, 
		string tau_mapname="tau_map_"+name
		string tauerr_mapname="tauerr_map_"+name
		string beta_mapname="beta_map_"+name
		string betaerror_mapname="betaerr_map_"+name
		string A_mapname="A_map_"+name
		string Aerror_mapname="Aerr_map_"+name
		//string g2tNormmapName="g2_norm_map_"+name
		string rsname="rs_map"+name
		string fitmapName="fit_map_"+name
		//string fitXmapName="fitX_map_"+name
		string t_rsName="t_rs_"+name
		//string fitNormmapName="fit_norm_map_"+name
		string rangeName="fitRange_"+name
		
		make/o/n=(N) t_map, tauerr_map, beta_map, betaerror_map, A_map, Aerror_map, range
		make/o/n=(M) g2t_temp
		make/o/n=(samplepts) g2t_rs, t_rs
		//SetScale/P x 0, 0.948,"", g2t_temp //STEM image series, 0.84 s per frame, total 0.948 s per frame.
 		//SetScale/P x 0, 0.0025,"", g2t_temp
		//SetScale/P x 0, 0.2276,"", g2t_temp
		SetScale/P x 0, 0.27716,"", g2t_temp
		//duplicate/o g2t_map fit_map//, g2tNorm_map, fitNorm_map
		//g2tNorm_map=0
		make/o/n=(200, N) fit_map//, fitX_map
		fit_map=0
		//fitX_map=0
		make/o/n=(samplepts, N) rs_map=0
		//fitNorm_map=0
		
		for (i=0; i<N; i+=1)
			g2t_temp[]=g2t_map[p][i]
			
			resampleg2(g2t_temp, "temp", samplepts)
			wave g2t_rs_temp=$"g2t_rs_temp"
			g2t_rs=g2t_rs_temp
			rs_map[][i]=g2t_rs[p]
			
			findlevel /edge=2/p/q g2t_rs, 1
			//variable V_levelX
			Make/D/N=3/O W_coef
			W_coef[0] = {g2t_rs[0]-1, t_rs[V_levelX],1}
			FuncFit/NTHR=0/q relaxationfunc3 W_coef  g2t_rs[0, round(V_levelX)] /X=t_rs/D
			wave fit_g2t_rs=$"fit_g2t_rs"
			//wave fitX_g2t_rs=$"fitX_g2t_rs"
			wave W_sigma=$"W_sigma"
			A_map[i]=W_coef[0]
			Aerror_map[i]=W_sigma[0]
			t_map[i]=W_coef[1]
			tauerr_map[i]=W_sigma[1]
			beta_map[i]=W_coef[2]
			betaerror_map[i]=W_sigma[2]
			fit_map[][i]=fit_g2t_rs[p]
			//fitX_map[][i]=fitX_g2t_rs[p]
			range[i]=round(V_levelX)
		endfor
			
			//j=0
			//do
				//j=j+1
			//while (g2t_temp[j]>g2t_temp[j+1]*0.993)
			//range[i]=j

			//if (j<10||g2t_temp[j]>1)
				//t_map[i]=nan
				//terror_map[i]=nan
				//beta_map[i]=nan
				//betaerror_map[i]=nan
				//A_map[i]=nan
				//Aerror_map[i]=nan
				//g2tNorm_map[][i]=nan
				//fit_map[][i]=nan
				//fitNorm_map[][i]=nan
			//else
				//fitname=num2istr(i)
				//g2NormName="g2_norm"+fitname
				//g2FitName="fit_"+fitname
				//fitNormName="fit_norm"+fitname
				//parameterName="parameter_"+fitname
				
				//fitg2_noGraph(g2t_temp, range, fitname)
				
				//wave par=$parameterName
				//t_map[i]=par[0]
				//terror_map[i]=par[1]
				//beta_map[i]=par[2]
				//betaerror_map[i]=par[3]
				//A_map[i]=par[4]
				//Aerror_map[i]=par[5]
				
				//wave g2tNorm=$g2NormName
				//g2tNorm_map[][i]=g2tNorm[p]
				
				//wave g2fit=$g2FitName
				//fit_map[][i]=g2fit[p]
				
				//wave fitNorm=$fitNormName
				//fitNorm_map[][i]=fitNorm[p]		
			//endif
			
			//killwaves par, g2tNorm, g2fit, fitNorm, $parameterName, $g2NormName, $g2FitName, $fitNormName
		//endfor
		
		duplicate/o t_map $tau_mapname
		duplicate/o tauerr_map $tauerr_mapname
		duplicate/o beta_map $beta_mapname
		duplicate/o betaerror_map $betaerror_mapname
		duplicate/o A_map $A_mapname
		duplicate/o Aerror_map $Aerror_mapname
		//duplicate/o g2tNorm_map $g2tNormmapName
		duplicate/o fit_map $fitmapName
		//duplicate/o fitX_map $fitXmapName
		duplicate/o t_rs $t_rsName
		//duplicate/o fitNorm_map $fitNormmapName
		duplicate/o range $rangeName
		killwaves g2t_temp, t_map, beta_map, betaerror_map, A_map, Aerror_map, fit_map, tauerr_map, W_coef, W_sigma, g2t_rs_temp, g2t_rs,  t_rs, fit_g2t_rs, t_rs, range//,  fitX_map,  fitX_g2, par, g2tNorm, g2fit, fitNorm, g2tNorm_map, fit_map, fitNorm_map,range
end

Function fitg2_noGraph(g2t, endpoint, name)
		wave g2t
		variable endpoint
		string name
		
		variable A, dA, a1, b1, da1, db1, tau1, dtau1, b_1, db_1, tau, dtau, ptau
		string fitname="fit_"+name
		//string fit1name="fit1_"+name
		string parametername="parameter_"+name
		string normname="g2_norm"+name
		string fitnormname="fit_norm"+name
		
		duplicate/o g2t temp
		redimension/n=(endpoint+1) temp
		
		Make/D/N=3/O W_coef
		W_coef[0] = {1,10,1}
		FuncFit/Q/NTHR=0 relaxationfunc3 W_coef  temp
		A=W_coef[0]
		tau1=W_coef[1]
		b_1=W_coef[2]
		wave W_sigma=$"W_sigma"
		dA=W_sigma[0]
		dtau1=W_sigma[1]
		db_1=W_sigma[2]
		//print A, dA
		//duplicate/o $"fit_temp" $fit1name
		
		duplicate/o temp, y1
		y1[]=ln(-ln((temp[p]-1)/A)/2)
		make/o/n=(endpoint+1) x1
		//x1[1,endpoint]=ln(3.07*p)
		//x1[1,endpoint]=ln(1.574*p)
		//x1[1,endpoint]=ln(0.2276*p)
		x1[1,endpoint]=ln(0.2363*p)
		//x1[1,endpoint]=ln(2.1414*p)
		//Display y1 vs x1
		CurveFit/Q/M=2/W=0 line, y1[0,endpoint]/X=x1[0,endpoint]/D
		a1=W_coef[0]
		b1=W_coef[1]
		da1=W_sigma[0]
		db1=W_sigma[1]
		
		tau=exp(-a1/b1)
		dtau=tau/b1^2*(a1*db1-da1*b1)
		
		duplicate/o temp fit
		//fit[]=1+A*exp(-2*(3.07*p/tau)^b1)
		//fit[]=1+A*exp(-2*(1.574*p/tau)^b1)
		fit[]=1+A*exp(-2*(0.2276*p/tau)^b1)
		//fit[]=1+A*exp(-2*(2.1414*p/tau)^b1)
		
		//display temp
		//ModifyGraph log(bottom)=1
		//appendtograph fit_temp
		//AppendToGraph fit
		//print tau, dtau
		//print b1, db1
		
		duplicate/o g2t temp, temp_norm
		temp_norm[]=(temp[p]-1)/A
		duplicate/o fit fit_norm
		fit_norm[]=(fit[p]-1)/A
		
		//display temp_norm
		//ModifyGraph log(bottom)=1
		//appendtograph fit_norm
		
		//ptau=round(tau/3.07)
		//ptau=round(tau/1.574)
		ptau=round(tau/0.2276)
		//ptau=round(tau/2.1414)
		//print temp_norm[ptau]
		//if (temp_norm[ptau]>exp(-2)*1.01)
			//print "suspicisous fitting"
			//else
				//print "good fitting"
		//endif
		
		duplicate/o fit $fitname
		duplicate/o temp_norm $normname
		duplicate/o fit_norm $fitnormname
		make/o/n=(10) parameter
		parameter[0]=tau
		parameter[1]=dtau
		parameter[2]=b1
		parameter[3]=db1
		parameter[4]=A
		parameter[5]=dA
		parameter[6]=tau1
		parameter[7]=dtau1
		parameter[8]=b_1
		parameter[9]=db_1
		duplicate/o parameter $parametername
		killwaves temp, x1, y1, W_coef, W_sigma, fit, fit_temp, temp_norm, fit_norm, parameter
end

Function displayG2t(g2t_map, name)
		wave g2t_map
		string name
		
		variable i
		variable N=dimsize(g2t_map, 1)

		duplicate/o g2t_map g2t_0
		redimension/n=(-1,0) g2t_0
		g2t_0[]=g2t_map[p][0]
		string g2t_name="g2t_"+name+"_"+num2istr(0)
		duplicate/o g2t_0 $g2t_name
		Display  $g2t_name
		ModifyGraph log(bottom)=1
		
		for (i=1; i<N; i+=1)
			duplicate/o g2t_0 g2t_temp
			g2t_temp[]=g2t_map[p][i]
			g2t_name="g2t_"+name+"_"+num2istr(i)
			duplicate/o g2t_temp $g2t_name
			AppendToGraph $g2t_name
		endfor
		
		killwaves g2t_temp
end
			
Function tau_map(g2t_map, A0, B0, tau0,fitstart,fitend)       // under construction.
                wave g2t_map
                variable A0, B0, tau0,fitstart,fitend
                
                variable pixelsize=dimsize(g2t_map, 1), tp=dimsize(g2t_map,0)
                variable i
                string Aname="A_"+nameofwave(data)         // name is too long.
                string Bname="B_"+nameofwave(data)
                string tauname="tau_"+nameofwave(data)
                
                make/o/n=(pixelsize) tau_temp=0, A_temp=0, B_temp=0
                make/o/n=(tp) g2t_temp=0
                Make/o/D/N=3/O W_coef
                for (i=0; i<pixelsize; i+=1)
                                g2t_temp[]=g2t_map[p][i]
                                W_coef[0] = {A0, B0, tau0}
                                Funcfit/NTHR=0/Q relaxationfunc W_coef, g2t_temp[fitstart,fitend]
                                tau_temp[i]=W_coef[2]
                                A_temp[i]=W_coef[0]
                                B_temp[i]=W_coef[1]
                endfor
                
                tau_temp=tau_temp/3
                duplicate/o tau_temp $tauname
                duplicate/o A_temp $Aname
                duplicate/o B_temp $Bname
                killwaves tau_temp, A_temp, B_temp, g2t_temp
end

Function correlatefunc_k_map_BNLOrius(data, k, bin, xcenter, ycenter, xa, ya, xb, xc, yc, yd, xe, ye)    // <I(t')>^2 from I(t') only. binning pixels of numbers bin^2. 
// BNL Orius data
	       wave data
	       variable k, bin, xcenter, ycenter, xa, ya, xb, xc, yc, yd, xe, ye
	       
	       // Hardware bin4, software bin2. Camera length 380mm.
	       variable r=round(k/0.00352388/2), tp = DimSize(data, 2), neighbour=round((bin-1)/2)
	       variable xmin, xmax, y1, y2, temp, pixelsize=1000
	       variable kbc=(yc-ya)/(xc-xb), kde=(ye-yd)/(xe-xc)
	
               make/o/n=(tp) add_td=0
               make/o/n=(tp,pixelsize) g2t_map=0
               make/o/n=(pixelsize) x_map, y_map, I_ave
	       variable x, i, j
	       variable pcount=0
	       string name="_"+nameofwave(data)+"_"+num2str(k)+"_bin"+num2str(bin)
	       string g2t_mapname="g2t_map"+name
	       string g2tname="g2t"+name
	       //string g2t_rsname="g2t_rs"+name
	       string x_mapname="x_map"+name
	       string y_mapname="y_map"+name
	       string I_avename="I_ave"+name
                                                                                                                                                                         // Find the circular points on diffraction patterns with the reciprocal space vector k.
	       xmin=xcenter-r+1
	       xmax=xcenter+r
	       for(x=xmin; x<xmax; x+=1)
		        temp=round(sqrt(r^2-(xcenter-x)^2))
		        y1=ycenter+temp
		        y2=ycenter-temp
		        // Remove beam stop shadow.
		        if (x>xe && x<xa && y1>ye && y1<ya)
		              pcount=pcount
		       elseif (x>xb && x<xe && y1>yd && y1<ya)
		       	      pcount=pcount
		       elseif (x>xc && x<xb && y1>yd && y1<yc)
		       	      pcount=pcount
		       elseif (x>xc && x<xb && y1>yc && y1<kbc*(x-xc)+yc)
		              pcount=pcount
		        elseif (x>xc && x<xe && y1>kde*(x-xc)+yd && y1<yd)
		              pcount=pcount
		        else
		              make/o/n=(tp) bin_average=0
		              for(i=0-neighbour; i<neighbour+1; i+=1)
		                     for(j=0-neighbour; j<neighbour+1; j+=1)
			                        Imagetransform/beam={(x+i), (y1+j)} getbeam data
			                        wave W_beam = $"W_beam"
			                        bin_average=bin_average+W_beam
			            endfor
			      endfor
			      W_beam=bin_average/bin^2
			      //ave=mean(W_beam)
			      add_td=0
			      for (j=0;j<tp;j+=1)
			            for (i=0;i<(tp-j);i+=1)
			                      add_td[j]=add_td[j]+W_beam[i]*W_beam[i+j]
			            endfor
			      endfor

			      g2t_map[][pcount]=(tp-p)*add_td[p]/(sum(W_beam,0,tp-p-1))^2
			      x_map[pcount]=x
			      y_map[pcount]=y1
			      I_ave[pcount]=mean(W_beam)
			      
                               pcount=pcount+1
                    endif
                        
                        
                    if (x>xe && x<xa && y2>ye && y2<ya)
		              pcount=pcount
		       elseif (x>xb && x<xe && y2>yd && y2<ya)
		       	 pcount=pcount
		       elseif (x>xc && x<xb && y2>yd && y2<yc)
		       	 pcount=pcount
		       elseif (x>xc && x<xb && y2>yc && y2<kbc*(x-xc)+yc)
		              pcount=pcount
		        elseif (x>xc && x<xe && y2>kde*(x-xc)+yd && y2<yd)
		              pcount=pcount
		       else
		             make/o/n=(tp) bin_average=0
		             for(i=0-neighbour; i<neighbour+1; i+=1)
		                     for(j=0-neighbour; j<neighbour+1; j+=1)
			                        Imagetransform/beam={(x+i), (y2+j)} getbeam data
			                        wave W_beam = $"W_beam"
			                        bin_average=bin_average+W_beam
			            endfor
			      endfor
			      W_beam=bin_average/bin^2
			      //ave=mean(W_beam)
			      add_td=0
			      for (j=0;j<tp;j+=1)
			            for (i=0;i<(tp-j);i+=1)
			                      add_td[j]=add_td[j]+W_beam[i]*W_beam[i+j]
			            endfor
			      endfor

			      g2t_map[][pcount]=(tp-p)*add_td[p]/(sum(W_beam,0,tp-p-1))^2
			      x_map[pcount]=x
			      y_map[pcount]=y2
			      I_ave[pcount]=mean(W_beam)
			     
                               pcount=pcount+1
                      endif
                endfor
                
                redimension/n=(-1, pcount) g2t_map
                redimension/n=(pcount) x_map, y_map, I_ave
                SetScale/P x 0, 0.0584,"", g2t_map        //Orius, hardware bin4, 400*410 pixels, 0.05s
                //SetScale/P x 0,0.34,"", g2t_map            time-series spectrum imaging, 0.1s per frame
                //SetScale/P x 0,1.07,"", g2t_map          //CCD serach, 1s per frame
                //SetScale/P x 0,3.07,"", g2t_map
                //SetScale/P x 0,0.171,"", g2t_map
                //SetScale/P x 0,2.075,"", g2t_map
                //SetScale/P x 0,1.574,"", g2t_map        //CCD search, 1.5 s per frame
                //SetScale/P x 0, 0.2276,"", g2t_map       //CCD search, ERC, 0.1 s per frame, Nov. 2013
                //SetScale/P x 0, 0.2363,"", g2t_map       //CCD search, ERC, 0.1 s per frame, March 2013
                //SetScale/P x 0, 2.1414,"", g2t_map       //CCD search, ERC, 2 s per frame
                
                imagetransform sumallrows g2t_map
                duplicate/o W_sumRows g2t
                g2t=g2t/pcount
                SetScale/P x 0,0.0584,"", g2t
                //SetScale/P x 0,0.34,"", g2t
                //SetScale/P x 0,1.07,"", g2t
                //SetScale/P x 0,3.07,"", g2t                 //CCD search, 3 s per frame.
                //SetScale/P x 0,0.171,"", g2t               //CCD search, 0.1 s per frame
                //SetScale/P x 0,2.075,"", g2t                 //CCD search, 2 s per frame.
                //SetScale/P x 0,1.574,"", g2t 			 //CCD search, 1.5 s per frame.
                //SetScale/P x 0,0.2276,"", g2t
                //SetScale/P x 0,0.2363,"", g2t
                //SetScale/P x 0, 2.1414,"", g2t
                
                //make/o/n=(100) g2t_rs, t_rs
                //t_rs=alog(log(0.0584)+p/99*(log(0.0584*2100)-log(0.0584)))
                //g2t_rs[]=g2t[t_rs[p]/0.0584]
                
                duplicate/o g2t_map, $g2t_mapname
                duplicate/o g2t, $g2tname
                //duplicate/o g2t_rs, $g2t_rsname
                duplicate/o x_map, $x_mapname
                duplicate/o y_map, $y_mapname
                duplicate/o I_ave, $I_avename
                
                killwaves g2t_map, g2t, x_map, y_map, I_ave, add_td, bin_average
                print pcount
end

Function correlatefunc_k_map_BNL(data, k, bin, xcenter, ycenter, xa, ya, yb, xc, yc, xd, xe, ye)    // <I(t')>^2 from I(t'), I(t')I(t'+t), and I all. Binning pixels of numbers bin^2. 
// BNL K2 data
	       wave data
	       variable k, bin, xcenter, ycenter, xa, ya, yb, xc, yc, xd, xe, ye
	      
	       variable r=k/0.008052248 // Bin 8. Camera length 380mm.
	       //variable r=k/0.01006531 // Bin 10. Camera length 380mm.
	       variable tp = DimSize(data, 2), neighbour=round((bin-1)/2)
	       variable samplesize=round(2*pi*r)
	       variable kbc=(yc-yb)/(xc-xa), kde=(ye-yc)/(xe-xd)
	
               make/o/n=(tp) add_td=0
               make/o/n=(tp,samplesize) g2t_tp_map=0, g2t_tpt_map=0, g2t_ta_map=0
               make/o/n=(samplesize) x_map, y_map, I_ave
               make/o/n=(samplesize) theta, xp, yp
               theta[]=p*2*pi/samplesize
               xp=round(xcenter+r*cos(theta))
               yp=round(ycenter+r*sin(theta))
               
	       variable i, j, ave
	       variable pcount=0
	       string name="_"+nameofwave(data)+"_"+num2str(k)+"_bin"+num2str(bin)
	       string g2t_tp_mapname="g2tmap_tp"+name
	       string g2t_tpt_mapname="g2tmap_tpt"+name
	       string g2t_ta_mapname="g2tmap_ta"+name
	       string g2t_tp_name="g2t_tp"+name
	       string g2t_tpt_name="g2t_tpt"+name
	       string g2t_ta_name="g2t_ta"+name
	       string x_mapname="x_map"+name
	       string y_mapname="y_map"+name
	       string I_avename="I_ave"+name
	                                                                                                                                                               // Find the circular points on diffraction patterns with the reciprocal space vector k.

	       for(k=0; k<samplesize; k+=1)
		        // Remove beam stop shadow
		        if (xp[k]>xc && xp[k]<xd && yp[k]>ya && yp[k]<yc)
		              pcount=pcount
		       elseif (xp[k]>xa && xp[k]<xc && yp[k]>ya && yp[k]<yb)
		       	      pcount=pcount
		       elseif (xp[k]>xa && xp[k]<xc && yp[k]>yb && yp[k]<kbc*(xp[k]-xa)+yb)
		       	      pcount=pcount
		       elseif (xp[k]>xd && xp[k]<xe && yp[k]>ya && yp[k]<ye)
		              pcount=pcount
		        elseif (xp[k]>xd && xp[k]<xe && yp[k]>ye && yp[k]<kde*(xp[k]-xd)+yc)
		              pcount=pcount
		       elseif(xp[k]>xa && xp[k]<xe && yp[k]>ya && yp[k]<yb)
		        		pcount=pcount
		        else
		              make/o/n=(tp) bin_average=0
		              for(i=0-neighbour; i<neighbour+1; i+=1)
		                     for(j=0-neighbour; j<neighbour+1; j+=1)
			                        Imagetransform/beam={(xp[k]+i), (yp[k]+j)} getbeam data
			                        wave W_beam = $"W_beam"
			                        bin_average=bin_average+W_beam
			            endfor
			      endfor
			      W_beam=bin_average/bin^2
			      ave=mean(W_beam)
			      add_td=0
			      for (j=0;j<tp;j+=1)
			            for (i=0;i<(tp-j);i+=1)
			                      add_td[j]=add_td[j]+W_beam[i]*W_beam[i+j]
			            endfor
			      endfor

			      g2t_tp_map[][pcount]=(tp-p)*add_td[p]/(sum(W_beam,0,tp-p-1))^2
			      g2t_tpt_map[][pcount]=(tp-p)*add_td[p]/sum(W_beam,0,tp-p-1)/sum(W_beam,p,tp-1)
			      g2t_ta_map[][pcount]=add_td[p]/(tp-p)/ave^2
			      x_map[pcount]=xp[k]
			      y_map[pcount]=yp[k]
			      I_ave[pcount]=ave
			      
                               pcount=pcount+1
                        endif
                 endfor
                
                redimension/n=(-1, pcount) g2t_tp_map, g2t_tpt_map, g2t_ta_map
                redimension/n=(pcount) x_map, y_map, I_ave
                SetScale/P x 0, 0.005,"", g2t_tp_map, g2t_tpt_map, g2t_ta_map        //K2, 400 frames/s, adding 2 frames together
                
                imagetransform sumallrows g2t_tp_map
                duplicate/o W_sumRows g2t
                g2t=g2t/pcount
                SetScale/P x 0,0.005,"", g2t
                
                duplicate/o g2t_tp_map, $g2t_tp_mapname
                duplicate/o g2t, $g2t_tp_name
                
                imagetransform sumallrows g2t_tpt_map
                duplicate/o W_sumRows g2t
                g2t=g2t/pcount
                SetScale/P x 0,0.005,"", g2t
                
                duplicate/o g2t_tpt_map, $g2t_tpt_mapname
                duplicate/o g2t, $g2t_tpt_name
                
                 imagetransform sumallrows g2t_ta_map
                duplicate/o W_sumRows g2t
                g2t=g2t/pcount
                SetScale/P x 0,0.005,"", g2t
                
                duplicate/o g2t_ta_map, $g2t_ta_mapname
                duplicate/o g2t, $g2t_ta_name
                
                duplicate/o x_map, $x_mapname
                duplicate/o y_map, $y_mapname
                duplicate/o I_ave, $I_avename
                
                killwaves g2t_tp_map, g2t_tpt_map, g2t_ta_map, g2t, x_map, y_map, I_ave, add_td, bin_average, W_beam, W_sumrows, xp, yp, theta
                print pcount
end


Function correlatefunc_k_map_BNL1(data, k, bin, xcenter, ycenter, xa, ya, yb, xc, yc, xd, xe, ye)    // <I(t')>^2 from I(t') only. binning pixels of numbers bin^2. 
// BNL K2 data
	       wave data
	       variable k, bin, xcenter, ycenter, xa, ya, yb, xc, yc, xd, xe, ye
	       
	       // Bin 8. Camera length 380mm.
	       variable r=round(k/0.007700351), tp = DimSize(data, 2), neighbour=round((bin-1)/2)
	       variable xmin, xmax, y1, y2, temp, pixelsize=1000
	       variable kbc=(yc-yb)/(xc-xa), kde=(ye-yc)/(xe-xd)
	
               make/o/n=(tp) add_td=0
               make/o/n=(tp,pixelsize) g2t_map=0
               make/o/n=(pixelsize) x_map, y_map, I_ave
	       variable x, i, j
	       variable pcount=0
	       string name="_"+nameofwave(data)+"_"+num2str(k)+"_bin"+num2str(bin)
	       string g2t_mapname="g2t_map"+name
	       string g2tname="g2t"+name
	       string x_mapname="x_map"+name
	       string y_mapname="y_map"+name
	       string I_avename="I_ave"+name
	                                                                                                                                                               // Find the circular points on diffraction patterns with the reciprocal space vector k.
	       xmin=xcenter-r+1
	       xmax=xcenter+r
	       for(x=xmin; x<xmax; x+=1)
		        temp=round(sqrt(r^2-(xcenter-x)^2))
		        y1=ycenter+temp
		        y2=ycenter-temp
		        // Remove beam stop shadow.
		        if (x>xc && x<xd && y1>ya && y1<yc)
		              pcount=pcount
		       elseif (x>xa && x<xc && y1>ya && y1<yb)
		       	      pcount=pcount
		       elseif (x>xa && x<xc && y1>yb && y1<kbc*(x-xa)+yb)
		       	      pcount=pcount
		       elseif (x>xd && x<xe && y1>ya && y1<ye)
		              pcount=pcount
		        elseif (x>xd && x<xe && y1>ye && y1<kde*(x-xd)+yc)
		              pcount=pcount
		        else
		              make/o/n=(tp) bin_average=0
		              for(i=0-neighbour; i<neighbour+1; i+=1)
		                     for(j=0-neighbour; j<neighbour+1; j+=1)
			                        Imagetransform/beam={(x+i), (y1+j)} getbeam data
			                        wave W_beam = $"W_beam"
			                        bin_average=bin_average+W_beam
			            endfor
			      endfor
			      W_beam=bin_average/bin^2
			      //ave=mean(W_beam)
			      add_td=0
			      for (j=0;j<tp;j+=1)
			            for (i=0;i<(tp-j);i+=1)
			                      add_td[j]=add_td[j]+W_beam[i]*W_beam[i+j]
			            endfor
			      endfor

			      g2t_map[][pcount]=(tp-p)*add_td[p]/(sum(W_beam,0,tp-p-1))^2
			      x_map[pcount]=x
			      y_map[pcount]=y1
			      I_ave[pcount]=mean(W_beam)
			      
                               pcount=pcount+1
                        endif
                        
        	 	 	if (x>xc && x<xd && y2>ya && y2<yc)
		              pcount=pcount
		       elseif (x>xa && x<xc && y2>ya && y2<yb)
		       	      pcount=pcount
		       elseif (x>xa && x<xc && y2>yb && y2<kbc*(x-xa)+yb)
		       	      pcount=pcount
		       elseif (x>xd && x<xe && y2>ya && y2<ye)
		              pcount=pcount
		        elseif (x>xd && x<xe && y2>ye && y2<kde*(x-xd)+yc)
		              pcount=pcount
		       else
		             make/o/n=(tp) bin_average=0
		             for(i=0-neighbour; i<neighbour+1; i+=1)
		                     for(j=0-neighbour; j<neighbour+1; j+=1)
			                        Imagetransform/beam={(x+i), (y2+j)} getbeam data
			                        wave W_beam = $"W_beam"
			                        bin_average=bin_average+W_beam
			            endfor
			      endfor
			      W_beam=bin_average/bin^2
			      //ave=mean(W_beam)
			      add_td=0
			      for (j=0;j<tp;j+=1)
			            for (i=0;i<(tp-j);i+=1)
			                      add_td[j]=add_td[j]+W_beam[i]*W_beam[i+j]
			            endfor
			      endfor

			      g2t_map[][pcount]=(tp-p)*add_td[p]/(sum(W_beam,0,tp-p-1))^2
			      x_map[pcount]=x
			      y_map[pcount]=y2
			      I_ave[pcount]=mean(W_beam)
			     
                               pcount=pcount+1
                      endif
                endfor
                
                redimension/n=(-1, pcount) g2t_map
                redimension/n=(pcount) x_map, y_map, I_ave
                SetScale/P x 0, 0.005,"", g2t_map        //K2, 400 frames/s, adding 2 frames together
                //SetScale/P x 0,0.34,"", g2t_map            time-series spectrum imaging, 0.1s per frame
                //SetScale/P x 0,1.07,"", g2t_map          //CCD serach, 1s per frame
                //SetScale/P x 0,3.07,"", g2t_map
                //SetScale/P x 0,0.171,"", g2t_map
                //SetScale/P x 0,2.075,"", g2t_map
                //SetScale/P x 0,1.574,"", g2t_map        //CCD search, 1.5 s per frame
                //SetScale/P x 0, 0.2276,"", g2t_map       //CCD search, ERC, 0.1 s per frame, Nov. 2013
                //SetScale/P x 0, 0.2363,"", g2t_map       //CCD search, ERC, 0.1 s per frame, March 2013
                //SetScale/P x 0, 2.1414,"", g2t_map       //CCD search, ERC, 2 s per frame
                
                imagetransform sumallrows g2t_map
                duplicate/o W_sumRows g2t
                g2t=g2t/pcount
                SetScale/P x 0,0.005,"", g2t
                //SetScale/P x 0,0.34,"", g2t
                //SetScale/P x 0,1.07,"", g2t
                //SetScale/P x 0,3.07,"", g2t                 //CCD search, 3 s per frame.
                //SetScale/P x 0,0.171,"", g2t               //CCD search, 0.1 s per frame
                //SetScale/P x 0,2.075,"", g2t                 //CCD search, 2 s per frame.
                //SetScale/P x 0,1.574,"", g2t 			 //CCD search, 1.5 s per frame.
                //SetScale/P x 0,0.2276,"", g2t
                //SetScale/P x 0,0.2363,"", g2t
                //SetScale/P x 0, 2.1414,"", g2t
                
                duplicate/o g2t_map, $g2t_mapname
                duplicate/o g2t, $g2tname
                duplicate/o x_map, $x_mapname
                duplicate/o y_map, $y_mapname
                duplicate/o I_ave, $I_avename
                
                killwaves g2t_map, g2t, x_map, y_map, I_ave, add_td, bin_average
                print pcount
end

Function correlatefunc_k_map_ERC(data,k,bin,xcenter,ycenter,xa,ya,xb,yb,xc,xd,yd,xe,xf)    // <I(t')>^2 from I(t') only. binning pixels of numbers bin^2.
	       wave data
	       variable k, bin, xcenter, ycenter,xa,ya,xb,yb,xc,xd,yd,xe,xf
	       
	       variable r=round(k/0.006967), tp = DimSize(data, 2), neighbour=round((bin-1)/2)
	       variable xmin, xmax, y1, y2, temp, pixelsize=1000
	       variable kab=(yb-ya)/(xb-xa), kcd=(yd-yb)/(xd-xc), kef=(ya-yd)/(xf-xe)
	
               make/o/n=(tp) add_td=0
               make/o/n=(tp,pixelsize) g2t_map=0
               make/o/n=(pixelsize) x_map, y_map, I_ave
	       variable x, i, j
	       variable pcount=0
	       string name="_"+nameofwave(data)+"_"+num2str(k)+"_bin"+num2str(bin)
	       string g2t_mapname="g2t_map"+name
	       string g2tname="g2t"+name
	       string x_mapname="x_map"+name
	       string y_mapname="y_map"+name
	       string I_avename="I_ave"+name
	                                                                                                                                                                   // Find the circular points on diffraction patterns with the reciprocal space vector k.
	       xmin=xcenter-r+1
	       xmax=xcenter+r
	       for(x=xmin; x<xmax; x+=1)
		        temp=round(sqrt(r^2-(xcenter-x)^2))
		        y1=ycenter+temp
		        y2=ycenter-temp
		        if (x>xd && x<xe && y1>yd && y1<ya)
		              pcount=pcount
		       elseif (x>xc && x<xd && y1>kcd*(x-xc)+yb)
		       	      pcount=pcount
		       elseif (x>xb && x<xc && y1>yb)
		       	      pcount=pcount
		       elseif (x>xa && x<xb && y1>kab*(x-xa)+ya)
		              pcount=pcount
		        elseif (x>xe && x<xf && y1>kef*(x-xe)+yd)
		              pcount=pcount
		        else
		              make/o/n=(tp) bin_average=0
		              for(i=0-neighbour; i<neighbour+1; i+=1)
		                     for(j=0-neighbour; j<neighbour+1; j+=1)
			                        Imagetransform/beam={(x+i), (y1+j)} getbeam data
			                        wave W_beam = $"W_beam"
			                        bin_average=bin_average+W_beam
			            endfor
			      endfor
			      W_beam=bin_average/bin^2
			      //ave=mean(W_beam)
			      add_td=0
			      for (j=0;j<tp;j+=1)
			            for (i=0;i<(tp-j);i+=1)
			                      add_td[j]=add_td[j]+W_beam[i]*W_beam[i+j]
			            endfor
			      endfor

			      g2t_map[][pcount]=(tp-p)*add_td[p]/(sum(W_beam,0,tp-p-1))^2
			      x_map[pcount]=x
			      y_map[pcount]=y1
			      I_ave[pcount]=mean(W_beam)
			      
                               pcount=pcount+1
                        endif
                        
        	 	if (x>xd && x<xe && y1>yd && y1<ya)   // WRONG?
		              pcount=pcount
		       elseif (x>xc && x<xd && y1>kcd*(x-xc)+yb)
		       	      pcount=pcount
		       elseif (x>xb && x<xc && y1>yb)
		       	      pcount=pcount
		       elseif (x>xa && x<xb && y1>kab*(x-xa)+ya)
		              pcount=pcount
		        elseif (x>xe && x<xf && y1>kef*(x-xe)+yd)
		              pcount=pcount
		       else
		             make/o/n=(tp) bin_average=0
		             for(i=0-neighbour; i<neighbour+1; i+=1)
		                     for(j=0-neighbour; j<neighbour+1; j+=1)
			                        Imagetransform/beam={(x+i), (y2+j)} getbeam data
			                        wave W_beam = $"W_beam"
			                        bin_average=bin_average+W_beam
			            endfor
			      endfor
			      W_beam=bin_average/bin^2
			      //ave=mean(W_beam)
			      add_td=0
			      for (j=0;j<tp;j+=1)
			            for (i=0;i<(tp-j);i+=1)
			                      add_td[j]=add_td[j]+W_beam[i]*W_beam[i+j]
			            endfor
			      endfor

			      g2t_map[][pcount]=(tp-p)*add_td[p]/(sum(W_beam,0,tp-p-1))^2
			      x_map[pcount]=x
			      y_map[pcount]=y2
			      I_ave[pcount]=mean(W_beam)
			     
                               pcount=pcount+1
                      endif
                endfor
                
                redimension/n=(-1, pcount) g2t_map
                redimension/n=(pcount) x_map, y_map, I_ave
                //SetScale/P x 0,0.34,"", g2t_map            time-series spectrum imaging, 0.1s per frame
                //SetScale/P x 0,1.07,"", g2t_map          //CCD serach, 1s per frame
                //SetScale/P x 0,3.07,"", g2t_map
                //SetScale/P x 0,0.171,"", g2t_map
                //SetScale/P x 0,2.075,"", g2t_map
                //SetScale/P x 0,1.574,"", g2t_map        //CCD search, 1.5 s per frame
                //SetScale/P x 0, 0.2276,"", g2t_map       //CCD search, ERC, 0.1 s per frame, Nov. 2013
                SetScale/P x 0, 0.2363,"", g2t_map       //CCD search, ERC, 0.1 s per frame, March 2013
                //SetScale/P x 0, 2.1414,"", g2t_map       //CCD search, ERC, 2 s per frame
                
                imagetransform sumallrows g2t_map
                duplicate/o W_sumRows g2t
                g2t=g2t/pcount
                //SetScale/P x 0,0.34,"", g2t
                //SetScale/P x 0,1.07,"", g2t
                //SetScale/P x 0,3.07,"", g2t                 //CCD search, 3 s per frame.
                //SetScale/P x 0,0.171,"", g2t               //CCD search, 0.1 s per frame
                //SetScale/P x 0,2.075,"", g2t                 //CCD search, 2 s per frame.
                //SetScale/P x 0,1.574,"", g2t 			 //CCD search, 1.5 s per frame.
                //SetScale/P x 0,0.2276,"", g2t
                SetScale/P x 0,0.2363,"", g2t
                //SetScale/P x 0, 2.1414,"", g2t
                
                duplicate/o g2t_map, $g2t_mapname
                duplicate/o g2t, $g2tname
                duplicate/o x_map, $x_mapname
                duplicate/o y_map, $y_mapname
                duplicate/o I_ave, $I_avename
                
                killwaves g2t_map, g2t, x_map, y_map, I_ave, add_td, bin_average
                print pcount
end

Function correlatefunc_k_map(data,k,bin,xcenter,ycenter,xblockstart,xblockend,yblockstart,yblockend)    // Four ways of g2 calculation. binning pixels of numbers bin^2.
	       wave data
	       variable k, bin, xcenter, ycenter,xblockstart,xblockend,yblockstart,yblockend	       
	       
	       variable r=k/0.00122961	// 4.4m, bin8.
	       //variable r=k/(0.00545597*2) // 510 mm, Bin 8.
	       //variable r=k/0.00545597 // 510 mm, Bin 4.
	       //variable r=k/(0.000487225*2) // Bin 8. Camera length 5.5 m.
	       variable tp = DimSize(data, 2), neighbour=round((bin-1)/2)
	       variable samplesize=round(2*pi*r)
	
              make/o/n=(tp) add_td=0//, add_4n=0, add_4di=0, add_4dinp=0
              make/o/n=(tp,samplesize) g2t_tpt_map=0//, g2t_tp_map=0, g2t_ta_map=0//, g2t4_map
              make/o/n=(samplesize) x_map, y_map, I_ave
              make/o/n=(samplesize) theta, xp, yp
              theta[]=p*2*pi/samplesize
              xp=round(xcenter+r*cos(theta))
              yp=round(ycenter+r*sin(theta))
               
	       variable i, j, ave, ave_i, ave_inp
	       variable pcount=0
	       string name="_"+nameofwave(data)+"_"+num2str(k)+"_bin"+num2str(bin)
	       //string g2t_tp_mapname="g2tmap_tp"+name
	       string g2t_tpt_mapname="g2tmap_tpt"+name
	       //string g2t_ta_mapname="g2tmap_ta"+name
	       //string g2t4_mapname="g2tmap4"+name
	       //string g2t_tp_name="g2t_tp"+name
	       string g2t_tpt_name="g2t_tpt"+name
	       //string g2t_ta_name="g2t_ta"+name
	       //string g2t4_name="g2t4"+name
	       string x_mapname="x_map"+name
	       string y_mapname="y_map"+name
	       string I_avename="I_ave"+name
	                                                                                                                                                               // Find the circular points on diffraction patterns with the reciprocal space vector k.

	       for(k=0; k<samplesize; k+=1)
		        // Remove beam stop shadow
		        if (xp[k]>xblockstart && xp[k]<xblockend && yp[k]>yblockstart && yp[k]<yblockend)
		              pcount=pcount
		        else
		              make/o/n=(tp) bin_average=0
		              for(i=0-neighbour; i<neighbour+1; i+=1)
		                     for(j=0-neighbour; j<neighbour+1; j+=1)
			                        Imagetransform/beam={(xp[k]+i), (yp[k]+j)} getbeam data
			                        wave W_beam = $"W_beam"
			                        bin_average=bin_average+W_beam
			            endfor
			      endfor
			      W_beam=bin_average/bin^2
			      ave=mean(W_beam)
			      add_td=0
			      //add_4n=0
			      //add_4di=0
			      //add_4dinp=0
			      for (j=0;j<tp;j+=1)
			      		//ave_i=mean(W_beam, 0, tp-j-1)
			      		//ave_inp=mean(W_beam, j, tp-1)
			            for (i=0;i<(tp-j);i+=1)
			                      add_td[j]=add_td[j]+W_beam[i]*W_beam[i+j]
			                      //add_4n[j]=add_4n[j]+(W_beam[i]-ave_i)*(W_beam[i+j]-ave_inp)
			                      //add_4di[j]=add_4di[j]+(W_beam[i]-ave_i)^2
			                      //add_4dinp[j]=add_4dinp[j]+(W_beam[i+j]-ave_inp)^2
			            endfor
			      endfor

			      //g2t_tp_map[][pcount]=(tp-p)*add_td[p]/(sum(W_beam,0,tp-p-1))^2
			      g2t_tpt_map[][pcount]=(tp-p)*add_td[p]/sum(W_beam,0,tp-p-1)/sum(W_beam,p,tp-1)
			      //g2t_ta_map[][pcount]=add_td[p]/(tp-p)/ave^2
			      //g2t4_map[][pcount]=add_4n[p]/sqrt(add_4di[p])/sqrt(add_4dinp[p])
			      x_map[pcount]=xp[k]
			      y_map[pcount]=yp[k]
			      I_ave[pcount]=ave
			      
                         pcount=pcount+1
                    endif
             endfor
                
             redimension/n=(-1, pcount) g2t_tpt_map//, g2t_tp_map, g2t_ta_map//, g2t4_map
             redimension/n=(pcount) x_map, y_map, I_ave
             
             SetScale/P x 0, 1.528,"", g2t_tpt_map//, g2t_tp_map, g2t_ta_map	// cluster size map, 1.528 nm interpixel distance. 
             //SetScale/P x 0, 0.248587083,"", g2t_tp_map, g2t_tpt_map, g2t_ta_map//, g2t4_map        //200 ^2 pixels, bin8, 0.18 s per frame, search.
             //SetScale/P x 0, 0.267647,"", g2t_tp_map, g2t_tpt_map, g2t_ta_map        //200 ^2 pixels, bin8, 0.2 s per frame, search.
             //SetScale/P x 0, 0.26547,"", g2t_tp_map, g2t_tpt_map, g2t_ta_map        //350 ^2 pixels, bin4, 0.2 s per frame, search.
             //SetScale/P x 0, 2,"", g2t_tp_map, g2t_tpt_map, g2t_ta_map        //cluster size map, 2 nm interpixel distance.
             //SetScale/P x 0, 0.005,"", g2t_tp_map, g2t_tpt_map, g2t_ta_map        //K2, 400 frames/s, adding 2 frames together
                
             //imagetransform sumallrows g2t_tp_map
             //duplicate/o W_sumRows g2t
             //g2t=g2t/pcount
             //SetScale/P x 0,1.528,"", g2t
             //SetScale/P x 0,0.248587083,"", g2t
             //SetScale/P x 0,0.267647,"", g2t
             //SetScale/P x 0,0.26547,"", g2t
             //SetScale/P x 0,2,"", g2t
             //SetScale/P x 0,0.005,"", g2t
                
             //duplicate/o g2t_tp_map, $g2t_tp_mapname
             //duplicate/o g2t, $g2t_tp_name
                
             imagetransform sumallrows g2t_tpt_map
             duplicate/o W_sumRows g2t
             g2t=g2t/pcount
             SetScale/P x 0,1.528,"", g2t
             //SetScale/P x 0,0.248587083,"", g2t
             //SetScale/P x 0,0.267647,"", g2t
             //SetScale/P x 0,0.26547,"", g2t
             //SetScale/P x 0,2,"", g2t
             //SetScale/P x 0,0.005,"", g2t
                
             duplicate/o g2t_tpt_map, $g2t_tpt_mapname
             duplicate/o g2t, $g2t_tpt_name
                
             //imagetransform sumallrows g2t_ta_map
             //duplicate/o W_sumRows g2t
             //g2t=g2t/pcount
             //SetScale/P x 0,1.528,"", g2t
             //SetScale/P x 0,0.248587083,"", g2t
             //SetScale/P x 0,0.267647,"", g2t
             //SetScale/P x 0,0.26547,"", g2t
             //SetScale/P x 0,2,"", g2t
             //SetScale/P x 0,0.005,"", g2t
                
             //duplicate/o g2t_ta_map, $g2t_ta_mapname
             //duplicate/o g2t, $g2t_ta_name
             
            	//imagetransform sumallrows g2t4_map
             //duplicate/o W_sumRows g2t
             //g2t=g2t/pcount
             //SetScale/P x 0,0.248587083,"", g2t
             //SetScale/P x 0,0.267647,"", g2t
             //SetScale/P x 0,0.26547,"", g2t
             //SetScale/P x 0,2,"", g2t
             //SetScale/P x 0,0.005,"", g2t
                
             //duplicate/o g2t4_map, $g2t4_mapname
             //duplicate/o g2t, $g2t4_name
                
             duplicate/o x_map, $x_mapname
             duplicate/o y_map, $y_mapname
             duplicate/o I_ave, $I_avename
                
             killwaves g2t_tpt_map//, g2t_tp_map, g2t_ta_map//, g2t4_map, 
             killwaves g2t, x_map, y_map, I_ave, add_td//, add_4n, add_4di, add_4dinp
             killwaves bin_average, W_beam, W_sumrows, xp, yp, theta
             print pcount
end

Function correlatefunc_half_t(data,k,bin,xcenter,ycenter,xblockstart,xblockend,yblockstart,yblockend)    // Three ways of <I(t')>^2. binning pixels of numbers bin^2.
	       wave data
	       variable k, bin, xcenter, ycenter,xblockstart,xblockend,yblockstart,yblockend	       
	       
	       variable r=k/(0.00545597*2) // 510 mm, Bin 8..
	       variable tp = DimSize(data, 2)/2, neighbour=round((bin-1)/2)
	       variable samplesize=round(2*pi*r)
	
              make/o/n=(tp) add_td=0
              make/o/n=(tp,samplesize) g2t_tp_map=0, g2t_tpt_map=0, g2t_ta_map=0
              make/o/n=(samplesize) x_map, y_map, I_ave
              make/o/n=(samplesize) theta, xp, yp
              theta[]=p*2*pi/samplesize
              xp=round(xcenter+r*cos(theta))
              yp=round(ycenter+r*sin(theta))
               
	       variable i, j, ave
	       variable pcount=0
	       string name="_"+nameofwave(data)+"_"+num2str(k)+"_bin"+num2str(bin)
	       string g2t_tp_mapname="g2tmap_ht_tp"+name
	       string g2t_tpt_mapname="g2tmap_ht_tpt"+name
	       string g2t_ta_mapname="g2tmap_ht_ta"+name
	       string g2t_tp_name="g2t_ht_tp"+name
	       string g2t_tpt_name="g2t_ht_tpt"+name
	       string g2t_ta_name="g2t_ht_ta"+name
	       string x_mapname="x_map"+name
	       string y_mapname="y_map"+name
	       string I_avename="I_ave"+name
	                                                                                                                                                               // Find the circular points on diffraction patterns with the reciprocal space vector k.

	       for(k=0; k<samplesize; k+=1)
		        // Remove beam stop shadow
		        if (xp[k]>xblockstart && xp[k]<xblockend && yp[k]>yblockstart && yp[k]<yblockend)
		              pcount=pcount
		        else
		              make/o/n=(tp) bin_average=0
		              for(i=0-neighbour; i<neighbour+1; i+=1)
		                     for(j=0-neighbour; j<neighbour+1; j+=1)
			                        Imagetransform/beam={(xp[k]+i), (yp[k]+j)} getbeam data
			                        wave W_beam = $"W_beam"
			                        bin_average=bin_average+W_beam
			            endfor
			      endfor
			      W_beam=bin_average/bin^2
			      ave=mean(W_beam)
			      add_td=0
			      for (j=0;j<tp;j+=1)
			            for (i=0;i<tp;i+=1)
			                      add_td[j]=add_td[j]+W_beam[i]*W_beam[i+j]
			            endfor
			      endfor

			      g2t_tp_map[][pcount]=tp*add_td[p]/(sum(W_beam,0,tp-1))^2
			      g2t_tpt_map[][pcount]=tp*add_td[p]/sum(W_beam,0,tp-1)/sum(W_beam,p,tp+p-1)
			      g2t_ta_map[][pcount]=add_td[p]/tp/ave^2
			      x_map[pcount]=xp[k]
			      y_map[pcount]=yp[k]
			      I_ave[pcount]=ave
			      
                         pcount=pcount+1
                    endif
             endfor
                
             redimension/n=(-1, pcount) g2t_tp_map, g2t_tpt_map, g2t_ta_map
             redimension/n=(pcount) x_map, y_map, I_ave
             
             SetScale/P x 0, 0.267647,"", g2t_tp_map, g2t_tpt_map, g2t_ta_map        //200 ^2 pixels, bin8, 0.2 s per frame, search.
                
             imagetransform sumallrows g2t_tp_map
             duplicate/o W_sumRows g2t
             g2t=g2t/pcount
             SetScale/P x 0,0.267647,"", g2t
                
             duplicate/o g2t_tp_map, $g2t_tp_mapname
             duplicate/o g2t, $g2t_tp_name
                
             imagetransform sumallrows g2t_tpt_map
             duplicate/o W_sumRows g2t
             g2t=g2t/pcount
             SetScale/P x 0,0.267647,"", g2t
                
             duplicate/o g2t_tpt_map, $g2t_tpt_mapname
             duplicate/o g2t, $g2t_tpt_name
                
             imagetransform sumallrows g2t_ta_map
             duplicate/o W_sumRows g2t
             g2t=g2t/pcount
             SetScale/P x 0,0.267647,"", g2t
                
             duplicate/o g2t_ta_map, $g2t_ta_mapname
             duplicate/o g2t, $g2t_ta_name
                
             duplicate/o x_map, $x_mapname
             duplicate/o y_map, $y_mapname
             duplicate/o I_ave, $I_avename
                
             killwaves g2t_tp_map, g2t_tpt_map, g2t_ta_map, g2t, x_map, y_map, I_ave, add_td, bin_average, W_beam, W_sumrows, xp, yp, theta
             print pcount
end

Function correlatefunc_k_map1(data,k,bin,xcenter,ycenter,xblockstart,xblockend,yblockstart,yblockend)    // <I(t')>^2 from I(t') only. binning pixels of numbers bin^2.
	       wave data
	       variable k, bin, xcenter, ycenter,xblockstart,xblockend,yblockstart,yblockend
	       
	       //variable r=round(k/(0.001318073*2)), tp = DimSize(data, 2), neighbour=round((bin-1)/2)  // 2.15 m, bin 8.
	       variable r=round(k/0.01115), tp = DimSize(data, 2), neighbour=round((bin-1)/2) // 510mm
	       variable xmin, xmax, y1, y2, temp, pixelsize=1000
	
               make/o/n=(tp) add_td=0
               make/o/n=(tp,pixelsize) g2t_map=0
               make/o/n=(pixelsize) x_map, y_map, I_ave
	       variable x, i, j
	       variable pcount=0
	       string name="_"+nameofwave(data)+"_"+num2str(k)+"_bin"+num2str(bin)
	       string g2t_mapname="g2t_map"+name
	       string g2tname="g2t"+name
	       string x_mapname="x_map"+name
	       string y_mapname="y_map"+name
	       string I_avename="I_ave"+name
	                                                                                                                                                                   // Find the circular points on diffraction patterns with the reciprocal space vector k.
	       xmin=xcenter-r+1
	       xmax=xcenter+r
	       for(x=xmin; x<xmax; x+=1)
		        temp=round(sqrt(r^2-(xcenter-x)^2))
		        y1=ycenter+temp
		        y2=ycenter-temp
		        if (x>xblockstart && x<xblockend && y1>yblockstart && y1<yblockend)
		              pcount=pcount
		       else
		              make/o/n=(tp) bin_average=0
		              for(i=0-neighbour; i<neighbour+1; i+=1)
		                     for(j=0-neighbour; j<neighbour+1; j+=1)
			                        Imagetransform/beam={(x+i), (y1+j)} getbeam data
			                        wave W_beam = $"W_beam"
			                        bin_average=bin_average+W_beam
			            endfor
			      endfor
			      W_beam=bin_average/bin^2
			      //ave=mean(W_beam)
			      add_td=0
			      for (j=0;j<tp;j+=1)
			            for (i=0;i<(tp-j);i+=1)
			                      add_td[j]=add_td[j]+W_beam[i]*W_beam[i+j]
			            endfor
			      endfor

			      g2t_map[][pcount]=(tp-p)*add_td[p]/(sum(W_beam,0,tp-p-1))^2
			      x_map[pcount]=x
			      y_map[pcount]=y1
			      I_ave[pcount]=mean(W_beam)
			      
                               pcount=pcount+1
                        endif
                        
                        if (x>xblockstart && x<xblockend && y2>yblockstart && y2<yblockend)
		             pcount=pcount
		       else
		             make/o/n=(tp) bin_average=0
		             for(i=0-neighbour; i<neighbour+1; i+=1)
		                     for(j=0-neighbour; j<neighbour+1; j+=1)
			                        Imagetransform/beam={(x+i), (y2+j)} getbeam data
			                        wave W_beam = $"W_beam"
			                        bin_average=bin_average+W_beam
			            endfor
			      endfor
			      W_beam=bin_average/bin^2
			      //ave=mean(W_beam)
			      add_td=0
			      for (j=0;j<tp;j+=1)
			            for (i=0;i<(tp-j);i+=1)
			                      add_td[j]=add_td[j]+W_beam[i]*W_beam[i+j]
			            endfor
			      endfor

			      g2t_map[][pcount]=(tp-p)*add_td[p]/(sum(W_beam,0,tp-p-1))^2
			      x_map[pcount]=x
			      y_map[pcount]=y2
			      I_ave[pcount]=mean(W_beam)
			     
                               pcount=pcount+1
                      endif
                endfor
                
                redimension/n=(-1, pcount) g2t_map
                redimension/n=(pcount) x_map, y_map, I_ave
                //SetScale/P x 0,0.2192825,"", g2t_map    //CCD search, 0.15 s.
                //SetScale/P x 0,0.382878,"", g2t_map    //CCD search, 0.3 s.
                //SetScale/P x 0,0.67233,"", g2t_map    //CCD search, 0.6 s.
                //SetScale/P x 0,0.34,"", g2t_map            //time-series spectrum imaging, 0.1s per frame
                SetScale/P x 0,1.07,"", g2t_map          //CCD search, 1s per frame
                //SetScale/P x 0,3.07,"", g2t_map         //CCD search, 3 s.
                //SetScale/P x 0,0.171,"", g2t_map        //CCD search, 0.1 s.
                //SetScale/P x 0,2.075,"", g2t_map
                //SetScale/P x 0,1.574,"", g2t_map        //CCD search, 1.5 s per frame
                
                imagetransform sumallrows g2t_map
                duplicate/o W_sumRows g2t
                g2t=g2t/pcount
                //SetScale/P x 0,0.2192825,"", g2t
                //SetScale/P x 0,0.382878,"", g2t
                //SetScale/P x 0,0.67233,"", g2t
                //SetScale/P x 0,0.34,"", g2t
                SetScale/P x 0,1.07,"", g2t
                //SetScale/P x 0,3.07,"", g2t                 //CCD search, 3 s per frame.
                //SetScale/P x 0,0.171,"", g2t               //CCD search, 0.1 s per frame
                //SetScale/P x 0,2.075,"", g2t                 //CCD search, 2 s per frame.
                //SetScale/P x 0,1.574,"", g2t 			 //CCD search, 1.5 s per frame.
                
                duplicate/o g2t_map, $g2t_mapname
                duplicate/o g2t, $g2tname
                duplicate/o x_map, $x_mapname
                duplicate/o y_map, $y_mapname
                duplicate/o I_ave, $I_avename
                
                killwaves g2t_map, g2t, x_map, y_map, I_ave, add_td, bin_average
                print pcount
end

// Fitting resampled g2 along the log t axis.
Function fitg2_rs(g2t, endpoint, name)
		wave g2t
		variable endpoint
		string name
		
		variable A, dA, a1, b1, da1, db1, tau1, dtau1, b_1, db_1, tau, dtau, ptau
		string fitname="fit_"+name
		string fit1name="fit1_"+name
		string parametername="parameter_"+name
		string normname="g2_"+name+"_norm"
		string fitnormname="fit_"+name+"_norm"
		
		duplicate/o g2t temp
		redimension/n=(endpoint+1) temp
		
		Make/D/N=3/O W_coef
		W_coef[0] = {1,1,1}
		FuncFit/Q/NTHR=0 relaxationfunc3 W_coef  temp[0,endpoint]/X=t_rs[0, endpoint]
		A=W_coef[0]
		tau1=W_coef[1]
		b_1=W_coef[2]
		wave W_sigma=$"W_sigma"
		dA=W_sigma[0]
		dtau1=W_sigma[1]
		db_1=W_sigma[2]
		print A, dA
		duplicate/o $"fit_temp" $fit1name
		
		duplicate/o temp, y1
		y1[]=ln(-ln((temp[p]-1)/A)/2)
		make/o/n=(endpoint+1) x1
		x1[0,endpoint]=ln(0.248587083*p)  // Bin8, 200 pixels^2, 0.18s, 1000 frames.
		//x1[0,endpoint]=ln(0.0584*p)  // BNL, Orius, 0.0584 s
		//x1[1,endpoint]=ln(0.005*p)  // BNL, 400 f/s, adding 2 frames together.
		//x1[1,endpoint]=ln(0.171*p)   // CCD search, 0.1 s.
		//x1[1,endpoint]=ln(0.34*p)   // drift corrected time-series, 0.1 s.
		//x1[1,endpoint]=ln(0.2192825*p)   // 0.15 s
		//x1[1,endpoint]=ln(0.67233*p)   // 0.6 s.
		//x1[1,endpoint]=ln(3.07*p)         // CCD search, 3 s.
		//x1[1,endpoint]=ln(1.574*p)
		//x1[1,endpoint]=ln(0.2276*p)
		//x1[1,endpoint]=ln(0.2363*p)  // ERC, bin 4, 0.1 s per frame.
		//x1[1,endpoint]=ln(2.1414*p)
		Display y1 vs x1
		CurveFit/Q/M=2/W=0 line, y1[0,endpoint]/X=x1[0,endpoint]/D
		a1=W_coef[0]
		b1=W_coef[1]
		da1=W_sigma[0]
		db1=W_sigma[1]
		
		tau=exp(-a1/b1)
		dtau=tau/b1^2*(a1*db1-da1*b1)
		
		duplicate/o temp fit
		fit[]=1+A*exp(-2*(0.248587083*p/tau)^b1)
		//fit[]=1+A*exp(-2*(0.0584*p/tau)^b1)
		//fit[]=1+A*exp(-2*(0.005*p/tau)^b1)
		//fit[]=1+A*exp(-2*(0.171*p/tau)^b1)
		//fit[]=1+A*exp(-2*(0.34*p/tau)^b1)
		//fit[]=1+A*exp(-2*(0.2192825*p/tau)^b1)
		//fit[]=1+A*exp(-2*(0.67233*p/tau)^b1)
		//fit[]=1+A*exp(-2*(3.07*p/tau)^b1)
		//fit[]=1+A*exp(-2*(1.574*p/tau)^b1)
		//fit[]=1+A*exp(-2*(0.2276*p/tau)^b1)
		//fit[]=1+A*exp(-2*(0.2363*p/tau)^b1)
		//fit[]=1+A*exp(-2*(2.1414*p/tau)^b1)
		
		display temp
		ModifyGraph log(bottom)=1
		appendtograph fit_temp
		AppendToGraph fit
		print tau, dtau
		print b1, db1
		
		duplicate/o g2t temp, temp_norm
		temp_norm[]=(temp[p]-1)/A
		duplicate/o fit fit_norm
		fit_norm[]=(fit[p]-1)/A
		
		display temp_norm
		ModifyGraph log(bottom)=1
		appendtograph fit_norm
		
		ptau=round(tau/0.248587083)
		//ptau=round(tau/0.0584)
		//ptau=round(tau/0.005)
		//ptau=round(tau/0.171)
		//ptau=round(tau/0.34)
		//ptau=round(tau/0.2192825)
		//ptau=round(tau/0.67233)
		//ptau=round(tau/3.07)
		//ptau=round(tau/1.574)
		//ptau=round(tau/0.2276)
		//ptau=round(tau/0.2363)
		//ptau=round(tau/2.1414)
		print temp_norm[ptau]
		if (temp_norm[ptau]>exp(-2)*1.01)
			print "suspicisous fitting"
			else
				print "good fitting"
		endif
		
		duplicate/o fit $fitname
		duplicate/o temp_norm $normname
		duplicate/o fit_norm $fitnormname
		make/o/n=(10) parameter
		parameter[0]=tau
		parameter[1]=dtau
		parameter[2]=b1
		parameter[3]=db1
		parameter[4]=A
		parameter[5]=dA
		parameter[6]=tau1
		parameter[7]=dtau1
		parameter[8]=b_1
		parameter[9]=db_1
		duplicate/o parameter $parametername
		killwaves temp, x1, y1, W_coef, W_sigma, fit, fit_temp, temp_norm, fit_norm
end

Function fitg2(g2t, endpoint, name)
		wave g2t
		variable endpoint
		string name
		
		variable A, dA, a1, b1, da1, db1, tau1, dtau1, b_1, db_1, tau, dtau, ptau
		string fitname="fit_"+name
		string fit1name="fit1_"+name
		string parametername="parameter_"+name
		string normname="g2_"+name+"_norm"
		string fitnormname="fit_"+name+"_norm"
		
		duplicate/o g2t temp
		redimension/n=(endpoint+1) temp
		
		Make/D/N=3/O W_coef
		W_coef[0] = {1,10,1}
		FuncFit/Q/NTHR=0 relaxationfunc3 W_coef  temp[1,endpoint]
		A=W_coef[0]
		tau1=W_coef[1]
		b_1=W_coef[2]
		wave W_sigma=$"W_sigma"
		dA=W_sigma[0]
		dtau1=W_sigma[1]
		db_1=W_sigma[2]
		print A, dA
		duplicate/o $"fit_temp" $fit1name
		
		duplicate/o temp, y1
		y1[]=ln(-ln((temp[p]-1)/A)/2)
		make/o/n=(endpoint+1) x1
		x1[1,endpoint]=ln(0.0584*p)  // BNL, Orius, 0.05s.
		//x1[1,endpoint]=ln(0.005*p)  // BNL, 400 f/s, adding 2 frames together.
		//x1[1,endpoint]=ln(0.171*p)   // CCD search, 0.1 s.
		//x1[1,endpoint]=ln(0.34*p)   // drift corrected time-series, 0.1 s.
		//x1[1,endpoint]=ln(0.2192825*p)   // 0.15 s
		//x1[1,endpoint]=ln(0.67233*p)   // 0.6 s.
		//x1[1,endpoint]=ln(3.07*p)         // CCD search, 3 s.
		//x1[1,endpoint]=ln(1.574*p)
		//x1[1,endpoint]=ln(0.2276*p)
		//x1[1,endpoint]=ln(0.2363*p)  // ERC, bin 4, 0.1 s per frame.
		//x1[1,endpoint]=ln(2.1414*p)
		Display y1 vs x1
		CurveFit/Q/M=2/W=0 line, y1[1,endpoint]/X=x1[1,endpoint]/D
		a1=W_coef[0]
		b1=W_coef[1]
		da1=W_sigma[0]
		db1=W_sigma[1]
		
		tau=exp(-a1/b1)
		dtau=tau/b1^2*(a1*db1-da1*b1)
		
		duplicate/o temp fit
		fit[]=1+A*exp(-2*(0.0584*p/tau)^b1)
		//fit[]=1+A*exp(-2*(0.005*p/tau)^b1)
		//fit[]=1+A*exp(-2*(0.171*p/tau)^b1)
		//fit[]=1+A*exp(-2*(0.34*p/tau)^b1)
		//fit[]=1+A*exp(-2*(0.2192825*p/tau)^b1)
		//fit[]=1+A*exp(-2*(0.67233*p/tau)^b1)
		//fit[]=1+A*exp(-2*(3.07*p/tau)^b1)
		//fit[]=1+A*exp(-2*(1.574*p/tau)^b1)
		//fit[]=1+A*exp(-2*(0.2276*p/tau)^b1)
		//fit[]=1+A*exp(-2*(0.2363*p/tau)^b1)
		//fit[]=1+A*exp(-2*(2.1414*p/tau)^b1)
		
		display temp
		ModifyGraph log(bottom)=1
		appendtograph fit_temp
		AppendToGraph fit
		print tau, dtau
		print b1, db1
		
		duplicate/o g2t temp, temp_norm
		temp_norm[]=(temp[p]-1)/A
		duplicate/o fit fit_norm
		fit_norm[]=(fit[p]-1)/A
		
		display temp_norm
		ModifyGraph log(bottom)=1
		appendtograph fit_norm
		
		ptau=round(tau/0.0584)
		//ptau=round(tau/0.005)
		//ptau=round(tau/0.171)
		//ptau=round(tau/0.34)
		//ptau=round(tau/0.2192825)
		//ptau=round(tau/0.67233)
		//ptau=round(tau/3.07)
		//ptau=round(tau/1.574)
		//ptau=round(tau/0.2276)
		//ptau=round(tau/0.2363)
		//ptau=round(tau/2.1414)
		print temp_norm[ptau]
		if (temp_norm[ptau]>exp(-2)*1.01)
			print "suspicisous fitting"
			else
				print "good fitting"
		endif
		
		duplicate/o fit $fitname
		duplicate/o temp_norm $normname
		duplicate/o fit_norm $fitnormname
		make/o/n=(10) parameter
		parameter[0]=tau
		parameter[1]=dtau
		parameter[2]=b1
		parameter[3]=db1
		parameter[4]=A
		parameter[5]=dA
		parameter[6]=tau1
		parameter[7]=dtau1
		parameter[8]=b_1
		parameter[9]=db_1
		duplicate/o parameter $parametername
		killwaves temp, x1, y1, W_coef, W_sigma, fit, fit_temp, temp_norm, fit_norm
end
		
Function amorphous(g2t_map, I_ave, I_threshold, g2tname)
		wave g2t_map, I_ave
		variable I_threshold
		string g2tname
		
		variable n=numpnts(I_ave), count=0, i
		make/n=(1000) g2t=0
		for (i=0;i<n;i+=1)
			if (I_ave[i]<I_threshold)
				g2t[]=g2t[p]+g2t_map[p][i]
				count=count+1
			endif
		endfor
		
		g2t=g2t/count
		//SetScale/P x 0,3.07,"", g2t
		SetScale/P x 0,2.1414,"", g2t
		duplicate/o g2t, $g2tname
		killwaves g2t
		print count
end

Function correlatefunc(data,xstart,xend,ystart,yend)
	       wave data
	       variable xstart,xend,ystart,yend
	       variable xp = xend-xstart, yp = yend-ystart, zp = DimSize(data, 2)
	
               make/o/n=(xp, yp, zp) g2t
	       variable i, j, k, l, add, ave
	       string g2tname="g2t_"+nameofwave(data)+"_x"+num2istr(xstart)+"to"+num2istr(xend)+"_y"+num2istr(ystart)+"to"+num2istr(yend)
	       
	       for(i=0; i<xp; i+=1)
		         for(j=0; j<yp; j+=1)
			            Imagetransform/beam={(i+xstart), (j+ystart)} getbeam data
			            wave W_beam = $"W_beam"
			            ave=mean(W_beam)
			            for (k=0;k<zp;k+=1)
			                        add=0;
			                        for (l=0; l<zp-k; l+=1)
			                                     add=add+W_beam[l]*W_beam[l+k]
			                        endfor
                                                 g2t[i][j][k]=add/ave^2/(zp-k)
                                    endfor
                          endfor
                endfor
                
                rename g2t, $g2tname
end

// Calculate annular average of data acquired at BNL, cl 380mm, on Orius.
Function annular_average_BNLOrius(data, xcenter, ycenter, xa, ya, xb, xc, yc, yd, xe, ye, strip_width)
                 wave data
                 variable xcenter, ycenter, xa, ya, xb, xc, yc, yd, xe, ye, strip_width
                 
                 variable xp = DimSize(data, 0), yp = DimSize(data, 1), zp = Dimsize(data, 2)
                 variable i, radius
                 variable kbc=(yc-ya)/(xc-xb), kde=(ye-yd)/(xe-xc)
                 
                 string name_ave=NameofWave(data)+"_ave"
                 string name_av=NameofWave(data)+"_annular_av"
                 
                 make/o/n=(xp,yp) temp_image
                 
                 temp_image=data[p][q][0]
                 temp_image[xe, xa][ye, ya] = nan
                 temp_image[xb, xe][yd, ya]=nan
                 temp_image[xc, xb][yd, yc]=nan
                 temp_image[xc, xb][yc, kbc*(p-xc)+yc]=nan
                 temp_image[xc, xe][kde*(p-xc)+yd, yd]=nan
                 radius = AnnularAverage(temp_image, xcenter, ycenter, strip_width)
                 wave annular_av = $"annular_av"
                 
                 variable kp=dimsize(annular_av,0)
                 make/o/n=(kp,zp) annular_ave
                 annular_ave[][0]=annular_av[p]
                 
                 for (i=1; i<zp; i+=1)
                 	temp_image=data[p][q][0]
                 	temp_image[xe, xa][ye, ya] = nan
                 	temp_image[xb, xe][yd, ya]=nan
                 	temp_image[xc, xb][yd, yc]=nan
                 	temp_image[xc, xb][yc, kbc*(p-xc)+yc]=nan
                 	temp_image[xc, xe][kde*(p-xc)+yd, yd]=nan
                 	radius = AnnularAverage(temp_image, xcenter, ycenter, strip_width)
                 	annular_ave[][i] = annular_av[p]
                 endfor
                 
                 for (i=0;i<kp;i+=1)
                       matrixop/o temp=row(annular_ave,i)
                       annular_av[i]=mean(temp)
                 endfor
                 
                 SetScale/P x 0,0.00352388*2*strip_width,"", annular_ave // Orius, hardware bin 4, software bin 2.
                 SetScale/P x 0,0.00352388*2*strip_width,"", annular_av
                 duplicate/o annular_ave, $name_ave; killwaves annular_ave
                 duplicate/o annular_av, $name_av; killwaves annular_av
end

Function annular_average_BNL(data, xcenter, ycenter, xa, ya, yb, xc, yc, xd, xe, ye, strip_width)
                 wave data
                 variable xcenter, ycenter, xa, ya, yb, xc, yc, xd, xe, ye, strip_width
                 
                 variable xp = DimSize(data, 0), yp = DimSize(data, 1), zp = Dimsize(data, 2)
                 variable i, radius
                 variable kbc=(yc-yb)/(xc-xa), kde=(ye-yc)/(xe-xd)
                 
                 string name_ave=NameofWave(data)+"_ave"
                 string name_av=NameofWave(data)+"_annular_av"
                 
                 make/o/n=(xp,yp) temp_image
                 
                 temp_image=data[p][q][0]
                 temp_image[xc, xd][ya, yc] = nan
                 temp_image[xa, xc][ya,yb]=nan
                 temp_image[xa, xc][yb, kbc*(p-xa)+yb]=nan
                 temp_image[xd, xe][ya, ye]=nan
                 temp_image[xd, xe][ye, kde*(p-xd)+yc]= nan
                 radius = AnnularAverage(temp_image, xcenter, ycenter, strip_width)
                 wave annular_av = $"annular_av"
                 
                 variable kp=dimsize(annular_av,0)
                 make/o/n=(kp,zp) annular_ave
                 annular_ave[][0]=annular_av[p]
                 
                 for (i=1; i<zp; i+=1)
                       temp_image=data[p][q][i]
                       temp_image[xc, xd][ya,yc] = nan
                       temp_image[xa, xc][ya,yb]=nan
                       temp_image[xa, xc][yb, kbc*(p-xa)+yb]=nan
                       temp_image[xd, xe][ya, ye]=nan
                       temp_image[xd, xe][ye, kde*(p-xd)+yc]= nan
                       radius = AnnularAverage(temp_image, xcenter, ycenter, strip_width)
                       annular_ave[][i] = annular_av[p]
                 endfor
                 
                 for (i=0;i<kp;i+=1)
                       matrixop/o temp=row(annular_ave,i)
                       annular_av[i]=mean(temp)
                 endfor
                 
                 SetScale/P x 0,0.008052248*strip_width,"", annular_ave // K2, software binning 8.
                 SetScale/P x 0,0.008052248*strip_width,"", annular_av
                 //SetScale/P x 0,0.01006531*strip_width,"", annular_ave // K2, software binning 10.
                 //SetScale/P x 0,0.01006531*strip_width,"", annular_av
                 duplicate/o annular_ave, $name_ave; killwaves annular_ave
                 duplicate/o annular_av, $name_av; killwaves annular_av
end

Function annular_average_ERC(data, xcenter, ycenter, xa,ya,xb,yb,xc,xd,yd,xe,xf, strip_width)
                 wave data
                 variable xcenter, ycenter, xa,ya,xb,yb,xc,xd,yd,xe,xf,strip_width
                 
                 variable xp = DimSize(data, 0), yp = DimSize(data, 1), zp = Dimsize(data, 2)
                 variable i, radius
                 variable kab=(yb-ya)/(xb-xa), kcd=(yd-yb)/(xd-xc), kef=(ya-yd)/(xf-xe)
                 
                 string name_ave=NameofWave(data)+"_ave"
                 string name_av=NameofWave(data)+"_annular_av"
                 
                 make/o/n=(xp,yp) temp_image
                 
                 temp_image=data[p][q][0]
                 temp_image[xd,xe][yd,ya] = nan
                 temp_image[xc,xd][kcd*(p-xc)+yb,ya]=nan
                 temp_image[xb,xc][yb,ya]=nan
                 temp_image[xa,xb][kab*(p-xa)+ya, ya]=nan
                 temp_image[xe,xf][kef*(p-xe)+yd, ya]= nan
                 radius = AnnularAverage(temp_image, xcenter, ycenter, strip_width)
                 wave annular_av = $"annular_av"
                 
                 variable kp=dimsize(annular_av,0)
                 make/o/n=(kp,zp) annular_ave
                 annular_ave[][0]=annular_av[p]
                 
                 for (i=1; i<zp; i+=1)
                       temp_image=data[p][q][i]
                       temp_image[xd,xe][yd,ya] = nan
                       temp_image[xc,xd][kcd*(p-xc)+yb,ya]=nan
                       temp_image[xb,xc][yb,ya]=nan
                       temp_image[xa,xb][kab*(p-xa)+ya, ya]=nan
                       temp_image[xe,xf][kef*(p-xe)+yd, ya]= nan
                       radius = AnnularAverage(temp_image, xcenter, ycenter, strip_width)
                       annular_ave[][i] = annular_av[p]
                 endfor
                 
                 for (i=0;i<kp;i+=1)
                       matrixop/o temp=row(annular_ave,i)
                       annular_av[i]=mean(temp)
                 endfor
                 
                 SetScale/P x 0,0.006967*strip_width,"", annular_ave
                 SetScale/P x 0,0.006967*strip_width,"", annular_av
                 duplicate/o annular_ave, $name_ave; killwaves annular_ave
                 duplicate/o annular_av, $name_av; killwaves annular_av
end

Function annular_average(data, xcenter, ycenter, xblockstart, xblockend, yblockstart, yblockend, strip_width)
                 wave data
                 variable xcenter, ycenter, xblockstart, xblockend, yblockstart, yblockend,strip_width
                 
                 variable xp = DimSize(data, 0), yp = DimSize(data, 1), zp = Dimsize(data, 2)
                 variable i, radius
                 
                 string name_ave=NameofWave(data)+"_ave"
                 string name_av=NameofWave(data)+"_annular_av"
                 
                 make/o/n=(xp,yp) temp_image
                 
                 temp_image=data[p][q][0]
                 temp_image[xblockstart,xblockend][yblockstart,yblockend] = nan
                 radius = AnnularAverage(temp_image, xcenter, ycenter, strip_width)
                 wave annular_av = $"annular_av"
                 
                 variable kp=dimsize(annular_av,0)
                 make/o/n=(kp,zp) annular_ave
                 annular_ave[][0]=annular_av[p]
                 
                 for (i=1; i<zp; i+=1)
                       temp_image=data[p][q][i]
                       temp_image[xblockstart,xblockend][yblockstart,yblockend] = nan
                       radius = AnnularAverage(temp_image, xcenter, ycenter, strip_width)
                       annular_ave[][i] = annular_av[p]
                 endfor
                 
                 for (i=0;i<kp;i+=1)
                       matrixop/o temp=row(annular_ave,i)
                       annular_av[i]=mean(temp)
                 endfor
                 
                 SetScale/P x 0,0.000328214*strip_width,"", annular_ave   //bin 1, 2.15m. New gun. 06-15-15.
                 //SetScale/P x 0,0.00136495*strip_width,"", annular_ave   //bin 1, 510mm. New gun. 06-15-15.
                 //SetScale/P x 0,0.003044273*strip_width,"", $name_ave    //camera length 1.05 m, bin 4.
                 //SetScale/P x 0,0.003044273*2*strip_width,"", $name_ave    //camera length 1.05 m, bin 8.
                 //SetScale/P x 0,0.000487225*2*strip_width,"", $name_ave    //camera length 5.5 m, bin 8.
                 //SetScale/P x 0,0.000614807*2*strip_width,"", annular_ave  //bin 8, 4.4 m
                 //SetScale/P x 0,0.000614807*strip_width,"", annular_ave  //bin 4, 4.4 m
                 //SetScale/P x 0,0.000965769*2*strip_width,"", annular_ave  //bin 8, 2.7 m
                 //SetScale/P x 0,0.000965769*strip_width,"", annular_ave  //bin 4, 2.7 m
                 //SetScale/P x 0,0.001196894*strip_width,"", annular_ave  //bin 4, 2.15 m
                 //SetScale/P x 0,0.001318073*2*strip_width,"", annular_ave  //bin 8, 2.15 m
                 //SetScale/P x 0,0.003123587*strip_width,"", annular_ave
                 //SetScale/P x 0,0.00545597*2*strip_width,"", annular_ave     //bin 8, 510mm
                 //SetScale/P x 0,0.00545597*strip_width,"", annular_ave   //bin 4, 510mm
                 //SetScale/P x 0,0.00335707/2*strip_width,"", annular_ave   //bin 2, EFSTEM 840mm
                 //SetScale/P x 0,0.00335707*strip_width,"", annular_ave   //bin 4, EFSTEM 840mm
                 
                 SetScale/P x 0,0.000328214*strip_width,"", annular_av
                 //SetScale/P x 0,0.00136495*strip_width,"", annular_av
                 //SetScale/P x 0,0.003044273*strip_width,"", annular_av
                 //SetScale/P x 0,0.003044273*2*strip_width,"", annular_av
                 //SetScale/P x 0,0.000487225*2*strip_width,"", annular_av
                 //SetScale/P x 0,0.000614807*2*strip_width,"", annular_av
                 //SetScale/P x 0,0.000614807*strip_width,"", annular_av
                 //SetScale/P x 0,0.000965769*2*strip_width,"", annular_av
                 //SetScale/P x 0,0.000965769*strip_width,"", annular_av
                 //SetScale/P x 0,0.001196894*strip_width,"", annular_av
                 //SetScale/P x 0,0.001318073*2*strip_width,"", annular_av
                 //SetScale/P x 0,0.003123587*strip_width,"", annular_av
                 //SetScale/P x 0,0.00545597*2*strip_width,"", annular_av
                 //SetScale/P x 0,0.00545597*strip_width,"", annular_av
                 //SetScale/P x 0,0.00335707/2*strip_width,"", annular_av
                 //SetScale/P x 0,0.00335707*strip_width,"", annular_av
                 duplicate/o annular_ave, $name_ave; killwaves annular_ave
                 duplicate/o annular_av, $name_av; killwaves annular_av
end
                 
Function relaxation_annular(avedata,kmin,kmax,strip_width)
                 wave avedata
                 variable kmin,kmax,strip_width
                 
                 variable  kp=dimsize(avedata,0), zp = Dimsize(avedata, 1)
                 variable i,kpmin,kpmax
                 
                 string tauname="tau_"+nameofwave(avedata)+"_"+num2str(kmin)+" to "+num2str(kmax)
                 string betaname="beta_"+nameofwave(avedata)+"_"+num2str(kmin)+" to "+num2str(kmax)
                 string g2tname="g2t_"+nameofwave(avedata)+"_"+num2str(kmin)+" to "+num2str(kmax)
                              
                 kpmin=round(kmin/0.011150/strip_width)
                 kpmax=round(kmax/0.011150/strip_width)
                 kp=kpmax-kpmin
                 
                 make/o/n=(kp, zp) g2t
                 make/o/n=(kp) tau, beta_parameter
                 make/o/n=(zp) g2t_temp
                 
                 for (i=0;i<kp;i+=1)      
                       Imagetransform/G=(i+kpmin) getrow avedata
                       wave W_temp=$"W_ExtractedRow"
                       Correlate/auto W_temp, W_temp
                       g2t_temp[] =W_temp[p+zp-1]
                       
                       Setscale/P y, 0, 1, "", g2t_temp                        
                       //display g2t_temp
                       Make/o/D/N=2/O W_coef
                       W_coef[0] = {100,1}
                       Funcfit/N/Q=1/NTHR=0 relaxationfunc W_coef, g2t_temp
                       tau[i]=W_coef[0]
                       beta_parameter[i]=W_coef[1]
                       
                       g2t[i][]=g2t_temp[q]
                 endfor
                 
                 SetScale/P x kmin,0.011150*strip_width,"", g2t
                 SetScale/P x kmin,0.011150*strip_width,"", tau
                 SetScale/P x kmin,0.011150*strip_width,"", beta_parameter
                 duplicate/o tau, $tauname
                 duplicate/o beta_parameter, $betaname
                 duplicate/o g2t, $g2tname
                 killwaves tau, beta_parameter, g2t
                 
end

function C2tAnnularAv(c2t_st, cx, cy, xmin, xmax, ymin, ymax)
	wave c2t_st
	variable cx, cy, xmin, xmax, ymin, ymax
	
	variable zp = Dimsize(c2t_st, 2)
	
	ImageTransform/P = 0 getplane c2t_st
	wave M_imageplane = $"M_imageplane"
	variable radius = AnnularAverageCorners(M_ImagePlane, cx, cy, 2)
	make/o/n=(radius, zp) c2t_aav
	wave annular_av = $"annular_av"
	setscale/P x, 0, deltax(annular_av), c2t_aav
	Setscale/P y, 0, DimDelta(c2t_st, 2), c2t_aav
	
	variable i
	for(i=0; i<zp; i+=1)
		ImageTransform/P=(i) getplane c2t_st
		M_ImagePlane[xmin, xmax][ymin, ymax] = nan
		AnnularAverageCorners(M_ImagePlane, cx, cy, 2)
		c2t_aav[][i] = annular_av[p]
	endfor
	
end

function relaxation_k(data, k, allowance)
               wave data
               variable k, allowance
               
               variable ps=round((k-allowance)/0.01115), pe=round((k+allowance)/0.01115)
               variable kp=dimsize(data,0), tp=dimsize(data,1), i, j, add
               string g2tname=nameofwave(data)+"_"+num2str(k)+"_"+num2str(allowance)
               
               make/o/n=(tp) g2t
               for (i=0;i<tp;i+=1)
                       add=0
                       for (j=ps; j<pe; j+=1)
                                  add=add+data[j][i]
                       endfor
                       g2t[i]=add/(pe-ps)
               endfor
               
               duplicate/o g2t, $g2tname
               killwaves g2t
end

Function tau_k_map3(g2t_map, name)
		wave g2t_map
		string name
		
		variable i, j
		variable M=dimsize(g2t_map,0), N=dimsize(g2t_map, 1)
		string fitname, g2NormName, g2FitName, fitNormname, parameterName
		string tau_mapname="tau_map_"+name
		string terror_mapname="tauerror_map_"+name
		string beta_mapname="beta_map_"+name
		string betaerror_mapname="betaerror_map_"+name
		string A_mapname="A_map_"+name
		string Aerror_mapname="Aerror_map_"+name
		string g2tNormmapName="g2_norm_map_"+name
		string fitmapName="fit_map_"+name
		string fitNormmapName="fit_norm_map_"+name
		string rangeName="fit_range_map_"+name
		
		make/o/n=(N) t_map, terror_map, beta_map, betaerror_map, A_map, Aerror_map, range
		make/o/n=(M) g2t_temp
		SetScale/P x 0, 0.2276,"", g2t_temp
		duplicate/o g2t_map g2tNorm_map, fit_map, fitNorm_map
		g2tNorm_map=0
		fit_map=0
		fitNorm_map=0
		
		for (i=0; i<N; i+=1)
			g2t_temp[]=g2t_map[p][i]
			
			j=0
			do
				j=j+1
			while (g2t_temp[j]>g2t_temp[j+1]*0.997)
			range[i]=j

			if (range[i]<10)
				t_map[i]=nan
				terror_map[i]=nan
				beta_map[i]=nan
				betaerror_map[i]=nan
				A_map[i]=nan
				Aerror_map[i]=nan
				g2tNorm_map[][i]=nan
				fit_map[][i]=nan
				fitNorm_map[][i]=nan
			else
				fitname=num2istr(i)
				g2NormName="g2_norm"+fitname
				g2FitName="fit_"+fitname
				fitNormName="fit_norm"+fitname
				parameterName="parameter_"+fitname
				
				fitg2_noGraph(g2t_temp, range, fitname)
				
				wave par=$parameterName
				if (par[0]<range[i]*0.2276)
					t_map[i]=par[0]
					terror_map[i]=par[1]
					beta_map[i]=par[2]
					betaerror_map[i]=par[3]
					A_map[i]=par[4]
					Aerror_map[i]=par[5]
					
					wave g2tNorm=$g2NormName
					g2tNorm_map[][i]=g2tNorm[p]
					
					wave g2fit=$g2FitName
					fit_map[][i]=g2fit[p]
					
					wave fitNorm=$fitNormName
					fitNorm_map[][i]=fitNorm[p]
					
					killwaves par, g2tNorm, g2fit, fitNorm, $parameterName, $g2NormName, $g2FitName, $fitNormName
				else
					t_map[i]=nan
					terror_map[i]=nan
					beta_map[i]=nan
					betaerror_map[i]=nan
					A_map[i]=nan
					Aerror_map[i]=nan
					g2tNorm_map[][i]=nan
					fit_map[][i]=nan
					fitNorm_map[][i]=nan
				endif
			endif
		endfor
		
		duplicate/o t_map $tau_mapname
		duplicate/o terrror_map $terror_mapname
		duplicate/o beta_map $beta_mapname
		duplicate/o betaerror_map $betaerror_mapname
		duplicate/o A_map $A_mapname
		duplicate/o Aerror_map $Aerror_mapname
		duplicate/o g2tNorm_map $g2tNormmapName
		duplicate/o fit_map $fitmapName
		duplicate/o fitNorm_map $fitNormmapName
		duplicate/o range $rangeName
		killwaves g2t_temp, t_map, terror_map, beta_map, betaerror_map, A_map, Aerror_map, g2tNorm_map, fit_map, fitNorm_map,range
end

Function tau_k_map2(g2t_map, name)
		wave g2t_map
		string name
		
		variable i
		variable M=dimsize(g2t_map,0), N=dimsize(g2t_map, 1)
		string fitname, g2NormName, g2FitName, fitNormname, parameterName
		string tau_mapname="tau_map_"+name
		string terror_mapname="tauerror_map_"+name
		string beta_mapname="beta_map_"+name
		string betaerror_mapname="betaerror_map_"+name
		string A_mapname="A_map_"+name
		string Aerror_mapname="Aerror_map_"+name
		string g2tNormmapName="g2_norm_map_"+name
		string fitmapName="fit_map_"+name
		string fitNormmapName="fit_norm_map_"+name
		string rangeName="fit_range_map_"+name
		
		make/o/n=(N) t_map, terror_map, beta_map, betaerror_map, A_map, Aerror_map, range
		make/o/n=(M) g2t_temp
		SetScale/P x 0, 0.2276,"", g2t_temp
		duplicate/o g2t_map g2tNorm_map, fit_map, fitNorm_map
		g2tNorm_map=0
		fit_map=0
		fitNorm_map=0
		
		for (i=0; i<N; i+=1)
			g2t_temp[]=g2t_map[p][i]
			
			wavestats/Q/W/Z g2t_temp
			wave M_stats=$"M_WaveStats"
			range[i]=M_stats[13]

			if (range[i]<10)
				t_map[i]=nan
				terror_map[i]=nan
				beta_map[i]=nan
				betaerror_map[i]=nan
				A_map[i]=nan
				Aerror_map[i]=nan
				g2tNorm_map[][i]=nan
				fit_map[][i]=nan
				fitNorm_map[][i]=nan
			else
				fitname=num2istr(i)
				g2NormName="g2_norm"+fitname
				g2FitName="fit_"+fitname
				fitNormName="fit_norm"+fitname
				parameterName="parameter_"+fitname
				
				fitg2_noGraph(g2t_temp, range, fitname)
				
				wave par=$parameterName
				t_map[i]=par[0]
				terror_map[i]=par[1]
				beta_map[i]=par[2]
				betaerror_map[i]=par[3]
				A_map[i]=par[4]
				Aerror_map[i]=par[5]
				
				wave g2tNorm=$g2NormName
				g2tNorm_map[][i]=g2tNorm[p]
				
				wave g2fit=$g2FitName
				fit_map[][i]=g2fit[p]
				
				wave fitNorm=$fitNormName
				fitNorm_map[][i]=fitNorm[p]		
			endif
			
			killwaves par, g2tNorm, g2fit, fitNorm, $parameterName, $g2NormName, $g2FitName, $fitNormName
		endfor
		
		duplicate/o t_map $tau_mapname
		duplicate/o terrror_map $terror_mapname
		duplicate/o beta_map $beta_mapname
		duplicate/o betaerror_map $betaerror_mapname
		duplicate/o A_map $A_mapname
		duplicate/o Aerror_map $Aerror_mapname
		duplicate/o g2tNorm_map $g2tNormmapName
		duplicate/o fit_map $fitmapName
		duplicate/o fitNorm_map $fitNormmapName
		duplicate/o range $rangeName
		killwaves g2t_temp, t_map, tauerror_map, beta_map, betaerror_map, A_map, Aerror_map, par, g2tNorm, g2fit, fitNorm, g2tNorm_map, fit_map, fitNorm_map,range
end

Function fitg2_noGraph2(g2t, endpoint, name)
		wave g2t
		variable endpoint
		string name
		
		variable A, dA, a1, b1, da1, db1, tau1, dtau1, b_1, db_1, tau, dtau, ptau
		string fitname="fit_"+name
		//string fit1name="fit1_"+name
		string parametername="parameter_"+name
		string normname="g2_norm"+name
		string fitnormname="fit_norm"+name
		
		duplicate/o g2t temp
		redimension/n=(endpoint+1) temp
		
		Make/D/N=3/O W_coef
		W_coef[0] = {1,10,1}
		FuncFit/Q/NTHR=0 relaxationfunc3 W_coef  temp
		A=W_coef[0]
		tau1=W_coef[1]
		b_1=W_coef[2]
		wave W_sigma=$"W_sigma"
		dA=W_sigma[0]
		dtau1=W_sigma[1]
		db_1=W_sigma[2]
		//print A, dA
		//duplicate/o $"fit_temp" $fit1name
		
		duplicate/o temp, y1
		y1[]=ln(-ln((temp[p]-1)/A)/2)
		make/o/n=(endpoint+1) x1
		//x1[1,endpoint]=ln(3.07*p)
		//x1[1,endpoint]=ln(1.574*p)
		x1[1,endpoint]=ln(0.2276*p)
		//x1[1,endpoint]=ln(2.1414*p)
		//Display y1 vs x1
		CurveFit/Q/M=2/W=0 line, y1[1,endpoint]/X=x1[1,endpoint]/D
		a1=W_coef[0]
		b1=W_coef[1]
		da1=W_sigma[0]
		db1=W_sigma[1]
		
		tau=exp(-a1/b1)
		dtau=tau/b1^2*(a1*db1-da1*b1)
		
		duplicate/o temp fit
		//fit[]=1+A*exp(-2*(3.07*p/tau)^b1)
		//fit[]=1+A*exp(-2*(1.574*p/tau)^b1)
		fit[]=1+A*exp(-2*(0.2276*p/tau)^b1)
		//fit[]=1+A*exp(-2*(2.1414*p/tau)^b1)
		
		//display temp
		//ModifyGraph log(bottom)=1
		//appendtograph fit_temp
		//AppendToGraph fit
		//print tau, dtau
		//print b1, db1
		
		duplicate/o g2t temp, temp_norm
		temp_norm[]=(temp[p]-1)/A
		duplicate/o fit fit_norm
		fit_norm[]=(fit[p]-1)/A
		
		//display temp_norm
		//ModifyGraph log(bottom)=1
		//appendtograph fit_norm
		
		//ptau=round(tau/3.07)
		//ptau=round(tau/1.574)
		ptau=round(tau/0.2276)
		//ptau=round(tau/2.1414)
		//print temp_norm[ptau]
		//if (temp_norm[ptau]>exp(-2)*1.01)
			//print "suspicisous fitting"
			//else
				//print "good fitting"
		//endif
		
		duplicate/o fit $fitname
		duplicate/o temp_norm $normname
		duplicate/o fit_norm $fitnormname
		make/o/n=(10) parameter
		parameter[0]=tau
		parameter[1]=dtau
		parameter[2]=b1
		parameter[3]=db1
		parameter[4]=A
		parameter[5]=dA
		parameter[6]=tau1
		parameter[7]=dtau1
		parameter[8]=b_1
		parameter[9]=db_1
		duplicate/o parameter $parametername
		killwaves temp, x1, y1, W_coef, W_sigma, fit, fit_temp, temp_norm, fit_norm, parameter
end

Function fitg21(g2t, endpoint, name)
		wave g2t
		variable endpoint
		string name
		
		variable A, B, dA, dB, a1, b1, da1, db1, tau1, dtau1, b_1, db_1, tau, dtau, ptau
		string fitname="fit_"+name
		string fit1name="fit1_"+name
		string parametername="parameter_"+name
		string normname="g2_"+name+"_norm"
		string fitnormname="fit_"+name+"_norm"
		
		duplicate/o g2t temp
		redimension/n=(endpoint+1) temp
		
		Make/D/N=4/O W_coef
		W_coef[0] = {1,1,500,1}
		FuncFit/Q/NTHR=0 relaxationfunc W_coef  temp
		A=W_coef[0]
		B=W_coef[1]
		tau1=W_coef[2]
		b_1=W_coef[3]
		wave W_sigma=$"W_sigma"
		dA=W_sigma[0]
		dB=W_sigma[1]
		dtau1=W_sigma[2]
		db_1=W_sigma[3]
		print A, dA
		print B, dB
		duplicate/o $"fit_temp" $fit1name
		
		duplicate/o temp, y1
		y1[]=ln(-ln((temp[p]-A)/B)/2)
		make/o/n=(endpoint+1) x1
		x1[1,endpoint]=ln(3.07*p)
		Display y1 vs x1
		CurveFit/Q/M=2/W=0 line, y1[1,endpoint]/X=x1[1,endpoint]/D
		a1=W_coef[0]
		b1=W_coef[1]
		da1=W_sigma[0]
		db1=W_sigma[1]
		
		tau=exp(-a1/b1)
		dtau=tau/b1^2*(a1*db1-da1*b1)
		
		duplicate/o temp fit
		fit[]=A+B*exp(-2*(3.07*p/tau)^b1)
		
		display temp
		ModifyGraph log(bottom)=1
		appendtograph fit_temp
		AppendToGraph fit
		print tau, dtau
		print b1, db1
		
		duplicate/o g2t temp, temp_norm
		temp_norm[]=(temp[p]-A)/B
		duplicate/o fit fit_norm
		fit_norm[]=(fit[p]-A)/B
		
		display temp_norm
		ModifyGraph log(bottom)=1
		appendtograph fit_norm
		
		ptau=round(tau/3.07)
		print temp_norm[ptau]
		if (temp_norm[ptau]>exp(-2)*1.01)
			print "suspicisous fitting"
			else
				print "good fitting"
		endif
		
		duplicate/o fit $fitname
		duplicate/o temp_norm $normname
		duplicate/o fit_norm $fitnormname
		make/o/n=(12) parameter
		parameter[0]=tau
		parameter[1]=dtau
		parameter[2]=b1
		parameter[3]=db1
		parameter[4]=A
		parameter[5]=dA
		parameter[6]=B
		parameter[7]=dB
		parameter[8]=tau1
		parameter[9]=dtau1
		parameter[10]=b_1
		parameter[11]=db_1
		duplicate/o parameter $parametername
		killwaves temp, x1, y1, W_coef, W_sigma, fit, fit_temp, parameter, temp_norm, fit_norm
end

Function tau_k_map1(g2t_map, name)
		wave g2t_map
		string name
		
		variable i, range, V_levelX
		variable N=dimsize(g2t_map, 1)
		string fitname
		string tau_mapname="tau_map"+name
		string tauerror_mapname="tauerror_map"+name
		string beta_mapname="beta_map"+name
		string betaerror_mapname="betaerror_map"+name
		string A_mapname="A_map"+name
		string Aerror_mapname="Aerror_map"+name
		make/n=(N) t_map, tauerror_map, beta_map, betaerror_map, A_map, Aerror_map
		
		duplicate/o g2t_map g2t_temp
		redimension/n=(-1,0) g2t_temp
		
		for (i=0; i<N; i+=1)
			g2t_temp[]=g2t_map[p][i]
			wavestats/Q/W/Z g2t_temp
			wave M_stats=$"M_WaveStats"
			//Findlevel /EDGE=2/P/Q g2t_temp 1
			//range=min(M_stats[9], round(V_levelX))
			range=M_stats[9]
			//differentiate
			fitname=name+"_"+num2istr(i)
			fitg2(g2t_temp, range, fitname)
			wave parameter=$"parameter"
			t_map[i]=parameter[0]
			tauerror_map[i]=parameter[1]
			beta_map[i]=parameter[2]
			betaerror_map=parameter[3]
			A_map=parameter[4]
			Aerror_map=parameter[5]
		endfor
		
		duplicate/o t_map $tau_mapname
		duplicate/o tauerrror_map $tauerror_mapname
		duplicate/o beta_map $beta_mapname
		duplicate/o betaerror_map $betaerror_mapname
		duplicate/o A_map $A_mapname
		duplicate/o Aerror_map $Aerror_mapname
		killwaves g2t_temp, t_map, tauerror_map, beta_map, betaerror_map, A_map, Aerror_map, M_stats
end

Function correlatefunc_k_map2(data,k,bin,xcenter,ycenter,xblockstart,xblockend,yblockstart,yblockend)    // <I(t')>^2 from I(t') and I(t+t'). binning pixels of numbers bin^2.
	       wave data
	       variable k, bin, xcenter, ycenter,xblockstart,xblockend,yblockstart,yblockend
	       
	       variable r=round(k/0.01115), tp = DimSize(data, 2), neighbour=round((bin-1)/2)
	       variable xmin, xmax, y1, y2, temp, pixelsize=1000
	
               make/o/n=(tp) add_td=0
               make/o/n=(tp,pixelsize) g2t_map=0
               make/o/n=(pixelsize) x_map, y_map, I_ave
	       variable x, i, j
	       variable pcount=0
	       string name="_"+nameofwave(data)+"_"+num2str(k)+"_bin"+num2str(bin)
	       string g2t_mapname="g2t_map"+name
	       string g2tname="g2t"+name
	       string x_mapname="x_map"+name
	       string y_mapname="y_map"+name
	       string I_avename="I_ave"+name
	                                                                                                                                                                   // Find the circular points on diffraction patterns with the reciprocal space vector k.
	       xmin=xcenter-r+1
	       xmax=xcenter+r
	       for(x=xmin; x<xmax; x+=1)
		        temp=round(sqrt(r^2-(xcenter-x)^2))
		        y1=ycenter+temp
		        y2=ycenter-temp
		        if (x>xblockstart && x<xblockend && y1>yblockstart && y1<yblockend)
		              pcount=pcount
		       else
		              make/o/n=(tp) bin_average=0
		              for(i=0-neighbour; i<neighbour+1; i+=1)
		                     for(j=0-neighbour; j<neighbour+1; j+=1)
			                        Imagetransform/beam={(x+i), (y1+j)} getbeam data
			                        wave W_beam = $"W_beam"
			                        bin_average=bin_average+W_beam
			            endfor
			      endfor
			      W_beam=bin_average/bin^2
			      //ave=mean(W_beam)
			      add_td=0
			      for (j=0;j<tp;j+=1)
			            for (i=0;i<(tp-j);i+=1)
			                      add_td[j]=add_td[j]+W_beam[i]*W_beam[i+j]
			            endfor
			      endfor

			      g2t_map[][pcount]=(tp-p)*add_td[p]/sum(W_beam,0,tp-p-1)/sum(W_beam, p, tp-1)
			      x_map[pcount]=x
			      y_map[pcount]=y1
			      I_ave[pcount]=mean(W_beam)
			      
                               pcount=pcount+1
                        endif
                        
                        if (x>xblockstart && x<xblockend && y2>yblockstart && y2<yblockend)
		             pcount=pcount
		       else
		             make/o/n=(tp) bin_average=0
		             for(i=0-neighbour; i<neighbour+1; i+=1)
		                     for(j=0-neighbour; j<neighbour+1; j+=1)
			                        Imagetransform/beam={(x+i), (y2+j)} getbeam data
			                        wave W_beam = $"W_beam"
			                        bin_average=bin_average+W_beam
			            endfor
			      endfor
			      W_beam=bin_average/bin^2
			      //ave=mean(W_beam)
			      add_td=0
			      for (j=0;j<tp;j+=1)
			            for (i=0;i<(tp-j);i+=1)
			                      add_td[j]=add_td[j]+W_beam[i]*W_beam[i+j]
			            endfor
			      endfor

			      g2t_map[][pcount]=(tp-p)*add_td[p]/sum(W_beam,0,tp-p-1)/sum(W_beam, p, tp-1)
			      x_map[pcount]=x
			      y_map[pcount]=y2
			      I_ave[pcount]=mean(W_beam)
			     
                               pcount=pcount+1
                      endif
                endfor
                
                redimension/n=(-1, pcount) g2t_map
                redimension/n=(pcount) x_map, y_map, I_ave
                SetScale/P x 0,0.34*1000/tp,"", g2t_map
                
                imagetransform sumallrows g2t_map
                duplicate/o W_sumRows g2t
                g2t=g2t/pcount
                SetScale/P x 0,0.34*1000/tp,"", g2t
                
                duplicate/o g2t_map, $g2t_mapname
                duplicate/o g2t, $g2tname
                duplicate/o x_map, $x_mapname
                duplicate/o y_map, $y_mapname
                duplicate/o I_ave, $I_avename
                
                killwaves g2t_map, g2t, x_map, y_map, I_ave, add_td, bin_average
                print pcount
end

Function relaxation_2D(data,xstart,xend,ystart,yend)         //under construction.
                wave data
                variable xstart,xend,ystart,yend
                variable xp=xend-xstart+1, yp=yend-ystart+1, zp = DimSize(data, 2)
                variable i, j,k,l
                string outputName=nameofwave(data)+"_"+num2istr(xstart)+"_"+num2istr(ystart)+"_"+num2istr(xp)+"by"+num2istr(yp)
                string g2tname="g2t_"+outputName
                string tauname="tau_"+outputName
                string Aname="A_"+outputName
                
                make/o/n=(xp, yp, zp) g2t=0
                make/o/n=(xp,yp) tau=0
                make/o/n=(xp,yp) A=0
                make/o/n=(zp) add_td=0, g2t_temp=0
                
                for (i=0; i<xp; i+=1)
                        for (j=0; j<yp; j+=1)
                                imagetransform/beam={(i+xstart), (j+ystart)} getbeam data
                                wave W_beam = $"W_beam"
                                add_td=0
                                g2t_temp=0
			       for (k=0;k<zp;k+=1)
			            for (l=0;l<(zp-k);l+=1)
			                      add_td[k]=add_td[k]+W_beam[l]*W_beam[k+l]
			            endfor
			       endfor
			       g2t_temp[]=(zp-p)*add_td[p]/sum(W_beam,0,zp-p-1)/sum(W_beam, p, zp-1)
			       
			       SetScale/P x 0,0.34,"", g2t_temp
			       g2t[i][j][]=g2t_temp[r]
			       
                                Make/o/D/N=2/O W_coef
                                W_coef[0] = {1,100}
                                Funcfit/NTHR=0/Q relaxationfunc W_coef, g2t_temp
                                tau[i][j]=W_coef[1]
                                A[i][j]=W_coef[0]
                         endfor
                endfor
                 
                 duplicate/o g2t, $g2tname
                 duplicate/o tau, $tauname
                 duplicate/o A, $Aname
end

Function correlatefunc_k(data,k,bin,xcenter,ycenter,xblockstart,xblockend,yblockstart,yblockend)    // <I(t')>^2 only from I(t') and I(t+t'). binning pixels of numbers bin^2.
	       wave data
	       variable k, bin, xcenter, ycenter,xblockstart,xblockend,yblockstart,yblockend
	       
	       variable r=round(k/0.01115), tp = DimSize(data, 2), neighbour=round((bin-1)/2)
	       variable xmin, xmax, y1, y2, temp
	
               make/o/n=(tp) g2t=0, add_td=0
	       variable x, i, j
	       variable pcount=0
	       string g2tname="g2t_"+nameofwave(data)+"_"+num2str(k)+"_bin_"+num2str(bin)
	                                                                                                                                                                   // Find the circular points on diffraction patterns with the reciprocal space vector k.
	       xmin=xcenter-r+1
	       xmax=xcenter+r
	       for(x=xmin; x<xmax; x+=1)
		        temp=round(sqrt(r^2-(xcenter-x)^2))
		        y1=ycenter+temp
		        y2=ycenter-temp
		        if (x>xblockstart && x<xblockend && y1>yblockstart && y1<yblockend)
		              g2t=g2t
		       else
		              make/o/n=(tp) bin_average=0
		              for(i=0-neighbour; i<neighbour+1; i+=1)
		                     for(j=0-neighbour; j<neighbour+1; j+=1)
			                        Imagetransform/beam={(x+i), (y1+j)} getbeam data
			                        wave W_beam = $"W_beam"
			                        bin_average=bin_average+W_beam
			            endfor
			      endfor
			      W_beam=bin_average/bin^2
			      //ave=mean(W_beam)
			      add_td=0
			      for (j=0;j<tp;j+=1)
			            for (i=0;i<(tp-j);i+=1)
			                      add_td[j]=add_td[j]+W_beam[i]*W_beam[i+j]
			            endfor
			      endfor

			      g2t[]=g2t[p]+(tp-p)*add_td[p]/sum(W_beam,0,tp-p-1)/sum(W_beam, p, tp-1)

                               pcount=pcount+1
                        endif
                        
                        if (x>xblockstart && x<xblockend && y2>yblockstart && y2<yblockend)
		             g2t=g2t
		       else
		             make/o/n=(tp) bin_average=0
		             for(i=0-neighbour; i<neighbour+1; i+=1)
		                     for(j=0-neighbour; j<neighbour+1; j+=1)
			                        Imagetransform/beam={(x+i), (y2+j)} getbeam data
			                        wave W_beam = $"W_beam"
			                        bin_average=bin_average+W_beam
			            endfor
			      endfor
			      W_beam=bin_average/bin^2
			      //ave=mean(W_beam)
			      add_td=0
			      for (j=0;j<tp;j+=1)
			            for (i=0;i<(tp-j);i+=1)
			                      add_td[j]=add_td[j]+W_beam[i]*W_beam[i+j]
			            endfor
			      endfor

			      g2t[]=g2t[p]+(tp-p)*add_td[p]/sum(W_beam,0,tp-p-1)/sum(W_beam, p, tp-1)
			     
                               pcount=pcount+1
                      endif
                endfor
                
                g2t=g2t/pcount
                SetScale/P x 0,0.34,"", g2t
                duplicate/o g2t, $g2tname
                killwaves g2t
                print pcount
end

Function correlatefunc_k6(data,k,bin,xcenter,ycenter,xblockstart,xblockend,yblockstart,yblockend)    // <I(t')>^2 only from I(t') and I(t+t'). binning pixels of numbers bin^2.
	       wave data
	       variable k, bin, xcenter, ycenter,xblockstart,xblockend,yblockstart,yblockend
	       
	       variable r=round(k/0.01115), tp = DimSize(data, 2), neighbour=round((bin-1)/2)
	       variable xmin, xmax, y1, y2, temp
	
               make/o/n=(tp) g2t=0
	       variable x, i, j, add_td, add_t, ave
	       variable pcount=0
	       string g2tname="g2t_"+nameofwave(data)+"_"+num2str(k)+"_bin_"+num2str(bin)
	                                                                                                                                                                   // Find the circular points on diffraction patterns with the reciprocal space vector k.
	       xmin=xcenter-r+1
	       xmax=xcenter+r
	       for(x=xmin; x<xmax; x+=1)
		        temp=round(sqrt(r^2-(xcenter-x)^2))
		        y1=ycenter+temp
		        y2=ycenter-temp
		        if (x>xblockstart && x<xblockend && y1>yblockstart && y1<yblockend)
		              g2t=g2t
		       else
		              make/o/n=(tp) bin_average=0
		              for(i=0-neighbour; i<neighbour+1; i+=1)
		                     for(j=0-neighbour; j<neighbour+1; j+=1)
			                        Imagetransform/beam={(x+i), (y1+j)} getbeam data
			                        wave W_beam = $"W_beam"
			                        bin_average=bin_average+W_beam
			            endfor
			      endfor
			      W_beam=bin_average/bin^2
			      //ave=mean(W_beam)
			      for (i=0;i<tp;i+=1)
			            add_td=0
			            add_t=0
			            for (j=0; j<(tp-i); j+=1)
			                        add_td=add_td+W_beam[j]*W_beam[j+i]
			                        add_t=add_t+W_beam[j]+W_beam[j+i]
			            endfor
			            ave=add_t/(tp-i)/2
                                     g2t[i]=g2t[i]+add_td/(tp-i)/ave^2
                              endfor
                              pcount=pcount+1
                        endif
                        if (x>xblockstart && x<xblockend && y2>yblockstart && y2<yblockend)
		             g2t=g2t
		       else
		             make/o/n=(tp) bin_average=0
		             for(i=0-neighbour; i<neighbour+1; i+=1)
		                     for(j=0-neighbour; j<neighbour+1; j+=1)
			                        Imagetransform/beam={(x+i), (y2+j)} getbeam data
			                        wave W_beam = $"W_beam"
			                        bin_average=bin_average+W_beam
			            endfor
			      endfor
			      W_beam=bin_average/bin^2
			      //ave=mean(W_beam)
			      for (i=0;i<tp;i+=1)
			             add_td=0
			             add_t=0
			             for (j=0; j<(tp-i); j+=1)
			                       add_td=add_td+W_beam[j]*W_beam[j+i]
			                       add_t=add_t+W_beam[j]+W_beam[j+i]
			             endfor
			             ave=add_t/(tp-i)/2
                                      g2t[i]=g2t[i]+add_td/(tp-i)/ave^2
                              endfor
                              pcount=pcount+1
                      endif
                endfor
                
                g2t=g2t/pcount
                SetScale/P x 0,0.34,"", g2t
                duplicate/o g2t, $g2tname
                killwaves g2t
                print pcount
end

Function correlatefunc_k5(data,k,bin,xcenter,ycenter,xblockstart,xblockend,yblockstart,yblockend)    // <I(t')>^2 over all delay time points. binning pixels of numbers bin^2.
	       wave data
	       variable k, bin, xcenter, ycenter,xblockstart,xblockend,yblockstart,yblockend
	       
	       variable r=round(k/0.01115), tp = DimSize(data, 2), neighbour=round((bin-1)/2)
	       variable xmin, xmax, y1, y2, temp
	
               make/o/n=(tp) g2t=0
	       variable x, i, j, add_td, ave
	       variable pcount=0
	       string g2tname="g2t_"+nameofwave(data)+"_"+num2str(k)+"_bin_"+num2str(bin)
	                                                                                                                                                                   // Find the circular points on diffraction patterns with the reciprocal space vector k.
	       xmin=xcenter-r+1
	       xmax=xcenter+r
	       for(x=xmin; x<xmax; x+=1)
		        temp=round(sqrt(r^2-(xcenter-x)^2))
		        y1=ycenter+temp
		        y2=ycenter-temp
		        if (x>xblockstart && x<xblockend && y1>yblockstart && y1<yblockend)
		              g2t=g2t
		       else
		              make/o/n=(tp) bin_average=0
		              for(i=0-neighbour; i<neighbour+1; i+=1)
		                     for(j=0-neighbour; j<neighbour+1; j+=1)
			                        Imagetransform/beam={(x+i), (y1+j)} getbeam data
			                        wave W_beam = $"W_beam"
			                        bin_average=bin_average+W_beam
			            endfor
			      endfor
			      W_beam=bin_average/bin^2
			      ave=mean(W_beam)
			      for (i=0;i<tp;i+=1)
			            add_td=0
			            for (j=0; j<(tp-i); j+=1)
			                        add_td=add_td+W_beam[j]*W_beam[j+i]
			            endfor
                                     g2t[i]=g2t[i]+add_td/(tp-i)/ave^2
                              endfor
                              pcount=pcount+1
                        endif
                        if (x>xblockstart && x<xblockend && y2>yblockstart && y2<yblockend)
		             g2t=g2t
		       else
		             make/o/n=(tp) bin_average=0
		             for(i=0-neighbour; i<neighbour+1; i+=1)
		                     for(j=0-neighbour; j<neighbour+1; j+=1)
			                        Imagetransform/beam={(x+i), (y2+j)} getbeam data
			                        wave W_beam = $"W_beam"
			                        bin_average=bin_average+W_beam
			            endfor
			      endfor
			      W_beam=bin_average/bin^2
			      ave=mean(W_beam)
			      for (i=0;i<tp;i+=1)
			             add_td=0
			             for (j=0; j<(tp-i); j+=1)
			                       add_td=add_td+W_beam[j]*W_beam[j+i]
			             endfor
                                      g2t[i]=g2t[i]+add_td/(tp-i)/ave^2
                              endfor
                              pcount=pcount+1
                      endif
                endfor
                
                g2t=g2t/pcount
                SetScale/P x 0,0.34,"", g2t
                duplicate/o g2t, $g2tname
                killwaves g2t
                print pcount
end
.
Function correlatefunc_k4(data,k,allowance,xcenter,ycenter,xblockstart,xblockend,yblockstart,yblockend)    // <I(t')>^2 over all delay time points. no binning.
	       wave data
	       variable k, allowance, xcenter, ycenter,xblockstart,xblockend,yblockstart,yblockend
	       
	       variable rmin=round((k-allowance)/0.01115), rmax=round((k+allowance)/0.01115), tp = DimSize(data, 2)
	       variable xmin, xmax, y1, y2, temp
	
               make/o/n=(tp) g2t=0
	       variable r, x, i, j, add_td, ave
	       variable pcount=0
	       string g2tname="g2t_"+nameofwave(data)+"_"+num2str(k)+"_"+num2str(allowance)
	                                                                                                                                                                   // Find the circular points on diffraction patterns with the reciprocal space vector k.
	       for(r=rmin; r<rmax; r+=1)
	                 xmin=xcenter-r+1
	                 xmax=xcenter+r
		         for(x=xmin; x<xmax; x+=1)
		                     temp=round(sqrt(r^2-(xcenter-x)^2))
		                     y1=ycenter+temp
		                     y2=ycenter-temp
		                     if (x>xblockstart && x<xblockend && y1>yblockstart && y1<yblockend)
		                                 g2t=g2t
		                     else
			                        Imagetransform/beam={(x), (y1)} getbeam data
			                        wave W_beam = $"W_beam"
			                        ave=mean(W_beam)
			                        for (i=0;i<tp;i+=1)
			                                   add_td=0
			                                   //add_t=0
			                                   for (j=0; j<(tp-i); j+=1)
			                                             add_td=add_td+W_beam[j]*W_beam[j+i]
			                                             //add_t=add_t+W_beam[j]+W_beam[j+i]
			                                   endfor
			                                   //ave=add_t/(tp-i)/2
                                                            g2t[i]=g2t[i]+add_td/(tp-i)/ave^2
                                                endfor
                                                pcount=pcount+1
                                     endif
                                     if (x>xblockstart && x<xblockend && y2>yblockstart && y2<yblockend)
		                                 g2t=g2t
		                     else
			                        Imagetransform/beam={(x), (y2)} getbeam data
			                        wave W_beam = $"W_beam"
			                        ave=mean(W_beam)
			                        for (i=0;i<tp;i+=1)
			                                   add_td=0
			                                   //add_t=0
			                                   for (j=0; j<(tp-i); j+=1)
			                                             add_td=add_td+W_beam[j]*W_beam[j+i]
			                                             //add_t=add_t+W_beam[j]+W_beam[j+i]
			                                   endfor
			                                   //ave=add_t/(tp-i)/2
                                                            g2t[i]=g2t[i]+add_td/(tp-i)/ave^2
                                                endfor
                                                pcount=pcount+1
                                     endif
                          endfor
                endfor
                
                g2t=g2t/pcount
                SetScale/P x 0,0.34,"", g2t
                duplicate/o g2t, $g2tname
                killwaves g2t
                print pcount
end

Function relaxationfunc2(w,t) : FitFunc
	Wave w
	Variable t

	//CurveFitDialog/ These comments were created by the Curve Fitting dialog. Altering them will
	//CurveFitDialog/ make the function less convenient to work with in the Curve Fitting dialog.
	//CurveFitDialog/ Equation:
	//CurveFitDialog/ f(t) =1+A*exp(-2*t/tau)
	//CurveFitDialog/ End of Equation
	//CurveFitDialog/ Independent Variables 1
	//CurveFitDialog/ t
	//CurveFitDialog/ Coefficients 2
	//CurveFitDialog/ w[0] = A
	//CurveFitDialog/ w[1] = tau
	
	return 1+w[0]*exp(-2*t/w[1])
End

// Fitting with two f(t)=f(0)exp(-(t/tau)^beta) terms, a change to KWW equation.
Function relaxationfunc4(w,t) : FitFunc
	Wave w
	Variable t

	//CurveFitDialog/ These comments were created by the Curve Fitting dialog. Altering them will
	//CurveFitDialog/ make the function less convenient to work with in the Curve Fitting dialog.
	//CurveFitDialog/ Equation:
	//CurveFitDialog/ f(t) =1+A*exp(-2*(t/tau1)^beta1)+A*exp(-2*(t/tau2)^beta2)
	//CurveFitDialog/ End of Equation
	//CurveFitDialog/ Independent Variables 1
	//CurveFitDialog/ t
	//CurveFitDialog/ Coefficients 5
	//CurveFitDialog/ w[0] = A
	//CurveFitDialog/ w[1] = tau1
	//CurveFitDialog/ w[2] = beta1
	//CurveFitDialog/ w[3] = tau2
	//CurveFitDialog/ w[4] = beta2
	
	return 1+w[0]*exp(-2*(t/w[1])^w[2])+w[0]*exp(-2*(t/w[3])^w[4])
End

// Fitting with two expoential terms, a change to KWW equation.
Function relaxationfunc5(w,t) : FitFunc
	Wave w
	Variable t

	//CurveFitDialog/ These comments were created by the Curve Fitting dialog. Altering them will
	//CurveFitDialog/ make the function less convenient to work with in the Curve Fitting dialog.
	//CurveFitDialog/ Equation:
	//CurveFitDialog/ f(t) =1+A1*exp(-2*(t/tau1)^beta1)+A2*exp(-2*(t/tau2)^beta2)
	//CurveFitDialog/ End of Equation
	//CurveFitDialog/ Independent Variables 1
	//CurveFitDialog/ t
	//CurveFitDialog/ Coefficients 6
	//CurveFitDialog/ w[0] = A1
	//CurveFitDialog/ w[1] = tau1
	//CurveFitDialog/ w[2] = beta1
	//CurveFitDialog/ w[3] = A2
	//CurveFitDialog/ w[4] = tau2
	//CurveFitDialog/ w[5] = beta2
	
	return 1+w[0]*exp(-2*(t/w[1])^w[2])+w[3]*exp(-2*(t/w[4])^w[5])
End

Function relaxationfunc3(w,t) : FitFunc
	Wave w
	Variable t

	//CurveFitDialog/ These comments were created by the Curve Fitting dialog. Altering them will
	//CurveFitDialog/ make the function less convenient to work with in the Curve Fitting dialog.
	//CurveFitDialog/ Equation:
	//CurveFitDialog/ f(t) =1+A*exp(-2*(t/tau)^beta)
	//CurveFitDialog/ End of Equation
	//CurveFitDialog/ Independent Variables 1
	//CurveFitDialog/ t
	//CurveFitDialog/ Coefficients 3
	//CurveFitDialog/ w[0] = A
	//CurveFitDialog/ w[1] = tau
	//CurveFitDialog/ w[2] = beta
	
	return 1+w[0]*exp(-2*(t/w[1])^w[2])
End

Function relaxationfunc1(w,t) : FitFunc
	Wave w
	Variable t

	//CurveFitDialog/ These comments were created by the Curve Fitting dialog. Altering them will
	//CurveFitDialog/ make the function less convenient to work with in the Curve Fitting dialog.
	//CurveFitDialog/ Equation:
	//CurveFitDialog/ f(t) =A+B*exp(-2*t/tau)
	//CurveFitDialog/ End of Equation
	//CurveFitDialog/ Independent Variables 1
	//CurveFitDialog/ t
	//CurveFitDialog/ Coefficients 3
	//CurveFitDialog/ w[0] = A
	//CurveFitDialog/ w[1] = B
	//CurveFitDialog/ w[2] = tau
	
	return w[0]+w[1]*exp(-2*t/w[2])
End

Function relaxationfunc(w,t) : FitFunc
	Wave w
	Variable t

	//CurveFitDialog/ These comments were created by the Curve Fitting dialog. Altering them will
	//CurveFitDialog/ make the function less convenient to work with in the Curve Fitting dialog.
	//CurveFitDialog/ Equation:
	//CurveFitDialog/ f(t) =A+B*exp(-2*(t/tau)^beta)
	//CurveFitDialog/ End of Equation
	//CurveFitDialog/ Independent Variables 1
	//CurveFitDialog/ t
	//CurveFitDialog/ Coefficients 4
	//CurveFitDialog/ w[0] = A
	//CurveFitDialog/ w[1] = B
	//CurveFitDialog/ w[2] = tau
	//CurveFitDialog/ w[3] = beta
	
	return w[0]+w[1]*exp(-2*(t/w[2])^w[3])
End

Function VTF(w,T) : FitFunc
	Wave w
	Variable T

	//CurveFitDialog/ These comments were created by the Curve Fitting dialog. Altering them will
	//CurveFitDialog/ make the function less convenient to work with in the Curve Fitting dialog.
	//CurveFitDialog/ Equation:
	//CurveFitDialog/ f(t) =A*exp(B/(T-570))
	//CurveFitDialog/ End of Equation
	//CurveFitDialog/ Independent Variables 1
	//CurveFitDialog/ T
	//CurveFitDialog/ Coefficients 2
	//CurveFitDialog/ w[0] = A
	//CurveFitDialog/ w[1] = B
	
	return w[0]*exp(W[1]/(T-570))
End

Function VTF1(w,T) : FitFunc
	Wave w
	Variable T

	//CurveFitDialog/ These comments were created by the Curve Fitting dialog. Altering them will
	//CurveFitDialog/ make the function less convenient to work with in the Curve Fitting dialog.
	//CurveFitDialog/ Equation:
	//CurveFitDialog/ f(t) =A*exp(B/(T-T0))
	//CurveFitDialog/ End of Equation
	//CurveFitDialog/ Independent Variables 1
	//CurveFitDialog/ T
	//CurveFitDialog/ Coefficients 3
	//CurveFitDialog/ w[0] = A
	//CurveFitDialog/ w[1] = B
	//CurveFitDialog/ w[2] = T0
	
	return w[0]*exp(W[1]/(T-w[2]))
End

Function average_k(data,k,allowance,xcenter,ycenter,xblockstart,xblockend,yblockstart,yblockend)
	       wave data
	       variable k, allowance, xcenter, ycenter,xblockstart,xblockend,yblockstart,yblockend
	       
	       variable rmin=round((k-allowance)/0.01115), rmax=round((k+allowance)/0.01115), tp = DimSize(data, 2)
	       variable xmin, xmax, y1, y2, temp
	
               make/o/n=(tp) I_ave=0
	       variable r, x
	       variable pcount=0
	       string I_avename="I_ave_"+nameofwave(data)+"_"+num2str(k)+"_"+num2str(allowance)
	       
	       for(r=rmin; r<rmax; r+=1)
	                 xmin=xcenter-r+1
	                 xmax=xcenter+r
		         for(x=xmin; x<xmax; x+=1)
		                     temp=round(sqrt(r^2-(xcenter-x)^2))
		                     y1=ycenter+temp
		                     y2=ycenter-temp
		                     if (x>xblockstart && x<xblockend && y1>yblockstart && y1<yblockend)
		                                 I_ave=I_ave
		                     else
			                        Imagetransform/beam={(x), (y1)} getbeam data
			                        wave W_beam = $"W_beam"
			                        I_ave=I_ave+W_beam
			                        pcount=pcount+1
			            endif
			            if (x>xblockstart && x<xblockend && y2>yblockstart && y2<yblockend)
		                                 I_ave=I_ave
		                     else
			                        Imagetransform/beam={(x), (y1)} getbeam data
			                        wave W_beam = $"W_beam"
			                        I_ave=I_ave+W_beam
			                        pcount=pcount+1
			            endif
                          endfor
                endfor
                
                I_ave=I_ave/pcount
                SetScale/P x 0,0.34,"", I_ave
                duplicate/o I_ave, $I_avename
                killwaves I_ave
                print pcount
end

Function correlatefunc_k2(data,k,allowance,xcenter,ycenter,xblockstart,xblockend,yblockstart,yblockend) // aveage <I(t')> only from I(t') and I(t+t')
	       wave data
	       variable k, allowance, xcenter, ycenter,xblockstart,xblockend,yblockstart,yblockend
	       
	       variable rmin=round((k-allowance)/0.01115), rmax=round((k+allowance)/0.01115), tp = DimSize(data, 2)
	       variable xmin, xmax, y1, y2, temp
	
               make/o/n=(tp) g2t=0
	       variable r, x, i, j, add_td, add_t, ave
	       variable pcount=0
	       string g2tname="g2t_"+nameofwave(data)+"_"+num2str(k)+"_"+num2str(allowance)+"_2"
	       
	       for(r=rmin; r<rmax; r+=1)
	                 xmin=xcenter-r+1
	                 xmax=xcenter+r
		         for(x=xmin; x<xmax; x+=1)
		                     temp=round(sqrt(r^2-(xcenter-x)^2))
		                     y1=ycenter+temp
		                     y2=ycenter-temp
		                     if (x>xblockstart && x<xblockend && y1>yblockstart && y1<yblockend)
		                                 g2t=g2t
		                     else
			                        Imagetransform/beam={(x), (y1)} getbeam data
			                        wave W_beam = $"W_beam"
			                        //ave=mean(W_beam)
			                        for (i=0;i<tp;i+=1)
			                                   add_td=0
			                                   add_t=0
			                                   for (j=0; j<(tp-i); j+=1)
			                                             add_td=add_td+W_beam[j]*W_beam[j+i]
			                                             add_t=add_t+W_beam[j]+W_beam[j+i]
			                                   endfor
			                                   ave=add_t/(tp-i)/2
                                                            g2t[i]=g2t[i]+add_td/(tp-i)/ave^2
                                                endfor
                                                pcount=pcount+1
                                     endif
                                     if (x>xblockstart && x<xblockend && y2>yblockstart && y2<yblockend)
		                                 g2t=g2t
		                     else
			                        Imagetransform/beam={(x), (y2)} getbeam data
			                        wave W_beam = $"W_beam"
			                        //ave=mean(W_beam)
			                        for (i=0;i<tp;i+=1)
			                                   add_td=0
			                                   add_t=0
			                                   for (j=0; j<(tp-i); j+=1)
			                                             add_td=add_td+W_beam[j]*W_beam[j+i]
			                                             add_t=add_t+W_beam[j]+W_beam[j+i]
			                                   endfor
			                                   ave=add_t/(tp-i)/2
                                                            g2t[i]=g2t[i]+add_td/(tp-i)/ave^2
                                                endfor
                                                pcount=pcount+1
                                     endif
                          endfor
                endfor
                
                g2t=g2t/pcount
                SetScale/P x 0,0.34,"", g2t
                duplicate/o g2t, $g2tname
                killwaves g2t
                print pcount
end

Function correlatefunc_1Dwave(data, dt)    // Take annular average intnesity I(k), then calculate g2t.
	       wave data
	       variable dt
	       
	       variable tp = DimSize(data, 0)
	
               make/o/n=(tp) g2t=0
	       variable i, j, add_td, add_t, add_t_td
	       string g2tname="g2t_"+nameofwave(data)
	       
	       for(i=0; i<tp; i+=1)
	                 add_td=0
	                 add_t=0
	                 add_t_td=0
		         for(j=0; j<tp-i; j+=1)
			                add_td=add_td+data[j]*data[j+i]
			                add_t=add_t+data[j]
			                add_t_td=add_t_td+data[j+i]
			endfor
                         g2t[i]=(tp-i)*add_td/add_t/add_t_td
               endfor
               
               SetScale/P x 0,dt,"", g2t 
               duplicate/o g2t, $g2tname
               killwaves g2t
end

Function correlatefunc_k3(data)    // Take annular average intnesity I(k), then calculate g2t.
	       wave data
	       
	       variable tp = DimSize(data, 0)
	
               make/o/n=(tp) g2t=0
	       variable i, j, add_td, ave=mean(data)
	       string g2tname="g2t_"+nameofwave(data)
	       
	       for(i=0; i<tp; i+=1)
	                 add_td=0
		         for(j=0; j<tp-i; j+=1)
			                add_td=add_td+data[j]*data[j+i]
			endfor
                         g2t[i]=add_td/(tp-i)/ave^2
               endfor
               
               SetScale/P x 0,0.34,"", g2t 
               duplicate/o g2t, $g2tname
               killwaves g2t
end

Function correlatefunc_test(data)    // under construction. test on a data[tp], average over t' to t'+t. Denominator <I(t')>*<I(t'+t)>
	       wave data
	       
	       variable tp = DimSize(data, 0)
	
               make/o/n=(tp) g2t=0, add_td=0
	       variable i, j
	       string g2tname="g2t_"+nameofwave(data)
	      
	       for(i=0; i<tp; i+=1)
		         for(j=0; j<(tp-i); j+=1)
			                add_td[i]=add_td[i]+data[j]*data[j+i]
			endfor
               endfor
               
               g2t[]=(tp-p)*add_td[p]/sum(data,0,tp-p-1)/sum(data, p, tp-1) 
               
               duplicate/o g2t, $g2tname
               killwaves g2t, add_td
end

Function correlatefunc_test3(data)    // test on a data[tp], average over t' to t'+t.
	       wave data
	       
	       variable tp = DimSize(data, 0)
	
               make/o/n=(tp) g2t=0
	       variable i, j, add_td, add_t, ave
	       string g2tname="g2t_"+nameofwave(data)
	       
	       for(i=0; i<tp; i+=1)
	                 add_td=0
	                 add_t=0
		         for(j=0; j<tp-i; j+=1)
			                add_td=add_td+data[j]*data[j+i]
			                add_t=add_t+data[j]+data[j+i]
			endfor
			ave=add_t/(tp-i)/2
                         g2t[i]=add_td/(tp-i)/ave^2
               endfor
                
               duplicate/o g2t, $g2tname
               killwaves g2t
end	       

Function correlatefunc_test2(data)    // test on a data[tp], average over all frames.
	       wave data
	       
	       variable tp = DimSize(data, 0)
	
               make/o/n=(tp) g2t=0
	       variable i, j, add_td, ave=mean(data)
	       string g2tname="g2t_"+nameofwave(data)
	       
	       for(i=0; i<tp; i+=1)
	                 add_td=0
		         for(j=0; j<tp-i; j+=1)
			                add_td=add_td+data[j]*data[j+i]
			endfor
                         g2t[i]=add_td/(tp-i)/ave^2
               endfor
                
               duplicate/o g2t, $g2tname
               killwaves g2t
end	       
	                             
Function relaxation_points(data, x, y)
               wave data
               variable x, y
               variable zp = DimSize(data, 2)
               string dataname=nameofwave(data)
               string g2tname="g2t_"+dataname+"_"+num2istr(x)+"_"+num2istr(y)
               string parameterName="parameter"+dataname+"_"+num2istr(x)+"_"+num2istr(y)
               string intensityname="intensity"+dataname+"_"+num2istr(x)+"_"+num2istr(y)
               
               make/o/n=(zp) g2t
               
               redimension/s data
               Imagetransform/beam={(x), (y)} getbeam data
               wave W_beam = $"W_beam"
               duplicate/o W_beam intensity
               
               Correlate/auto W_beam, W_beam
               g2t[] = W_beam[p+zp-1]
               Setscale/P x, 0, 1, "", g2t
               Display g2t
               rename g2t, $g2tname
               
               Make/o/D/N=2/O W_coef
               W_coef[0] = {100,1}
               Funcfit/NTHR=0/Q relaxationfunc W_coef, $g2tname/D
               rename W_coef, $parameterName
               rename intensity, $intensityname
               
               printf "tau, beta: %g, %g\r" W_coef[0], W_coef[1]
              
end

Function relaxation_2D1(data,xstart,xend,ystart,yend)
                wave data
                variable xstart,xend,ystart,yend
                variable xp=xend-xstart+1, yp=yend-ystart+1, zp = DimSize(data, 2)
                variable i, j
                string outputName=nameofwave(data)+"_"+num2istr(xstart)+"_"+num2istr(ystart)+"_"+num2istr(xp)+"by"+num2istr(yp)
                string tauname="tau_"+outputName
                string betaname="beta_"+outputName
                
                make/o/n=(xp, yp, zp) g2t
                make/o/n=(xp,yp) tau
                make/o/n=(xp,yp) beta_parameter
                
                for (i=0; i<xp; i+=1)
                        for (j=0; j<yp; j+=1)
                                Imagetransform/beam={(i+xstart), (j+ystart)} getbeam data
                                wave W_beam = $"W_beam"
                                Correlate/auto W_beam, W_beam
                                g2t[i][j][] = W_beam[r+zp-1]
                                Setscale/P x, 0, 1, "", g2t
                                
                                Make/o/D/N=2/O W_coef
                                W_coef[0] = {100,1}
                                Funcfit/NTHR=0/Q relaxationfunc W_coef, g2t[i][j][]
                                tau[i][j]=W_coef[0]
                                beta_parameter[i][j]=W_coef[1]
                         endfor
                 endfor
                 
                 rename tau, $tauname
                 rename beta_parameter, $betaname
end