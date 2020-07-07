;*************************************************
; clmrun.ncl
;************************************************
load "/$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "/$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "/$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
;************************************************
begin
;************************************************
; open file and read in data
;************************************************

year_begin = 1980
year_end   = 2009
nyear = year_end - year_begin + 1
output_name = "noSAI_upd_GPP_diff"+year_begin +"_"+year_end
;title_name = "GPP year "+year_begin+ " to year " + year_end +" (new_fsun-default)"
title_name = "GPP difference (No SAI - Corrected Nitrogen Scaler)"

GPP_control         = new((/96,144,nyear/),double)
GPP_control_monthly = new((/96,144,12/),double)
GPP_control_AVG     = new((/96,144/),double)
GPP_compare            = new((/96,144,nyear/),double)
GPP_compare_monthly    = new((/96,144,12/),double)
GPP_compare_AVG        = new((/96,144/),double)
diff_AVG_GPP        = new((/96,144/),double)
Areasum_control_GPP = new((/96,144/),double)
Areasum_compare_GPP    = new((/96,144/),double)

do year=year_begin,year_end,1
  n=year-year_begin

  do month=1,12
     fname1 = "/home/reichpb/eebutler/github/E3SM/run/fsun_g_v1_f19_f19_ICB20TRCNPRDCTCBC/run/fsun_g_v1_f19_f19_ICB20TRCNPRDCTCBC.clm2.h0."+ sprinti("%0.4d",year) + "-" + sprinti("%0.2d",month) +".nc"
     f1 = addfile(fname1,"r")
     GPP_control_monthly(:,:,(month-1))=f1->GPP(0,:,:)

     fname2 = "/home/reichpb/eebutler/github/E3SM/run/fsun_nosai_v1_f19_f19_ICB20TRCNPRDCTCBC/run/fsun_nosai_v1_f19_f19_ICB20TRCNPRDCTCBC.clm2.h0."+ sprinti("%0.4d",year) + "-" + sprinti("%0.2d",month) +".nc"
     f2 = addfile(fname2,"r")
     GPP_compare_monthly(:,:,(month-1))=f2->GPP(0,:,:)

  end do

  GPP_control(:,:,n)=dim_avg_n(GPP_control_monthly,2)
  GPP_compare(:,:,n)=dim_avg_n(GPP_compare_monthly,2)

end do

GPP_control=GPP_control*3600*24*365   ;365 days per year gC/m2
GPP_control@units     = "(gC/m2 per year)"
GPP_compare=GPP_compare*3600*24*365   ;365 days per year gC/m2
GPP_compare@units     = "(gC/m2 per year)"

GPP_control_AVG = dim_avg_n(GPP_control,2)  ; averged GPP over the years
GPP_compare_AVG    = dim_avg_n(GPP_compare,2)  ; averged GPP over the years
landmask = f1->landmask(:,:)
area     = f1->area(:,:)*10e6 ; from km^2 to m^2
landfrac = f1->landfrac(:,:); 
GPP_control_AVG = mask(GPP_control_AVG,landmask,1)
n1=num(.not.ismissing(GPP_control_AVG))
print(n1)
GPP_compare_AVG    = mask(GPP_compare_AVG,landmask,1)

do lat=0,95
   do lon=0,143
       Areasum_control_GPP(lat,lon) = GPP_control_AVG(lat,lon)*area(lat,lon)*landfrac(lat,lon)
       Areasum_compare_GPP(lat,lon)    = GPP_compare_AVG(lat,lon)*area(lat,lon)*landfrac(lat,lon)
   end do
end do

Areasum_control_GPP@_FillValue = -999
Areasum_compare_GPP@_FillValue = -999
Total_control_GPP = sum(Areasum_control_GPP)/10e15;
Total_compare_GPP    = sum(Areasum_compare_GPP)/10e15;
print(Total_control_GPP);
print(Total_compare_GPP);
diff_AVG_GPP(:,:) = f1->GPP(0,:,:)   ;trick to get the attribution
diff_AVG_GPP@units     = "[gC/m~S~2~N~yr]"
diff_AVG_GPP@long_name = "gross primary production"

diff_AVG_GPP = GPP_compare_AVG - GPP_control_AVG

diff_AVG_GPP@_FillValue           = -999.
opt = True
opt@fout = "avg_gpp_noSAI.txt"
write_matrix(diff_AVG_GPP,"144f9.2",opt)

opt = True
opt@fout = "land_frac.txt"
write_matrix(landfrac,"144f9.2",opt)

opt = True
opt@fout = "land_area.txt"
write_matrix(area,"144f9.2",opt)
 
;************************************************
; create 1 plot
;************************************************
;**************************************************
; plot
;**************************************************
    wks  = gsn_open_wks("eps",output_name)         ; send graphics to PNG file
    gsn_define_colormap(wks,"BlAqGrYeOrRe")        ; choose colormap   
    
    res                     = True
    res@gsnDraw             = False
    res@gsnFrame            = False
    res@gsnSpreadColors     = True
    res@lbLabelBarOn        = False            ; turn off individual label bars
    res@mpFillOn            = False            ; no need   
    res@cnLevelSelectionMode= "ManualLevels"   ; manual set levels
    res@cnMinLevelValF      = -300.
    res@cnMaxLevelValF      = 300.
    res@cnLevelSpacingF     = 30              	; 15 contour levels        
    res@cnFillOn            = True            	; color fill plot
    res@cnLinesOn           = False
    res@cnLineLabelsOn      = False
    res@cnInfoLabelOn       = False
    res@gsnZonalMean        = True             	; put zonal on side
    res@gsnZonalMeanXMinF   = -100.		; set minimum X-axis value for zonal mean plot  
    res@gsnZonalMeanXMaxF   = 300.		; set maximum X-axis value for zonal mean plot  
    res@gsnAddCyclic          = True    	; data already has cyclic point

    plot = new(1,graphic)    

    diff_AVG_GPP@_FillValue           = -999.
    plot = gsn_csm_contour_map_ce(wks,diff_AVG_GPP,res)
    
    resPanel                  = True
    resPanel@txString         = title_name                 ; add common title
    resPanel@gsnMaximize      = True                       ; use full page
    resPanel@gsnPaperOrientation = "landscape"
    resPanel@gsnPanelLabelBar = True                       ; add color bar
    resPanel@lbLabelAutoStride = True
  
    gsn_panel(wks,plot,(/1,1/),resPanel)               ; now draw as one 
end
