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
output_name = "fsun_Tveg_diff"+year_begin +"_"+year_end
;title_name = "TV year "+year_begin+ " to year " + year_end +" (new_fsun-default)"
title_name = "Vegetation Temperature difference (New Fsun - Default)"

TV_control         = new((/96,144,nyear/),double)
TV_control_monthly = new((/96,144,12/),double)
TV_control_AVG     = new((/96,144/),double)
TV_compare            = new((/96,144,nyear/),double)
TV_compare_monthly    = new((/96,144,12/),double)
TV_compare_AVG        = new((/96,144/),double)
diff_AVG_TV        = new((/96,144/),double)
Areasum_control_TV = new((/96,144/),double)
Areasum_compare_TV    = new((/96,144/),double)

do year=year_begin,year_end,1
  n=year-year_begin

  do month=1,12

     fname1 = "/home/reichpb/eebutler/github/E3SM/run/global2_480_v1_f19_f19_ICB20TRCNPRDCTCBC/run/global2_480_v1_f19_f19_ICB20TRCNPRDCTCBC.clm2.h0."+ sprinti("%0.4d",year) + "-" + sprinti("%0.2d",month) +".nc"
     f1 = addfile(fname1,"r")
     TV_control_monthly(:,:,(month-1))=f1->TV(0,:,:)

     fname2 = "/home/reichpb/eebutler/github/E3SM/run/fsun_g_v1_f19_f19_ICB20TRCNPRDCTCBC/run/fsun_g_v1_f19_f19_ICB20TRCNPRDCTCBC.clm2.h0."+ sprinti("%0.4d",year) + "-" + sprinti("%0.2d",month) +".nc"
     f2 = addfile(fname2,"r")
     TV_compare_monthly(:,:,(month-1))=f2->TV(0,:,:)

  end do

  TV_control(:,:,n)=dim_avg_n(TV_control_monthly,2)
  TV_compare(:,:,n)=dim_avg_n(TV_compare_monthly,2)

end do

;TV_control=TV_control*3600*24*365   ;365 days per year gC/m2
TV_control@units     = "C"
;TV_compare=TV_compare*3600*24*365   ;365 days per year gC/m2
TV_compare@units     = "C"

TV_control_AVG = dim_avg_n(TV_control,2)  ; averged TV over the years
TV_compare_AVG    = dim_avg_n(TV_compare,2)  ; averged TV over the years
landmask = f1->landmask(:,:)
area     = f1->area(:,:)*10e6 ; from km^2 to m^2
landfrac = f1->landfrac(:,:); 
TV_control_AVG = mask(TV_control_AVG,landmask,1)
n1=num(.not.ismissing(TV_control_AVG))
print(n1)
TV_compare_AVG    = mask(TV_compare_AVG,landmask,1)

do lat=0,95
   do lon=0,143
       Areasum_control_TV(lat,lon) = TV_control_AVG(lat,lon)*area(lat,lon)*landfrac(lat,lon)
       Areasum_compare_TV(lat,lon)    = TV_compare_AVG(lat,lon)*area(lat,lon)*landfrac(lat,lon)
   end do
end do

Areasum_control_TV@_FillValue = -999
Areasum_compare_TV@_FillValue = -999
Total_control_TV = sum(Areasum_control_TV)/10e15;
Total_compare_TV    = sum(Areasum_compare_TV)/10e15;
print(Total_control_TV);
print(Total_compare_TV);
diff_AVG_TV(:,:) = f1->TV(0,:,:)   ;trick to get the attribution
diff_AVG_TV@units     = "[~S~o~N~C]"
diff_AVG_TV@long_name = "vegetation temperature"

diff_AVG_TV = TV_compare_AVG - TV_control_AVG

diff_AVG_TV@_FillValue           = -999.

opt = True
opt@fout = "avg_tveg.txt"
write_matrix(diff_AVG_TV,"144f9.2",opt)
 
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
    res@cnMinLevelValF      = -1.
    res@cnMaxLevelValF      = 1.
    res@cnLevelSpacingF     = 0.1              	; 15 contour levels        
    res@cnFillOn            = True            	; color fill plot
    res@cnLinesOn           = False
    res@cnLineLabelsOn      = False
    res@cnInfoLabelOn       = False
    res@gsnZonalMean        = False             	; put zonal on side
    res@gsnZonalMeanXMinF   = -0.1		; set minimum X-axis value for zonal mean plot  
    res@gsnZonalMeanXMaxF   = 0.1		; set maximum X-axis value for zonal mean plot  
    res@gsnAddCyclic          = True    	; data already has cyclic point

    plot = new(1,graphic)    

    diff_AVG_TV@_FillValue           = -999.
    plot = gsn_csm_contour_map_ce(wks,diff_AVG_TV,res)
    
    resPanel                  = True
    resPanel@txString         = title_name                 ; add common title
    resPanel@gsnMaximize      = True                       ; use full page
    resPanel@gsnPaperOrientation = "landscape"
    resPanel@gsnPanelLabelBar = True                       ; add color bar
    resPanel@lbLabelAutoStride = True
  
    gsn_panel(wks,plot,(/1,1/),resPanel)               ; now draw as one 
end
