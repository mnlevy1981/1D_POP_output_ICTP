load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"

begin

  out_type = "X11" ; plot to screen (change to "pdf" to make plots)
  exercise = 3

  if (exercise.eq.3) then
  casename = "C.1D_1D.r8_cvmix"
    nlevs_in_plot = 40
    min_sst       = 14.8  ; deg C
    max_sst       = 15.   ; deg C
    max_diff      = 0.025 ; m^2/s
    diff_contour  = 0.001 ; m^2/s
  end if
  if (exercise.eq.4) then
    casename = "C.1D_1D.no_taux"
    nlevs_in_plot = 100
    min_sst       = 14.  ; deg C
    max_sst       = 15.  ; deg C
    max_diff      = 0.15 ; m^2/s
    diff_contour  = 0.01 ; m^2/s
  end if
  if (exercise.eq.5) then
    casename = "C.1D_1D.taux_and_windstress"
    nlevs_in_plot = 25
    min_sst       = 14.8  ; deg C
    max_sst       = 17.   ; deg C
    max_diff      = 0.024 ; m^2/s
    diff_contour  = 0.001 ; m^2/s
  end if
  max_depth = 1.d*nlevs_in_plot    ; units are meters
  filename = "~/scratch/"+casename+"/run/"+casename+".pop.h.0001-01-02-01200.nc"

  ; Suffix is used for file names when making plots
  suffix = ".ex"+exercise

  ; Read data from netcdf file
  f = addfile(filename,"r")
  time = f->time - 366.d
  time@long_name = "time (days)"
  z_iface = -f->z_w_bot(0:nlevs_in_plot-1)/100.d   ; Convert from cm -> m
  z_cntr  = -f->z_t(0:nlevs_in_plot-1)/100.d       ; Convert from cm -> m
  z_iface@long_name = "depth (m)"
  z_cntr@long_name  = "depth (m)"

  h  = -f->HBLT(:,0,0)/100.d    ; Convert from cm -> m
  temp = f->TEMP(:,:,0,0)
  sst = temp(:,0)
  u = transpose(f->UVEL(:,0:nlevs_in_plot-1,0,0))/100.d   ; Convert from cm/s -> m/s
  u@long_name = "Zonal Velocity"
  u@units = "m/s"
  v = transpose(f->VVEL(:,0:nlevs_in_plot-1,0,0))/100.d   ; Convert from cm/s -> m/s
  v@long_name = "Meridonal Velocity"
  v@units = "m/s"
  Kt = transpose(f->VDC_T(:,0:nlevs_in_plot-1,0,0)/10000.d)   ; Convert from cm^2/s -> m^2/s
  Kt@long_name = "Diffusivity Coefficient for Temperature"
  Kt@units = "m^2/s"
  Ku = transpose(f->VVC(:,0:nlevs_in_plot-1,0,0)/10000.d)     ; Convert from cm^2/s -> m^2/s
  Ku@long_name = "Diffusivity Coefficient for Momentum"
  Ku@units = "m^2/s"
  delete(f)

  ; Set up u, v, Kt, and Ku  to plot correctly
  u!0 = "z"
  u!1 = "t"
  u&t = time
  u&z = z_cntr
  v!0 = "z"
  v!1 = "t"
  v&t = time
  v&z = z_cntr
  Kt!0 = "z"
  Kt!1 = "t"
  Kt&t = time
  Kt&z = z_cntr
  Ku!0 = "z"
  Ku!1 = "t"
  Ku&t = time
  Ku&z = z_cntr


  ; Print some statistics to sreen
  print((/"Deepest boundary layer depth: "+min(h)/))
  print((/"SST varies between "+min(sst)+" and "+max(sst)/))
  print((/"u varies between "+min(u)+" and "+max(u)/))
  print((/"Max temp diffusivity: "+max(Kt)/))
  print((/"Max momentum diffusivity: "+max(Ku)/))
  print((/"Ratio of max diffusivities: "+max(Kt)/max(Ku)/))

  ; Plot settings that are the same for HBLT and SST plots
  pres = True
  pres@gsnMaximize = True
  pres@gsnDraw = False
  pres@gsnFrame = False
  pres@tiMainFuncCode = ":"
  pres@tiXAxisFuncCode = ":"
  pres@tiYAxisFuncCode = ":"
  pres@trXMinF = 0.
  pres@trXMaxF = 15.
  pres@tiXAxisString = "time (days)"

  ; Contour settings that are the same for u, v, kappa, and nu
  cres = True
  cres@gsnMaximize = True
  cres@gsnDraw = False
  cres@gsnFrame = False
  cres@cnFillOn = True
  cres@gsnSpreadColors = True
  cres@cnLinesOn = False
  cres@cnLevelSelectionMode = "ManualLevels"

  ; Plot boundary layer depth
  wks1 = gsn_open_wks(out_type, "OBL_Depth"+suffix)
  pres@trYMinF = -max_depth
  pres@trYMaxF = 0.
  pres@tiMainString = "Boundary Layer Depth"
  pres@tiYAxisString = "depth (meters)"
  plot = gsn_csm_xy(wks1, time, h, pres)
  draw(plot)
  frame(wks1)

  ; Plot SST
  wks2 = gsn_open_wks(out_type, "SST"+suffix)
  pres@trYMinF = min_sst
  pres@trYMaxF = max_sst
  pres@tiMainString = "Sea Surface Temperature"
  pres@tiYAxisString = "temperature (deg C)"
  plot = gsn_csm_xy(wks2, time, sst, pres)
  draw(plot)
  frame(wks2)

  ; Contour map for velocity fields
  wks3 = gsn_open_wks(out_type, "UVEL"+suffix)
  cres@cnMinLevelValF = -0.1
  cres@cnMaxLevelValF = 0.15
  plot = gsn_csm_contour(wks3, u,cres)
  draw(plot)
  frame(wks3)

  wks4 = gsn_open_wks(out_type, "VVEL"+suffix)
  cres@cnMinLevelValF = -0.15
  cres@cnMaxLevelValF = 0.1
  plot = gsn_csm_contour(wks4, v,cres)
  draw(plot)
  frame(wks4)

  ; Contour map for diffusivity fields
  wks5 = gsn_open_wks(out_type, "Kt"+suffix)
  cres@cnMinLevelValF = 0.
  cres@cnMaxLevelValF = max_diff
  plot = gsn_csm_contour(wks5, Kt,cres)
  draw(plot)
  frame(wks5)

  wks6 = gsn_open_wks(out_type, "Ku"+suffix)
  cres@cnMinLevelValF = 0.
  cres@cnMaxLevelValF = max_diff
  plot = gsn_csm_contour(wks6, Ku,cres)
  draw(plot)
  frame(wks6)

  ; Close all plots
  delete(wks1)
  delete(wks2)
  delete(wks3)
  delete(wks4)
  delete(wks5)
  delete(wks6)


end
