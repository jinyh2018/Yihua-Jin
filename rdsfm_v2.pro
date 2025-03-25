;---------------------------------------------------------------------------------
;              Residual Distribution based Spatio-temporal data Fusion Method
;        
;                   Developed by Yihua Jin, email: jinyh2018@ybu.edu.cn
;  This program was created by referencing and  modifying the mad_run5.pro provided by Mort Canty (2013), 
;  and FSDAF.pro provided by Xiaolin Zhu (2016). 
;  Thank you to Mort Canty and Xiaolin Zhu for providing the IDL code for the program.               
;              
;;---------------------------------------------------------------------------------

Function mad_iter, fid1, fid2, dims1, dims2, pos1, pos2, m_fid, niter=niter, rho=rho, $
  sigma=sigma, A=A, B=B, means1=means1, means2=means2, $
  verbose=verbose, lam=lam

  COMPILE_OPT idl2

  if n_elements(niter) eq 0 then niter = 50
  if n_elements(verbose) eq 0 then verbose = 0
  if n_elements(lam) eq 0 then lam = 0.0

  envi_file_query, fid1, interleave=interleave1
  envi_file_query, fid2, interleave=interleave2

  num_cols = dims1[2]-dims1[1]+1
  num_rows = dims1[4]-dims1[3]+1
  num_pixels = (num_cols*num_rows)
  num_bands = n_elements(pos1)

  if num_bands eq 1 then interleave1 = (interleave2 = 2)

  ; column vector of 1s
  ones = fltarr(num_cols,1)+1.0

  if m_fid ne -1 then mask = envi_get_data(fid=m_fid,dims=dims1,pos=0) $
  else mask = bytarr(num_cols,num_rows)+1B
  if niter gt 0 then rhos=fltarr(num_bands,niter)-1

  if niter gt 1 then begin
    window, 12, xsize=600, ysize=400, title='IR-MAD'
    wset,12
    plot, indgen(niter), rhos[0,*], psym=4, xrange = [0,niter+1], yrange = [0,1.1], $
      xtitle='iteration', ytitle='Correlations',color=0,background='FFFFFF'XL
  endif

  ; get tile ids
  tile_id1 = envi_init_tile(fid1,pos1,num_tiles=num_tiles, $
    interleave=interleave1,xs=dims1[1],xe=dims1[2],ys=dims1[3],ye=dims1[4])
  tile_id2 = envi_init_tile(fid2,pos2, $
    interleave=interleave2,xs=dims2[1],xe=dims2[2],ys=dims2[3],ye=dims2[4])

  sigMADs = fltarr(num_bands,num_cols)

  ; length penalization
  Omega_L = Identity(num_bands)

  ; begin iteration
  iter=0L
  interrupting = 0
  delta = 1.0
  old_rho = fltarr(num_bands)
  repeat begin
    progressbar = Obj_New('cgprogressbar', Text=' ',$
      title='IR-MAD',xsize=300,ysize=20,/cancel)
    progressbar->start
    cpm = obj_new('CPM',2*num_bands)
    txt = 'Iter '+strtrim(iter,2)
    ; spectral tiling
    for tile_index=0L,num_tiles-1 do begin
      if progressbar->CheckCancel() then begin
        if (iter gt 0) and not interrupting then begin
          iter = niter
          interrupting = 1
          txt = 'Interrupting...'
        end else begin
          print,'Calculation aborted'
          obj_destroy, cpm
          envi_tile_done, tile_id1
          envi_tile_done, tile_id2
          wdelete,12
          message, 'Calculation aborted'
        endelse
      endif
      if tile_index mod 10 eq 0 then begin
        pct = (tile_index)*100/(num_tiles)
        progressbar->Update,fix(pct)
      endif
      tile1 = envi_get_tile(tile_id1,tile_index)
      tile2 = envi_get_tile(tile_id2,tile_index)
      if interleave1 eq 1 then tile1 = transpose(tile1)
      if interleave2 eq 1 then tile2 = transpose(tile2)
      if iter gt 0 then begin
        mads = ((tile1-means1)##A-(tile2-means2)##B)
        chi_sqr = total((mads/sigMADs)^2,1)
        weights=1-chisqr_pdf(chi_sqr,num_bands)
      end else weights = ones
      ;       only sample pixels under the mask
      indices = where(mask[*,tile_index],count)
      if count gt 0 then $
        cpm->update,[tile1[*,indices],tile2[*,indices]],weights=weights[indices]
    endfor
    progressbar->destroy

    ; canonical correlation
    SS = cpm->Covariance()
    means = cpm->Means()
    Obj_Destroy, cpm

    S11 = SS[0:num_bands-1,0:num_bands-1]
    S11 = (1-lam)*S11 + lam*omega_L
    S22 = SS[num_bands:*,num_bands:*]
    S22 = (1-lam)*S22 + lam*omega_L
    S12 = SS[num_bands:*,0:num_bands-1]
    S21 = SS[0:num_bands-1,num_bands:*]

    C1 = S12 ## invert(S22,/double) ## S21
    B1 = S11
    C2 = S21 ## invert(S11,/double) ## S12
    B2 = S22

    if num_bands gt 1 then begin
      gen_eigenproblem, C1, B1, A, mu2
      gen_eigenproblem, C2, B2, B, mu2
    end else begin
      mu2 = [C1/B1]
      A = [1/sqrt(B1)]
      B = [1/sqrt(B2)]
    endelse

    mu = sqrt(mu2)
    a2=diag_matrix(transpose(A)##A)
    b2=diag_matrix(transpose(B)##B)
    sigma = sqrt( (2-lam*(a2+b2))/(1-lam)-2*mu )
    rho=mu*(1-lam)/sqrt( (1-lam*a2)*(1-lam*b2) )

    sigMads = ones##sigma
    means1  = ones##means[0:num_bands-1]
    means2  = ones##means[num_bands:*]

    ; stopping criterion
    delta = max(abs(rho-old_rho))
    old_rho = rho

    ; ensure sum of positive correlations between X and U is positive
    ; their covariance matrix is S11##A
    invsig_x = diag_matrix(1/sqrt(diag_matrix(S11)))
    sum = total(invsig_x##S11##A,2)
    A = A##diag_matrix(sum/abs(sum))

    ; ensure positive correlation between each pair of canonical variates
    cov = diag_matrix(transpose(A)##S12##B)
    B = B##diag_matrix(cov/abs(cov))

    if iter gt 0 and iter eq niter then goto, done

    ;    print sigmas and plot canonical correlations
    if verbose then begin
      print, 'delta = '+strtrim(delta,2)
      print, 'iteration '+strtrim(iter,2)
      print, reverse(sigMADs[*,0]), format='("Sigma MADs",'+strtrim(num_bands,2)+'f10.5)'
    endif
    if (niter gt 1) and (iter lt niter) then begin
      rhos[*,iter] = rho
      wset,12
      plot, indgen(iter+1)+1, rhos[0,0:iter], psym=-4, xrange=[0,niter+1],  yrange = [0,1.1], $
        xtitle='iteration', ytitle='Correlations',color=0,background='FFFFFF'XL
      for i=1,num_bands-1 do $
        oplot, indgen(iter+1)+1, rhos[i,0:iter], psym=-4, color=0
    endif
    done:
    iter=iter+1
  endrep until (iter gt niter) or (delta lt 0.001)
  envi_tile_done, tile_id1
  envi_tile_done, tile_id2

  ; successfully finished, so
  return, 0

End

Pro GetData,ImgData = ImgData,ns = ns,nl = nl,nb = nb,Data_Type = Data_Type,$
  FileName = FileName,Map_info = map_Info, Fid = Fid
  Filter = ['all file;*.*']
  Envi_Open_File,FileName,R_Fid = Fid
  Envi_File_Query,Fid,ns = ns,nl = nl,nb = nb,Data_Type = Data_Type
  map_info = envi_get_map_info(fid=Fid)
  dims = [-1,0,ns - 1 ,0,nl - 1]
  case Data_Type Of
    1:ImgData = BytArr(ns,nl,nb)    ;  BYTE  Byte
    2:ImgData = IntArr(ns,nl,nb)    ;  INT  Integer
    3:ImgData = LonArr(ns,nl,nb)    ;  LONG  Longword integer
    4:ImgData = FltArr(ns,nl,nb)    ;  FLOAT  Floating point
    5:ImgData = DblArr(ns,nl,nb)    ;  DOUBLE  Double-precision floating
    6:ImgData = COMPLEXARR(ns,nl,nb); complex, single-precision, floating-point
    9:ImgData = DCOMPLEXARR(ns,nl,nb);complex, double-precision, floating-point
    12:ImgData = UINTARR(ns,nl,nb)   ; unsigned integer vector or array
    13:ImgData = ULONARR(ns,nl,nb)   ;  unsigned longword integer vector or array
    14:ImgData = LON64ARR(ns,nl,nb)   ;a 64-bit integer vector or array
    15:ImgData = ULON64ARR(ns,nl,nb)   ;an unsigned 64-bit integer vector or array
  EndCase
  For i = 0,nb-1 Do Begin
    Dt = Envi_Get_Data(Fid = Fid,dims = dims,pos=i)
    ImgData[*,*,i] = Dt[*,*]
  EndFor
End

;function of multiple linear regression without interception
Function P_OLS, ind_v,dep_v,k,low,high    ;format: ind_v:k*n,dep_v:1*n

  common V_PUB, matrix
  common V_PUB2, y

  nb=n_elements(dep_v)
  x=reform(dep_v)
  varend=fltarr(k)
  yend=fltarr(k,nb)
  for ii=0,k-1 do begin
    varend[ii]=sqrt(total((ind_v[ii,*]-mean(ind_v[ii,*]))^2))
    yend[ii,*]=ind_v[ii,*]
  endfor
  var=sqrt(total((dep_v-mean(dep_v))^2))

  y=dep_v
  matrix=yend
  y=transpose(double(y))

  glow=fltarr(k)+low
  ghigh=fltarr(k)+high
  gintial=fltarr(k)+1.0/float(k)
  gbnd   = [glow, ghigh]
  Lbnd =[0,100]
  nobj   = 1
  g      = gintial
  Lcomp = 'HMBL11'
  nobj=0
  CONSTRAINED_MIN, g, gbnd, Lbnd, nobj, Lcomp, inform
  L =total((matrix ## g- y)^2)
  return,g
END

FUNCTION HMBL11, g
  common V_PUB
  common V_PUB2
  L=total((matrix ## g-y)^2)
  RETURN, L
END
;///////////////////////////////////////////////////////////////////////////////////////////
;...................................Main Program............................................
;///////////////////////////////////////////////////////////////////////////////////////////

Pro RDSFM_v2

  t0=systime(1)                  ;the initial time of program running

  ;please set the following parameters
  ;----------------------------------------------------------------------
  iter=25              ; set the iteration of IR-MAD
  penalization=0                ; set the penalization, the default value is 0
  w=20                 ;set the half window size, if 25, the window size is 25*2+1=51
  num_similar_pixel=20         ;set number of similar pixels
  min_class=4.0                ;set the estimated minimum and maximum number of classes
  max_class=10.0
  num_pure=100                 ;number of most purest coarse pixels in each class selected fro change calculation
  DN_min=0.0                   ;set the range of DN value of the image,If byte, 0 and 255
  DN_max=10000.0
  scale_factor=16              ;set the scale factor, it is integer=coarse resolution/fine resolution, e.g., 480/30=16
  block_size=20               ;set the size of block, e.g., 20 means 20*20 coarse pixels, if process whole ETM scene, set 30~50
  background=0                 ;set the value of background pixels. 0 means that pixels will be considered as background if one of its bands= 0
  background_band=1            ;which band with value = background indicating background pixels. Sometimes, background pixels have different values in different bands
  temp_file='D:\temppp'          ;set path of a folder to store temporary files
  file_mkdir,temp_file
  ;------------------------------------------------------------------------
e5 = envi(/current)
   if e5 eq !null then $
      message, 'This extension requires an interactive ENVI session   

;open the fine image of the first pair
envi_select, title='Choose fine image image of the first pair; F1 (and mask if desired)', fid=fid1, dims=dims1, pos=pos1, /mask, m_fid=m_fid
if (fid1 eq -1) then  begin
   print,'Cancelled'
   return
endif
envi_file_query, fid1, fname=fname1, xstart=xstart1, ystart=ystart1, interleave=interleave1, nb=nb, ns=ns, nl=nl, sname=sname
print,'first image: ',fname1

if (interleave1 eq 0) and (nb gt 1) then begin
    answer = dialog_message(file_basename(fname1)+' will be coverted to BIP. Continue?',/question)
    if answer eq 'No' then begin
       print,'Cancelled'
       return
    endif
    if strmid(sname,0,3) eq '[Me' then begin
       void = dialog_message('Conversion not possible, data in memory',/error)
       print,'Conversion to BIP not possible, cancelled
       return
    endif
    dims =  [-1L,0,ns-1,0,nl-1]
    pos = lindgen(nb)
    ENVI_DOIT, 'CONVERT_INPLACE_DOIT', fid=fid1, o_interleave=2, dims=dims, pos=pos, r_fid=r_fid
    fid1 = r_fid
    interleave1 = 2
endif

if nb eq 1 then interleave1 = 2

map_info=envi_get_map_info(fid=fid1)    
envi_convert_file_coordinates, fid1, dims1[1], dims1[3], ee, n, /to_map   ; 이게 왜 필요할가? 
map_info.mc=[0D, 0D, ee, n]     ; 이건 또 뭘가?

;open the coarse image of the first pair
envi_select, title='Choose coarse image of the first pair(C1)', fid=fid3, dims=dims3,pos=pos3
envi_file_query, fid3, fname=fname3

;open the Coarse image of the prediction time
envi_select, title='Choose coarse image of the prediction time(C2)', fid=fid2, dims=dims2,pos=pos2
;envi_file_query, fid2, fname=fname3

if (fid2 eq -1) then begin
   print,'Cancelled'
   return
endif

envi_file_query, fid2, fname=fname2, xstart=xstart2, ystart=ystart2, interleave=interleave2, $
                 nb=nb, ns=ns, nl=nl, sname=sname
print,'second image: ',fname2

if (n_elements(pos1) ne n_elements(pos2)) or $
  ((dims1[2]-dims1[1]) ne (dims2[2]-dims2[1])) or $
  ((dims1[4]-dims1[3]) ne (dims2[4]-dims2[3])) then begin
  print, 'Spectral/spatial subset sizes are different. Aborting.'
  Message, 'Spectral/spatial subsets different'
endif

if (interleave2 eq 0) and (nb gt 1) then begin
    answer = dialog_message(file_basename(fname2)+' will be coverted to BIP. Continue?',/question)
    if answer eq 'No' then begin
       print,'cancelled'
       return
    endif
    if strmid(sname,0,3) eq '[Me' then begin
       void = dialog_message('Conversion not possible, data in memory',/error)
       print,'Conversion to BIP not possible, cancelled
       return
    endif
    dims =  [-1L,0,ns-1,0,nl-1]
    pos = lindgen(nb)
    ENVI_DOIT, 'CONVERT_INPLACE_DOIT', fid=fid2, o_interleave=2, dims=dims, pos=pos, r_fid=r_fid
    fid2 = r_fid
    interleave2 = 2
endif

if nb eq 1 then interleave2 = 2

num_cols = dims1[2]-dims1[1]+1
num_rows = dims1[4]-dims1[3]+1
num_pixels = (num_cols*num_rows)
num_bands = n_elements(pos1)

chi_sqr=fltarr(num_cols)

niter=iter  ; number of iteration setting
lam=penalization   ; number of penalization

; MAD output
base=widget_auto_base(title='MAD Output')
sb=widget_base(base, /row, /frame)
wp=widget_outfm(sb, uvalue='outf', /auto)
result_mad=auto_wid_mng(base)
if not result_mad.accept then $
  print, 'MAD output cancelled'
  
; CV output
base = widget_auto_base(title='CV Output (root name)')
sb = widget_base(base, /row, /frame)
wp = widget_outfm(sb, uvalue='outf', /auto)
result_cv = auto_wid_mng(base)
if not result_cv.accept then print, 'CV output cancelled'

; Stats output
outfile_stats=dialog_pickfile(filter='*.mad',default_extension='mad',/write,/overwrite_prompt,title='Save MAD stats to disk')

; iterated MAD ***********************************************************************
if mad_iter(fid1, fid2, dims1, dims2, pos1, pos2, m_fid, niter=niter, rho=rho, $
       sigma=sigma, A=A, B=B, lam=lam, means1=means1, means2=means2,/verbose) eq -1 $
       then Message, 'mad_iter was cancelled or failed'   
; ************************************************************************************
print, 'Canonical correlations:'
print, reverse(rho)

vMs = reverse(sigma^2)
sigMads=(1+fltarr(1,num_cols))##reverse(sigma)

; output to file or memory
txt1=''
txt=''
if result_mad.accept then begin
  if result_mad.outf.in_memory then begin
    MAD_array = fltarr(num_cols,num_rows,num_bands+1)
    txt = 'MADs -> memory,'
  end else begin
    ;envi_write_envi_file, result_mad, out_name=txt1
    openw, unit_mad, result_mad.outf.name, /get_lun
    txt1 = 'MADs -> file,'
  endelse
endif
if result_cv.accept then begin
  if result_cv.outf.in_memory then begin
    CV1_array = fltarr(num_cols,num_rows,num_bands)
    CV2_array = fltarr(num_cols,num_rows,num_bands)
    txt = txt+' CVs -> memory'
  end else begin
    openw, unit_cv1, result_cv.outf.name+'_1', /get_lun
    openw, unit_cv2, result_cv.outf.name+'_2', /get_lun
    txt = txt+' CVs -> file'
  endelse
endif
txt=txt+' ...'

tile_id1 = envi_init_tile(fid1,pos1,num_tiles=num_tiles, interleave=interleave1,xs=dims1[1],xe=dims1[2],ys=dims1[3],ye=dims1[4])
tile_id2 = envi_init_tile(fid2,pos2, interleave=interleave2,xs=dims2[1],xe=dims2[2],ys=dims2[3],ye=dims2[4])
progressbar = Obj_New('cgprogressbar', Text='0', title=txt,xsize=250,ysize=20)
progressbar->start
abort=0
for tile_index=0L,num_tiles-1 do begin
  tile1 = envi_get_tile(tile_id1,tile_index)
  tile2 = envi_get_tile(tile_id2,tile_index)
  if interleave1 eq 1 then tile1 = transpose(tile1)
  if interleave2 eq 1 then tile2 = transpose(tile2)
  CV1s = (tile1-means1)##A
  CV2s = (tile2-means2)##B
  MADs = reverse(CV1s - CV2s,1)
  chi_sqr = float(total((MADs/sigMADs)^2,1))
  if result_mad.accept then begin
    if result_mad.outf.in_memory then begin
      for i=0,num_bands-1 do MAD_array[*,tile_index,i] = float(MADs[i,*])
      MAD_array[*,tile_index,num_bands] = chi_sqr
    end else begin
      writeu,unit_mad,[[float(MADs)],float(transpose(chi_sqr))]
    endelse
  endif
  if result_cv.accept then begin
    if result_cv.outf.in_memory then begin
      for i=0,num_bands-1 do begin
        CV1_array[*,tile_index,i] = float(CV1s[i,*])
        CV2_array[*,tile_index,i] = float(CV2s[i,*])
      endfor
    end else begin
      writeu,unit_cv1,float(CV1s)
      writeu,unit_cv2,float(CV2s)
    endelse
  endif
  if progressbar->CheckCancel() then begin
    print,'output aborted'
    abort=1
    tile_index=num_tiles
  endif
  pct=tile_index*100/num_tiles
  progressbar->Update,pct
endfor
progressbar->destroy
envi_tile_done, tile_id1
envi_tile_done, tile_id2
if result_mad.accept then if not result_mad.outf.in_memory then free_lun, unit_mad
if result_cv.accept then if not result_cv.outf.in_memory then begin
  free_lun, unit_cv1
  free_lun, unit_cv2
endif
if abort then return

if result_mad.accept then begin
  bnames = strarr(num_bands+1)
  bnames[0:num_bands-1] = 'MAD('+file_basename(fname1)+':'+file_basename(fname2)+') '+strtrim(lindgen(num_bands)+1,2)
  ;bnames[num_bands] = 'CHI_SQR('+file_basename(fname1)+':'+file_basename(fname2)+')'
  if result_mad.outf.in_memory then begin
    envi_enter_data, MAD_array, $
      r_fid=r_fid, $
      map_info=map_info, $
      bnames=bnames, $
      xstart=xstart1+dims1[1], ystart=ystart1+dims1[3], $
      descrip='MAD: t1= '+file_basename(fname1)+' t2= '+file_basename(fname2)
    print, 'MADs written to memory'
    envi_assign_header_value, fid=r_fid, keyword='varMADs',value=vMs
  end else begin
    envi_setup_head ,fname=result_mad.outf.name, ns=num_cols, $
      r_fid=r_fid, $
      nl=num_rows, nb=num_bands+1, $
      data_type=4, interleave=2, /write, /open, $
      bnames=bnames, $
      map_info=map_info, $
      xstart=xstart1+dims1[1], ystart=ystart1+dims1[3], $
      descrip='MAD: t1= '+fname1+' t2= '+fname2
    print, 'File created ', result_mad.outf.name
    envi_assign_header_value, fid=r_fid, keyword='varMADs',value=vMs
    envi_write_file_header,r_fid
  endelse
endif

envi_assign_header_value, fid=r_fid, keyword='varMADs', $
  value=2*(1-rho)

if result_cv.accept then begin
  bnames1 = strarr(num_bands)
  bnames2 = strarr(num_bands)
  bnames1= 'CV1('+file_basename(fname1)+':'+file_basename(fname2)+') '+strtrim(lindgen(num_bands)+1,2)
  bnames2= 'CV2('+file_basename(fname1)+':'+file_basename(fname2)+') '+strtrim(lindgen(num_bands)+1,2)
  if result_cv.outf.in_memory then begin
    envi_enter_data, CV1_array, $
      map_info=map_info, $
      bnames=bnames1, $
      xstart=xstart1+dims1[1], ystart=ystart1+dims1[3], $
      descrip='CV1: t1= '+file_basename(fname1)+' t2= '+file_basename(fname2)
    envi_enter_data, CV2_array, $
      map_info=map_info, $
      bnames=bnames2, $
      xstart=xstart1+dims1[1], ystart=ystart1+dims1[3], $
      descrip='CV2: t1= '+file_basename(fname1)+' t2= '+file_basename(fname2)
    print, 'CVs written to memory'
  end else begin
    envi_setup_head ,fname=result_cv.outf.name+'_1', ns=num_cols, $
      nl=num_rows, nb=num_bands, $
      data_type=4, interleave=2, /write, /open, $
      bnames=bnames1, $
      map_info=map_info, $
      xstart=xstart1+dims1[1], ystart=ystart1+dims1[3], $
      descrip='CV1: t1= '+fname1+' t2= '+fname2
    envi_setup_head ,fname=result_cv.outf.name+'_2', ns=num_cols, $
      nl=num_rows, nb=num_bands, $
      data_type=4, interleave=2, /write, /open, $
      bnames=bnames2, $
      map_info=map_info, $
      xstart=xstart1+dims1[1], ystart=ystart1+dims1[3], $
      descrip='CV2: t1= '+fname1+' t2= '+fname2
    print, 'File created ', result_cv.outf.name+'_1'
    print, 'File created ', result_cv.outf.name+'_2'
  endelse
endif

; output statistics

if (outfile_stats eq '') then print,'output of MAD stats was cancelled' $
else begin
  openw,lun,outfile_stats,/get_lun
  printf,lun,'; MAD statistics for '+'t1= '+file_basename(fname1)+' t2= '+file_basename(fname2)
  printf,lun,'; '+systime(0)
  printf,lun,strtrim(num_bands,2)
  printf,lun, [[means1[*,0]],[means2[*,0]],[A],[B],[rho],[sigma]],format='('+strtrim(num_bands,2)+'E26.16)'
  free_lun,lun
  print, 'MAD statistics written to '+outfile_stats
endelse

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;complete the MAD calculation and spatiotemporal data fusion start
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


map_info = envi_get_map_info(fid = fid1)
patch_long=block_size*scale_factor
orig_ns=ns
orig_nl=nl
n_ns=ceil(float(orig_ns)/patch_long)
n_nl=ceil(float(orig_nl)/patch_long)

ind_patch1=intarr(4,n_ns*n_nl)           ;divide the whole scene into 1000*1000 block
ind_patch=intarr(4,n_ns*n_nl)
location=intarr(4,n_ns*n_nl)

for i_ns=0,n_ns-1,1 do begin
  for i_nl=0,n_nl-1,1 do begin
    ind_patch1[0,n_ns*i_nl+i_ns]=i_ns*patch_long
    ind_patch[0,n_ns*i_nl+i_ns]=max([0,ind_patch1[0,n_ns*i_nl+i_ns]-scale_factor])
    location[0,n_ns*i_nl+i_ns]=ind_patch1[0,n_ns*i_nl+i_ns]-ind_patch[0,n_ns*i_nl+i_ns]

    ind_patch1[1,n_ns*i_nl+i_ns]=min([ns-1,(i_ns+1)*patch_long-1])
    ind_patch[1,n_ns*i_nl+i_ns]=min([ns-1,ind_patch1[1,n_ns*i_nl+i_ns]+scale_factor])
    location[1,n_ns*i_nl+i_ns]=ind_patch1[1,n_ns*i_nl+i_ns]-ind_patch1[0,n_ns*i_nl+i_ns]+location[0,n_ns*i_nl+i_ns]

    ind_patch1[2,n_ns*i_nl+i_ns]=i_nl*patch_long
    ind_patch[2,n_ns*i_nl+i_ns]=max([0,ind_patch1[2,n_ns*i_nl+i_ns]-scale_factor])
    location[2,n_ns*i_nl+i_ns]=ind_patch1[2,n_ns*i_nl+i_ns]-ind_patch[2,n_ns*i_nl+i_ns]

    ind_patch1[3,n_ns*i_nl+i_ns]=min([nl-1,(i_nl+1)*patch_long-1])
    ind_patch[3,n_ns*i_nl+i_ns]=min([nl-1,ind_patch1[3,n_ns*i_nl+i_ns]+scale_factor])
    location[3,n_ns*i_nl+i_ns]=ind_patch1[3,n_ns*i_nl+i_ns]-ind_patch1[2,n_ns*i_nl+i_ns]+location[2,n_ns*i_nl+i_ns]
  endfor
endfor

tempoutname=temp_file+'\temp_F1'

pos=indgen(nb)
for isub=0,n_ns*n_nl-1,1 do begin
  dims=[-1,ind_patch[0,isub],ind_patch[1,isub],ind_patch[2,isub],ind_patch[3,isub]]
  envi_doit, 'resize_doit', fid=fid1, pos=pos, dims=dims, interp=0, rfact=[1,1], $
    out_name=tempoutname+strtrim(isub+1,1), r_fid=r_fid1
  envi_file_mng, id=r_fid1, /remove
endfor

;resize coarse1 file
tempoutname=temp_file+'\temp_C1'
pos=indgen(nb)
for isub=0,n_ns*n_nl-1,1 do begin
  dims=[-1,ind_patch[0,isub],ind_patch[1,isub],ind_patch[2,isub],ind_patch[3,isub]]
  envi_doit, 'resize_doit', fid=fid3, pos=pos, dims=dims, interp=0, rfact=[1,1], $
    out_name=tempoutname+strtrim(isub+1,1), r_fid=r_fid3
  envi_file_mng, id=r_fid3, /remove
endfor
envi_file_mng, id=fid3, /remove

;; resize coarse2 file
tempoutname=temp_file+'\temp_C2'
pos=indgen(nb)
for isub=0,n_ns*n_nl-1,1 do begin
  dims=[-1,ind_patch[0,isub],ind_patch[1,isub],ind_patch[2,isub],ind_patch[3,isub]]
  envi_doit, 'resize_doit', fid=fid2, pos=pos, dims=dims, interp=0, rfact=[1,1], $
    out_name=tempoutname+strtrim(isub+1,1), r_fid=r_fid2
  envi_file_mng, id=r_fid2, /remove
endfor
envi_file_mng, id=fid2, /remove

; resize mad file
tempoutname=temp_file+'\temp_mad'
pos=indgen(nb)
for isub=0, n_ns*n_nl-1, 1 do begin
  dims=[-1, ind_patch[0,isub], ind_patch[1,isub], ind_patch[2,isub], ind_patch[3,isub]]
  envi_doit, 'resize_doit', fid=r_fid, pos=pos, dims=dims, interp=0, rfact=[1,1],$
    out_name=tempoutname+strtrim(isub+1,1), r_fid=r_fid4
  envi_file_mng, id=r_fid4, /remove
endfor
;envi_file_mng, id=r_fid, /remove

GetData, ImgData=fine1_whole, Filename=fname1, fid=fid5
GetData, ImgData=coarse1_whole, Filename=fname3, fid=fid6
GetData, ImgData=coarse2_whole, Filename=fname2, fid=fid7

change_mean=fltarr(ns, nl, nb)
for ib=0, nb-1, 1 do begin
  mean=mean(coarse2_whole[*,*,ib]-fine1_whole[*,*,ib])
  change_mean[*,*,ib]=make_array(ns, nl, value=mean)
endfor

;replace background pixels by mean of non-background to avoid its effects on classification
background_whole=bytarr(ns,nl)
ind_back=where(fine1_whole[*,*,background_band-1] eq background, num_back)
if (num_back gt 0) then begin
  background_whole[ind_back]=1
  for iband=0, nb-1, 1 do begin
    temp=fine1_whole[*,*,iband]
    temp[ind_back]=mean(temp[where(background_whole eq 0)])
    fine1_whole[*,*,iband]=temp
  endfor
endif
tempoutname11=temp_file+'\fine1_nobackground'
Envi_Write_Envi_File,fine1_whole,Out_Name = tempoutname11
ind_back=0;clear this variable
temp=0;clear this variable
fine1_whole=0 ;clear this variable
background_whole=0; clear this variable

envi_open_file,tempoutname11,r_fid=fid00

;step1: get spectral classes from fine resolution image at t1 by isodata
;parameter of isodata
CHANGE_THRESH = .05
NUM_CLASSES = max_class
ITERATIONS = 20
ISO_MERGE_DIST = 0.05*DN_max
ISO_MERGE_PAIRS = 2
ISO_MIN_PIXELS = 200
ISO_SPLIT_SMULT = 1
ISO_SPLIT_STD = 0.05*DN_max
MIN_CLASSES = min_class
out_bname = 'IsoData'
out_name=temp_file+'\class_ISODATA'
ENVI_File_Query, fid00, DIMS=dims, NB=nb1
ENVI_DOIT, 'class_doit', fid=fid00, pos=indgen(nb1), dims=dims, $
  out_bname=out_bname, out_name=out_name, method=4, $
  r_fid=r_fid, $
  NUM_CLASSES = NUM_CLASSES, $
  ITERATIONS = ITERATIONS, $
  CHANGE_THRESH = CHANGE_THRESH, $
  ISO_MERGE_DIST = ISO_MERGE_DIST, $
  ISO_MERGE_PAIRS = ISO_MERGE_PAIRS, $
  ISO_MIN_PIXELS = ISO_MIN_PIXELS, $
  ISO_SPLIT_SMULT = ISO_SPLIT_SMULT, $
  ISO_SPLIT_STD = ISO_SPLIT_STD, $
  MIN_CLASSES = MIN_CLASSES

envi_open_file,out_name,r_fid=fid
envi_file_query, fid, dims=dims

tempoutname=temp_file+'\class'
pos=indgen(1)
for isub=0,n_ns*n_nl-1,1 do begin
  dims=[-1,ind_patch[0,isub],ind_patch[1,isub],ind_patch[2,isub],ind_patch[3,isub]]
  envi_doit, 'resize_doit', fid=fid, pos=pos, dims=dims, interp=0, rfact=[1,1], $
    out_name=tempoutname+strtrim(isub+1,1), r_fid=r_fid1
  envi_file_mng, id=r_fid1, /remove
endfor

envi_file_mng, id=fid, /remove
envi_file_mng, id=fid00, /remove
envi_file_mng, id=fid5, /remove
envi_file_mng, id=fid6, /remove
envi_file_mng, id=fid7, /remove

;------------------------------------------------------------------
;process  each block
;-------------------------------------------------------------------

print,'there are total',n_ns*n_nl,' blocks'

for isub=0,n_ns*n_nl-1,1 do begin

  ;open each block image

  FileName = temp_file+'\temp_F1'
  GetData,ImgData=fine1,ns = ns,nl = nl,nb = nb,Data_Type = Data_Type,FileName = FileName+strtrim(isub+1,1),Fid = Fid11
  fine1=float(fine1)

  FileName = temp_file+'\temp_C1'
  GetData,ImgData=coarse1,FileName = FileName+strtrim(isub+1,1),Fid = Fid12
  coarse1=FLOAT(coarse1)

  FileName = temp_file+'\temp_C2'
  GetData,ImgData=coarse2,FileName = FileName+strtrim(isub+1,1),Fid = Fid13
  coarse2=FLOAT(coarse2)

  FileName = temp_file+'\class'
  GetData,ImgData=L1_class0,FileName = FileName+strtrim(isub+1,1),Fid = Fid14

  FileName = temp_file+'\temp_mad'
  GetData, ImgData=mad, FileName=FileName+strtrim(isub+1, 1), Fid=Fid16

  num_class=max(L1_class0)
  ;recode the classification map if the subset does not have all classes
  i_new_c=0
  L1_class=intarr(ns,nl)
  for iclass=0, num_class-1,1 do begin
    ind_ic=where(L1_class0 eq iclass+1 and fine1[*,*, background_band-1] ne background, num_ic)
    if (num_ic gt 0) then begin
      L1_class[ind_ic]=i_new_c+1
      i_new_c=i_new_c+1
    endif
  endfor

  num_class=max(L1_class)

  if (num_class gt 0) then begin   ;do not process if the whole subset is background
    ;correct extreme noise in fine1 becase extreme values will affect the allowed data range
    for ib=0,nb-1, 1 do begin
      sortIndex = Sort(fine1[*,*,ib])
      sortIndices = (Findgen(float(ns)*nl+1))/(float(ns)*nl)
      Percentiles=[0.0001, 0.9999]
      dataIndices = Value_Locate(sortIndices, Percentiles)
      data_1_4= (fine1[*,*,ib])[sortIndex[dataIndices]]
      ;correct too small values
      ind_small=where(fine1[*,*,ib] le data_1_4[0] or fine1[*,*,ib] lt DN_min)
      temp=fine1[*,*,ib]
      temp[ind_small]=min((fine1[*,*,ib])[where(fine1[*,*,ib] gt data_1_4[0] and fine1[*,*,ib] ge DN_min)])
      fine1[*,*,ib]=temp
      ;correct too large values
      ind_large=where(fine1[*,*,ib] ge data_1_4[1] or fine1[*,*,ib] gt DN_max)
      temp=fine1[*,*,ib]
      temp[ind_large]=max((fine1[*,*,ib])[where(fine1[*,*,ib] lt data_1_4[1] and fine1[*,*,ib] le DN_max)])
      fine1[*,*,ib]=temp
    endfor

    ;get index image between coarse and fine resolutions
    ii=0
    ns_c=floor(ns/scale_factor)
    nl_c=floor(nl/scale_factor)
    index_f=intarr(ns,nl)
    index_c=intarr(ns_c,nl_c)
    for i=0, ns_c-1, 1 do begin
      for j=0,nl_c-1,1 do begin
        index_f[i*scale_factor:(i+1)*scale_factor-1, j*scale_factor:(j+1)*scale_factor-1]=ii
        index_c[i,j]=ii
        ii=ii+1.0
      endfor
    endfor

    ;col and row index
    row_ind=intarr(ns,nl)
    col_ind=intarr(ns,nl)
    for i=0,ns-1,1 do begin
      col_ind[i,*]=i
    endfor
    for i=0,nl-1,1 do begin
      row_ind[*,i]=i
    endfor

    ;resample coarse image to coarse resolution
    fine_c1=fltarr(ns_c,nl_c,nb)
    coarse_c1=fltarr(ns_c,nl_c,nb)
    coarse_c2=fltarr(ns_c,nl_c,nb)
    row_c=fltarr(ns_c,nl_c)
    col_c=fltarr(ns_c,nl_c)
    for ic=0,ns_c-1, 1 do begin
      for jc=0,nl_c-1, 1 do begin
        ind_c=where(index_f eq index_c[ic,jc])
        row_c[ic,jc]= mean(row_ind[ind_c])
        col_c[ic,jc]= mean(col_ind[ind_c])
        for ib=0,nb-1,1 do begin
          fine_c1[ic,jc,ib]=mean((fine1[*,*,ib])[ind_c])
          coarse_c1[ic,jc,ib]=mean((coarse1[*,*,ib])[ind_c])
          coarse_c2[ic,jc,ib]=mean((coarse2[*,*,ib])[ind_c])
        endfor
      endfor
    endfor

    ;step 2: get fracture of each class within each coarse pixel at t1
    Fraction1=fltarr(ns_c,nl_c,num_class)
    for ic=0,ns_c-1, 1 do begin
      for jc=0,nl_c-1, 1 do begin
        ind_c=where(index_f eq index_c[ic,jc], num_c)
        L1_class_c=L1_class[ind_c]
        for iclass=0, num_class-1,1 do begin
          ind_ic=where(L1_class_c eq iclass+1, num_ic)
          Fraction1[ic,jc,iclass]=float(num_ic)/float(num_c)
        endfor
        if (total(Fraction1[ic,jc,*]) le 0.999) then begin   ;avoild pixels have background fine pixels
          Fraction1[ic,jc,*]=0
        endif
      endfor
    endfor

    ;get the heterogenity of each fine pixel
    het_index=fltarr(ns,nl)
    ;het_index1=fltarr(ns, nl, nb)
    scale_d=w

    for i=0,ns-1, 1 do begin
      for j=0,nl-1, 1 do begin
        ai=max([0,i-scale_d])       ; the window location
        bi=min([ns-1,i+scale_d])
        aj=max([0,j-scale_d])
        bj=min([nl-1,j+scale_d])
        class_t=L1_class[i,j]
        ;select same-class pixels
        ind_same_class=where(L1_class[ai:bi,aj:bj] eq class_t, num_sameclass)
        het_index[i,j]=float(num_sameclass)/((bi-ai+1.0)*(bj-aj+1.0))

      endfor
    endfor

    c_rate=fltarr(num_class,nb)

    ;allowed change value for each band
    min_allow=fltarr(nb)
    max_allow=fltarr(nb)
    for ib=0,nb-1,1 do begin
      min_allow[ib]=min(coarse_c2[*,*,ib]-coarse_c1[*,*,ib])-stddev(coarse_c2[*,*,ib]-coarse_c1[*,*,ib])
      max_allow[ib]=max(coarse_c2[*,*,ib]-coarse_c1[*,*,ib])+stddev(coarse_c2[*,*,ib]-coarse_c1[*,*,ib])
    endfor

    for ib=0,nb-1, 1 do begin
      X_matrix=fltarr(num_class,num_pure*num_class)
      y_matrix=fltarr(1,num_pure*num_class)
      ii=0
      for ic=1, num_class, 1 do begin
        order=sort(Fraction1[*,*,ic-1])
        order=reverse(order)
        ind_f=where(Fraction1[*,*,ic-1] gt 0.01, num_f) ;make sure all selected modis pixel contain class i
        num_pure1=min([num_f,num_pure])
        change_c=(coarse_c2[*,*,ib])[order[0:(num_pure1-1)]]-(coarse_c1[*,*,ib])[order[0:(num_pure1-1)]]

        ;only use 0.1-0.9 samples to exclude the land cover change pixels
        sortIndex = Sort(change_c)
        sortIndices = (Findgen(num_pure1+1))/float(num_pure1)
        Percentiles=[0.25, 0.75]
        ;Percentiles=[0.1, 0.9]
        dataIndices = Value_Locate(sortIndices, Percentiles)
        data_1_4= change_c[sortIndex[dataIndices]]
        ind_nonchange=where(change_c ge data_1_4[0] and change_c le data_1_4[1],num_nonc)
        y_matrix[0,ii:ii+num_nonc-1]=change_c[ind_nonchange]
        for icc=0, num_class-1,1 do begin
          f_c=(Fraction1[*,*,icc])[order[0:(num_pure1-1)]]
          X_matrix[icc,ii:ii+num_nonc-1]=f_c[ind_nonchange]
        endfor
        ii=ii+num_nonc
      endfor
      X_matrix=X_matrix[*,0:ii-1]
      y_matrix=y_matrix[0,0:ii-1]
      result=P_OLS(X_matrix,y_matrix,num_class,min_allow[ib],max_allow[ib])  ;format: ind_v:k*n,dep_v:1*n
      c_rate[*,ib]=result[*]
    endfor

    ;;step3: predict L2 assuming no land cover change
    L2_1=fine1
    for ic=1, num_class, 1 do begin
      ind_L1_class=where(L1_class eq ic)
      for ib=0,nb-1, 1 do begin
        temp=L2_1[*,*,ib]
        temp[ind_L1_class]=(fine1[*,*,ib])[ind_L1_class]+c_rate[ic-1,ib]
        L2_1[*,*,ib]=temp
      endfor
    endfor

    ;allowed minmum value for each band at t2
    min_allow=fltarr(nb)
    max_allow=fltarr(nb)
    for ib=0,nb-1,1 do begin
      min_allow0=min([min(coarse2[*,*,ib]),min(L2_1[*,*,ib])])
      min_allow[ib]=max([min_allow0,DN_min])
      max_allow0=max([max(coarse2[*,*,ib]),max(L2_1[*,*,ib])])
      max_allow[ib]=min([max_allow0,DN_max])
    endfor

    ;; step5: distribute changed values and residuals
    L2_1c=fltarr(ns_c, nl_c, nb)
    for ic=0, ns_c-1, 1 do begin
      for jc=0, nl_c-1, 1 do begin
        ind_c=where(index_f eq index_c[ic, jc])
        for ib=0, nb-1, 1 do begin
          L2_1c[ic, jc, ib]=mean((L2_1[*,*,ib])[ind_c])
        endfor
      endfor
    endfor

    w_mad=mad/change_mean

    real_change=coarse_c2-coarse_c1
    real_change2=coarse_c2-fine_c1
    predict_change=L2_1c-fine_c1
    residual=real_change-predict_change

    change_21=fltarr(ns, nl, nb)
    for ic=0, ns_c-1, 1 do begin
      for jc=0, nl_c-1, 1 do begin
        ind_c=where(index_f eq index_c[ic, jc], num_ii)

        for ib=0, nb-1, 1 do begin
          r_ch=real_change2[ic,jc,ib]    ;;;;;;;add
          diff_change=residual[ic, jc,ib]
          w_unform=fltarr(num_ii)
          w_unform[*]=abs(diff_change)
          ;w_change=predict_change[ind_c]+w_unform*(1.0-het_index[ind_c])*m_nochange[ind_c]+real_change*((w_mad[*,*,ib])[ind_c])*m_changed[ind_c]+0.000001
          ;w_change=w_unform*(1.0-het_index[ind_c])*m_nochange[ind_c]+r_ch*((w_mad[*,*,ib])[ind_c])*m_changed[ind_c]+0.000001 ;;;;
          w_change=w_unform*(1-het_index[ind_c])+r_ch*((w_mad[*,*,ib])[ind_c])
          w_change=w_change/(mean(w_change))

          ;avoid extreme weights
          ind_extrem=where(w_change gt 10, num_extrem)
          if (num_extrem gt 0) then begin
            w_change[ind_extrem]=mean(w_change)
          endif
          w_change=w_change/(mean(w_change))

          temp=change_21[*,*,ib]
          temp[ind_c]=w_change*diff_change  ; or diff_change
          change_21[*,*,ib]=temp
        endfor
      endfor
    endfor

    fine2_2=L2_1+change_21

    ;correct abnormal detected change
    for ib=0,nb-1, 1 do begin
      temp=fine2_2[*,*,ib]
      ind_min=where(temp lt min_allow[ib], num_min)
      if (num_min gt 0) then begin
        temp[ind_min]=min_allow[ib]
      endif
      ind_max=where(temp gt max_allow[ib], num_max)
      if (num_max gt 0) then begin
        temp[ind_max]=max_allow[ib]
      endif
      fine2_2[*,*,ib]=temp
    endfor

    change_21=fine2_2-fine1

  endif else begin
    change_21=fine1-fine1
  endelse
  ;
  change_21=change_21[location[0,isub]:location[1,isub],location[2,isub]:location[3,isub],*]
  change_21=float(change_21)

  print,'finish change prediction step ',isub+1,' block'

  tempoutname1=temp_file+'\temp_change'
  Envi_Write_Envi_File,change_21,Out_Name = tempoutname1+strtrim(isub+1,1)
  envi_file_mng, id=Fid11, /remove
  envi_file_mng, id=Fid12, /remove
  envi_file_mng, id=Fid13, /remove
  envi_file_mng, id=Fid14, /remove
  ENVI_FILE_MNG, ID=FID16, /REMOVE
  ;envi_file_mng, id=fid8, /remove

endfor
;
;;--------------------------------------------------------------------------------------
;mosiac all the change patch

mfid=intarr(n_ns*n_nl)
mdims=intarr(5,n_ns*n_nl)
mpos=intarr(nb,n_ns*n_nl)
pos=indgen(nb)
x0=intarr(n_ns*n_nl)
y0=intarr(n_ns*n_nl)

for isub=0,n_ns*n_nl-1,1 do begin
  envi_open_file, tempoutname1+strtrim(isub+1,1), r_fid= sub_fid
  if (sub_fid eq -1) then begin
    envi_batch_exit
    return
  endif
  envi_file_query,  sub_fid, ns=sub_ns, nl=sub_nl
  mfid[isub] = sub_fid
  mpos[*,isub] = indgen(nb)
  mdims[*,isub] = [-1,0, sub_ns-1,0, sub_nl-1]
  x0[isub]=ind_patch1[0,isub]
  y0[isub]=ind_patch1[2,isub]
endfor

xsize = orig_ns
ysize = orig_nl
pixel_size = [1.,1.]

use_see_through = replicate(1L,n_ns*n_nl)
see_through_val = replicate(0L,n_ns*n_nl)

out_name=temp_file+'_change'
envi_doit, 'mosaic_doit', fid=mfid, pos=mpos, $
  dims=mdims, out_name=out_name, xsize=xsize, $
  ysize=ysize, x0=x0, y0=y0, georef=0,MAP_INFO=map_info, $
  out_dt=4, pixel_size=pixel_size, $
  background=0, see_through_val=see_through_val, $
  use_see_through=use_see_through

for i=0,n_ns*n_nl-1,1 do begin
  envi_file_mng, id=mfid[i], /remove, /delete
endfor

;;##############step 6: final prediction

FileName6 = out_name
envi_open_file,FileName6,r_fid=fid
tempoutname=temp_file+'\temp_change'
pos=indgen(nb)
for isub=0,n_ns*n_nl-1,1 do begin
  dims=[-1,ind_patch[0,isub],ind_patch[1,isub],ind_patch[2,isub],ind_patch[3,isub]]
  envi_doit, 'resize_doit', fid=fid, pos=pos, dims=dims, interp=0, rfact=[1,1], $
    out_name=tempoutname+strtrim(isub+1,1), r_fid=r_fid1
  envi_file_mng, id=r_fid1, /remove
endfor
envi_file_mng, id=fid, /remove,/delete


for isub=0,n_ns*n_nl-1,1 do begin

  ;open each block image

  FileName = temp_file+'\temp_F1'
  GetData,ImgData=fine1,ns = ns,nl = nl,nb = nb,Data_Type = Data_Type,FileName = FileName+strtrim(isub+1,1),Fid = Fid11
  fine1=float(fine1)

  FileName = temp_file+'\temp_C1'
  GetData,ImgData=coarse1,FileName = FileName+strtrim(isub+1,1),Fid = Fid12
  coarse1=FLOAT(coarse1)

  FileName = temp_file+'\temp_C2'
  GetData,ImgData=coarse2,FileName = FileName+strtrim(isub+1,1),Fid = Fid13
  coarse2=FLOAT(coarse2)

  FileName = temp_file+'\class'
  GetData,ImgData=L1_class,FileName = FileName+strtrim(isub+1,1),Fid = Fid14

  FileName = temp_file+'\temp_change'
  GetData,ImgData=change_21,FileName = FileName+strtrim(isub+1,1),Fid = Fid15
  change_21=FLOAT(change_21)
  ;place the blended result
  fine2=fine1[location[0,isub]:location[1,isub],location[2,isub]:location[3,isub],*]


  ;compute the distance of each pixel in the window with the target pixel (integrate window)
  D_D_all=fltarr((w*2+1),(w*2+1))
  for jw=0,w*2,1 do begin
    for iw=0,w*2,1 do begin
      D_D_all[iw,jw]=((w-iw)^2+(w-jw)^2)^0.5+0.00001
    endfor
  endfor
  D_D_all=reform(D_D_all,(w*2+1)*(w*2+1))

  for i=location[0,isub],location[1,isub],1 do begin           ;retieve each target pixel
    for j=location[2,isub],location[3,isub],1 do begin
      if (fine1[i,j,background_band-1] ne background  ) then begin    ;do not process the background

        ai=max([0,i-w])       ; the window location
        bi=min([ns-1,i+w])
        aj=max([0,j-w])
        bj=min([nl-1,j+w])

        ci=i-ai      ;location of target pixel
        cj=j-aj

        class_t=L1_class[i,j]
        ;select same-class pixels
        ind_same_class=where(L1_class[ai:bi,aj:bj] eq class_t, num_sameclass)
        fine1_wind=fine1[ai:bi,aj:bj,*]

        ;compute weight for each simialr pixel
        D_D_cand=fltarr(num_sameclass)        ;spatial distance
        if ((bi-ai+1)*(bj-aj+1) lt (w*2.0+1)*(w*2.0+1)) then begin   ;not integrate window
          for icand=0,num_sameclass-1,1 do begin
            iw=ai+(ind_same_class[icand] mod (bi-ai+1))
            jw=aj+(ind_same_class[icand]/(bi-ai+1))
            D_D_cand[icand]=((i-iw)^2+(j-jw)^2)^0.5+0.00001                           ;smaller value indicates closer
          endfor
        endif else begin
          D_D_cand=D_D_all[ind_same_class]      ;integrate window
        endelse

        ;normalize these distances
        D_D_cand=1.0+D_D_cand/w
        C_D=1.0/D_D_cand
        weight=C_D/total(C_D)

        for iband=0,nb-1,1 do begin
          ;predict the value
          change_cand=(change_21[ai:bi,aj:bj,iband])[ind_same_class]
          fine2[i-location[0,isub],j-location[2,isub],iband]=fine1[i,j,iband]+total(weight*change_cand)
          ;fine2[i-location[0,isub],j-location[2,isub],iband]=fine1[i,j,iband]+change_cand)

          ;revise the abnormal prediction
          if (fine2[i-location[0,isub],j-location[2,isub],iband] lt DN_min ) then begin ;correct abnomal prediction
            another_predict=max([DN_min, fine1[i,j,iband]+(coarse2[i,j,iband]-coarse1[i,j,iband])])
            fine2[i-location[0,isub],j-location[2,isub],iband]=min([DN_max,another_predict])
          endif
          if (fine2[i-location[0,isub],j-location[2,isub],iband] gt DN_max ) then begin ;correct abnomal prediction
            another_predict=min([DN_max, fine1[i,j,iband]+(coarse2[i,j,iband]-coarse1[i,j,iband])])
            fine2[i-location[0,isub],j-location[2,isub],iband]=max([DN_min,another_predict])
          endif
        endfor

      endif
    endfor

  endfor

  ; change the type of prediction into the type same as the input image
  case Data_Type Of
    1:fine2 = Byte(fine2)    ;  BYTE  Byte
    2:fine2 = FIX(fine2)     ;  INT  Integer
    3:fine2 = LONG(fine2)    ;  LONG  Longword integer
    4:fine2 = FLOAT(fine2)   ;  FLOAT  Floating point
    5:fine2 = DOUBLE(fine2)  ;  DOUBLE  Double-precision floating
    6:fine2 = COMPLEX(fine2); complex, single-precision, floating-point
    9:fine2 = DCOMPLEX(fine2);complex, double-precision, floating-point
    12:fine2 = UINT(fine2)   ; unsigned integer vector or array
    13:fine2 = ULONG(fine2)   ;  unsigned longword integer vector or array
    14:fine2 = LONG64(fine2)   ;a 64-bit integer vector or array
    15:fine2 = ULONG64(fine2)   ;an unsigned 64-bit integer vector or array
  EndCase

  print,'finish final prediction ',isub+1,' block'
  tempoutname1=temp_file+'\temp_blended'
  Envi_Write_Envi_File,fine2,Out_Name = tempoutname1+strtrim(isub+1,1)
  envi_file_mng, id=Fid11, /remove, /delete
  envi_file_mng, id=Fid12, /remove, /delete
  envi_file_mng, id=Fid13, /remove, /delete
  envi_file_mng, id=Fid14, /remove, /delete
  envi_file_mng, id=Fid15, /remove, /delete

endfor

print, 'time used:', floor((systime(1)-t0)/3600), 'h',floor(((systime(1)-t0) mod 3600)/60),'m',(systime(1)-t0) mod 60,'s'
;#########################################
;mosiac all the blended patch

mfid=intarr(n_ns*n_nl)
mdims=intarr(5,n_ns*n_nl)
mpos=intarr(nb,n_ns*n_nl)
pos=indgen(nb)
x0=intarr(n_ns*n_nl)
y0=intarr(n_ns*n_nl)

for isub=0,n_ns*n_nl-1,1 do begin
  envi_open_file, tempoutname1+strtrim(isub+1,1), r_fid= sub_fid
  if (sub_fid eq -1) then begin
    envi_batch_exit
    return
  endif
  envi_file_query,  sub_fid, ns=sub_ns, nl=sub_nl
  mfid[isub] = sub_fid
  mpos[*,isub] = indgen(nb)
  mdims[*,isub] = [-1,0, sub_ns-1,0, sub_nl-1]
  x0[isub]=ind_patch1[0,isub]
  y0[isub]=ind_patch1[2,isub]
endfor

xsize = orig_ns
ysize = orig_nl
pixel_size = [1.,1.]

use_see_through = replicate(1L,n_ns*n_nl)
see_through_val = replicate(0L,n_ns*n_nl)

out_name=Dialog_PickFile(Title = 'Enter the filename of the blended image')

envi_doit, 'mosaic_doit', fid=mfid, pos=mpos, $
  dims=mdims, out_name=out_name, xsize=xsize, $
  ysize=ysize, x0=x0, y0=y0, georef=0,MAP_INFO=map_info, $
  out_dt=Data_Type, pixel_size=pixel_size, $
  background=0, see_through_val=see_through_val, $
  use_see_through=use_see_through

for i=0,n_ns*n_nl-1,1 do begin
  envi_file_mng, id=mfid[i], /remove, /delete
endfor
 
END  