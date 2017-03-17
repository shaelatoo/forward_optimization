function forward_comparison_optimization,images,hdrs,magfile, $
    scale,maxits=maxits,simplex=simplex, $
    tims=tims,min_values=min_values


;;;;;   UNTESTED   ;;;;;;;;;

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;                                                            ;;
;; Purpose: This program optimizes a PFSS magnetic field model;;
;;            using comparison between an input image series  ;;
;;            and forward-modeled images based on the model.  ;;
;;            The comparison is quantitative, consisting of a ;;
;;            calculation of the correlation between angular  ;;
;;            and radial strips from each image.              ;;
;;                                                            ;;
;; Inputs: images - a series of coronal images to which the   ;;
;;            model should be compared                        ;;
;;         hdrs - array of fits header files corresponding to ;;
;;           the images                                       ;;
;;         magfile - magnetogram file to use to create the    ;;
;;            initial model                                   ;;
;;         scale - parameter to set scale size of initial     ;;
;;            simplex - set relative to the expected range of ;;
;;            the magnetogram coefficients                    ;;
;;         maxlvar - maximum l value of spherical harmonic    ;;
;;            coefficients to be optimized                    ;;
;;         maxits - maximum number of times to run through the;;
;;            update loop before giving up on convergence     ;;
;;                                                            ;;
;; Returns: magnetogram corresponding to optimized field      ;;
;;                                                            ;;
;; Keyword Outputs:                                           ;;
;;   simplex: (nvert x nvert+1) array of vertices.  Each      ;;
;;      vertex represents a point in the solution space at    ;;
;;      which the penalty function was evaluated.  if the     ;;
;;      optimization converged, the first image,              ;;
;;      simplex[*,*,0], will contain the solution correspond- ;;
;;      ing to the minimum penalty function.  otherwise, it   ;;
;;      gives some idea where in the solution space the algor-;;
;;      ithm was looking.  nvert is determined by the maxlvar ;;
;;      keyword, nvert=(maxlvert+1)*(maxlvert+2); as simplex  ;;
;;      grows the computation becames more time-consuming and ;;
;;      less accurate                                         ;;
;;   tims - array giving the time (in seconds) since beginning;; 
;;      optimization, with a time corresponding to each       ;;
;;      element of min_values
;;   min_values - array giving the lowest value in the func-  ;;
;;      tion_value array at the end of each run through the   ;;
;;      optimization loop
;;                                                            ;;
;; Keywords: ;;
;;             ;; 
;;                                                            ;;
;; Dependencies: ;;
;;                                                            ;;
;; Created: 02/08/2017                                        ;;
;;                                                            ;;
;; Created by: Shaela Jones                                   ;;
;;                                                            ;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

; initialize common block variables
@pfss_opt_parameters
@pfss_data_block


; parameters
defaultscale=3.
defaultmaxits=5000L
expansion=2.0          ; expansion, reflection,contraction,reduction are parameters of the amoeba optimization algorithm
reflection=-1.0
contraction=0.5
reduction=0.5
tc_default=2.*7.*24.*60.*60.   ; default time cutoff: two weeks!
minrange=0.001     ; when y values become very similar, algorithm is terminated
;lats=[]
;lons=[] 


; check that pfss_opt_parameters have been defined
if N_ELEMENTS(rss) eq 0 or N_ELEMENTS(rgrid) eq 0 or $
  N_ELEMENTS(magtype) eq 0 or N_ELEMENTS(nlat0) eq 0 or $
  N_ELEMENTS(noreset) eq 0 then begin
  print,'Important variables from pfss_opt_parameters not defined'
  return,-1
endif
lmax=nlat0    ;; make lmax a keyword?
if N_ELEMENTS(maxlvar) eq 0 then begin
  print,'Maxlvar not set.
  return,-1
endif


; check that user has provided header files
if N_ELEMENTS(hdrs) eq 0 then begin
  print,'Please provide header file data.'
  return,-1
endif


; initializations
if n_elements(maxits) eq 0 then maxits = defaultmaxits
if N_ELEMENTS(scale) eq 0 then scale=defaultscale
min_values=DBLARR(maxits)
tims=DBLARR(maxits)
if NOT(KEYWORD_SET(time_cutoff)) then time_cutoff=tc_default
nangles=N_ELEMENTS(angles)


; read magfile and initialize some elements of data block
PFSS_MAG_CREATE_SJ,magnetogram,magtype,nlat0,file=magfile,/quiet
magnetogram=magnetogram-mean(magnetogram)
PFSS_GET_POTL_COEFFS,magnetogram,rtop=rss
nlat=n_elements(theta)
nlon=2*nlat
cth=cos(theta)
magmov=magnetogram


; determine observer position, image center, and resolution for each image
EXTRACT_HEADER_DATA,hdrs,ctrs=ctrs,res=res,obspos=obspos, $
     occultr=occultr
if obspos[0] eq -1 then begin    ; observer position is not in header file - should be user-provided
  obspos=[[lats],[lons]]
endif


; fill "now" variable in PFSS structure
image=MRDFITS(magfile,0,hdr)
magwcs=FITSHEAD2WCS(hdr)
now=magwcs.time.fits_date


; initialize simplex & magt
magt=SPHERICAL_TRANSFORM(magnetogram,cth,lmax=lmax)
realmagt=REAL_PART(magt)
imagmagt=IMAGINARY(magt)
smagt=REFORM(realmagt[1,0:1])
for i=2,maxlvar do smagt=[smagt,REFORM(realmagt[i,0:i])]
for i=1,maxlvar do smagt=[smagt,REFORM(imagmagt[i,1:i])]
sim0=smagt
nvert=N_ELEMENTS(sim0)
simplex=sim0#REPLICATE(1.0,nvert+1)
for i=0,nvert-1 do simplex[i,i+1]=sim0[i]+ $
  scale[i<(N_ELEMENTS(scale)-1)]


; initialize vector of penalty function values
psum=0
y=FLTARR(nvert+1)
for i=0,nvert do y[i]=HARMONIC_FORWARD_TRYPOINT(simplex,y, $
  psum,i,1.,images,ctrs,obspos,res,occultr, $
  /penalty_only)
ncalls=LONG(nvert)    ; number of penalty function calculations




; optimization loop
cnt=0L
psum = TOTAL(simplex,2)
TIC
while cnt lt maxits do begin   ;Each iteration
  ord = SORT(y)
  lowest = ord[0]    ;Lowest point
  min_values[cnt]=y[lowest]   ; keeps track of min value each iteration
  highest = ord[nvert]   ;Highest point
  next_highest = ord[nvert-1]  ;Next highest point
  range = abs(y[highest]) + abs(y[lowest]) ;Denominator = interval
  if range ge minrange then rtol = 2.0 * abs(y[highest]-y[lowest])/range $
  else rtol = ftol / 2.   ;Terminate if interval is 0
  
  if rtol lt ftol then begin ;Done?
    t = y[0] & y[0] = y[lowest] & y[lowest] = t ;Sort so fcn min is 0th elem
    t = simplex[*,lowest] & simplex[*,lowest] = simplex[*,0] & simplex[*,0] = t
    print,'Simplex minimum size has been reached'
    newmagt=EXTRACT_TRANSFORM(reform(t),magt,maxlvar)
    result=INV_SPHERICAL_TRANSFORM(newmagt,cth)
    noreset=0
    return,result
  endif
  
  if TOC() gt time_cutoff then begin   ; time cutoff reached?
    t = y[0] & y[0] = y[lowest] & y[lowest] = t ;Sort so fcn min is 0th elem
    t = simplex[*,lowest] & simplex[*,lowest] = simplex[*,0] & simplex[*,0] = t
    print,'Optimization time cutoff has been reached'
    newmagt=EXTRACT_TRANSFORM(REFORM(t),magt,maxlvar)
    result=INV_SPHERICAL_TRANSFORM(newmagt,cth)
    noreset=0
    return, result
  endif
  
  ; try a reflection
  ytry=HARMONIC_FORWARD_TRYPOINT(simplex,y,psum,highest, $
     reflection,images,spclatlon,ctrs,obspos,res, $
     occultr,/penalty_only)
  ncalls++
  
  ; if ytry is better than the best point, expand
  if ytry le y[lowest] then begin
    ytry=HARMONIC_FORWARD_TRYPOINT(simplex,y,psum,highest, $
       expansion,images,spclatlon,ctrs,obspos,res, $
       occultr,/penalty_only)
    ncalls++
  endif else if ytry ge y[next_highest] then begin
    ; if ytry is the new worst point, contract
    ysave = y[highest]
    ytry=HARMONIC_FORWARD_TRYPOINT(simplex,y,psum,highest, $
      contraction,images,spclatlon,ctrs,obspos,res, $
      occultr,/penalty_only)
    ncalls++
    if ytry ge ysave then begin
      ; if the contracted point is still the new worst, reduce
      for i=0,nvert do begin
        if i ne lowest then begin
          psum = reduction * (simplex[*,i] + simplex[*,lowest])
          simplex[*,i] = psum
          y[i]=HARMONIC_FORWARD_TRYPOINT(simplex,y,psum,i, $
            1.0,images,spclatlon,ctrs,obspos,res, $
            occultr,/penalty_only)
        endif
      endfor
      ncalls = ncalls + nvert
      psum = TOTAL(simplex,2)
    endif   ;ytry ge ysave
  endif
  tims[cnt]=TOC()
  if (cnt mod 100) eq 0 then begin
    print,"Elapsed time = ",tims[cnt]/60.,' mins'
    print,'Penalty function at best vertex = ',min_values[cnt]
  endif
  cnt++
endwhile


; if reached this point, algorithm failed to converge
;   return best known point
print,'phase_varying_amoeba failed to converge to a solution in ', $
  'the specified number of iterations.'
t = y[0] & y[0] = y[lowest] & y[lowest] = t ;Sort so fcn min is 0th elem
t = simplex[*,lowest] & simplex[*,lowest] = simplex[*,0] & simplex[*,0] = t
newmagt=EXTRACT_TRANSFORM(REFORM(t),magt,maxlvar)
result=INV_SPHERICAL_TRANSFORM(newmagt,cth)
return, result  



end
