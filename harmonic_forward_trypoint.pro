function HARMONIC_FORWARD_TRYPOINT,simplex,y,psum,windex, $
  fac,images,ctrs,obspos,res,maxlvar,occultr, $
  penalty_only=penalty_only


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;                                                            ;;
;; Purpose: Alters one of the vertices in the simplex var-    ;;
;;            iable.  If the altered vertex gives a lower pen-;;
;;            alty function value it replaces the old vertex  ;;
;;            in the simplex.  Alternatively, this routine can;;
;;            be used to simply evaluate the objective        ;;
;;            function at the specified vertex. Called by     ;;
;;            forward_comparison_optimization.  The objective ;;
;;            function in this case is given by the inverse of;;
;;            the sums of the correlations between radial and ;;
;;            azimuthal strips from the input images and      ;;
;;            forward-modeled images.                         ;;
;;                                                            ;;
;; Inputs: simplex - nvert x (nvert+1) array of locations in  ;;
;;           the solution space where the objective function  ;;
;;           is known - used to find optimal solution         ;;
;;         psum - "center of mass" of the vertices in the     ;;
;;           simplex array = total(simplex,2) before perturb. ;;
;;         y - value of penalty function corresponding to each;;
;;           vertex in simplex                                ;;
;;         windex - index of the vertex of interest; where    ;;
;;           objective function should be calculated or vertex;;
;;           should be improved                               ;;
;;         fac - factor by which to alter position of worst   ;;
;;           vertex; or, set to 1 if penalty_only keyword is  ;;
;;           to be used                                       ;;
;;         images - series of images to compare with the field;;
;;           model                                            ;;
;;         ctrs - locations of sun center (in pixels) in each ;;
;;           image                                            ;;
;;         obspos - 2xn array giving latitude, longitude of   ;;
;;           the imaging instrument that produced each image, ;;
;;           in degrees                                       ;;
;;         res - 2xn array giving the x-,y-resolution of each ;;
;;           image in arcsec/pixel                            ;; 
;;         maxlvar - maximum l value of harmonic coefficients ;;
;;           being optimized; needed to convert simplex vert- ;;
;;           ices back into transforms                        ;;
;;         occultr - approximate height (in pixels) of        ;;
;;            occulter in each image                          ;;
;;                                                            ;;
;;                                                            ;;
;; Returns: Value of objective function associated with the   ;;
;;            vertex in simplex[*,fac]                        ;;
;;                                                            ;;
;; Keywords: penalty_only - calculate objective function asso-;;
;;             ciated with the vertex in simplex[*,fac], with ;;
;;             no alteration                                  ;;
;;                                                            ;;
;; Dependencies: extract_transform, calc_phis,pfss software   ;;
;;                 library from SolarSoft                     ;;
;;                                                            ;;
;; Created: 02/08/2017                                        ;;
;;                                                            ;;
;; Created by: Shaela Jones                                   ;;
;;                                                            ;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

; grab common block variables
@pfss_opt_parameters
@pfss_data_block
@pfss_forward_vars


; initializations
nlat=N_ELEMENTS(lat)
nrix=N_ELEMENTS(rix)
szsimplex=SIZE(simplex)
szimages=SIZE(images)
nimages=szimages[3]     ; check this
nx=szimages[1]
ny=szimage[2]


; calculate minimum, maximum x,y of input images
xrange=res*nx
yrange=res*ny
xxmin=-xrange/2.
xxmax=xrange/2.
yymin=-yrange/2.
yymax=yrange/2.



; find new vertex
fac1=(1-fac)/N_ELEMENTS(psum)
fac2=fac1-fac
new_vertex=psum*fac1-simplex[*,windex]*fac2


; re-configure transform, calculate phiat,phibt in pfss common block
newmagt=EXTRACT_TRANSFORM(new_vertex,magt,maxlvar)
CALC_PHIS,newmagt


; extrapolate magnetic field of vertex, forward model images
PFSS_POTL_FIELD,rss,rgrid,/trunc,/quiet
savfile='field_temp.sav'
save,BPH,BR,BTH,I,LAT,LON,NLAT,NLON,NR,PHI,PHIAT,PHIBT, $
     RIX,THETA,filename=savfile


; convert forward-modeled images to polar images
model_polarims=FLTARR(nimages,nx,ny)
for i=0,nimages-1 do begin
  FOR_DRIVE,'pfssmod',inst='WL',cmer=obspos[1,i], $
       bang=obspos[0,i],ngrid=nx,ngy=ny, $
       rindex=rix,quantmap=modelimi,occult=occultr[i], $
       xxmin=xxmin[i],xxmax=xxmax[i],yymin=yymin[i], $
       yymax=yymax[i]
  model_polarims[i,*,*]=POLAR_IMAGE(modelimi,[nx/2.,ny/2.], $
       nrix,nlat)
 ; polar_image will create an evenly spaced r grid - okay with this?       
; i think there's a way to do the model calculation once and then 
;     rotate it several times to calculate the image from different
;     perspectives; how to do this?  how much time would it save?
endfor



; calculate correlation coefficients
; faster way to do this calculation?
; may have a problem with fully or partialy occulted pixels
rcorr=FLTARR(nrix)
thetacorr=FLTARR(nlat)
for i=0,nrix-1 do begin
  rcorr[i]=CORRELATION(REFORM(model_polarims[*,i]), $
        REFORM(constraint_polarims[*,i]))
endfor
for j=0,nlat-1 do begin
  thetacorr[j]=CORRELATION(REFORM(model_polarims[j,*]), $
       REFORM(constraint_polarims[j,*]))
endfor



; calculate objective function
objective_func=1./TOTAL(1+rcorr)+1./TOTAL(1+thetacorr)
objective_func=normalization*objective_func


; add more terms later?


return,objective_func
end
