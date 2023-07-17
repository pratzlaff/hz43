pro print_to_file, filename, lambda, model, ratio, err
  openw, lun, filename, /get_lun
  printf, lun, format='(%"# lambda\tmodel_flux\tratio\terr")'
  for j=0, n_elements(lambda)-1 do printf, lun, lambda[j], model[j], ratio[j], err[j], format='(%"%0.1f\t%g\t%g\t%g")'
  free_lun, lun
end

pro print_ratios
;; From Nick's 2014-11-14 email:
;;
;; Hi Pete,
;; I think that this should do it for you-
;;
;; restore, '/data/hal9000/HZ43/SAV/data_model_ratio.sav'
;;
;; ;Neg wvl, ratio, and error arrays
;; *wvlmb[i]
;; *rat_m[i]
;; *rat_me[i]
;;
;; ;Pos wvl, ratio, and error arrays

;; *wvlpb[i]
;; *rat_p[i]
;; *rat_pe[i]
;;
;; ;Key
;; i=0   Black curve
;; i=19  Green curve
;; i=20  Red curve
;; i=21  Blue curve

  restore, 'data_model_ratio.sav'

  index = [15, 16, 19, 20, 21]
  years = ['2010', '2011', '2012', '2013', '2014']

  for j=0,n_elements(index)-1 do begin
    i = index[j]
    print_to_file, years[j]+'_pos.ratio', *wvlpb[i], *sedspec_pb[i], *rat_p[i], *rat_pe[i]
    print_to_file, years[j]+'_neg.ratio', *wvlmb[i], *sedspec_mb[i], *rat_m[i], *rat_me[i]
  endfor
end
