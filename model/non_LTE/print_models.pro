pro print_model_to_file, filename, lambda, flux, T_eff, log_g

  openw, lun, filename, /get_lun
  printf, lun, T_eff, format='(%"# T_eff=%d")'
  printf, lun, log_g, format='(%"# log_g=%0.1f")'
  printf, lun, format='(%"# lambda\tmodel_flux")'
  for i=0, n_elements(lambda)-1 do printf, lun, lambda[i], flux[i], format='(%"%0.1f\t%g")'
  free_lun, lun

end

pro print_models
  files=[ $
'0050200_7.80_00005-55000.WFP', $
'0050200_8.00_00005-55000.WFP', $
'0050200_8.20_00005-55000.WFP', $
'0050500_7.80_00005-55000.WFP', $
'0050500_8.00_00005-55000.WFP', $
'0050500_8.20_00005-55000.WFP', $
'0050800_7.80_00005-55000.WFP', $
'0050800_8.00_00005-55000.WFP', $
'0050800_8.20_00005-55000.WFP', $
'0051100_7.80_00005-55000.WFP', $
'0051100_8.00_00005-55000.WFP', $
'0051100_8.20_00005-55000.WFP', $
'0051400_7.80_00005-55000.WFP', $
'0051400_8.00_00005-55000.WFP', $
'0051400_8.20_00005-55000.WFP', $
'sed510_79' $
]

  T_eff = [ $
           50200, 50200, 50200, $
           50500, 50500, 50500, $
           50800, 50800, 50800, $
           51100, 51100, 51100, $
           51400, 51400, 51400, $
           51000 $
          ]

  log_g = [ $
           7.8, 8.0, 8.2, $
           7.8, 8.0, 8.2, $
           7.8, 8.0, 8.2, $
           7.8, 8.0, 8.2, $
           7.8, 8.0, 8.2, $
           7.9 $
          ]

  restore,'non_LTE_models.sav'

  for i=0,n_elements(T_eff)-1 do begin
    filename = string(i, format='(%"%02d")') + '.model'
    print_model_to_file, filename, sedwvl, *sedspec[i], T_eff[i], log_g[i]
  endfor
end
