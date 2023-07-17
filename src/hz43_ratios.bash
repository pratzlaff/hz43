#! /bin/bash

export PYTHONPATH=/data/legs/rpete/flight/analysis_functions

for iqe in N0008 N0009 N0010 N0011 N0012
do
  export ARFPATH=/data/legs/rpete/flight/hz43/arfs/I_qe_${iqe}:/data/legs/rpete/flight/hrcs_qeu/ARD/v14/arfs/S_qeu_N0014
 /usr/bin/python3 /data/legs/rpete/flight/hrcs_qe/N0015/ECR/src/hz43_s_rates.py  -r -o plots/hz43_hrc_ratios_i_qe_${iqe}.pdf
done
