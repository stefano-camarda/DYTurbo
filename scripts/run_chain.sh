#!/bin/bash



#./scripts/merge.sh --proc wp,wm --term RES2P,CT2P --qtymerge qt010y-11 --outdir CT10nnlo_RESCT2P_160512
#./scripts/merge.sh --proc wp,wm,z0 --term RES2P,CT2P --qtymerge qt010y-11 --outdir CT10nnlo_RESCT2P_160512 --find_missing

#RERUN:
# ./scripts/submit_DYTURBO.sh --mogon --rerun dyturbo_wp_lhc7_CT10nnlo_array_o2qt3040y-11tRES2P_seed 130,131,132,133,134,135,136,137 RUN
# ./scripts/submit_DYTURBO.sh --mogon --rerun dyturbo_wp_lhc7_CT10nnlo_array_o2qt4050y-5-3tRES2P_seed 119 RUN
# ./scripts/submit_DYTURBO.sh --mogon --rerun dyturbo_wp_lhc7_CT10nnlo_array_o2qt4050y-11tRES2P_seed 142 RUN
# ./scripts/submit_DYTURBO.sh --mogon --rerun dyturbo_wp_lhc7_CT10nnlo_array_o2qt4050y13tRES2P_seed 106,122,136,148 RUN
# ./scripts/submit_DYTURBO.sh --mogon --rerun dyturbo_wp_lhc7_CT10nnlo_array_o2qt4050y35tRES2P_seed 116,120 RUN
# ./scripts/submit_DYTURBO.sh --mogon --rerun dyturbo_wp_lhc7_CT10nnlo_array_o2qt8090y-11tRES2P_seed 118,120 RUN

# ./scripts/submit_DYTURBO.sh --mogon --rerun dyturbo_wp_lhc7_CT10nnlo_array_o2qt010y-3-1tCT2P_seed 121 RUN
# ./scripts/submit_DYTURBO.sh --mogon --rerun dyturbo_wp_lhc7_CT10nnlo_array_o2qt5060y-3-1tCT2P_seed 120 RUN

# ./scripts/submit_DYTURBO.sh --mogon --rerun dyturbo_wm_lhc7_CT10nnlo_array_o2qt1020y-11tRES2P_seed 108,109,112 RUN
# ./scripts/submit_DYTURBO.sh --mogon --rerun dyturbo_wm_lhc7_CT10nnlo_array_o2qt5060y13tRES2P_seed 150 RUN
# ./scripts/submit_DYTURBO.sh --mogon --rerun dyturbo_wm_lhc7_CT10nnlo_array_o2qt6070y35tRES2P_seed 103,138,139 RUN
# ./scripts/submit_DYTURBO.sh --mogon --rerun dyturbo_wm_lhc7_CT10nnlo_array_o2qt8090y-11tRES2P_seed 136 RUN
# ./scripts/submit_DYTURBO.sh --mogon --rerun dyturbo_wm_lhc7_CT10nnlo_array_o2qt8090y13tRES2P_seed 101,128 RUN
# ./scripts/submit_DYTURBO.sh --mogon --rerun dyturbo_wm_lhc7_CT10nnlo_array_o2qt90100y-3-1tRES2P_seed 143 RUN
# ./scripts/submit_DYTURBO.sh --mogon --rerun dyturbo_wm_lhc7_CT10nnlo_array_o2qt90100y35tRES2P_seed 113 RUN

# ./scripts/submit_DYTURBO.sh --mogon --rerun dyturbo_wm_lhc7_CT10nnlo_array_o2qt4050y35tCT2P_seed 120,121,122,123,124,125,126,127,128,129,130,131,132,133,134,135,136,137,138,139,140,141,142,143,144,145,146,147,148,149,150 RUN
# ./scripts/submit_DYTURBO.sh --mogon --rerun dyturbo_wm_lhc7_CT10nnlo_array_o2qt5060y-5-3tCT2P_seed 100 RUN

# RERUN:
#./scripts/submit_DYTURBO.sh --mogon --rerun dyturbo_wm_lhc7_CT10nnlo_array_o2qt6070y35tRES2P_seed 103 RUN
#./scripts/submit_DYTURBO.sh --mogon --rerun dyturbo_wm_lhc7_CT10nnlo_array_o2qt7080y35tRES2P_seed 138,139 RUN
#./scripts/submit_DYTURBO.sh --mogon --rerun dyturbo_wm_lhc7_CT10nnlo_array_o2qt90100y35tRES2P_seed 113 RUN
#./scripts/submit_DYTURBO.sh --mogon --rerun dyturbo_wm_lhc7_CT10nnlo_array_o2qt4050y35tCT2P_seed 121,136,148,149 RUN

#./scripts/submit_DYTURBO.sh --mogon --proc wp --term RES2D,CT2D --pdfset MMHT2014nnlo68cl --pdfvar array
#./scripts/submit_DYTURBO.sh --mogon --proc wm --term RES2D,CT2D --pdfset MMHT2014nnlo68cl --pdfvar array

#./scripts/merge.sh --proc wp,wm --term RES2P,CT2P --pdfset CT10nnlo --qtymerge qt010y-11 --outdir CT10nnlo_RESCT2D_160523 RUN

# DONE -- 2P wp,wm CT10nnlo

# 2D
## CT10nnlo
#./scripts/submit_DYTURBO.sh --mogon --proc z0 --term RES2D,CT2D --pdfset CT10nnlo --pdfvar array RUN
#./scripts/merge.sh --proc z0 --term RES2D,CT2D --pdfset CT10nnlo --qtymerge qt05y-11 --outdir CT10nnlo_RESCT2D_160523 RUN
####./scripts/submit_DYTURBO.sh --mogon --proc wp,wm --term RES2D,CT2D --pdfset CT10nnlo --pdfvar array RUN
####./scripts/merge.sh --proc wp,wm --term RES2D,CT2D --pdfset CT10nnlo --qtymerge qt05y-11 --outdir CT10nnlo_RESCT2D_160523 RUN

#./scripts/merge.sh --proc z0 --term RES2D,CT2D --pdfset CT10nnlo --qtymerge qt510y-11 --outdir CT10nnlo_RESCT2D_160523
#./scripts/merge.sh --proc wp,wm --term RES2P,CT2P --pdfset CT10nnlo --qtymerge qt010y-11 --outdir CT10nnlo_RESCT2D_160523

# DONE -- 2D z0 CT10nnlo

## MMHT
#./scripts/submit_DYTURBO.sh --mogon --proc z0 --term RES2D,CT2D --pdfset MMHT2014nnlo68cl --pdfvar array RUN
#./scripts/merge.sh --proc z0 --term RES2D,CT2D --pdfset MMHT2014nnlo68cl --qtymerge qt05y-11 --outdir MMHT14_RESCT2D_160523
# DONE z0 MMHT14 

#./scripts/submit_DYTURBO.sh --mogon --proc wp,wm --term RES2D,CT2D --pdfset MMHT2014nnlo68cl --pdfvar array RUN
#./scripts/merge.sh --proc wp,wm --term RES2D,CT2D --pdfset MMHT2014nnlo68cl --qtymerge qt05y-11 --outdir MMHT14_RESCT2D_160523 RUN

## CT14NNLO
#./scripts/submit_DYTURBO.sh --mogon --proc z0 --term RES2D,CT2D --pdfset CT14nnlo --pdfvar array RUN
#./scripts/merge.sh --proc z0 --term RES2D,CT2D --pdfset CT14nnlo --qtymerge qt05y-11 --outdir CT14nnlo_RESCT2D_160523 RUN
#./scripts/submit_DYTURBO.sh --mogon --proc wp,wm --term RES2D,CT2D --pdfset CT14nnlo --pdfvar array RUN
#./scripts/merge.sh --proc wp,wm,z0 --term RES2D,CT2D --pdfset CT14nnlo --qtymerge qt05y-11 --outdir CT14nnlo_RESCT2D_160523 #RUN

#RERUN:
#./scripts/submit_DYTURBO.sh --mogon --proc z0 --term RES2P,CT2P --pdfset CT10nnlo --pdfvar array RUN
#./scripts/merge.sh --proc z0 --term RES2P,CT2P --qtymerge qt05y-11 --outdir CT10nnlo_RESCT2P_160523 RUN

#./scripts/submit_DYTURBO.sh --mogon --proc wp,wm --term RES2P,CT2P --pdfset CT10nnlo --pdfvar array
#./scripts/merge.sh --proc wp,wm --term RES2P,CT2P --qtymerge qt010y-11 --outdir CT10nnlo_RESCT2P_160512

#============================================
# MERGE GRID
#============================================
#./scripts/merge.sh --proc z0 --term REAL --indir  results_grid/group.perf-jets  --gridmerge v146.*_results_merge --outdir grid_fabrice_160524 RUN
#./scripts/merge.sh --proc z0 --term VIRT --indir  results_grid/group.perf-jets  --gridmerge v146.*_results_merge --outdir grid_fabrice_160524  RUN

#./scripts/merge.sh --proc wp,wm --term REAL --indir  results_grid/group.perf-jets  --gridmerge v146.*_results_merge --outdir grid_fabrice_160524
#./scripts/merge.sh --proc wp,wm --term VIRT --indir  results_grid/group.perf-jets  --gridmerge v146.*_results_merge --outdir grid_fabrice_160524  RUN

# ./scripts/merge.sh --proc wp,wm,z0 --term REAL,VIRT --indir  results_grid/group.perf-jets  --gridmerge v146.*_results_merge --outdir grid_fabrice_160524 RUN

# MMHT
#./scripts/merge.sh --proc wp,wm,z0 --term REAL,VIRT --indir  results_grid/group.perf-jets --pdfset MMHT2014nnlo68cl  --gridmerge v146.*_results_merge --outdir grid_fabrice_160524 RUN
#./scripts/merge.sh --proc wp --term VIRT --indir  results_grid/group.perf-jets --pdfset MMHT2014nnlo68cl  --gridmerge v146.*_results_merge --outdir grid_fabrice_160524 RUN
#./scripts/merge.sh --proc wm --term REAL --indir  results_grid/group.perf-jets --pdfset MMHT2014nnlo68cl  --gridmerge v146.*_results_merge --outdir grid_fabrice_160524 RUN
#./scripts/merge.sh --proc wp --term REAL --indir  results_grid/group.perf-jets --pdfset MMHT2014nnlo68cl  --gridmerge v146.*_results_merge --outdir grid_fabrice_160524 RUN


#CT14
#./scripts/merge.sh --proc wp,wm,z0 --term REAL,VIRT --indir  results_grid/group.perf-jets --pdfset CT14nnlo  --gridmerge v146.*_results_merge --outdir grid_fabrice_160524 RUN




# FIXED ORDER

./scripts/submit_DYTURBO.sh --grid --proc wp,wm,z0 --pdfset CT10nnlo --pdfvar all --term VV --seeds 500
./scripts/submit_DYTURBO.sh --grid --proc wp,wm,z0 --pdfset CT10nnlo --pdfvar array --term FIXCT2D 

# submit with niter=0 -- same seed
