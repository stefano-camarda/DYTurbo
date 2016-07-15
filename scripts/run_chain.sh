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
#DONE ./scripts/submit_DYTURBO.sh --mogon --proc z0 --term RES2D,CT2D --pdfset MMHT2014nnlo68cl --pdfvar array RUN
#DONE ./scripts/merge.sh --proc z0 --term RES2D,CT2D --pdfset MMHT2014nnlo68cl --qtymerge qt05y-11 --outdir MMHT14_RESCT2D_160523

# DONE ./scripts/submit_DYTURBO.sh --mogon --proc wp,wm --term RES2D,CT2D --pdfset MMHT2014nnlo68cl --pdfvar array RUN
# DONE./scripts/merge.sh --proc wp,wm --term RES2D,CT2D --pdfset MMHT2014nnlo68cl --qtymerge qt05y-11 --outdir MMHT14_RESCT2D_160523 RUN

## CT14NNLO
#DONE ./scripts/submit_DYTURBO.sh --mogon --proc wm,wp,z0 --term RES2D,CT2D --pdfset CT14nnlo --pdfvar array RUN
#DONE ./scripts/merge.sh --proc wp,wm,z0 --term RES2D,CT2D --pdfset CT14nnlo --qtymerge qt05y-11 --outdir CT14_RESCT2D_160523 RUN

#./scripts/merge.sh --proc z0 --term RES2D,CT2D --pdfset CT14nnlo --qtymerge qt05y-11 --outdir CT14nnlo_RESCT2D_160523 RUN
#./scripts/submit_DYTURBO.sh --mogon --proc wp,wm --term RES2D,CT2D --pdfset CT14nnlo --pdfvar array RUN
#./scripts/merge.sh --proc wp,wm,z0 --term RES2D,CT2D --pdfset CT14nnlo --qtymerge qt05y-11 --outdir CT14nnlo_RESCT2D_160523 #RUN

#RERUN:
#./scripts/submit_DYTURBO.sh --mogon --proc z0 --term RES2P,CT2P --pdfset CT10nnlo --pdfvar array RUN
#./scripts/merge.sh --proc z0 --term RES2P,CT2P --qtymerge qt05y-11 --outdir CT10nnlo_RESCT2P_160523 RUN

#./scripts/submit_DYTURBO.sh --mogon --proc wp,wm --term RES2P,CT2P --pdfset CT10nnlo --pdfvar array
#./scripts/merge.sh --proc wp,wm --term RES2P,CT2P --qtymerge qt010y-11 --outdir CT10nnlo_RESCT2P_160512

#TEST LO
#./scripts/submit_DYTURBO.sh --mogon --proc z0 --pdfset CT10nnlo --pdfvar all --term LO --order 1 --seeds 100
#./scripts/submit_DYTURBO.sh --mogon --proc z0 --pdfset CT10nnlo --pdfvar 0 --term REAL --seeds 100 --qtlo 30 --qthi 90



##  -- PROFILED PDFS

## PROFILLED MMHT
#./scripts/submit_DYTURBO.sh --lxbatch --proc z0 --term RES2D,CT2D --pdfset MMHTProf68cl --pdfvar array RUN
#./scripts/merge.sh --proc z0 --term RES2D,CT2D --pdfset MMHTProf68cl --qtymerge qt05y-11 --outdir MMHT14Prof_RESCT2D_160523

#./scripts/submit_DYTURBO.sh --lxbatch --proc wp,wm --term RES2D,CT2D --pdfset MMHTProf68cl --pdfvar array RUN
#./scripts/merge.sh --proc wp,wm --term RES2D,CT2D --pdfset MMHTProf68cl --qtymerge qt05y-11 --outdir MMHT14Prof_RESCT2D_160523

# rerun
#./scripts/submit_DYTURBO.sh --lxbatch --proc z0 --term CT2D  --pdfset MMHTProf68cl --pdfvar array RUN
#./scripts/submit_DYTURBO.sh --lxbatch --proc wp --term RES2D --pdfset MMHTProf68cl --pdfvar array RUN

#./scripts/merge.sh --proc wp,wm,z0 --term RES2D,CT2D --pdfset MMHTProf68cl --qtymerge qt05y-11 --outdir MMHT14Prof_RESCT2D_160523 #RUN

# PROFILLED CT10
#./scripts/submit_DYTURBO.sh --mogon --proc wp,wm,z0 --term RES2D,CT2D  --pdfset CT10nnlo68clProfiled --pdfvar array RUN
#./scripts/submit_DYTURBO.sh --mogon --proc wp --term RES2D,CT2D  --pdfset CT10nnlo68clProfiled --pdfvar array RUN
#./scripts/submit_DYTURBO.sh --mogon --proc wm --term RES2D       --pdfset CT10nnlo68clProfiled --pdfvar array RUN

#./scripts/merge.sh --proc wp,wm,z0 --term RES2D,CT2D --pdfset CT10nnlo68clProfiled --qtymerge qt05y-11 --outdir CT10Prof_RESCT2D_160523
#./scripts/merge.sh --proc z0 --term RES2D,CT2D --pdfset CT10nnlo68clProfiled --qtymerge qt05y-11 --outdir CT10Prof_RESCT2D_160523 RUN
#./scripts/merge.sh --proc wm --term CT2D --pdfset CT10nnlo68clProfiled --qtymerge qt05y-11 --outdir CT10Prof_RESCT2D_160523 RUN
#./scripts/merge.sh --proc wp --term RES2D,CT2D --pdfset CT10nnlo68clProfiled --qtymerge qt05y-11 --outdir CT10Prof_RESCT2D_160523 RUN
#./scripts/merge.sh --proc wm --term RES2D --pdfset CT10nnlo68clProfiled --qtymerge qt05y-11 --outdir CT10Prof_RESCT2D_160523 RUN



# PROFILLED CT14
#./scripts/submit_DYTURBO.sh --mogon --proc wp,wm,z0 --term RES2D,CT2D  --pdfset CT14nnloProf68cl --pdfvar array #RUN
#./scripts/merge.sh --proc wp,wm,z0 --term RES2D,CT2D --pdfset CT14nnloProf68cl --qtymerge qt05y-11 --outdir CT14Prof_RESCT2D_160523 RUN

#============================================
# MERGE GRID DONE
#============================================
#./scripts/merge.sh --proc z0 --term REAL --indir  results_grid/group.perf-jets  --gridmerge v146.*_results_merge --outdir grid_fabrice_160524 RUN
#./scripts/merge.sh --proc z0 --term VIRT --indir  results_grid/group.perf-jets  --gridmerge v146.*_results_merge --outdir grid_fabrice_160524  RUN

#./scripts/merge.sh --proc wp,wm --term REAL --indir  results_grid/group.perf-jets  --gridmerge v146.*_results_merge --outdir grid_fabrice_160524
#./scripts/merge.sh --proc wp,wm --term VIRT --indir  results_grid/group.perf-jets  --gridmerge v146.*_results_merge --outdir grid_fabrice_160524  RUN

# ./scripts/merge.sh --proc wp,wm,z0 --term REAL,VIRT --indir  results_grid/group.perf-jets  --gridmerge v146.*_results_merge --outdir grid_fabrice_160524 RUN

# MMHT
#./scripts/merge.sh --proc wp,wm,z0 --term REAL,VIRT --indir  results_grid/group.perf-jets --pdfset MMHT2014nnlo68cl  --gridmerge v146.*_results_merge --outdir grid_fabrice_160524 #RUN

#./scripts/merge.sh --proc wp --term VIRT --indir  results_grid/group.perf-jets --pdfset MMHT2014nnlo68cl  --gridmerge v146.*_results_merge --outdir grid_fabrice_160524 RUN
#./scripts/merge.sh --proc wm --term REAL --indir  results_grid/group.perf-jets --pdfset MMHT2014nnlo68cl  --gridmerge v146.*_results_merge --outdir grid_fabrice_160524 RUN
#./scripts/merge.sh --proc wp --term REAL --indir  results_grid/group.perf-jets --pdfset MMHT2014nnlo68cl  --gridmerge v146.*_results_merge --outdir grid_fabrice_160524 RUN


#CT14
#./scripts/merge.sh --proc wp,wm,z0 --term REAL,VIRT --indir  results_grid/group.perf-jets --pdfset CT14nnlo  --gridmerge v146.*_results_merge --outdir grid_fabrice_160524

# download profilled
#rucio ls group.perf-jets.dyturbo_*_lhc7_*_all_o2qt050y-55t*_seed_v1463677142_00_results_merge.root

#rucio ls group.perf-jets.dyturbo_*_lhc7_MMHTProf68cl_all_o2qt050y-55t*_seed_v1463677142_00_results_merge.root
#./scripts/merge.sh --proc wp,wm,z0 --term REAL,VIRT --indir  results_grid/group.perf-jets --pdfset MMHTProf68cl  --gridmerge v146.*_results_merge --outdir grid_fabrice_160524 RUN

## REDOWNLOAD REAL

#rucio download group.phys-sm.dyturbo_*_lhc7_*_all_o2qt050y-55tREAL_seed_v1463677145_results_merge.root

#./scripts/merge.sh --proc wp,wm,z0 --term REAL --indir results_grid/group.phys-sm --pdfset CT10nnlo,CT10nnlo68clProfiled,MMHT2014nnlo68cl,MMHTProf68cl,CT14nnlo,CT14nnloProf68cl --gridmerge v146.*_results_merge --outdir grid_fabrice_160603

# ATLAS2
#./scripts/merge.sh --proc wp,wm,z0 --term REAL --indir results_grid/group.phys-sm --pdfset CT10nnlo             --gridmerge v146.*_results_merge --outdir grid_fabrice_160604 RUN
#./scripts/merge.sh --proc wp,wm,z0 --term REAL --indir results_grid/group.phys-sm --pdfset CT10nnlo68clProfiled --gridmerge v146.*_results_merge --outdir grid_fabrice_160604 RUN
# ATLAS3
#./scripts/merge.sh --proc wp,wm,z0 --term REAL --indir results_grid/group.phys-sm --pdfset MMHT2014nnlo68cl     --gridmerge v146.*_results_merge --outdir grid_fabrice_160604 RUN
#./scripts/merge.sh --proc wp,wm,z0 --term REAL --indir results_grid/group.phys-sm --pdfset MMHTProf68cl         --gridmerge v146.*_results_merge --outdir grid_fabrice_160604 RUN
# ATLAS4
#./scripts/merge.sh --proc wp,wm,z0 --term REAL --indir results_grid/group.phys-sm --pdfset CT14nnlo             --gridmerge v146.*_results_merge --outdir grid_fabrice_160604 RUN
#./scripts/merge.sh --proc wp,wm,z0 --term REAL --indir results_grid/group.phys-sm --pdfset CT14nnloProf68cl     --gridmerge v146.*_results_merge --outdir grid_fabrice_160604 RUN


#============================================
# FIXED ORDER
#============================================
#./scripts/submit_DYTURBO.sh --lxbatch --proc wp,wm,z0 --pdfset CT10nnlo --pdfvar array --term FIXCT2D RUN
#./scripts/submit_DYTURBO.sh --lxbatch --proc wp,wm,z0 --pdfset CT14nnlo --pdfvar array --term FIXCT2D RUN
#./scripts/submit_DYTURBO.sh --lxbatch --proc wp,wm,z0 --pdfset MMHT2014nnlo68cl --pdfvar array --term FIXCT2D RUN

# PROFILLED
#./scripts/submit_DYTURBO.sh --lxbatch --proc wp,wm,z0 --pdfset CT10nnlo68clProfiled --pdfvar array --term FIXCT2D RUN
#./scripts/submit_DYTURBO.sh --lxbatch --proc wp,wm,z0 --pdfset CT14nnloProf68cl --pdfvar array --term FIXCT2D RUN
#./scripts/submit_DYTURBO.sh --lxbatch --proc wp,wm,z0 --pdfset MMHTProf68cl --pdfvar array --term FIXCT2D RUN

#SKIP./scripts/submit_DYTURBO.sh --grid --proc wp,wm,z0 --pdfset CT10nnlo --pdfvar all --term VV --seeds 500
# SKIP ./scripts/submit_DYTURBO.sh --mogon --proc wp,wm,z0 --pdfset CT10nnlo --pdfvar array --term FIXCT2D  RUN
# SKIP ./scripts/submit_DYTURBO.sh --mogon --proc wp,wm,z0 --pdfset CT14nnlo --pdfvar array --term FIXCT2D  RUN
# SKIP ./scripts/submit_DYTURBO.sh --mogon --proc wp,wm,z0 --pdfset MMHT2014nnlo68cl --pdfvar array --term FIXCT2D  RUN

#  
#  gitshashort=`git rev-parse --short HEAD`
#  outGridName=GRID_FIXCT_CT10_$gitshashort
#  
#  
#  ./scripts/submit_DYTURBO.sh --grid --proc wp,wm,z0 --pdfset CT10nnlo --pdfvar array --term FIXCT2D
#  outGridName=GRID_FIXCT_CT10_$gitshashort
#  rm -rf $outGridName
#  mv GRID $outGridName && tar czvf ${outGridName}.tgz $outGridName
#  
#  ./scripts/submit_DYTURBO.sh --grid --proc wp,wm,z0 --pdfset CT14nnlo --pdfvar array --term FIXCT2D
#  outGridName=GRID_FIXCT_CT14_$gitshashort
#  rm -rf $outGridName
#  mv GRID $outGridName && tar czvf ${outGridName}.tgz $outGridName
#  
#  ./scripts/submit_DYTURBO.sh --grid --proc wp,wm,z0 --pdfset MMHT2014nnlo68cl --pdfvar array --term FIXCT2D
#  outGridName=GRID_FIXCT_MMHT14_$gitshashort
#  rm -rf $outGridName
#  mv GRID $outGridName && tar czvf ${outGridName}.tgz $outGridName

# submit with niter=0 -- same seed



#============================================
## NEW RUN v01
#============================================

done_samples_019="
dyturbo_z0_lhc7_CT10nnlo_all_o2tFIXCT_seed_v1466990784_results_merge.root
dyturbo_z0_lhc7_CT10nnlo_all_o2tREAL_seed_v1466990784_results_merge.root
dyturbo_z0_lhc7_CT10nnlo_all_o2tVIRT_seed_v1466990784_results_merge.root
dyturbo_z0_lhc7_CT10nnlo_all_o2tVV_seed_v1466990784_results_merge.root

dyturbo_wp_lhc7_CT10nnlo_all_o2tFIXCT_seed_v1466990784_results_merge.root
dyturbo_wp_lhc7_CT10nnlo_all_o2tVV_seed_v1466990784_results_merge.root

dyturbo_wm_lhc7_CT10nnlo_all_o2tFIXCT_seed_v1466990784_results_merge.root
dyturbo_wm_lhc7_CT10nnlo_all_o2tVV_seed_v1466990784_results_merge.root
"
samples_019="
dyturbo_wp_lhc7_CT10nnlo_all_o2tREAL_seed_v1466990784_results_merge.root
dyturbo_wp_lhc7_CT10nnlo_all_o2tVIRT_seed_v1466990784_results_merge.root
dyturbo_wm_lhc7_CT10nnlo_all_o2tREAL_seed_v1466990784_results_merge.root
dyturbo_wm_lhc7_CT10nnlo_all_o2tVIRT_seed_v1466990784_results_merge.root
"

done_samples_006="
dyturbo_wp_lhc7_CT14nnlo_all_o2tREAL_seed_v1466990840_results_merge.root
dyturbo_wm_lhc7_CT14nnlo_all_o2tREAL_seed_v1466990840_results_merge.root
dyturbo_z0_lhc7_CT14nnlo_all_o2tREAL_seed_v1466990840_results_merge.root

dyturbo_wp_lhc7_CT14nnlo_all_o2tFIXCT_seed_v1466990840_results_merge.root
dyturbo_wp_lhc7_CT14nnlo_all_o2tVV_seed_v1466990840_results_merge.root
dyturbo_wp_lhc7_CT14nnlo_all_o2tVIRT_seed_v1466990840_results_merge.root

dyturbo_wm_lhc7_CT14nnlo_all_o2tFIXCT_seed_v1466990840_results_merge.root
dyturbo_wm_lhc7_CT14nnlo_all_o2tVV_seed_v1466990840_results_merge.root
dyturbo_wm_lhc7_CT14nnlo_all_o2tVIRT_seed_v1466990840_results_merge.root

dyturbo_z0_lhc7_CT14nnlo_all_o2tFIXCT_seed_v1466990840_results_merge.root
dyturbo_z0_lhc7_CT14nnlo_all_o2tVV_seed_v1466990840_results_merge.root
dyturbo_z0_lhc7_CT14nnlo_all_o2tVIRT_seed_v1466990840_results_merge.root
"

done_samples_078="
dyturbo_wm_lhc7_MMHT2014nnlo68cl_all_o2tFIXCT_seed_v1466990932_results_merge.root
dyturbo_wm_lhc7_MMHT2014nnlo68cl_all_o2tVV_seed_v1466990932_results_merge.root
dyturbo_wm_lhc7_MMHT2014nnlo68cl_all_o2tVIRT_seed_v1466990932_results_merge.root

dyturbo_wp_lhc7_MMHT2014nnlo68cl_all_o2tFIXCT_seed_v1466990932_results_merge.root
dyturbo_wp_lhc7_MMHT2014nnlo68cl_all_o2tVV_seed_v1466990932_results_merge.root
dyturbo_wp_lhc7_MMHT2014nnlo68cl_all_o2tVIRT_seed_v1466990932_results_merge.root

dyturbo_z0_lhc7_MMHT2014nnlo68cl_all_o2tFIXCT_seed_v1466990932_results_merge.root
dyturbo_z0_lhc7_MMHT2014nnlo68cl_all_o2tVV_seed_v1466990932_results_merge.root
dyturbo_z0_lhc7_MMHT2014nnlo68cl_all_o2tREAL_seed_v1466990932_results_merge.root
dyturbo_z0_lhc7_MMHT2014nnlo68cl_all_o2tVIRT_seed_v1466990932_results_merge.root
"

samples_078="
dyturbo_wm_lhc7_MMHT2014nnlo68cl_all_o2tREAL_seed_v1466990932_results_merge.root
dyturbo_wp_lhc7_MMHT2014nnlo68cl_all_o2tREAL_seed_v1466990932_results_merge.root
"

real="
dyturbo_wp_lhc7_CT10nnlo_all_o2tREAL_seed_v1466990784_results_merge.root
dyturbo_wm_lhc7_CT10nnlo_all_o2tREAL_seed_v1466990784_results_merge.root
dyturbo_z0_lhc7_CT10nnlo_all_o2tREAL_seed_v1466990784_results_merge.root

dyturbo_wp_lhc7_CT14nnlo_all_o2tREAL_seed_v1466990840_results_merge.root
dyturbo_wm_lhc7_CT14nnlo_all_o2tREAL_seed_v1466990840_results_merge.root
dyturbo_z0_lhc7_CT14nnlo_all_o2tREAL_seed_v1466990840_results_merge.root

dyturbo_wp_lhc7_MMHT2014nnlo68cl_all_o2tREAL_seed_v1466990932_results_merge.root
dyturbo_wm_lhc7_MMHT2014nnlo68cl_all_o2tREAL_seed_v1466990932_results_merge.root
dyturbo_z0_lhc7_MMHT2014nnlo68cl_all_o2tREAL_seed_v1466990932_results_merge.root
"

test_real="dyturbo_wp_lhc7_CT10nnlo_all_o2tREAL_seed_v1466990784_results_merge.root
dyturbo_wm_lhc7_CT10nnlo_all_o2tREAL_seed_v1466990784_results_merge.root
"

test_real14="
dyturbo_z0_lhc7_CT14nnlo_all_o2tREAL_seed_v1466990840_results_merge.root
"

#samples=$samples_019
#samples=$done_samples_006
#samples=$samples_078
#samples=$test_real
samples=$test_real14

for job in $samples
do

    outlier=
    [[ $job =~ REAL ]] && outlier=o
    echo $job
    #rucio ls group.phys-sm.$job
    #echo
    #rucio download --dir results_grid group.phys-sm.$job
    echo
    ls results_grid/group.phys-sm.$job/*.root* | wc -l
    #ls results_merge/v01/${job} 
    #/usr/bin/time -v ./bin/dyturbo-merger -d results_merge/v01/average_${job} results_grid/group.phys-sm.$job/group.phys-sm.*.root 2>&1
    /usr/bin/time -v ./bin/dyturbo-merger -Td$outlier results_merge/v01/${job} results_grid/group.phys-sm.$job/group.phys-sm.*.root 2>&1
    #/usr/bin/time -v ./../DYTURBO_CLIdev/bin/dyturbo-merger -d$outlier results_merge/v01/old_${job} results_grid/group.phys-sm.$job/group.phys-sm.*.root 2>&1
    echo
    echo
done | tee merge-`hostname`.log
