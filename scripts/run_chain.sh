#!/bin/bash

thiswd=`pwd`
AS_LIST="0.1150 0.1170 0.1180 0.1182 0.1194 0.1200"
#GPAR_LIST="0.5 0.8 1.1 1.2 1.5"
GPAR_LIST="0.5 0.8 1.1 1.2 1.5"

rmmkcd(){
    rm -rf $1 && mkdir -p $1 && cd $1 || exit 2
}

prepare_lhapdf(){
    for as in $AS_LIST
    do
        asstring=`echo $as | sed "s|\.||"`
        newpdf=CT10nnlo_AS_${asstring}
        rmmkcd `lhapdf-config --prefix`/share/LHAPDF/$newpdf
        cat ../CT10nnlo/CT10nnlo.info |  sed "
           s|^AlphaS_MZ:.*$|AlphaS_MZ: ${as}|;
           s|^AlphaS_Type:.*$|AlphaS_Type: ode|;
           s|^AlphaS_Vals:.*$|AlphaS_Reference: ${as}|;
           s|^AlphaS_Qs:.*$|AlphaS_MassReference: 91.1876|;
           s|^NumMembers:.*$|NumMembers: 1|;" > ${newpdf}.info
        cp ../CT10nnlo/CT10nnlo_0000.dat ${newpdf}_0000.dat
    done
}

run_alphaS_scan(){
    for gparam in $GPAR_LIST
    do
        for as in $AS_LIST
        do
            asstring=`echo $as | sed "s|\.||"`
            output=$thiswd/results/CT10nnlo_AS_${asstring}_G_${gparam}.root
            $thiswd/bin/dyturbo $thiswd/input/alphaS_scan.in --gparam $gparam --pdfset CT10nnlo_AS_$asstring | tee out.log 
            cp results.root $output
        done
    done
}

submit_as_scan(){
    pdflist=""
    for as in $AS_LIST
    do
        asstring=`echo $as | sed "s|\.||"`
        pdflist=$pdflist"CT10nnlo_AS_$asstring,"
    done
    pdflist=${pdflist%,}
    basecommand="./scripts/submit_jobs.sh --mogon --infile input/alphaS_scan.in --proc z0 --mbins 1,66,116 --pdfset $pdflist --seeds 2 --yes"
    # finite order
    #$DRY_RUN $basecommand --order 1 --term CT
    #$DRY_RUN $basecommand --order 2 --term CT
    #$DRY_RUN $basecommand --order 1 --term VJLO
    #$DRY_RUN $basecommand --order 2 --term VJVIRT
    $DRY_RUN $basecommand --order 2 --term VJREAL
    #
    #for gparam in $GPAR_LIST
    #do
    #    $DRY_RUN $basecommand --order 1 --term BORN --gparam $gparam
    #    $DRY_RUN $basecommand --order 2 --term BORN --gparam $gparam
    #done
}

check_as_scan_out(){
    order=$1
    term=$2
    outfile=${dirbase}_${term}/${filebase}_o${order}t${term}_${filesuffix}
    ls $outfile  > /dev/null 2>&1 || filelist="$filelist
    $outfile"
    filelist="$filelist`ack -L "Successfully completed." $outfile 2> /dev/null`"
}

check_as_scan(){
    echo "Checking as scan output ... "
    filelist=""
    for as in $AS_LIST
    do
        asstring=`echo $as | sed "s|\.||"`
        #results_z0_CT10nnlo_AS_01200_VJVIRT/
        pdfset="CT10nnlo_AS_${asstring}"
        dirbase="results_z0_${pdfset}"
        #dyturbo_z0_CT10nnlo_AS_01200_0_o2tVJVIRT_seed_1.out
        filesuffix=seed_1.out
        filebase="dyturbo_z0_${pdfset}_0"
        check_as_scan_out 1 CT
        check_as_scan_out 1 VJLO
        check_as_scan_out 2 CT
        check_as_scan_out 2 VJREAL
        check_as_scan_out 2 VJVIRT
        for gparam in $GPAR_LIST
        do
            filebase="dyturbo_z0_CT10nnlo_AS_*_0_g${gparam}"
            check_as_scan_out 1 BORN
            check_as_scan_out 2 BORN
        done
    done
    [[ $filelist != "" ]] && echo -e "$filelist\n There were failed jobs" && exit 2
    echo "DONE"
}

merge_as_scan(){
    outdir=/etapfs03/atlashpc/cuth/LGNTuple/Base-2.4.20/LGNTupleMaker/LGAnalysis/Zpt13TeVAnalysis/share/templates
    order=2
    mkdir -p $outdir
    for gparam in $GPAR_LIST
    do
        for as in $AS_LIST
        do
            asstring=`echo $as | sed "s|\.||"`
            output=$outdir/CT10nnlo_AS_${asstring}_G_${gparam}.root
            inputs=""
            if [[ $order == 1 ]]
            then
                inputs+=" results_z0_CT10nnlo_AS_${asstring}_BORN/dyturbo_z0_CT10nnlo_AS_${asstring}_0_g${gparam}_o1t*_seed_1.root"
                inputs+=" results_z0_CT10nnlo_AS_${asstring}_CT/dyturbo_z0_CT10nnlo_AS_${asstring}_0_o1t*_seed_1.root"
                inputs+=" results_z0_CT10nnlo_AS_${asstring}_VJLO/dyturbo_z0_CT10nnlo_AS_${asstring}_0_o1t*_seed_1.root"
            fi
            if [[ $order == 2 ]]
            then
                inputs+=" results_z0_CT10nnlo_AS_${asstring}_BORN/dyturbo_z0_CT10nnlo_AS_${asstring}_0_g${gparam}_o2t*_seed_1.root"
                #inputs+=" results_z0_CT10nnlo_AS_${asstring}_CT/dyturbo_z0_CT10nnlo_AS_${asstring}_0_o2t*_seed_1.root"
                #inputs+=" results_z0_CT10nnlo_AS_${asstring}_VJREAL/dyturbo_z0_CT10nnlo_AS_${asstring}_0_o2t*_seed_1.root"
                #inputs+=" results_z0_CT10nnlo_AS_${asstring}_VJVIRT/dyturbo_z0_CT10nnlo_AS_${asstring}_0_o2t*_seed_1.root"
            fi
            # test input files
            ls $inputs > /dev/null || exit 1
            # merge
            $DRY_RUN hadd -f $output $inputs || exit 1
        done
    done
}

copy_to_accuracy(){
    rsync -avP results_z0_CT10* etap-accuracy.physik.uni-mainz.de:/home/cuth/work/unimainz/ATLAS/workarea/LGNtupleMaker/LGAnalysis/Zpt13TeVAnalysis/share/templates/DYTURBO/
}

DRY_RUN=
#DRY_RUN="echo "
main(){
    #prepare_lhapdf
    #run_alphaS_scan
    #submit_as_scan
    #check_as_scan
    #merge_as_scan
    copy_to_accuracy
}

main
exit 0

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

# submit v02 
# add version and change kinematics
# only nominal

#SUBM="./scripts/submit_jobs_wmass.sh --grid --griduser jcuth --voms phys-sm"
# central
#$SUBM --pdfset CT10nnlo         --version v02 --pdfvar 0 --infile input/wmass.in  --seeds 1000 --proc wp,wm --mbins 96,20,500 --order 2 --term CT,VV,REAL,VIRT
#$SUBM --pdfset CT14nnlo         --version v02 --pdfvar 0 --infile input/wmass.in  --seeds 1000 --proc wp,wm --mbins 96,20,500 --order 2 --term CT,VV,REAL,VIRT
#$SUBM --pdfset MMHT2014nnlo68cl --version v02 --pdfvar 0 --infile input/wmass.in  --seeds 1000 --proc wp,wm --mbins 96,20,500 --order 2 --term CT,VV,REAL,VIRT
#$SUBM --pdfset CT10nnlo         --version v02 --pdfvar 0 --infile input/wmass.in  --seeds 1000 --proc z0    --mbins 10,66,116 --order 2 --term CT,VV,REAL,VIRT
#$SUBM --pdfset CT14nnlo         --version v02 --pdfvar 0 --infile input/wmass.in  --seeds 1000 --proc z0    --mbins 10,66,116 --order 2 --term CT,VV,REAL,VIRT
#$SUBM --pdfset MMHT2014nnlo68cl --version v02 --pdfvar 0 --infile input/wmass.in  --seeds 1000 --proc z0    --mbins 10,66,116 --order 2 --term CT,VV,REAL,VIRT
# PDFvar
#$SUBM --pdfset CT10nnlo         --version v02 --pdfvar all --infile input/wmass_lowstat.in  --seeds 1000 --proc wp,wm --mbins 96,20,500 --order 2 --term CT,VV,REAL,VIRT
#$SUBM --pdfset CT14nnlo         --version v02 --pdfvar all --infile input/wmass_lowstat.in  --seeds 1000 --proc wp,wm --mbins 96,20,500 --order 2 --term CT,VV,REAL,VIRT
#$SUBM --pdfset MMHT2014nnlo68cl --version v02 --pdfvar all --infile input/wmass_lowstat.in  --seeds 1000 --proc wp,wm --mbins 96,20,500 --order 2 --term CT,VV,REAL,VIRT
#$SUBM --pdfset CT10nnlo         --version v02 --pdfvar all --infile input/wmass_lowstat.in  --seeds 1000 --proc z0    --mbins 10,66,116 --order 2 --term CT,VV,REAL,VIRT
#$SUBM --pdfset CT14nnlo         --version v02 --pdfvar all --infile input/wmass_lowstat.in  --seeds 1000 --proc z0    --mbins 10,66,116 --order 2 --term CT,VV,REAL,VIRT
#$SUBM --pdfset MMHT2014nnlo68cl --version v02 --pdfvar all --infile input/wmass_lowstat.in  --seeds 1000 --proc z0    --mbins 10,66,116 --order 2 --term CT,VV,REAL,VIRT


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

# for job in $samples
# do
# 
#     outlier=
#     [[ $job =~ REAL ]] && outlier=o
#     echo $job
#     #rucio ls group.phys-sm.$job
#     #echo
#     #rucio download --dir results_grid group.phys-sm.$job
#     echo
#     ls results_grid/group.phys-sm.$job/*.root* | wc -l
#     #ls results_merge/v01/${job} 
#     #/usr/bin/time -v ./bin/dyturbo-merger -d results_merge/v01/average_${job} results_grid/group.phys-sm.$job/group.phys-sm.*.root 2>&1
#     /usr/bin/time -v ./bin/dyturbo-merger -Td$outlier results_merge/v01/${job} results_grid/group.phys-sm.$job/group.phys-sm.*.root 2>&1
#     #/usr/bin/time -v ./../DYTURBO_CLIdev/bin/dyturbo-merger -d$outlier results_merge/v01/old_${job} results_grid/group.phys-sm.$job/group.phys-sm.*.root 2>&1
#     echo
#     echo
# done | tee merge-`hostname`.log

#job=dyturbo_z0_lhc7_CT14nnlo_all_o2tREAL_seed_v1466990840_results_merge.root
#/usr/bin/time -v ./bin/dyturbo-merger -d  results_merge/v01/average_$job results_grid/group.phys-sm.$job/*root*
#/usr/bin/time -v ./bin/dyturbo-merger -do results_merge/v01/$job         results_grid/group.phys-sm.$job/*root*
