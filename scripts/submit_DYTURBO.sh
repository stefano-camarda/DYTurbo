#!/bin/bash

##
## @file submit_resbos.sh
##
## Description of the script file
##
## @author cuto <Jakub.Cuth@cern.ch>
## @date 2014-12-02

DRYRUN=echo 
DRYRUN=


# define directories
dyturbo_project=/etapfs03/atlashpc/cuth/DYTURBO
batch_script_dir=$dyturbo_project/scripts/batch_scripts
in_files_dir=$dyturbo_project/scripts/infiles
result_dir=$dyturbo_project/results
# templates
batch_template=$dyturbo_project/scripts/run_DYTURBO_TMPL.sh
dyturbo_in_tmpl=$dyturbo_project/scripts/DYTURBO_TMPL.in


job_name=unset
in_file=unset
sh_file=unset
process=unset
collider=unset
random_seed=unset
lobin=unset
hibin=unset
sample=unset
member=unset

prepare_script(){
    mkdir -p $batch_script_dir
    mkdir -p $result_dir
    sample=${pdfset}_${variation}
    job_name=dyturbo_${process}_${collider}_${sample}_${qtregion}_${random_seed}
    sh_file=$batch_script_dir/$job_name.sh
    echo $job_name
    # prepare in
    prepare_in
    # prepare script
    cp $batch_template tmp
    sed -i "s|JOBNAME|$job_name|g               " tmp
    sed -i "s|OUTDIR|$result_dir|g              " tmp
    sed -i "s|DYTURBOROOTDIR|$dyturbo_project|g " tmp
    sed -i "s|DYTURBOINPUTFILE|$in_file|g       " tmp
    mv tmp $sh_file
    chmod +x $sh_file
}

prepare_in(){
    mkdir -p $in_files_dir
    in_file=$in_files_dir/$job_name.in
    # order check
    order=1
    if [[ $pdfset == CT10nlo  ]]; then order=1; fi;
    if [[ $pdfset == CT10nnlo ]]; then order=2; fi;
    # variations
    gpar=1e0
    member=0
    re='^[0-9]+$'
    pdfsetname=$pdfset
    if [[ $variation == g_05  ]]; then member=0; gpar=0.5e0;                        fi;
    if [[ $variation == g_15  ]]; then member=0; gpar=1.5e0;                        fi;
    if [[ $variation == as_*  ]]; then member=0; pdfsetname=${pdfset}_${variation}; fi;
    if [[ $variation =~ $re   ]]; then member=$variation;                           fi;
    cp $dyturbo_in_tmpl tmp
    sed -i "s|LOQTBIN|$lobin|g        " tmp
    sed -i "s|HIQTBIN|$hibin|g        " tmp
    sed -i "s|SETPDFSET|$pdfsetname|g " tmp
    sed -i "s|SETMEMBER|$member|g     " tmp
    sed -i "s|SETORDER|$order|g       " tmp
    sed -i "s|SETGPAR|$gpar|g         " tmp
    mv tmp $in_file

}

submit_job(){
    bsub < $sh_file
}

jobsDone(){
    bsub < scripts/send_mail.sh
}

submit_Z_dyturbo(){
    process=z0
    pdfset=CT10nnlo
    random_seed=100101
    lobin=unset
    hibin=unset
    collider=lhc7
    for ibin in 0 2 4 6 8 10 12 14 16 18 22 26 30 34 38 42 46 50 54 60 70 80 100 150 200 300 800
    do
        if [[ $lobin == unset ]]
        then
            lobin=$ibin
            continue
        fi
        hibin=$ibin
        qtregion=qt${lobin}${hibin}
        for variation in `seq 0 50` g_05 g_15 as_0117 as_0119
        do
            prepare_script
            $DRYRUN submit_job
        done
        lobin=$ibin
    done
    $DRYRUN jobsDone
}

clear_files(){
    rm -rf scripts/batch_scripts/*.sh
    rm -rf scripts/infiles/*.in
}

# MAIN
clear_files
submit_Z_dyturbo
