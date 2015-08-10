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

prepare_script(){
    mkdir -p $batch_script_dir
    mkdir -p $result_dir
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
    order=1
    if [[ $sample == CT10nlo  ]]; then order=1; fi;
    if [[ $sample == CT10nnlo ]]; then order=2; fi;

    cp $dyturbo_in_tmpl tmp
    sed -i "s|LOQTBIN|$lobin|g " tmp
    sed -i "s|HIQTBIN|$hibin|g " tmp
    sed -i "s|SAMPLENAME|$sample|g " tmp
    sed -i "s|SETORDER|$order|g " tmp
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
    sample=CT10nnlo
    random_seed=100101
    lobin=unset
    hibin=unset
    for sample in CT10nnlo CT10nlo
    do
        for ibin in 0 2 4 6 8 10 12 14 16 18 22 26 30 34 38 42 46 50 54 60 70 80 100 150 200 300 800
        do
            if [[ $lobin == unset ]]
            then
                lobin=$ibin
                continue
            fi
            hibin=$ibin
            qtregion=qt${lobin}${hibin}
            prepare_script
            $DRYRUN submit_job
            lobin=$ibin
        done
    done
    $DRYRUN jobsDone
}



# MAIN
submit_Z_dyturbo
