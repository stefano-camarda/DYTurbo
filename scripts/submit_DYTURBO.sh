#!/bin/bash

##
## @file submit_resbos.sh
##
## Description of the script file
##
## @author cuto <Jakub.Cuth@cern.ch>
## @date 2014-12-02

DRYRUN=echo 
#DRYRUN=


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
loqtbin=unset
hiqtbin=unset
loybin=unset
hiybin=unset
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
    if [[ $variation == g_05  ]]; then gpar=0.5e0;                        fi;
    if [[ $variation == g_15  ]]; then gpar=1.5e0;                        fi;
    if [[ $variation == as_*  ]]; then pdfsetname=${pdfset}_${variation}; fi;
    if [[ $variation =~ $re   ]]; then member=$variation;                 fi;
    cp $dyturbo_in_tmpl tmp
    sed -i "s|LOQTBIN|$loqtbin|g        " tmp
    sed -i "s|HIQTBIN|$hiqtbin|g        " tmp
    sed -i "s|LOYBIN|$loybin|g        " tmp
    sed -i "s|HIYBIN|$hiybin|g        " tmp
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
    DRYRUN=echo 
    read -p "Do you want to submit results ? " -n 1 -r
    echo    # (optional) move to a new line
    if [[ $REPLY =~ ^[Yy]$ ]]
    then
        DRYRUN=
    fi
    process=z0
    pdfset=CT10nnlo
    random_seed=100101
    loqtbin=unset
    hiqtbin=unset
    collider=lhc7
    for iqtbin in 0 2 4 6 8 10 12 14 16 18 22 26 30 34 38 42 46 50 54 60 70 80 100 150 200 300 800
    do
        if [[ $loqtbin == unset ]]; then loqtbin=$iqtbin; continue; fi; hiqtbin=$iqtbin
        for iybin in 0 1 2 2.4
        do
            if [[ $loybin == unset ]]; then loybin=$iybin; continue; fi; hiybin=$iybin
            qtregion=qt${loqtbin}${hiqtbin}y${loybin}${hiybin}
            for variation in `seq 0 50` g_05 g_15 as_0117 as_0119
            do
                prepare_script
                $DRYRUN submit_job
            done
            loybin=$iybin
        done
        loybin=unset
        hiybin=unset
        loqtbin=$iqtbin
    done
    $DRYRUN jobsDone
}

clear_files(){
    echo "Clearing setup files"
    rm -f scripts/batch_scripts/*.sh
    rm -f scripts/infiles/*.in
    echo "Done"
}

clear_results(){
    read -p "Are you sure you want to delete all current results ? " -n 1 -r
    echo    # (optional) move to a new line
    if [[ $REPLY =~ ^[Yy]$ ]]
    then
        rm -f results/*.{root,out,err}
    fi
    echo "Done"
}

# MAIN
clear_files
#clear_results
submit_Z_dyturbo
