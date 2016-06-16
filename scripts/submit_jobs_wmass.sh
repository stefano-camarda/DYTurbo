#!/bin/bash

##
## @file submit_jobs_wmass.sh
##
## Description of the script file
##
## @author Jakub Cuth <Jakub.Cuth@cern.ch>
## @date 2016-06-16

help(){
    echo " Help is on TODO list :), But you can read parse_input function to check possibilities"
}

main(){
    #parse_input $*
    parse_input --lxbatch --proc z0 --pdfset CT10nnlo,CT10nnlo68clProfiled --pdfvar all
    clear_files

    submit_jobs_wmass
    if [[ $batchsystem =~ lxplus|mogon|localrun  ]] 
    then
        ask_submit
        $DRYRUN submit_jobs_wmass
    fi

}

parse_input(){
    batchsystem=unset
    proclist=z0
    #
    order=1
    termlist="LO VV FIXCT"
    order=2
    termlist="REAL VIRT FIXCT VV"
    #
    pdflist="CT10nnlo CT10nnlo68clProfiled"
    pdfvarlist=all
    pdfvarlist="0 1 2 3"
    pdfvarlist=0
    #
    infile=infile/wmass.in
    #
    seedlist=100,201 # run 100 and 201
    seedlist=100-201 # run form 100 to 201
    seedlist=800    # run from 1 to 800

    while [[ $# > 0 ]]
    do
        key=$1
        #echo debug key $key
        case $key in
            # TARGET
            --local)
                target=localrun
                ;;
            --grid)
                target=grid
                ;;
            --mogon)
                target=mogon
                ;;
            --lxbatch)
                target=lxbatch
                ;;
            --proc)
                proclist="`echo $2|sed 's|,| |g'`"
                shift
                ;;
            --order)
                order=$2
                shift
                ;;
            --term)
                termlist="`echo $2|sed 's|,| |g'`"
                shift
                ;;
            --pdfset)
                pdflist="`echo $2|sed 's|,| |g'`"
                shift
                ;;
            --pdfvar)
                pdfvarlist="`echo $2|sed 's|,| |g'`"
                shift
                ;;
            --infile)
                infile=$2
                shift
                ;;
            --seeds)
                seedlist=$2
                shift
                ;;
            #  HELP and OTHER
            -h|--help)
                help
                exit 0
                ;;
            *)
                echo Uknown option: $key
                help
                exit 3
                ;;
        esac
        shift
    done
}

submit_jobs_wmass(){
    program=dyturbo
    random_seed=seed
    # check seed
    [[ $target =~ grid ]] || [[ $seedlist =~ - ]]  || seedlist=1-$seedlist
    # check order term
    [[ $order == 1 ]] && [[ $termlist =~ REAL|VIRT ]] && echo "WRONG ORDER $order TO TERM $termlist" && return 3
    [[ $order == 2 ]] && [[ $termlist =~ LO        ]] && echo "WRONG ORDER $order TO TERM $termlist" && return 3 
    # loops, lopps and loops
    for pdfset in pdflist 
    do
        if [[ $target =~ grid ]]
        then
            cubacores=0
            prepare_tarbal
        fi
        for process in $proclist
        do
            for terms in $termlist
            do
                for variation in $pdfvarlist
                do
                    qtregion=`echo o${order}t${terms} | sed "s/\.//g;s/ //g"`
                    prepare_script
                    [[ $target =~ mogon|lxbatch|localrun ]] &&  submit_job
                    [[ $target =~ grid ]]  && add_to_tarbal && submit_job2grid
                done # variation
            done #term
        done #proc
        # finish tarbal
    done #pdfsets
}

#directories 
dyturbo_project=`pwd -P`
batch_script_dir=$dyturbo_project/scripts/batch_scripts
queue=atlasshort
batch_template=$dyturbo_project/scripts/run_on_batch_tmpl.sh

# this goes away
in_files_dir=unset
cubacores=unset
job_name=unset
in_file=unset
sh_file=unset
process=unset
collider=unset
random_seed=unset
loqtbin=unset
terms=ALL
hiqtbin=unset
loybin=unset
hiybin=unset
sample=unset
member=unset
program=unset
fiducial=unset
detfiducial=0
makelepcuts=false
gpar=1
tarbalfile=unset
seedlist=unset

prepare_script(){
    # temporary add lhapdf to path
    PATH=$PATH:lhapdf6/bin
    if ! lhapdf-config --prefix > /dev/null
    then
        echo "Cannot reach lhapdf-config. Please setup the PATH to LHAPDF"
        exit 3
    fi
    result_dir=$dyturbo_project/results_${process}_${pdfset}_${terms}
    mkdir -p $batch_script_dir
    mkdir -p $result_dir
    #
    sample=${pdfset}_${variation}
    job_name=${program}_${process}_${collider}_${sample}_${qtregion}_${random_seed}
    sh_file=$batch_script_dir/$job_name.sh
    echo $job_name $seedlist
    #
    mlo=50 && mhi=1000
    [[ $proces =~ z0 ]] && mlo=66 && mhi=116
    # argumets
    argumets= "--proc $process \
               --mtbins 100,$mlo,$mhi
               --pdfset $pdfset \
               --pdfvar $variation \
               --order $order \
               --term $terms \
               $infile
    "
    # job queue
    nprocessors=1
    walltime=5:00
    if [[ $program =~ dyres ]] || [[ $variation =~ all ]] || [[ $terms =~ [23]P ]]
    then
        walltime=20:00
        queue=atlaslong
    fi
    [[ $target =~ lxbatch ]] && queue=8nh && walltime=8:00
    #a prepare script
    cp $batch_template tmp
    sed -i "s|JOBNAME|$job_name|g               " tmp
    sed -i "s|SEEDLIST|$seedlist|g              " tmp
    sed -i "s|OUTDIR|$result_dir|g              " tmp
    sed -i "s|SETQUEUE|$queue|g                 " tmp
    sed -i "s|SETWALLTIME|$walltime|g           " tmp
    sed -i "s|DYTURBOROOTDIR|$dyturbo_project|g " tmp
    sed -i "s|DYTURBOINPUTFILE|$in_file|g       " tmp
    sed -i "s|SETNPROCESSORS|$nprocessors|g     " tmp
    sed -i "s|SETPROGARGUMETS|$argumets|g     " tmp
    sed -i "s|SETLHAPDFLIB|`lhapdf-config --prefix`/lib|g           " tmp
    sed -i "s|SETLHAPDFDATA|`lhapdf-config --prefix`/share/LHAPDF|g " tmp
    [[ $target =~ lxbatch ]] && sed -i "s|^#BSUB -R.*$||g      "  tmp
    [[ $target =~ lxbatch ]] && sed -i "s|^#BSUB -app.*$|#BSUB -M 2000000|g " tmp
    [[ $target =~ lxbatch ]] && sed -i "s|/jobdir/|/pool/|g      "  tmp
    mv tmp $sh_file
    chmod +x $sh_file
}

submit_job(){
    if [[ $target == localrun ]]
    then
        # pickup random number
        export LSB_JOBINDEX=133
        $DRYRUN $sh_file
    else
        if [ -a $result_dir/$job_name.root ]
        then
            echo "Output already exists, skipping submission."
        else
            $DRYRUN bsub < $sh_file
        fi
    fi
}


DRYRUN=echo 
ask_submit(){
    read -p "Do you want to submit jobs ? " -n 1 -r
    echo    # (optional) move to a new line
    if [[ $REPLY =~ ^[Yy]$ ]]
    then
        DRYRUN=
    fi
}


## run MAIN
main $*


