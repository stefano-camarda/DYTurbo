#!/bin/bash

##
## @file submit_jobs_wmass.sh
##
## Description of the script file
##
## @author Jakub Cuth <Jakub.Cuth@cern.ch>
## @date 2016-06-16

#set -x

help(){
    #echo " Help is on TODO list :), But you can read parse_input function to check possibilities"

    echo "
USAGE: ./scripts/submit_jobs_wmass.sh TARGET OPTIONS --seeds N [--griduser USER --gridvoms VOMS]

TARGETS : (Need to specify one of those)
    --local         Prepare for run locally.
    --grid          Prepare grid submission directory (compilation on local)
    --mogon         Prepare standard scripts for mogon (Mainz cluster).
    --lxbatch       Prepare standard scripts for lxbatch (CERN LXPLUS batch system).
    --help          Print this help and die.

OPTIONS :
               [default]    {available}            Description

    --infile   [wmass.in]   {input file path}      Common input file with all settings.
    --proc     [z0]         {z0,wp,wm}             Set the process and mass integration (can be comma separated list)
    --pdfset   [CT10nnlo]   {LHAPDFname}           Set LHAPDF set name (can be comma separated list)
    --pdfvar   [0]          {int or all or array}  Set member number or run all
    --order    [1]          {1,2}                  1: NLL+NLO 2: NNLL+NNLO
    --term     [LO]         {RES,CT}               Monte-Carlo integration with order as above
                            {RES3D,CT3D}           Cubature 3D integration with order set above
                            {RES2D,CT2D}           Cubature 2D integration with order set above
                            {REAL,VIRT,LO}         Real, virt and V+J LO  with MC
                            {VV,FIXCT,FIXCT2D}     Fixed terms
    --seeds                 {int or range or list} MANDATORY: Set range (for batch) or Njobs (grid)
    --griduser              {GRID username}        MANDATORY IF GRID
    --gridvoms              {voms settings}        If you want to run with group privileges.

    "


}

main(){
    parse_input $*
    clear_files
    #
    submit_jobs_wmass
    #
    if [[ $target =~ lxplus|mogon|localrun  ]]
    then
        ask_submit
        $DRYRUN submit_jobs_wmass
    elif [[ $target =~ grid  ]]
    then
        echo
        echo " Your GRID username: $cernuser"
        echo " Your GRID group: $cerngroup"
        echo " Your GRID voms: $gridofficial"
        ask_submit
        # check panda voms proxy
        if ! voms-proxy-info || ! which prun
        then
            echo "Please setup panda and voms : "
            echo "setupATLAS "
            echo "lsetup rucio panda "
            echo "voms-proxy-init -voms atlas -valid 96:00"
            exit 7
        fi
        if [[ $DRYRUN =~ echo ]] 
        then
            echo "To submit to grid just type:"
        fi
        for gridpdf in $pdflist
        do
            $DRYRUN cd $dyturbo_project/GRID_$gridpdf
            $DRYRUN ./subm.sh
        done
    else
        echo "UKNOWN TARGET: $target"
    fi
}

parse_input(){
    target=unset
    proclist=z0
    #
    order=1
    #termlist="LO VV FIXCT"
    termlist=LO
    #order=2
    #termlist="REAL VIRT FIXCT VV"
    #
    pdflist="CT10nnlo CT10nnlo68clProfiled"
    #pdfvarlist=all
    #pdfvarlist="0 1 2 3"
    pdfvarlist=0
    #
    infile=input/wmass.in
    #
    #seedlist=100,201 # run 100 and 201
    #seedlist=100-201 # run form 100 to 201
    #seedlist=800    # run from 1 to 800
    seedlist=unset
    #
    voms=unset
    cernuser=unset

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
            --mbins)
                mbins=$2
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
            --griduser)
                cernuser=$2
                shift
                ;;
            --voms)
                voms=$2
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
    # check infile
    if ! ls $infile > /dev/null
    then
        echo " Input file not found pleae set correct with '--infile name'  "
        exit 5
    fi
    # check seed
    if [[ $seedlist == unset ]]
    then
        echo " Seed argument is mandatory please set '--seeds [integer/range/list]'"
        exit 5
    fi
    if [[ $target =~ grid ]]
    then
        [[ $cernuser == unset ]] \
            && echo " GRID usernaname is mandatory please set '--griduser NAME'" \
            && echo " You can also specify your voms by '--voms VOMS'" \
            &&   exit 6
    else # not grid
        [[ $seedlist =~ -|, ]]  || seedlist=1-$seedlist
    fi
    # check order term
    [[ $order == 1 ]] && [[ $termlist =~ REAL|VIRT ]] && echo "WRONG ORDER $order TO TERM $termlist" && return 3
    [[ $order == 2 ]] && [[ $termlist =~ LO        ]] && echo "WRONG ORDER $order TO TERM $termlist" && return 3 
    # loops, lopps and loops
    for pdfset in $pdflist 
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
                    [[ $target =~ grid ]]  && submit_job2grid
                done # variation
            done #term
        done #proc
        [[ $target =~ grid ]] && finalize_grid_submission
    done #pdfsets
}

#directories 
dyturbo_project=`pwd -P`
batch_script_dir=$dyturbo_project/scripts/batch_scripts
queue=atlasshort
batch_template=$dyturbo_project/scripts/run_on_batch_tmpl.sh
collider=lhc7

# this goes away
in_files_dir=unset
cubacores=unset
job_name=unset
in_file=unset
sh_file=unset
process=unset
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
    mbins=32,40,200
    [[ $process =~ z0 ]] && mbins=10,66,116
    # arguments
    arguments="input.in --proc $process --mbins $mbins --pdfset $pdfset --pdfvar $variation --order $order --term $terms"
    [[ $target =~ lxbatch|mogon|localrun ]] && arguments=$arguments" --seed $LSB_JOBINDEX "
    # make sure we make some noise on grid
    [[ $target =~ mogon ]] && arguments=$arguments" --verbose "
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
    sed -i "s|DYTURBOINPUTFILE|$dyturbo_project/$infile|g       " tmp
    sed -i "s|SETNPROCESSORS|$nprocessors|g     " tmp
    sed -i "s|SETPROGARGUMETS|$arguments|g     " tmp
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
        export LSB_JOBINDEX=1
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

clear_files(){
    echo -n "Clearing setup files .. "
    if [[ $(bjobs -o name 2> /dev/null | grep dyturbo) ]] 
    then
        echo -n "There are jobs running, skip clearing .. "
        #rm -f scripts/batch_scripts/*.sh
        #rm -f scripts/infiles/*.in
    else
        rm -f scripts/batch_scripts/*.sh
        #rm -f scripts/infiles/*.in
    fi
    rm -f scripts/infiles/*.tar
    rm -f scripts/infiles/*.tar.gz
    rm -f scripts/cmd_list
    echo "Done"
}

submit_job2grid(){
    echo "PRUN $job_name $seedlist $arguments" >> scripts/grid_submit.cmd
}

prepare_tarbal(){
    griddir=GRID_$pdfset
    gridv=v`date +%s`
    DYTURBOVERSION=`grep PACKAGE_VERSION config.h | cut -d\" -f2`
    # check you have cvmfs rootconfig
    rootprefix=`root-config --prefix`
    if [[ $rootprefix =~ /cvmfs/atlas.cern.ch/.*/[56].[0-9][0-9].[0-9][0-9]-(x86_64|i686)-slc[56]-gcc ]]
    then
        ROOTVERSION=`basename $rootprefix | cut -d- -f1`
        CMTVERSION=`basename $rootprefix | cut -d- -f2-`
    else
        echo "Current ROOT path is '$rootprefix'"
        echo "Your root setup is not comming from cvmfs, please do 'lsetup root'."
        exit 4
    fi
    # prepare executables or tarbal
    if [[ $target =~ compile ]]
    then
        echo "NOT TESTED $target" && exit 6
        echo "Making a tarbal.. please wait"
        if ! make dist > /dev/null
        then
            echo "Compilation problem.. Try to 'make install'. I am exiting."
            exit 3
        fi
    else
        echo "Making a executable... please wait"
        if ! make install > /dev/null
        then
            echo "Compilation problem.. Try to 'make install'. I am exiting."
            exit 3
        fi
        # test the lhapdf
        PATH=$PATH:lhapdf6/bin
        if ! lhapdf-config --prefix > /dev/null
        then
            echo "Cannot reach lhapdf-config. Please setup the PATH to LHAPDF"
            exit 3
        fi
    fi
    # official mode
    gridofficial=" " # space at the end
    cerngroup=user
    if [[ $voms != unset ]]
    then
        gridofficial=" --official --voms $voms"
        cerngroup=group
    fi
    # prepare PRUN command
    cat  scripts/run_prun.sh              >  scripts/grid_submit.cmd
    echo ""                               >> scripts/grid_submit.cmd
    echo "gridv=$gridv"                   >> scripts/grid_submit.cmd
    echo "CERNUSER=$cernuser"             >> scripts/grid_submit.cmd
    echo "CERNGROUP=$cerngroup"           >> scripts/grid_submit.cmd
    echo "OFFICIAL=$gridofficial"         >> scripts/grid_submit.cmd
    echo "DYTURBOVERSION=$DYTURBOVERSION" >> scripts/grid_submit.cmd
    echo "ROOTVERSION=$ROOTVERSION"       >> scripts/grid_submit.cmd
    echo "CMTVERSION=$CMTVERSION"         >> scripts/grid_submit.cmd
    echo "target=$target"                 >> scripts/grid_submit.cmd
    echo ""                               >> scripts/grid_submit.cmd
    echo "# Setup rucio and panda "       >> scripts/grid_submit.cmd
    echo "# lsetup rucio panda"           >> scripts/grid_submit.cmd
    echo "# voms-proxy-init -voms atlas -valid 96:00" >> scripts/grid_submit.cmd
}

finalize_grid_submission(){
    CP="rsync -a"
    #echo "compressing... $tarbalfile"
    #gzip $tarbalfile
    if [[ $target =~ compile ]]
    then
        # put source code
        $CP dyturbo-${DYTURBOVERSION}.tar.gz $griddir/
    else
        # prepare folders
        mkdir -p $griddir/bin
        mkdir -p $griddir/lib
        mkdir -p $griddir/LHAPDF
        # copy exec
        $CP bin/* $griddir/bin/
        $CP lib/*.so* $griddir/lib/
        $CP lib/*.la  $griddir/lib/
        $CP `lhapdf-config --libdir`/libLHAPDF.so $griddir/lib/
        # copy pdfset
        lhapdfdir=`lhapdf-config --datadir`
        $CP $lhapdfdir/$pdfset $griddir/LHAPDF
        $CP $lhapdfdir/lhapdf.conf $griddir/LHAPDF
        $CP $lhapdfdir/pdfsets.index $griddir/LHAPDF
    fi
    # submision scripts
    cat scripts/grid_submit.cmd > $griddir/subm.sh
    chmod +x $griddir/subm.sh
    # on-site scripts
    sed "s|SETTARGET|$target|g" scripts/compile_grid.sh > $griddir/compile_grid.sh
    sed "s|SETTARGET|$target|g" scripts/run_grid.sh     > $griddir/run_grid.sh
    # add inputfile
    $CP input/default.in $griddir/
    $CP $infile $griddir/input.in
    #
    #ls -hla --color=auto $griddir
    #echo
    #echo "Go to GRID folder 'cd GRID' edit subm.sh (change user name, role) and run it './subm.sh' "
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


