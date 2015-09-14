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
dyturbo_project=`pwd -P`
batch_script_dir=$dyturbo_project/scripts/batch_scripts
in_files_dir=$dyturbo_project/scripts/infiles
result_dir=$dyturbo_project/results
# templates
batch_template=$dyturbo_project/scripts/run_DYTURBO_TMPL.sh
dyturbo_in_tmpl=$dyturbo_project/scripts/DYTURBO_TMPL.in


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
    nprocessors=$(($cubacores+1))
    #
    cp $batch_template tmp
    sed -i "s|JOBNAME|$job_name|g               " tmp
    sed -i "s|OUTDIR|$result_dir|g              " tmp
    sed -i "s|DYTURBOROOTDIR|$dyturbo_project|g " tmp
    sed -i "s|DYTURBOINPUTFILE|$in_file|g       " tmp
    sed -i "s|SETNPROCESSORS|$nprocessors|g       " tmp
    mv tmp $sh_file
    chmod +x $sh_file
}

prepare_in(){
    mkdir -p $in_files_dir
    in_file=$in_files_dir/$job_name.in
    # terms
    doRES="false"
    doCTM="false"
    doREA="false"
    doVIR="false"
    doLOR="false"
    if [[ $terms =~ ^ALL$ ]]
    then
        doRES="true"
        doCTM="true"
        doREA="true"
        doVIR="true"
        doLOR="true"
    fi;
    if [[ $terms =~ RES  ]]; then doRES="true"; fi;
    if [[ $terms =~ CT   ]]; then doCTM="true"; fi;
    if [[ $terms =~ REAL ]]; then doREA="true"; fi;
    if [[ $terms =~ VIRT ]]; then doVIR="true"; fi;
    if [[ $terms =~ LO   ]]; then doLOR="true"; fi;
    # special 3d
    resumdim=4
    if [[ $terms =~ RES3D  ]]; then resumdim=3; fi;
    # order check
    order=1
    if [[ $pdfset == CT10nlo   ]]; then order=1; fi;
    if [[ $pdfset == CT10nnlo  ]]; then order=2; fi;
    if [[ $pdfset == CTZPT2    ]]; then order=2; fi;
    if [[ $pdfset == ZPT-CT10  ]]; then order=2; fi;
    # process (default z0)
    lomass=66.
    himass=116.
    rmass=91.1876
    width=2.495
    nproc=3
    if [[ $process =~ ^w[pm]$ ]]
    then
        lomass=10.
        himass=1000.
        rmass=80.385
        width=2.091
    fi
    if [[ $process =~ ^wp$ ]]; then nproc=1; fi;
    if [[ $process =~ ^wm$ ]]; then nproc=2; fi;
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
    sed -i "s|LOQTBIN|$loqtbin|g      " tmp
    sed -i "s|HIQTBIN|$hiqtbin|g      " tmp
    sed -i "s|LOYBIN|$loybin|g        " tmp
    sed -i "s|HIYBIN|$hiybin|g        " tmp
    sed -i "s|SETSEED|$random_seed|g  " tmp
    sed -i "s|SETMLO|$lomass|g        " tmp
    sed -i "s|SETMHI|$himass|g        " tmp
    sed -i "s|SETRMASS|$rmass|g       " tmp
    sed -i "s|SETWIDTH|$width|g       " tmp
    sed -i "s|SETNPROC|$nproc|g       " tmp
    sed -i "s|SETPDFSET|$pdfsetname|g " tmp
    sed -i "s|SETMEMBER|$member|g     " tmp
    sed -i "s|SETORDER|$order|g       " tmp
    sed -i "s|SETGPAR|$gpar|g         " tmp
    sed -i "s|SETDORES|$doRES|g         " tmp
    sed -i "s|SETDOCTM|$doCTM|g         " tmp
    sed -i "s|SETDOREA|$doREA|g         " tmp
    sed -i "s|SETDOVIR|$doVIR|g         " tmp
    sed -i "s|SETDOLOR|$doLOR|g         " tmp
    sed -i "s|SETCUBACORES|$cubacores|g " tmp
    sed -i "s|SETRESUMDIM|$resumdim|g   " tmp
    mv tmp $in_file
}

submit_job(){
    if [[ `hostname` =~ cuth-dell  ]]
    then
        bash -x $sh_file
    else
        bsub < $sh_file
    fi
}

jobsDone(){
    bsub < scripts/send_mail.sh
}

make_range(){
    lobin=$1
    hibin=$2
    Nbins=$3
    make_range=""
    for i in `seq 0 $Nbins`
    do
        step=` echo "scale=1; ($hibin - $lobin) / $Nbins" | bc ` #$((   ))
        val=` echo "scale=1; $lobin + $i * $step" | bc -l | sed -r "s/^\./0./g" ` # bc <<< "$lobin + $i * $step"` #$((   )) $(( $lobin + $i * $step ))
        make_range="$make_range $val"
    done
    echo -n $make_range
}

submit_Z_dyturbo(){
    DRYRUN=echo 
    read -p "Do you want to submit jobs ? " -n 1 -r
    echo    # (optional) move to a new line
    if [[ $REPLY =~ ^[Yy]$ ]]
    then
        DRYRUN=
    fi
    # PROC SWITCH
    process=z0
    process=wp
    # TERM SWITCH
    terms=RESCT
    terms=VIRT
    terms=REAL
    # PDF SWITCH
    pdfset=CT10nnlo
    pdfset=CTZPT2
    pdfset=ZPT-CT10
    # defaults
    cubacores=8
    random_seed=100101
    loqtbin=unset
    hiqtbin=unset
    collider=lhc7
    seedlist=100101
    for process in z0 # z0 wp
    do
        for terms in RES CT REAL VIRT # REAL # RES CT VIRT  # RES3D CT REAL VIRT
        do
            # Z0
            #qtbinlist="0 2 4 6 8 10 12 14 16 18 22 26 30 34 38 42 46 50 54 60 70 80 100 150 200 300 800"
            #ybinlist="0 1 2 2.4"
            #variationList=" `seq 0 50` g_05 g_15 as_0117 as_0119"
            # WPM
            #if [[ $process =~ ^w[pm]$ ]]
            #then
            qtbinlist=`make_range 0.0 10.0 20.0` #`make_range 0.0 100.0 200.0`
            jqtbinmax=5
            ybinlist=`make_range 0.0 2.0 10.0` #`make_range 0.0 5.0 25.0`
            jybinmax=5
            variationList="0"
            #fi
            if [[ $terms =~ RES|CT|REAL|VIRT ]]
            then
                qtbinlist=`make_range 0.0 100.0 50.0` #`make_range 0.0 100.0 50.0`
                ybinlist="0 1 2 3 4 5" #"0 1 2 3 4 5"
                seedlist=`seq 100101 100101`
                cubacores=8
                jqtbinmax=1
                jybinmax=1
            fi
            # special run full space
            qtbinlist="0 100"
            ybinlist="0 5" #"0 1 2 3 4 5"
            seedlist=`seq 100101 100101`
            cubacores=8
            jqtbinmax=1
            jybinmax=1
            for random_seed in $seedlist
            do
                loqtbin=unset
                jqtbin=1
                allqtbins=""
                for iqtbin in $qtbinlist
                do
                    if [[ $loqtbin == unset ]]; then
                        loqtbin=$iqtbin;
                        continue;
                    fi;
                    # merge pt bins
                    allqtbins="$allqtbins $iqtbin"
                    if [[ $jqtbin != "$jqtbinmax" ]]; then
                        jqtbin=$(( $jqtbin + 1 ))
                        continue;
                    fi;
                    hiqtbin=$allqtbins
                    jqtbin=1
                    allqtbins=""
                    # end merge pt bins
                    #
                    loybin=unset
                    jybin=1
                    allybins=""
                    for iybin in $ybinlist
                    do
                        if [[ $loybin == unset ]]; then
                            loybin=$iybin;
                            continue;
                        fi;
                        # merge y bins
                        allybins="$allybins $iybin"
                        if [[ $jybin != "$jybinmax" ]]; then
                            jybin=$(( $jybin + 1 ))
                            continue;
                        fi;
                        hiybin=$allybins
                        jybin=1
                        allybins=""
                        # end merge y bins
                        #
                        qtregion=`echo qt${loqtbin}${iqtbin}y${loybin}${iybin}t${terms} | sed "s/\.//g;s/ //g"`
                        for variation in $variationList
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
            done
        done
    done
    #$DRYRUN jobsDone
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
        rm -f results/*.{root,log,out,err}
    fi
    echo "Done"
}

# MAIN
clear_files
clear_results
submit_Z_dyturbo
