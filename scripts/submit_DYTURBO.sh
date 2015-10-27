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
program=unset
fiducial=unset
detfiducial=0
makelepcuts=false
gpar=1

prepare_script(){
    mkdir -p $batch_script_dir
    mkdir -p $result_dir
    sample=${pdfset}_${variation}
    job_name=${program}_${process}_${collider}_${sample}_${qtregion}_${random_seed}
    sh_file=$batch_script_dir/$job_name.sh
    echo $job_name
    # prepare in
    prepare_in
    # prepare script
    nprocessors=$(($cubacores+1))
    #
    cp $batch_template tmp
    sed -i "s|JOBNAME|$job_name|g               " tmp
    sed -i "s|SEEDLIST|$seedlist|g              " tmp
    sed -i "s|OUTDIR|$result_dir|g              " tmp
    sed -i "s|DYTURBOROOTDIR|$dyturbo_project|g " tmp
    sed -i "s|DYTURBOINPUTFILE|$in_file|g       " tmp
    sed -i "s|SETNPROCESSORS|$nprocessors|g     " tmp
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
        termstring="tota"
    fi;
    if [[ $terms =~ RES  ]]; then doRES="true"; termstring="none"; fi;
    if [[ $terms =~ CT   ]]; then doCTM="true"; termstring="none"; fi;
    if [[ $terms =~ REAL ]]; then doREA="true"; termstring="real"; fi;
    if [[ $terms =~ VIRT ]]; then doVIR="true"; termstring="virt"; fi;
    if [[ $terms =~ LO   ]]; then doLOR="true"; termstring="lord"; fi;
    # special 3d
    resumdim=4
    ctdim=8
    if [[ $terms =~ RES3D  ]]; then resumdim=3; fi;
    if [[ $terms =~ CT3D   ]]; then    ctdim=3; fi;
    if [[ $terms =~ CT2D   ]]; then    ctdim=2; fi;
    # order check
    order=1
    if [[ $pdfset == CT10nlo    ]]; then order=1; fi;
    if [[ $pdfset == CT10nnlo   ]]; then order=2; fi;
    if [[ $pdfset == CTZPT2     ]]; then order=2; fi;
    if [[ $pdfset == ZPT-CT10   ]]; then order=2; fi;
    if [[ $pdfset == WZZPT-CT10 ]]; then order=2; fi;
    if [[ $pdfset == CT14       ]]; then order=2; fi;
    # process (default z0)
    lomass=66.
    himass=116.
    rmass=91.1876
    width=2.495
    nproc=3 
    nprocmcfm=41
    #nprocmcfm=41
    #if [[ $order == 2 ]]; then nprocmcfm=44; fi;
    if [[ $process =~ ^w[pm]$ ]]
    then
        lomass=10
        himass=1000
        rmass=80.385
        width=2.091
    fi
    if [[ $process =~ ^wp$ ]];
    then
        nproc=1;
        nprocmcfm=11;
        #if [[ $order == 1 ]]; then nprocmcfm=11; fi;
    fi;
    if [[ $process =~ ^wm$ ]];
    then
        nproc=2;
        nprocmcfm=16;
        #if [[ $order == 1 ]]; then nprocmcfm=16; fi;
    fi;
    # fiducial
    detfiducial=0
    if [[ $fiducial == D0    ]]; then detfiducial=1; fi;
    if [[ $fiducial == CDF   ]]; then detfiducial=2; fi;
    if [[ $fiducial == ATLAS ]]; then detfiducial=3; fi;
    if [[ $fiducial == CMS   ]]; then detfiducial=4; fi;
    # collider: default lhc7
    ih1=1
    ih2=1
    sroot=7e3
    if [[ $collider == tev1 ]]
    then
        ih2=-1
        sroot=1.8e3
        if [[ $process == z0     ]]; then lomass=30; himass=150; fi;
        if [[ $process =~ ^w[pm] ]]; then lomass=40; himass=120; fi;
    fi;
    if [[ $collider == tev2 ]]; then ih2=-1; sroot=1.96e3; fi;
    if [[ $collider == lhc8 ]]; then sroot=8e3; fi;
    # variations
    #gpar=1e0
    member=0
    re='^[0-9]+$'
    pdfsetname=$pdfset
    if [[ $variation == g_05  ]]; then gpar=0.5e0;                        fi;
    if [[ $variation == g_15  ]]; then gpar=1.5e0;                        fi;
    if [[ $variation == as_*  ]]; then pdfsetname=${pdfset}_${variation}; fi;
    if [[ $variation =~ $re   ]]; then member=$variation;                 fi;
    # set correct input template
    in_tmpl=$dyturbo_in_tmpl
    if [[ $job_name =~ ^dyturbo_ ]]; then in_tmpl=$dyturbo_project/scripts/DYTURBO_TMPL.in; fi;
    if [[ $job_name =~ ^dyres_   ]]; then in_tmpl=$dyturbo_project/scripts/DYRES_TMPL.in;   fi;
    if [[ $job_name =~ ^mcfm_    ]]; then in_tmpl=$dyturbo_project/scripts/MCFM_TMPL.in;    fi;
    cp $in_tmpl tmp
    sed -i "s|LOQTBIN|$loqtbin|g           ;
            s|HIQTBIN|$hiqtbin|g           ;
            s|LOYBIN|$loybin|g             ;
            s|HIYBIN|$hiybin|g             ;
            s|SETSEED|$random_seed|g       ;
            s|SETMLO|$lomass|g             ;
            s|SETMHI|$himass|g             ;
            s|SETRMASS|$rmass|g            ;
            s|SETWIDTH|$width|g            ;
            s|SETNPROC|$nproc|g            ;
            s|SETPDFSET|$pdfsetname|g      ;
            s|SETMEMBER|$member|g          ;
            s|SETORDER|$order|g            ;
            s|SETGPAR|$gpar|g              ;
            s|SETDORES|$doRES|g            ;
            s|SETDOCTM|$doCTM|g            ;
            s|SETDOREA|$doREA|g            ;
            s|SETDOVIR|$doVIR|g            ;
            s|SETDOLOR|$doLOR|g            ;
            s|SETCUBACORES|$cubacores|g    ;
            s|SETRESUMDIM|$resumdim|g      ;
            s|SETCTDIM|$ctdim|g            ;
            s|SETMCFMPROC|$nprocmcfm|g     ;
            s|SETTERMSTRING|$termstring|g  ;
            s|SETLEPCUTS|$makelepcuts|g    ;
            s|SETDETFIDUCIAL|$detfiducial|g;
            s|SETIH1|$ih1|g                ;
            s|SETIH2|$ih2|g                ;
            s|SETSROOT|$sroot|g            ; " tmp
    mv tmp $in_file
}

submit_job(){
    if [[ `hostname` =~ cuth-dell  ]]
    then
        $DRYRUN bash -x $sh_file
    else
        if [ -a results/$job_name.root ]
        then
            echo "Output already exists, skipping submission."
        else
            $DRYRUN bsub < $sh_file
        fi
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
    program=dyturbo
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
                            submit_job
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

submit_Wwidth(){
    # ask about running/submitting
    DRYRUN=echo 
    read -p "Do you want to submit jobs ? " -n 1 -r
    echo    # (optional) move to a new line
    if [[ $REPLY =~ ^[Yy]$ ]]
    then
        DRYRUN=
    fi
    # general setup
    program=dyturbo
    loqtbin=0
    hiqtbin=1000
    loybin=-5
    hiybin=5
    collider=lhc7
    cubacores=8
    variation=0
    #
    seedlist=10100
    random_seed=seed
    batch_template=$dyturbo_project/scripts/run_DYTURBO_Array_TMPL.sh
    # submit
    for process in z0 # wp wm
    do
        for order in 2 # 1 2
        do
            pdfset=CT10nlo
            termlist="RES CT LO"
            if [[ $order == 2 ]];
            then
                pdfset=CT10nnlo
                termlist="RES CT REAL VIRT"
            fi
            #fiducialsList="FULL D0 CDF ATLAS CMS"
            fiducialsList="FULL D0"
            for fiducial in $fiducialsList
            do
                #colliderlist="tev1 tev2 lhc7 lhc8"
                colliderlist="tev1"
                makelepcuts=false
                if [[ $fiducial == D0    ]]; then colliderlist=tev1; makelepcuts=true; fi;
                if [[ $fiducial == CDF   ]]; then colliderlist=tev2; makelepcuts=true; fi;
                if [[ $fiducial == ATLAS ]]; then colliderlist=lhc7; makelepcuts=true; fi;
                if [[ $fiducial == CMS   ]]; then colliderlist=lhc8; makelepcuts=true; fi;
                for collider in $colliderlist
                do
                    for terms in $termlist
                    do
                        qtregion=`echo f${fiducial}qt${loqtbin}${hiqtbin}y${loybin}${hiybin}t${terms} | sed "s/\.//g;s/ //g"`
                        prepare_script
                        submit_job
                    done
                done
            done
        done
    done
}

submit_allProg(){
    # ask about running/submitting
    DRYRUN=echo 
    read -p "Do you want to submit jobs ? " -n 1 -r
    echo    # (optional) move to a new line
    if [[ $REPLY =~ ^[Yy]$ ]]
    then
        DRYRUN=
    fi
    # full phase space
    loqtbin=0
    hiqtbin=100
    loybin=0
    hiybin=5
    collider=lhc7
    random_seed=100101
    startSeed=100201
    cubacores=8
    variation=0
    gpar=.83175
    # testing the array submission
    batch_template=$dyturbo_project/scripts/run_DYTURBO_Array_TMPL.sh
    #for program in dyres
    for program in dyturbo # mcfm # dyturbo dyres mcfm
    do
        for process in wp # wp wm z0
        do
            makelepcuts=false
            #if [[ $process =~ z0 ]]; then makelepcuts=true; fi
            for order in 2 # 3
            do
                # set PDF ?
                pdfset=CT10nlo
                if [[ $order == 2 ]]; then pdfset=ZPT-CT10; fi;
                if [[ $order == 3 ]]; then pdfset=WZZPT-CT10; order=2; fi;
                # set terms
                termlist="ALL"
                if [[ $program =~ ^dyturbo ]] 
                then
                    termlist="RES CT LO"
                    #if [[ $order == 2 ]]; then termlist="RES CT REAL VIRT"; fi;
                    if [[ $order == 2 ]]; then termlist="REAL"; fi;
                fi
                if [[ $program =~ ^dyres ]] 
                then
                    cubacores=0
                    termlist="ALL"
                    if [[ $order == 2 ]]; then termlist="ALL REAL VIRT"; fi;
                fi
                if [[ $program =~ ^mcfm ]] 
                then
                    cubacores=0
                    termlist="LO"
                    if [[ $order == 2 ]]; then termlist="REAL VIRT"; fi;
                fi
                #
                for terms in $termlist
                do
                    #NSeeds=20
                    #if [[ terms == RES ]]; then NSeeds=100; fi;
                    #if [[ terms == REAL ]]; then NSeeds=500; fi;
                    NSeeds=50
                    NSeeds=100
                    endSeed=$(( $startSeed + $NSeeds - 1 ))
                    seedlist="1000-1999"
                    #seedlist="5000"
                    random_seed=seed
                    #for random_seed in `seq $startSeed $endSeed`
                    #do
                        #if [[ $terms != REAL ]]; then continue; fi;
                        qtregion=`echo qt${loqtbin}${hiqtbin}y${loybin}${hiybin}t${terms} | sed "s/\.//g;s/ //g"`
                        prepare_script
                        submit_job
                    #done
                done
            done
        done
    done
}


clear_files(){
    echo "Clearing setup files"
    if [[ $(bjobs 2> /dev/null) ]] 
    then
        echo There are jobs running, skip clearing
    else
        rm -f scripts/batch_scripts/*.sh
        rm -f scripts/infiles/*.in
    fi
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
#submit_Z_dyturbo
submit_allProg
#submit_Wwidth
