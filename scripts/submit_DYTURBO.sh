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
tarbalfile=unset
seedlist=unset

prepare_script(){
    mkdir -p $batch_script_dir
    mkdir -p $result_dir
    sample=${pdfset}_${variation}
    job_name=${program}_${process}_${collider}_${sample}_${qtregion}_${random_seed}
    sh_file=$batch_script_dir/$job_name.sh
    echo $job_name $seedlist
    # prepare in
    prepare_in
    # prepare script
    nprocessors=$(($cubacores+1))
    walltime=5:00
    queue=atlasshort
    if [[ $program =~ dyres ]]
    then
        walltime=20:00
        queue=atlaslong
    fi
    #
    cp $batch_template tmp
    sed -i "s|JOBNAME|$job_name|g               " tmp
    sed -i "s|SEEDLIST|$seedlist|g              " tmp
    sed -i "s|OUTDIR|$result_dir|g              " tmp
    sed -i "s|SETQUEUE|$queue|g              " tmp
    sed -i "s|SETWALLTIME|$walltime|g              " tmp
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
    # order check -- trust user
    # order=1
    # if [[ $pdfset == CT10nlo    ]]; then order=1; fi;
    # if [[ $pdfset == CT10nnlo         ]]; then order=2; fi;
    # if [[ $pdfset == CTZPT2           ]]; then order=2; fi;
    # if [[ $pdfset == ZPT-CT10         ]]; then order=2; fi;
    # if [[ $pdfset == WZZPT-CT10       ]]; then order=2; fi;
    # if [[ $pdfset == CT14             ]]; then order=2; fi;
    # if [[ $pdfset == MMHT2014nnlo68cl ]]; then order=2; fi;
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
    if [[ $fiducial == CMS7   ]]; then detfiducial=4; fi;
    if [[ $fiducial == CMS8   ]]; then detfiducial=5; fi;
    # collider: default lhc7
    ih1=1
    ih2=1
    sroot=7e3
    # TEV1
    if [[ $collider == tev1 ]]
    then
        ih2=-1
        sroot=1.8e3
        if [[ $process == z0     ]]; then lomass=75; himass=105; fi;
        if [[ $process =~ ^w[pm] ]]; then lomass=40; himass=120; fi;
    fi;
    # TEV2
    if [[ $collider == tev2 ]];
    then
        ih2=-1;
        sroot=1.96e3;
        if [[ $process == z0     ]]; then lomass=30; himass=500; fi;
    fi
    # TEV2
    if [[ $collider == lhc8 ]]
    then
        sroot=8e3
        if [[ $process == z0     ]] && [[ $fiducial =~ CMS ]]; then lomass=60; himass=120; fi;
    fi
    # variations
    #gpar=1e0
    member=0
    re='^[0-9]+$'
    pdfsetname=$pdfset
    if [[ $variation == g_05  ]]; then gpar=0.5e0;                        fi;
    if [[ $variation == g_15  ]]; then gpar=1.5e0;                        fi;
    if [[ $variation == as_*  ]]; then pdfsetname=${pdfset}_${variation}; fi;
    if [[ $variation =~ $re   ]]; then member=$variation;                 fi;
    # gpar
    if [[ $pdfset == WZZPT-CT10 ]] 
    then
        gpar=0.9097
        if [[ $variation == 53    ]]; then gpar=0.97330;                      fi;
        if [[ $variation == 54    ]]; then gpar=0.84610;                      fi;
    fi;
    if [[ $pdfset == ZPT-CT10 ]] 
    then
        gpar=0.83175
        if [[ $variation == 53    ]]; then gpar=0.88990;                      fi;
        if [[ $variation == 54    ]]; then gpar=0.77360;                      fi;
    fi;
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

prepare_tarbal(){
    tarbalfile=$in_files_dir/${program}_${pdfset}_${gridv}.tar
    exclude="-X scripts/excl" #'--exclude="*.o" --exclude="*.lo" --exclude="*.Po" --exclude="*.Plo" --exclude="*.deps*" --exclude="*.a" --exclude="*.pdf"'
    # add scripts and default input
    tar cf $tarbalfile --transform='s|.*/||g' scripts/run_grid.sh scripts/compile_grid.sh input/default.in
    # add libraries
    tar rf $tarbalfile lib/lib* bin/dyturbo
    #  # add autotools config
    #  tar rf $tarbalfile configure.ac install-cuba dyturbo-config.in Makefile.am input/
    #  # add dyturbo source code
    #  tar rf $tarbalfile $exclude src/ dyres/ mcfm/ dynnlo/ dyres/ cernlib/ Cuba-4.2/
    # add wanted PDFset
    tar rhf $tarbalfile -C lhapdf6/share/LHAPDF/ $pdfset/ lhapdf.conf pdfsets.index
}

finalize_grid_submission(){
    echo "compressing... $tarbalfile"
    gzip $tarbalfile
    rm -rf GRID/*
    rsync -avP $tarbalfile.gz GRID/
    cat scripts/grid_submit.cmd | column -t > GRID/subm.sh
    chmod +x GRID/subm.sh
    ls -l GRID
}


add_to_tarbal(){
    # hack the interations
    sed -i "s|^cubaverbosity   *=.*$|cubaverbosity    = 2  |g" $in_files_dir/$job_name.in
    sed -i "s|^vegasncallsRES  *=.*$|vegasncallsRES   = 5e5|g" $in_files_dir/$job_name.in
    sed -i "s|^vegasncallsCT   *=.*$|vegasncallsCT    = 5e7|g" $in_files_dir/$job_name.in
    sed -i "s|^vegasncallsLO   *=.*$|vegasncallsLO    = 1e8|g" $in_files_dir/$job_name.in
    sed -i "s|^vegasncallsREAL *=.*$|vegasncallsREAL  = 2e8|g" $in_files_dir/$job_name.in
    sed -i "s|^vegasncallsVIRT *=.*$|vegasncallsVIRT  = 1e8|g" $in_files_dir/$job_name.in
    # add input file
    tar rf $tarbalfile scripts/infiles/$job_name.in
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

submit_job2grid(){
    $DRYRUN 
    intargz=`basename $tarbalfile`.gz
    echo prun --exec \". run_grid.sh ${job_name} %RNDM:1 \" \
    --outDS user.\$GRIDUSER.${job_name}.out \
    --outputs=HIST:results_merge.root \
    --noCompile \
    --nJobs $seedlist \
    --inTarBall=$intargz \
    >> scripts/grid_submit.cmd
    #--site ANALY_CERN_SHORT \
    #--excludeFile="out_*"
    #--nJobs 1 \
    #--long \
    #--bexec \"./compile_grid.sh\" \
    #--rootVer=6.02/12 --cmtConfig=x86_64-slc6-gcc48-opt \
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
    for process in z0 wp wm
    do
        for order in 2 # 1 2
        do
            termlist="RES CT LO"
            if [[ $order == 2 ]];
            then
                termlist="RES CT REAL VIRT"
            fi
            #fiducialsList="FULL D0 CDF ATLAS CMS"
            fiducialsList="D0 CDF ATLAS CMS7 CMS8"
            for fiducial in $fiducialsList
            do
                #colliderlist="tev1 tev2 lhc7 lhc8"
                colliderlist="tev1"
                makelepcuts=false
                if [[ $fiducial == D0    ]]; then colliderlist=tev1; makelepcuts=false; pdfset=CTEQ4M;           fi;
                if [[ $fiducial == CDF   ]]; then colliderlist=tev2; makelepcuts=false; pdfset=CTEQ5L;           fi;
                if [[ $fiducial == ATLAS ]]; then colliderlist=lhc7; makelepcuts=false; pdfset=CTEQ66;           fi;
                if [[ $fiducial == CMS7  ]]; then colliderlist=lhc7; makelepcuts=false; pdfset=CTEQ6L; fi;
                if [[ $fiducial == CMS8  ]]; then colliderlist=lhc8; makelepcuts=false; pdfset=MSTW2008nnlo68cl; fi;
                pdfset=MMHT2014nnlo68cl
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
    variation=0
    gpar=.83175
    # testing the array submission
    batch_template=$dyturbo_project/scripts/run_DYTURBO_Array_TMPL.sh
    for program in dyres #  dyturbo dyres mcfm
    do
        for process in wm # wp wm z0
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
                    cubacores=8
                    termlist="RES CT LO"
                    #if [[ $order == 2 ]]; then termlist="RES CT REAL VIRT"; fi;
                    if [[ $order == 2 ]]; then termlist="REAL"; fi;
                    seedlist="1000-1999"
                fi
                if [[ $program =~ ^dyres ]] 
                then
                    cubacores=0
                    termlist="ALL"
                    if [[ $order == 2 ]]; then termlist="ALL REAL VIRT"; fi;
                    seedlist="5000"
                fi
                if [[ $program =~ ^mcfm ]] 
                then
                    cubacores=0
                    termlist="LO"
                    if [[ $order == 2 ]]; then termlist="REAL VIRT"; fi;
                fi
                for terms in $termlist
                do
                    random_seed=seed
                    qtregion=`echo qt${loqtbin}${hiqtbin}y${loybin}${hiybin}t${terms} | sed "s/\.//g;s/ //g"`
                    prepare_script
                    submit_job
                done
            done
        done
    done
}

submit_grid(){
    # ask about running/submitting == turned off due to separete file for prun
    #  DRYRUN=echo 
    #  read -p "Do you want to submit jobs ? " -n 1 -r
    #  echo    # (optional) move to a new line
    #  if [[ $REPLY =~ ^[Yy]$ ]]
    #  then
    #      DRYRUN=
    #  fi
    # lsetup rucio panda
    # voms-proxy-init atlas
    # full phase space
    loqtbin=0
    hiqtbin=100
    loybin=0
    hiybin=5
    collider=lhc7
    random_seed=100101
    startSeed=100201
    cubacores=0
    variation=0
    gpar=.83175
    program=dyturbo
    pdfset=WZZPT-CT10; order=2;
    #pdfset=ZPT-CT10; order=2;
    #
    gridv=v`date +%s`
    prepare_tarbal
    echo -e "#!/bin/bash\nGRIDUSER=jcuth\n" > scripts/grid_submit.cmd
    #
    for process in z0 # wp wm z0
    do
        makelepcuts=false
        #if [[ $process =~ z0 ]]; then makelepcuts=true; fi
        for variation in 0 # `seq 0 54`
        do
            #terms
            termlist="RES CT LO"
            if [[ $order == 2 ]]; then termlist="RES CT REAL VIRT"; fi;
            #termlist="RES"
            for terms in $termlist
            do
                seedlist=100
                if [[ $terms == REAL ]]; then seedlist=2000; fi;
                random_seed=seed
                seedlist=1
                # prepare config
                qtregion=`echo ${gridv}qt${loqtbin}${hiqtbin}y${loybin}${hiybin}t${terms} | sed "s/\.//g;s/ //g"`
                prepare_script
                add_to_tarbal
                submit_job2grid
            done
        done
    done
    finalize_grid_submission
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
    rm -f scripts/infiles/*.tar
    rm -f scripts/infiles/*.tar.gz
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
#submit_allProg
#submit_Wwidth
submit_grid
