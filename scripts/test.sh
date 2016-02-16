


echo Starting test
#set -x

#make -C Cuba-4.2

#  ./configure --enable-debug &&
#  make install &&
#  tdir=run_000 &&
#  rm -rf $tdir &&
#  mkdir $tdir &&
#  cd $tdir &&
#  
#  /usr/bin/time -v ./../bin/dyturbo ../input/test.in &&
#  cd .. &&
#  
####

#./configure --enable-root --enable-Ofast &&
make install &&

#./bin/dyturbo input/test.in &&
#hadd -f merge.root results*.root &&
#echo -e " #include<iostream> \n void integr(){ double err; std::cout << qt_y_total->IntegralAndError(-1,-1,-1,-1,err) << std::endl; cout << err << endl; std::cout << h_qtVy->IntegralAndError(-1,-1,-1,-1,err) << std::endl; cout << err << endl; }" > scripts/integr.C &&
#root -l -q merge.root scripts/integr.C && #$SCRIPTSYS/misc/Browse.C &&

#export LSB_JOBINDEX=133 && rm -rf results/dyturbo_wm_lhc7_WZZPT-CT10_0_qt0100y05tRES_seed_$LSB_JOBINDEX.root &&

#./scripts/submit_DYTURBO.sh &&
#./scripts/batch_scripts/dyturbo_wm_lhc7_WZZPT-CT10_0_qt0100y05tRES_seed.sh &&
#tdir=run_000 &&
#rm -rf $tdir &&
#mkdir $tdir &&
#cd $tdir &&


#/usr/bin/time -v ./../bin/dyturbo ../input/test.in &&
#gdb --args ./../bin/dyturbo ../input/test.in &&

#hadd results_merge.root results*.root &&


#python scripts/pt_plot.py &&

/usr/bin/time -v ./bin/merger -X real.root results_Stefano/AiMoments/z-8tev-nnll-real/{1,2,3,4}*/results.root &&
/usr/bin/time -v ./bin/merger -X virt.root results_Stefano/AiMoments/z-8tev-nnll-virt/{1,2,3,4}*/results.root &&

python scripts/pt_plot.py &&
#set +x
echo Finished test || echo TEST FAILURE

