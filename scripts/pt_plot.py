#!/usr/bin/python2.7
# -*- coding: utf-8 -*- 

## Documentation for file
#
# More details. 
#
# @file pt_plot.py
# @author cuto <Jakub.Cuth@cern.ch>
# @date 2015-08-07

# lsetup root "sft releases/pyanalysis/1.5_python2.7-0dd7c"

import sys,os,re
from copy import deepcopy
# for batch add -b before first ROOT calling 
sys.argv.append('-b') # run in batch

#from PlotTools import *
from ROOT  import TGraphErrors, TFile, TGraphAsymmErrors
import ROOT as RT
#import asciitable
from array import array

from PlotTools import *
pl=PlotTools()

#import array

def plot_pt(fname):

    TABLE = asciitable.read(fname, delimiter=",")
    gr = TGraphErrors()

    print TABLE
    i=0
    for d in TABLE :
        xval=float(d[1]+d[0])/2
        xerr=float(d[1]-d[0])/2
        tot_val=d[-4]
        tot_err=d[-3]
        print xval, xerr, yval, yerr
        gr.SetPoint(i,xval,yval)
        gr.SetPointError(i,xerr,yerr)
        i+=1
        pass
    #pl = PlotTools()
    pl.NewCanvas("dyturbo_qt_CT10nnlo");
    pl.SetFrameStyle1D([gr],minX=0,maxX=100)
    #gPad.SetLogy()
    gr.Draw("PLSAME");
    pl.Save()
    pl = PlotTools()
    pl.NewCanvas("dyturbo_qt_CT10nnlo40");
    pl.SetFrameStyle1D([gr],maxX=40)
    #gPad.SetLogy()
    gr.Draw("PLSAME");
    pl.Save()
    pass

def print_results():
    for variation in [str(x) for x in range(0,51)] + ["g_05", "g_15", "as_0117", "as_0119"]:
        fname = "results_merge/dyturbo_z0_lhc7_CT10nnlo_{0}_qtMerge_100101.root".format(variation)
        tf = TFile.Open(fname,"read")
        for term in ["resum" , "ct"]:
            gr = tf.Get("qt_"+term)
            print
            print " {0:^13} : {1:^29}".format("bin",term+" "+variation)
            print "-"*(1+ 5 +3+ 5 +3+ 12 +5+ 12)
            lobin = -1
            hibin = -1
            #for ibin in [0, 2, 4, 6, 8, 10, 12, 14, 16, 18, 22, 26, 30, 34, 38, 42, 46, 50, 54, 60, 70, 80, 100, 150, 200, 300, 800]:
            #if lobin < 0 :
            #lobin=ibin
            #continue
            for i in range(gr.GetN()):

                lobin = gr.GetX()[i] - gr.GetEX()[i]
                hibin = gr.GetX()[i] + gr.GetEX()[i]
                val = gr.GetY ()[i]
                err = gr.GetEY()[i]
                print " {0:5} - {1:5} : {2:12e} +/- {3:12e}".format(lobin,hibin,val,err)
                #lobin=ibin
                pass
            print
            pass
        tf.Close()
    pass



def print_table() :
    #pl = PlotTools()
    fname = "results_merge/merge_RESUM_CT.root"
    hname="qt_y_{}_{}"
    htmp = pl.GetHistSetName("dummy",fname,hname.format("resum","0"))
    xbinlist = range(1,htmp.GetNbinsX() +1)
    ybinlist = range(1,htmp.GetNbinsY() +1)
    #header
    line = "* lobin hibin"
    for variation in [str(x) for x in range(0,51)] + ["g_05", "g_15", "as_0117", "as_0119"]:
        line+=" "
        line+=variation
        pass
    print line
    #table
    for xbin in xbinlist :
        line = "";
        binlo = htmp.GetXaxis().GetBinLowEdge(xbin)
        binhi = htmp.GetXaxis().GetBinUpEdge(xbin)
        line += str(binlo) + " " + str(binhi)
        for variation in [str(x) for x in range(0,51)] + ["g_05", "g_15", "as_0117", "as_0119"]:
            val=0
            err2=0
            for term in ["resum" , "ct"] :
                h=pl.GetHist(fname,hname.format(term,variation))
                for ybin in ybinlist :
                    ibin = h.GetBin(xbin,ybin)
                    val+=      h.GetBinContent(ibin)
                    err2+= pow(h.GetBinError  (ibin),2)
                pass
            line += " " + str(val)
            #line += " " + str(sqrt(err2))
            pass
        print line
        pass
    pass

def draw_projections(hlist,proj):
    resum = 0
    ct    = 0
    real  = 0
    virt  = 0
    varname="qt"
    MAXX=3
    #
    proc = hlist[0].GetName().split("_")[0]
    #
    if "X" in proj :
        varname="qt"
        MAXX=30
        resum = hlist[0].ProjectionX("{}_{}_{}".format (proc, "resum", varname ))
        ct    = hlist[1].ProjectionX("{}_{}_{}".format (proc, "ct",    varname ))
        real  = hlist[2].ProjectionX("{}_{}_{}".format (proc, "real",  varname ))
        virt  = hlist[3].ProjectionX("{}_{}_{}".format (proc, "virt",  varname ))
    elif "Y" in proj :
        varname="y"
        MAXX=4
        resum = hlist[0].ProjectionY("{}_{}_{}".format (proc, "resum", varname ))
        ct    = hlist[1].ProjectionY("{}_{}_{}".format (proc, "ct",    varname ))
        real  = hlist[2].ProjectionY("{}_{}_{}".format (proc, "real",  varname ))
        virt  = hlist[3].ProjectionY("{}_{}_{}".format (proc, "virt",  varname ))
    else :
        resum = hlist[0]
        ct    = hlist[1]
        real  = hlist[2]
        virt  = hlist[3]
    # sum them
    total = resum.Clone("{}_{}_{}".format (proc, "tot",  varname ))
    total.Add(ct   )
    total.Add(real )
    total.Add(virt )
    finite = ct.Clone("{}_{}_{}".format (proc, "fin",  varname ))
    finite.Add(real)
    finite.Add(virt)
    # set colors
    #f = TFile(proc+"qt.root","RECREATE")
    hists = [ total, resum , finite, ct    , real  , virt  ]
    N = len(hists)
    for i,h in enumerate(hists):
        h.SetTitle(h.GetName())
        h.SetLineColor( pl.AutoCompareColor(i,N) )
        h.SetMarkerColor( pl.AutoCompareColor(i,N) )
        h.SetMarkerStyle(20+i)
        h.SetMarkerSize(1)
        h.SetDrawOption("E0")
        h.GetYAxis().SetTitle("d#sigma / d"+varname+"[pb]")
        my_integral = 0
        my_interr2 = 0
        print h.GetName()
        for ibin,binval in enumerate(h):
            errval = h.GetBinError(ibin)
            perc = 0. if binval==0 else abs(errval/binval*100)
            print ibin,binval,errval, "{:.2}%".format(perc)
            if ibin == 5 : break
            my_integral+=binval
            my_interr2+=pow(errval,2)
            pass
        errval = pow(my_interr2,.5)
        perc = 0. if my_integral==0. else abs(errval/my_integral*100)
        print " integral {} +- {} ({}%)".format(my_integral,errval,perc)
        print
        #h.Write()
        pass
    # compare all
    pl.CompareHistsInList("{}_{}_{}".format (proc, varname, "allterms"  ) , hists     , compareType="ratio", doStyle=False, maxX=MAXX)
    pl.CompareHistsInList("{}_{}_{}".format (proc, varname, "resct"     ) , hists[1:4], compareType="ratio", doStyle=False, maxX=MAXX)
    pl.CompareHistsInList("{}_{}_{}".format (proc, varname, "ctrealvirt") , hists[2:] , compareType="ratio", doStyle=False, maxX=MAXX)
    pl.CompareHistsInList("{}_{}_{}".format (proc, varname, "totres"    ) , hists[0:2], compareType="ratio", doStyle=False, maxX=MAXX)
    pl.CompareHistsInList("{}_{}_{}".format (proc, varname, "total"     ) , hists[0:1], compareType="ratio", doStyle=False, maxX=MAXX)
    #pl.CompareHistsInList(proc+"_qt_realorig"  , [realorig], doStyle=False, maxX=MAXX)
    pass

def w_pt(proc="wp"):
    #proc="wp"
    #proc="z0"
    #pl = PlotTools()
    fname = "results_merge/dyturbo_{}_lhc7_CTZPT2_0_qtyMerget{}_100101.root"
    #fname = "results/dyturbo_wp_lhc7_CTZPT2_0_qt0020y05t{}_100101.root"
    ftmpres = fname.format(proc,"RESCT")
    ftmpct  = fname.format(proc,"RESCT")
    ftmpreal= fname.format(proc,"REAL") # "dyturbo_wp_lhc7_CTZPT2_0_qt0020y05tREAL_100101.root"
    ftmpvirt= fname.format(proc,"VIRT") # "dyturbo_wp_lhc7_CTZPT2_0_qt0020y05tVIRT_100101.root"
    # get histograms
    resum = pl.GetHistSetName( "res2D"      ,   ftmpres , "qt_y_resum" ).ProjectionX( proc+"_resum_pt", 10, 15)
    ct    = pl.GetHistSetName( "ct2D"       ,   ftmpct  , "qt_y_ct"    ).ProjectionX( proc+"_ct_pt",    10, 15)
    real  = pl.GetHistSetName( proc+"_real_pt" ,ftmpreal,  "h_qt"       )
    virt  = pl.GetHistSetName( proc+"_virt_pt" ,ftmpvirt,  "h_qt"       )
    # fix them :D
    #realorig = real.Clone("wp_real_pt_NOCORR")
    #realorig.SetTitle("real before correction")
    # real from 4 - 6 GeV
    # bin_start = 10
    # bin_end = 11
    # val_start = real.GetBinContent(bin_start-1)
    # val_end = real.GetBinContent(bin_end+1)
    # slope = val_end-val_start
    # binlist = range (bin_start,bin_end+1)
    # N = len(binlist)
    # for i,ibin in enumerate(binlist):
    #     newval = val_start + slope * i / N
    #     real.SetBinContent(ibin,newval)
    #     pass
    # fix the normalization
    real.Scale(1)
    virt.Scale(1)
    hists= [ resum, ct, real, virt]
    draw_projections(hists,"none")
    pl.MakePreviewFromList(0,proc+"_qt")
    pass

def getFromFiles(proc, term):
    MIN_QT=0
    MIN_Y =0
    MAX_QT=10
    MAX_Y =2
    hout=0
    term_def = {
            #"resum" : "RES3D",
            "resum" : "RES",
            "ct"    : "CT",
            "real"  : "REAL",
            "virt"  : "VIRT",
            }
    hname    ="h_qtVy"
    hnameNorm="qt_y_"+term
    name="{}_{}_{}".format(proc,"qty",term)
    dir="results/"
    fterm =  term_def[term]
    matchstring=r'.*dyturbo_'+proc+'_.*t'+fterm+'_.*\.root'
    matchstring=r'.*dyturbo_'+proc+'_.*qt0100y05t'+fterm+'_.*\.root'
    files =  [dir+f for f in os.listdir(dir) if re.match(matchstring, f) ]
    files.sort()
    print files
    for i,fname in enumerate(files) :
        si = "_" +str(i)
        # get all files histograms
        h2      = pl.GetHistSetName(name+si,fname,hname)
        h2_norm = pl.GetHistSetName(name+"_norm"+si,fname,hnameNorm)
        if "blablabla" in term :
            h2 = h2_norm
        else :
            #normalise and check
            Ih2,err = getIntegralError(h2)
            print_res(h2.GetName(),Ih2,err)
            Inorm,err = getIntegralError(h2_norm)
            print_res(h2_norm.GetName(),Inorm,err)
            h2.Scale(Inorm/Ih2)
            Ih2,err = getIntegralError(h2)
            print_res(h2.GetName(),Ih2,err)
        if hout==0 :
            hout = h2.Clone(name)
        else :
            hout.Add(h2)
        pass
    return hout

def w_pt_y():
    proc="wp"
    #proc="z0"
    # 
    h_resum = getFromFiles(proc,"resum")
    h_ct    = getFromFiles(proc,"ct")
    h_real  = getFromFiles(proc,"real")
    h_virt  = getFromFiles(proc,"virt")
    #
    # h_fin = h_resum.Clone("fin")
    # h_fin .Add(h_ct   )
    # h_fin .Add(h_real )
    # h_fin .Add(h_virt )
    # h_tot = h_ct   .Clone("tot")
    # h_tot .Add(h_real )
    # h_tot .Add(h_virt )
    #
    hists2D = [h_resum, h_ct, h_real, h_virt]
    #
    draw_projections(hists2D,"X")
    draw_projections(hists2D,"Y")
    pl.MakePreviewFromList(0,proc+"_qty")
    pass


def getDataHistZ(hMC):
    DATA= [
            [  0       , 2     ,    0.02822    ,       0.2715   ,        0.3694  ,        0.3639   ,       0.5853 ],
            [  2       , 4     ,     0.0584    ,       0.1667   ,        0.3191  ,        0.3459   ,       0.4993 ],
            [  4       , 6     ,    0.05805    ,       0.1698   ,         0.228  ,        0.3614   ,       0.4598 ],
            [  6       , 8     ,    0.04917    ,       0.1847   ,        0.2151  ,        0.3602   ,       0.4584 ],
            [  8       , 10    ,    0.04076    ,        0.205   ,        0.2364  ,        0.3433   ,       0.4645 ],
            [  10      , 12    ,     0.0338    ,       0.2277   ,        0.2562  ,        0.3426   ,       0.4846 ],
            [  12      , 14    ,    0.02815    ,       0.2451   ,        0.2603  ,        0.3442   ,       0.4963 ],
            [  14      , 16    ,    0.02375    ,       0.2691   ,        0.2597  ,        0.3431   ,       0.5075 ],
            [  16      , 18    ,    0.02012    ,       0.2973   ,        0.2747  ,        0.3407   ,        0.529 ],
            [  18      , 22    ,    0.01595    ,        0.254   ,        0.2537  ,        0.3429   ,       0.4964 ],
            [  22      , 26    ,      0.012    ,       0.3024   ,        0.2815  ,        0.3557   ,       0.5452 ],
            [  26      , 30    ,   0.009166    ,       0.3406   ,        0.3084  ,        0.3566   ,       0.5816 ],
            [  30      , 34    ,   0.007242    ,       0.3902   ,        0.3275  ,        0.3525   ,       0.6195 ],
            [  34      , 38    ,   0.005802    ,       0.4375   ,        0.3548  ,        0.3473   ,       0.6618 ],
            [  38      , 42    ,   0.004641    ,       0.4889   ,        0.3896  ,        0.3499   ,       0.7163 ],
            [  42      , 46    ,   0.003777    ,         0.53   ,        0.4272  ,        0.3481   ,       0.7646 ],
            [  46      , 50    ,   0.003172    ,       0.5688   ,         0.433  ,        0.3673   ,       0.8037 ],
            [  50      , 54    ,   0.002593    ,       0.6426   ,        0.4601  ,        0.3673   ,       0.8715 ],
            [  54      , 60    ,   0.002104    ,        0.613   ,        0.4283  ,        0.3749   ,       0.8365 ],
            [  60      , 70    ,   0.001492    ,       0.5526   ,        0.4396  ,        0.3785   ,       0.8012 ],
            [  70      , 80    ,   0.0009851   ,       0.6894   ,        0.4946  ,        0.4282   ,       0.9505 ],
            [  80      , 100   ,   0.0005525   ,       0.6184   ,        0.4858  ,        0.4447   ,       0.9034 ],
            [  100     , 150   ,   0.0001918   ,       0.6264   ,        0.5307  ,        0.6531   ,        1.049 ],
            [  150     , 200   ,   4.891e-05   ,        1.258   ,        0.7213  ,         0.628   ,         1.58 ],
            [  200     , 300   ,   1.081e-05   ,        1.882   ,         1.402  ,         1.335   ,          2.7 ],
            [  300     , 800   ,   3.985e-07   ,        4.195   ,         2.045  ,         1.324   ,        4.851 ],
            ]
    binsx = [ x[0] for x in DATA ]
    binsx.append(DATA[-1][1])
    ar_bins = array('d',binsx)
    #for i,a in enumerate(binsx):
        #ar_bins[i]=binsx[i]
    N = len(binsx)-1
    hData = TH1D ("Zdata", "Z DATA;q_{T}[GeV]", N, ar_bins );
    for i,x in enumerate(DATA):
        hData.SetBinContent(i+1, x[2]);
        hData.SetBinError  (i+1, x[2]*x[6]/100);
    # rebin same MC
    hMC = hMC.Rebin(N,"Z MC",ar_bins)
    # divide by bin width
    for ibin in range(1,hMC.GetNbinsX()+1) :
        val = hMC.GetBinContent(ibin)
        width = hMC.GetBinWidth(ibin)
        hMC.SetBinContent(ibin,val/width)
        pass
    hData.Scale(1./hData.Integral())
    hMC.Scale(1./hMC.Integral())
    return [hData, hMC]










# add by hand
# dyturbo_z0_lhc7_CT10nnlo_9_qt1214y12_100101
# |    12 -    14 |     1 -     2 |   -2585.71 +/-0.000531072 (   1360.85s) | -2585.71 +/-0.000531072 (   1360.85s) |
# dyturbo_z0_lhc7_CT10nnlo_12_qt1416y22.4_100101
# |    14 -    16 |     2 -   2.4 |   -187.808 +/-1.76105e-06 (   1353.08s) | -187.808 +/-1.76105e-06 (   1353.08s) |

def getHistRename(name,fname,hname):
    ff = TFile(fname,"READ")
    o = ff.Get(hname)
    h = o.Clone(name)
    h.SetTitle(name)
    h.SetDirectory(0)
    return h

def fixJob(hist, binX, binY, val, err):
    bin = hist.FindBin(binX,binY)
    print "BEFORE ", hist.GetBinContent (bin), hist.GetBinError   (bin)
    hist.SetBinContent (bin, val)
    hist.SetBinError   (bin, err)
    print "AFTER ", hist.GetBinContent (bin), hist.GetBinError   (bin)
    return hist

def merge_all_hist():
    hlist = list()
    for merge in ["RESUM", "CT"] :
        for variation in [str(x) for x in range(0,51)] + ["g_05", "g_15", "as_0117", "as_0119"]:
            fname = "results_merge/dyturbo_z0_lhc7_CT10nnlo_{0}_qtyMerge{1}_100101.root".format(variation,merge)
            hname = "qt_y_{0}".format(merge.lower())
            name = hname+"_"+variation
            h =  getHistRename(name, fname, hname)
            if "CT" in merge :
                if "9" == variation :
                    h = fixJob(h,13,1.5,-2585.71,0.000531072)
                if "12" == variation :
                    h = fixJob(h,15,2.2,-187.808 , 1.76105e-06)
            hlist.append ( h )
            pass
        pass
    # write outfile
    ff = TFile.Open("results_merge/merge_RESUM_CT.root","RECREATE")
    for h in hlist :
        h.Write()
        pass
    ff.Write()
    ff.Close()
    pass


def print_res(name, var, err):
    errpc = 0 if err==0 else abs(err/var)*100
    print " {:10} : {:>+15.4f} +- {:>15.4f} ({:.2}%) ".format( name, float(var), float(err), float(errpc) )
    pass

def getIntegralError(h):
    # integral
    val=0
    err = Double(0.)
    doFull = False
    if h.GetDimension()==2 :
        xmin= -1 if doFull else 1
        xmax= -1 if doFull else h.GetNbinsX()+1
        ymin= -1 if doFull else 1
        ymax= -1 if doFull else h.GetNbinsY()+1
        val = h.IntegralAndError( xmin, xmax, ymin, ymax, err)
    else :
        xmin= -1 if doFull else 1
        xmax= -1 if doFull else h.GetNbinsX()+1
        val = h.IntegralAndError( xmin, xmax, err)
        pass
    return val,float(err)

def check_cancelation():
    pl = PlotTools()
    terms=[
       [ "RES"   , "resum" ] ,
       [ "CT"    , "ct"    ] ,
       [ "REAL"  , "real"  ] ,
       [ "VIRT"  , "virt"  ] ,
    ]
    processes= [ "z0" , "wp" ]
    filetmp="results/dyturbo_{}_lhc7_CTZPT2_0_qt0020y05t{}_100101.root"
    #histtmp="qt_y_virt"
    hists=list()
    for proc in processes :
        print proc
        fin = 0
        err2fin = 0
        tot = 0
        err2tot= 0
        htot=0
        hfin=0
        for term in terms :
            histname="h_qt"
            histname="h_qtVy"
            #histname="h_y"
            histres="qt_y_"+term[1]
            fname=filetmp.format(proc,term[0])
            if "resum" in term[1] :
                fname="results_merge/dyturbo_{}_lhc7_CTZPT2_0_qtyMerget{}_100101.root".format(proc,term[0])
                histname=histres
                pass
            h =  pl.GetHistSetName(proc+"_"+term[1], fname ,histname)
            hres =  pl.GetHistSetName(proc+"norm_"+term[1], fname ,histres)
            #print time.ctime(os.path.getmtime(fname))
            #
            val=0
            err=0
            #normalize properly
            val,err = getIntegralError(h)
            norm,normerr = getIntegralError(hres)
            h.Scale(norm/val)
            val,err = getIntegralError(h)
            err2 = float(err*err)
            #
            print_res(term[1], val, err)
            #
            hists.append(h)
            if htot != 0:
                htot .Add(h)
                if(term[1]!="resum") : hfin .Add(h)
            if htot == 0:
                htot = h.Clone(proc+"_"+"tot")
                hfin = h.Clone(proc+"_"+"fin")
            if term[1]!="resum" :
                fin+=val
                err2fin += err2
            tot += val
            err2tot+= err2
            pass
        hists.append(fin)
        hists.append(tot)
        #
        print_res("fin" , fin, sqrt(err2fin))
        print_res("tot" , tot, sqrt(err2tot))
        #
        w_pt(proc,hists)
    pass

def addIntegralError(v1,e1,v2,e2):
    return v1+v2, sqrt(e1*e1+e2*e2)

def scale_err(hist,scale):
    if scale!=scale : return hist
    for ibin,binval in enumerate(hist):
        olderr=hist.GetBinError(ibin)
        hist.SetBinError(ibin, olderr*scale)
    return hist

def root_file_integral():
    #filetmp="results/dyturbo_z0_lhc7_CT10nlo_0_qt0100y05t{}_100101.root"
    #filetmp="results/dyturbo_z0_lhc7_CT10nnlo_0_qt0100y05t{}_100101.root"
    filetmp="results/mcfm_z0_lhc7_CT10nlo_0_qt0100y05t{}_100101.root"
    #filetmp="results/mcfm_z0_lhc7_CT10nnlo_0_qt0100y05t{}_100101.root"
    #filetmp="run_dir/results4.root"
    #filetmp="results/dyturbo_z0_lhc7_CT10nnlo_0_qt010y01t{}_100101.root"
    #
    programs  = ["dyturbo"] #, "mcfm"]
    orders    = ["ZPT-CT10"] #, "CT10nlo", "CT10nnlo"]
    #
    #filetmp1="results/{}_{}_lhc7_{}_0_qt0100y05t{}_100101.root"
    #filetmp2="results_merge/dipole_plotting/{}_{}_lhc7_{}_0_qt0100y05t{}_100101.root"
    #filetmp="results_merge/{}_{}_lhc7_{}_0_qt0100y05t{}_merge.root"
    #filetmp="results_merge/{}_{}_lhc7_{}_0_qt0100y05t{}_outliers.root"
    #
    #processes = ["z0"] "wp","wm","z0"] , "z0"]
    #filetmp="results_merge/z0_fiducial_profilingClosure/{}_{}_lhc7_{}_0_qt0100y05t{}_outliers.root"
    #
    programs  = ["dyturbo"] #, "mcfm"]
    processes = [ "wp" ] #,"wm" ]
    filetmp="results_merge/wpm_predictions_151012/{}_{}_lhc7_{}_0_qt0100y05t{}_outliers.root"
    filetmpREAL="results_merge/wp_real_1000seeds_151026/{}_{}_lhc7_{}_0_qt0100y05t{}_outliers.root" #results_merge/wpm_300seed_151019/{}_{}_lhc7_{}_0_qt0100y05t{}_seed_outliers.root" results_merge/wpm_300seeds_151019/dyturbo_wp_lhc7_ZPT-CT10_0_qt0100y05tREAL_seed_outliers.root
    #            
    #
    hist="h_qtVy"
    hist_integr="qt_y_{}"
    TERMS=[ 
            ( "TOT" , "tot"  ),
            ( "FIN" , "fin"  ),
            ( "RES"  , "resum" ),
            #( "RES3D", "resum" ),
            ( "CT"   , "ct"    ),
            #( "CT3D" , "ct"    ),
            ( "LO"   , "lo"    ),
            ( "REAL" , "real"  ),
            #( "REALOUT" , "real"  ),
            ( "VIRT" , "virt"  ),
            #( "ALL"  , "total" ),
            ] 
    #
    #programs  = ["dyres"] #, "mcfm"]
    #processes = [ "wp" ] #,"wm" ]
    #filetmp="results_merge/dyres_closure/{}_{}_lhc7_{}_0_qt0100y05t{}_seed_outliers.root"
    #TERMS= [
            #( "ALL", "total" ),
            #]
    #
    for proc in processes :
        for prog in programs :
            for pdf in orders:
                print " {} {} {}".format(prog,proc,pdf)
                tot=0
                toterr=0
                htot=0
                fin=0
                finerr=0
                hfin=0
                hlist=list()
                titltmp="{}_{}_{}_{}_{}"
                for term,term_lowcas in TERMS :
                    hist_res_term = term_lowcas
                    filename=filetmp
                    Nfiles = 1
                    if "mcfm" in prog : hist_res_term = "total"
                    if "FIN"  in term : hist_res_term = "total"
                    if "TOT"  in term : hist_res_term = "total"
                    if "REAL"  in term : filename=filetmpREAL
                    try :
                        h = pl.GetHist(filename.format(prog,proc,pdf,term),hist)
                    except ValueError:
                        continue
                    except :
                        raise
                    titl=titltmp.format(prog,proc,pdf,term,"qty")
                    h.SetName(titl)
                    h.SetTitle(proc+" "+term)
                    h.Scale(1./Nfiles)
                    h.SetContour(200)
                    # clone grandtotal hist
                    if htot==0 :
                        titl=titltmp.format(prog,proc,pdf,"total","qty")
                        htot = h.Clone(titl)
                        htot.SetTitle(titl)
                        htot.Reset()
                    if hfin==0 :
                        titl=titltmp.format(prog,proc,pdf,"finite","qty")
                        hfin = h.Clone(titl)
                        hfin.SetTitle(titl)
                        hfin.Reset()
                    # zero the total (mcfm summing only real+virt)
                    if "mcfm" in prog and "REAL" in term :
                        tot =0
                        toterr=0
                        htot.Reset()
                        fin =0
                        finerr=0
                        hfin.Reset()
                        pass
                    integ , inerr = getIntegralError(h)
                    print_res(term,integ,inerr)
                    # from full phase space
                    hres = pl.GetHist(filetmp.format(prog,proc,pdf,term),hist_integr.format(hist_res_term))
                    integr , ineer = getIntegralError(hres)
                    print_res(term_lowcas,integr,ineer)
                    if not "3D" in term:
                        hlist.append(h)
                        if not "TOT" in term and not "FIN" in term: 
                            tot,toterr = addIntegralError(tot,toterr, integr,ineer)
                            htot.Add(h)
                            if not "RES" in term:
                                fin,finerr = addIntegralError(fin,finerr, integr,ineer)
                                hfin.Add(h)
                                pass
                            pass
                        pass
                    pass
                # after all terms
                integ , inerr = getIntegralError(hfin)
                print_res("hFIN",integ,inerr)
                print_res("hfin",fin,finerr)
                #
                integ , inerr = getIntegralError(htot)
                print_res("hTOT",integ,inerr)
                print_res("htot",tot,toterr)
                #
                hlist.insert(0,hfin)
                hlist.insert(0,htot)
                # 0   1   2       3       4   5  6    7
                # tot fin haddtot haddfin res ct real virt
                #hlistProject = [ hlist[4], hlist[2], ]
                hlistProject = list()
                hlistProject.append(pl.GetHistSetTitNam("real 1000j",filetmpREAL.format("dyturbo","wp","ZPT-CT10","REAL"), hist))
                hlistProject.append(pl.GetHistSetTitNam("real 300j", "results_merge/wpm_real_300seeds_151019/dyturbo_wp_lhc7_ZPT-CT10_0_qt0100y05tREAL_seed_outliers.root", hist))
                hlistProject.append(pl.GetHistSetTitNam("real 100j", "results_merge/wpm_real_100seeds_151016/dyturbo_wp_lhc7_ZPT-CT10_0_qt0100y05tREAL_outliers.root", hist))
                hlistProject .append(  hlist[6] ); hlistProject[-1].SetTitle("real 50j")
                #
                #
                hlistProject = [ hlist[4], hlist[2], ]
                hlistProject.append(pl.GetHistSetTitNam("dyres", "results_merge/dyres_closure/dyres_wp_lhc7_ZPT-CT10_0_qt0100y05tALL_seed_outliers.root", hist))
                #
                hlqt = [ x.ProjectionX(x.GetName()+"_qt") for x in hlistProject ]
                hly  = [ x.ProjectionY(x.GetName()+"_y" ) for x in hlistProject ]
                # set axis title
                hlistProject[0].GetXaxis().SetTitle("q_{T}[GeV]")
                hlistProject[0].GetYaxis().SetTitle("y")
                hlqt[0].GetXaxis().SetTitle("q_{T}[GeV]")
                hly[0] .GetXaxis().SetTitle("y")
                # set 
                pl.CompareHistsInList( titltmp.format(prog,proc,pdf,"all","qt" ) , hlqt , maxX=30, compareType="ratio0" )
                pl.CompareHistsInList( titltmp.format(prog,proc,pdf,"all","y"  ) , hly  , maxX=30, compareType="ratio0" )
                #gStyle.SetPalette(RT.kDarkBodyRadiator)
                pl.CompareHistsInList( titltmp.format(prog,proc,pdf,"tot","qty" ) , hlistProject[0:1])
                # data MC comparison
                if "z0" in proc :
                    pl.CompareHistsInList( titltmp.format(prog,proc,pdf,"datamc","qt" ) , getDataHistZ(hlqt[1]), maxX=100, compareType="ratio0")
                continue
                Pallets= [
                    ["kDeepSea",                  RT.kDeepSea                  ],
                    ["kGreyScale",                RT.kGreyScale                ],
                    ["kDarkBodyRadiator",         RT.kDarkBodyRadiator         ],
                    ["kBlueYellow",               RT.kBlueYellow               ],
                    ["kRainBow",                  RT.kRainBow                  ],
                    ["kInvertedDarkBodyRadiator", RT.kInvertedDarkBodyRadiator ],
                    ["kBird",                     RT.kBird                     ],
                    ["kCubehelix",                RT.kCubehelix                ],
                    ["kGreenRedViolet",           RT.kGreenRedViolet           ],
                    ["kBlueRedYellow",            RT.kBlueRedYellow            ],
                    ["kOcean",                    RT.kOcean                    ],
                    ["kColorPrintableOnGrey",     RT.kColorPrintableOnGrey     ],
                    ["kAlpine",                   RT.kAlpine                   ],
                    ["kAquamarine",               RT.kAquamarine               ],
                    ["kArmy",                     RT.kArmy                     ],
                    ["kAtlantic",                 RT.kAtlantic                 ],
                    ["kAurora",                   RT.kAurora                   ],
                    ["kAvocado",                  RT.kAvocado                  ],
                    ["kBeach",                    RT.kBeach                    ],
                    ["kBlackBody",                RT.kBlackBody                ],
                    ["kBlueGreenYellow",          RT.kBlueGreenYellow          ],
                    ["kBrownCyan",                RT.kBrownCyan                ],
                    ["kCMYK",                     RT.kCMYK                     ],
                    ["kCandy",                    RT.kCandy                    ],
                    ["kCherry",                   RT.kCherry                   ],
                    ["kCoffee",                   RT.kCoffee                   ],
                    ["kDarkRainBow",              RT.kDarkRainBow              ],
                    ["kDarkTerrain",              RT.kDarkTerrain              ],
                    ["kFall",                     RT.kFall                     ],
                    ["kFruitPunch",               RT.kFruitPunch               ],
                    ["kFuchsia",                  RT.kFuchsia                  ],
                    ["kGreyYellow",               RT.kGreyYellow               ],
                    ["kGreenBrownTerrain",        RT.kGreenBrownTerrain        ],
                    ["kGreenPink",                RT.kGreenPink                ],
                    ["kIsland",                   RT.kIsland                   ],
                    ["kLake",                     RT.kLake                     ],
                    ["kLightTemperature",         RT.kLightTemperature         ],
                    ["kLightTerrain",             RT.kLightTerrain             ],
                    ["kMint",                     RT.kMint                     ],
                    ["kNeon",                     RT.kNeon                     ],
                    ["kPastel",                   RT.kPastel                   ],
                    ["kPearl",                    RT.kPearl                    ],
                    ["kPigeon",                   RT.kPigeon                   ],
                    ["kPlum",                     RT.kPlum                     ],
                    ["kRedBlue",                  RT.kRedBlue                  ],
                    ["kRose",                     RT.kRose                     ],
                    ["kRust",                     RT.kRust                     ],
                    ["kSandyTerrain",             RT.kSandyTerrain             ],
                    ["kSienna",                   RT.kSienna                   ],
                    ["kSolar",                    RT.kSolar                    ],
                    ["kSouthWest",                RT.kSouthWest                ],
                    ["kStarryNight",              RT.kStarryNight              ],
                    ["kSunset",                   RT.kSunset                   ],
                    ["kTemperatureMap",           RT.kTemperatureMap           ],
                    ["kThermometer",              RT.kThermometer              ],
                    ["kValentine",                RT.kValentine                ],
                    ["kVisibleSpectrum",          RT.kVisibleSpectrum          ],
                    ["kWaterMelon",               RT.kWaterMelon               ],
                    ["kCool",                     RT.kCool                     ],
                    ["kCopper",                   RT.kCopper                   ],
                    ["kGistEarth",                RT.kGistEarth                ],
                    #["kViridis",                  RT.kViridis                  ]
                ]
                for name, code in Pallets :
                    gStyle.SetPalette(code)
                    pl.CompareHistsInList( titltmp.format(prog,proc,pdf,"tot","qty"+name ) , hlistProject[0:1])
                pass
            pass
        pass
    pl.MakePreviewFromList(0,"all_compare")
    pass

def quick_calc():
    pb="pb"
    print "nnlo"
    print "z0 mcfm fb"
    print_res("real", -846121.25, 1568.0395)
    print_res("virt", 2129394.3 , 775.93249)
    mytot,myerr=addIntegralError(-846121.25, 1568.0395, 2129394.3 , 775.93249 )
    print_res("myto",  mytot, myerr)
    print_res("tota",  1283273.0886 ,  1749.5196 )
    print "z0 dyturbo"
    mytot,myerr= 0,0
    print_res( "ct"   , -288420 , 180.135 ); mytot,myerr = addIntegralError(mytot,myerr, -288420 , 180.135 )
    print_res( "real" , -412256 , 3149.67 ); mytot,myerr = addIntegralError(mytot,myerr, -412256 , 3149.67 )
    print_res( "virt" , 710744  , 230.544 ); mytot,myerr = addIntegralError(mytot,myerr, 710744  , 230.544 )
    print_res( "fin"  , mytot,myerr )
    print_res( "res"  , 450026  , 153.251 ); mytot,myerr = addIntegralError(mytot,myerr, 450026  , 153.251 )
    print_res( "tot"  , mytot,myerr )
    print "nlo"
    print "z0 mcfm fb"
    print_res("real", -846121.25, 1568.0395)
    print_res("virt", 2129394.3 , 775.93249)
    mytot,myerr=addIntegralError(-846121.25, 1568.0395, 2129394.3 , 775.93249 )
    print_res("myto",  mytot, myerr)
    print_res("tota",  1283273.0886 ,  1749.5196 )
    print "z0 dyturbo"
    mytot,myerr= 0,0
    print_res( "ct"   , -421637 ,   138.485 ); mytot,myerr = addIntegralError(mytot,myerr, -288420 , 180.135 )
    print_res( "lord" , 435281 ,    128.38 ); mytot,myerr = addIntegralError(mytot,myerr, -412256 , 3149.67 )
    print_res( "fin"  , mytot,myerr )
    print_res( "res"  ,  442687 ,   163.995); mytot,myerr = addIntegralError(mytot,myerr, 450026  , 153.251 )
    print_res( "tot"  , mytot,myerr )
    pass


def GetValError(fiducial, col, proc,full=False):
    #file_tmpl="results_merge/dyturbo_{}_{}_CT10nnlo_0_f{}qt01000y-55t{}_100101.root"
    #file_tmpl="results/dyturbo_{}_{}_CT10nnlo_0_f{}qt01000y-55t{}_100101.root"
    #file_tmpl="results_merge/wwidth_properD0massCut/dyturbo_{}_{}_CT10nnlo_0_f{}qt01000y-55t{}_100101.root"
    #file_tmpl_D0="results_merge/wwidth_D0_zmumu/dyturbo_{}_{}_CT10nnlo_0_f{}qt01000y-55t{}_seed_10100.root"
    file_tmpl="results_merge/wwidth_151029_MMHT/dyturbo_{}_{}_MMHT2014nnlo68cl_0_f{}qt01000y-55t{}_seed_10100.root"
    term="TOT"
    hist="h_qtVy"
    if full :
        hist="qt_y_total"
    totInt,totIne=0,0
    for term in [ "RES", "CT", "REAL", "VIRT" ]:
        ft = file_tmpl
        #if "z0" in proc and "tev1" in col and "D0" in fiducial : ft = file_tmpl_D0
        #if "z0" in proc and "tev1" in col : ft = file_tmpl_D0
        #print ft ,proc,col,fiducial,term
        h = pl.GetHist(ft.format(proc,col,fiducial,term),hist)
        integ , inerr = getIntegralError(h)
        print term, integ, inerr, inerr/integ
        totInt,totIne = addIntegralError(integ,inerr,totInt,totIne)
        pass
    print "TOT ", totInt, totIne, totIne/totInt
    print
    return totInt, totIne/totInt

def print_wwidth_res(int,err):
    print int,err,

def wwidth_table():
    experiments = [
            [ "D0"    , "tev1" ],
            [ "CDF"   , "tev2" ],
            [ "ATLAS" , "lhc7" ],
            [ "CMS7"  , "lhc7" ],
            [ "CMS8"  , "lhc8" ],
        ]
    processes = [ "wp", "wm", "z0"]
    for exp,col in experiments :
        #new header
        print " table for exp", exp, "@", col, " in order W+ W- Z"
        print " full  | fiducial | eff "
        for proc in processes :
            # full
            #full_v, full_e_r = GetValError("FULL",col,proc)
            full_v, full_e_r = GetValError(exp,col,proc,True)
            print_wwidth_res(full_v, full_e_r)
            # fiducial
            fid_v, fid_e_r = GetValError(exp,col,proc)
            print_wwidth_res(fid_v, fid_e_r)
            # efficiency
            eff_v = fid_v / full_v
            eff_e_r = sqrt ( full_e_r**2 + fid_e_r**2 )
            print_wwidth_res(eff_v, eff_e_r)
            # newline
            print
            pass
        pass
        # ratio
    pass


def find_fluctuations():
    outnametmp="{}_{}_{}"
    titletmp = "{} {} {}"
    #filetmp = "results_merge/wpm_predictions_151012/dyturbo_{}_lhc7_ZPT-CT10_0_qt0100y05t{}_outliers.root"
    #filetmpREAL = "results_merge/wp_real_1000seeds_151026/dyturbo_{}_lhc7_ZPT-CT10_0_qt0100y05t{}_outliers.root"
    filetmp = "results_merge/grid_151116/dyturbo_{}_lhc7_WZZPT-CT10_0_v1447428851qt0100y05t{}_outliers.root"
    proc = "wp"
    hname="h_qtVy"
    terms = [
            "CT",
            "RES",
            "REAL",
            "VIRT",
            #"REAL 1000"
            ]
    for proc in [ "wp", "wm", "z0" ]:
        for var in [ "qt", "y"] :
            fluct_list=list()
            for term in terms :
                title = titletmp.format(proc,term,var)
                filename = filetmp.format(proc,term)
                #if "REAL 1000" == term :
                    #filename = filetmpREAL.format(proc,"REAL")
                hist2d = pl.GetHistSetTitNam(titletmp.format(proc,term,"qty"),filename,hname)
                hist = 0
                if var == "qt" :
                    hist = hist2d.ProjectionX(title)
                elif var == "y" :
                    hist = hist2d.ProjectionY(title)
                else : 
                    raise ValueError(" unknown var ")
                hist_moving = pl.CreateMovingAverageHist(hist,1)
                hist_ratio = pl.CreateRatio0Hists(hist_moving,hist)
                fluct_list.append(hist_ratio)
                fluct_list[-1].SetTitle(title)
                pass
            pl.CompareHistsInList(outnametmp.format(proc,var,"MovAvg"),fluct_list, compareType=None)
            # zooming
            MINY=0.8
            MINX=0
            MAXX=3.5
            if "qt" == var :
                MINY=-300
                MINX=5
                MAXX=50
            pl.CompareHistsInList(outnametmp.format(proc,var,"MovAvgZoom"), fluct_list, minX=MINX, maxX=MAXX,minY=MINY,compareType=None)
            pass
    pl.MakePreviewFromList(0,"mov_avg")
    pass


def PlotBand(name,hlist):
    cent=hlist[0]
    pl.AutoSetStyle(hlist) #,"c")
    band = TGraphAsymmErrors()
    for i in range(1,cent.GetNbinsX()+1) :
        x  = cent.GetBinCenter (i)
        y  = cent.GetBinContent(i)
        ey1 = hlist[1].GetBinContent(i)-y
        ey2 = hlist[2].GetBinContent(i)-y
        ex = cent.GetBinWidth(i)/2.
        band.SetPoint(i-1,x,y)
        eyh=abs(ey1)
        eyl=abs(ey2)
        if ey1*ey2 < 0 and ey1 < 0 :
                eyl = abs(ey1)
                eyh = abs(ey2)
        band.SetPointError(i-1,ex,ex,eyl,eyh)
        pass
    band.Draw()
    band.SetName(name+"_bandG")
    band.SetTitle(name+"_bandG")
    band.SetDrawOption("3")
    band.SetFillColor(860-9)
    mainHists = [ band,cent]
    subHists = pl.CreateRatio0Hists(cent,hlist)
    pl.NewCanvas(name+"band")
    # var
    var = name.split("_")[2]
    maxx=40
    if var=="y" :
        maxx=3.8
    pl.DrawHistCompareSubPlot(mainHists,subHists, drawOpt="3", compareType="ratio0", maxX=maxx, cdiv=0.3)
    pl.Save()
    pass

def PlotUnc2D(name,hlist):
    #cent
    cent=hlis[0]
    hists = [cent]
    #pos=var-cent
    hists .append(hlist[1])
    hists[-1].Add(cent,-1)
    hists[-1].SetTitle("pos")
    #neg=var-cent
    hists .append(hlist[1])
    hists[-1].Add(cent,-1)
    hists[-1].SetTitle("neg")
    #pos+neg
    hists .append(hists[-2].Clone("pos+neg"))
    hists[-1].Add(hists[-2])
    hists[-1].SetTitle("posPneg")
    # plot all
    for h in hists:
        pl.NewCanvas(name+h.Title)
        pl.SetFrameStyle2D([h])
        pl.Save()
    pass

def CreateUncertPlot(name,variations,hname,fbase):
    centvar = ["0"]
    centvar += [ str (x) for x in variations]
    #hlist = list()
    #hlist.append(pl.GetHistSetTitNam(name,hname,fbase.format())
    hlist = [pl.GetHistSetTitNam(name,fbase.format(x),hname) for x in centvar]
    plotclass = hlist[0].ClassName()
    if plotclass == "TH1D" :
        PlotBand(name,hlist)
    elif plotclass == "TH2D" :
        PlotUnc2D(name,hlist)
    else :
        raise ValueError("not know plotclass: "+plotclass)
    pass


def uncert_as_g():
    processes = [
            #"wp",
            "wm",
            #"z0",
            ]
    plots = [
            "h_qt",
            "h_y",
            ]
    filetmp = "results_merge/grid_151116/dyturbo_{}_lhc7_WZZPT-CT10_{}_v1447428851qt0100y05t{}_outliers.root"
    for proc in processes:
        filebase= filetmp.format(proc,"{}","TOT")
        for plot in plots :
            CreateUncertPlot( "_".join([proc,plot,"alphaS" ]), [51,52] , plot, filebase)
            #CreateUncertPlot( "_".join(proc,plot,"g"      ), [53,54] , plot, filebase)
        pass
    pl.MakePreviewFromList(0,"unc")
    pass


class  makeInfo:
    def __init__(s,procD=0,plotD=0,termD=0,fileD=0):
        # NULL
        s.proc_name = ""
        s.proc_titl = ""
        #
        s.plot_name = ""
        s.plot_titl = ""
        s.plot_maxx = ""
        #
        s.term_name  = ""
        s.term_title = ""
        s.term_cline = ""
        s.term_cfill = ""
        #
        s.filebase   = ""
        # DEFINE
        if procD!=0 :
            s.proc_name = procD[0]
            s.proc_titl = procD[1]
            pass
        if plotD!=0 :
            s.plot_name = plotD[0]
            s.plot_titl = plotD[1]
            s.plot_maxx = plotD[2]
            pass
        if termD!=0 :
            s.term_name  = termD[0]
            s.term_title = termD[1]
            s.term_cline = pl.ColorHTML(termD[2])
            s.term_cfill = pl.ColorHTML(termD[2]+"88")
            pass
        if termD!=0 :
            s.filebase   = fileD
            pass
        s.name = "_".join([
                s.proc_name,
                s.term_name,
                s.plot_name
                ])
        pass

class makeUncInfo():
    def __init__(s,uncD):
        s.name = uncD[0]
        s.titl = uncD[1]
        s.varN = [ str(x) for x in uncD[2] ]
        s.types = uncD[3].split(",")
        s.cfill = pl.ColorHTML(uncD[5]+"88")
        s.cline = pl.ColorHTML(uncD[5])
    pass

class  makeUncDescrInfo:
    def __init__(s,uncDList):
        s.Get=dict()
        for uncD in uncDList :
            info=makeUncInfo(uncD)
            s.Get[info.name] = info
            pass
        s.size = len(uncDList)
        pass

    def __getitem__(s,name):
        return s.Get[name]

    def getNameByVarN(s,varN):
        for name,uncI in s.Get.iteritems():
            if varN in uncI.varN :
                return name
        return ""

    def iterkeys(s):
        return s.Get.iteritems()


class TheoUncStudy:
    def __init__(s):
        s.processesDesc = [
                [ "wp" , "W^{+}#rightarrowl^{+}#nu" ],
                [ "wm" , "W^{-}#rightarrowl^{-}#nu" ],
                [ "z0" , "Z#rightarrowll"           ],
                ]
        s.plotsDesc = [
                [  "h_qt"         , "q_{T}[GeV];#sigma[pb]" , 100  ],
                #[  "h_y"          , "y;#sigma[pb]"          , 5   ],
                #[ "h_qtVy"       , "q_{T}[GeV];y"          , 1e8 ],
                #[ "p_qtVy_A4"    , "q_{T}[GeV];y;A4"       , 1e8 ],
                #[  "p_qtVy_A4_prfx" , "q_{T}[GeV];A4"         , 40  ],
                #[  "p_qtVy_A4_prfy" , "y;A4"                  , 3   ],
                #[  "p_qtVy_A4_outliers_px" , "q_{T}[GeV];A4"         , 40  ],
                #[  "p_qtVy_A4_outliers_py" , "y;A4"                  , 3   ],
                ]
        s.uncDescr =  [
                ["stat"  , "stat."           , [0]         , "error"       , "#BEC4C9", "#88919A" ],
                ["alphas", "#alpha_{S} var." , [54,53]     , "pos,neg,sym" , "#FFAB91", "#FF6737" ],
                ["gpar"  , "g var."          , [52,51]     , "pos,neg,sym" , "#D1FF91", "#ABFE37" ],
                #["pdf"   , "PDF"             , range(1,51) , "pos,neg,sym" , "#99CDFF", "#45A2FC" ],
                ]
        s.infoUnc = makeUncDescrInfo(s.uncDescr)
        s.termDesc = [
                [ "TOT"  , ""         , "#88919a" ] ,
                [ "FIN"  , "(fin.)"   , "#d95f02" ] ,
                [ "RES"  , "(res.)"   , "#7570b3" ] ,
                [ "CT"   , "(ct.)"    , "#e7298a" ] ,
                [ "REAL" , "(real.)"  , "#66a61e" ] ,
                [ "VIRT" , "(virt.)"  , "#e6ab02" ] ,
                ]
        s.file_template = "results_merge/grid_151123/dyturbo_{}_lhc7_WZZPT-CT10_{}_v1447428851qt0100y05t{}_outliers.root"
        s.MaxUnc=0.015
        pl.MakeNiceGradient()
        pass

    def GetPlot(s,info,var,name=""):
        # if tot or fin load smartly :
        #
        h=0
        match={
             #"TOT"  : ["RES","CT","REAL","VIRT"] ,
             #"FIN"  : ["CT","REAL","VIRT"]       ,
             "TOT"  : ["TOT"]  ,
             "FIN"  : ["FIN"]  ,
             "RES"  : ["RES"]  ,
             "CT"   : ["CT"]   ,
             "REAL" : ["REAL"] ,
             "VIRT" : ["VIRT"] ,
            }
        for term in match[info.term_name] :
            tinfo=deepcopy(info)
            tinfo.term_name=term
            hist = s.getplot(tinfo,var,term)
            if not "TH1D" in hist.ClassName() or not "TH2D" in hist.ClassName() :
                h = hist.Clone(name)
                break;
            if h == 0 :
                h = pl.EmptyClone(hist,name)
            xsec =  s.getplot(tinfo,"qt_y_total","xsec")
            c = xsec.Integral()  / hist.Integral() 
            h.Add(hist,c)
        return h

    def getplot(s,info,var,name=""):
        if name == "" : name=s.infoUnc.getNameByVarN(var)
        fname=info.filebase.format(info.proc_name,var,info.term_name)
        hname=info.plot_name
        title=name+";"+info.plot_titl
        h = 0
        try :
           h = pl.GetHistSetName(name,fname,hname)
           pass
        except ValueError as detail:
           #print detail
           # try to remove last term of name and then make projection
           splitname=hname.split("_")
           hname = "_".join( splitname[:-1] )
           proj=splitname[-1]
           h2d = pl.GetHistSetName(name+"2d",fname,hname)
           h = pl.GetProjection(h2d,"_"+proj)
       # set axis
        h.SetTitle(title)
        return h

    def getallbands(s,info,central):
        bands=dict()
        for u_name,uInfo in s.infoUnc.iterkeys() :
            # load all variations
            htmp = [ s.GetPlot(info,var,uInfo.name+str(i)) for i,var in enumerate(uInfo.varN)]
            htmp = [central] + htmp
            for type in uInfo.types :
                bname = "_".join([uInfo.name,type])
                bands[bname]=pl.MakeUncBand(bname,htmp,band=type,rel=True)
                # add decomposition of statistical uncertainty (band per each term relative to total)
                if info.term_name=="TOT" and type == "error" :
                    alltrms= [ "FIN" , "RES"  , "CT"   , "REAL" , "VIRT" ]
                    for trm in alltrms :
                        centr = htmp[0]
                        infocp = deepcopy(info)
                        infocp.term_name=trm
                        term_hist = s.GetPlot(infocp,"0",trm)
                        bname="stat_"+trm
                        bands[bname] = pl.MakeUncBand(bname,[centr,term_hist],band="errOth",rel=True)
                #bands[bname].Print()
                pass
            pass
        # create stack
        stacks = [
            ["alphas","gpar"],
            ["stat", "alphas"],
            ["stat","alphas+gpar"],
            #["stat+alphas+gpar","pdf"],
            ]
        bands["stat_sym"] = bands["stat_error"]
        if info.term_name == "TOT": bands["stat_TOT"] = bands["stat_error"]
        type="sym"
        for st1,st2 in stacks :
            res="+".join([st1,st2])
            resName = "_".join([res , type])
            st1name = "_".join([st1 , type])
            st2name = "_".join([st2 , type])
            htmp = [bands[st1name], bands[st2name] ]
            CorMatrix = [[0]] if "alphas" == st1 and "gpar" == st2 else [[0.20]]
            bands[resName]=pl.CombUncBand(resName,htmp,correl=CorMatrix, opts="rel" )
            #bands[resName].Print()
        return bands

    def PlotCentralWithBand(s,ctrl,allbands,info):
        minx=0
        maxx=info.plot_maxx
        if info.term_name=="REAL"  :
            minx=50
            maxx=80
        name = info.name+"_CentBand"
        central=ctrl.Clone("hh")
        central.SetLineWidth(2)
        central.SetLineColor(kAzure)
        central_totband = pl.MakeBandGraph( "centralband", central, [allbands["stat+alphas+gpar_sym"]], "band,rel" )
        central_totband.SetTitle(" Full unc. band " + info.term_title)
        central_totband.SetFillColor(kAzure-9)
        for i,bin in enumerate(central):
            central.SetBinError(i,0);
        # plot
        pl.NewCanvas(name)
        hlist=[central,central_totband]
        pl.SetFrameStyle1D(hlist,minX=minx, maxX=maxx)
        central_totband.Draw("SAMEE3")
        central.Draw("LSAME")
        pl.DrawLegend([central_totband],"fl",legx=0.65,legy=0.80,scale=1.1)
        pl.WriteText(info.proc_titl,0.68,0.85,tsize=0.05)
        pl.Save()
        pass

    def PlotStackUncertainty(s,allbands,info):
        minx=0
        maxx=info.plot_maxx
        uncmax= s.MaxUnc*2 
        forcrng=True
        if info.term_name=="REAL"  :
            minx=50
            maxx=80
            uncmax=1e8
            forcrng=False
        # systematic + statistic
        name = info.name+"_UncStack"
        stackList=[
                [ "stat_sym"             , "stat"   ] ,
                [ "stat+alphas_sym"      , "alphas" ] ,
                [ "stat+alphas+gpar_sym" , "gpar"   ] ,
                #[ "alphas+gpar_sym" , "gpar"   ] ,
                #[ "alphas_sym" , "alphas"   ] ,
                #[ "alphas+gpar_sym" , "alphas"   ] ,
                #[ "gpar_sym" , "gpar"   ] ,
                #[ "gpar_sym" , "gpar"   ] ,
                #[ "stat+alphas+gpar_sym" , "gpar"   ] ,
                ]
        hlist=list()
        titl=""
        for bname,uncName in stackList :
            uncD = s.infoUnc[uncName]
            hlist.append(allbands[bname])
            hlist[-1].SetFillColor(uncD.cfill)
            hlist[-1].SetLineColor(uncD.cline)
            hlist[-1].SetLineWidth(2)
            titl += uncD.titl
            hlist[-1].SetTitle(titl+" "+info.term_title)
            titl += "#oplus"
            #hlist[-1].Print()
            hlist[-1].GetYaxis().SetTitle("rel. unc.")
        pl.NewCanvas(name)
        pl.SetFrameStyle1D(hlist,minX=minx, maxX=maxx,maxY=uncmax,forceRange=forcrng)
        #pl.DrawHistCompare([central_totband,central])
        for h in reversed(hlist) :
            h.Draw("same,f")
        pl.DrawLegend(hlist,"f",legx=0.2,scale=0.8)
        pl.WriteText(info.proc_titl,0.8,0.8,tsize=0.06)
        pl.Save()
        #
        # statistic per term
        if info.term_name == "TOT" :
            uncList=[
                    [ "TOT"    , "FIN"    , "RES"                                     ],
                    [            "FIN"    ,            "CT"     , "REAL"   , "VIRT"   ],
                    [                       "RES"    , "CT"     , "REAL"   , "VIRT"   ],
                    [ "TOT"    , "FIN"    , "RES"    , "CT"     , "REAL"   , "VIRT"   ],
                    ]
            for i,termList in enumerate(uncList) :
                name = info.name+"_UncStatPerTerm"+str(i)
                hlist=list()
                for trm in termList:
                    #uncD = s.infoUnc[uncName]
                    bname = "stat_"+trm
                    hlist.append(allbands[bname])
                    trmD = [ x for x in s.termDesc if x[0]==trm][0]
                    print trmD
                    tmp_info = makeInfo(0,0,trmD,0)
                    hlist[-1].SetFillColor(tmp_info.term_cfill)
                    hlist[-1].SetLineColor(tmp_info.term_cline)
                    hlist[-1].SetLineWidth(2)
                    hlist[-1].SetTitle("stat "+tmp_info.term_title)
                    #hlist[-1].Print()
                    hlist[-1].GetYaxis().SetTitle("rel. unc. wrt. total")
                    pass
                pl.NewCanvas(name)
                pl.SetFrameStyle1D(hlist,maxX=info.plot_maxx,maxY=s.MaxUnc*2,forceRange=True)
                #pl.DrawHistCompare([central_totband,central])
                for h in hlist :
                    h.Draw("same,f")
                pl.DrawLegend(hlist,"f",legx=0.2,scale=0.8)
                pl.WriteText(info.proc_titl,0.8,0.8,tsize=0.06)
                pl.Save()
                pass
            pass
        pass

    def PlotPosNegEnvelopes(s,allbands,info):
        tmp=pl.canvasSettings
        pl.canvasSettings = [0,0,1200,400]
        uncList=[
                "gpar"   ,
                "alphas" ,
                "stat"   ,
                ]
        for uncName in uncList :
            maxiy=s.MaxUnc
            miniy=-maxiy
            hlist=list()
            btypes=[ "pos", "neg" ]
            if "stat" in uncName : 
                btypes = ["sym"]
                miniy=0
            for btype in btypes:
                uncD = s.infoUnc[uncName]
                bname="_".join([uncName,btype])
                hlist.append(allbands[bname])
                hlist[-1].SetFillColor(uncD.cfill)
                hlist[-1].SetLineColor(uncD.cline)
                hlist[-1].SetLineWidth(2)
                hlist[-1].SetTitle(uncD.titl+" "+info.term_title)
                hlist[-1].GetYaxis().SetTitle("rel. unc.")
            name = info.name+"_UncEnv_" + uncName
            pl.NewCanvas(name)
            pl.SetFrameStyle1D(hlist,maxX=info.plot_maxx,scale=2,minY=miniy,maxY=maxiy,forceRange=True)
            pl.axis.GetYaxis().SetNdivisions(5)
            gPad.SetLeftMargin  ( 0.120 )
            gPad.SetTopMargin  ( 0.050 )
            #pl.axis.GetListOfPrimitives().Print()
            #pl.DrawHistCompare([central_totband,central])
            for h in hlist :
                h.Draw("same")
            pl.DrawLegend(hlist[0:1],"f",legx=0.15,legy=0.87)
            pl.WriteText(info.proc_titl,0.8,0.8,tsize=0.09)
            pl.Save()
            pass
        #
        # decompose statistic per term
        if info.term_name == "TOT" :
            uncList=[
                    "TOT"  ,
                    "FIN"  ,
                    "RES"  ,
                    "CT"   ,
                    "REAL" ,
                    "VIRT" ,
                    ]
            for trm in uncList :
                hlist=list()
                maxiy=s.MaxUnc
                miniy=0
                bname="stat_"+trm
                hlist.append(allbands[bname])
                trmD = [ x for x in s.termDesc if x[0]==trm][0]
                tmp_info = makeInfo(0,0,trmD)
                hlist[-1].SetFillColor(tmp_info.term_cfill)
                hlist[-1].SetLineColor(tmp_info.term_cline)
                hlist[-1].SetLineWidth(2)
                hlist[-1].SetTitle("stat "+tmp_info.term_title)
                hlist[-1].GetYaxis().SetTitle("rel. unc. wrt. total")
                name = info.name+"_UncEnvStatTerm_" + trm
                pl.NewCanvas(name)
                pl.SetFrameStyle1D(hlist,maxX=info.plot_maxx,scale=2,minY=miniy,maxY=maxiy,forceRange=True)
                pl.axis.GetYaxis().SetNdivisions(5)
                gPad.SetLeftMargin  ( 0.120 )
                gPad.SetTopMargin  ( 0.050 )
                #pl.axis.GetListOfPrimitives().Print()
                #pl.DrawHistCompare([central_totband,central])
                for h in hlist :
                    h.Draw("same")
                pl.DrawLegend(hlist[0:1],"f",legx=0.15,legy=0.87)
                pl.WriteText(info.proc_titl,0.8,0.8,tsize=0.09)
                pl.Save()
                pass
        pl.canvasSettings=tmp
        pass

    def PlotCentral2D(s,central,info):
        name=info.name+"Cent2D"
        pl.NewCanvas(name)
        pl.SetFrameStyle2D([central])
        central.Draw("same,colz")
        central.SetContour(99)
        pl.WriteText("shape "+info.term_name,0.15,0.87,tsize=0.035,tcol=pl.ColorHTML("#e0e0e0"))
        pl.WriteText(info.proc_titl,0.68,0.85,tsize=0.05,tcol=pl.ColorHTML("#e0e0e0"))
        pl.Save()
        pass

    def PlotUnc2D(s,central,allbands,info):
        blist = [ ]
        tlist = [ "TOT" , "FIN" , "RES" , "CT" , "REAL", "VIRT" ]
        if info.term_name == "TOT" : blist += [ "stat_"+x for x in tlist ]
        for i,bname in enumerate(blist) :
            name=info.name+"_Unc2D_"+bname
            termtit = [x for x in s.termDesc if x[0]==tlist[i]][0][1]
            pl.NewCanvas(name)
            pl.SetFrameStyle2D([allbands[bname]])
            allbands[bname].Draw("same,colz")
            allbands[bname].GetZaxis().SetRangeUser(0,s.MaxUnc*4)
            allbands[bname].SetContour(99)
            pl.WriteText("stat unc "+ termtit + " rel. wrt. total",0.15,0.87,tsize=0.035,tcol=pl.ColorHTML("#303030"))
            pl.WriteText(info.proc_titl,0.68,0.688,tsize=0.05,tcol=pl.ColorHTML("#303030"))
            pl.Save()
        pass


    def DoStudy(s):
        for proc in s.processesDesc :
            for term in s.termDesc :
                for plot in s.plotsDesc :
                    info = makeInfo(proc , plot , term , s.file_template)
                    # get central
                    central = s.GetPlot(info,"0","central")
                    dim = central.GetDimension()
                    # create uncertainty bands
                    allbands = s.getallbands(info,central)
                    if dim == 1 :
                        # create central plot with total unc band
                        s.PlotCentralWithBand(central,allbands, info )
                        # create stack plot with correct uncertainty combination
                        s.PlotStackUncertainty(allbands, info )
                        # create pos-neg envelope plot
                        s.PlotPosNegEnvelopes(allbands, info )
                        # @todo: PDF
                        pass
                    elif dim == 2 :
                        # create central plot
                        s.PlotCentral2D(central,info)
                        # create total unc plot with correct uncertainty combination
                        s.PlotUnc2D(central,allbands,info)
                        # create sym unc per each
                        pass
                    else :
                        raise ImplementationError("Plots for dim {}.".format(dim))
                    pass
                pass
            pass
        pl.MakePreviewFromList(0,"unc_study")
        pass

    def DoPDFQuadStudy(s):
        fname1 = "results_merge/quad_151210/dyturbo_wm_lhc7_WZZPT-CT10_{}_qt0100y05t{}_seed_merge.root"
        fname2 = "results_merge/grid_151201/dyturbo_wm_lhc7_WZZPT-CT10_{}_v1447428851qt0100y05t{}_outliers.root"
        fname3 = "results_merge/allPDF_20151216/dyturbo_wm_lhc7_WZZPT-CT10_all_qt0100y05t{}_seed_1010.root"
        hname1 = "qt_y_total"
        hname2 = "h_qtVy"
        # get
        centrals = [
            pl.GetHistSetTitNam("RES3D" , fname1.format(0, "RES3D") , hname1),
            pl.GetHistSetTitNam("CT3D"  , fname1.format(0, "CT3D")  , hname1),
            pl.GetHistSetTitNam("REAL"  , fname2.format(0, "REAL")  , hname2),
            pl.GetHistSetTitNam("VIRT"  , fname2.format(0, "VIRT")  , hname2),
            #pl.GetHistSetTitNam("REALall", fname3.format(   "REAL") , hname2),
            #pl.GetHistSetTitNam("VIRTall", fname3.format(   "VIRT") , hname2),
            ]
        central = centrals[1].Clone("cent"); centTerms=centrals[1].GetName();
        #central.Add(centrals[1]); centTerms+=centrals[1].GetName();
        #central.Add(centrals[2]); centTerms+=centrals[2].GetName();
        #central.Add(centrals[3]); centTerms+=centrals[3].GetName();
        central.Print()
        # do it for total
        pdfvars=list()
        for i in range(1,3) :
            termvars = [
                pl.GetHistSetTitNam("RES3D" +str(i), fname1.format(i, "RES3D" ) , hname1        ),
                pl.GetHistSetTitNam("CT3D"  +str(i), fname1.format(i, "CT3D"  ) , hname1        ),
                pl.GetHistSetTitNam("REAL"  +str(i), fname3.format(   "REAL"  ) , hname2+str(i) ),
                pl.GetHistSetTitNam("VIRT"  +str(i), fname3.format(   "VIRT"  ) , hname2+str(i) ),
            ]
            pdfvar=termvars[1];
            #pdfvar.Add(termvars[1])
            #pdfvar.Add(termvars[2])
            #pdfvar.Add(termvars[3])
            pdfvars.append(pdfvar)
            pdfvars[-1].Print()
        pdfTerms=""
        #pdfTerms+="RES3D"
        pdfTerms+="CT3D"
        #pdfTerms+="REAl"
        #pdfTerms+="VIRT"
        # prepare
        band=pl.MakeUncBand("pdf",[central]+pdfvars,band="pdf",rel=True)
        band.Print()
        band.SetMaximum(1e-1)
        zoom=""
        zoom="_zoom2"
        # plot
        pl.NewCanvas("QuadPDF_"+hname1+"_pdf"+pdfTerms+"_rel"+centTerms+zoom)
        pl.SetFrameStyle2D([band],maxX=10)
        band.Draw("same,colz")
        #central.Draw("same,colz")
        pl.Save()
        pl.MakePreviewFromFolder(pl.imgDir+"/QuadPDF_*.pdf")
        pass

def makeStatPlot():
    name="res_stat_quad"
    file="results_merge/quad_151207/dyturbo_wm_lhc7_WZZPT-CT10_0_qt0100y05tRES3D_seed_merge.root"
    hfin="qt_y_total"
    hin = pl.GetHistSetTitNam(name,file,hfin)
    h = pl.MakeUncBand("stat",[hin],band="error",rel=True)
    h.SetTitle("stat rel unc")
    pl.MakeNiceGradient()
    pl.NewCanvas(name)
    pl.SetFrameStyle2D([h])
    h.SetMaximum(0.001)
    h.Draw("COLZ")
    pl.Save()
    pass

def plot_profile():
    fname = "run_dir_133/results_merge.root"
    pname = "p_qtVy_A0"
    pO = pl.GetHist(fname,pname)
    p = pl.SwitchTH2Axes(pO,TProfile2D)
    #p.RebinX(5)
    hxy = pl.GetProjection(p   ,ax="_pxy"  )
    px  = pl.GetProjection(p   ,ax="_prfx" )
    py  = pl.GetProjection(p   ,ax="_prfy" )
    hx  = pl.GetProjection(hxy ,ax="_px"   )
    hy  = pl.GetProjection(hxy ,ax="_py"   )
    #
    for h in [p,hxy,px,py,hx,hy]:
        MINY=0.
        MAXY=2.
        dim = h.GetDimension()
        if dim == 2 :
            MINY="-inf"
            MAXY="inf"
            pass
        pl.NewCanvas(h.GetName())
        pl.SetFrameStyle(h,minY=MINY,maxY=MAXY)
        h.Draw("COLZ,SAME")
        pl.Save()
        pass
    #
    pl.MakePreviewFromList(0,"profiles")
    pass

def benchmark():
    samples=[
            #[ "DYTURBO-v0.9.6"      , "results_merge/benchmark_v0_160125"   , "dyturbo_{}_lhc7_CT10nnlo_0_bm0qt0100ym55t{}_seed_outliers.root"   , ["REAL","VIRT","CT","RES"] , "h_qt" , 1. ],
            #[ "DYTURBO-v0.9.6.1" , "results_merge/benchmark_v0.1_160125" , "dyturbo_{}_lhc7_CT10nnlo_0_bm0.1qt0100ym55t{}_seed_outliers.root" , ["REAL","VIRT","CT","RES"] , "h_qt" , 1. ],
            #[ "DYRES-v1.0"          , "results_merge/Stefano_dyturbo_v1"    , "{}{}.root"                                                        , ["r","v"]                  , "pt"   , 1  ],
            #[ "DYTURBO-v0.9.6_PDF" , "results_merge/benchmark_v1_160125" , "dyturbo_{}_lhc7_CT10nnlo_0_bm1qt0100ym55t{}_seed_outliers.root" , ["REAL","VIRT","CT","RES"] , "h_qt" , 1. ],

            [ "DYTURBO-v0.9.6.1"   , "results_merge/benchmark_v0.1_160125" , "dyturbo_{}_lhc7_CT10nnlo_0_bm0.1qt0100ym55t{}_seed_outliers.root" , ["REAL","VIRT","CT","RES"] , "h_qtVy"   , 1. ],
            [ "DYRES-v1.0"         , "results_merge/Stefano_dyturbo_v1"    , "{}{}.root"                                                        , ["r","v"]                  , "yvspt"  , 1  ],
            #[ "DYTURBO-v0.9.6"     , "results_merge/benchmark_v0_160125"   , "dyturbo_{}_lhc7_CT10nnlo_0_bm0qt0100ym55t{}_seed_outliers.root"   , ["REAL","VIRT","CT","RES"] , "h_qtVy" , 1. ],
            #[ "DYTURBO-v0.9.6_PDF" , "results_merge/benchmark_v1_160125"   , "dyturbo_{}_lhc7_CT10nnlo_0_bm1qt0100ym55t{}_seed_outliers.root"   , ["REAL","VIRT","CT","RES"] , "h_qtVy" , 1. ],

            #[ "DYRES-v1.0Res"         , "results_merge/Stefano_dyturbo_v1"  , "{}{}.root"                                                      , ["v"]                      , "pt"     , 1       ],
            #[ "DYTURBO-v0.9.6Res"     , "results_merge/benchmark_v0_160125"   , "dyturbo_{}_lhc7_CT10nnlo_0_bm0qt0100ym55t{}_seed_outliers.root" , ["RES"] , "h_qt" , 1. ],
            #[ "DYTURBO-v0.9.6.1Res"   , "results_merge/benchmark_v0.1_160125" , "dyturbo_{}_lhc7_CT10nnlo_0_bm0.1qt0100ym55t{}_seed_outliers.root" , ["RES"] , "h_qt" , 1. ],
            #[ "DYTURBO-v0.9.6_PDFRes" , "results_merge/benchmark_v1_160125"   , "dyturbo_{}_lhc7_CT10nnlo_0_bm1qt0100ym55t{}_seed_outliers.root" , ["RES"] , "h_qt" , 1. ],

            #[ "DYRES-v1.0Res2D"       , "results_merge/Stefano_dyturbo_v1"  , "{}{}.root"                                                      , ["v"]                      , "yvspt"  , 1       ],
            #[  "DYTURBO-v0.9.6Res2D"   , "results_merge/benchmark_v0_160125" , "dyturbo_{}_lhc7_CT10nnlo_0_bm0qt0100ym55t{}_seed_outliers.root" , ["RES"]                    , "h_qtVy" , 1.      ],
            #[  "DYTURBO-v0.9.6_PDFRes" , "results_merge/benchmark_v1_160125" , "dyturbo_{}_lhc7_CT10nnlo_0_bm1qt0100ym55t{}_seed_outliers.root" , ["RES"]                    , "h_qtVy" , 1.      ],

            #[ "DYRES-v1.0Fin"         , "results_merge/Stefano_dyturbo_v1"  , "{}{}.root"                                                      , ["r"]                , "pt"   , 1  ],
            #[ "DYTURBO-v0.9.6Fin"     , "results_merge/benchmark_v0_160125" , "dyturbo_{}_lhc7_CT10nnlo_0_bm0qt0100ym55t{}_seed_outliers.root" , ["REAL","VIRT","CT"] , "h_qt" , 1. ],
            #[ "DYTURBO-v0.9.6.1Fin" , "results_merge/benchmark_v0.1_160125" , "dyturbo_{}_lhc7_CT10nnlo_0_bm0.1qt0100ym55t{}_seed_outliers.root" , ["REAL","VIRT","CT","RES"] , "h_qt" , 1. ],
            #[ "DYTURBO-v0.9.6_PDFFin" , "results_merge/benchmark_v1_160125" , "dyturbo_{}_lhc7_CT10nnlo_0_bm1qt0100ym55t{}_seed_outliers.root" , ["REAL","VIRT","CT"] , "h_qt" , 1. ],

            #[ "DYRES-v1.0Fin"         , "results_merge/Stefano_dyturbo_v1"  , "{}{}.root"                                                      , ["r"]                , "yvspt"   , 1  ],
            #[ "DYTURBO-v0.9.6Fin"     , "results_merge/benchmark_v0_160125" , "dyturbo_{}_lhc7_CT10nnlo_0_bm0qt0100ym55t{}_seed_outliers.root" , ["REAL","VIRT","CT"] , "h_qtVy" , 1. ],
            #[ "DYTURBO-v0.9.6_PDFFin" , "results_merge/benchmark_v1_160125" , "dyturbo_{}_lhc7_CT10nnlo_0_bm1qt0100ym55t{}_seed_outliers.root" , ["REAL","VIRT","CT"] , "h_qtVy" , 1. ],

            #  merge
            #[ "DYRES-v1.0Fin"         , "results_merge/Stefano_dyturbo_v1"  , "{}{}.root"                                                      , ["r"]                      , "pt"     , 1       ],
            #[ "DYTURBO-v0.9.6Fin"     , "results_merge/benchmark_v0_160125" , "dyturbo_{}_lhc7_CT10nnlo_0_bm0qt0100ym55t{}_seed_merge.root"    , ["REAL","VIRT","CT"]       , "h_qt"   , 1./101. ],
            #[ "DYRES-v1.0Res"         , "results_merge/Stefano_dyturbo_v1"  , "{}{}.root"                                                      , ["v"]                      , "pt"     , 1       ],
            #[ "DYTURBO-v0.9.6Res"     , "results_merge/benchmark_v0_160125" , "dyturbo_{}_lhc7_CT10nnlo_0_bm0qt0100ym55t{}_seed_merge.root"    , ["RES"]                    , "h_qt"   , 1./101. ],
            #[ "DYRES-v1.0"            , "results_merge/Stefano_dyturbo_v1"  , "{}{}.root"                                                      , ["r","v"]                  , "pt"     , 1       ],
            #[ "DYTURBO-v0.9.6"        , "results_merge/benchmark_v0_160125" , "dyturbo_{}_lhc7_CT10nnlo_0_bm0qt0100ym55t{}_seed_merge.root"    , ["REAL","VIRT","CT","RES"] , "h_qt"   , 1./101. ],
        ]
    #
    #blueH   = pl.ColorHTML("45A2FC") #"99CDFF")
    #violetH = pl.ColorHTML("FD37E4") #"FF91F1")
    #greenH  = pl.ColorHTML("D6FF37") #"E9FF91")
    #orangeH = pl.ColorHTML("FFAE37") #"FFD291")
    #
    redH     = pl.ColorHTML("FF6737")
    greenH   = pl.ColorHTML("ABFE37")
    blueH    = pl.ColorHTML("45A2FC")
    torquaH  = pl.ColorHTML("37FD98")
    CONF = [
         [   "bm0",
             [ # samples
                 #[ "DYTURBO (pol. inter. PDF)" ,   "results_merge/benchmark_v0.1_160125" , "dyturbo_{}_lhc7_CT10nnlo_0_bm0.1qt0100ym55t{}_seed_outliers.root" , "h_qt" , 6],
                 #[ "DYTURBO (pol. inter. PDF) average" , "results_merge/benchmark_v0.2_160129" , "dyturbo_{}_lhc7_CT10nnlo_0_bm0qt0100ym55t{}_seed_outliers.root" , "h_qt_rebin_average" , 6],
                 #[ "DYTURBO (pol. inter. PDF)"         , "results_merge/benchmark_v0.2_160129" , "dyturbo_{}_lhc7_CT10nnlo_0_bm0qt0100ym55t{}_seed_outliers.root" , "h_qt_rebin" , 6],
                 #[ "DYRES-v1.0 average"                , "results_merge/Stefano_dyturbo_v1_160201"    , "{}{}.root"                                                        , "pt_rebin_average"   , torquaH],
                 #[ "DYRES-v1.0"                        , "results_merge/Stefano_dyturbo_v1_160201"    , "{}{}.root"                                                        , "pt_rebin"   , 4],

                 #[ "DYTURBO (pol. inter. PDF) average" , "results_merge/benchmark_v0.2_160201_o2" , "dyturbo_{}_lhc7_CT10nnlo_0_bm0qt0100ym55t{}_seed_outliers.root" , "h_qt_rebin_average" , 6],
                 #[ "DYTURBO (pol. inter. PDF)" , "results_merge/benchmark_v0.2_160201_o2" , "dyturbo_{}_lhc7_CT10nnlo_0_bm0qt0100ym55t{}_seed_outliers.root" , "h_qt_rebin" , redH],
                 #[ "DYRES-v1.0 average"                , "results_merge/Stefano_dyturbo_v1_160201_o2"    , "{}{}.root"                                                        , "pt_rebin_average"   , torquaH],
                 [ "DYTURBO (pol. inter. PDF)", "results_merge/benchmark_v0_160204_WZ"      , "dyturbo_{}_lhc7_CT10nnlo_0_bm0qt0100ym55t{}_seed_outliers.root", "h_qt_rebin", redH]   ,
                 [ "DYRES-v1.0"               , "results_merge/Stefano_dyturbo_v1_160201_o2", "{}{}.root"                                                     , "pt_rebin"  , torquaH],
             ],
             [ # terms
                 [ "RES" , ["RES"],                    ["v"]     ],
                 [ "FIN" , ["REAL","VIRT","CT"],       ["r"]     ],
                 [ ""    , ["REAL","VIRT","CT","RES"], ["r","v"] ],

                 #[ "RES"   , ["RES"],    ["RES"]       ],
                 #[ "CT"   , ["CT"],    ["CT"]       ],
                 #[ "REAL" , ["REAL"],    ["REAL"]       ],
                 #[ "VIRT" , ["VIRT"],    ["VIRT"]       ],
                 #[ "FIN" , ["REAL","VIRT","CT"],  ["REAL","VIRT","CT"] ],
                 #[ "" , ["REAL","VIRT","CT","RES"],  ["REAL","VIRT","CT","RES"] ],

                 #[ ""    , ["REAL","VIRT","CT","RES"], ["REAL","VIRT","CT","RES"], ["r","v"], ["r","v"] ],
                 #[ "FIN" , ["REAL","VIRT","CT"],  ["REAL","VIRT","CT"],    ["r"],    ["r"]     ],
                 #[ "RES" , ["RES"], ["RES"],  ["v"],                    ["v"]     ],

                 #[ "RES" , ["v"],                    ["v"]     ],
             ]
         ],
         [   "bm1",
             [ # samples 
                 #[ "DYTURBO (num. integr. PDF)"  , "results_merge/benchmark_v1_160125"   , "dyturbo_{}_lhc7_CT10nnlo_0_bm1qt0100ym55t{}_seed_outliers.root"   , "h_qt" , blueH   ],
                 #[ "DYTURBO (pol. inter. PDF)" , "results_merge/benchmark_v0_160125" , "dyturbo_{}_lhc7_CT10nnlo_0_bm0qt0100ym55t{}_seed_outliers.root" , "h_qt" , redH],
                 #[ "DYTURBO (pol. inter. PDF)"   , "results_merge/benchmark_v0.1_160125" , "dyturbo_{}_lhc7_CT10nnlo_0_bm0.1qt0100ym55t{}_seed_outliers.root" , "h_qt" , redH ],
                 [ "DYTURBO (num. integr. PDF)" , "results_merge/benchmark_v1_160202"      , "dyturbo_{}_lhc7_CT10nnlo_0_bm0qt0100ym55t{}_seed_outliers.root" , "h_qt" , blueH ],
                 [ "DYTURBO (pol. inter.  PDF)" , "results_merge/benchmark_v0.2_160201_o2" , "dyturbo_{}_lhc7_CT10nnlo_0_bm0qt0100ym55t{}_seed_outliers.root" , "h_qt" , redH ],
             ],
             [ # terms 
                 #[ ""    , ["REAL","VIRT","CT","RES"] , ["REAL","VIRT","CT","RES"] ],
                 #[ "FIN" , ["REAL","VIRT","CT"]       , ["REAL","VIRT","CT"]       ],
                 [ "RES" , ["RES"]                    , ["RES"]                    ],
             ]
         ],
          [   "bm2",
              [ # samples
                 [ "DYTURBO (QUADRATURE)"  , "results_merge/benchmark_v2_160202"   , "dyturbo_{}_lhc7_CT10nnlo_0_bm2qt0100ym55t{}_seed_outliers.root"   , "h_qt" , greenH   ],
                 [ "DYTURBO (VEGAS)"  , "results_merge/benchmark_v1_160202"   , "dyturbo_{}_lhc7_CT10nnlo_0_bm0qt0100ym55t{}_seed_outliers.root"   , "h_qt" , blueH   ],
              ],
              [ # terms
                  #[ ""    , ["REAL","VIRT","CT3D","RES3D"] , ["REAL","VIRT","CT","RES"]  ],
                  #[ "FIN" , ["REAL","VIRT","CT3D"]       , ["REAL","VIRT","CT"]        ],
                  [ "RES" , ["RES3D"]                    , ["RES"]                     ],
              ]
          ]
    ]
    #
    procs = ["z0"] # [ "wp", "wm", "z0" ]
    #procs = [ "z0" , "wp", "wm" ]
    projs = ["_px", "_py"] # # ["_px","_py"]
    for cfg in CONF :
        bname   = cfg[0]
        samples = cfg[1]
        terms   = cfg[2]
        for proc in procs :
            proctit = ""
            if "wp" in proc : proctit= "W^{+}"
            if "wm" in proc : proctit= "W^{-}"
            if "z0" in proc : proctit= "Z"
            for term in terms :
                term_name = term[0]
                term_list = term[1:]
                for proj in projs:
                    hists=list()
                    for isampl,sampl_conf in enumerate(samples) :
                        title=sampl_conf[0]
                        hname=sampl_conf[3]
                        scol=sampl_conf[4]
                        # merger creates projection from TH2D
                        if "Vy" in hname or "yvs" in hname : hname+=proj
                        fname=sampl_conf[1]+"/"+sampl_conf[2]
                        # start with empty histogram
                        tmpterm=term_list[isampl][0]
                        h= pl.EmptyClone(pl.GetHistSetTitNam(title,fname.format(proc,tmpterm),hname),title)
                        scale=1
                        for trm in term_list[isampl] :
                            hnametmp=hname
                            if "v" == trm :
                                hnametmp+="_average"
                                #hname+=" average"
                            #if "3D" in trm :
                                #hnametmp=
                            tmp=pl.GetHist(fname.format(proc,trm),hnametmp)
                            c=scale
                            # reweight by X section histogram
                            try:
                                htot=pl.GetHist(fname.format(proc,trm),"qt_y_total")
                                xsec=htot.Integral()
                                intprime=tmp.Integral()
                                inttilde=0
                                intnow=0
                                for ibin,val in enumerate(tmp) :
                                    inttilde+=val
                                    intnow+=val*tmp.GetBinWidth(ibin)
                                #c = inttilde/intnow
                                #if "REAL" in trm : c*=101./1001.
                                print c
                                pass
                            except ValueError:
                                pass
                            #tmp.Scale(c)
                            h.Add(tmp)
                            pass # sum over terms
                        if "Vy" in hname or "yvs" in hname :
                            # do projection by hand
                            #h = pl.GetProjection(h,proj)
                            pass
                        else :
                            #skip for 1D histograms
                            if "_py" in proj : continue
                            pass
                        if "_px" in proj:
                            # its probably pt so rebin and scale by bin width (equidistant)
                            h.Print("range")
                            #h.Rebin(5)
                            #h.Scale(1./h.GetBinWidth(1))
                            h.GetXaxis().SetTitle("p_{T}[GeV]")
                            h.GetYaxis().SetTitle("1/p_{T}d#sigma/dp_{T}[fb.GeV^{-2}]")
                        if "_py" in proj:
                            h.GetXaxis().SetTitle("y")
                            h.GetYaxis().SetTitle("#frac{d#sigma}{dy}[fb]")
                        h.SetLineColor(scol)
                        h.SetLineWidth(2)
                        h.SetMarkerColor(scol)
                        hists.append(h)
                        pass # loop over samples
                    if len(hists)==0: continue
                    pl.CompareHistsInList(bname+"_"+term_name+"_"+proc+proj,hists,compareType="ratio0",legx=0.5,doStyle=False,doSave=False)
                    # cosmetics
                    pl.c1.cd(0)
                    pl.WriteText(proctit+" "+term_name,0.7,0.8,tsize=0.04)
                    pl.c1.cd(2)
                    # find maximal deviation
                    maxdev=0.
                    for bin in list(gPad.GetListOfPrimitives())[3]:
                        if bin == 0: continue
                        maxdev = max([maxdev,abs(1-bin)])
                        print bin, maxdev*100
                    if "FIN" not in term_name : pl.WriteText(" max dev {:.2g}%".format(maxdev*100),0.7,0.8,tsize=0.08)
                    pl.c1.cd()
                    pl.Save()
                    pass # loop over proj
                pass # loop over terms
            pass # loop over processes
        pass # loop over benchmark tests
    pl.MakePreviewFromList(0,"bm_all")
    pass

def plot_y():
    pl.CompareHistsInFiles(
            "proj_y_test",
            [
                [ "noNorm y"            , "rebin_noNorm.root", "h_qt"]              ,
                #[ "noNorm y aver"      , "rebin_noNorm.root", "h_y_average"]      ,
                [ "noNorm proj y"       , "rebin_noNorm.root", "h_qtVy_px"]        ,
                #[ "noNorm proj y errors", "rebin_noNorm.root", "h_qtVy_pu"]        ,
                #[ "noNorm proj y noOver", "rebin_noNorm.root", "h_qtVy_pv"]        ,
                #[ "noNorm proj y aver" , "rebin_noNorm.root", "h_qtVy_py_average"],
                #[ "XNorm y"             , "rebin_XNorm.root" , "h_y"]              ,
                #[ "XNorm y aver"       , "rebin_XNorm.root" , "h_y_average"]      ,
                #[ "XNorm proj y"        , "rebin_XNorm.root" , "h_qtVy_py"]        ,
                #[ "XNorm proj y aver"  , "rebin_XNorm.root" , "h_qtVy_py_average"],
            ],
            compareType="ratio0",
            legx=0.6
            )

## Documentation for main
#
# More details. 
if __name__ == '__main__' :
    #print_results();
    #print_table();
    #check_cancelation()
    #w_pt();
    #w_pt_y();
    #merge_all_hist()
    #plot_pt("results/pt_table_CT10nnlo.txt")
    #quick_calc()
    #root_file_integral()
    #wwidth_table()
    #uncert_as_g()
    #find_fluctuations()
    #DY = TheoUncStudy()
    #DY.DoStudy()
    #DY.DoPDFQuadStudy()
    #plot_profile()
    #benchmark()
    #makeStatPlot()
    plot_y()
    pass



