#!/usr/bin/python
# -*- coding: utf-8 -*- 

## Documentation for file
#
# More details. 
#
# @file pt_plot.py
# @author cuto <Jakub.Cuth@cern.ch>
# @date 2015-08-07

import sys,os,re
# for batch add -b before first ROOT calling 
sys.argv.append('-b') # run in batch

#from PlotTools import *
from ROOT  import TGraphErrors, TFile
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
    filetmpREAL="results_merge/wpm_300seeds_151019/{}_{}_lhc7_{}_0_qt0100y05t{}_seed_outliers.root" #results_merge/wpm_300seed_151019/{}_{}_lhc7_{}_0_qt0100y05t{}_seed_outliers.root" results_merge/wpm_300seeds_151019/dyturbo_wp_lhc7_ZPT-CT10_0_qt0100y05tREAL_seed_outliers.root
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
                hlistProject.append(pl.GetHistSetTitNam("real 300j", "results_merge/wpm_300seeds_151019/dyturbo_wp_lhc7_ZPT-CT10_0_qt0100y05tREAL_seed_outliers.root", hist))
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


def GetValError(fiducial, col, proc):
    file_tmpl="results_merge/dyturbo_{}_{}_CT10nnlo_0_f{}qt01000y-55t{}_100101.root"
    #file_tmpl="results/dyturbo_{}_{}_CT10nnlo_0_f{}qt01000y-55t{}_100101.root"
    file_tmpl="results_merge/wwidth_properD0massCut/dyturbo_{}_{}_CT10nnlo_0_f{}qt01000y-55t{}_100101.root"
    term="TOT"
    hist="qt_y_total"
    totInt,totIne=0,0
    for term in [ "RES", "CT", "REAL", "VIRT" ]:
        h = pl.GetHist(file_tmpl.format(proc,col,fiducial,term),hist)
        integ , inerr = getIntegralError(h)
        #print term, integ, inerr, inerr/integ
        totInt,totIne = addIntegralError(integ,inerr,totInt,totIne)
        pass
    #print "tot ", totInt, totIne, totIne/totInt
    return totInt, totIne/totInt

def print_wwidth_res(int,err):
    print int,err,

def wwidth_table():
    experiments = [
            [ "D0"    , "tev1" ],
            #[ "CDF"   , "tev2" ],
            #[ "ATLAS" , "lhc7" ],
            #[ "CMS"   , "lhc8" ],
        ]
    processes = [ "wp", "wm", "z0"]
    for exp,col in experiments :
        #new header
        print " table for exp", exp, "@", col, " in order W+ W- Z"
        print " full  | fiducial | eff "
        for proc in processes :
            # full
            full_v, full_e_r = GetValError("FULL",col,proc)
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
    root_file_integral()
    #wwidth_table()
    pass



