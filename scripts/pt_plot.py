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
#import asciitable

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
    pl = PlotTools()
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



## Documentation for main
#
# More details. 
if __name__ == '__main__' :
    #print_results();
    merge_all_hist()
    #plot_pt("results/pt_table_CT10nnlo.txt")
    pass



