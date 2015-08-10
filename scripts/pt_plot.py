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


from PlotTools import *
from ROOT  import TGraphErrors

import asciitable
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


## Documentation for main
#
# More details. 
if __name__ == '__main__' :
    plot_pt("results/pt_table_CT10nnlo.txt")
    pass



