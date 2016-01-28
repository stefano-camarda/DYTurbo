#!/usr/bin/python
# -*- coding: utf-8 -*- 

## Documentation for file
#
# More details. 
#
# @file PlotTools.py
# @author cuto <Jakub.Cuth@cern.ch>
# @date 2013-09-13



# sys
# import sys,os
# for batch add -b before first ROOT calling 
# sys.argv.append('-b') # run in batch


import os,time
from ROOT import ROOT, TObject, TTree, TH2F, TH1D, TProfile, TProfile2D, TH1F, TH1I, TH1C, TGraph, TGraphAsymmErrors, TF1, TMath, TFile, TCanvas, TBox, TLegend, TColor, gPad, gStyle, gROOT, Double, TLatex, TMarker, TLine

# very handy for debugging -- missing backtrace
import logging as MSG
FORMAT = '%(levelname)s %(asctime)-15s from %(filename)s:%(lineno)d `%(funcName)s` : %(message)s'
MSG.basicConfig(
        format=FORMAT,
        #filename="last.log",
        level=MSG.DEBUG
        )

def doMSG(lvl):
    return MSG.getLogger().isEnabledFor(lvl)


import colorsys

def log10(fl):
    return TMath.Log10(fl)

def sqrt(fl):
    return TMath.Sqrt(fl)


class HISTDEF:
    def __init__ (s,
            _var = "var"      ,
            _name = "name"    ,
            _nbins = 1        ,
            _xmin = 0.        ,
            _xmax = 1.        ,
            _title = "title"  ,
            _xtitle = "x"     ,
            _weight = "1"
            ) :
        s.var    = _var
        s.name   = _name   
        s.nbins  = _nbins  
        s.xmin   = _xmin   
        s.xmax   = _xmax   
        s.title  = _title  
        s.xtitle = _xtitle 
        s.weight = _weight 
        pass

class TreeHist():
    def __init__(s, _filename, _treename="physics") :
        s.filename = _filename
        s.treename = _treename
        s.histoDef = list()
        s.hists = list()
        s.looplimit = 0
        s.doTreeLooping=True
        pass

    def VarList(s):
        return [ h.name for h in s.histoDef]

    def loadtree(s):
        s.inf=TFile.Open(s.filename, "READ")
        s.tree= s.inf.Get(s.treename)

    def makehistos_byDraw(s):
        for hdf in s.histoDef :
            s.tree.Draw(hdf.var.replace("ev.",""),"weight_evt*"+hdf.weight.replace("ev.",""))
            h = s.tree.GetHistogram().Clone(hdf.name)
            #h.SetName(hdf.name)
            h.SetTitle(hdf.title+";"+hdf.xtitle+";entries")
            s.hists.append(h)
        pass

    def makehistos_byLooping(s):
        for hdf in s.histoDef :
            h = TH1F (hdf.name, hdf.title+";"+hdf.xtitle+";entries", hdf.nbins, hdf.xmin, hdf.xmax)
            s.hists.append(h)
            pass

        print "looping the tree {0} from file {1}".format(s.treename,s.filename)
        i_evt = 0
        start = time.time()
        for ev in s.tree :
            i_evt+=1
            if s.looplimit==i_evt : break
            if ( ( TMath.Log( i_evt ) / TMath.Log (2) ) % 1 == 0):
                now = time.time()
                print "evt: ", i_evt, "({0} kHz)".format(i_evt/( 1000*(now-start) ) )
            for hdf in s.histoDef :
                i_h = s.histoDef.index(hdf)
                var = eval (hdf.var)
                weight = eval (hdf.weight)
                s.hists[i_h].Fill(var, weight)
                pass
            pass

    def AddHistoDef(s,
            _var = "var"      ,
            _name = "name"    ,
            _nbins = 1        ,
            _xmin = 0.        ,
            _xmax = 1.        ,
            _title = "title"  ,
            _xtitle = "x"     ,
            _weight = "1"
            ) :
        s.histoDef.append(HISTDEF(
            _var    ,
            _name   ,
            _nbins  ,
            _xmin   ,
            _xmax   ,
            _title  ,
            _xtitle ,
            _weight
            ))
        pass




    def MakeHistos(s):
        s.loadtree()
        if s.doTreeLooping :
            s.makehistos_byLooping()
        else:
            s.makehistos_byDraw()
        pass

    def Save(s,outname="hists.root"):
        of = TFile.Open(outname,"RECREATE")
        of.cd()
        for h in s.hists:
            h.Write()
            pass
        of.Write()
        of.Close()



## Documentation for a class.
#
# More details.
class PlotTools:

    ## The constructor.
    def __init__(s):
        s.rootFile = 0;
        s.imgDir = "share/img/"
        #                  name         HistName, axis prefix, parameter index
        s.GausVarSettings=dict()
        s.GausVarSettings["mean"  ] = [ "Mean"  , "mean"     , 0                ]
        s.GausVarSettings["rms"   ] = [ "RMS"   , "RMS"      , 0                ]
        s.GausVarSettings["mu"    ] = [ "Mu"    , "#mu"      , 1                ]
        s.GausVarSettings["sigma" ] = [ "Sigm"  , "#sigma"   , 2                ]

        #                   x, y, width, heigh
        s.canvasSettings = [0, 0, 800  , 600   ]
        s.c1 = 0;
        s.axis = 0;
        s.latex = TLatex();
        s.boxes = list();
        s.axes = list();
        s.measurements = list();
        s.legend = 0;
        #s.c_white_a10 = s.make_color_transparent(ROOT.kBlue-10,50)
        s.outImgFormats=["png","eps","pdf"]
        s.outImgFormats=["root","pdf","eps"]

        s.forceRange=False

        s.updated_plots=list()

        s.contourHists = list()
        # colors
        s.new_col = list()
        # color set 1
        s.nice2_purple  = s.ColorHTML("8B3A96")
        s.nice2_blue    = s.ColorHTML("405CB5")
        s.nice2_turque  = s.ColorHTML("4BAD8B")
        s.nice2_green   = s.ColorHTML("587F13")
        s.nice2_Dorange = s.ColorHTML("C2942F")
        # color  set 2 -- lighter
        s.nice_purple = s.ColorHTML("7A2396")
        s.nice_blue   = s.ColorHTML("2E7EB0")
        s.nice_green  = s.ColorHTML("41943C")
        s.nice_orange = s.ColorHTML("CCA11A")
        s.nice_red    = s.ColorHTML("BD3529")


        s.l_TFileWithTree = list()
        s.TFileWithTree = 0
        s.l_TTree = list()
        s.TTree = 0

    def NewCanvas(s, name) :
        gROOT.Reset()
        if s.c1 != None and s.c1 != 0 : s.c1.Close();
        s.c1 = TCanvas(name,name.replace("_"," "), s.canvasSettings[0], s.canvasSettings[1], s.canvasSettings[2], s.canvasSettings[3]);
        s.axes[:]=[]
        s.measurements[:]=[]
        s.c1.cd()

    def Save(s) :
        s.c1.cd()
        s.axis.Draw("SAME")
        gPad.Update()
        # ensure the img dir exists
        if not os.path.exists(s.imgDir):
            os.makedirs(s.imgDir)
        # save with differnt types
        for imgType in s.outImgFormats :
            gPad.SaveAs(s.imgDir+"/"+s.c1.GetName()+"."+imgType)
        # add name to update list
        s.updated_plots.append(s.imgDir+"/"+s.c1.GetName())
        pass

    def ls(s,filename):
        fsplit =  filename.split(":")
        f = TFile(fsplit[0],"r")
        if len(fsplit)>1 :
            f.cd(fsplit[1])
        f.ls()
        f.Close()
        pass

    def GausTH2toTH1(s, inHist,var="mean", plotSlices=False ):
        debug=False
        nbinx = inHist.GetXaxis().GetNbins();
        nbiny = inHist.GetYaxis().GetNbins();
        var_title = s.GausVarSettings[var][1]+" ("+inHist.GetYaxis().GetTitle()+")"
        #outHist = TH1F("G"+var+"_"+inHist.GetName(), 
                       #";"+inHist.GetXaxis().GetTitle()+";"+var_title, 
                       #nbinx, inHist.GetXaxis().GetXmin(), inHist.GetXaxis().GetXmax()
                       #)
        outHist = inHist.ProjectionX("G"+var+"_"+inHist.GetName(), 1,1)
        outHist.SetTitle(";"+inHist.GetXaxis().GetTitle()+";"+var_title)
        outHist.Reset()
        #sliceHist = TH1F("slice","slice", nbiny, inHist.GetXaxis().GetXmin(), inHist.GetXaxis().GetXmax())
        sliceHist = inHist.ProjectionY("slice", 1,1)
        sliceHist.Reset()
        gausFun = TF1("gausPlusConst","gaus(0)+[3]")
        for i_x in range(1,nbinx+1) :
            integral = 0
            for i_y in range(1,nbiny+1) :
                sliceHist.SetBinContent(i_y, inHist.GetBinContent(i_x,i_y))
                integral = integral + inHist.GetBinContent(i_x,i_y)
            if integral > 50 :
                mean    = sliceHist.GetMean();
                #meanErr = sliceHist.GetMeanError();
                rms     = sliceHist.GetRMS();
                #rmsErr  = sliceHist.GetRMSError();
                #integral = sliceHist.Integrate("width")
                gausFun.SetRange(mean-2.0*rms, mean+2.0*rms)
                #gausFun.SetParameter(1, integral)
                gausFun.SetParameter(0, integral)
                gausFun.SetParLimits(0, 0, 10*integral)
                gausFun.SetParameter(1, mean)
                gausFun.SetParLimits(1, mean-3.*rms, mean+3*rms)
                gausFun.SetParameter(2, rms)
                gausFun.SetParLimits(2, 0., 5*rms)
                #gausFun.SetParameter(3, 0)
                #sliceHist.Fit(gausFun, "WQR")
                sliceHist.Fit(gausFun, "WRQ")
                if gausFun.GetParameter(2) < 0. : 
                    print "SIGMA LOWER THAN 0:",gausFun.GetParameter(2)
                    outHist.SetBinContent(i_x, 0.)
                    outHist.SetBinError  (i_x, 0.)
                if debug :
                    print "{10}:{11}: A {0:f}+-{1:f}, mu {2:f}+-{3:f}, sigma {4:f}+-{5:f}, C {6:f}+-{7:f}, Mean {8:f}, RMS{9:f}".format(
                            gausFun.GetParameter(0),
                            gausFun.GetParError (0),
                            gausFun.GetParameter(1),
                            gausFun.GetParError (1),
                            gausFun.GetParameter(2),
                            gausFun.GetParError (2),
                            gausFun.GetParameter(3),
                            gausFun.GetParError (3),
                            mean,
                            rms,
                            inHist.GetName(),i_x
                            )
                if var=="mean":
                    outHist.SetBinContent(i_x, mean)
                    outHist.SetBinError  (i_x, 0)
                elif var=="rms":
                    outHist.SetBinContent(i_x, rms)
                    outHist.SetBinError  (i_x, 0)
                else :
                    outHist.SetBinContent(i_x, gausFun.GetParameter(s.GausVarSettings[var][2]))
                    outHist.SetBinError  (i_x, gausFun.GetParError (s.GausVarSettings[var][2]))
            else :
                outHist.SetBinContent(i_x, 0.)
                outHist.SetBinError  (i_x, 0.)
                #print " Low statistics", inHist.GetName(), "in x bin", i_x
                pass
        return outHist;

    def EmptyClone(s, h, name):
        o = h.Clone(name)
        o.Reset()
        o.SetDirectory(0)
        return o

    def RatioHist(s, h1, h2 ) :
        h3 = s.EmptyClone(h1, "R_"+h1.GetName()+"__"+h2.GetName())
        h3.GetYaxis().SetTitle(h1.GetTitle()+"/"+h2.GetTitle())
        for i_x in range(1,h3.GetNbinsX()) :
            v2 = h2.GetBinContent(i_x)
            if v2 == 0 : continue;
            v1 = h1.GetBinContent(i_x)
            e1 = h1.GetBinError  (i_x)
            e2 = h2.GetBinError  (i_x)
            h3.SetBinContent(i_x, v1/v2)
            h3.SetBinError  (i_x, TMath.Sqrt(TMath.Power(e1/v1,2) + TMath.Power(e2/v2,2) ))
        return h3

    def DiffHist(s, h1, h2 ) :
        h3 = s.EmptyClone(h1, "D_"+h1.GetName()+"__"+h2.GetName() )
        h3.Reset()
        h3.GetYaxis().SetTitle(h1.GetTitle()+"-"+h2.GetTitle())
        for i_x in range(1,h3.GetNbinsX()) :
            v2 = h2.GetBinContent(i_x)
            if v2 == 0 : continue;
            v1 = h1.GetBinContent(i_x)
            e1 = h1.GetBinError  (i_x)
            e2 = h2.GetBinError  (i_x)
            h3.SetBinContent(i_x, v1-v2)
            h3.SetBinError  (i_x, TMath.Sqrt(TMath.Power(e1*v2,2) + TMath.Power(e2*v1,2) ))
        return h3

    def MakeUncBand(s,name,hlist,**kwargs):
        relative = kwargs["rel"]  if "rel"  in kwargs else True
        utype    = kwargs["band"] if "band" in kwargs else "error"
        h = s.EmptyClone(hlist[0],name)
        dim = h.GetDimension()
        # loop bins
        ist=0
        for xbin in range(0,h.GetNbinsX()+2):
            ybinlist = [0] if dim<2 else range(0,h.GetNbinsY()+2)
            for ybin in ybinlist:
                zbinlist = [0] if dim<3 else range(0,h.GetNbinsZ()+2)
                for zbin in zbinlist:
                    ibin = h.GetBin(xbin,ybin,zbin)
                    err =  hlist[0].GetBinError(ibin)
                    val =  hlist[0].GetBinContent(ibin)
                    unc = 1
                    pos = neg = sym = err_other = 0
                    try :
                        err_other = hlist[1].GetBinError(ibin)
                        pos = hlist[1].GetBinContent(ibin)
                        neg = hlist[2].GetBinContent(ibin)
                        sym= (abs(pos - neg))/2.
                        pos = pos - val
                        neg = neg - val
                        # pdf
                        CL99to68=1/2.705
                        N = len(hlist)-1
                        allind = range(N)
                        pdfVar = [ hlist[i+1].GetBinContent(ibin) for i in allind]
                        pdfVarSq = [ (x-val)**2 for x in pdfVar]
                        pdfDWSq = [ CL99to68*(pdfVarSq[i] + pdfVarSq[i+1])/4. for i in range(N/2)]
                        pdfsigma = sqrt(sum(pdfDWSq))
                        #print "    PDFERROR ", ibin, pdfsigma
                        #print "    ", len(pdfVar), pdfVar
                        #print "    ", len(pdfVarSq), pdfVarSq
                        #print "    ", len(pdfDWSq), pdfDWSq
                        #print "    PDFERROR END"
                        #if ist>100 : return
                        #ist+=1
                    except IndexError:
                        pass
                    if relative and val!=0 : unc/=abs(val)
                    if   utype=="stat"   : unc *= sqrt(abs(val))
                    elif utype=="error"  : unc *= err
                    elif utype=="errOth" : unc *= err_other
                    elif utype=="pos"    : unc *= pos
                    elif utype=="neg"    : unc *= neg
                    elif utype=="sym"    : unc *= sym
                    elif utype=="pdf"    : unc *= pdfsigma
                    else : raise ValueError("Uknown error type {}".format(utype))
                    #print ibin,unc
                    h.SetBinContent(ibin,unc)
                    pass
                pass
            pass
        return h

    def CombUncBand(s,name,hlist,**kwargs):
        CorMatrix = kwargs["correl"] if "correl" in kwargs else [[0]]
        opts      = kwargs["opts"]   if "opts"   in kwargs else "rel"
        h = s.EmptyClone(hlist[0],name)
        dim = h.GetDimension()
        N = len(hlist)
        # loop bins
        for xbin in range(0,h.GetNbinsX()+2):
            ybinlist = [0] if dim<2 else range(0,h.GetNbinsY()+2)
            for ybin in ybinlist:
                zbinlist = [0] if dim<3 else range(0,h.GetNbinsZ()+2)
                for zbin in zbinlist:
                    ibin = h.GetBin(xbin,ybin,zbin)
                    #vals=list()
                    vals = sorted([ x.GetBinContent(ibin) for x in hlist])
                    unc = 0
                    if "rel" in opts :
                        # uncor
                        for v in vals : unc += v*v
                        # cor
                        for i in range (0,N) :
                            for j in range(i+1,N) :
                                unc+= -2.*CorMatrix[i][j-1]*abs(vals[i]*vals[j])
                    else : raise ValueError("Uknown option {}".format(opts))
                    #print ibin, unc, vals
                    #hlist[0].Print()
                    #hlist[1].Print()
                    h.SetBinContent(ibin,sqrt(unc))
                    pass
                pass
            pass
        return h

    def MakeBandGraph(s,name,cent,hlist,opt="diff"):
        #band = TGraphAsymmErrors()
        band = s.EmptyClone(cent,name)
        for i in range(1,cent.GetNbinsX()+1) :
            x  = cent.GetBinCenter (i)
            y  = cent.GetBinContent(i)
            #band.SetPoint(i-1,x,y)
            band.SetBinContent(i,y)
            eyl=eyh=exl=exh=0.
            exl=exh= cent.GetBinWidth(i)/2.
            if "diff" in opt :
                ey1 =  hlist[1].GetBinContent(i)-y
                ey2 = -hlist[2].GetBinContent(i)+y
                eyh=abs(ey1)
                eyl=abs(ey2)
                if ey1*ey2 < 0 and ey1 < 0 :
                    eyl = abs(ey1)
                    eyh = abs(ey2)
                    pass
                if "sym" in opt :
                    eyl=eyh=(eyh+eyl)/2
                pass
            elif "band" in opt :
                ey=hlist[0].GetBinContent(i)
                #ey=cent.GetBinError(i)/y
                if "rel" in opt and y!=0 :
                    ey*=abs(y)
                eyl=eyh=ey
                #print i,eyl,eyh,ey,y,hlist[0].GetBinContent(i)
                pass
            else : raise ValueError("Uknown option {}".format(opt))
            #band.SetPointError(i-1,exl,exh,eyl,eyh)
            band.SetBinError(i,ey)
            pass
        return band

    def GetObject(s,filename,objname,classname) :
        f = TFile(filename, "r")
        if not f.IsOpen() :
            raise ValueError("Can not open file '"+filename+"'")
        h = f.Get(objname)
        if type(h) == TObject : 
            raise ValueError("No "+classname+" "+"'"+objname+"'"+" in file '"+"'"+filename+"'"+"'")
        if classname == "histogram" :
            h.SetDirectory(0)
            if type(h) == TProfile :
                name  = "prf_"+h.GetName()
                title = h.GetTitle()
                title += ";"
                title += h.GetXaxis().GetTitle()
                title += ";"
                title += h.GetYaxis().GetTitle()
                nbins = h.GetNbinsX()
                xmin = h.GetXaxis().GetXmin()
                xmax = h.GetXaxis().GetXmax()
                hh = TH1D(name, title, nbins, xmin, xmax)
                hh.SetDirectory(0)
                for ibin in range(0,h.GetNbinsX()+2):
                    v = h.GetBinContent ( ibin )
                    e = h.GetBinError   ( ibin )
                    hh.SetBinContent ( ibin , v )#if v==v else 0. )
                    hh.SetBinError   ( ibin , e )#if e==e else 0. )
                return hh
        return h

    def GetProjection(s,h,ax="_px") :
        if "TH1" in h.ClassName() : return h.Clone(h.GetName()+ax)
        name=h.GetName()+ax
        ytit  = h.GetZaxis().GetTitle()
        hp=0
        if   ax == "_px": # out of 2D
            n = h.GetNbinsX()
            #hp=h.ProjectionX(name,1,n,"e")
            hp=h.ProjectionX(name) #,1,n,"e")
            hp.GetYaxis().SetTitle(ytit)
            #hp=h.ProfileX(name,1,n,"e")
            #for i in range(n-10,n+1):
                #hp.SetBinContent(i,0)
        elif ax == "_py":
            n = h.GetNbinsY()
            hp=h.ProjectionY(name,1,n,"e")
            hp.GetYaxis().SetTitle(ytit)
            #hp=h.ProfileY(name,1,n-10,"e")
        elif ax == "_pxy":
            hp=h.ProjectionXY(name,"w")
            x_tit = h.GetXaxis().GetTitle()
            y_tit = h.GetYaxis().GetTitle()
            z_tit = h.GetZaxis().GetTitle()
            tit = "{};{};{};{}".format(h.GetTitle(),x_tit,y_tit,z_tit)
            hp.SetTitle(tit)#.GetTitle())
            #hp.GetXaxis().SetTitle(x_tit)
            #hp.GetYaxis().SetTitle(x_tit)
            #hp.GetZaxis().SetTitle(x_tit)
        elif ax == "_prfy":
            n = h.GetNbinsY()
            hp=h.ProfileY(name+"prof",1,n,"w").ProjectionX(name)
            #hp=h.ProjectionX(name,1,n,"w")
            hp.GetYaxis().SetTitle(ytit)
        elif ax == "_prfx":
            n = h.GetNbinsX()
            hp=h.ProfileX(name+"prof",1,n,"w").ProjectionX(name)
            #hp=h.ProjectionX(name,1,n,"w")
            hp.GetYaxis().SetTitle(ytit)
            # correct for nan in error
            #  for xbin in range(1,hp.GetNbinsX()+1):
            #      v1 = hp.GetBinContent (xbin)
            #      e1 = hp.GetBinError   (xbin)
            #      if v1!=v1 : MSG.warning(" nan in value: hist {} bin {}".format(hp.GetName(),xbin))
            #      if e1!=e1 :
            #          MSG.warning(" nan in error: hist {} bin {}".format(hp.GetName(),xbin))
            #          # sum the entries
            #          ybinlist = range(1,h.GetNbinsY()+1)
            #          entries = list()
            #          for ybin in ybinlist:
            #              ibin = h.GetBin(xbin,ybin)
            #              n1 = h.GetBinEntries(ibin)
            #              entries.append(abs(n1))
            #              pass
            #          tot=np.sum(entries)
            #          hp.SetBinError(xbin, sqrt(tot))
            #hp=h.ProfileX(name,1,n,"e")
            #for i in range(n-10,n+1):
                #hp.SetBinContent(i,0)
        else : raise RuntimeError("Uknown type of projection `{}`".format(ax) )
        return hp

    def GetGraph(s, filename, objname) :
        h = s.GetObject(filename,objname,"graph")
        return h

    def GetHist(s, filename, objname) :
        h = s.GetObject(filename,objname,"histogram")
        return h

    def GetFunction(s, filename, objname) :
        h = s.GetObject(filename,objname,"function")
        return h

    def GetCanvas(s, filename, objname) :
        h = s.GetObject(filename,objname,"canvas")
        return h

    def GetHistSetName(s, newname, filename, histname) :
        h = s.GetHist(filename,histname)
        h.SetName(newname)
        #f.Close()
        return h

    def GetHistSetTitle(s, newtitle, filename, histname) :
        h = s.GetHist(filename,histname)
        h.SetTitle(newtitle)
        #f.Close()
        return h

    def GetHistSetTitNam(s, newtitle, filename, histname) :
        h = s.GetHist(filename,histname)
        h.SetName(newtitle)
        h.SetTitle(newtitle)
        #f.Close()
        return h

    def GetTree(s, filename, treename) :
        # get tree
        s.l_TFileWithTree . append( TFile(filename, "r")                )
        s.l_TTree         . append( s.l_TFileWithTree[-1].Get(treename) )
        # set tmp variables
        s.TFileWithTree = s.l_TFileWithTree [-1]
        s.TTree         = s.l_TTree         [-1]
        # test and print
        s.TTree.Print()
        if type(s.TTree) == TObject : 
            print "No tree ", "'"+treename+"'"," in file '", "'"+filename+"'", "'"
            return 0
        #f.Close()
        #return t
        pass

    def GetHistogramFromLastTree(s,name, variable, weight, Nevents=int(1e10)):
        s.TTree.Draw(variable,weight,"",Nevents)
        h = s.TTree.GetHistogram().Clone(name)
        h.SetTitle(name)
        return h

    def SwitchTH2Axes(s,h_in,H2) :
        outname = h_in.GetName()+"_yx"
        z_tit = h_in.GetZaxis().GetTitle()
        x_tit = h_in.GetYaxis().GetTitle()
        x_N   = h_in.GetYaxis().GetNbins()
        x_lo  = h_in.GetYaxis().GetXmin()
        x_hi  = h_in.GetYaxis().GetXmax()
        y_tit = h_in.GetXaxis().GetTitle()
        y_N   = h_in.GetXaxis().GetNbins()
        y_lo  = h_in.GetXaxis().GetXmin()
        y_hi  = h_in.GetXaxis().GetXmax()
        tit = "{};{};{};{}".format(h_in.GetTitle(),x_tit,y_tit,z_tit)
        h_out = H2(outname,tit ,x_N,x_lo,x_hi ,y_N,y_lo,y_hi )
        for xbin in range(0,x_N+2) :
            for ybin in range(0,y_N+2) :
                ibinT = h_in.GetBin(ybin,xbin)
                ibin  = h_out.GetBin(xbin,ybin)
                v = h_in.GetBinContent(ybin,xbin)
                e = h_in.GetBinError  (ybin,xbin)
                n = 0
                if H2 == TProfile2D :
                    n = h_in.GetBinEntries(ibinT)
                    v*=n
                h_out.SetBinContent(xbin,ybin,v)
                h_out.SetBinError  (xbin,ybin,e)
                if H2 == TProfile2D : h_out.SetBinEntries  (ibin,n)
        return h_out


    def GetRangesFromHists(s, hlist, **kwargs): 
        # old parse:

        # supported options: dim=1, logY=False, logX=False, maxY="inf", minY="-inf", maxX="inf", minX="-inf") :
        dim  = kwargs[ "dim"  ] if "dim"  in kwargs else 1
        logX = kwargs[ "logX" ] if "logX" in kwargs else False
        logY = kwargs[ "logY" ] if "logY" in kwargs else False

        minX = kwargs["minX"] if "minX" in kwargs else "-inf"
        maxX = kwargs["maxX"] if "maxX" in kwargs else "inf"
        minY = kwargs["minY"] if "minY" in kwargs else "-inf"
        maxY = kwargs["maxY"] if "maxY" in kwargs else "inf"

        xmin = (float(minX) if minX !="-inf" else hlist[0].GetXaxis().GetXmin())
        xmax = (float(maxX) if maxX !="inf"  else hlist[0].GetXaxis().GetXmax())
        ymin = (float(minY) if minY !="-inf" else hlist[0].GetYaxis().GetXmin())
        ymax = (float(maxY) if minY !="inf"  else hlist[0].GetYaxis().GetXmax())

        if dim==1 :
            # get range
            y_mins = list()
            y_maxs = list()
            #print
            for h in hlist :
                if maxX != "inf"  or minX != "-inf" : h.GetXaxis().SetRangeUser(float(minX), float(maxX))
                if logY :
                    if float(maxY) < 0 : maxY="inf"
                    if float(minY) < 0 : minY=0
                if "TGraph" in h.ClassName() :
                    #h.Draw("APL")
                    #y_mins.append(h.GetMinimum())
                    #y_maxs.append(h.GetMaximum())
                    lower = list()
                    upper = list()
                    for i in range(h.GetN()) :
                        ce = h.GetY()[i]
                        lo = 0
                        hi = 0
                        if "TGraphErrors" == h.ClassName() :
                            lo = hi = h.GetEY()[i]
                        elif "TGraphAsymmErrors" == h.ClassName() :
                            lo = h.GetEYlow ()[i]
                            hi = h.GetEYhigh()[i]
                        if float(minY) < ce - lo : lower .append (ce-lo)
                        if float(maxY) > ce + hi : upper .append (ce+hi)
                        pass
                    y_mins.append( min(lower) )
                    y_maxs.append( max(upper) )
                else :
                    y_mins.append(h.GetMinimum(float(minY)))
                    y_maxs.append(h.GetMaximum(float(maxY)))
                #print y_mins[-1], y_maxs[-1]
            #print min(y_mins) ,  max(y_maxs)
            if logY :
                dy = TMath.Power(10, ( TMath.Log10(max(y_maxs)) - TMath.Log10(min(y_mins)) )*0.1)
                if dy == 0 : 
                    ymin = 1e-1
                    ymin = 1e+1
                else:
                    #print dy
                    ymax = max(y_maxs)*dy
                    ymin = min(y_mins)/dy
            else:
                dy = (max(y_maxs) - min(y_mins))*0.1
                ymax =  max(y_maxs)+dy
                ymin = (min(y_mins)-dy if min(y_mins)!=0 else 0 )
            #print s.axis.GetMinimum(), s.axis.GetMaximum()
        if dim==2 :
            raise NotImplementedError("range not implemented for dim=2");
            pass
        if dim==3 :
            raise NotImplementedError("range not implemented for dim=2");
            pass
        return xmin,xmax,ymin,ymax

    def SetFrameStyle(s,horlist,**kwargs):
        H = horlist
        if not type(H) is list :
            H = [horlist]
        dim = H[0].GetDimension()
        if dim == 1 :
            s.SetFrameStyle1D(H,**kwargs)
        if dim == 2 :
            s.SetFrameStyle2D(H,**kwargs)
        pass

    def SetFrameStyle1D(s,hlist, **kwargs):
        # possible options: scale = 1.0, logY=False, logX=False, maxY="inf", minY="-inf", maxX="inf", minX="-inf", forceRange=False) :
        #parse args
        scale      = kwargs[ "scale"      ] if "scale"      in kwargs else 1.0
        fixRange   = kwargs[ "fixRange"   ] if "fixRange"   in kwargs else True
        forceRange = kwargs[ "forceRange" ] if "forceRange" in kwargs else False
        logX       = kwargs[ "logX"       ] if "logX"       in kwargs else False
        logY       = kwargs[ "logY"       ] if "logY"       in kwargs else False
        # Graph issues
        hlist0 = hlist[0]
        if "TGraph" in hlist0.ClassName() :
            gPad.SetLogx((1 if logX else 0))
            gPad.SetLogy((1 if logY else 0))
            hlist[0].Draw("APL")
            hlist0=hlist0.GetHistogram()
            hlist0.SetTitle(hlist[0].GetTitle())
        # get range
        kwargs["dim"] = 1
        xmin, xmax, ymin, ymax = s.GetRangesFromHists(hlist, **kwargs) #dim=1, logY=logY, logX=logX, maxY=maxY, minY=minY, maxX=maxX, minX=minX )
        if forceRange :
            #print "Force Range", kwargs["maxY"], kwargs["minY"]
            if "maxY" in kwargs : ymax = kwargs["maxY"]
            if "minY" in kwargs : ymin = kwargs["minY"]
            if "maxX" in kwargs : xmax = kwargs["maxX"]
            if "minX" in kwargs : xmin = kwargs["minX"]
            pass
        #MSG.debug(" x min max: {} {} y min max: {} {}".format( xmin, xmax, ymin,ymax) )
        title = ";{};{}".format( hlist0.GetXaxis().GetTitle(), hlist0.GetYaxis().GetTitle())
        xaxis_orig = hlist0.GetXaxis()
        if xaxis_orig.GetBinLabel != "" :
            # seems the histogram has labels, copy them
            s.axes.append(s.EmptyClone(hlist0,"axis"))
            s.axes[-1].SetTitle(title)
        else:
            # not labels, axis will be one-bin histogram
            s.axes.append(TH1F("axis",title ,1,xmin,xmax))
        s.axis = s.axes[-1]
        s.axis.SetDirectory(0)
        s.axis.SetLineColor(0)
        s.axis.SetLineWidth(0)
        s.axis.SetTitle(";{};{}".format(
            hlist0.GetXaxis().GetTitle(),
            hlist0.GetYaxis().GetTitle())
            )
        if fixRange :
            s.axis.SetMaximum(ymax)
            s.axis.SetMinimum(ymin)
            s.axis.GetXaxis().SetRangeUser(xmin,xmax)
            pass
        # global pad style settings
        #s.c1.cd()
        gStyle.SetOptStat(0)
        gPad.SetTicks()
        gPad.SetLeftMargin  ( 0.120 )
        gPad.SetRightMargin ( 0.040 )
        gPad.SetBottomMargin( 0.120*scale  )
        gPad.SetTopMargin( 0.050*scale  )
        #  axis label title offset
        s.axis.SetLabelOffset( 0.012, "X" ) # label offset on x axis
        s.axis.SetLabelOffset( 0.012, "Y" ) # label offset on x axis
        s.axis.GetXaxis().SetTitleOffset( 1.33 )
        s.axis.GetYaxis().SetTitleOffset( 1.25/scale  )
        #  axis label/title size
        labelSize = 0.045*scale
        s.axis.GetXaxis().SetTitleSize( labelSize )
        s.axis.GetYaxis().SetTitleSize( labelSize )
        labelSize = 0.04*scale
        s.axis.GetXaxis().SetLabelSize( labelSize )
        s.axis.GetYaxis().SetLabelSize( labelSize )
        # draw
        gPad.Clear()
        s.axis.Draw()
        gPad.SetLogx((1 if logX else 0))
        gPad.SetLogy((1 if logY else 0))
        pass

    def SetFrameStyle2D(s,hlist,  **kwargs) :
        #parse args 
        scale      = kwargs[ "scale"      ] if "scale"      in kwargs else 1.0
        #fixRange   = kwargs[ "fixRange"   ] if "fixRange"   in kwargs else True
        #forceRange = kwargs[ "forceRange" ] if "forceRange" in kwargs else False
        logX       = kwargs[ "logX"       ] if "logX"       in kwargs else False
        logY       = kwargs[ "logY"       ] if "logY"       in kwargs else False
        logZ       = kwargs[ "logZ"       ] if "logZ"       in kwargs else False
        frameColor = kwargs[ "frameCol"   ] if "frameCol"   in kwargs else 0
        kwargs["dim"] = 2
        # set view port
        xmax = "inf"  if not "maxX" in kwargs else kwargs["maxX"]
        xmin = "-inf" if not "minX" in kwargs else kwargs["minX"]
        ymax = "inf"  if not "maxY" in kwargs else kwargs["maxY"]
        ymin = "-inf" if not "minY" in kwargs else kwargs["minY"]
        # create axis
        s.axes.append(s.EmptyClone(hlist[0],"axis"))
        s.axis = s.axes[-1]
        s.axis.SetTitle("")
        s.axis.GetXaxis().SetRangeUser(float(xmin),float(xmax))
        s.axis.GetYaxis().SetRangeUser(float(ymin),float(ymax))
        # global pad style settings
        # s.c1.cd()
        gStyle.SetOptStat(0)
        gPad.SetFrameFillColor(frameColor)
        gPad.SetTicks()
        gPad.SetLeftMargin   ( 0.120*TMath.Sqrt(scale) )
        gPad.SetRightMargin  ( 0.140                   )
        gPad.SetBottomMargin ( 0.120*scale             )
        gPad.SetTopMargin    ( 0.050*scale             )
        #  axis label title offset
        s.axis.SetLabelOffset( 0.012, "X" ) # label offset on x axis
        s.axis.SetLabelOffset( 0.012, "Y" ) # label offset on x axis
        s.axis.GetXaxis().SetTitleOffset( 1.33 )
        s.axis.GetYaxis().SetTitleOffset( 1.25 )
        #  axis label/title size
        labelSize = 0.045*scale
        s.axis.GetXaxis().SetTitleSize( labelSize )
        s.axis.GetYaxis().SetTitleSize( labelSize )
        labelSize = 0.04*scale
        s.axis.GetXaxis().SetLabelSize( labelSize )
        s.axis.GetYaxis().SetLabelSize( labelSize )
        # draw
        s.axis.Draw()
        gPad.SetLogx((1 if logX else 0))
        gPad.SetLogy((1 if logY else 0))
        gPad.SetLogz((1 if logZ else 0))

    def DrawContour(s, histlist, contours=0) :
        # define contours
        if contours == 0 : 
            contours = [ # sigma levels
                3.17310522791349303e-01,
                4.55002597802489639e-02,
                2.69979614651120765e-03,
                6.33424887837306439e-05,
                5.73303137155808246e-07
                ]
        if isinstance(contours,int) :
            maxi = histlist[0].GetMinimum()
            mini = histlist[0].GetMaximum()
            step = float(maxi-mini)/contours
            contours = range(mini, maxi , step)
        # test contours
        N_cont=0;
        N_cont=len(contours)
        if N_cont == 0: 
            MSG.warning( "Bad contours definition: "+str(contours) )
        # do pseudo overlay
        for drawTechnique in ["samecont0","samecont3"]:
            for hist in histlist :
                # for each contour draw
                for i_cont in range(N_cont) :
                    col_line = s.MakeTColorBrighter(hist.GetLineColor(), i_cont, N_cont)
                    col_fill = s.MakeTColorBrighter(hist.GetFillColor(), i_cont, N_cont)
                    s.contourHists.append (hist.Clone());
                    s.contourHists[-1].SetFillColor(col_fill)
                    s.contourHists[-1].SetLineColor(col_line)
                    s.contourHists[-1].SetContour(1)
                    s.contourHists[-1].SetContourLevel(0, contours[i_cont])
                    s.contourHists[-1].Draw(drawTechnique)

    def DrawLegend(s, objList, drawOptions, **kwargs) :
        scale      = kwargs[ "scale"      ] if "scale"      in kwargs else 1.0
        x          = kwargs[ "legx"          ] if "legx"          in kwargs else 0.7
        y          = kwargs[ "legy"          ] if "legy"          in kwargs else 0.88
        cols       = kwargs[ "legCol"          ] if "legCol"          in kwargs else 1
        ytop =  y
        dy = 0.093*scale * (len(objList))
        dx = 0.48*scale*cols
        s.legend = TLegend(x, ytop-dy, x+dx, ytop)
        s.legend.SetNColumns(cols)
        s.legend.SetLineWidth(0)
        s.legend.SetLineStyle(0)
        s.legend.SetLineColor(0)
        s.legend.SetFillStyle(0)
        #s.DrawBox(x,ytop-dy,dx,dy,fcol=s.c_white_a10)
        #s.legend.SetFillColor(s.c_white_a10)

        if type(drawOptions) != list :
            drawOptions = [ drawOptions for x in objList ]
        for obj in objList :
            i_obj = objList.index(obj)
            title = obj.GetTitle().split(";")[0]
            s.legend.AddEntry(obj, title,drawOptions[i_obj])
        s.legend.Draw()

    # def TColorGetRGB(s,color):
    #     c = gROOT.GetColor(color)
    #     return [c.GetRed(),c.GetGreen(),c.GetBlue()]

    # def MakeBrighterRGB(s, inColorRGB, q=0.3) :
    #     #HEX = 255 # if your insput is 0..1=>1, if 0..256=>1
    #     ## inRGB = map(ord, inColorRGB.decode('hex') ) # makee 3-tuple of numbers 0..256
    #     #inRGB = inColorRGB # 
    #     #print inRGB, q
    #     ## keep hue H = TMath.Sqrt(3)*(R-B)/(2*(R-G-B))
    #     ## luminance Y := 0.2126 R + 0.7152 G + 0.0722 B
    #     #aRGB = [0.2126, 0.7152, 0.0722] # 3-tuple of coeficients
    #     ##maxY = 252 # our maximal luminosity, white is 256
    #     #maxY = 0.8
    #     #minY = sum(a*c for a,c in zip (aRGB, inRGB) ) # current luminosity
    #     #stepY = maxY-minY #lumi step
    #     #outCol = [(a*stepY*q+c) for a,c in zip(aRGB, inRGB)]
    #     #print outCol
    #     #outCol = [int(a*HEX) for a in outCol]
    #     #print outCol
    #     #return "".join(map(chr, outCol)).encode('hex') # code back to hex
    #     pass

    def MakeTColorGrayish(s, inColor):
        col = gROOT.GetColor(inColor)
        H = col.GetHue       ()/360
        S = col.GetSaturation()
        L = col.GetLight     ()
        R = col.GetRed       ()
        G = col.GetGreen     ()
        B = col.GetBlue      ()
        S = S/3
        R,G,B = colorsys.hls_to_rgb(H,L,S)
        return TColor.GetColor(R,G,B)
        pass

    def MakeTColorBrighter(s, inColor, index, N_col):
        col = gROOT.GetColor(inColor)
        H = col.GetHue       ()/360
        S = col.GetSaturation()
        L = col.GetLight     ()
        R = col.GetRed       ()
        G = col.GetGreen     ()
        B = col.GetBlue      ()
        #print inColor, H ,S ,L ,R ,G ,B 
        Lmax = 1
        Lmin = L
        Lstep = (Lmax-Lmin)/N_col 
        L = L + index*Lstep
        R,G,B = colorsys.hls_to_rgb(H,L,S)
        #print inColor, H ,S ,L ,R ,G ,B 
        #print

        return TColor.GetColor(R,G,B)



    def AutoCompareColor(s, i_col, N_col) :
        step = 1./(N_col)
        H = step*i_col
        S = 0.8
        L = 0.4
        R,G,B = colorsys.hls_to_rgb(H,L,S)
        ocol =  TColor.GetColor(R,G,B)
        #MSG.debug( "RGB: {},{},{} ROOT index {}".format( R,G,B,ocol) )
        return ocol

    def AutoCompareColorLight(s, i_col, N_col) :
        return s.MakeTColorBrighter(s.AutoCompareColor(i_col, N_col), 1, 3)

    def AutoCompareFromList(s, i, l):
        N_list = len(l)
        return l[i%N_list]

    def AutoCompareLine(s, i):
        line_list = range(11)
        return s.AutoCompareFromList(i,line_list)
    
    def AutoCompareMarker(s,i):
        marker_list = [
                20,21,22,23,29,33,34, # filled
                24,25,26,32,30,27,28, # empty
                1,2,3,5,6,7           # crosses and dots
                ]
        return s.AutoCompareFromList(i,marker_list)

    def AutoCompareTest(s, N_test=20) :
        # there are 20 different markers
        hlist=list()
        for i in range(N_test) :
            hlist.append(TH1F(str(i), str(i), 1, 0, 1 ))

            hlist[-1].SetLineColor(s.AutoCompareColor(i,N_test))
            hlist[-1].SetLineStyle(s.AutoCompareLine(i))
            hlist[-1].SetLineWidth(2)

            hlist[-1].SetFillColor(s.AutoCompareColorLight(i, N_test))
            hlist[-1].SetFillStyle(1001)

            hlist[-1].SetMarkerColor(s.AutoCompareColor(i,N_test))
            hlist[-1].SetMarkerStyle(s.AutoCompareMarker(i))
            hlist[-1].SetMarkerSize (3)

        s.ShowStyle("AutoCompare"+str(N_test), hlist)

    def AutoSetStyle(s, hists , opt="clmf") :
        N = len(hists)
        for h in hists :
            i = hists.index(h)
            if "c" in opt: h.SetLineColor(s.AutoCompareColor(i,N))
            #if "l" in opt: h.SetLineStyle(s.AutoCompareLine(i))
            #h.SetLineWidth(2)

            #if "c" in opt: h.SetFillColor(s.AutoCompareColorLight(i, N))
            #if "f" in opt: h.SetFillStyle(1001)

            #if "c" in opt: h.SetMarkerColor(s.AutoCompareColor(i,N))
            #if "m" in opt: h.SetMarkerStyle(s.AutoCompareMarker(i))
            h.SetMarkerSize (1)

    def CreateMovingAverageHist(s,H,nbins=2):
        # -- try smooth
        H_MA = H.Clone("average")
        H_MA.Smooth()
        return H_MA
        #
        # take nbins before and nbins after and calculate a mean
        H_MA = s.EmptyClone(H,"average")
        binlist = range(nbins+1,H.GetNbinsX()-nbins+1)
        for ibin in binlist :
            unc = H.GetBinError(ibin)
            # get bin for averaging
            valrange=range (ibin-nbins, ibin+nbins+1) 
            vals = [ H.GetBinContent(ii) for ii in valrange ]
            # make an average
            avg = sum(vals)/len(vals) 
            # set new bin val
            H_MA.SetBinContent(ibin, avg )
            H_MA.SetBinError  (ibin, unc )
            print ibin, valrange, vals, avg
        return H_MA


    def CreateChiHists(s, A, B):
        C = s.CreateHistsFromFun(A, B, s.dif_chi)
        return C

    def CreateRatioHists(s, A, B):
        C = s.CreateHistsFromFun(A, B, s.divide_bins)
        return C

    def CreateRatio0Hists(s, A, B):
        C = s.CreateHistsFromFun(A, B, s.divide_bins0)
        return C

    def CreateSubtrctHists(s, A, B):
        C = s.CreateHistsFromFun(A, B, s.norm_subtract_bins)
        #C[0].GetYaxis().SetTitle( "1-"+A.GetTitle()+"/other")
        return C

    def CreateHistsFromFun(s, A, B, function) :
        # create C=A/B
        # B could be list
        C = list()
        Blist = B
        if not isinstance(B, list): # if is not list
            Blist = [B]
        # create 
        if "TGraph" in A.ClassName() : # graphs
            for b in Blist:
                c = b.Clone(b.GetName()+"_fun")
                for i in range(b.GetN()) :
                    a_val = A.GetY()[i]
                    a_err = 0
                    b_val = b.GetY()[i]
                    b_err = 0
                    c_val = 0
                    c_err2 = 0
                    if "TGraphErrors" == A.ClassName() :
                        a_err = A.GetEY()[i]
                        b_err = b.GetEY()[i]
                    elif "TGraphAsymmErrors" == A.ClassName() :
                        a_err = A.GetEYlow()[i]
                        b_err = b.GetEYlow()[i]
                        [c_val, c_err2] = function(a_val, a_err, b_val, b_err)
                        a_err = A.GetEYhigh()[i]
                        b_err = b.GetEYhigh()[i]
                    [c_val, c_err] = function(a_val, a_err, b_val, b_err)
                    c.SetPoint(i,A.GetX()[i], c_val)
                    if "TGraphErrors" == A.ClassName() :
                        c.SetPointError(i,0, c_err)
                    elif "TGraphAsymmErrors" == A.ClassName() :
                        c.SetPointError(i,0,0, c_err2, c_err)
                    pass
                C.append(c)
            pass
        else : # histograms
            dim=A.GetDimension()
            for b in Blist:
                c = s.EmptyClone(b, b.GetName()+"_fun")
                for xbin in range(1, c.GetNbinsX()+1):
                    ybinlist = [0] if dim<2 else range(1,A.GetNbinsY()+1)
                    for ybin in ybinlist:
                        zbinlist = [0] if dim<3 else range(1,A.GetNbinsZ()+1)
                        for zbin in zbinlist:
                            a_val = A.GetBinContent (xbin,ybin,zbin)
                            a_err = A.GetBinError   (xbin,ybin,zbin)
                            b_val = b.GetBinContent (xbin,ybin,zbin)
                            b_err = b.GetBinError   (xbin,ybin,zbin)
                            [c_val, c_err] = function(a_val, a_err, b_val, b_err)
                            c.SetBinContent(xbin,ybin,zbin, c_val);
                            c.SetBinError  (xbin,ybin,zbin, c_err);
                            #print "a {} b {} c {}".format(a_val,b_val,c_val)
                #print "C"
                C.append(c)
                pass
            pass
        if isinstance(B, list) :
            return C
        else :
            return C[0]

    def DrawHistCompareSubPlot(s, mainHists, subHists, **kwargs): 
        #supported kwargs: logy=False, logx=False , ymax="inf", ymin="-inf", xmax="inf", xmin="-inf", forcerange=False) :
        compareType = kwargs[ "compareType" ] if "compareType" in kwargs else "subtract"
        cdiv        = kwargs[ "cdiv"        ] if "cdiv"        in kwargs else 0.55
        drawOpt     = kwargs[ "drawOpt"     ] if "drawOpt"     in kwargs else "same"

        s.c1.Divide(1,2)
        # draw main
        s.c1.cd(1)
        kwargs["scale"] = 1./cdiv
        s.SetFrameStyle1D(mainHists, **kwargs) # scale = 1./cdiv, logY=logy, logX=logx, maxY=ymax, minY=ymin, maxX=xmax, minX=xmin, forceRange=forcerange)
        s.DrawHistCompare(mainHists,drawOpt)
        # draw subplot
        s.c1.cd(2)
        ratioArgs = {
                "scale" : 1./(1-cdiv) ,
                }
        if "logX" in kwargs : ratioArgs["logX"] = kwargs["logX"]
        if "minX" in kwargs : ratioArgs["minX"] = kwargs["minX"]
        if "maxX" in kwargs : ratioArgs["maxX"] = kwargs["maxX"]
        if compareType == "ratio" :
            #ratioArgs["logY"]=False
            ratioArgs["minY"]= -25
            ratioArgs["maxY"]= 25
            s.SetFrameStyle1D(subHists, **ratioArgs) # scale = 1./(1-cdiv), logY=True, minY=0.04, maxY=25, logX=logx)
        elif compareType=="ratio0" :
            #ratioArgs["logY"]=False
            ratioArgs["minY"]= 0.955
            ratioArgs["maxY"]= 1.045
            ratioArgs["forceRange"]= True
            s.SetFrameStyle1D(subHists, **ratioArgs) # scale = 1./(1-cdiv), logY=True, minY=0.04, maxY=25, logX=logx)
            s.axis.GetYaxis().SetNdivisions(5,2,0,True)
        elif compareType=="subtract" :
            ratioArgs["minY"]=-10.
            ratioArgs["maxY"]= 10.
            s.SetFrameStyle1D(subHists, **ratioArgs) # scale = 1./(1-cdiv), minY=-10, maxY=10, logX=logx)
        elif compareType=="chi" :
            ratioArgs["minY"]= -3.0
            ratioArgs["maxY"]= 3.0
            #ratioArgs["fixRange"]=False
            ratioArgs["forceRange"]=True
            s.SetFrameStyle1D(subHists, **ratioArgs) # scale = 1./(1-cdiv), minY=-10, maxY=10, logX=logx)
        else :
            raise NotImplementedError("Uknown compare type: "+compareType)
        s.DrawHistCompare(subHists,drawOpt)
        # edit the position
        s.AdjustSubplot(cdiv)

    def AdjustSubplot(s, cdiv=0.8 ):
        # cdiv is canvas division
        #mainplot
        s.c1.cd(1)
        gPad.SetPad(0., 1.0-cdiv, 1.0, 1.0)
        gPad.SetBottomMargin(0)
        #subplot
        s.c1.cd(2)
        gPad.SetPad(0., 0.0, 1.0, 1.0-cdiv)
        gPad.SetTopMargin(0)
        s.c1.cd()

    def DrawMeasurement(s, x="inf", xerr="-inf", xerrP="-inf", xerrN="-inf", 
                           y="inf", yerr="-inf", yerrP="-inf", yerrN="-inf", 
                        type="band", 
                        cfill=TColor.GetColor("#AAAAAA"), cline=TColor.GetColor("#000000"), spoint=20 ):
        wline=1.
        wmarker=1.5
        wpoint=2.0
        rel_bar_size=0.02
        if xerr != "-inf": 
            xerrP=xerr
            xerrN=xerr
        if yerr != "-inf": 
            yerrP=yerr
            yerrN=yerr
        if type == "band" :
            # get ranges for bands and lines
            xmin = s.axis.GetXaxis().GetXmin()
            xmax = s.axis.GetXaxis().GetXmax()
            ymin = s.axis.GetYaxis().GetXmin()
            ymax = s.axis.GetYaxis().GetXmax()
            # plot y measurement
            if y != "inf" :
                if yerrP != "-inf" and yerrN != "-inf" :
                    s.DrawMeasureBand(xmin,y-abs(yerrN),xmax,y+abs(yerrP), cfill, 1001 )
                    # this?
                    #s.DrawMeasureLine(xmin,y-abs(yerrN),xmax,y+abs(yerrP), cline, 1, 1 )
                    #s.DrawMeasureLine(xmin,y-abs(yerrN),xmax,y+abs(yerrP), cline, 1, 1 )
                    # or this?
                    # s.DrawMeasureLine(xmin,y-abs(yerrN),xmax,y-abs(yerrN), cline, 1, 1 )
                    # s.DrawMeasureLine(xmin,y+abs(yerrP),xmax,y+abs(yerrP), cline, 1, 1 )
                s.DrawMeasureLine(xmin,y,xmax,y, cline, 1, 2)
            # plot x measurement
            if x != "inf" :
                if xerrP != "-inf" and xerrN != "-inf" :
                    s.DrawMeasureBand(x-abs(xerrN),ymin,x+abs(xerrP),ymax, cfill, 1001 )
                    s.DrawMeasureLine(x-abs(xerrN),ymin,x+abs(xerrP),ymax, cline, 1, 1 )
                    s.DrawMeasureLine(x-abs(xerrN),ymin,x+abs(xerrP),ymax, cline, 1, 1 )
                s.DrawMeasureLine(x,ymin,x,ymax, cline, 1, 2)
        if type == "point" :
            DrawMeasurePoint(x,xerrP,xerrN,y,yerrP,yerrN,col,20)
            pass

    def DrawMeasureLine(s, x1,y1,x2,y2,col,width,style):
        s.measurements.append(TGraph(2))
        s.measurements[-1].SetPoint(0, x1, y1)
        s.measurements[-1].SetPoint(1, x2, y2)
        s.measurements[-1].SetLineColor(col)
        s.measurements[-1].SetLineStyle(style)
        s.measurements[-1].SetLineWidth(width)
        s.measurements[-1].Draw("L")

    def DrawMeasureBand(s, x1,y1,x2,y2,col,style) :
        s.measurements.append(TGraph(4))
        s.measurements[-1].SetPoint(0, x1, y1)
        s.measurements[-1].SetPoint(1, x1, y2)
        s.measurements[-1].SetPoint(2, x2, y2)
        s.measurements[-1].SetPoint(3, x2, y1)
        s.measurements[-1].SetFillColor(col)
        #s.measurements[-1].SetLineColor(ROOT.kRed)
        s.measurements[-1].SetFillStyle(style)
        s.measurements[-1].Draw("F")

    def DrawMeasurePoint(s, x,xerrP,xerrN,y,yerrP,yerrN,cpoint,spoint,wpoint) :
        pass

    def DrawHistStack(s, listOfHistos ) :
        pass

    def DrawHistCompare(s, listOfHistos,opt="" ) :
        for h in listOfHistos :
            dim = 1
            plotopt=str(h.GetDrawOption()).replace("same","")+","+opt.replace("same","")
            try :
                dim = h.GetDimension()
            except AttributeError:
                pass
            if dim ==2 : 
                plotopt+=",same,COLZ"
            else :
                plotopt+=",same"
            h.Draw(plotopt)
        gPad.Update();

    def WriteText(s, text, x=0.2,y=0.8, **kwargs):
        useNDC =  kwargs ["useNDC"] if "useNDC" in kwargs else True
        col   = kwargs ["tcol"   ] if "tcol"   in kwargs else 1
        align = kwargs ["talign" ] if "talign" in kwargs else 11
        size  = kwargs ["tsize"  ] if "tsize"   in kwargs else 0.025
        s.latex.SetTextColor(col)
        s.latex.SetTextSize(size)
        s.latex.SetTextAlign(align)
        if useNDC: s.latex.SetNDC();
        else : s.latex.SetNDC(False)
        s.latex.DrawLatex(x,y,text)
        pass

    def WriteTextWithBox(s, text, x=0.2,y=0.8, **kwargs) : 
        tsize  = kwargs ["tsize"  ] if "tsize"   in kwargs else 0.025
        #tcol   = kwargs ["tcolor" ] if "tcolor"  in kwargs else 1
        w = tsize*1.6
        h = tsize*1.2
        x_toff = tsize*2
        y_toff = tsize*0.1
        s.DrawBox(x,y,w,h, **kwargs)
        s.WriteText(text, x+x_toff, y+y_toff, **kwargs)
        pass

    def DrawMarker(s,x,y,**kwargs):
        useNDC =  kwargs ["useNDC"] if "useNDC" in kwargs else True
        c_mark = kwargs ["mcolor"] if "mcolor" in kwargs else 1
        s_mark = kwargs ["mstyle"] if "mstyle" in kwargs else 20
        w_mark = kwargs ["msize"]  if "msize"  in kwargs else 2.
        # uncert
        y_err = kwargs ["yerr"] if "yerr" in kwargs else 0
        x_err = kwargs ["xerr"] if "xerr" in kwargs else 0
        if not "lcolor"  in kwargs : kwargs["lcolor"] = c_mark
        if not "lwidth"  in kwargs : kwargs["lwidth"] = 2
        # empty markers
        empty_markers={ 24:20, 25:21, 26:22, 27:33, 28:34, 30:29, 32:23 }
        is_empty_marker = s_mark in empty_markers
        thickness=0.25
        if is_empty_marker :
            s_mark = empty_markers[s_mark]
        # set style
        s.mark = TMarker()
        s.mark.SetMarkerColor (c_mark)
        s.mark.SetMarkerStyle (s_mark)
        s.mark.SetMarkerSize  (w_mark)
        if y_err!=0 :
            s.DrawErrorBar(x,y-yerr,x,y+yerr, **kwargs)
        if x_err!=0 :
            s.DrawErrorBar(x-x_err,y,x+x_err,y, **kwargs)
            pass
        if useNDC :
            x = s.NDC2val_x(x)
            y = s.NDC2val_y(y)
        if is_empty_marker :
            s.mark.DrawMarker(x,y)
            s.mark.SetMarkerColor(0) # can be change to GetFrameFillColor
            s.mark.SetMarkerSize(w_mark*(1-thickness))
        s.mark.DrawMarker(x,y)
        pass

    def DrawErrorBar(s,x1,y1,x2,y2, **kwargs):
        useNDC =  kwargs ["useNDC"] if "useNDC" in kwargs else True
        errBar_style = kwargs ["estyle"]  if "estyle"  in kwargs else "-"
        if "-" in errBar_style :
            s.DrawLine(x1,y1,x2,y2, **kwargs)
            pass
        if "=" in errBar_style :
            raise NotImplementedError("Please make the bar error.")
        if "|" in errBar_style :
            # maker size 1 corresponds to 8pix
            pix8_NDC=8./gPad.GetWw()
            dd = kwargs["msize"] * 0.5 * pix8_NDC
            if x1==x2 : # y error -- bars arw in x
                if not useNDC : dd = s.NDC2val_x(dd,True)
                s.DrawLine(x1-dd,y1,x1+dd,y1, **kwargs)
                s.DrawLine(x1-dd,y2,x1+dd,y2, **kwargs)
                pass
            if y1==y2 : # x error -- bars are in y
                if not useNDC : dd = s.NDC2val_y(dd,True)
                s.DrawLine(x1,y1-dd,x1,y1+dd, **kwargs)
                s.DrawLine(x2,y1-dd,x2,y1+dd, **kwargs)
                pass
            pass
        pass

    def DrawLine(s,x1,y1,x2,y2, **kwargs) :
        useNDC =  kwargs ["useNDC"] if "useNDC" in kwargs else True
        c_line = kwargs ["lcolor"]  if "lcolor"  in kwargs else 1
        s_line = kwargs ["lstyle"]  if "lstyle"  in kwargs else 1
        w_line = kwargs ["lwidth"]  if "lwidth"  in kwargs else 2
        s.line = TLine()
        s.line.SetLineColor (c_line)
        s.line.SetLineStyle (s_line)
        s.line.SetLineWidth (w_line)
        if useNDC :
            s.line.DrawLineNDC(x1,y1,x2,y2)
        else:
            s.line.DrawLine(x1,y1,x2,y2)
        pass



    def DrawBox(s, x, y, width, height, **kwargs) :
        c_fill = kwargs ["fcol"] if "fcol" in kwargs else 1
        #c_trans = kwargs ["fcol"] if "fcol" in kwargs else 1
        c_line = kwargs ["lcol"] if "lcol" in kwargs else -1
        useNDC =  kwargs ["useNDC"] if "useNDC" in kwargs else True
        #print "box ",x,y,width,height,c_fill,c_line
        # coordinates
        x_min = x
        y_min = y
        x_max = x+width
        y_max = y+height
        if useNDC : 
            x_min = s.NDC2val_x( x_min ) 
            y_min = s.NDC2val_y( y_min ) 
            x_max = s.NDC2val_x( x_max ) 
            y_max = s.NDC2val_y( y_max ) 
            #print "NDC " , x_min, x_max, y_min, y_max
        # box printing
        s.boxes.append(TBox( x_min , y_min , x_max , y_max))
        s.boxes[-1].SetFillColor(c_fill)
        #s.boxes[-1].    SetFillColorAlpha (c_fill,0.5)
        #s.boxes[-1].SetFillTran(c_fill)
        s.boxes[-1].SetFillStyle(1001)
        if c_line > 0 :
            s.boxes[-1].SetLineColor(c_line);
            s.boxes[-1].SetLineStyle(ROOT.kSolid);
            s.boxes[-1].SetLineWidth(2);
            s.boxes[-1].Draw("fl");
        else :
            s.boxes[-1].Draw("f");
            pass
        pass


    def CopyStyle(s, fromH, toH) :
        styles=[
                "LineColor",
                "LineWidth",
                "LineStyle",
                "FillColor",
                "FillStyle",
                "MarkerColor",
                "MarkerSize",
                "MarkerStyle",
                ]
        for st in styles :
            exec("toH.Set"+st+"(fromH.Get"+st+"())")

    def ShowStyle(s, name, listOfHistos) :
        s.NewCanvas(name)
        N_hist=len(listOfHistos)
        hBar = list()
        hDot = list()
        names=[h.GetTitle() for h in listOfHistos]
        for h in listOfHistos :
            i_h = listOfHistos.index(h)
            hBar.append(TH1I("bar_"+h.GetName(), "", N_hist, 0, N_hist))
            s.CopyStyle(h,hBar[-1])
            hBar[-1].Fill(i_h,1)

            hDot.append(TH1I("dot_"+h.GetName(), "", N_hist, 0, N_hist))
            s.CopyStyle(h,hDot[-1])
            hDot[-1].Fill(i_h,2)
        s.SetFrameStyle1D(hDot)
        [s.axis.GetXaxis().SetBinLabel(i+1, listOfHistos[i].GetTitle()) for i in range(N_hist)]
        for h in hBar:
            h.Draw("SAMEB")
        for h in hDot:
            h.Draw("SAMEPE")
        s.Save()


    def ColorHTML(s,col):
        coltext=str(col).replace("#","")
        if len(coltext) == 8 :
            nonA = TColor.GetColor("#"+coltext[0:-2])
            alpha = int(coltext[-2:],16)/256. # hex to int and divide by max
            print "chtml", nonA, alpha
            return s.make_color_transparent(nonA,alpha)
        return TColor.GetColor("#"+coltext)

    def MakeNiceGradient(s):
        import numpy as np
        # red -- light yellow -- blue
        Red    = [ 165 , 215 , 244 , 253 , 254 , 224 , 171 , 116 , 69  , 49  ][::-1]
        Green  = [ 0   , 48  , 109 , 174 , 224 , 243 , 217 , 173 , 117 , 54  ][::-1]
        Blue   = [ 38  , 39  , 67  , 97  , 144 , 248 , 233 , 209 , 180 , 149 ][::-1]
        Length = [ 0   , 1   , 2   , 3   , 4   , 5   , 6   , 7   , 8   , 9   ] 
        ln = len(Length)
        TColor.CreateGradientColorTable(
            len(Length),
            np.array(Length )/(ln-1.) ,
            np.array(Red    )/254.    ,
            np.array(Green  )/254.    ,
            np.array(Blue   )/254.    ,
            50
            )
        pass

    def make_color_transparent(s,old_icol,val):
       new_icol=1000+old_icol
       old_col = gROOT.GetColor(old_icol)
       new_name=old_col.GetName()+"_"+str(val)
       if gROOT.GetColor(new_icol) != None : return new_icol;
       s.new_col .append(  TColor(new_icol, old_col.GetRed(),old_col.GetGreen(),old_col.GetBlue(),new_name,val) )
       #print "new transparent color ", new_icol
       print "transparent color ", new_icol, gROOT.GetColor(new_icol)
       print "new color ", s.new_col[-1].GetNumber(), gROOT.GetColor(new_icol).GetName()
       return s.new_col[-1].GetNumber();

    def norm_subtract_bins(s, a_val, a_err, b_val, b_err) :
        # C = 1-A/B
        if b_val==0 : return [0,0]
        c_val = 1-float(a_val)/b_val
        c_err = TMath.Sqrt(
                TMath.Power(float(a_err)/b_val, 2) +
                TMath.Power(float(b_err)*a_val/(b_val*b_val), 2)
                )
        return [c_val, c_err]

    #def norm_subtract_bins(s, a_val, a_err, b_val, b_err) :
        ## C = A/B
        #if b_val==0 : return [0,0]
        #c_val = 1-float(a_val)/b_val
        #c_err = TMath.Sqrt(
                #TMath.Power(float(a_err)/b_val, 2) +
                #TMath.Power(float(b_err)*a_val/(b_val*b_val), 2)
                #)
        #return [c_val, c_err]

    def dif_chi(s, a_val,a_err,b_val,b_err):
        sigma = TMath.Sqrt(a_err**2 + b_err**2)
        if sigma==0 : sigma = TMath.Sqrt(a_val+b_val)
        if sigma==0 : return [0,0]
        c_val = (a_val - b_val) / sigma
        c_err = 1.
        return [c_val, c_err]

    def divide_bins(s, a_val, a_err, b_val, b_err) :
        # C = A/B
        if b_val==0 : return [0,0]
        c_val = float(a_val)/b_val
        #c_err = TMath.Sqrt(
                #TMath.Power(float(a_err)/b_val, 2) +
                #TMath.Power(float(b_err)*a_val/(b_val*b_val), 2)
                #)
        a_rel=a_err/a_val if a_val!=0 else 0.
        b_rel=b_err/b_val
        c_err = c_val * TMath.Sqrt(a_rel*a_rel + b_rel*b_rel)
        return [c_val, c_err]

    def divide_bins0(s, a_val, a_err, b_val, b_err) :
        # C = B/A
        return s.divide_bins(b_val,b_err,a_val,a_err)


    def is_a_number(s, inpt) :
        try:
            float(inpt)
            return True
        except ValueError:
            return False

    # helpers for coordinates on canvas
    def NDC2val_x (s, x, rel = False) : 
        gPad.Update()
        x_min = gPad.GetX1() #  frame->GetXaxis()->GetXmax();
        x_max = gPad.GetX2() #  frame->GetXaxis()->GetXmin();
        if rel : return x*(x_max-x_min)
        return x_min+x*(x_max-x_min)

    def NDC2val_y (s, y, rel = False) :
        gPad.Update()
        y_min = gPad.GetY1() # frame->GetYaxis()->GetXmax();
        y_max = gPad.GetY2() # frame->GetYaxis()->GetXmin();
        if rel : return y*(y_max-y_min)
        return y_min+y*(y_max-y_min)

    def val_y2NDC (s, y, rel = False) :
        gPad.Update()
        y_min = gPad.GetY1() # frame->GetYaxis()->GetXmax();
        y_max = gPad.GetY2() # frame->GetYaxis()->GetXmin();
        width = (y_max-y_min)
        if rel : return y/width
        return (y-y_min)/width

    def CompareHistsInFilesOld(s, name, hflist, normalise=False, logy=False, logx=False , ymax="inf", ymin="-inf", xmax="inf", xmin="-inf", forcerange=False) :
        # you can use list CompareHistsInFiles
        pass

    def CompareHistsInFiles(s, name, hflist, **kwargs) :
        # hflist = title, filename, histname
        hists = [ s.GetHistSetTitle(x[0], x[1], x[2]) for x in hflist ]
        s.CompareHistsInList(name, hists, **kwargs)

    def CompareHistsInList(s, name, hists, **kwargs) :
        # supported kwargs:  normalise=False, 
        # old: logy=False, logx=False , ymax="inf", ymin="-inf", xmax="inf", xmin="-inf", forcerange=False

        # parse options
        normalise= kwargs["normalise"] if "normalise" in kwargs else False
        setXlabel= kwargs["setXlabel"] if "setXlabel" in kwargs else False
        compareType = kwargs[ "compareType" ] if "compareType" in kwargs else "subtract"
        doSave = kwargs[ "doSave" ] if "doSave" in kwargs else True
        doStyle = kwargs[ "doStyle" ] if "doStyle" in kwargs else True

        s.NewCanvas(name)
        if (doStyle) : s.AutoSetStyle(hists)
        # for compare no fill, unless we develop transparent colors
        [ h.SetFillStyle(0) for h in hists ]

        if setXlabel :
            [ h.GetXaxis().SetTitle(setXlabel) for h in hists]

        if normalise :
            [ h.Scale ( 1./ h.Integral()) for h in hists]

        if len(hists)<1 :
            compareType=None

        try :
            s.hRatio = list()
            kwargsRatio = kwargs
            kwargsRatio = kwargs
            if compareType == "ratio":
                s.hRatio = s.CreateRatioHists(hists[0],hists)
                s.hRatio[0].GetYaxis().SetTitle( hists[0].GetTitle()+"/other")
            elif compareType == "ratio0":
                s.hRatio = s.CreateRatio0Hists(hists[0],hists)
                s.hRatio[0].GetYaxis().SetTitle( "ratio" ) # "other/"+hists[0].GetTitle())
            elif compareType == "subtract":
                s.hRatio = s.CreateSubtrctHists(hists[0],hists)
                s.hRatio[0].GetYaxis().SetTitle( "1-"+hists[0].GetTitle()+"/other")
            elif compareType == "chi":
                s.hRatio = s.CreateChiHists(hists[0],hists)
                s.hRatio[0].GetYaxis().SetTitle( "#chi")
                kwargsRatio["drawOpt"]= "p,same"
            else :
                raise NotImplementedError("Uknown compare type: "+str(compareType))
            s.DrawHistCompareSubPlot(hists, s.hRatio, **kwargs)  # ratiology, logx, ymax, ymin, xmax, xmin, forcerange
            s.c1.cd(1)
            s.DrawLegend(hists,"pl",**kwargs)
        except NotImplementedError :
            s.SetFrameStyle1D(hists,**kwargs)
            s.DrawHistCompare(hists)
            s.DrawLegend(hists,"pl", **kwargs)

        if doSave : s.Save()
        pass


    def MakePreviewFromFolder(s, path) :
        hasPDFTK=False
        if hasPDFTK :
            os.system(" rm -f {1}/preview.pdf; pdftk {0}/*.pdf cat output {1}/preview.pdf".format(path, s.imgDir))
        else :
            # create list (with wildcards)
            import glob
            hlist = glob.glob(path)
            hlist = [ x.replace(".pdf","") for x in hlist]
            # use prev, from list
            s.MakePreviewFromList(hlist)

    def MakePreviewFromList(s, figlist = 0, fname="preview") :
        if figlist == 0 : figlist = s.updated_plots
        hasPDFTK=False
        if hasPDFTK :
            q = "rm -f {0}/{1}.pdf; pdftk "
            for f in figlist :
                q+=f+".pdf "
                pass
            q+= "cat output {0}/{1}.pdf"
            qqdo = q.format(s.imgDir,fname)
            #print qqdo
            os.system(qqdo)
        else :
            # use tex
            # create tex file
            # create header
            textxt = r"""
\documentclass{beamer}
% add page numbers
\addtobeamertemplate{navigation symbols}{}{%
\usebeamerfont{footline}%
\usebeamercolor[fg]{footline}%
\hspace{1em}%
\insertframenumber/\inserttotalframenumber
}
\setbeamercolor{background canvas}{bg=black}
\begin{document}

"""
            # make body
            tmpl=r"""
\begin{frame}{TITLE}
\begin{center}
\includegraphics[width=0.9\textwidth]{FILE}
\end{center}
\end{frame}

"""
            for fig in figlist:
                figname = fig.split("/")[-1]
                imgfile=figname+".pdf"
                title=figname.replace("_", " ")
                textxt+=tmpl.replace("TITLE",title).replace("FILE",imgfile)
            # close file
            textxt+=r"""
\end{document}
            """
            # write file
            texfile = open(s.imgDir+"/"+fname+".tex",'w')
            texfile.write(textxt)
            texfile.close()
            # compile tex file
            os.system("cd {}; pdflatex -interaction=batchmode {}.tex  > /dev/null".format(s.imgDir, fname))
            os.system("cd {}; pdflatex -interaction=batchmode {}.tex  > /dev/null".format(s.imgDir, fname))



## Documentation for main
#
# More details. 
if __name__ == '__main__' :
    print "Just use this class."
    pass



kAzure=860

# 0 - fill
# 4 - line
# 3 - darker fill
# 1 - lighter fill

#####  Color Palette by Paletton.com
######  Palette URL: http://paletton.com/#uid=40c0-0kdP++00++70++jv+Yp2+2
# colorscheme RGB (where blue is kAzure-9)
# *** Red
# 
#    shade 0 = #FFAB91 = rgb(255,171,145) = rgba(255,171,145,1) = rgb0(1,0.671,0.569)
#    shade 1 = #FFFFFF = rgb(255,255,255) = rgba(255,255,255,1) = rgb0(1,1,1)
#    shade 2 = #FFD5C7 = rgb(255,213,199) = rgba(255,213,199,1) = rgb0(1,0.835,0.78)
#    shade 3 = #FF8964 = rgb(255,137,100) = rgba(255,137,100,1) = rgb0(1,0.537,0.392)
#    shade 4 = #FF6737 = rgb(255,103, 55) = rgba(255,103, 55,1) = rgb0(1,0.404,0.216)
# 
# *** Blue
# 
#    shade 0 = #99CDFF = rgb(153,205,255) = rgba(153,205,255,1) = rgb0(0.6,0.804,1)
#    shade 1 = #FFFFFF = rgb(255,255,255) = rgba(255,255,255,1) = rgb0(1,1,1)
#    shade 2 = #CBE6FF = rgb(203,230,255) = rgba(203,230,255,1) = rgb0(0.796,0.902,1)
#    shade 3 = #6FB8FF = rgb(111,184,255) = rgba(111,184,255,1) = rgb0(0.435,0.722,1)
#    shade 4 = #45A2FC = rgb( 69,162,252) = rgba( 69,162,252,1) = rgb0(0.271,0.635,0.988)
# 
# *** Green
# 
#    shade 0 = #D1FF91 = rgb(209,255,145) = rgba(209,255,145,1) = rgb0(0.82,1,0.569)
#    shade 1 = #FFFFFF = rgb(255,255,255) = rgba(255,255,255,1) = rgb0(1,1,1)
#    shade 2 = #E8FFC7 = rgb(232,255,199) = rgba(232,255,199,1) = rgb0(0.91,1,0.78)
#    shade 3 = #BEFF64 = rgb(190,255,100) = rgba(190,255,100,1) = rgb0(0.745,1,0.392)
#    shade 4 = #ABFE37 = rgb(171,254, 55) = rgba(171,254, 55,1) = rgb0(0.671,0.996,0.216)
# 
# *** Torqua
# 
#    shade 0 = #91FFC7 = rgb(145,255,199) = rgba(145,255,199,1) = rgb0(0.569,1,0.78)
#    shade 1 = #FFFFFF = rgb(255,255,255) = rgba(255,255,255,1) = rgb0(1,1,1)
#    shade 2 = #C7FFE3 = rgb(199,255,227) = rgba(199,255,227,1) = rgb0(0.78,1,0.89)
#    shade 3 = #64FFB0 = rgb(100,255,176) = rgba(100,255,176,1) = rgb0(0.392,1,0.69)
#    shade 4 = #37FD98 = rgb( 55,253,152) = rgba( 55,253,152,1) = rgb0(0.216,0.992,0.596)
#
#

######  Color Palette by Paletton.com
######  Palette URL: http://paletton.com/#uid=73-1q0kdP++00++70++jv+Yp2+2
# 
# 
# *** Blue
# 
#    shade 0 = #99CDFF = rgb(153,205,255) = rgba(153,205,255,1) = rgb0(0.6,0.804,1)
#    shade 1 = #FFFFFF = rgb(255,255,255) = rgba(255,255,255,1) = rgb0(1,1,1)
#    shade 2 = #CBE6FF = rgb(203,230,255) = rgba(203,230,255,1) = rgb0(0.796,0.902,1)
#    shade 3 = #6FB8FF = rgb(111,184,255) = rgba(111,184,255,1) = rgb0(0.435,0.722,1)
#    shade 4 = #45A2FC = rgb( 69,162,252) = rgba( 69,162,252,1) = rgb0(0.271,0.635,0.988)
# 
# *** Violete
# 
#    shade 0 = #FF91F1 = rgb(255,145,241) = rgba(255,145,241,1) = rgb0(1,0.569,0.945)
#    shade 1 = #FFFFFF = rgb(255,255,255) = rgba(255,255,255,1) = rgb0(1,1,1)
#    shade 2 = #FFC7F8 = rgb(255,199,248) = rgba(255,199,248,1) = rgb0(1,0.78,0.973)
#    shade 3 = #FF64EC = rgb(255,100,236) = rgba(255,100,236,1) = rgb0(1,0.392,0.925)
#    shade 4 = #FD37E4 = rgb(253, 55,228) = rgba(253, 55,228,1) = rgb0(0.992,0.216,0.894)
# 
# *** Green
# 
#    shade 0 = #E9FF91 = rgb(233,255,145) = rgba(233,255,145,1) = rgb0(0.914,1,0.569)
#    shade 1 = #FFFFFF = rgb(255,255,255) = rgba(255,255,255,1) = rgb0(1,1,1)
#    shade 2 = #F4FFC7 = rgb(244,255,199) = rgba(244,255,199,1) = rgb0(0.957,1,0.78)
#    shade 3 = #DFFF64 = rgb(223,255,100) = rgba(223,255,100,1) = rgb0(0.875,1,0.392)
#    shade 4 = #D6FF37 = rgb(214,255, 55) = rgba(214,255, 55,1) = rgb0(0.839,1,0.216)
# 
# *** Orange
# 
#    shade 0 = #FFD291 = rgb(255,210,145) = rgba(255,210,145,1) = rgb0(1,0.824,0.569)
#    shade 1 = #FFFFFF = rgb(255,255,255) = rgba(255,255,255,1) = rgb0(1,1,1)
#    shade 2 = #FFE8C7 = rgb(255,232,199) = rgba(255,232,199,1) = rgb0(1,0.91,0.78)
#    shade 3 = #FFC064 = rgb(255,192,100) = rgba(255,192,100,1) = rgb0(1,0.753,0.392)
#    shade 4 = #FFAE37 = rgb(255,174, 55) = rgba(255,174, 55,1) = rgb0(1,0.682,0.216)
# 
# #####  Generated by Paletton.com (c) 2002-2014

#### Special kAzure-9: almost grey
#
#   shade 0 = #BEC4C9
#   shade 1 = #FFFFFF
#   shade 2 = #DEE1E4
#   shade 3 = #A4ABB2
#   shade 4 = #88919A
#





### Color Brewer 2.0

