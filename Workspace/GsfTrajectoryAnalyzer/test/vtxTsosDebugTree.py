import sys
from math import sqrt,exp,log,pi
from PyMultiGaussianState1D import *
import ROOT

varAxisTitles = [ "q/p [GeV^{-1}]", "local dx/dz", "local dy/dz", "local dx [cm]", "local dy [cm]" ]

gobjects = [ ]

def drawMixture(pad,title,mgs,xmin=None,xmax=None,superimpose=False):
    g = mgs.graph()
    color = 4 if not superimpose else 2
    g.SetLineColor(color)
    g.SetLineWidth(2)
    gobjects.append(g)

    x1,x2 = mgs.xminmax()
    mode = mgs.mode()
    y2 = mgs.mode()[1]
    if xmin==None:
        xmin = x1-0.05*(x2-x1)
    if xmax==None:
        xmax = x2+0.05*(x2-x1)
    if not superimpose:
        hframe = pad.DrawFrame(xmin,y2/5000.,xmax,2.*y2)
        hframe.GetXaxis().SetTitle(title)
    g.Draw("C")

    graphs = mgs.singleGraphs()
    gobjects.extend(graphs)
    for g in graphs:
        g.SetLineColor(color)
        g.SetLineStyle(2)
        g.Draw("C")
    pad.SetLogy(1)
    pad.SetGridx(1)
    pad.SetGridy(1)
    pad.Update()

    line = "Mean : {0:10.2e}#pm{1:9.2e}  Mode : {2:10.2e}".format(mgs.mean(),mgs.sigma(),mode[0])
    if superimpose:
        latex = ROOT.TLatex(0.1,0.91,line)
        latex.SetTextColor(2)
    else:
        latex = ROOT.TLatex(0.1,0.95,line)
        latex.SetTextColor(4)
    latex.SetNDC(1)
#        pave = ROOT.TPaveText(0.1,0.80,0.5,0.9,"NDC NB")
    latex.SetTextSize(0.04)
    latex.SetTextAlign(11)
    latex.Draw()
    gobjects.append(latex)

class DrawSingleTrack:

    def __init__(self,fname,index=None,otherCanvas=None,event=(None,None,None)):
        self.fname = fname
        self.index = index
        self.otherCanvas = otherCanvas
        self.canvasTitle = None

        self.event = event
        self.runLumiEvt = None
        self.itrack = 0
        self.canvas = None
        self.icnv = 0

        self.vtxMode = None
        self.hits = ""
        self.coords = [ ]
        self.ival = None

        self.gobjects = [ ]
        

    def drawTrack(self):

        if self.otherCanvas==None:
            self.canvas = ROOT.TCanvas("cnv","cnv",900,900)
            self.canvas.Divide(3,3)
        else:
            self.canvas = self.otherCanvas
        self.icnv = 0

        self.vtxMode = None
        self.hits = ""
        self.coords = [ ]
        self.ival = None
        showEvent = True
        tf = ROOT.TFile(self.fname)
        tfd = tf.Get("trackAnalyzer")
        tree = tfd.Get("GsfTree")
        self.hits = ""
        mgses = { }
        for ev in tree:

            if ev.ic<0:
                px, py, pz = None, None, None
                self.canvasTitle = "Run "+str(ev.run)+" Lumi "+str(ev.lumi)+" Evt "+str(ev.evt)
                self.runLumiEvt = ( ev.run, ev.lumi, ev.evt )
                showEvent = ( self.event[0]==None or self.event[0]==self.runLumiEvt[0] ) and \
                    ( self.event[1]==None or self.event[1]==self.runLumiEvt[1] ) and \
                    ( self.event[2]==None or self.event[2]==self.runLumiEvt[2] )
                self.vtxMode = ""
#                fields = l[:-1].split()
#                coord = None
#                if "from" in fields[2]:
#                    self.vtxMode = fields[3][:-1]
#                else:
#                    self.vtxMode = "cartesian"
#                print "vtxMode = ",self.vtxMode
                for c in [ "px", "py", "pz", "qp", "dxdz", "dydz" ]:
                    self.coords.append(c)
                    mgses[c] = PyMultiGaussanState1D()
                continue

            if self.vtxMode!=None and ev.ic>=0:
                mgses["px"].addState(ev.wgt,ev.gsfPx,gsfSigPx)
                mgses["py"].addState(ev.wgt,ev.gsfPy,gsfSigPy)
                mgses["pz"].addState(ev.wgt,ev.gsfPz,gsfSigPz)
                mgses["qp"].addState(ev.wgt,ev.locPars(0),sqrt(locCov(0,0)))
                mgses["dxdz"].addState(ev.wgt,ev.locPars(1),sqrt(locCov(1,1)))
                mgses["dydz"].addState(ev.wgt,ev.locPars(2),sqrt(locCov(2,2)))
                if ev.ic<(ev.nc-1):
                    continue

                self.icnv = 0
                for c in [ 'px', 'py', 'pz', 'qp', 'dxdz', 'dydz' ]:
                    self.icnv += 1
                    self.canvas.cd(self.icnv)
                    if ( self.index==None or self.itrack==self.index ) and showEvent:
                        if c.startswith("p"):
                            drawMixture(ROOT.gPad,c,mgses[c],xmin=-200.,xmax=200., \
                                        superimpose=(self.otherCanvas!=None))
                            if c=='px':
                                px = mgses[c].mode()[0]
                            elif c=='py':
                                py = mgses[c].mode()[0]
                            elif c=='pz':
                                pz = mgses[c].mode()[0]
                        elif c=="qp":
                            drawMixture(ROOT.gPad,c,mgses[c],xmin=-1./2.5,xmax=1./2.5, \
                                            superimpose=(self.otherCanvas!=None))
                        else:
                            drawMixture(ROOT.gPad,c,mgses[c], \
                                            superimpose=(self.otherCanvas!=None))
                    ROOT.gPad.Update()
                if self.canvasTitle!=None:
                    title = "Track "+str(self.itrack)
                    p3 = ROOT.TVector3(px,py,pz)
                    title += " ; Pt {0:4.1f} Eta {1:5.2f}".format(p3.Perp(),p3.Eta())
                    title += " ; "+self.canvasTitle
                    self.canvas.SetTitle(title)
                if ( self.index==None or self.itrack==self.index ) and showEvent:
                    self.canvas.Update()
                    print "hits: ",self.hits
                    self.hits = ""
                    raw_input("Enter")
                if self.itrack==self.index:
                    return self.canvas
                else:
                    self.itrack += 1
                self.icnv = 0
                self.gobjects = [ ]



if __name__=="__main__":

    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument('--itrack', '-i', type=int, default=None)
    parser.add_argument('--run', '-r', type=int, default=None)
    parser.add_argument('--lumi', '-l', type=int, default=None)
    parser.add_argument('--evt', '-e', type=int, default=None)
    parser.add_argument('files',nargs='+',default=None)
    args = parser.parse_args()
    assert len(args.files)<=2 and (len(args.files)==1 or args.itrack!=None)

    drawClass = DrawSingleTrack(fname=args.files[0],index=args.itrack,event=(args.run,args.lumi,args.evt))
    cnv = drawClass.drawTrack()

    if len(args.files)>1:    
        drawClass2 = DrawSingleTrack(fname=args.files[1],index=args.itrack,otherCanvas=cnv)
        cnv2 = drawClass2.drawTrack()

