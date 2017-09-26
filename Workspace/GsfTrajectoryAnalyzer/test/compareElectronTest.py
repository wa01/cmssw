import sys,os

class Electron:

    def __init__(self,run,lumi,evt,pt,eta,phi):
        self.event = ( run, lumi, evt )
        self.pt = pt
        self.eta = eta
        self.phi = phi

        self.scE = None
        self.scEcorr = None
        self.scEta = None

        self.lostHits = None
        self.validHits = None
        self.hitPattern = None

    def __eq__(self,other):
#        print "Compare"
#        print "  ",self.event,other.event
#        print "  ",self.lostHits,other.lostHits
#        print "  ",self.validHits,other.validHits
#        print "  ",self.hitPattern,other.hitPattern
        return self.event==other.event and \
            self.lostHits==other.lostHits and \
            self.validHits==other.validHits and \
            self.hitPattern==other.hitPattern

        
def readElectrons(fname):
    event = None
    electron = None
    result = { }
    for l in open(fname):
        fields = l[:-1].strip().split()
        if len(fields)==3 and fields[0].isdigit() and fields[1].isdigit() and fields[2].isdigit():
            run = int(fields[0])
            lumi = int(fields[1])
            evt = int(fields[2])
            event = ( run, lumi, evt )
            if not event in result:
                result[event] = [ ]
            electron = None
        elif l.startswith("electron pt,eta,phi"):
            electron = Electron(run,lumi,evt,float(fields[3]),float(fields[4]),float(fields[5]))
            result[event].append(electron)
        elif l.startswith("gsfTrack lostHits"):
            electron.lostHits = int(fields[7])
            electron.validHits = int(fields[8])
        elif l.startswith("GsfTrack hits"):
            electron.hitPattern = fields[2:]
        elif l.find(" std /  corrected cluster energy")>=0:
            electron.scE = float(fields[7])
            electron.scEcorr = float(fields[8])
            electron.scEta = float(fields[9])
    return result

res1 = readElectrons(sys.argv[1])
res2 = readElectrons(sys.argv[2])

#assert len(res1)==len(res2)

evts1 = sorted(res1.keys())
evts2 = sorted(res2.keys())
assert evts1==evts2

for evt in evts1:
    print "Event ",evt
    eles1 = res1[evt]
    eles2 = res2[evt]
    for e1 in eles1:
        nmatch = 0
        for e2 in eles2:
            if e2==e1:
                nmatch += 1
                if abs(e1.eta-e2.eta)>0.1:
                    line = "* "
                else:
                    line = "  "
                line += "  pts   :  {0:5.1f} {1:5.1f}".format(e1.pt,e2.pt)
                line += "  etas  :  {0:5.2f} {1:5.2f}".format(e1.eta,e2.eta)
                line += "  phis  :  {0:5.2f} {1:5.2f}".format(e1.phi,e2.phi)
                line += "  pts   :  {0:5.1f} {1:5.1f}".format(e1.pt,e2.pt)
                line += "  scE   :  {0:5.1f} {1:5.1f}".format(e1.scE,e2.scE)
                line += "  scEta :  {0:5.2f} {1:5.2f}".format(e1.scEta,e2.scEta)
                print line
                print "     ",e1.hitPattern
        if nmatch==0:
            print "*** no match ***"
        elif nmatch>1:
            print "*** multiple matches ***"

