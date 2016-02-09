# import ROOT in batch mode
import sys
oldargv = sys.argv[:]
sys.argv = [ '-b-' ]
import ROOT
ROOT.gROOT.SetBatch(True)
sys.argv = oldargv
import operator

# load FWLite C++ libraries
ROOT.gSystem.Load("libFWCoreFWLite.so");
ROOT.gSystem.Load("libDataFormatsFWLite.so");
ROOT.AutoLibraryLoader.enable()

# load FWlite python libraries
from DataFormats.FWLite import Handle, Events

triggerBits, triggerBitLabel = Handle("edm::TriggerResults"), ("TriggerResults","","HLT")
triggerPrescales, triggerPrescaleLabel  = Handle("pat::PackedTriggerPrescales"), "patTrigger"

# open file (you can use 'edmFileUtil -d /store/whatever.root' to get the physical file name)
#events = Events("file:/data/users/hart/condor/AMSB_chargino500GeV_ctau10cm_step4_User/AMSB_chargino_step4_24.root")
events = Events("file:toymodel_250gev_10cm.root")

maxNEvts = -1

#dict {'HLT_a': 1, 'HLT_b': 2, ...}
triggerBinNumbers = {}
numTriggers = 0
triggerFireCounts = {}
            
# set up histogram
output = ROOT.TFile("triggerFires.root", "RECREATE")
histogram = ROOT.TH1D("h_hlt", "HLT fires", 500, 1, 501)
            
for iev,event in enumerate(events):
    event.getByLabel(triggerBitLabel, triggerBits)
    names = event.object().triggerNames(triggerBits.product())

    for i in xrange(triggerBits.product().size()):

        # if this is the first encounter of a trigger path,
        # count it in numTriggers and add its name to the list with 0 fires
        if names.triggerName(i) not in triggerBinNumbers:
            numTriggers += 1
            triggerBinNumbers[names.triggerName(i)] = numTriggers # bin number
            triggerFireCounts[names.triggerName(i)] = 0 # times fired
            histogram.GetXaxis().SetBinLabel(numTriggers, names.triggerName(i))

        # find 
        binNumber = triggerBinNumbers[names.triggerName(i)]

        if triggerBits.product().accept(i):
            triggerFireCounts[names.triggerName(i)] += 1
            histogram.Fill(binNumber)
                        
    if maxNEvts > 0 and iev >= maxNEvts-1:
        break


for x, y in sorted(triggerFireCounts.items(), key=operator.itemgetter(1), reverse=True):
    print x, y
    
output.Write()
output.Close()
