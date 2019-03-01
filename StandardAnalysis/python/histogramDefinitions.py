import FWCore.ParameterSet.Config as cms

###############################################
##### Set up the histograms to be plotted #####
###############################################

newMuonHistograms = cms.PSet(
    inputCollection = cms.vstring("muons"),
    histograms = cms.VPSet(
        # promptFinalState
        cms.PSet(
            binsX = cms.untracked.vdouble(100, 0, 100),
            inputVariables = cms.vstring('abs(genMatchedParticle.promptFinalState.pdgId) - 1000000'),
            name = cms.string('genMatchedPromptFinalStateID'),
            title = cms.string(';PDG ID - 1e6')
        ),
        cms.PSet(
            binsX = cms.untracked.vdouble(30, 0.0, 0.3),
            inputVariables = cms.vstring('dRToGenMatchedParticle.promptFinalState'),
            name = cms.string('dRToGenMatchedPromptFinalStateID'),
            title = cms.string(';#DeltaR')
        ),
        cms.PSet(
            binsX = cms.untracked.vdouble(200, 0, 1000),
            inputVariables = cms.vstring('genMatchedParticle.promptFinalState.pt'),
            name = cms.string('genMatchedPromptFinalStatePt'),
            title = cms.string(';genMatchedPromptFinalState Pt')
        ),
        # hardProcessFinalState
        cms.PSet(
            binsX = cms.untracked.vdouble(100, 0, 100),
            inputVariables = cms.vstring('abs(genMatchedParticle.hardProcessFinalState.pdgId) - 1000000'),
            name = cms.string('genMatchedHardProcessFinalStateID'),
            title = cms.string(';PDG ID - 1e6')
        ),
        cms.PSet(
            binsX = cms.untracked.vdouble(30, 0.0, 0.3),
            inputVariables = cms.vstring('dRToGenMatchedParticle.hardProcessFinalState'),
            name = cms.string('dRToGenMatchedHardProcessFinalStateID'),
            title = cms.string(';#DeltaR')
        ),
        cms.PSet(
            binsX = cms.untracked.vdouble(200, 0, 1000),
            inputVariables = cms.vstring('genMatchedParticle.hardProcessFinalState.pt'),
            name = cms.string('genMatchedHardProcessFinalStatePt'),
            title = cms.string(';genMatchedHardProcessFinalState Pt')
        ),
        cms.PSet(
            binsX = cms.untracked.vdouble(210, -10, 200),
            inputVariables = cms.vstring('(pt - genMatchedParticle.promptFinalState.pt) / genMatchedParticle.promptFinalState.pt'),
            name = cms.string('genMatchedPtDifference'),
            title = cms.string(';(muon pt - genMatched pt) / genMatched pt')
        ),
        cms.PSet(
            binsX = cms.untracked.vdouble(2000, -1000, 1000),
            inputVariables = cms.vstring('pt - genMatchedParticle.promptFinalState.pt'),
            name = cms.string('genMatchedPtDifference'),
            title = cms.string(';muon pt - genMatched pt')
        ),
    )
)

muonD0Histograms = cms.PSet(
    inputCollection = cms.vstring("muons", "beamspots"),
    histograms = cms.VPSet (
        cms.PSet(
            name = cms.string("muonD0Beamspot"),
            title = cms.string("Muon d_{0} wrt Beamspot; muon d_{0} [cm]"),
            binsX = cms.untracked.vdouble(1000, -20, 20),
            inputVariables = cms.vstring("(-(muon.vx - beamspot.x0)*muon.py + (muon.vy - beamspot.y0)*muon.px)/muon.pt"),
        ),
        cms.PSet(
            name = cms.string("muonD0SigBeamspot"),
            title = cms.string("Muon d_{0}/#sigma(d_{0}) wrt Beamspot; d_{0}/#sigma(d_{0})"),
            binsX = cms.untracked.vdouble(1000, -20, 20),
            inputVariables = cms.vstring("((-(muon.vx - beamspot.x0)*muon.py + (muon.vy - beamspot.y0)*muon.px)/muon.pt)/hypot(muon.innerTrack.d0Error, hypot(beamspot.x0Error, beamspot.y0Error))"),
        ),
        cms.PSet(
            name = cms.string("muonAbsD0BeamspotVsAbsDz"),
            title = cms.string("Muon |d_{0}| wrt Beamspot vs. Muon |d_{z}|; muon |d_{z}| [cm]; |muon d_{0}| [cm]"),
            binsX = cms.untracked.vdouble(1000, 0, 20),
            binsY = cms.untracked.vdouble(1000, 0, 20),
            inputVariables = cms.vstring("abs((muon.vz - beamspot.z0) - ((muon.vx - beamspot.x0)*muon.px + (muon.vy - beamspot.y0)*muon.py)/muon.pt*(muon.pz/muon.pt))","abs((-(muon.vx - beamspot.x0)*muon.py + (muon.vy - beamspot.y0)*muon.px)/muon.pt)"),
        ),
    )
)

muonInnerTrackHistograms = cms.PSet(
    inputCollection = cms.vstring("muons"),
    histograms = cms.VPSet(
        cms.PSet(
            name = cms.string("nLostStripHits"),
            title = cms.string("Number of lost strip hits; nLost"),
            binsX = cms.untracked.vdouble(30, 0, 30),
            inputVariables = cms.vstring("innerTrack.hitPattern_.numberOfLostStripHits"),
        ),
        cms.PSet(
            name = cms.string("nLostPixelHits"),
            title = cms.string("Number of lost pixel hits; nLost"),
            binsX = cms.untracked.vdouble(30, 0, 30),
            inputVariables = cms.vstring("innerTrack.hitPattern_.numberOfLostPixelHits"),
        ),
        cms.PSet(
            name = cms.string("nValidPixelHits"),
            title = cms.string("Number of valid pixel hits; nValid"),
            binsX = cms.untracked.vdouble(30, 0, 30),
            inputVariables = cms.vstring("innerTrack.hitPattern_.numberOfValidPixelHits"),
        ),
        cms.PSet(
            name = cms.string("nValidStripHits"),
            title = cms.string("Number of valid strip hits; nValid"),
            binsX = cms.untracked.vdouble(30, 0, 30),
            inputVariables = cms.vstring("innerTrack.hitPattern_.numberOfValidStripHits"),
        ),
        cms.PSet (
            name = cms.string("trackerLayersWithMeasurement"),
            title = cms.string("Muon Number of Tracker Layer with Measurement;muon trackerLayersWithMeasurement"),
            binsX = cms.untracked.vdouble(30, 0, 30),
            inputVariables = cms.vstring("muon.innerTrack.hitPattern_.trackerLayersWithMeasurement"),
        ),
    )
)

muonOuterTrackHistograms = cms.PSet(
    inputCollection = cms.vstring("muons"),
    histograms = cms.VPSet(
        cms.PSet (
            name = cms.string("nValidMuonHits"),
            title = cms.string("Muon Number of Valid Muon Hits;muon numberOfValidMuonHits"),
            binsX = cms.untracked.vdouble(30, 0, 30),
            inputVariables = cms.vstring("muon.outerTrack.hitPattern_.numberOfValidMuonHits"),
        ),
        cms.PSet (
            name = cms.string("nLostMuonHits"),
            title = cms.string("Muon Number of Lost Muon Hits;nLost"),
            binsX = cms.untracked.vdouble(30, 0, 30),
            inputVariables = cms.vstring("muon.outerTrack.hitPattern_.numberOfLostMuonHits"),
        ),
    )
)

muonHistograms = cms.PSet(
    inputCollection = cms.vstring("muons","beamspots"),
    histograms = cms.VPSet (
        cms.PSet (
            name = cms.string("muonPt"),
            title = cms.string("Muon Transverse Momentum;muon p_{T} [GeV]"),
            binsX = cms.untracked.vdouble(1000, 0, 1000),
            inputVariables = cms.vstring("muon.pt"),
        ),
        cms.PSet (
            name = cms.string("muonEta"),
            title = cms.string("Muon Pseudorapidity;muon #eta"),
            binsX = cms.untracked.vdouble(50, -2.5, 2.5),
            inputVariables = cms.vstring("muon.eta"),
        ),
        cms.PSet (
            name = cms.string("muonPhi"),
            title = cms.string("Muon Azimuthal Angle;muon #phi"),
            binsX = cms.untracked.vdouble(64, -3.2, 3.2),
            inputVariables = cms.vstring("muon.phi"),
        ),
        cms.PSet (
            name = cms.string("muonPFMuonFlag"),
            title = cms.string("Muon PFMuonFlas;muon isPFMuon"),
            binsX = cms.untracked.vdouble(4, -2, 2),
            inputVariables = cms.vstring("muon.isPFMuon"),
        ),
        cms.PSet (
            name = cms.string("muonNumberOfValidMuonHits"),
            title = cms.string("Muon Number of Valid Muon Hits;muon numberOfValidMuonHits"),
            binsX = cms.untracked.vdouble(30, 0, 30),
            inputVariables = cms.vstring("muon.globalTrack.hitPattern_.numberOfValidMuonHits"),
        ),
        cms.PSet (
            name = cms.string("muonNormalizedChi2"),
            title = cms.string("Muon normalizedChi2; muon #Chi^{2}_{Norm}"),
            binsX = cms.untracked.vdouble(30, 0, 30),
            inputVariables = cms.vstring("muon.globalTrack.normalizedChi2"),
        ),
        cms.PSet (
            name = cms.string("muonNumberOfMatchedStations"),
            title = cms.string("Muon Number of Matched Stations;muon numberOfMatchedStations"),
            binsX = cms.untracked.vdouble(30, 0, 30),
            inputVariables = cms.vstring("muon.numberOfMatchedStations"),
        ),
        cms.PSet (
            name = cms.string("muonSumChargedHadronPt"),
            title = cms.string("Muon sumChargedHadronPt;muon sumChargedHadronPt [GeV]"),
            binsX = cms.untracked.vdouble(50, 0, 50),
            inputVariables = cms.vstring("muon.pfIsolationR04_.sumChargedHadronPt"),
        ),
        cms.PSet (
            name = cms.string("muonSumNeutralHadronEt"),
            title = cms.string("Muon sumNeutralHadronEt;muon sumNeutralHadronEt [GeV]"),
            binsX = cms.untracked.vdouble(50, 0, 50),
            inputVariables = cms.vstring("muon.pfIsolationR04_.sumNeutralHadronEt"),
        ),
        cms.PSet (
            name = cms.string("muonSumPhotonEt"),
            title = cms.string("Muon sumPhotonEt;muon sumPhotonEt [GeV]"),
            binsX = cms.untracked.vdouble(50, 0, 50),
            inputVariables = cms.vstring("muon.pfIsolationR04_.sumPhotonEt"),
        ),
        cms.PSet (
            name = cms.string("muonSumPUPt"),
            title = cms.string("Muon sumPUPt;muon sumPUPt [GeV]"),
            binsX = cms.untracked.vdouble(50, 0, 50),
            inputVariables = cms.vstring("muon.pfIsolationR04_.sumPUPt"),
        ),
        cms.PSet (
            name = cms.string("muonDxy"),
            title = cms.string("Muon IP;muon d_{xy}"),
            binsX = cms.untracked.vdouble(200, -0.01, 0.01),
            inputVariables = cms.vstring("(-(muon.vx - beamspot.x0)*muon.py + (muon.vy - beamspot.y0)*muon.px)/muon.pt"),
        ),
        cms.PSet (
            name = cms.string("absMuonDxyPrompt"),
            title = cms.string("Muon IP;|muon d_{xy}| [cm]"),
            binsX = cms.untracked.vdouble(100, 0, 0.01),
            inputVariables = cms.vstring("abs((-(muon.vx - beamspot.x0)*muon.py + (muon.vy - beamspot.y0)*muon.px))/muon.pt"),
        ),
        cms.PSet (
            name = cms.string("absMuonDxyDisplaced"),
            title = cms.string("Muon IP;|muon d_{xy}| [cm]"),
            binsX = cms.untracked.vdouble(100, 0.01, 0.02),
            inputVariables = cms.vstring("abs((-(muon.vx - beamspot.x0)*muon.py + (muon.vy - beamspot.y0)*muon.px))/muon.pt"),
        ),
        cms.PSet (
            name = cms.string("absMuonDxyInclusive"),
            title = cms.string("Muon IP;|muon d_{xy}| [cm]"),
            binsX = cms.untracked.vdouble(500, 0, 0.5),
            inputVariables = cms.vstring("abs((-(muon.vx - beamspot.x0)*muon.py + (muon.vy - beamspot.y0)*muon.px))/muon.pt"),
        ),
        cms.PSet (
            name = cms.string("muonDz"),
            title = cms.string("Muon Dz;muon d_{z}"),
            binsX = cms.untracked.vdouble(200, -10, 10),
            inputVariables = cms.vstring("(muon.vz - beamspot.z0) - ((muon.vx - beamspot.x0)*muon.px + (muon.vy - beamspot.y0)*muon.py)/muon.pt*(muon.pz/muon.pt)"),
        ),
        cms.PSet (
            name = cms.string("muonDxySignificanceSmall"),
            title = cms.string("Muon Dxy Significance; muon d_{xy}/#delta_{d_{xy}}"),
            binsX = cms.untracked.vdouble(400, -2, 2),
            inputVariables = cms.vstring('((-(muon.vx - beamspot.x0)*muon.py + (muon.vy - beamspot.y0)*muon.px)/muon.pt)/hypot(muon.innerTrack.d0Error, hypot(beamspot.x0Error, beamspot.y0Error))'),
        ),
        cms.PSet (
            name = cms.string("muonAbsDxySignificanceSmall"),
            title = cms.string("Muon Dxy Significance; muon d_{xy}/#delta_{d_{xy}}"),
            binsX = cms.untracked.vdouble(200, 0, 2),
            inputVariables = cms.vstring('abs(((-(muon.vx - beamspot.x0)*muon.py + (muon.vy - beamspot.y0)*muon.px)/muon.pt)/hypot(muon.innerTrack.d0Error, hypot(beamspot.x0Error, beamspot.y0Error)))'),
        ),
        cms.PSet (
            name = cms.string("muonDxySignificanceMedium"),
            title = cms.string("Muon Dxy Significance; muon d_{xy}/#delta_{d_{xy}}"),
            binsX = cms.untracked.vdouble(100, -5, 5),
            inputVariables = cms.vstring('((-(muon.vx - beamspot.x0)*muon.py + (muon.vy - beamspot.y0)*muon.px)/muon.pt)/hypot(muon.innerTrack.d0Error, hypot(beamspot.x0Error, beamspot.y0Error))'),
        ),
        cms.PSet (
            name = cms.string("muonAbsDxySignificanceMedium"),
            title = cms.string("Muon Dxy Significance; muon d_{xy}/#delta_{d_{xy}}"),
            binsX = cms.untracked.vdouble(50, 0, 5),
            inputVariables = cms.vstring('abs(((-(muon.vx - beamspot.x0)*muon.py + (muon.vy - beamspot.y0)*muon.px)/muon.pt)/hypot(muon.innerTrack.d0Error, hypot(beamspot.x0Error, beamspot.y0Error)))'),
        ),
        cms.PSet (
            name = cms.string("muonDxySignificanceLarge"),
            title = cms.string("Muon Dxy Significance; muon d_{xy}/#delta_{d_{xy}}"),
            binsX = cms.untracked.vdouble(100, -10, 10),
            inputVariables = cms.vstring('((-(muon.vx - beamspot.x0)*muon.py + (muon.vy - beamspot.y0)*muon.px)/muon.pt)/hypot(muon.innerTrack.d0Error, hypot(beamspot.x0Error, beamspot.y0Error))'),
        ),
        cms.PSet (
            name = cms.string("muonAbsDxySignificanceLarge"),
            title = cms.string("Muon Dxy Significance; muon d_{xy}/#delta_{d_{xy}}"),
            binsX = cms.untracked.vdouble(50, 0, 10),
            inputVariables = cms.vstring('abs(((-(muon.vx - beamspot.x0)*muon.py + (muon.vy - beamspot.y0)*muon.px)/muon.pt)/hypot(muon.innerTrack.d0Error, hypot(beamspot.x0Error, beamspot.y0Error)))'),
        ),
        cms.PSet (
            name = cms.string("muonDbetaIsolation"),
            title = cms.string("Muon Isolation; muon #Delta#beta Isolation"),
            binsX = cms.untracked.vdouble(150, 0, 1.5),
            inputVariables = cms.vstring("(muon.pfIsolationR04_.sumChargedHadronPt + max(0.0,muon.pfIsolationR04_.sumNeutralHadronEt + muon.pfIsolationR04_.sumPhotonEt - 0.5*muon.pfIsolationR04_.sumPUPt))/muon.pt"),
        ),
        cms.PSet (
            name = cms.string("muonDxyPhi"),
            title = cms.string("Muon IP vs Azimuthal Angle;muon #phi"),
            binsX = cms.untracked.vdouble(64, -3.2, 3.2),
            binsY = cms.untracked.vdouble(200, -0.01, 0.01),
            inputVariables = cms.vstring("muon.phi","(-(muon.vx - beamspot.x0)*muon.py + (muon.vy - beamspot.y0)*muon.px)/muon.pt"),
        ),
        cms.PSet (
            name = cms.string("muonEtaPt"),
            title = cms.string("Muon Pseudorapidity vs Transverse Momentum;muon p_{T} [GeV]"),
            binsX = cms.untracked.vdouble(100,0,500),
            binsY = cms.untracked.vdouble(5, -2.5, 2.5),
            inputVariables = cms.vstring("muon.pt","muon.eta"),
        ),
  )
)

dimuonHistograms = cms.PSet(
    inputCollection = cms.vstring("muons", "muons"),
    histograms = cms.VPSet(
        cms.PSet(
            name = cms.string("dimuonInvmass"),
            title = cms.string("Muon - Muon invariant mass"),
            binsX = cms.untracked.vdouble(1000, 0, 1000),
            inputVariables = cms.vstring("invMass(muon, muon)"),
        ),
    )
)
    
muonJetHistograms = cms.PSet(
    inputCollection = cms.vstring("muons","jets"),
    histograms = cms.VPSet (
        cms.PSet (
            name = cms.string("muonJetDeltaR"),
            title = cms.string("Muon - Jet #DeltaR"),
            binsX = cms.untracked.vdouble(300, 0, 3),
            inputVariables = cms.vstring("deltaR(muon,jet)"),
        ),
    )
)

jetHistograms = cms.PSet(
    inputCollection = cms.vstring("jets"),
    histograms = cms.VPSet (
        cms.PSet (
            name = cms.string("jetPt"),
            title = cms.string("Jet Transverse Momentum;electron p_{T} [GeV]"),
            binsX = cms.untracked.vdouble(300, 0, 300),
            inputVariables = cms.vstring("pt"),
        ),
        cms.PSet (
            name = cms.string("jetEta"),
            title = cms.string("Jet Pseudorapidity;electron #eta"),
            binsX = cms.untracked.vdouble(50, -2.5, 2.5),
            inputVariables = cms.vstring("eta"),
        ),
        cms.PSet (
            name = cms.string("jetPhi"),
            title = cms.string("Jet Azimuthal Angle;electron #phi"),
            binsX = cms.untracked.vdouble(64, -3.2, 3.2),
            inputVariables = cms.vstring("phi"),
        ),
        cms.PSet (
            name = cms.string("jetCSV"),
            title = cms.string("Jet CSV;jet CSV"),
            binsX = cms.untracked.vdouble(1000, 0, 1),
            inputVariables = cms.vstring("pfCombinedInclusiveSecondaryVertexV2BJetTags"),
        )
    )
)

metHistograms = cms.PSet(
    inputCollection = cms.vstring("mets"),
    histograms = cms.VPSet (
        cms.PSet (
            name = cms.string("metPt"),
            title = cms.string("Met Transverse Momentum;p_{T}_{jet} [GeV]"),
            binsX = cms.untracked.vdouble(500, 0, 500),
            inputVariables = cms.vstring("pt"),
        ),
        cms.PSet (
            name = cms.string("metPhi"),
            title = cms.string("Met Azumithal Angle; #phi"),
            binsX = cms.untracked.vdouble(64, -3.2, 3.2),
            inputVariables = cms.vstring("phi"),
        ),
    )
)


