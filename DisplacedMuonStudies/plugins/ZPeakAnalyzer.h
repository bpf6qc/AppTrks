#include <string>

#include <iostream>
#include <algorithm>
#include <regex>
#include "math.h"
#include "TH1.h"
#include "TH2.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TTree.h"
#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/HLTReco/interface/TriggerEvent.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "DataFormats/HLTReco/interface/TriggerObject.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "DataFormats/Math/interface/deltaPhi.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "TVector3.h"
#include "TLorentzVector.h"

using namespace std;

class ZPeakAnalyzer : public edm::EDAnalyzer
{
 public:
  ZPeakAnalyzer(const edm::ParameterSet &);
  ~ZPeakAnalyzer();

  void analyze(const edm::Event &, const edm::EventSetup &);
  bool IsGoodMuon(const reco::Track recoTrack);

  edm::Service<TFileService> fs;

 private:
  edm::Handle<reco::TrackCollection> TrackCollection;
  edm::Handle<reco::TrackCollection> SecondaryTrackCollection;
  edm::Handle<trigger::TriggerEvent> triggerSummary;

  edm::InputTag recoTrack_;
  edm::InputTag recoTrack2_;
  edm::InputTag triggerResultsLabel_;

  edm::Handle<edm::TriggerResults> triggerResults;

  TH1D * h_invmass;
};
