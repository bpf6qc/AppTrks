#include "ZPeakAnalyzer.h"

using namespace std;
using namespace reco;
using namespace edm;

ZPeakAnalyzer::ZPeakAnalyzer(const edm::ParameterSet &cfg) :
  recoTrack_(cfg.getParameter<edm::InputTag>("recoTrack")),
  recoTrack2_(cfg.getParameter<edm::InputTag>("recoTrack2")),
  triggerResultsLabel_(cfg.getParameter<edm::InputTag>("triggerResultsLabel"))
{
  TH1D::SetDefaultSumw2();
  h_invmass = fs->make<TH1D>("invmass", "Di-muon invariant mass;m(#mu, #mu) [GeV];Events", 75, 0, 150);
}

ZPeakAnalyzer::~ZPeakAnalyzer() {

}

ZPeakAnalyzer::analyze(const edm::Event &event, const edm::EventSetup &setup) {

  if(event.getByLabel(recoTrack_, TrackCollection)) {
    for(reco::TrackCollection::const_iterator recoTrack = TrackCollection->begin(); recoTrack != TrackCollection->end(); recoTrack++) {

      if(event.getByLabel(recoTrack2_, SecondaryTrackCollection)) {
	for(reco::TrackCollection::const_iterator recoTrack2 = SecondaryTrackCollection->begin(); recoTrack2 != SecondaryTrackCollection->end(); recoTrack2++) {

	  if(recoTrack->phi() != recoTrack2->phi() && IsGoodMuon(*recoTrack) && IsGoodMuon(*recoTrack2)) {

	    TLorentzVector muon1(recoTrack->px(), recoTrack->py(), recoTrack->pz(), recoTrack->p());
	    TLorentzVector muon2(recoTrack2->px(), recoTrack2->py(), recoTrack2->pz(), recoTrack2->p());

	    double invmass = (muon1 + muon2).M();
	    h_invmass->Fill(invmass);
	  }

	}
      }

    }
  }

}

bool ZPeakAnalyzer::IsGoodMuon(const reco::Track recoTrack) {

  bool passes = (recoTrack.pt() > 10 &&
		 fabs(recoTrack.eta()) < 2.4 &&
		 recoTrack.hitPattern().numberOfValidMuonHits() > 0 &&
		 recoTrack.hitPattern().muonStationsWithValidHits() >= 2);

  return passes;

}

DEFINE_FWK_MODULE(ZPeakAnalyzer);
