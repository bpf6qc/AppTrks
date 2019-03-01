#include "ZPeakAnalyzer.h"

using namespace std;
using namespace reco;
using namespace edm;

ZPeakAnalyzer::ZPeakAnalyzer(const edm::ParameterSet &cfg) :
  trackCollectionLabel_(cfg.getParameter<edm::InputTag>("trackCollectionLabel")),
  genCollectionLabel_(cfg.getParameter<edm::InputTag>("genCollectionLabel")),
  triggerResultsLabel_(cfg.getParameter<edm::InputTag>("triggerResultsLabel"))
{
  TH1D::SetDefaultSumw2();

  h_chargino_pt = fs->make<TH1D>("chargino_pt", "Gen chargino PT", 200, 0, 2000);
  h_chargino_eta = fs->make<TH1D>("chargino_eta", "Gen chargino eta", 100, -5, 5);

  h_sa_pt = fs->make<TH1D>("sa_pt", "StandAlone muon PT", 200, 0, 2000);
  h_sa_eta = fs->make<TH1D>("sa_eta", "StandAlone muon eta", 100, -5, 5);
  
}

ZPeakAnalyzer::~ZPeakAnalyzer() {

}

void ZPeakAnalyzer::analyze(const edm::Event &event, const edm::EventSetup &setup) {

  if(event.getByLabel(trackCollectionLabel_, TrackCollection)) {
    for(reco::TrackCollection::const_iterator recoTrack = TrackCollection->begin(); recoTrack != TrackCollection->end(); recoTrack++) {

      TLorentzVector muon(recoTrack->px(), recoTrack->py(), recoTrack->pz(), recoTrack->p());
      
      h_sa_pt->Fill(muon.Pt());
      h_sa_eta->Fill(muon.Eta());
    }
  }

  if(event.getByLabel(genCollectionLabel_, GenParticlesCollection)) {
    for(reco::GenParticleCollection::const_iterator part = GenParticlesCollection->begin(); part != GenParticlesCollection->end(); part++) {

      if(abs(part->pdgId()) != 1000024) continue;

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
