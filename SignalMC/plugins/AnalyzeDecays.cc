// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/Candidate/interface/Candidate.h"

#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"

#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/Candidate/interface/CandMatchMap.h"
#include "DataFormats/Common/interface/AssociationVector.h"

#include "SimDataFormats/Track/interface/SimTrackContainer.h"
#include "SimDataFormats/Vertex/interface/SimVertexContainer.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "TH2D.h"
#include "TH1D.h"
#include "TLorentzVector.h"

#include <ext/algorithm>  // For sorting?

using namespace edm;
using namespace std;
using namespace reco;

class AnalyzeDecays : public edm::EDAnalyzer {
public:
  explicit AnalyzeDecays(const edm::ParameterSet&);
  ~AnalyzeDecays();

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
  
private:
  virtual void beginJob() ;
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;

  virtual void beginRun(edm::Run const&, edm::EventSetup const&);
  virtual void endRun(edm::Run const&, edm::EventSetup const&);
  virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);
  virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);

  // ----------member data ---------------------------

  TH1D * h_neutralinoSumPt;

  TH1D * h_nVtxNeutralinoToChargino;
  TH1D * h_nVtxNeutralinoParent;
  TH1D * h_nGenNeutralinoTot;
  TH1D * h_nGenNeutralinoSel;

  TH1D * h_nNeutralinoNoDecay;
  TH1D * h_missingVtx;
  TH1D * h_numdgt;
  
  TH1D * h_vxy;
  TH1D * h_vx;
  TH1D * h_vy;
  TH1D * h_vz;
  TH2D * h_vxyz;
  TH1D * h_decayVxy;
  TH1D * h_decayVx;
  TH1D * h_decayVy;
  TH1D * h_decayVz;
  TH2D * h_decayVxyz;
  TH2D * h_decayVxyzWide;
  TH1D * h_decayLength;
  TH1D * h_ctauSmall;
  TH1D * h_ctauMedium;
  TH1D * h_ctauLarge;
  TH1D * h_ctauTruncatedSmall;
  TH1D * h_ctauTruncatedMedium;
  TH1D * h_ctauTruncatedLarge;

  TH1D * h_daughter0Id;
  TH1D * h_daughter0E;
  TH1D * h_daughter1Id;
  TH1D * h_daughter1E;
  TH1D * h_EDiff;

  TH1D * h_neutralinoPt;

  TH1D * h_GenEta;
  TH1D * h_GenEtaSel;
  TH1D * h_GenEtaFoundVtx;

  TH1D * h_dPhi_decay;
  TH1D * h_dR_decay;
  
  TH1D * h_reco_pt;
  TH1D * h_reco_eta;
  
  edm::InputTag genParticleTag_;
  edm::InputTag recoParticleTag_;

  bool isParticleGun_;
  double MaxEta_;
  bool quiet_;

  struct LessById {
    bool operator()(const SimTrack &tk1, const SimTrack &tk2) const { return tk1.trackId() < tk2.trackId(); }
    bool operator()(const SimTrack &tk1, unsigned int    id ) const { return tk1.trackId() < id;            }
    bool operator()(unsigned int     id, const SimTrack &tk2) const { return id            < tk2.trackId(); }
  };

};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
AnalyzeDecays::AnalyzeDecays(const edm::ParameterSet& iConfig) {

  genParticleTag_ = iConfig.getParameter<edm::InputTag>("genParticleTag");
  recoParticleTag_ = iConfig.getParameter<edm::InputTag>("recoParticleTag");
  isParticleGun_  = iConfig.getUntrackedParameter<bool>("isParticleGun", false);
  MaxEta_         = iConfig.getUntrackedParameter<double>("MaxEta", 2.5);
  quiet_          = iConfig.getUntrackedParameter<bool>("quiet", false);

  //now do what ever initialization is needed
  edm::Service<TFileService> fs;

  h_neutralinoSumPt = fs->make<TH1D>("neutralinoSumPt", ";p_{T} (#chi^{0}#chi^{0})", 100, 0, 500);
  h_neutralinoPt = fs->make<TH1D>("neutralinoPt", ";p_{T} (#chi^{0})", 100, 0, 500);
  
  h_nVtxNeutralinoToChargino = fs->make<TH1D>("nVtxNeutralinoToChargino", ";nVtx with #chi^{0}#rightarrow#chi^{#pm}", 10, -0.5, 9.5);
  h_nVtxNeutralinoParent = fs->make<TH1D>("nVtxNeutralinoParent", ";nVtx with #chi^{0} parent", 10, -0.5, 9.5);
  h_nGenNeutralinoTot = fs->make<TH1D>("nGenNeutralinoTot", "nGen #chi^{0}", 10, -0.5, 9.5);
  h_nGenNeutralinoSel = fs->make<TH1D>("nGenNeutralinoSel", "nSelected #chi^{0}", 10, -0.5, 9.5);

  h_nNeutralinoNoDecay = fs->make<TH1D>("nNeutralinoNoDecay", ";nGen #chi^{0} with no decay", 10, -0.5, 9.5);
  h_missingVtx = fs->make<TH1D>("missingVtx", ";# missing #chi^{0} vertices", 10, -0.5, 9.5);
  h_numdgt = fs->make<TH1D>("numdgt", ";nDaughters of #chi^{0}", 10, -0.5, 9.5);

  h_GenEta = fs->make<TH1D>("GenEta" , ";#eta, generated #chi^{#pm}" , 100 , -7, 7);
  h_GenEtaSel = fs->make<TH1D>("GenEtaSel" , ";#eta, selected #chi^{#pm}" , 100 , -7, 7);
  h_GenEtaFoundVtx = fs->make<TH1D>("GenEtaFoundVtx" ,    ";#eta, generated #chi^{#pm} with decay vertex" , 100 , -7, 7);

  h_vxy = fs->make<TH1D>("neutralino_vxy", "neutralino vxy", 50, 0, 0.5);
  h_vx = fs->make<TH1D>("neutralino_vx", "neutralino_vx", 50, 0, 0.5);
  h_vy = fs->make<TH1D>("neutralino_vy", "neutralino_vy", 50, 0, 0.5);
  h_vz = fs->make<TH1D>("neutralino_vz", "neutralino_vz", 50, 0, 0.5);
  h_vxyz = fs->make<TH2D>("neutralino_vxyz", ";neutralino_vz;neutralino_vxy", 50, 0, 0.5, 50, 0, 0.5);
  h_decayVxy = fs->make<TH1D>("chargino_Vxy", "chargino_Vxy", 100, 0, 1000);
  h_decayVx = fs->make<TH1D>("chargino_Vx", "chargino_Vx", 100, 0, 1000);
  h_decayVy = fs->make<TH1D>("chargino_Vy", "chargino_Vy", 100, 0, 1000);
  h_decayVz = fs->make<TH1D>("chargino_Vz", "chargino_Vz", 100, 0, 1500);
  h_decayVxyz = fs->make<TH2D>("chargino_Vxyz", "Position of #chi^{0}#rightarrow#chi^{#pm}#pi^{#mp} decay;|z| [cm];|#rho| [cm]" , 100, 0, 1500, 100, 0, 1000);
  h_decayVxyzWide = fs->make<TH2D>("chargino_VxyzWide",  "Position of #chi^{0}#rightarrow#chi^{#pm}#pi^{#mp} decay;|z| [cm];|#rho| [cm]" , 1000, 0, 3000, 1000, 0, 2000);
  h_decayLength = fs->make<TH1D>("decayLength", "decayLength", 100, 0, 1000);
  h_ctauSmall = fs->make<TH1D>("ctauSmall", ";c#tau [cm]", 100, 0, 100);
  h_ctauMedium = fs->make<TH1D>("ctauMedium", ";c#tau [cm]", 100, 0, 1000);
  h_ctauLarge = fs->make<TH1D>("ctauLarge", ";c#tau [cm]", 100, 0, 10000);
  h_ctauTruncatedSmall = fs->make<TH1D>("ctauTruncatedSmall", ";c#tau [cm]", 100, 0, 100);
  h_ctauTruncatedMedium = fs->make<TH1D>("ctauTruncatedMedium", ";c#tau [cm]", 100, 0, 1000);
  h_ctauTruncatedLarge = fs->make<TH1D>("ctauTruncatedLarge", ";c#tau [cm]", 100, 0, 10000);

  h_daughter0Id = fs->make<TH1D>("daughter0Id", ";daughter 0 PDG ID", 100, 1000020, 1000030);
  h_daughter1Id = fs->make<TH1D>("daughter1Id", ";daughter 1 PDG ID", 100, 200, 220);
  h_daughter0E  = fs->make<TH1D>("daughter0E",  ";daughter 0 energy [GeV]", 100, 0, 1000);
  h_daughter1E  = fs->make<TH1D>("daughter1E", ";daughter 1 energy [GeV]", 100, 0, 500);
  h_EDiff       = fs->make<TH1D>("EDiff", ";(total daughters' energy) - (chargino energy) [GeV]", 100, -150, 150);

  h_dPhi_decay = fs->make<TH1D>("dPhi_decay", ";#Delta#phi(#chi^{0}, #chi^{#pm})", 20, -0.005, 0.005);
  h_dR_decay = fs->make<TH1D>("dR_decay", ";#Delta R(#chi^{0}, #chi^{#pm})", 10, 0, 0.01);
  
  // RECO quantities
  h_reco_pt = fs->make<TH1D>("reco_pt", "RECO muon PT", 200, 0, 2000);
  h_reco_eta = fs->make<TH1D>("reco_eta", "RECO muon eta", 100, -5, 5);
  
}


AnalyzeDecays::~AnalyzeDecays()
{

  // do anything here that needs to be done at desctruction time
  // (e.g. close files, deallocate resources etc.)

}




//
// member functions
//

// ------------ method called for each event  ------------
void
AnalyzeDecays::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  using namespace edm;

  Handle<reco::GenParticleCollection > genParticles;
  iEvent.getByLabel(genParticleTag_,genParticles);

  TLorentzVector pTot(0,0,0,0);

  int nGenNeutralinoTot = 0;
  int nNeutralinoNoDecay = 0;
  int nGenNeutralinoSel = 0;
  int nVtxNeutralinoToChargino = 0;
   
  double vxy;
  double vx;
  double vy;
  double vz;
  double decayVxy;
  double decayVx;
  double decayVy;
  double decayVz;
  double decayLength;
  double ctau;

  quiet_ || cout << "checking gen particles, size = " << genParticles->size() << endl;

  if (genParticles.isValid()) {
    for( size_t k = 0; k < genParticles->size(); k++ ){
      const reco::Candidate & neutralino = (*genParticles)[k];

      int status = neutralino.status();
      int pdgId  = neutralino.pdgId();
      int numdgt = neutralino.numberOfDaughters();

      // since status 3 particles are always first in the list,
      // we can break once we're through these
      if(!isParticleGun_ && neutralino.status() != 3) break;
      if(abs(neutralino.pdgId()) != 1000022) continue; // only consider neutralino parents

      nGenNeutralinoTot++;

      h_GenEta->Fill(neutralino.eta());
      if (fabs(neutralino.eta()) > MaxEta_) continue;

      nGenNeutralinoSel++;
      h_GenEtaSel->Fill(neutralino.eta());

      quiet_ || cout << "Found gen neutralino with:  "
		     << " pt = " << neutralino.pt()
		     << ", eta = " << neutralino.eta()
		     << ", status=" << status
		     << "; pdgId=" << pdgId
		     << "; numdgt=" << numdgt << endl;
       
      TLorentzVector p4(neutralino.px(),
			neutralino.py(),
			neutralino.pz(),
			neutralino.energy());
      pTot += p4;

      h_neutralinoPt->Fill(p4.Pt());

      vx = neutralino.vx();
      vy = neutralino.vy();
      vz = neutralino.vz();
      vxy = sqrt(vx*vx + vy*vy);

      const reco::Candidate *mother   = &neutralino;
      const reco::Candidate *daughter = &neutralino;

      // Descend the decay chain until no daughters have the same PDG ID as neutralino.
      while(true) {
	bool foundDauSamePdgId = false;
	for(uint i = 0; i < daughter->numberOfDaughters(); i++) {
	  if(daughter->daughter(i)->pdgId() == neutralino.pdgId()) {
	    foundDauSamePdgId = true;
	    mother = daughter;
	    daughter = daughter->daughter(i);
	    break;
	  }
	}
	if (!foundDauSamePdgId) break;
      }

      const reco::Candidate * pion = 0;
      const reco::Candidate * chargino = 0;
      
      for(uint i = 0; i < daughter->numberOfDaughters(); i++) {
	if(abs(daughter->daughter(i)->pdgId()) == 211 && !pion) pion = daughter->daughter(i);
	if(abs(daughter->daughter(i)->pdgId()) == 1000024 && !chargino) chargino = daughter->daughter(i);
      }

      if(!pion || !chargino) continue;
      
      h_GenEtaFoundVtx->Fill(neutralino.eta());
      nVtxNeutralinoToChargino++;

      double daughter0Id = -99;
      double daughter0E  = -99;
      double daughter1Id = -99;
      double daughter1E  = -99;
      double daughterTotE  = 0;

      decayVx = chargino->vx ();
      decayVy = chargino->vy ();
      decayVz = chargino->vz ();
      decayVxy = sqrt(decayVx*decayVx + decayVy*decayVy);

      TVector3 source(neutralino.vx (), neutralino.vy (), neutralino.vz ());
      TVector3 sink(chargino->vx (), chargino->vy (), chargino->vz ());

      quiet_ || cout << "Debug:  BNstop particle pdg id = " << neutralino.pdgId()
		     << ", status = " << neutralino.status()
		     << ", pt = " << neutralino.pt()
		     << ", eta = " << neutralino.eta()
		     << ", num charginos = " << neutralino.numberOfDaughters()
		     << endl
		     << "  production:  "
		     << "  vx = " << neutralino.vx ()
		     << "  vy = " << neutralino.vy ()
		     << "  vz = " << neutralino.vz ()
		     << endl
		     << "  decay:  "
		     << "  vx = " << chargino->vx ()
		     << "  vy = " << chargino->vy ()
		     << "  vz = " << chargino->vz ()
		     << endl
		     << "  chargino PDG ID:  " << chargino->pdgId()
		     << ", status = " << chargino->status()
		     << ", pt = " << chargino->pt()
		     << ", eta = " << chargino->eta()
		     << endl
		     << "  mother PDG ID:  " << mother->pdgId()
		     << ", status:  " << mother->status()
		     << ", numdau:  " << mother->numberOfDaughters()
		     << endl;
	 
      for (uint i=0; i<mother->numberOfDaughters(); i++) {
	quiet_ || cout << "  Mother has daughter " << i << ": " << mother->daughter(i)->pdgId() << endl;
      }

      if (source == sink) {
	// Set to non-physical values if no daughters are found
	decayLength  = -99;
	ctau         = -99;
	//          betaAtDecay  = -99;
	//          gammaAtDecay = -99;
	decayVx  = -99;
	decayVy  = -99;
	decayVz  = -99;
	decayVxy = -99;
      }
      else {
	decayLength = (sink - source).Mag ();
	ctau = (sink - source).Mag () / (neutralino.p4 ().Beta () * neutralino.p4 ().Gamma ());
      }

      h_vxy->Fill(fabs(vxy));
      h_vx->Fill(fabs(vx));
      h_vy->Fill(fabs(vy));
      h_vz->Fill(fabs(vz));
      h_vxyz->Fill(fabs(vz), fabs(vxy));
      h_decayVxy->Fill(fabs(decayVxy));
      h_decayVx->Fill(fabs(decayVx));
      h_decayVy->Fill(fabs(decayVy));
      h_decayVz->Fill(fabs(decayVz));
      h_decayVxyz->Fill(fabs(decayVz), fabs(decayVxy));
      h_decayVxyzWide->Fill(fabs(decayVz), fabs(decayVxy));
      h_decayLength->Fill(decayLength);
      h_ctauSmall->Fill(ctau);
      h_ctauMedium->Fill(ctau);
      h_ctauLarge->Fill(ctau);
      (decayVxy < 17.5 && decayVz < 27.0) && h_ctauTruncatedSmall->Fill(ctau);
      (decayVxy < 175.0 && decayVz < 270.0) && h_ctauTruncatedMedium->Fill(ctau);
      h_ctauTruncatedLarge->Fill(ctau);
      h_numdgt->Fill(numdgt);

      if (decayLength == 0) {
	quiet_ || cout << "Warning:  found particle " << pdgId
		       << " with decayLength = " << decayLength
		       << ", status = " << status
		       << ", numdgt = " << numdgt
		       << endl;
	nNeutralinoNoDecay++;
      }

      if(neutralino.numberOfDaughters() >= 1){
	daughter0Id = neutralino.daughter(0)->pdgId();
	daughter0E  = neutralino.daughter(0)->energy();
	daughterTotE += daughter0E;
	quiet_ || cout << "Debug:  daughter0Id = " << neutralino.daughter(0)->pdgId()
		       << ", status = " << neutralino.daughter(0)->status()
		       << ", pt = " <<     neutralino.daughter(0)->pt()
		       << ", eta = " <<    neutralino.daughter(0)->eta()
		       << ", vx = " << neutralino.daughter(0)->vx()
		       << ", vy = " << neutralino.daughter(0)->vy()
		       << ", vz = " << neutralino.daughter(0)->vz()
		       << endl;


	if(neutralino.numberOfDaughters() >= 2){
	  daughter1Id = neutralino.daughter(1)->pdgId();
	  daughter1E  = neutralino.daughter(1)->energy();
	  daughterTotE += daughter1E;

	  quiet_ || cout << "Debug:  daughter1Id = " << neutralino.daughter(1)->pdgId()
			 << ", status = " << neutralino.daughter(1)->status()
			 << ", pt = " <<        neutralino.daughter(1)->pt()
			 << ", eta = " <<       neutralino.daughter(1)->eta()
			 << ", vx = " << neutralino.daughter(1)->vx()
			 << ", vy = " << neutralino.daughter(1)->vy()
			 << ", vz = " << neutralino.daughter(1)->vz()
			 << endl;

	}
      }

      h_daughter0Id->Fill(daughter0Id);
      h_daughter0E ->Fill(daughter0E);
      h_daughter1Id->Fill(fabs(daughter1Id));
      h_daughter1E ->Fill(daughter1E);
      if (decayLength!=0) {
	h_EDiff      ->Fill(daughterTotE - neutralino.energy());
      }

      h_dPhi_decay->Fill(TVector2::Phi_mpi_pi(neutralino.phi() - chargino->phi()));
      h_dR_decay->Fill(TMath::Sqrt((neutralino.eta() - chargino->eta())*(neutralino.eta() - chargino->eta()) +
				   (neutralino.phi() - chargino->phi())*(neutralino.phi() - chargino->phi())));
       
      quiet_ || cout << "Debug:  daughter0E = " << daughter0E
		     << ", daughter1E = " << daughter1E
		     << ", mother EE = " << neutralino.energy()
		     << ", mother status = " << neutralino.status()
		     << ", |daughter0Id| = 1e6 + " << fabs(daughter0Id) - 1e6
		     << ", daughter0 status = " << neutralino.daughter(0)->status()
		     << ", daughter1Id = " << daughter1Id
		     << ", Ediff = " << daughterTotE - neutralino.energy()
		     << endl;
    }

  }

  h_neutralinoSumPt->Fill(pTot.Pt());
  h_nVtxNeutralinoToChargino->Fill(nVtxNeutralinoToChargino);
  h_nGenNeutralinoTot->Fill(nGenNeutralinoTot);
  h_nGenNeutralinoSel->Fill(nGenNeutralinoSel);
  h_nNeutralinoNoDecay->Fill(nNeutralinoNoDecay);
  h_missingVtx->Fill(nGenNeutralinoSel - nVtxNeutralinoToChargino);
   
  quiet_ || cout << "Found"
		 << " nVtxCharginoToNeutralino = " << nVtxNeutralinoToChargino
		 << ", missingVtx = " << nGenNeutralinoSel - nVtxNeutralinoToChargino
		 << ", nGenCharginoTot = " << nGenNeutralinoTot
		 << ", nGenCharginoSel = " << nGenNeutralinoSel
		 << endl;

  quiet_ || cout << "Debugging:  missingVtx = " << nGenNeutralinoSel - nVtxNeutralinoToChargino << endl;

  Handle<reco::TrackCollection > recoParticles;
  iEvent.getByLabel(recoParticleTag_, recoParticles);
   
  if(recoParticles.isValid()) {
    
    for(size_t k = 0; k < recoParticles->size(); k++){
      const reco::Track & recoTrack = (*recoParticles)[k];

      TLorentzVector muon(recoTrack.px(), recoTrack.py(), recoTrack.pz(), recoTrack.p());
       
      h_reco_pt->Fill(muon.Pt());
      h_reco_eta->Fill(muon.Eta());

    }

  }

#ifdef THIS_IS_AN_EVENT_EXAMPLE
  Handle<ExampleData> pIn;
  iEvent.getByLabel("example",pIn);
#endif

#ifdef THIS_IS_AN_EVENTSETUP_EXAMPLE
  ESHandle<SetupData> pSetup;
  iSetup.get<SetupRecord>().get(pSetup);
#endif
}


// ------------ method called once each job just before starting event loop  ------------
void
AnalyzeDecays::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void
AnalyzeDecays::endJob()
{
}

// ------------ method called when starting to processes a run  ------------
void
AnalyzeDecays::beginRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a run  ------------
void
AnalyzeDecays::endRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when starting to processes a luminosity block  ------------
void
AnalyzeDecays::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a luminosity block  ------------
void
AnalyzeDecays::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
AnalyzeDecays::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(AnalyzeDecays);
