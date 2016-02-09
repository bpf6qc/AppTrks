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

#include "SimDataFormats/Track/interface/SimTrackContainer.h"
#include "SimDataFormats/Vertex/interface/SimVertexContainer.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "TH2D.h"
#include "TH1D.h"
#include "TLorentzVector.h"

#include <ext/algorithm>  // For sorting?

using namespace std;

class DecayAnalyzer : public edm::EDAnalyzer {
   public:
      explicit DecayAnalyzer(const edm::ParameterSet&);
      ~DecayAnalyzer();

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


  edm::InputTag genParticleTag_;

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
DecayAnalyzer::DecayAnalyzer(const edm::ParameterSet& iConfig) {

  genParticleTag_ = iConfig.getParameter<edm::InputTag>("genParticleTag");
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

  h_vxy = fs->make<TH1D>("vxy", "vxy", 100, 0, 5);
  h_vx = fs->make<TH1D>("vx", "vx", 100, 0, 5);
  h_vy = fs->make<TH1D>("vy", "vy", 100, 0, 5);
  h_vz = fs->make<TH1D>("vz", "vz", 100, 0, 50);
  h_vxyz = fs->make<TH2D>("vxyz", ";vz;vxy", 100, 0, 50, 100, 0, 5);
  h_decayVxy = fs->make<TH1D>("decayVxy", "decayVxy", 100, 0, 1000);
  h_decayVx = fs->make<TH1D>("decayVx", "decayVx", 100, 0, 1000);
  h_decayVy = fs->make<TH1D>("decayVy", "decayVy", 100, 0, 1000);
  h_decayVz = fs->make<TH1D>("decayVz", "decayVz", 100, 0, 1500);
  h_decayVxyz = fs->make<TH2D>("decayVxyz", "Position of #chi^{0}#rightarrow#chi^{#pm}#pi^{#mp} decay;|z| [cm];|#rho| [cm]" , 100, 0, 1500, 100, 0, 1000);
  h_decayVxyzWide = fs->make<TH2D>("decayVxyzWide",  "Position of #chi^{0}#rightarrow#chi^{#pm}#pi^{#mp} decay;|z| [cm];|#rho| [cm]" , 1000, 0, 3000, 1000, 0, 2000);
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

}


DecayAnalyzer::~DecayAnalyzer()
{

   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}




//
// member functions
//

// ------------ method called for each event  ------------
void
DecayAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;

   Handle<reco::GenParticleCollection > genParticles;
   iEvent.getByLabel(genParticleTag_,genParticles);
   //   iEvent.getByLabel("genParticles",genParticles);

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
       const reco::Candidate & mcParticle = (*genParticles)[k];

       int status = mcParticle.status();
       int pdgId  = mcParticle.pdgId();
       int numdgt = mcParticle.numberOfDaughters();

       // since status 3 particles are always first in the list,
       // we can break once we're through these
       if(!isParticleGun_ && mcParticle.status() != 3) break;
       
       if(abs(mcParticle.pdgId()) != 1000022) continue; // only consider neutralino parents

       nGenNeutralinoTot++;

       h_GenEta->Fill(mcParticle.eta());
       if (fabs(mcParticle.eta()) > MaxEta_) continue;

       nGenNeutralinoSel++;
       h_GenEtaSel->Fill(mcParticle.eta());

       quiet_ || cout << "Found gen particle with:  "
                      << " pt = " << mcParticle.pt()
                      << ", eta = " << mcParticle.eta()
                      << ", status=" << status
                      << "; pdgId=" << pdgId
                      << "; numdgt=" << numdgt << endl;
       
       TLorentzVector p4(mcParticle.px(),
                         mcParticle.py(),
                         mcParticle.pz(),
                         mcParticle.energy());
       pTot += p4;

       h_neutralinoPt->Fill(p4.Pt());

       vx = mcParticle.vx();
       vy = mcParticle.vy();
       vz = mcParticle.vz();
       vxy = sqrt(vx*vx + vy*vy);

       const reco::Candidate *mother   = &mcParticle;
       const reco::Candidate *daughter = &mcParticle;

       // Descend the decay chain until no daughters have the same PDG ID as mcParticle.
       while (true) {
         bool foundDauSamePdgId = false;
         for (uint i=0; i<daughter->numberOfDaughters(); i++) {
           if (daughter->daughter(i)->pdgId() == mcParticle.pdgId()) {
             foundDauSamePdgId = true;
             mother = daughter;
             daughter = daughter->daughter(i);
             break;
           }
         }
         if (!foundDauSamePdgId) break;
       }
       // Now daughter has no daughters with the same PDG ID as mcParticle.
       // Next choose the daughter with the outermost production vertex, in case there are multiple vertices
       // (e.g., an electron delta ray can produce a vertex before the decay vertex)
       double radiusLastVtx = -99;
       int idxDauLastVtx = -99;
       for (uint i=0; i<daughter->numberOfDaughters(); i++) {
         double testVx = daughter->daughter(i)->vx();
         double testVy = daughter->daughter(i)->vy();
         double testVz = daughter->daughter(i)->vz();
         double radius = sqrt(testVx*testVx + testVy*testVy + testVz*testVz);
         if (radius > radiusLastVtx) {
           radiusLastVtx = radius;
           idxDauLastVtx = i;
         }
       }
       if (idxDauLastVtx>=0) {
         mother = daughter;
         daughter = daughter->daughter(idxDauLastVtx);
       }

       // old:
//        while (daughter->numberOfDaughters () &&          // Find the daughter that is
//               (daughter->status () == 3 ||               // not status 3
//                daughter->pdgId() == mcParticle.pdgId())  // and has different pdgId
//               )  daughter = daughter->daughter (0);


       if (daughter->pdgId() != 1000024 &&
           abs(daughter->pdgId()) != 211) {
         quiet_ || cout << "DebuggingWarning:  will skip event with daughter PDG ID = " << daughter->pdgId()
                        << " and chargino ID = " << mcParticle.pdgId() << endl;
         continue;
       }

       h_GenEtaFoundVtx->Fill(mcParticle.eta());
       nVtxNeutralinoToChargino++;

       double daughter0Id = -99;
       double daughter0E  = -99;
       double daughter1Id = -99;
       double daughter1E  = -99;
       double daughterTotE  = 0;

       decayVx = daughter->vx ();
       decayVy = daughter->vy ();
       decayVz = daughter->vz ();
       decayVxy = sqrt(decayVx*decayVx + decayVy*decayVy);

       TVector3 source (mcParticle.vx (), mcParticle.vy (), mcParticle.vz ()),
         sink (daughter->vx (), daughter->vy (), daughter->vz ());

       quiet_ || cout << "Debug:  BNstop particle pdg id = " << mcParticle.pdgId()
                      << ", status = " << mcParticle.status()
                      << ", pt = " << mcParticle.pt()
                      << ", eta = " << mcParticle.eta()
                      << ", num daughters = " << mcParticle.numberOfDaughters()
                      << endl
                      << "  production:  "
                      << "  vx = " << mcParticle.vx ()
                      << "  vy = " << mcParticle.vy ()
                      << "  vz = " << mcParticle.vz ()
                      << endl
                      << "  decay:  "
                      << "  vx = " << daughter->vx ()
                      << "  vy = " << daughter->vy ()
                      << "  vz = " << daughter->vz ()
                      << endl
                      << "  daughter PDG ID:  " << daughter->pdgId()
                      << ", status = " << daughter->status()
                      << ", pt = " << daughter->pt()
                      << ", eta = " << daughter->eta()
                      << endl
                      << "  mother PDG ID:  " << mother->pdgId()
                      << ", status:  " << mother->status()
                      << ", numdau:  " << mother->numberOfDaughters()
                      << endl;
       for (uint i=0; i<mother->numberOfDaughters(); i++) {
          quiet_ || cout << "  Mother has daughter " << i << ": " << mother->daughter(i)->pdgId() << endl;
       }

       if (source==sink) {
         // Set to non-physical values if no daughters are found
         decayLength  = -99;
         ctau         = -99;
//          betaAtDecay  = -99;
//          gammaAtDecay = -99;
         decayVx  = -99;
         decayVy  = -99;
         decayVz  = -99;
         decayVxy = -99;
       } else {
         decayLength = (sink - source).Mag ();
         ctau = (sink - source).Mag () / (mcParticle.p4 ().Beta () * mcParticle.p4 ().Gamma ());
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

       if (decayLength==0) {
         quiet_ || cout << "Warning:  found particle " << pdgId
                        << " with decayLength = " << decayLength
                        << ", status = " << status
                        << ", numdgt = " << numdgt
                        << endl;
         nNeutralinoNoDecay++;
       }

       if( (mcParticle.numberOfDaughters()>=1) ){
         daughter0Id = mcParticle.daughter(0)->pdgId();
         daughter0E  = mcParticle.daughter(0)->energy();
         daughterTotE += daughter0E;
         quiet_ || cout << "Debug:  daughter0Id = " << mcParticle.daughter(0)->pdgId()
                        << ", status = " << mcParticle.daughter(0)->status()
                        << ", pt = " <<     mcParticle.daughter(0)->pt()
                        << ", eta = " <<    mcParticle.daughter(0)->eta()
                        << ", vx = " << mcParticle.daughter(0)->vx()
                        << ", vy = " << mcParticle.daughter(0)->vy()
                        << ", vz = " << mcParticle.daughter(0)->vz()
                        << endl;


         if( (mcParticle.numberOfDaughters()>=2) ){
           daughter1Id = mcParticle.daughter(1)->pdgId();
           daughter1E  = mcParticle.daughter(1)->energy();
           daughterTotE += daughter1E;

           quiet_ || cout << "Debug:  daughter1Id = " << mcParticle.daughter(1)->pdgId()
                          << ", status = " << mcParticle.daughter(1)->status()
                          << ", pt = " <<        mcParticle.daughter(1)->pt()
                          << ", eta = " <<       mcParticle.daughter(1)->eta()
                          << ", vx = " << mcParticle.daughter(1)->vx()
                          << ", vy = " << mcParticle.daughter(1)->vy()
                          << ", vz = " << mcParticle.daughter(1)->vz()
                          << endl;

         }
       }

       h_daughter0Id->Fill(daughter0Id);
       h_daughter0E ->Fill(daughter0E);
       h_daughter1Id->Fill(fabs(daughter1Id));
       h_daughter1E ->Fill(daughter1E);
       if (decayLength!=0) {
         h_EDiff      ->Fill(daughterTotE - mcParticle.energy());
       }
       quiet_ || cout << "Debug:  daughter0E = " << daughter0E
                      << ", daughter1E = " << daughter1E
                      << ", mother EE = " << mcParticle.energy()
                      << ", mother status = " << mcParticle.status()
                      << ", |daughter0Id| = 1e6 + " << fabs(daughter0Id) - 1e6
                      << ", daughter0 status = " << mcParticle.daughter(0)->status()
                      << ", daughter1Id = " << daughter1Id
                      << ", Ediff = " << daughterTotE - mcParticle.energy()
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
DecayAnalyzer::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void
DecayAnalyzer::endJob()
{
}

// ------------ method called when starting to processes a run  ------------
void
DecayAnalyzer::beginRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a run  ------------
void
DecayAnalyzer::endRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when starting to processes a luminosity block  ------------
void
DecayAnalyzer::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a luminosity block  ------------
void
DecayAnalyzer::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
DecayAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(DecayAnalyzer);
