#include "DisplacedMuonProducer.h"

#if IS_VALID(muons)
#include "OSUT3Analysis/AnaTools/interface/CommonUtils.h"

DisplacedMuonProducer::DisplacedMuonProducer (const edm::ParameterSet &cfg) :
  muons_ (cfg.getParameter<edm::InputTag> ("muons")),
  cfg_ (cfg)
{
  produces<vector<osu::Muon> > ();
}

DisplacedMuonProducer::~DisplacedMuonProducer ()
{
}

void
DisplacedMuonProducer::produce (edm::Event &event, const edm::EventSetup &setup)
{
  edm::Handle<vector<osu::Muon>> collection;
  if (!anatools::getCollection (muons_, collection, event, false))
    return;
  edm::Handle<vector<osu::Mcparticle> > particles;
  anatools::getCollection (edm::InputTag ("", ""), particles, event);

  pl_ = auto_ptr<vector<osu::Muon> > (new vector<osu::Muon> ());
  for (const auto &object : *collection)
    {
      osu::Muon muon1 (object, particles, cfg_);
      for (const auto &object : *collection)
	{
	  osu::Muon muon2 (object, particles, cfg_);
	  if(durp


	  if(muon1.isGlobalMuon() && muon2.isGlobalMuon() && PassCosmicSelection(muon1,muon2) && muon1.phi() != muon2.phi() && muon1.pt() > 30)
	    {
	      pl_->push_back (muon1);
	      break;
	    }
	}
    }

  event.put (pl_);
  pl_.reset ();
}

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(DisplacedMuonProducer);

#endif
