#ifndef AppTrksEVENT_VARIABLE_PRODUCER
#define AppTrksEVENT_VARIABLE_PRODUCER

#include "OSUT3Analysis/AnaTools/interface/EventVariableProducer.h"
#include "OSUT3Analysis/AnaTools/interface/DataFormat.h"
#include "OSUT3Analysis/AnaTools/interface/ValueLookupTree.h"
struct OriginalCollections
{
  edm::Handle<vector<pat::Electron> >       electrons;
  edm::Handle<vector<pat::Jet> >            jets;
  edm::Handle<vector<pat::Muon> >           muons;
  edm::Handle<vector<reco::Vertex> >        primaryvertexs;
  edm::Handle<vector<PileupSummaryInfo>>    pileupinfos;
  edm::Handle<edm::TriggerResults>          triggers;
};

class AppTrksEventVariableProducer : public EventVariableProducer
  {
    public:
	AppTrksEventVariableProducer (const edm::ParameterSet &);
        void getOriginalCollections (const unordered_set<string> &objectsToGet, const edm::ParameterSet &collections, OriginalCollections &handles, const edm::Event &event);
        bool passCleaning (double eta, double phi, OriginalCollections &handles);
        ~AppTrksEventVariableProducer ();
        OriginalCollections handles_;

    private:
	void AddVariables(const edm::Event &);
        string type_;
        string triggerPath_;
        double triggerScalingFactor_;
  };
#endif
