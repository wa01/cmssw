// -*- C++ -*-
//
// Package:    Workspace/GsfTrackAnalyzer
// Class:      GsfTrackAnalyzer
// 
/**\class GsfTrackAnalyzer GsfTrackAnalyzer.cc Workspace/GsfTrackAnalyzer/plugins/GsfTrackAnalyzer.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Wolfgang Adam
//         Created:  Fri, 15 Sep 2017 10:06:29 GMT
//
//


// system include files
#include <memory>
#include <iostream>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/ESHandle.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectronFwd.h"

#include "DataFormats/GsfTrackReco/interface/GsfTrack.h"
#include "DataFormats/GsfTrackReco/interface/GsfTrackFwd.h"

#include "Geometry/Records/interface/GlobalTrackingGeometryRecord.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"

#include "Geometry/CommonDetUnit/interface/GlobalTrackingGeometry.h"
#include "MagneticField/Engine/interface/MagneticField.h"

#include "TrackingTools/TrajectoryState/interface/TrajectoryStateOnSurface.h"
#include "TrackingTools/GsfTools/interface/MultiTrajectoryStateTransform.h"
#include "TrackingTools/GsfTools/interface/MultiTrajectoryStateMode.h"
#include "TrackingTools/GsfTools/interface/MultiGaussianStateTransform.h"
#include "TrackingTools/GsfTools/interface/MultiGaussianState1D.h"
#include "TrackingTools/GsfTools/interface/GaussianSumUtilities1D.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "TTree.h"
#include "TVector3.h"

#include <map>
#include <vector>

//
// class declaration
//

// If the analyzer does not use TFileService, please remove
// the template argument to the base class so the class inherits
// from  edm::one::EDAnalyzer<> and also remove the line from
// constructor "usesResource("TFileService");"
// This will improve performance in multithreaded jobs.

class GsfTrackAnalyzer : public edm::one::EDAnalyzer<edm::one::SharedResources>  {
   public:
      explicit GsfTrackAnalyzer(const edm::ParameterSet&);
      ~GsfTrackAnalyzer();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


   private:
      virtual void beginJob() override;
      virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
      virtual void endJob() override;

      // ----------member data ---------------------------
private:
  edm::EDGetTokenT<reco::GsfElectronCollection> gsfElectronToken_;
  edm::EDGetTokenT<reco::GsfTrackCollection> gsfTrackToken_;
  edm::ESHandle<GlobalTrackingGeometry> trackingGeometryHandle_;
  edm::ESHandle<MagneticField> magneticFieldHandle_;

  MultiTrajectoryStateMode mtsMode_;

  TTree* tree_;
  int run_;
  int lumi_;
  long long evt_;
  short int ics_[2];
  short int isEle_;
  float scPars_[3];
  TVector3 eleMom_;
  TVector3 gsfMomMean_;
  TVector3 gsfMomModeCart_;
  TVector3 gsfMomModeLoc_;
  
  double wgt_;
  AlgebraicVector5 locPars_;
  AlgebraicSymMatrix55 locCov_;
  double gsfMomCart_[3];
  double gsfErrMomCart_[3];
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
GsfTrackAnalyzer::GsfTrackAnalyzer(const edm::ParameterSet& iConfig)

{
   //now do what ever initialization is needed
   usesResource("TFileService");

   gsfElectronToken_ = consumes< reco::GsfElectronCollection >(edm::InputTag("gedGsfElectrons"));
   gsfTrackToken_ = consumes<reco::GsfTrackCollection >(edm::InputTag("electronGsfTracks"));
}


GsfTrackAnalyzer::~GsfTrackAnalyzer()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
GsfTrackAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;

   iSetup.get<GlobalTrackingGeometryRecord>().get(trackingGeometryHandle_); 
   iSetup.get<IdealMagneticFieldRecord>().get(magneticFieldHandle_);
   MultiTrajectoryStateTransform mtst(&*trackingGeometryHandle_,&*magneticFieldHandle_);

   edm::Handle<reco::GsfElectronCollection> gsfElectronHandle;
   iEvent.getByToken(gsfElectronToken_, gsfElectronHandle);

   // std::vector<reco::GsfTrackRef> electronTrackRefs;
   std::map<const reco::GsfTrackRef, const reco::GsfElectron*> electronTrackRefs;
   for ( reco::GsfElectronCollection::const_iterator iele=gsfElectronHandle->begin();
	 iele!=gsfElectronHandle->end(); ++iele ) {
     if ( !(iele->gsfTrack().isNull()) ) 
       electronTrackRefs[iele->gsfTrack()] = &*iele;
       // electronTrackRefs.push_back(iele->gsfTrack());
   }

   edm::Handle<reco::GsfTrackCollection > gsfTrackHandle;
   iEvent.getByToken(gsfTrackToken_, gsfTrackHandle);

   // for ( reco::GsfTrackCollection::const_iterator igsf=gsfTrackHandle->begin();
	 // igsf!=gsfTrackHandle->end(); ++igsf ) {
   for ( size_t jgsf=0; jgsf<gsfTrackHandle->size(); ++jgsf ) {
     const reco::GsfTrack* igsf(&(*gsfTrackHandle)[jgsf]);
     reco::GsfTrackRef gsfRef(gsfTrackHandle,jgsf);
     const reco::GsfElectron* electron(0);
     std::map<const reco::GsfTrackRef, const reco::GsfElectron*>::const_iterator imap = 
       electronTrackRefs.find(gsfRef);
       // std::find(electronTrackRefs.begin(),electronTrackRefs.end(),gsfRef);
     if ( imap!=electronTrackRefs.end() ) {
       isEle_ = 1;
       electron = (*imap).second;
       reco::SuperClusterRef scRef = electron->superCluster();
       scPars_[0] = scRef->correctedEnergy();
       scPars_[1] = scRef->eta();
       scPars_[2] = scRef->phi();
       eleMom_.SetPtEtaPhi(electron->pt(),electron->eta(),electron->phi());
     }
     else {
       isEle_ = 0;
       scPars_[0] = 0.;
       scPars_[1] = 0.;
       scPars_[2] = 0.;
       eleMom_.SetXYZ(0.,0.,0.);
     }

     TrajectoryStateOnSurface innTSOS = mtst.innerStateOnSurface(*igsf);
     TrajectoryStateOnSurface vtxTSOS = mtst.extrapolatedState(innTSOS,GlobalPoint(0.,0.,0.));

     ics_[0] = vtxTSOS.components().size();
     ics_[1] = -1;

     GlobalVector mom(vtxTSOS.globalMomentum());
     gsfMomMean_.SetXYZ(mom.x(),mom.y(),mom.z());
     mtsMode_.momentumFromModeCartesian(vtxTSOS,mom);
     gsfMomModeCart_.SetXYZ(mom.x(),mom.y(),mom.z());
     mtsMode_.momentumFromModeLocal(vtxTSOS,mom);
     gsfMomModeLoc_.SetXYZ(mom.x(),mom.y(),mom.z());

     wgt_ = vtxTSOS.weight();
     locPars_ = vtxTSOS.localParameters().vector();
     locCov_ = vtxTSOS.localError().matrix();

     for ( size_t i=0; i<3; ++i ) {
       gsfMomCart_[i] = 0.;
       gsfErrMomCart_[i] = 0.;
     }

     tree_->Fill();

     for ( size_t ic=0; ic<vtxTSOS.components().size(); ++ic ) {
       TrajectoryStateOnSurface tsos = vtxTSOS.components()[ic];
       
       ics_[1] = ic;
       wgt_ = tsos.weight();
       locPars_ = tsos.localParameters().vector();
       locCov_ = tsos.localError().matrix();

       mom = tsos.globalMomentum();
       const AlgebraicSymMatrix66& errCart = tsos.cartesianError().matrix();
       gsfMomCart_[0] = mom.x();
       gsfMomCart_[1] = mom.y();
       gsfMomCart_[2] = mom.z();
       for ( size_t i=0; i<3; ++i )  gsfErrMomCart_[i] = sqrt(errCart(i+3,i+3));

       tree_->Fill();
     }
     // // std::cout << "Gsf pt = " << igsf->pt() 
     // // 	       << " innTSOS " << " " << innTSOS.components().size()
     // // 	       << " vtxTSOS " << " " << vtxTSOS.components().size() << std::endl;
     // // // LocalPoint lp(vtxTSOS.surface().toLocal(GlobalPoint(0.,0.,0.)));
     // // // std::cout << "   " << lp.x() << " " << lp.y() << " " << lp.z() << std::endl;
     // // // GlobalVector gy(vtxTSOS.surface().toGlobal(LocalVector(0.,1.,0.)));
     // // // std::cout << "   " << gy.x() << " " << gy.y() << " " << gy.z() << std::endl;

     // AlgebraicVector5 modeParameters;
     // AlgebraicSymMatrix55 modeCovariance;
     // // set parameters and variances for "mode" state (local parameters)
     // for ( unsigned int iv=0; iv<5; ++iv ) {
     //   MultiGaussianState1D state1D = MultiGaussianStateTransform::multiState1D(vtxTSOS,iv);
     //   GaussianSumUtilities1D utils(state1D);
     //   modeParameters(iv) = utils.mode().mean();
     //   modeCovariance(iv,iv) = utils.mode().variance();
     //   if ( !utils.modeIsValid() ) {
     // 	 // if mode calculation fails: use mean
     // 	 modeParameters(iv) = utils.mean();
     // 	 modeCovariance(iv,iv) = utils.variance();
     //   }
     // }

   }
}


// ------------ method called once each job just before starting event loop  ------------
void 
GsfTrackAnalyzer::beginJob()
{
  edm::Service<TFileService> fs;

  tree_ = fs->make<TTree>("GsfTree","GsfTree");
  tree_->Branch("run",&run_,"run/I");
  tree_->Branch("lumi",&lumi_,"lumi/I");
  tree_->Branch("evt",&evt_,"evt/L");
  tree_->Branch("isEle",&isEle_,"isEle/S");
  tree_->Branch("scPars",scPars_,"scE/F:scEta/F:scPhi/F");
  tree_->Branch("eleMom",&eleMom_);
  tree_->Branch("gsfMomMean",&gsfMomMean_);
  tree_->Branch("gsfMomModeCart",&gsfMomModeCart_);
  tree_->Branch("gsfMomModeLoc",&gsfMomModeLoc_);
  tree_->Branch("ics",ics_,"nc/S:ic/S");
  tree_->Branch("wgt",&wgt_,"wgt/D");
  tree_->Branch("locPars",&locPars_);
  tree_->Branch("locCov",&locCov_);
  tree_->Branch("gsfMomCart",gsfMomCart_,"gsfPx/D:gsfPy/D:gsfPz/D");
  tree_->Branch("gsfErrMomCart",gsfErrMomCart_,"gsfSigPx/D:gsfSigPy/D:gsfSigPz/D");
  
  
}

// ------------ method called once each job just after ending the event loop  ------------
void 
GsfTrackAnalyzer::endJob() 
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
GsfTrackAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(GsfTrackAnalyzer);
