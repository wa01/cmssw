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

#include "DataFormats/GsfTrackReco/interface/GsfTrack.h"
#include "DataFormats/GsfTrackReco/interface/GsfTrackFwd.h"

#include "Geometry/Records/interface/GlobalTrackingGeometryRecord.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"

#include "Geometry/CommonDetUnit/interface/GlobalTrackingGeometry.h"
#include "MagneticField/Engine/interface/MagneticField.h"

#include "TrackingTools/TrajectoryState/interface/TrajectoryStateOnSurface.h"
#include "TrackingTools/GsfTools/interface/MultiTrajectoryStateTransform.h"
#include "TrackingTools/GsfTools/interface/MultiGaussianStateTransform.h"
#include "TrackingTools/GsfTools/interface/MultiGaussianState1D.h"
#include "TrackingTools/GsfTools/interface/GaussianSumUtilities1D.h"


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
  edm::EDGetTokenT<edm::View<reco::GsfTrack> > gsfTrackToken_;
  edm::ESHandle<GlobalTrackingGeometry> trackingGeometryHandle_;
  edm::ESHandle<MagneticField> magneticFieldHandle_;
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

   gsfTrackToken_ = consumes<edm::View<reco::GsfTrack> >(edm::InputTag("electronGsfTracks"));
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

   edm::Handle<edm::View<reco::GsfTrack> > gsfTrackHandle;
   iEvent.getByToken(gsfTrackToken_, gsfTrackHandle);

   for ( edm::View<reco::GsfTrack>::const_iterator igsf=gsfTrackHandle->begin();
	 igsf!=gsfTrackHandle->end(); ++igsf ) {
     TrajectoryStateOnSurface innTSOS = mtst.innerStateOnSurface(*igsf);
     TrajectoryStateOnSurface vtxTSOS = mtst.extrapolatedState(innTSOS,GlobalPoint(0.,0.,0.));
     std::cout << "Gsf pt = " << igsf->pt() 
	       << " innTSOS " << " " << innTSOS.components().size()
	       << " vtxTSOS " << " " << vtxTSOS.components().size() << std::endl;
     // LocalPoint lp(vtxTSOS.surface().toLocal(GlobalPoint(0.,0.,0.)));
     // std::cout << "   " << lp.x() << " " << lp.y() << " " << lp.z() << std::endl;
     // GlobalVector gy(vtxTSOS.surface().toGlobal(LocalVector(0.,1.,0.)));
     // std::cout << "   " << gy.x() << " " << gy.y() << " " << gy.z() << std::endl;

     AlgebraicVector5 modeParameters;
     AlgebraicSymMatrix55 modeCovariance;
     // set parameters and variances for "mode" state (local parameters)
     for ( unsigned int iv=0; iv<5; ++iv ) {
       MultiGaussianState1D state1D = MultiGaussianStateTransform::multiState1D(vtxTSOS,iv);
       GaussianSumUtilities1D utils(state1D);
       modeParameters(iv) = utils.mode().mean();
       modeCovariance(iv,iv) = utils.mode().variance();
       if ( !utils.modeIsValid() ) {
	 // if mode calculation fails: use mean
	 modeParameters(iv) = utils.mean();
	 modeCovariance(iv,iv) = utils.variance();
       }
     }

   }
}


// ------------ method called once each job just before starting event loop  ------------
void 
GsfTrackAnalyzer::beginJob()
{
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
