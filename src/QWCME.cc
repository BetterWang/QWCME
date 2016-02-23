// -*- C++ -*-
//
// Package:    QWCME
// Class:      QWCME
// 
/**\class QWCME QWCME.cc QWAna/QWCME/src/QWCME.cc

*/
//
// Original Author:  Quan Wang
//         Created:  02/23/2016
//
//


// system include files
#include <memory>
#include <algorithm>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "DataFormats/Common/interface/TriggerResults.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include <DataFormats/ParticleFlowCandidate/interface/PFCandidateFwd.h>
#include <DataFormats/ParticleFlowCandidate/interface/PFCandidate.h>
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/HeavyIonEvent/interface/EvtPlane.h"
#include "TH1.h"
#include "TH2.h"
#include "TNtuple.h"
#include "TComplex.h"
#include <complex>


#include "QWAna/QWCME/interface/QWCME.h"


using namespace std;

//#ifdef QW_DEBUG
//
// constructors and destructor
//
QWCME::QWCME(const edm::ParameterSet& iConfig)
	:
		trackTag_( iConfig.getUntrackedParameter<edm::InputTag>("tracks_") )
	,	centralityToken_( consumes<int>(iConfig.getParameter<edm::InputTag>("centrality_")) )
	,	vertexToken_( consumes<reco::VertexCollection>(iConfig.getUntrackedParameter<edm::InputTag>("vertexSrc_")) )
	,	algoParameters_(iConfig.getParameter<std::vector<int> >("algoParameters"))
	 ,	hpp_2p({1, 1, -2})
	 ,	hnn_2n({1, 1, -2})
	 ,	hpp({1, 1})
	 ,	hnn({1, 1})
	 ,	h_2p({-2})
	 ,	h_2n({-2})
	 ,	hp_2p({1, -2})
	 ,	hn_2n({1, -2})
{
	//now do what ever initialization is needed
	minvz_ = iConfig.getUntrackedParameter<double>("minvz_", -15.);
	maxvz_ = iConfig.getUntrackedParameter<double>("maxvz_", 15.);
	dzdzerror_ = iConfig.getUntrackedParameter<double>("dzdzerror_", 3.);
	d0d0error_ = iConfig.getUntrackedParameter<double>("d0d0error_", 3.);
	chi2_ = iConfig.getUntrackedParameter<double>("chi2_", 40);
	pterrorpt_ = iConfig.getUntrackedParameter<double>("pterrorpt_", 0.1);

	fweight_ = iConfig.getUntrackedParameter<edm::InputTag>("fweight_", string("NA"));
	facceptance_ = iConfig.getUntrackedParameter<edm::InputTag>("facceptance_", string("NA"));
	bEff = iConfig.getUntrackedParameter<bool>("bEff_", false);
	bPhiEta = iConfig.getUntrackedParameter<bool>("bPhiEta_", false);
	bCentNoff = iConfig.getUntrackedParameter<bool>("bCentNoff_", false);
	bCaloMatching_ = iConfig.getUntrackedParameter<bool>("bCaloMaching", false);
	Noffmin_ = iConfig.getUntrackedParameter<int>("Noffmin_", 0);
	Noffmax_ = iConfig.getUntrackedParameter<int>("Noffmax_", 5000);
	cmode_ = iConfig.getUntrackedParameter<int>("cmode_", 1);
	bGen_ = iConfig.getUntrackedParameter<bool>("bGen_", false);
	nvtx_ = iConfig.getUntrackedParameter<int>("nvtx_", 100);
	bFlipEta_ = iConfig.getUntrackedParameter<bool>("bFlipEta_", false);
	reso_ = iConfig.getUntrackedParameter<double>("reso", 0.2);

	if ( bGen_ ) {
		trackGenToken_ = consumes<reco::GenParticle>(trackTag_);
	} else {
		trackToken_ = consumes<reco::TrackCollection>(trackTag_);
	}

	string streff = fweight_.label();
	if ( streff == string("NA") ) {
		cout << "!!! eff NA" << endl;
		bEff = false;
		fEffFak = 0;
	} else {
		fEffFak = new TFile(streff.c_str());
		if ( !fEffFak->IsOpen() ) {
			bEff = false;
		} else {
			cout << "!!! Using particle weight " << streff << endl;
			if ( bEff ) {
				cout << "!!! Apply Eff correction" << endl;
				for ( int i = 0; i < 20; i++ ) {
					if ( streff == string("PbPb_MB_TT_5TeV_v2.root") or streff == string("PbPb_dijet_TT_5TeV_v2.root") ) {
						TH2D * h = (TH2D*) fEffFak->Get("rTotalEff3D_0_5");
						for ( int c = 0; c < 10; c++ ) {
							hEff_cbin[c] = h;
						}
						h = (TH2D*) fEffFak->Get("rTotalEff3D_5_10");
						for ( int c = 10; c < 20; c++ ) {
							hEff_cbin[c] = h;
						}
						h = (TH2D*) fEffFak->Get("rTotalEff3D_10_30");
						for ( int c = 20; c < 60; c++ ) {
							hEff_cbin[c] = h;
						}
						h = (TH2D*) fEffFak->Get("rTotalEff3D_30_50");
						for ( int c = 60; c < 100; c++ ) {
							hEff_cbin[c] = h;
						}
						h = (TH2D*) fEffFak->Get("rTotalEff3D_50_100");
						for ( int c = 100; c < 200; c++ ) {
							hEff_cbin[c] = h;
						}
					}
				}
				cout << "!!! eff histo done" << endl;
			}
		}
	}
	string stracc = facceptance_.label();
	if ( stracc == string("NA") ) {
		cout << "!!! acc NA" << endl;
		bacc = false;
		facc = 0;
	} else {
		facc = new TFile(stracc.c_str());
		if ( !facc->IsOpen() ) {
			bacc = false;
		} else {
			cout << "!!! Using acceptance weight " << stracc << endl;
			bacc = true;
			for ( int cent = 0; cent < nCentBins; cent++ ) {
				for ( int ipt = 0; ipt < nPtBins; ipt++ ) {
					hacc[cent][ipt][0] = (TH2D*) facc->Get(Form("hPhiEta_%i_%i_0", cent, ipt));
					hacc[cent][ipt][1] = (TH2D*) facc->Get(Form("hPhiEta_%i_%i_1", cent, ipt));
				}
			}
		}
	}


	t = new QWEvent;
	memset(t, 0, sizeof(QWEvent));
	//
	edm::Service<TFileService> fs;
	for ( int cent = 0; cent < nCentBins; cent++ ) {
		hPt[cent] = fs->make<TH1D>(Form("hPt_%i", cent), "", nPtBins, ptbins);
		if ( bPhiEta ) {
			for ( int i = 0; i < nPtBins; i++ ) {
				//cout << "!! new histo cent = " << cent << " of " << nCentBins << "\t ipt = " << i  << " of " << nPtBins << endl;
				hPhiEta[cent][i][0] = fs->make<TH2D>(Form("hPhiEta_%i_%i_0", cent, i), "", 64, -Pi, Pi, 48, -2.4, 2.4);
				hPhiEta[cent][i][1] = fs->make<TH2D>(Form("hPhiEta_%i_%i_1", cent, i), "", 64, -Pi, Pi, 48, -2.4, 2.4);
			}
		}
	}

	for ( int cent = 0; cent < nCentBins; cent++ ) {
		//cout << "!! new histo cent = " << cent << endl;
		hdNdPtdEta[cent] = fs->make<TH2D>(Form("hdNdPtdEta_%i", cent), Form("hdNdPtdEta_%i", cent), nEtaBins, etabins, nPtBins, ptbins);
		hdNdPtdEtaPt[cent] = fs->make<TH2D>(Form("hdNdPtdEtaPt_%i", cent), Form("hdNdPtdEta_%i", cent), nEtaBins, etabins, nPtBins, ptbins);
	}


	if ( bCaloMatching_ ) {
		pfToken_ = consumes<reco::PFCandidateCollection>(iConfig.getUntrackedParameter<edm::InputTag>("pfTag"));
	}

	qpp_2p.resize(hpp_2p);
	qnn_2n.resize(hnn_2n);
	qpp.resize(hpp);
	qnn.resize(hnn);
	q_2p.resize(h_2p);
	q_2n.resize(h_2n);
	qp_2p.resize(hp_2p);
	qn_2n.resize(hn_2n);

	switch ( cmode_ ) {
		case 1:
			cqpp_2p = new correlations::closed::FromQVector( qpp_2p );
			cqnn_2n = new correlations::closed::FromQVector( qnn_2n );
			cqpp    = new correlations::closed::FromQVector( qpp    );
			cqnn    = new correlations::closed::FromQVector( qnn    );
			cq_2p   = new correlations::closed::FromQVector( q_2p   );
			cq_2n   = new correlations::closed::FromQVector( q_2n   );
			cqp_2p  = new correlations::closed::FromQVector( qp_2p  );
			cqn_2n  = new correlations::closed::FromQVector( qn_2n  );
			break;
		case 2:
			cqpp_2p = new correlations::recurrence::FromQVector( qpp_2p );
			cqnn_2n = new correlations::recurrence::FromQVector( qnn_2n );
			cqpp    = new correlations::recurrence::FromQVector( qpp    );
			cqnn    = new correlations::recurrence::FromQVector( qnn    );
			cq_2p   = new correlations::recurrence::FromQVector( q_2p   );
			cq_2n   = new correlations::recurrence::FromQVector( q_2n   );
			cqp_2p  = new correlations::recurrence::FromQVector( qp_2p  );
			cqn_2n  = new correlations::recurrence::FromQVector( qn_2n  );
			break;
		case 3:
			cqpp_2p = new correlations::recursive::FromQVector( qpp_2p );
			cqnn_2n = new correlations::recursive::FromQVector( qnn_2n );
			cqpp    = new correlations::recursive::FromQVector( qpp    );
			cqnn    = new correlations::recursive::FromQVector( qnn    );
			cq_2p   = new correlations::recursive::FromQVector( q_2p   );
			cq_2n   = new correlations::recursive::FromQVector( q_2n   );
			cqp_2p  = new correlations::recursive::FromQVector( qp_2p  );
			cqn_2n  = new correlations::recursive::FromQVector( qn_2n  );
			break;
	}

	trV = fs->make<TTree>("trV", "trV");
	trV->Branch("Noff", &gNoff, "Noff/I");
	trV->Branch("Mult", &gMult, "Mult/I");

	trV->Branch("dpp_2p", &dpp_2p, "dpp_2p/D");
	trV->Branch("dnn_2n", &dnn_2n, "dnn_2n/D");
	trV->Branch("dpp",    &dpp   , "dpp/D");
	trV->Branch("dnn",    &dnn   , "dnn/D");
	trV->Branch("d_2p",   &d_2p  , "d_2p/D");
	trV->Branch("d_2n",   &d_2n  , "d_2n/D");
	trV->Branch("dp_2p",  &dp_2p , "dp_2p/D");
	trV->Branch("dn_2n",  &dn_2n , "dn_2n/D");
}

bool
QWCME::CaloMatch(const reco::Track & track, const edm::Event & iEvent, unsigned int idx)
{
	if ( !bCaloMatching_ ) return true;
	edm::Handle<reco::PFCandidateCollection> pfCand;
	iEvent.getByToken( pfToken_, pfCand );
	double energy = 0;
	for ( reco::PFCandidateCollection::const_iterator it = pfCand->begin();
			it != pfCand->end();
			++it ) {
		if ( (it->particleId() != reco::PFCandidate::h) and
				(it->particleId() != reco::PFCandidate::e) and
				(it->particleId() != reco::PFCandidate::mu) ) continue;
		if ( idx == it->trackRef().key() ) {
			energy = it->ecalEnergy() + it->hcalEnergy();
			break;
		}
	}

	if( track.pt() < 20 || ( energy/( track.pt()*TMath::CosH(track.eta() ) ) > reso_ && (energy)/(TMath::CosH(track.eta())) > (track.pt() - 80.0) )  ) return true;
	else return false;
}

QWCME::~QWCME()
{

	// do anything here that needs to be done at desctruction time
	// (e.g. close files, deallocate resources etc.)

}

int
QWCumuV3::getNoffCent(const edm::Event& iEvent, const edm::EventSetup& iSetup, int& Noff)
{
        // very hard coded Noff track centrality cut
        using namespace edm;
        using namespace reco;
//      int Noff = 0;

        Handle<VertexCollection> vertexCollection;
        iEvent.getByToken(vertexToken_, vertexCollection);
        const VertexCollection * recoVertices = vertexCollection.product();

        int primaryvtx = 0;
        math::XYZPoint v1( (*recoVertices)[primaryvtx].position().x(), (*recoVertices)[primaryvtx].position().y(), (*recoVertices)[primaryvtx].position().z() );
        double vxError = (*recoVertices)[primaryvtx].xError();
        double vyError = (*recoVertices)[primaryvtx].yError();
        double vzError = (*recoVertices)[primaryvtx].zError();


        Handle<TrackCollection> tracks;
        iEvent.getByToken(trackToken_,tracks);
        for(TrackCollection::const_iterator itTrack = tracks->begin();
                itTrack != tracks->end();
                ++itTrack) {

                if ( !itTrack->quality(reco::TrackBase::highPurity) ) continue;
                if ( itTrack->charge() == 0 ) continue;
                if ( itTrack->pt() < 0.4 ) continue;

                double d0 = -1.* itTrack->dxy(v1);
                double derror=sqrt(itTrack->dxyError()*itTrack->dxyError()+vxError*vyError);
                double dz=itTrack->dz(v1);
                double dzerror=sqrt(itTrack->dzError()*itTrack->dzError()+vzError*vzError);
                if ( fabs(itTrack->eta()) > 2.4 ) continue;
                if ( fabs( dz/dzerror ) > 3. ) continue;
                if ( fabs( d0/derror ) > 3. ) continue;
                if ( itTrack->ptError()/itTrack->pt() > 0.1 ) continue;
                if ( itTrack->numberOfValidHits() < 11 ) continue;
                Noff++;
        }

        int cent = nCentNoff-1;
        while ( CentNoffCut[cent] <= Noff ) cent--;
        return cent;
}

void
QWCME::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
	if ( bGen_ ) analyzeGen(iEvent, iSetup);
	else analyzeData(iEvent, iSetup);
	if ( t->Mult == 0 ) return;

	for ( int i = 0; i < t->Mult; i++ ) {
		if ( t->Charge[i] > 0 ) {
			qpp_2p.fill(t->Phi[i], t->weight[i]);
			qpp.fill(t->Phi[i], t->weight[i]);
			q_2p.fill(t->Phi[i], t->weight[i]);
			qp_2p.fill(t->Phi[i], t->weight[i]);
		} else {
			qnn_2n.fill(t->Phi[i], t->weight[i]);
			qnn.fill(t->Phi[i], t->weight[i]);
			q_2n.fill(t->Phi[i], t->weight[i]);
			qn_2n.fill(t->Phi[i], t->weight[i]);
		}
	}

	dpp_2p = cqpp_2p->calculate( 3 , hpp_2p);
	dnn_2n = cqnn_2n->calculate( 3 , hnn_2n);
	dpp    = cqpp   ->calculate( 2 , hpp   );
	dnn    = cqnn   ->calculate( 2 , hnn   );
	d_2p   = cq_2p  ->calculate( 1 , h_2p  );
	d_2n   = cq_2n  ->calculate( 1 , h_2n  );
	dp_2p  = cqp_2p ->calculate( 2 , hp_2p );
	dn_2n  = cqn_2n ->calculate( 2 , hn_2n );

	trV->Fill();

	qpp_2p.reset();
	qnn_2n.reset();
	qpp   .reset();
	qnn   .reset();
	q_2p  .reset();
	q_2n  .reset();
	qp_2p .reset();
	qn_2n .reset();
}


void
QWCME::analyzeGen(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
}


// ------------ method called for each event  ------------
	void
QWCME::analyzeData(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
	//cout << "==> do event" << endl;
	using namespace edm;
	using namespace reco;

	t->Mult = 0;
	// vertex
	Handle<VertexCollection> vertexCollection;
	iEvent.getByToken(vertexToken_, vertexCollection);
	VertexCollection recoVertices = *vertexCollection;
	if ( recoVertices.size() > nvtx_ ) return;
	sort(recoVertices.begin(), recoVertices.end(), [](const reco::Vertex &a, const reco::Vertex &b){
			if ( a.tracksSize() == b.tracksSize() ) return a.chi2() < b.chi2() ? true:false;
			return a.tracksSize() > b.tracksSize() ? true:false;
			});

	int primaryvtx = 0;
	math::XYZPoint v1( recoVertices[primaryvtx].position().x(), recoVertices[primaryvtx].position().y(), recoVertices[primaryvtx].position().z() );
	double vxError = recoVertices[primaryvtx].xError();
	double vyError = recoVertices[primaryvtx].yError();
	double vzError = recoVertices[primaryvtx].zError();

	double vz = recoVertices[primaryvtx].z();
	if (fabs(vz) < minvz_ || fabs(vz) > maxvz_) {
		//cout << __LINE__ << " vz = " << vz << endl;
		return;
	}

	// centrality
	int bin = 0;
	int cbin = 0;
	t->Noff = 0;

	if ( bCentNoff ) {
		cbin = getNoffCent( iEvent, iSetup, t->Noff);
		if ( (t->Noff < Noffmin_) or (t->Noff >= Noffmax_) ) {
			return;
		}
	} else {
		edm::Handle<int> ch;
		iEvent.getByToken(centralityToken_,ch);
		bin = *(ch.product());
		if ( bin < 0 or bin >= 200 ) {
			//cout << __LINE__ << " bin = " << bin << endl;
			return;
		}
		while ( centbins[cbin+1] < bin*.5+0.1 ) cbin++;
		t->Noff = bin;
	}
	bin = cbin;

	// track
	Handle<TrackCollection> tracks;
	iEvent.getByToken(trackToken_,tracks);
	t->Cent = bin;
	t->vz = vz;

	for(TrackCollection::const_iterator itTrack = tracks->begin();
			itTrack != tracks->end();
			++itTrack) {
		if ( itTrack->charge() == 0 ) continue;
		if ( !itTrack->quality(reco::TrackBase::highPurity) ) continue;

		double d0 = -1.* itTrack->dxy(v1);
		double derror=sqrt(itTrack->dxyError()*itTrack->dxyError()+vxError*vyError);
		double dz=itTrack->dz(v1);
		double dzerror=sqrt(itTrack->dzError()*itTrack->dzError()+vzError*vzError);

		if ( fabs(itTrack->eta()) > 2.4 ) continue;
		if ( fabs( dz/dzerror ) > dzdzerror_ ) continue;
		if ( fabs( d0/derror ) > d0d0error_ ) continue;
		if ( itTrack->ptError()/itTrack->pt() > pterrorpt_ ) continue;
		if ( itTrack->numberOfValidHits() < 11 ) continue;
		if ( itTrack->normalizedChi2() / itTrack->hitPattern().trackerLayersWithMeasurement() > 0.15 ) continue;
		if ( find( algoParameters_.begin(), algoParameters_.end(), itTrack->originalAlgo() ) == algoParameters_.end() ) continue;
		if ( !CaloMatch(*itTrack, iEvent, itTrack - tracks->begin()) ) continue;

		t->RFP[t->Mult] = 1;
		t->Charge[t->Mult] = itTrack->charge();
		if ( (charge_ == 1) && (t->Charge[t->Mult]<0) ) {
			t->RFP[t->Mult] = 0;
		}
		if ( (charge_ == -1) && (t->Charge[t->Mult]>0) ) {
			t->RFP[t->Mult] = 0;
		}

		t->Pt[t->Mult] = itTrack->pt();
		if ( t->Pt[t->Mult] >= ptbins[nPtBins] || t->Pt[t->Mult] <= ptbins[0] ) {
			t->RFP[t->Mult] = 0;
		}
		t->Eta[t->Mult] = itTrack->eta();
		if (bFlipEta_) t->Eta[t->Mult] = - t->Eta[t->Mult];

		if ( bEff ) {
			t->rEff[t->Mult] = hEff_cbin[t->Noff]->GetBinContent( hEff_cbin[t->Noff]->FindBin(t->Eta[t->Mult], t->Pt[t->Mult] ) );
		} else {
			t->rEff[t->Mult] = 1.;
		}

		if ( t->rEff[t->Mult] <= 0.1 or TMath::IsNaN(t->rEff[t->Mult]) ) {
			t->RFP[t->Mult] = 0;
		}
		double weight = 1./t->rEff[t->Mult];

		double phi = itTrack->phi();

		double wacc = 1.;
		int ipt=0;

		while ( t->Pt[t->Mult] > ptbins[ipt+1] ) ipt++;
		if ( bacc ) {
			wacc = 1./hacc[bin][ipt][t->Charge[t->Mult]>0]->GetBinContent(hacc[bin][ipt][t->Charge[t->Mult]>0]->FindBin(phi, t->Eta[t->Mult]));
		}
		if ( bPhiEta ) hPhiEta[bin][ipt][t->Charge[t->Mult]>0]->Fill(phi, t->Eta[t->Mult], wacc);

		weight *= wacc;

		if ( (t->Pt[t->Mult] < rfpptmin_) || (t->Pt[t->Mult] > rfpptmax_) || t->Eta[t->Mult] < rfpmineta_ || t->Eta[t->Mult] > rfpmaxeta_ ) {
			t->RFP[t->Mult] = 0;
		}

		t->weight[t->Mult] = weight;

		hdNdPtdEta[bin]->Fill(t->Eta[t->Mult], t->Pt[t->Mult]);
		hdNdPtdEtaPt[bin]->Fill(t->Eta[t->Mult], t->Pt[t->Mult], t->Pt[t->Mult]);

		t->Phi[t->Mult] = phi;
		hPt[bin]->Fill(t->Pt[t->Mult]);

		t->Mult++;
	}
}




// ------------ method called once each job just before starting event loop  ------------
	void 
QWCME::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
	void 
QWCME::endJob() 
{
}

// ------------ method called when starting to processes a run  ------------
	void 
QWCME::beginRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a run  ------------
	void 
QWCME::endRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when starting to processes a luminosity block  ------------
	void 
QWCME::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a luminosity block  ------------
	void 
QWCME::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
QWCME::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
	//The following says we do not know what parameters are allowed so do no validation
	// Please change this to state exactly what you do use, even if it is no parameters
	edm::ParameterSetDescription desc;
	desc.setUnknown();
	descriptions.addDefault(desc);

	//Specify that only 'tracks' is allowed
	//To use, remove the default given above and uncomment below
	//ParameterSetDescription desc;
	//desc.addUntracked<edm::InputTag>("tracks","ctfWithMaterialTracks");
	//descriptions.addDefault(desc);
}

//////////////////////////////////////////


//define this as a plug-in
DEFINE_FWK_MODULE(QWCME);
