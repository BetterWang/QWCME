#include <correlations/Types.hh>
#include <correlations/Result.hh>
#include <correlations/QVector.hh>
#include <correlations/recursive/FromQVector.hh>
#include <correlations/recurrence/FromQVector.hh>
#include <correlations/closed/FromQVector.hh>
#include <TComplex.h>
#include <TH1.h>
#include <TH2.h>
#include <TTree.h>
#include <TNtupleD.h>
#include <TRandom3.h>
#include <TFile.h>
#include "QWConstV3.h"
#include <RecoHI/HiEvtPlaneAlgos/interface/HiEvtPlaneList.h>
//
// constants, enums and typedefs
//

//#define QW_DEBUG 1
//#define QW_PEREVENT 1

#define PRD(x) cout << "!!QW!! " << __LINE__ << " DEBUG OUTPUT " << (#x) << " = " << (x) << endl;
#define PR(x) cout << "!!QW!! " << __LINE__ << " DEBUG OUTPUT " << (#x) << endl;
//
// class declaration
//

const int NMAX_TRK = 10000;
typedef struct QWEvent_ {
	int     Cent;
	int     Mult;
	double  vz;
	int	Noff;
	double  Pt[NMAX_TRK];
	double  Eta[NMAX_TRK];
	double  Phi[NMAX_TRK];
	int     Charge[NMAX_TRK];
	double	rEff[NMAX_TRK];
	double	weight[NMAX_TRK];
	int	RFP[NMAX_TRK];
	int     RunId;
	int     EventId;
} QWEvent;

///////////////// Class ////////////////////////////

class QWCME : public edm::EDAnalyzer {
	public:
		explicit QWCME(const edm::ParameterSet&);
		~QWCME();

		static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

	private:
		virtual void beginJob() ;
		virtual void analyze(const edm::Event&, const edm::EventSetup&);
		void analyzeData(const edm::Event&, const edm::EventSetup&);
		void analyzeGen(const edm::Event&, const edm::EventSetup&);
		virtual void endJob() ;

		virtual void beginRun(edm::Run const&, edm::EventSetup const&);
		virtual void endRun(edm::Run const&, edm::EventSetup const&);
		virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);
		virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);

	/////////////////////////////////////////////
		//TRandom3 * gRandom;
		// ----------member data ---------------------------
		edm::InputTag					trackTag_;
		edm::EDGetTokenT<reco::TrackCollection>		trackToken_;
		edm::EDGetTokenT<reco::GenParticle>		trackGenToken_;
		edm::EDGetTokenT<int>				centralityToken_;
		edm::EDGetTokenT<reco::VertexCollection>	vertexToken_;
		edm::EDGetTokenT<reco::PFCandidateCollection>	pfToken_;

		edm::InputTag fweight_;
		edm::InputTag facceptance_;
	/////////////////////////////////////////////
		double	minvz_, maxvz_;
		double	dzdzerror_;
		double	d0d0error_;
		double	chi2_;
		double	pterrorpt_;

		bool	bacc;
		bool	bEff;
		bool	bPhiEta;
		bool	bCentNoff;
		bool	bCaloMatching_;
		int	Noffmin_;
		int	Noffmax_;
		int	cmode_;
		bool	bGen_;
		bool	bFlipEta_;
		std::vector<int> algoParameters_;

		unsigned int	nvtx_;

		double	reso_;
		QWEvent * t;
		TFile	* facc;
	/////////////////////////////////////////////
		TH1D * hPt[nCentBins];
		TH2D * hPhiEta[nCentBins][nPtBins][2];

		TH2D	* hdNdPtdEta[nCentBins];
		TH2D	* hdNdPtdEtaPt[nCentBins];

		TH2D * hEff_cbin[200];
		TH2D * hFak_cbin[200];

		TH2D * hacc[nCentBins][nPtBins][2];


		int gNoff;
		int gMult;

		correlations::FromQVector *	cqpp_2p;
		correlations::FromQVector *	cqnn_2n;
		correlations::FromQVector *	cqpp;
		correlations::FromQVector *	cqnn;
		correlations::FromQVector *	cq_2p;
		correlations::FromQVector *	cq_2n;
		correlations::FromQVector *	cqp_2p;
		correlations::FromQVector *	cqn_2n;


		correlations::QVector		qpp_2p;
		correlations::QVector		qnn_2n;
		correlations::QVector		qpp;
		correlations::QVector		qnn;
		correlations::QVector		q_2p;
		correlations::QVector		q_2n;
		correlations::QVector		qp_2p;
		correlations::QVector		qn_2n;

		correlations::HarmonicVector	hpp_2p;
		correlations::HarmonicVector	hnn_2n;
		correlations::HarmonicVector	hpp;
		correlations::HarmonicVector	hnn;
		correlations::HarmonicVector	h_2p;
		correlations::HarmonicVector	h_2n;
		correlations::HarmonicVector	hp_2p;
		correlations::HarmonicVector	hn_2n;


		void Sim();

		bool CaloMatch(const reco::Track&, const edm::Event&, unsigned int idx);
		int getNoffCent(const edm::Event&, const edm::EventSetup&, int& Noff);
};



