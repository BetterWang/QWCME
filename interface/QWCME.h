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
		double	dpp_2p;
		double	dnn_2n;
		double	dpp   ;
		double	dnn   ;
		double	dp_2p ;
		double	dn_2n ;
		double	d_2p  ;
		double	d_2n  ;
		double	dp    ;
		double	dn    ;

		double ipp_2p;
		double inn_2n;
		double ipp;
		double inn;
		double ip_2p;
		double in_2n;
		double i_2p;
		double i_2n;
		double ip;
		double in;

		double wpp_2p;
		double wnn_2n;
		double wpp;
		double wnn;
		double wp_2p;
		double wn_2n;
		double w_2p;
		double w_2n;
		double wp;
		double wn;


		correlations::QVector		qpp_2p;
		correlations::QVector		qnn_2n;
		correlations::QVector		qpp;
		correlations::QVector		qnn;
		correlations::QVector		qp_2p;
		correlations::QVector		qn_2n;
		correlations::QVector		q_2p;
		correlations::QVector		q_2n;
		correlations::QVector		qp;
		correlations::QVector		qn;

		correlations::HarmonicVector	hpp_2p;
		correlations::HarmonicVector	hpp;
		correlations::HarmonicVector	hp_2p;
		correlations::HarmonicVector	h_2p;
		correlations::HarmonicVector	hp;

		TTree * trV;


		bool CaloMatch(const reco::Track&, const edm::Event&, unsigned int idx);
		int getNoffCent(const edm::Event&, const edm::EventSetup&, int& Noff);
};



