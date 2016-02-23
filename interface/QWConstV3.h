#include <TMath.h>

const double Pi = TMath::Pi();
const double Pi2 = 2*Pi;

//PbPb centrality binning
const double centbins[]={0,5,10,15,20,25,30,35,40,50,60,70,80,90,100}; // nCentBins = 13
const Int_t nCentBins = sizeof(centbins)/sizeof(double)-1;

const double ptbins[25] = {
                0.3,  0.5,  1.0,  1.25,  1.5,  2.0,  2.5,  3.0,  3.5,   4.0, 5.0, 6.0, 7.0, 8.0,
		10.0, 12.0, 14.0, 20.0, 26.0, 35.0, 45.0, 60.0, 80.0, 100.0, 1000000.0}; // pPb pt binning nPtBins = 24;
const Int_t nPtBins = sizeof(ptbins)/sizeof(double)-1;

const double etabins[] = {
	        -2.4, -2.2, -2, -1.8, -1.6, -1.4, -1.2, -1, -0.8, -0.6, -0.4, -0.2, 0, 0.2, 0.4, 0.6, 0.8, 1, 1.2, 1.4, 1.6, 1.8, 2, 2.2, 2.4
};
const Int_t nEtaBins = sizeof(etabins)/sizeof(double)-1;


//const Int_t CentNoffCut[] = {100000, 350, 320, 300, 260, 240, 220, 185, 150, 120, 100, 80, 60, 50, 40, 30, 20, 10, 0};
//const Int_t nCentNoff = sizeof(CentNoffCut)/sizeof(Int_t);
