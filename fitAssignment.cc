#include "PhysicsTools/KinFitter/interface/TAbsFitParticle.h"
#include "PhysicsTools/KinFitter/interface/TFitConstraintM.h"
#include "PhysicsTools/KinFitter/interface/TFitParticleEtEtaPhi.h"
#include "PhysicsTools/KinFitter/interface/TKinFitter.h"
#include "Math/VectorUtil.h"
#include "kinFit.h"
#include <TDecompLU.h>
#include <TFile.h>
#include <TLorentzVector.h>
#include <TTree.h>
#include <TMatrixD.h>
#include <TString.h>
#include <TTree.h>
#include <functional>
#include <array>
#include <cassert>
#include <unistd.h>
#include <iostream>
#include <vector>
#include <map>
using ROOT::Math::VectorUtil::DeltaR;
using namespace std;

void tryCombinations(TLorentzVector *B1, TLorentzVector *B2, TLorentzVector *J1, TLorentzVector *J2, TLorentzVector *J3, TLorentzVector *J4, TMatrixD mB1, TMatrixD mB2, TMatrixD mJ1, TMatrixD mJ2, TMatrixD mJ3, TMatrixD mJ4, struct Selection &bestSelection, double dRLimit, TString option) {

    applyKinFit({B1, B2, J1, J2, J3, J4}, {mB1, mB2, mJ1, mJ2, mJ3, mJ4}, bestSelection, dRLimit, option);
    applyKinFit({B1, B2, J1, J3, J2, J4}, {mB1, mB2, mJ1, mJ3, mJ2, mJ4}, bestSelection, dRLimit, option);
    applyKinFit({B1, B2, J1, J4, J2, J3}, {mB1, mB2, mJ1, mJ4, mJ2, mJ3}, bestSelection, dRLimit, option);
    applyKinFit({B1, B2, J3, J4, J1, J2}, {mB1, mB2, mJ3, mJ4, mJ1, mJ2}, bestSelection, dRLimit, option);
    applyKinFit({B1, B2, J2, J4, J1, J3}, {mB1, mB2, mJ2, mJ4, mJ1, mJ3}, bestSelection, dRLimit, option);
    applyKinFit({B1, B2, J2, J3, J1, J4}, {mB1, mB2, mJ2, mJ3, mJ1, mJ4}, bestSelection, dRLimit, option);
}

void setbestcombi(vector<vector<double>> inputdata, vector<vector<double>> gendata, double dRLimit, string option) {
    // inputdata has the following structure:
    // inputdata[0] = pt
    // inputdata[1] = eta
    // inputdata[2] = phi
    // inputdata[3] = M
    // inputdata[4] = btag
    // GEN LEVEL INFO IS NEEDED FOR MATCHING AND HAS SAME STRUCTURE

    vector<double> inputpt = inputdata[0];
    vector<double> inputeta = inputdata[1];
    vector<double> inputphi = inputdata[2];
    vector<double> inputM = inputdata[3];


    vector<double> genpt = gendata[0];
    vector<double> geneta = gendata[1];
    vector<double> genphi = gendata[2];
    vector<double> genM = gendata[3];


    	struct Selection bestSelection;
	    bestSelection.chi2 = 10000;
	int N_permutations = 0;
	vector<TLorentzVector *> jets;
	for (size_t i = 0; i < inputdata[0].size(); i++) {
        TLorentzVector *jet = new TLorentzVector();
        jet->SetPtEtaPhiE(inputpt[i], inputeta[i], inputphi[i], inputM[i]);
        jets.push_back(jet);
    }
//maybe two seperate vectors from cf one for two btagged and one for 4 other jets
//otherwise we need btags from cf and sort here
	multimap<double, TLorentzVector *, greater<double>> btags;
	for (size_t i = 0; i < jets.size(); i++) {
	    TLorentzVector *jet = jets[i];
	    double btag = inputdata[4][i];
	    jets.push_back(jet);
	    btags.insert({btag, jet});
	}
	/* THIS WILL BE MOVED TO COLUMNFLOW
	// sort jets by Pt and only take the 6 with highest Pt
	sort(jets.begin(), jets.end(), [](TLorentzVector *&a, TLorentzVector *&b) {
	    return a->Pt() > b->Pt();
	});
	*/

	//	jets.resize(6);


	// try all 6 physical combinations of jets
	auto it = btags.begin();
	TLorentzVector *B1 = it->second;
	jets.erase(std::find(jets.begin(), jets.end(), B1));
	it++;
	TLorentzVector *B2 = it->second;
	jets.erase(std::find(jets.begin(), jets.end(), B2));

	TLorentzVector* J1 = jets[0];
	TLorentzVector* J2 = jets[1];
	TLorentzVector* J3 = jets[2];
	TLorentzVector* J4 = jets[3];
	  /* THIS WOULD BE FOR EXTRA JETS
	TLorentzVector* J5 = 0;
	TMatrixD mJ5(3,3);
	if(jets.size() > 4){
		 J5 = jets[4];
	mJ5 = setupMatrix(J5, false);
	}

	TLorentzVector* J6 = nullptr;
	TMatrixD mJ6(3,3);
	if(jets.size() > 5){
		J6 = jets[5];
		mJ6 = setupMatrix(J6,false);
	}
*/	
      //check if btagged jets are close enough IDK WHY THIS IS HERE
    /*
    if (DeltaR(*B1, *B2) <= 2) {
	    continue;
	}
	*/

	TMatrixD mB1 = setupMatrix(B1, true);
	TMatrixD mB2 = setupMatrix(B2, true);
	TMatrixD mJ1 = setupMatrix(J1, false);
	TMatrixD mJ2 = setupMatrix(J2, false);
	TMatrixD mJ3 = setupMatrix(J3, false);
	TMatrixD mJ4 = setupMatrix(J4, false);

	tryCombinations(B1, B2, J1, J2, J3, J4, mB1, mB2, mJ1, mJ2, mJ3, mJ4, bestSelection, dRLimit, option);
/*	
	if (jets.size() > 4) {
	tryCombinations(B1, B2, J5, J2, J3, J4,  mB1,  mB2,  mJ5,  mJ2,  mJ3,  mJ4, bestSelection, dRLimit, option);
	tryCombinations(B1, B2, J1, J5, J3, J4,  mB1,  mB2,  mJ1,  mJ5,  mJ3,  mJ4, bestSelection, dRLimit, option);
	tryCombinations(B1, B2, J1, J2, J5, J4,  mB1,  mB2,  mJ1,  mJ2,  mJ5,  mJ4, bestSelection, dRLimit, option);
	tryCombinations(B1, B2, J1, J2, J3, J5,  mB1,  mB2,  mJ1,  mJ2,  mJ3,  mJ5, bestSelection, dRLimit, option);
	}
*/	
      //MATCHING
 vector<TLorentzVector *> genJets;
 	for (size_t i = 0; i < gendata[0].size(); i++) {
        TLorentzVector *genjet = new TLorentzVector();
        genjet->SetPtEtaPhiM(genpt[i], geneta[i], genphi[i], genM[i]);
        genJets.push_back(genjet);
    }
    //check if the selection is correct with gen Information
      //combinationType ==1 means that the selection is correctly matched 
      //combinationType ==0 means that the selection is incorrectly matched
    //combinationType == -1 means that matching didnt work
    /*
    bestSelection.combinationType = -1;
    if (option == "ambiguous") {
	if (allmatch(genJets, jets)) {
	    bestSelection.combinationType = 1;
	} else {
	    bestSelection.combinationType = 0;
	}
    } else if (option == "unambiguous") {
	if (allmatchunambiguous(genJets, jets)) {
	    bestSelection.combinationType = 1;
	} else {
	    bestSelection.combinationType = 0;
	}
    }
   */
	//	mB1.Print();
	//	mB2.Print();
	//	mJ1.Print();
	//	mJ2.Print();
	//	mJ3.Print();
	//	mJ4.Print();

	
	cout << "Best Comb Chi2: " << bestSelection.chi2;
	cout << " Best Selection Pt  ";
	for (TLorentzVector *bestjet : bestSelection.bestPermutation) {
	    cout << bestjet->Pt() << " ";
	}
	cout << endl;

      //return struct for further work in CF
//return bestSelection;
}
