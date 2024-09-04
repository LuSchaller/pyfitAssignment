#include "PhysicsTools/KinFitter/interface/TAbsFitParticle.h"
#include "PhysicsTools/KinFitter/interface/TFitConstraintM.h"
#include "PhysicsTools/KinFitter/interface/TFitParticleEtEtaPhi.h"
#include "PhysicsTools/KinFitter/interface/TKinFitter.h"
//#include "addExtraJets.h"
#include "config.h"
//#include "matching.h"
#include <array>
#include <cassert>
#include <iostream>
#include <unistd.h>
#include <vector>

#include <TFile.h>
#include <TObject.h>
#include <TLorentzVector.h>
#include <TString.h>
#include <TTree.h>
using std::vector;
struct Selection {
    double chi2;
    vector<TLorentzVector *> bestPermutation;
    vector<TLorentzVector > fitJets;
    int combinationType;
};

TMatrixD setupMatrix(const TLorentzVector *object, bool isB) {
    TMatrixD CovM(3, 3);
    CovM.Zero();
    const double et = object->E() * std::fabs(sin(object->Theta()));
    const double eta = object->Eta();
    CovM(0, 0) = pow(CalcEt(et, eta, isB), 2);
    CovM(0, 0) *= pow(getEtaDependentScaleFactor(*object), 2);
    CovM(1, 1) = pow(CalcEta(et, eta, isB), 2);
    CovM(2, 2) = pow(CalcPhi(et, eta, isB), 2);
    return CovM;
}
// now overloaded to account for adding Jets
TMatrixD setupMatrix(const TLorentzVector *object1, const TLorentzVector *object2, bool isB) {
    TMatrixD CovM(3, 3);
    CovM.Zero();
    const double et1 = object1->E() * std::fabs(sin(object1->Theta()));
    const double eta1 = object1->Eta();
    const double et2 = object2->E() * std::fabs(sin(object2->Theta()));
    const double eta2 = object2->Eta();

    CovM(0, 0) = pow(CalcEt(et1, eta1, isB) * getEtaDependentScaleFactor(*object1), 2) + pow(CalcEt(et2, eta2, isB) * getEtaDependentScaleFactor(*object2), 2);
    CovM(1, 1) = pow(CalcEta(et1, eta1, isB), 2) + pow(CalcEta(et2, eta2, isB), 2);
    CovM(2, 2) = pow(CalcPhi(et1, eta1, isB), 2) + pow(CalcPhi(et2, eta2, isB), 2);
    return CovM;
}
void applyKinFit(vector<TLorentzVector *> jetSelection, vector<TMatrixD> decayprodmats, struct Selection &currentSelection, double dRLimit, TString option) {
    // create new TKinFitter
    TKinFitter *fitter_ = new TKinFitter("TopKinFitter", "TopKinFitter");
    fitter_->setMaxNbIter(200);
    fitter_->setMaxDeltaS(5e-5);
    fitter_->setMaxF(1e-4);
    fitter_->setVerbosity(0);

    // definiton of covariance matrices
    TMatrixD mB1 = decayprodmats[0];
    TMatrixD mB2 = decayprodmats[1];
    TMatrixD mW1P1 = decayprodmats[2];
    TMatrixD mW1P2 = decayprodmats[3];
    TMatrixD mW2P1 = decayprodmats[4];
    TMatrixD mW2P2 = decayprodmats[5];

    /* this would be for dealing with additional (<6) jets
    std::map<TLorentzVector *, TMatrixD> decayprodmap;
    decayprodmap.insert({jetSelection[0], mB1});
    decayprodmap.insert({jetSelection[1], mB2});
    decayprodmap.insert({jetSelection[2], mW1P1});
    decayprodmap.insert({jetSelection[3], mW1P2});
    decayprodmap.insert({jetSelection[4], mW2P1});
    decayprodmap.insert({jetSelection[5], mW2P2});
    vector<TLorentzVector *> extraJets;
   TLorentzVector addedjettodel;
    bool addevent = false;

    //cout << "All Pt pre Add : ";
    
    for (TLorentzVector *jet : extraJets) {
	bool matched = false;
	for (const auto &decayprod : decayprodmap) {
	    //	cout <<"DecayProd : " <<  decayprod.first->Pt()<< " jEt : " << jet->Pt() ;
	    if (abs((jet->Pt()) - decayprod.first->Pt()) > 0.0001)
		continue;
	    matched = true;
	    break;
	}
	if (matched) continue;
	TLorentzVector *closestdecayprodvecptr = 0;
	//cout << "New Jet :" <<closestdecayprodvecptr->Pt() << endl;
	bool hasclosestdecayprod = addExtraJets(*jet, jetSelection, dRLimit, closestdecayprodvecptr);
	if (hasclosestdecayprod) {
	    bool btag = false;
	    if (closestdecayprodvecptr == jetSelection[0] || closestdecayprodvecptr == jetSelection[1]) btag = true;
	    decayprodmap.at(closestdecayprodvecptr) = setupMatrix(closestdecayprodvecptr, jet, btag);
	    *closestdecayprodvecptr += *jet;
	    addedjettodel = *jet;
	    addevent = true;
	}
    }
    */
    // setup jetSelection and add to fitter_
    auto B1 = new TFitParticleEtEtaPhi("B1", "B1", jetSelection[0], &mB1);
    auto B2 = new TFitParticleEtEtaPhi("B2", "B2", jetSelection[1], &mB2);
    auto W1Prod1 = new TFitParticleEtEtaPhi("W1Prod1", "W1Prod1 ", jetSelection[2], &mW1P1);
    auto W1Prod2 = new TFitParticleEtEtaPhi("W1Prod2", "W1Prod2", jetSelection[3], &mW1P2);
    auto W2Prod1 = new TFitParticleEtEtaPhi("W2Prod1", "W2Prod1 ", jetSelection[4], &mW2P1);
    auto W2Prod2 = new TFitParticleEtEtaPhi("W2Prod2", "W2Prod2", jetSelection[5], &mW2P2);

    for (auto p : {B1, B2, W1Prod1, W1Prod2, W2Prod1, W2Prod2})
	fitter_->addMeasParticle(p);

    // set constants:
    double mW_ = 80.4;

    // set up constraints and add to fitter_
    auto kW1Mass = new TFitConstraintM("W1Mass", "W1Mass", 0, 0, mW_);
    auto kW2Mass = new TFitConstraintM("W2Mass", "W2Mass", 0, 0, mW_);
    auto kEqualTopMasses = new TFitConstraintM("EqualTopMasses", "EqualTopMasses", 0, 0, 0);

    // add particles/jetSelection to constraints
    kW1Mass->addParticles1(W1Prod1, W1Prod2);
    kW2Mass->addParticles1(W2Prod1, W2Prod2);
    kEqualTopMasses->addParticles1(B1, W1Prod1, W1Prod2);
    kEqualTopMasses->addParticles2(B2, W2Prod1, W2Prod2);

    for (auto c : {kW1Mass, kW2Mass, kEqualTopMasses})
	fitter_->addConstraint(c);
    // perform fit
    fitter_->fit();

    // ourfit has the newly fitted TLorentzvectors of the particles as entries, the order is the same
    // order the particles have been added to the fitter (B1,B2, W1Prod1, W1Prod2, W2Prod1, W2Prod2)


    // if fit didnt converge empty all top.reco vectors
    double thisChi2 = fitter_->getS();
    if (thisChi2 < currentSelection.chi2 && fitter_->getStatus() == 0) {
	    currentSelection.chi2 = thisChi2;
	    currentSelection.bestPermutation = jetSelection;
	    // add fitted jets to vector in order (B1,B2,W1Prod1,W1Prod2, W2Prod1, W2Prod2, W1, W2, Top1, Top2, TTBar)
	    for(int i=0 ; i<6 ; i++){
	     currentSelection.fitJets[i] = *(fitter_->get4Vec(i));
	    }
	currentSelection.fitJets[6] = currentSelection.fitJets[2] + currentSelection.fitJets[3];
	    currentSelection.fitJets[7] = currentSelection.fitJets[4] + currentSelection.fitJets[5];
	    currentSelection.fitJets[8] = currentSelection.fitJets[6] + currentSelection.fitJets[0];
	    currentSelection.fitJets[9] = currentSelection.fitJets[7] + currentSelection.fitJets[1];
	    currentSelection.fitJets[10] = currentSelection.fitJets[8] + currentSelection.fitJets[9];
	    }
	// overwrite reco if we find a better combination
	/* since we no longer deal with root files this is commented out
	topevent->recoB1.resize(1);
	topevent->recoB2.resize(1);
	topevent->recoW1Prod1.resize(1);
	topevent->recoW1Prod2.resize(1);
	topevent->recoW2Prod1.resize(1);
	topevent->recoW2Prod2.resize(1);
	topevent->recoB1[0] = *jetSelection[0];
	topevent->recoB2[0] = *jetSelection[1];
	topevent->recoW1Prod1[0] = *jetSelection[2];
	topevent->recoW1Prod2[0] = *jetSelection[3];
	topevent->recoW2Prod1[0] = *jetSelection[4];
	topevent->recoW2Prod2[0] = *jetSelection[5];

	// overwrite all fit TLorentzvectors
	topevent->fitB1.resize(1);
	topevent->fitB2.resize(1);
	topevent->fitW1Prod1.resize(1);
	topevent->fitW1Prod2.resize(1);
	topevent->fitW2Prod1.resize(1);
	topevent->fitW2Prod2.resize(1);
	topevent->fitW1.resize(1);
	topevent->fitW2.resize(1);
	topevent->fitTop1.resize(1);
	topevent->fitTop2.resize(1);
	topevent->fitTTBar.resize(1);
	topevent->fitB1[0] = *(fitter_->get4Vec(0));
	topevent->fitB2[0] = *(fitter_->get4Vec(1));
	topevent->fitW1Prod1[0] = *(fitter_->get4Vec(2));
	topevent->fitW1Prod2[0] = *(fitter_->get4Vec(3));
	topevent->fitW2Prod1[0] = *(fitter_->get4Vec(4));
	topevent->fitW2Prod2[0] = *(fitter_->get4Vec(5));
	topevent->fitW1[0] = topevent->fitW1Prod1[0] + topevent->fitW1Prod2[0];
	topevent->fitW2[0] = topevent->fitW2Prod1[0] + topevent->fitW2Prod2[0];
	topevent->fitTop2[0] = topevent->fitB1[0] + topevent->fitW1[0];
	topevent->fitTop1[0] = topevent->fitB2[0] + topevent->fitW2[0];
	topevent->fitTTBar[0] = topevent->fitTop1[0] + topevent->fitTop2[0];

	// overwrite fitChi2 and fitProb
	topevent->fitChi2.resize(1);
	topevent->fitProb.resize(1);
	topevent->fitChi2[0] = fitter_->getS();
	topevent->fitProb[0] = TMath::Prob(fitter_->getS(), fitter_->getNDF());
    }
    vector<TLorentzVector *> genJets;
    vector<TLorentzVector *> Jets;
    for (size_t ip = 0; ip < jetevent->genParton.size(); ip++)
	genJets.push_back(&(jetevent->genParton[ip]));
    //check if the selection is correct with gen Information
    topevent->combinationType.clear();
    topevent->combinationType.resize(1);
    if (option == "ambiguous") {
	if (allmatch(genJets, jetSelection)) {
	    topevent->combinationType[0] = 1;
	} else {
	    topevent->combinationType[0] = 0;
	}
    } else if (option == "unambiguous") {
	if (allmatchunambiguous(genJets, jetSelection,topevent)) {
	    topevent->combinationType[0] = 1;
	} else {
	    topevent->combinationType[0] = 0;
	}
    }*/
    delete fitter_;
    for (auto p : {B1, B2, W1Prod1, W1Prod2, W2Prod1, W2Prod2})
	delete p;
    for (auto c : {kW1Mass, kW2Mass, kEqualTopMasses})
	delete c;
}

