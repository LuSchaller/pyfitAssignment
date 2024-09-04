
// Where to find link to resolutions and config files:
// https://gitlab.desy.de/patrick.connor/uhh-top-mass/-/blob/alljets_run2/Configuration/analysis/analyzeTopHypotheses_alljets_cfg.py
// https://gitlab.desy.de/patrick.connor/uhh-top-mass/-/blob/alljets_run2/Configuration/python/stringResolutions_etEtaPhi_modified_cff.py

// Where to find link to setup of covariance matrix:
// https://github.com/cms-sw/cmssw/blob/CMSSW_8_0_X/TopQuarkAnalysis/TopKinFitter/src/CovarianceMatrix.cc#L138
// https://github.com/cms-sw/cmssw/blob/CMSSW_8_0_X/TopQuarkAnalysis/TopObjectResolutions/interface/Jet.h 


double CalcEta(double et, double eta, bool isB )
{
  double res= 0.0;
  if(isB){
// eta b
		if(std::abs(eta)<0.087) res = sqrt(pow(0.00605,2) + pow((1.63/et),2));
		else if(std::abs(eta)<0.174) res = sqrt(pow(0.00592,2) + pow((1.64/et),2));
		else if(std::abs(eta)<0.261) res = sqrt(pow(0.00584,2) + pow((1.65/et),2));
		else if(std::abs(eta)<0.348) res = sqrt(pow(0.00593,2) + pow((1.65/et),2));
		else if(std::abs(eta)<0.435) res = sqrt(pow(0.00584,2) + pow((1.68/et),2));
		else if(std::abs(eta)<0.522) res = sqrt(pow(0.00646,2) + pow((1.67/et),2));
		else if(std::abs(eta)<0.609) res = sqrt(pow(0.00661,2) + pow((1.67/et),2));
		else if(std::abs(eta)<0.696) res = sqrt(pow(0.00724,2) + pow((1.65/et),2));
		else if(std::abs(eta)<0.783) res = sqrt(pow(0.00763,2) + pow((1.67/et),2));
		else if(std::abs(eta)<0.870) res = sqrt(pow(0.00746,2) + pow((1.7/et),2));
		else if(std::abs(eta)<0.957) res = sqrt(pow(0.00807,2) + pow((1.7/et),2));
		else if(std::abs(eta)<1.044) res = sqrt(pow(0.00843,2) + pow((1.72/et),2));
		else if(std::abs(eta)<1.131) res = sqrt(pow(0.00886,2) + pow((1.74/et),2));
		else if(std::abs(eta)<1.218) res = sqrt(pow(0.0101,2) + pow((1.76/et),2));
		else if(std::abs(eta)<1.305) res = sqrt(pow(0.0127,2) + pow((1.78/et),2));
		else if(std::abs(eta)<1.392) res = sqrt(pow(0.0161,2) + pow((1.73/et),2));
		else if(std::abs(eta)<1.479) res = sqrt(pow(0.0122,2) + pow((1.77/et),2));
		else if(std::abs(eta)<1.566) res = sqrt(pow(0.0123,2) + pow((1.79/et),2));
		else if(std::abs(eta)<1.653) res = sqrt(pow(0.0111,2) + pow((1.79/et),2));
		else if(std::abs(eta)<1.740) res = sqrt(pow(0.0106,2) + pow((1.8/et),2));
		else if(std::abs(eta)<1.830) res = sqrt(pow(0.0115,2) + pow((1.83/et),2));
		else if(std::abs(eta)<1.930) res = sqrt(pow(0.012,2) + pow((1.88/et),2));
		else if(std::abs(eta)<2.043) res = sqrt(pow(0.0131,2) + pow((1.91/et),2));
		else if(std::abs(eta)<2.172) res = sqrt(pow(0.0134,2) + pow((1.92/et),2));
		else if(std::abs(eta)<2.322) res = sqrt(pow(0.0132,2) + pow((1.89/et),2));
		else /*if(std::abs(eta)<2.500)*/ res = sqrt(pow(0.0121,2) + pow((2.28/et),2));
	}else{
// eta udsc
		if (std::abs(eta)<0.087) res = sqrt(pow(0.00915,2) + pow((1.51/et),2));
		else if (std::abs(eta)<0.174) res = sqrt(pow(0.00887,2) + pow((1.53/et),2));
		else if (std::abs(eta)<0.261) res = sqrt(pow(0.00865,2) + pow((1.54/et),2));
		else if (std::abs(eta)<0.348) res = sqrt(pow(0.00867,2) + pow((1.55/et),2));
		else if (std::abs(eta)<0.435) res = sqrt(pow(0.00907,2) + pow((1.55/et),2));
		else if (std::abs(eta)<0.522) res = sqrt(pow(0.00844,2) + pow((1.59/et),2));
		else if (std::abs(eta)<0.609) res = sqrt(pow(0.00915,2) + pow((1.57/et),2));
		else if (std::abs(eta)<0.696) res = sqrt(pow(0.00856,2) + pow((1.58/et),2));
		else if (std::abs(eta)<0.783) res = sqrt(pow(0.00897,2) + pow((1.58/et),2));
		else if (std::abs(eta)<0.870) res = sqrt(pow(0.0095,2) + pow((1.6/et),2));
		else if (std::abs(eta)<0.957) res = sqrt(pow(0.00836,2) + pow((1.65/et),2));
		else if (std::abs(eta)<1.044) res = sqrt(pow(0.00782,2) + pow((1.68/et),2));
		else if (std::abs(eta)<1.131) res = sqrt(pow(0.0093,2) + pow((1.65/et),2));
		else if (std::abs(eta)<1.218) res = sqrt(pow(0.00986,2) + pow((1.69/et),2));
		else if (std::abs(eta)<1.305) res = sqrt(pow(0.0124,2) + pow((1.72/et),2));
		else if (std::abs(eta)<1.392) res = sqrt(pow(0.0181,2) + pow((1.63/et),2));
		else if (std::abs(eta)<1.479) res = sqrt(pow(0.0121,2) + pow((1.69/et),2));
		else if (std::abs(eta)<1.566) res = sqrt(pow(0.0122,2) + pow((1.69/et),2));
		else if (std::abs(eta)<1.653) res = sqrt(pow(0.00975,2) + pow((1.69/et),2));
		else if (std::abs(eta)<1.740) res = sqrt(pow(0.00881,2) + pow((1.71/et),2));
		else if (std::abs(eta)<1.830) res = sqrt(pow(0.00938,2) + pow((1.75/et),2));
		else if (std::abs(eta)<1.930) res = sqrt(pow(0.00894,2) + pow((1.8/et),2));
		else if (std::abs(eta)<2.043) res = sqrt(pow(0.00893,2) + pow((1.83/et),2));
		else if (std::abs(eta)<2.172) res = sqrt(pow(0.0099,2) + pow((1.82/et),2));
		else if (std::abs(eta)<2.322) res = sqrt(pow(0.00944,2) + pow((1.8/et),2));
		else /*if (std::abs(eta)<2.500)*/ res = sqrt(pow(0.0103,2) + pow((2.15/et),2));
	}
  return res;
}

double CalcEt(double et, double eta, bool isB)
{
  double res = 0.0; 
  if(isB){ 
// et b
		if(std::abs(eta)<0.087) res = et * (sqrt(pow(0.0686,2) + pow((1.03/sqrt(et)),2) + pow((3.43837/et),2)));
		else if(std::abs(eta)<0.174) res = et * (sqrt(pow(0.0737,2) + pow((1.01/sqrt(et)),2) + pow((3.46808/et),2)));
		else if(std::abs(eta)<0.261) res = et * (sqrt(pow(0.0657,2) + pow((1.07/sqrt(et)),2) + pow((3/et),2)));
		else if(std::abs(eta)<0.348) res = et * (sqrt(pow(0.062,2) + pow((1.07/sqrt(et)),2) + pow((3/et),2)));
		else if(std::abs(eta)<0.435) res = et * (sqrt(pow(0.0605,2) + pow((1.07/sqrt(et)),2) + pow((3/et),2)));
		else if(std::abs(eta)<0.522) res = et * (sqrt(pow(0.059,2) + pow((1.08/sqrt(et)),2) + pow((3/et),2)));
		else if(std::abs(eta)<0.609) res = et * (sqrt(pow(0.0577,2) + pow((1.08/sqrt(et)),2) + pow((3/et),2)));
		else if(std::abs(eta)<0.696) res = et * (sqrt(pow(0.0525,2) + pow((1.09/sqrt(et)),2) + pow((3/et),2)));
		else if(std::abs(eta)<0.783) res = et * (sqrt(pow(0.0582,2) + pow((1.09/sqrt(et)),2) + pow((3/et),2)));
		else if(std::abs(eta)<0.870) res = et * (sqrt(pow(0.0649,2) + pow((1.08/sqrt(et)),2) + pow((3/et),2)));
		else if(std::abs(eta)<0.957) res = et * (sqrt(pow(0.0654,2) + pow((1.1/sqrt(et)),2) + pow((3/et),2)));
		else if(std::abs(eta)<1.044) res = et * (sqrt(pow(0.0669,2) + pow((1.11/sqrt(et)),2) + pow((3/et),2)));
		else if(std::abs(eta)<1.131) res = et * (sqrt(pow(0.0643,2) + pow((1.15/sqrt(et)),2) + pow((3/et),2)));
		else if(std::abs(eta)<1.218) res = et * (sqrt(pow(0.0645,2) + pow((1.16/sqrt(et)),2) + pow((3/et),2)));
		else if(std::abs(eta)<1.305) res = et * (sqrt(pow(0.0637,2) + pow((1.19/sqrt(et)),2) + pow((3/et),2)));
		else if(std::abs(eta)<1.392) res = et * (sqrt(pow(0.0695,2) + pow((1.21/sqrt(et)),2) + pow((3/et),2)));
		else if(std::abs(eta)<1.479) res = et * (sqrt(pow(0.0748,2) + pow((1.2/sqrt(et)),2) + pow((3/et),2)));
		else if(std::abs(eta)<1.566) res = et * (sqrt(pow(0.0624,2) + pow((1.23/sqrt(et)),2) + pow((3/et),2)));
		else if(std::abs(eta)<1.653) res = et * (sqrt(pow(0.0283,2) + pow((1.25/sqrt(et)),2) + pow((3/et),2)));
		else if(std::abs(eta)<1.740) res = et * (sqrt(pow(0.0316,2) + pow((1.21/sqrt(et)),2) + pow((3/et),2)));
		else if(std::abs(eta)<1.830) res = et * (sqrt(pow(2.29e-07,2) + pow((1.2/sqrt(et)),2) + pow((3/et),2)));
		else if(std::abs(eta)<1.930) res = et * (sqrt(pow(5.18e-09,2) + pow((1.14/sqrt(et)),2) + pow((3.44819/et),2)));
		else if(std::abs(eta)<2.043) res = et * (sqrt(pow(2.17e-07,2) + pow((1.09/sqrt(et)),2) + pow((3.65053/et),2)));
		else if(std::abs(eta)<2.172) res = et * (sqrt(pow(3.65e-07,2) + pow((1.09/sqrt(et)),2) + pow((3.41422/et),2)));
		else if(std::abs(eta)<2.322) res = et * (sqrt(pow(2.02e-07,2) + pow((1.09/sqrt(et)),2) + pow((3.43837/et),2)));
		else /*if(std::abs(eta)<2.500)*/ res = et * (sqrt(pow(5.27e-07,2) + pow((1.12/sqrt(et)),2) + pow((3.48832/et),2)));
	}else{
// et udsc
		if (std::abs(eta)<0.087) res = et * (sqrt(pow(0.0591,2) + pow((1/sqrt(et)),2) + pow((3.12952/et),2)));
		else if (std::abs(eta)<0.174) res = et * (sqrt(pow(0.0619,2) + pow((0.975/sqrt(et)),2) + pow((3.37218/et),2)));
		else if (std::abs(eta)<0.261) res = et * (sqrt(pow(0.0574,2) + pow((1/sqrt(et)),2) + pow((3/et),2)));
		else if (std::abs(eta)<0.348) res = et * (sqrt(pow(0.0569,2) + pow((1.01/sqrt(et)),2) + pow((3/et),2)));
		else if (std::abs(eta)<0.435) res = et * (sqrt(pow(0.057,2) + pow((1/sqrt(et)),2) + pow((3/et),2)));
		else if (std::abs(eta)<0.522) res = et * (sqrt(pow(0.0522,2) + pow((1.02/sqrt(et)),2) + pow((3/et),2)));
		else if (std::abs(eta)<0.609) res = et * (sqrt(pow(0.0502,2) + pow((1.02/sqrt(et)),2) + pow((3/et),2)));
		else if (std::abs(eta)<0.696) res = et * (sqrt(pow(0.053,2) + pow((1.03/sqrt(et)),2) + pow((3/et),2)));
		else if (std::abs(eta)<0.783) res = et * (sqrt(pow(0.051,2) + pow((1.03/sqrt(et)),2) + pow((3/et),2)));
		else if (std::abs(eta)<0.870) res = et * (sqrt(pow(0.0549,2) + pow((1.04/sqrt(et)),2) + pow((3/et),2)));
		else if (std::abs(eta)<0.957) res = et * (sqrt(pow(0.0544,2) + pow((1.06/sqrt(et)),2) + pow((3/et),2)));
		else if (std::abs(eta)<1.044) res = et * (sqrt(pow(0.0519,2) + pow((1.09/sqrt(et)),2) + pow((3/et),2)));
		else if (std::abs(eta)<1.131) res = et * (sqrt(pow(0.0539,2) + pow((1.12/sqrt(et)),2) + pow((3/et),2)));
		else if (std::abs(eta)<1.218) res = et * (sqrt(pow(0.0492,2) + pow((1.16/sqrt(et)),2) + pow((3/et),2)));
		else if (std::abs(eta)<1.305) res = et * (sqrt(pow(0.0489,2) + pow((1.18/sqrt(et)),2) + pow((3/et),2)));
		else if (std::abs(eta)<1.392) res = et * (sqrt(pow(0.0414,2) + pow((1.25/sqrt(et)),2) + pow((3/et),2)));
		else if (std::abs(eta)<1.479) res = et * (sqrt(pow(0.0373,2) + pow((1.26/sqrt(et)),2) + pow((3/et),2)));
		else if (std::abs(eta)<1.566) res = et * (sqrt(pow(0.0125,2) + pow((1.24/sqrt(et)),2) + pow((3/et),2)));
		else if (std::abs(eta)<1.653) res = et * (sqrt(pow(1.37e-07,2) + pow((1.08/sqrt(et)),2) + pow((4.28528/et),2)));
		else if (std::abs(eta)<1.740) res = et * (sqrt(pow(2.37e-07,2) + pow((1.04/sqrt(et)),2) + pow((4.24972/et),2)));
		else if (std::abs(eta)<1.830) res = et * (sqrt(pow(2.3e-07,2) + pow((1/sqrt(et)),2) + pow((4.31393/et),2)));
		else if (std::abs(eta)<1.930) res = et * (sqrt(pow(1.25e-07,2) + pow((0.965/sqrt(et)),2) + pow((4.34276/et),2)));
		else if (std::abs(eta)<2.043) res = et * (sqrt(pow(5.78e-08,2) + pow((0.924/sqrt(et)),2) + pow((4.34276/et),2)));
		else if (std::abs(eta)<2.172) res = et * (sqrt(pow(4.25e-08,2) + pow((0.923/sqrt(et)),2) + pow((4.13793/et),2)));
		else if (std::abs(eta)<2.322) res = et * (sqrt(pow(0.00601,2) + pow((0.881/sqrt(et)),2) + pow((4.40828/et),2)));
		else /*if (std::abs(eta)<2.500)*/ res = et * (sqrt(pow(4.94e-08,2) + pow((0.86/sqrt(et)),2) + pow((4.65549/et),2)));
	}
	return res;
}

double CalcPhi(double et, double eta,  bool isB)
{
  double res = 0.0;
  if(isB){ 
// phi b
		if(std::abs(eta)<0.087) res = sqrt(pow(0.00787,2) + pow((1.74/et),2));
		else if(std::abs(eta)<0.174) res = sqrt(pow(0.00766,2) + pow((1.74/et),2));
		else if(std::abs(eta)<0.261) res = sqrt(pow(0.00755,2) + pow((1.74/et),2));
		else if(std::abs(eta)<0.348) res = sqrt(pow(0.00734,2) + pow((1.74/et),2));
		else if(std::abs(eta)<0.435) res = sqrt(pow(0.00734,2) + pow((1.75/et),2));
		else if(std::abs(eta)<0.522) res = sqrt(pow(0.00767,2) + pow((1.74/et),2));
		else if(std::abs(eta)<0.609) res = sqrt(pow(0.00742,2) + pow((1.75/et),2));
		else if(std::abs(eta)<0.696) res = sqrt(pow(0.00771,2) + pow((1.73/et),2));
		else if(std::abs(eta)<0.783) res = sqrt(pow(0.00758,2) + pow((1.76/et),2));
		else if(std::abs(eta)<0.870) res = sqrt(pow(0.00789,2) + pow((1.75/et),2));
		else if(std::abs(eta)<0.957) res = sqrt(pow(0.00802,2) + pow((1.76/et),2));
		else if(std::abs(eta)<1.044) res = sqrt(pow(0.0078,2) + pow((1.79/et),2));
		else if(std::abs(eta)<1.131) res = sqrt(pow(0.00806,2) + pow((1.82/et),2));
		else if(std::abs(eta)<1.218) res = sqrt(pow(0.00784,2) + pow((1.86/et),2));
		else if(std::abs(eta)<1.305) res = sqrt(pow(0.00885,2) + pow((1.9/et),2));
		else if(std::abs(eta)<1.392) res = sqrt(pow(0.0108,2) + pow((1.93/et),2));
		else if(std::abs(eta)<1.479) res = sqrt(pow(0.0112,2) + pow((2/et),2));
		else if(std::abs(eta)<1.566) res = sqrt(pow(0.0102,2) + pow((2.02/et),2));
		else if(std::abs(eta)<1.653) res = sqrt(pow(0.00857,2) + pow((2.01/et),2));
		else if(std::abs(eta)<1.740) res = sqrt(pow(0.00856,2) + pow((1.97/et),2));
		else if(std::abs(eta)<1.830) res = sqrt(pow(0.00761,2) + pow((1.95/et),2));
		else if(std::abs(eta)<1.930) res = sqrt(pow(0.00721,2) + pow((1.92/et),2));
		else if(std::abs(eta)<2.043) res = sqrt(pow(0.00722,2) + pow((1.86/et),2));
		else if(std::abs(eta)<2.172) res = sqrt(pow(0.00703,2) + pow((1.86/et),2));
		else if(std::abs(eta)<2.322) res = sqrt(pow(0.00845,2) + pow((1.86/et),2));
		else /*if(std::abs(eta)<2.500)*/ res = sqrt(pow(0.00975,2) + pow((2/et),2));
	}else{
// phi udsc
		if(std::abs(eta)<0.087) res = sqrt(pow(0.01,2) + pow((1.6/et),2));
		else if(std::abs(eta)<0.174) res = sqrt(pow(0.00982,2) + pow((1.61/et),2));
		else if(std::abs(eta)<0.261) res = sqrt(pow(0.0101,2) + pow((1.59/et),2));
		else if(std::abs(eta)<0.348) res = sqrt(pow(0.00988,2) + pow((1.6/et),2));
		else if(std::abs(eta)<0.435) res = sqrt(pow(0.0102,2) + pow((1.59/et),2));
		else if(std::abs(eta)<0.522) res = sqrt(pow(0.00982,2) + pow((1.6/et),2));
		else if(std::abs(eta)<0.609) res = sqrt(pow(0.00979,2) + pow((1.6/et),2));
		else if(std::abs(eta)<0.696) res = sqrt(pow(0.00925,2) + pow((1.62/et),2));
		else if(std::abs(eta)<0.783) res = sqrt(pow(0.00973,2) + pow((1.61/et),2));
		else if(std::abs(eta)<0.870) res = sqrt(pow(0.00971,2) + pow((1.62/et),2));
		else if(std::abs(eta)<0.957) res = sqrt(pow(0.00916,2) + pow((1.64/et),2));
		else if(std::abs(eta)<1.044) res = sqrt(pow(0.00959,2) + pow((1.66/et),2));
		else if(std::abs(eta)<1.131) res = sqrt(pow(0.00964,2) + pow((1.67/et),2));
		else if(std::abs(eta)<1.218) res = sqrt(pow(0.00969,2) + pow((1.71/et),2));
		else if(std::abs(eta)<1.305) res = sqrt(pow(0.00992,2) + pow((1.76/et),2));
		else if(std::abs(eta)<1.392) res = sqrt(pow(0.0124,2) + pow((1.79/et),2));
		else if(std::abs(eta)<1.479) res = sqrt(pow(0.0135,2) + pow((1.8/et),2));
		else if(std::abs(eta)<1.566) res = sqrt(pow(0.0107,2) + pow((1.85/et),2));
		else if(std::abs(eta)<1.653) res = sqrt(pow(0.00895,2) + pow((1.84/et),2));
		else if(std::abs(eta)<1.740) res = sqrt(pow(0.00902,2) + pow((1.81/et),2));
		else if(std::abs(eta)<1.830) res = sqrt(pow(0.00861,2) + pow((1.79/et),2));
		else if(std::abs(eta)<1.930) res = sqrt(pow(0.00877,2) + pow((1.75/et),2));
		else if(std::abs(eta)<2.043) res = sqrt(pow(0.00791,2) + pow((1.73/et),2));
		else if(std::abs(eta)<2.172) res = sqrt(pow(0.00775,2) + pow((1.73/et),2));
		else if(std::abs(eta)<2.322) res = sqrt(pow(0.00807,2) + pow((1.71/et),2));
		else /*if(std::abs(eta)<2.500)*/ res = sqrt(pow(0.0103,2) + pow((1.81/et),2));
	}
	return res;
}

// where to find the eta dependantscalefactors: 
// https://gitlab.desy.de/patrick.connor/uhh-top-mass/-/blob/alljets_run2/Configuration/analysis/analyzeTopHypotheses_alljets_cfg.py#L540
// https://gitlab.desy.de/patrick.connor/uhh-top-mass/-/blob/alljets_run2/Configuration/analysis/analyzeTopHypotheses_alljets_cfg.py#L570
double getEtaDependentScaleFactor(const TLorentzVector& object){
	double sf = 1.0;
	const double eta = object.Eta();
	if(std::abs(eta)<0.5) sf = 1.109;
	else if(std::abs(eta)<0.8) sf = 1.138;
	else if(std::abs(eta)<1.1) sf = 1.114;
	else if(std::abs(eta)<1.3) sf = 1.123;
	else if(std::abs(eta)<1.7) sf = 1.084;
	else if(std::abs(eta)<1.9) sf = 1.082;
	else if(std::abs(eta)<2.1) sf = 1.140;
	else if(std::abs(eta)<2.3) sf = 1.067;
	else if(std::abs(eta)<2.5) sf = 1.177;
	else if(std::abs(eta)<2.8) sf = 1.364;
	else if(std::abs(eta)<3.0) sf = 1.857;
	else if(std::abs(eta)<3.2) sf = 1.328;
	else sf = 1.16;
	
	return sf;
}
