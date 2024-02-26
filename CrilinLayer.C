#define converter_cxx
#include "CrilinLayer.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <vector>
#include <iostream>
#include <set>
#include <TMath.h>

// A SET OF GLOBAL VARIABLES -----------------------------------------------

const std::vector<int> _bins_x = {79, 83, 85, 87, 89};
const std::vector<float> _halflength_x = {397.8, 418.2, 428.4, 438.6, 448.8};
const int _bins_z = 435;
const float _halflength_z = 2213.2;
std::vector<std::vector<double>> _pars;

// SOME HELPER FUNCTIONS ---------------------------------------------------

double calculateTheta(double x, double y){
   return TMath::ATan2(y,x);
}

double calculateMean(const std::vector<double>& data, const std::vector<double>& weights){
   double avg = 0;
   double norm = 0;
   for(int i=0; i<data.size(); i++){
      avg+=data[i]*weights[i];
      norm+=weights[i];
   }
   return avg/norm;
}

double calculateWeightedStd(const std::vector<double>& values, const std::vector<double>& weights) {
    // Step 1: Calculate the weighted mean
    double weightedSum = 0.0;
    double weightSum = 0.0;
    for (size_t i = 0; i < values.size(); ++i) {
        weightedSum += values[i] * weights[i];
        weightSum += weights[i];
    }
    double weightedMean = weightedSum / weightSum;

    // Step 2: Calculate the sum of weighted squared differences from the weighted mean
    double weightedSquaredDiffSum = 0.0;
    for (size_t i = 0; i < values.size(); ++i) {
        double diff = values[i] - weightedMean;
        weightedSquaredDiffSum += weights[i] * diff * diff;
    }

    // Step 3: Calculate the weighted variance and weighted standard deviation
    double weightedVariance = weightedSquaredDiffSum / weightSum;
    double weightedStandardDeviation = std::sqrt(weightedVariance);

    return weightedStandardDeviation;
}

int assignWedge( double x, double y ){
    //
    // Assigns an index from 0 to 11 according to which edge the hit belongs to. 0 starting from -175 deg and proceeding clockwise
    // 
    double angle_rad = TMath::ATan2( y,x )-TMath::Pi()/12; // Range goes from -180 to 180 deg
    double index = (angle_rad+TMath::Pi())/(TMath::Pi()/6);
    return index>=0? index : 11;
}



double calculateBIB(double z, int ilayer){
    z *= 1e-3; // switch mm to m
    std::vector<double> pars = _pars[ilayer];
    double zsquare = z*z;
    double bib_depo = pars[0]*TMath::Exp( -zsquare/(pars[1]*pars[1]) )+pars[2]*zsquare+pars[3]*z+pars[4];
    bib_depo/=788.7;
    return bib_depo;
}

std::vector<int> unique(const std::vector<int>& inputVector) {
    std::vector<int> uniqueVector;
    std::set<int> uniqueSet;

    for (const int& value : inputVector) {
        if (uniqueSet.insert(value).second) {
            uniqueVector.push_back(value);
        }
    }

    return uniqueVector;
}

std::vector<double> Rotate( double x, double y, double angle ){
    std::vector<double> new_coords;
    new_coords.push_back( x*TMath::Cos(angle) - y*TMath::Sin(angle) );
    new_coords.push_back( x*TMath::Sin(angle) + y*TMath::Cos(angle) );
    return new_coords;
}

// CODE STARTS HERE! ----------------------------------------------------------

void converter::Loop()
{
   if (fChain == nullptr) return;

   Long64_t nentries = fChain->GetEntriesFast();
   Long64_t nbytes = 0, nb = 0;

   _pars = {
       {4.306977147942193938e+02, 2.307000459016088634e+00, 2.135263279401148040e+01, 1.355287070784944481e+00, -2.174233985138803007e+02},
       {5.675413618161136401e+04, 1.244668165318582354e+01, 3.495415009520628473e+02, 5.053136358620166080e-01, -5.667741847020672139e+04},
       {4.319225913733152083e+03, 1.149909885732619585e+01, 2.658189961307489924e+01, 2.730465576907189496e-01, -4.277767894649226037e+03},
       {-3.111283195165126472e+02, 4.700709378041818098e+00, -1.644565416793840029e+01, 2.085226596876513239e-01, 3.410961506055436416e+02},
       {-3.772776440519966127e+04, 1.301875270987308930e+01, -2.233485924720556000e+02, 2.325673499469212979e-01, 3.775870371152652660e+04}
   };

   std::vector<TH2D*> layers(5);
   for (int i_layer = 0; i_layer < 5; i_layer++) {
       layers[i_layer] = new TH2D("hlayer", "h_layer", _bins_z, -_halflength_z, _halflength_z, _bins_x[i_layer], -_halflength_x[i_layer], _halflength_x[i_layer]);
   }

   // Initialize output vectors outside the loop
   std::vector<double> hit_x;
   std::vector<double> hit_y;
   std::vector<double> hit_z;
   std::vector<double> hit_dE;
   std::vector<double> evt_dE;
   std::vector<double> recHit_dE;
   std::vector<double> initial_energy;
   std::vector<int> isSignal;

   for (int i_file = 0; i_file < 50; i_file++) {

        string fname = "/lustre/cmsdata/fnardi/MuColl_tuples/" + std::to_string(G4energy) + "GeV/converted_" + std::to_string(G4energy) + "GeV_debug_" + std::to_string(i_file) + ".root";

        TFile* ofile = new TFile(fname.c_str(), "RECREATE");

        double radius_cut = 200;

        TTree* tree = new TTree("converted_photons", "converted photons");
        tree->Branch("hit_x", &hit_x);
        tree->Branch("hit_y", &hit_y);
        tree->Branch("hit_z", &hit_z);
        tree->Branch("hit_dE", &hit_dE);
        tree->Branch("recHit_dE", &recHit_dE);
        tree->Branch("evt_dE", &evt_dE);
        tree->Branch("isSignal", &isSignal);
        tree->Branch("photon_E", &initial_energy);

        for (Long64_t jentry=20*i_file; jentry<20*(i_file+1); jentry++){
            
            Long64_t ientry = LoadTree(jentry);
            if (ientry < 0) break;
            nb = fChain->GetEntry(jentry);
            nbytes += nb;

            std::vector<double> sig_x, sig_y, sig_z, sig_dE;
            std::vector<int> wedges;

            // Pre-compute theta_mean and theta_std
            std::vector<double> theta;
            std::vector<double> energies;

            for (int i_hit = 0; i_hit < nsch; ++i_hit) {
                if (scene[i_hit] > 0.0005) { // CUT at 500keV
                    theta.push_back(calculateTheta(scpox[i_hit], scpoy[i_hit]));
                    energies.push_back(scene[i_hit]);
                }
            }
        
            double theta_mean = calculateMean(theta, energies);
            double theta_std = calculateWeightedStd(theta, energies);

            for (int i_hit = 0; i_hit < nsch; ++i_hit) {
                double theta = calculateTheta(scpox[i_hit], scpoy[i_hit]);
                if (scene[i_hit] > 0.0005 && (TMath::Abs(theta - theta_mean) < theta_std)) {
                   sig_x.push_back(scpox[i_hit]);
                   sig_y.push_back(scpoy[i_hit]);
                   sig_z.push_back(scpoz[i_hit]);
                   sig_dE.push_back(scene[i_hit]);

                   wedges.push_back(assignWedge(scpox[i_hit], scpoy[i_hit]));
                }
            }

            auto unique_wedges = unique(wedges);

            for (const auto& w : unique_wedges) {
                double angle = (8. - w) * TMath::Pi() / 6;

                for (int i_hits = 0; i_hits < sig_x.size(); ++i_hits) {
                    if (wedges[i_hits] == w) {
                        std::vector<double> coords = Rotate(sig_x[i_hits], sig_y[i_hits], angle);

                        int layer_idx = int((coords[1] - 1497.5) / 45);

                        layers[layer_idx]->Fill(sig_z[i_hits], coords[0], sig_dE[i_hits]);
                    }
                }

                for (int i_layer = 0; i_layer < 5; ++i_layer) {
                    for (int z_idx = 0; z_idx < layers[i_layer]->GetNbinsX(); ++z_idx) {
                        for (int x_idx = 0; x_idx < layers[i_layer]->GetNbinsY(); ++x_idx) {

                            double dE_signal = layers[i_layer]->GetBinContent(z_idx, x_idx);
                            double z = layers[i_layer]->GetXaxis()->GetBinCenter(z_idx);
                            double dE_bib = calculateBIB(z, i_layer);
                            double dE_voxel = dE_signal + dE_bib;

                            if (dE_voxel < 0.0005) continue; // Cut at 500 keV

                            double inner_x = layers[i_layer]->GetYaxis()->GetBinCenter(x_idx);

                            std::vector<double> final_xy = Rotate(inner_x, 1510. + i_layer * 45., -angle);

                            if (TMath::Abs(calculateTheta(final_xy[0], final_xy[1]) - theta_mean) < theta_std && TMath::Abs(z) < 200) {
                                hit_x.push_back(final_xy[0]);
                                hit_y.push_back(final_xy[1]);
                                hit_z.push_back(z);
                                recHit_dE.push_back(dE_signal);
                                hit_dE.push_back(dE_voxel);
                                isSignal.push_back(dE_signal > dE_bib ? 0 : -1);
                                evt_dE.push_back(layers[i_layer]->Integral());
                                initial_energy.push_back(G4energy);
                            }
                        }
                    }
                    layers[i_layer]->Reset();
                }
            }
            tree->Fill();
            
            hit_x.clear();
            hit_y.clear();
            hit_z.clear();
            recHit_dE.clear();
            hit_dE.clear();
            isSignal.clear();
            evt_dE.clear();
            initial_energy.clear();
            std::cout << "...>Processed event " << jentry << "\n";
            
        }
        tree->Write();
        ofile->Close();
        std::cout << "--> Wrote file #" << i_file << " for " << G4energy << " GeV\n";
   }
}