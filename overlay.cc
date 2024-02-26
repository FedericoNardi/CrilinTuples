#define converter_new_cxx
#include "overlay.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TH1.h>
#include <TLegend.h>
#include <vector>
#include <memory>
#include <iostream>
#include <set>
#include <TMath.h>

// A SET OF GLOBAL VARIABLES -----------------------------------------------

const std::vector<int> _bins_x = {79, 83, 85, 87, 89};
const std::vector<float> _halflength_x = {397.8, 418.2, 428.4, 438.6, 448.8};
const int _bins_z = 435;
const float _halflength_z = 2213.2;

bool DEBUG = true;

// SOME HELPER FUNCTIONS ---------------------------------------------------

double calculateTheta(double x, double y)
{
    return TMath::ATan2(y, x);
}

double calculateMean(const std::vector<double> &data, const std::vector<double> &weights)
{
    double avg = 0;
    double norm = 0;
    for (int i = 0; i < data.size(); i++)
    {
        avg += data[i] * weights[i];
        norm += weights[i];
    }
    return avg / norm;
}

double calculateWeightedStd(const std::vector<double> &values, const std::vector<double> &weights)
{
    // Step 1: Calculate the weighted mean
    double weightedSum = 0.0;
    double weightSum = 0.0;
    for (size_t i = 0; i < values.size(); ++i)
    {
        weightedSum += values[i] * weights[i];
        weightSum += weights[i];
    }
    double weightedMean = weightedSum / weightSum;

    // Step 2: Calculate the sum of weighted squared differences from the weighted mean
    double weightedSquaredDiffSum = 0.0;
    for (size_t i = 0; i < values.size(); ++i)
    {
        double diff = values[i] - weightedMean;
        weightedSquaredDiffSum += weights[i] * diff * diff;
    }

    // Step 3: Calculate the weighted variance and weighted standard deviation
    double weightedVariance = weightedSquaredDiffSum / weightSum;
    double weightedStandardDeviation = std::sqrt(weightedVariance);

    return weightedStandardDeviation;
}

int assignWedge(double x, double y)
{
    //
    // Assigns an index from 0 to 11 according to which edge the hit belongs to. 0 starting from -175 deg and proceeding clockwise
    //
    double angle_rad = TMath::ATan2(y, x) - TMath::Pi() / 12; // Range goes from -180 to 180 deg
    double index = (angle_rad + TMath::Pi()) / (TMath::Pi() / 6);
    return index >= 0 ? index : 11;
}

std::vector<int> unique(const std::vector<int> &inputVector)
{
    std::vector<int> uniqueVector;
    std::set<int> uniqueSet;

    for (const int &value : inputVector)
    {
        if (uniqueSet.insert(value).second)
        {
            uniqueVector.push_back(value);
        }
    }

    return uniqueVector;
}

std::vector<double> Rotate(double x, double y, double angle)
{
    std::vector<double> new_coords;
    new_coords.push_back(x * TMath::Cos(angle) - y * TMath::Sin(angle));
    new_coords.push_back(x * TMath::Sin(angle) + y * TMath::Cos(angle));
    return new_coords;
}

// CODE STARTS HERE! ----------------------------------------------------------

void converter_new::Loop()
{
    if (fChain == nullptr)
        return;

    Long64_t nentries = fChain->GetEntriesFast();
    Long64_t nbytes = 0, nb = 0;

    std::vector<TH2D *> layers(5);
    for (int i_layer = 0; i_layer < 5; i_layer++)
    {
        layers[i_layer] = new TH2D("hlayer", "h_layer", _bins_z, -_halflength_z, _halflength_z, _bins_x[i_layer], -_halflength_x[i_layer], _halflength_x[i_layer]);
    }

    std::vector<TH2D *> layers_t(5);
    for (int i_layer = 0; i_layer < 5; i_layer++)
    {
        layers_t[i_layer] = new TH2D("hlayer_t", "h_layer_t", _bins_z, -_halflength_z, _halflength_z, _bins_x[i_layer], -_halflength_x[i_layer], _halflength_x[i_layer]);
    }

    // Initialize output vectors outside the loop
    std::vector<double> hit_x;
    std::vector<double> hit_y;
    std::vector<double> hit_z;
    std::vector<double> hit_t;
    std::vector<double> hit_t_sig;
    std::vector<double> hit_dE;
    std::vector<double> evt_dE;
    std::vector<double> recHit_dE;
    std::vector<double> initial_energy;
    std::vector<int> isSignal;
    std::vector<double> signalFraction;

    for (int i_file = 0; i_file < 100; i_file++)
    {
        if(DEBUG) std::cout << "-> Running on file: " << i_file << "\n";

        string fname = "/lustre/cmsdata/fnardi/MuColl_tuples/timing/photons_" + std::to_string(i_file) + ".root";

        TFile *ofile = new TFile(fname.c_str(), "RECREATE");

        double radius_cut = 200;

        TTree *tree = new TTree("converted_photons", "converted photons");
        tree->Branch("hit_x", &hit_x);              // mm
        tree->Branch("hit_y", &hit_y);              // mm
        tree->Branch("hit_z", &hit_z);              // mm
        tree->Branch("hit_t", &hit_t);              // ns
        tree->Branch("hit_t_sig", &hit_t_sig);      // ns
        tree->Branch("hit_dE", &hit_dE);            // GeV
        tree->Branch("recHit_dE", &recHit_dE);
        tree->Branch("evt_dE", &evt_dE);
        tree->Branch("isSignal", &isSignal);
        tree->Branch("primary_E", &initial_energy);
        tree->Branch("signalFraction", &signalFraction);

        for (Long64_t jentry = 100 * i_file; jentry < 100 * (i_file + 1); jentry++)
        {
            // if(jentry>0){continue;}  // Cut for DEBUG
            Long64_t ientry = LoadTree(jentry);
            if (ientry < 0)
                break;
            nb = fChain->GetEntry(jentry);
            nbytes += nb;

            if(DEBUG) std::cout << "--> Processing event " << jentry << " with energy: " << primary_energy << "\n";

            std::vector<double> sig_x, sig_y, sig_z, sig_dE, sig_t;
            std::vector<int> wedges;

            // Pre-compute theta_mean and theta_std
            std::vector<double> theta;
            std::vector<double> energies;

            for (int i_hit = 0; i_hit < nsch; ++i_hit)
            {
                if (scene[i_hit] > 5E-6) // GeV
                { // CUT 
                    theta.push_back(calculateTheta(scpox[i_hit], scpoy[i_hit]));
                    energies.push_back(scene[i_hit]);
                }
            }

            double theta_mean = calculateMean(theta, energies);
            double theta_std = calculateWeightedStd(theta, energies);

            if(DEBUG) std::cout << "--> Theta mean: " << theta_mean << " Theta std: " << theta_std << "\n";

            for (int i_hit = 0; i_hit < nsch; ++i_hit)
            {
                double theta = calculateTheta(scpox[i_hit], scpoy[i_hit]);
                if (scene[i_hit] > 0.0005 && (TMath::Abs(theta - theta_mean) < theta_std))
                {
                    sig_x.push_back(scpox[i_hit]);
                    sig_y.push_back(scpoy[i_hit]);
                    sig_z.push_back(scpoz[i_hit]);
                    sig_dE.push_back(scene[i_hit]);
                    sig_t.push_back(sctim[i_hit][0]);
                    wedges.push_back(assignWedge(scpox[i_hit], scpoy[i_hit]));
                }
            }

            auto unique_wedges = unique(wedges);

            if (unique_wedges.size() > 3)
            {
                continue;
            }

            TFile *bib_file = new TFile("/lustre/cmsdata/fnardi/MuColl_tuples/bib_layers/bib_e_layers.root", "READ");
            TFile *bib_t_file = new TFile("/lustre/cmsdata/fnardi/MuColl_tuples/bib_layers/bib_t_layers.root", "READ");

            for (const auto &w : unique_wedges)
            { // EDITED unique
                if(DEBUG) std::cout << "---> Processing wedge " << w << "\n";
                std::vector<TH2D *> bib_layers(5);
                std::vector<TH2D *> bib_t_layers(5);
                for (int i_layer = 0; i_layer < 5; i_layer++)
                {
                    if(DEBUG) std::cout << "Getting: " << Form("h_layer_%d_%d", w, i_layer) << "\n";
                    bib_layers[i_layer] = dynamic_cast<TH2D *>(bib_file->Get(Form("h_layer_%d_%d", w, i_layer))->Clone());
                    if(DEBUG) std::cout << "energy hist " << w << '_' << i_layer << '_' << bib_file->Get(Form("h_layer_%d_%d", w, i_layer))->Clone() << std::endl;
                    bib_t_layers[i_layer] = dynamic_cast<TH2D *>(bib_t_file->Get(Form("h_layer_%d_%d", w, i_layer))->Clone());
                    if(DEBUG) std::cout << "time hist " << w << '_' << i_layer << '_' << bib_t_layers[i_layer]->GetName() << std::endl;
                }
                if(DEBUG) std::cout << "---> Got bib layers\n";
                for (int i_hits = 0; i_hits < sig_x.size(); ++i_hits)
                {
                    if (wedges[i_hits] == w)
                    {
                        // Filling energies and time for each layer
                        double angle = (8. - w) * TMath::Pi() / 6;
                        std::vector<double> coords = Rotate(sig_x[i_hits], sig_y[i_hits], angle);

                        int layer_idx = int((coords[1] - 1497.5) / 45);

                        layers[layer_idx]->Fill(sig_z[i_hits], coords[0], sig_dE[i_hits]);
                        layers_t[layer_idx]->Fill(sig_z[i_hits], coords[0], sig_t[i_hits]);
                    }
                } 
                auto h_tim_signal   = new TH1D("h_tim_signal", "h_tim_signal", 100, -0.6, 1.5);
                auto h_tim_bib      = new TH1D("h_tim_bib", "h_tim_bib", 100, -0.6, 1.5);
                auto h_tim_overall  = new TH1D("h_tim_overall", "h_tim_overall", 100, -0.6, 1.5);
                

                for (int i_layer = 0; i_layer < 5; ++i_layer)
                {
                    for (int z_idx = 0; z_idx < layers[i_layer]->GetNbinsX(); ++z_idx)
                    {
                        for (int x_idx = 0; x_idx < layers[i_layer]->GetNbinsY(); ++x_idx)
                        {

                            double dE_signal = layers[i_layer]->GetBinContent(z_idx, x_idx);
                            double z = layers[i_layer]->GetXaxis()->GetBinCenter(z_idx);
                            double dE_bib = bib_layers[i_layer]->GetBinContent(z_idx, x_idx);
                            double dE_voxel = dE_signal + dE_bib;
                            double t_sig = layers_t[i_layer]->GetBinContent(z_idx, x_idx);
                            double t_bib = bib_t_layers[i_layer]->GetBinContent(z_idx, x_idx);

                            if (dE_voxel < 0.0000005)
                                continue; // Cut at 500 keV

                            double inner_x = layers[i_layer]->GetYaxis()->GetBinCenter(x_idx);

                            double angle = (8. - w) * TMath::Pi() / 6;
                            std::vector<double> final_xy = Rotate(inner_x, 1510. + i_layer * 45., -angle);

                            // if (TMath::Abs(z) < 250 && TMath::Abs(calculateTheta(final_xy[0], final_xy[1]) - theta_mean) < theta_std)
                            {
                                hit_x.push_back(final_xy[0]);
                                hit_y.push_back(final_xy[1]);
                                hit_z.push_back(z);
                                recHit_dE.push_back(dE_signal);
                                hit_dE.push_back(dE_voxel);
                                hit_t.push_back(t_sig < t_bib ? t_sig : t_bib);
                                hit_t_sig.push_back(t_sig);
                                isSignal.push_back(dE_signal > dE_bib ? 0 : -1);
                                evt_dE.push_back(layers[i_layer]->Integral());
                                initial_energy.push_back(primary_energy);
                                signalFraction.push_back(dE_signal / dE_voxel);
                                 if(jentry == 1){
                                    double r = sqrt( pow(final_xy[0],2)+pow(final_xy[1],2)+pow(z,2) );
                                    double shift = r/3.E2;
                                    h_tim_signal->Fill(t_sig-shift, dE_signal);
                                    h_tim_bib->Fill(t_bib, dE_bib);
                                    auto t_hit = t_sig*1e-2-shift < t_bib ? t_sig*1e-2 : t_bib;
                                    h_tim_overall->Fill(t_hit, dE_signal);
                                    //if(DEBUG){
                                        std::cout << t_sig-shift << " " << t_sig << " " << shift << " " << t_bib << std::endl;
                                    //}
                                }
                            }
                        }
                    }
                    layers[i_layer]->Reset();
                    layers_t[i_layer]->Reset();
                }

                // Print timing histogram for first event only
                if (jentry == 1)
                {
                    TCanvas *c1 = new TCanvas("c1", "c1", 800, 600);
                    // c1->SetLogy();
                    // Title
                    h_tim_bib->SetTitle("BIB timing");
                    //h_tim_signal->SetLineColor(858);
                    
                    h_tim_bib->SetLineColor(795);
                    h_tim_bib->Draw("hist");
                    /*
                    h_tim_overall->SetLineColor(kSpring-6);
                    h_tim_overall->SetFillColor(kSpring-5);
                    h_tim_overall->SetFillStyle(3001);
                    // no stats
                    h_tim_signal->SetStats(0);
                    h_tim_bib->SetStats(0);
                    h_tim_overall->SetStats(0);
                    // Draw histograms
                    h_tim_signal->Draw("hist");
                    h_tim_bib->Draw("hist same");
                    h_tim_overall->Draw("hist same");
                    
                    // Legend top right
                    TLegend *leg = new TLegend(0.7, 0.7, 0.9, 0.9);
                    leg->AddEntry(dynamic_cast<TH1D*>(h_tim_signal), "Signal", "l");
                    leg->AddEntry(dynamic_cast<TH1D*>(h_tim_bib), "BIB", "l");
                    leg->AddEntry(dynamic_cast<TH1D*>(h_tim_overall), "Overall", "f");
                    // smaller font
                    leg->SetTextSize(0.03);
                    // Adapt box size
                    leg->SetFillColor(0);
                    leg->SetFillStyle(0);
                    leg->SetBorderSize(0);
                    // legend location top right
                    leg->Draw();
                    */
                    
                    c1->SaveAs("timing_hist_bib.png");
                    std::cout << "Total signal deposited: " << h_tim_signal->Integral() << std::endl;
                }
                // delete histograms
                // delete h_tim_signal;
                // delete h_tim_bib;
                // delete h_tim_overall;

                if(DEBUG) std::cout << "---> Reset layers \n";
                for (auto &elem : bib_layers)
                {
                    delete elem;
                }
                for (auto &elem : bib_t_layers)
                {
                    delete elem;
                }
                bib_layers.clear();
                bib_t_layers.clear();
            }
            tree->Fill();
            hit_x.clear();
            hit_y.clear();
            hit_z.clear();
            hit_t.clear();
            hit_dE.clear();
            recHit_dE.clear();
            evt_dE.clear();
            initial_energy.clear();
            isSignal.clear();
            if(DEBUG) std::cout << "...>Processed event " << jentry << "\n";
            bib_file->Close();
            bib_t_file->Close();
        }
        ofile->cd();
        tree->Write();
        ofile->Close();

        std::cout << "--> Wrote file #" << i_file << "\n";
    }
}