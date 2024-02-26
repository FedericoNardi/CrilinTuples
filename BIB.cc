#define BIB_cxx
#include "BIB.h"
#include <TMath.h>
#include <TH2.h>
#include <vector>
#include <string>
#include <iostream>
#include <set>
#include <memory>
#include <TH2.h>
#include <TProfile.h>
#include <TStyle.h>
#include <TCanvas.h>
// #include <TApplication.h>

// A SET OF GLOBAL VARIABLES -----------------------------------------------

const std::vector<int> _bins_x = {79, 83, 85, 87, 89};
const std::vector<float> _halflength_x = {397.8, 418.2, 428.4, 438.6, 448.8};
const int _bins_z = 435;
const float _halflength_z = 2213.2;

bool DEBUG = false;
bool DRAW_HIST = true;
bool CUT_T = true;

// SOME HELPER FUNCTIONS ---------------------------------------------------

double calculateTheta(double x, double y)
{
   return TMath::ATan2(y, x);
}

double calculateMean(const std::vector<double> &data, const std::vector<double> &weights)
{
   double avg = 0;
   double norm = 0;
   for (unsigned int i = 0; i < data.size(); i++)
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

void BIB::Loop(int _wedge, bool DEBUG, bool DRAW_HIST, bool CUT_T)
{
   int wedge = _wedge;
   double t_cut = 10.;

   std::string fname_e = "/lustre/cmsdata/fnardi/MuColl_tuples/bib_layers/bib_e_layers.root";
   std::string fname_t = "/lustre/cmsdata/fnardi/MuColl_tuples/bib_layers/bib_t_layers.root";
   std::unique_ptr<TFile> ofile_e(new TFile(fname_e.c_str(), "UPDATE"));
   std::unique_ptr<TFile> ofile_t(new TFile(fname_t.c_str(), "UPDATE"));

   std::vector<TH2D*> bib_layers(5);
   for (int i_layer = 0; i_layer < 5; i_layer++)
   {
      std::string hname = "h_layer_" + std::to_string(wedge) + "_" + std::to_string(i_layer);
      bib_layers[i_layer] = new TH2D(hname.c_str(), hname.c_str(), _bins_z, -_halflength_z, _halflength_z, _bins_x[i_layer], -_halflength_x[i_layer], _halflength_x[i_layer]);
   }
   std::vector<TH2D*> tim_layers(5);
   for (int i_layer = 0; i_layer < 5; i_layer++)
   {
      std::string hname = "h_layer_t_" + std::to_string(wedge) + "_" + std::to_string(i_layer);
      tim_layers[i_layer] = new TH2D(hname.c_str(), hname.c_str(), _bins_z, -_halflength_z, _halflength_z, _bins_x[i_layer], -_halflength_x[i_layer], _halflength_x[i_layer]);
   }
   std::cout << "Allocated histograms\n";

   std::string hname0 = "h_deposit" + std::to_string(wedge);
   TProfile* deposit = new TProfile(
      hname0.c_str(), 
      hname0.c_str(), 
      _bins_z, -_halflength_z, _halflength_z, 
      -10, 10
   );

   std::string hname1 = "h_time" + std::to_string(wedge);
   TH1D* time = new TH1D(hname1.c_str(), hname1.c_str(), 100, -1.6, 1.6);

   if(DEBUG) std::cout << "Allocated histograms. Getting entries\n";

   Long64_t nentries = fChain->GetEntriesFast();
   if(DEBUG) std::cout << "Starting loop...\n";

   double radius_cut = 300;

   if(DEBUG) std::cout << "Starting loop\n";
   
   for (Long64_t jentry = 0; jentry < nentries; jentry++)
   {
      fChain->GetEntry(jentry);
      if(DEBUG) std::cout << "Entry #" << jentry << "\n";
      for (int i_hit = 0; i_hit < nsch; ++i_hit)
      {
         if(DEBUG) std::cout << "Hit #" << i_hit << "\n";
         if (scene[i_hit] > 0.000005) //&& (TMath::Abs(scpoz[i_hit]) < radius_cut))
         {
            if (assignWedge(scpox[i_hit], scpoy[i_hit]) != wedge) {continue;}
            double angle = (8. - wedge) * TMath::Pi() / 6;
            if(DEBUG) std::cout << "Angle: " << angle << "\n";
            std::vector<double> coords = Rotate(scpox[i_hit], scpoy[i_hit], angle);
            if(DEBUG) std::cout << "Coords: " << coords[0] << ", " << coords[1] << "\n";
            int layer_idx = int((coords[1] - 1497.5) / 45);
            if(DEBUG) std::cout << "Layer idx: " << layer_idx << "\n";
            double t_norm = (sctim[i_hit][0]-TMath::Sqrt(pow(scpox[i_hit],2)+pow(scpoy[i_hit],2)+pow(scpoz[i_hit],2))/299.792458);
            if(CUT_T) t_cut=0.250;
            if( TMath::Abs(t_norm)<t_cut )
            {
               time->Fill(t_norm, scene[i_hit]);
               deposit->Fill(scpoz[i_hit],t_norm);
               bib_layers[layer_idx]->Fill(scpoz[i_hit], coords[0], scene[i_hit]);
               // fill tim_layers if there is no hit in the same bin or if the time is smaller
               if(tim_layers[layer_idx]->GetBinContent(tim_layers[layer_idx]->FindBin(scpoz[i_hit], coords[0])) == 0 || tim_layers[layer_idx]->GetBinContent(tim_layers[layer_idx]->FindBin(scpoz[i_hit], coords[0])) > t_norm)
               {
                  tim_layers[layer_idx]->SetBinContent(tim_layers[layer_idx]->FindBin(scpoz[i_hit], coords[0]), t_norm);
               }
            }
         }
      }
   }
   if(DEBUG) std::cout << "Loop over\n";

   if(DRAW_HIST){
   
      TCanvas* c1 = new TCanvas("c1", "c1", 1200, 600);
      c1->Divide(2);
      c1->cd(1);
      deposit->SetStats(kFALSE);
      deposit->SetTitle("hit time");
      deposit->GetXaxis()->SetTitle("z [mm]");
      deposit->GetYaxis()->SetTitle("time [ns]");
      deposit->Draw("hist");
      c1->cd(2);
      time->SetStats(kFALSE);
      time->SetTitle("energy deposition vs time");
      time->GetXaxis()->SetTitle("normalized hit time [ns]");
      time->GetYaxis()->SetTitle("energy deposition [GeV]");
      time->Draw("hist");
      if(CUT_T) c1->Print("img/timing_data_cut.png");
      else c1->Print("img/timing_data.png");


      // Draw the timing distribution of each cell stored in tim_layers
      TH1D* h_time_dist = nullptr;
      if(CUT_T) h_time_dist = new TH1D("h_time_d", "h_time_d", 100, -0.5, 0.5);
      else h_time_dist = new TH1D("h_time_d", "h_time_d", 100, -5., 1.5);
      // Raise Error if h_time_dist is not created
      if(h_time_dist == nullptr) throw std::runtime_error("Error: h_time_dist is not created");
      TCanvas* c2 = new TCanvas("c2", "c2", 800, 600);
      for(int iLayer=0; iLayer<5; iLayer++){
         for(int iBin=0; iBin<tim_layers[iLayer]->GetNbinsX(); iBin++){
            for(int jBin=0; jBin<tim_layers[iLayer]->GetNbinsY(); jBin++){
               h_time_dist->Fill(tim_layers[iLayer]->GetBinContent(iBin, jBin), bib_layers[iLayer]->GetBinContent(iBin, jBin));
            }
         }
      }
      h_time_dist->SetStats(kFALSE);
      h_time_dist->SetTitle("timing distribution");
      h_time_dist->GetXaxis()->SetTitle("normalized hit time [ns]");
      h_time_dist->GetYaxis()->SetTitle("energy [GeV]");
      h_time_dist->Draw("hist");
      if(CUT_T) c2->Print("img/timing_voxels_cut.png");
      else c2->Print("img/timing_voxels.png");  
   }  

   if(!DRAW_HIST){ // if the histograms are not drawn, save them to file
      ofile_e->cd();
      std::cout << "Writing to file ENERGY \n";
      for(int i_layer=0; i_layer<5; i_layer++)
      {
         bib_layers[i_layer]->Write();
      }
      ofile_t->cd();
      std::cout << "Writing to file TIME \n";
      for(int i_layer=0; i_layer<5; i_layer++)
      {
         tim_layers[i_layer]->Write();
      }
      ofile_t->Close();
      ofile_e->Close();

      std::cout << "Files closed\n";
   }
   // Deleting the histograms
   for(int i_layer=0; i_layer<5; i_layer++)
   {
      if(DEBUG) std::cout << "Deleted layer " << i_layer << "\n";
      delete bib_layers[i_layer];
      delete tim_layers[i_layer];
   }
   bib_layers.clear();
   tim_layers.clear();
   std::cout << "Exiting\n";
   return;
}

int main()
{  
   BIB bib;
   for(int wedge = 0; wedge < 12; wedge++) bib.Loop(wedge,0,0,1);
   return 0;
}
