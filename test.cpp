//
// Author: Alexei Prokudin <prokudin@jlab.org>
//
//========================================================= includes
#include "./src/parameters.h"
#include "./src/unpolarised.h"
#include "./src/sivers.h"
#include "./src/transversity.h"
#include "./src/collins.h"
#include "./src/stfunctions.h"

using namespace std;

void HermesSivers()
{
   TCanvas *resultsivers = new TCanvas("Sivers Hermes", "Sivers Hermes",10,45,700,500);
   resultsivers->SetHighLightColor(2);
   resultsivers->Range(0,0,1,1);
   resultsivers->SetFillColor(0);
   resultsivers->SetFillStyle(4000);
   resultsivers->SetBorderMode(0);
   resultsivers->SetBorderSize(2);
   resultsivers->SetFrameFillStyle(0);
   resultsivers->SetFrameBorderMode(0);
   TLatex *   tex = new TLatex(0.08,0.6,"A_{UT}^{sin (#phi_{S} - #phi_{h})}");
   tex->SetTextSize(0.08);
   tex->SetTextAngle(90);
   tex->SetLineWidth(2);
   tex->Draw();
      tex = new TLatex(0.85,0.05,"x");
   tex->SetLineWidth(2);
   tex->Draw();
  
// ------------>Primitives in pad: nemo
   TPad *nemo = new TPad("nemo", "This is nemo",0.1,0.1,1,1);
   nemo->Draw();
   nemo->cd();
   nemo->Range(-0.0375,-0.1875,0.4275,0.1875);
   nemo->SetFillColor(0);
   nemo->SetFillStyle(4000);
   nemo->SetBorderMode(0);
   nemo->SetBorderSize(0);
   nemo->SetFrameFillStyle(0);
   nemo->SetFrameBorderMode(0);
   nemo->SetFrameFillStyle(0);
   nemo->SetFrameBorderMode(0);

   // Predictions for pi+
   TMD torino10;
   string target = "proton";
   string hadron = "pi+";

   double S  = pow2(7.25374);
   double z  = 0.36;
   double Q2 = 2.41;
   double PhT = 0.36;

   double x[7] = {0.036,0.056,0.076,0.098,0.133,0.186,0.275};
   double AUTpip[7];

   for (size_t i = 0; i < 7; i++)
   {
    AUTpip[i] = torino10.FUTSivers(target, hadron, S, x[i], z, Q2, PhT)/torino10.FUU(target, hadron, S, x[i], z, Q2, PhT);
    //cout << AUTpip[i] << endl;
   }

   TGraph* grtest = new TGraph(7,x,AUTpip);
   grtest->Draw("cp");
   grtest->SetLineColor(2);
   grtest->SetLineWidth(3);

   TH1F *Graph_Graph1 = new TH1F("Graph_Graph1","",100,0.009,0.381);
   Graph_Graph1->SetMinimum(-0.08);
   Graph_Graph1->SetMaximum(0.08);
   Graph_Graph1->SetDirectory(0);
   Graph_Graph1->SetStats(0);
   Graph_Graph1->GetXaxis()->SetNdivisions(8);
   Graph_Graph1->GetXaxis()->SetLabelFont(42);
   Graph_Graph1->GetXaxis()->SetLabelOffset(0.04);
   Graph_Graph1->GetXaxis()->SetLabelSize(0.05);
   Graph_Graph1->GetXaxis()->SetTitleSize(0.035);
   Graph_Graph1->GetXaxis()->SetTitleColor(0);
   Graph_Graph1->GetXaxis()->SetTitleFont(42);
   Graph_Graph1->GetYaxis()->SetNdivisions(8);
   Graph_Graph1->GetYaxis()->SetLabelFont(42);
   Graph_Graph1->GetYaxis()->SetLabelOffset(0.01);
   Graph_Graph1->GetYaxis()->SetLabelSize(0.05);
   Graph_Graph1->GetYaxis()->SetTitleSize(0.035);
   Graph_Graph1->GetYaxis()->SetTitleFont(42);
   Graph_Graph1->GetZaxis()->SetLabelFont(42);
   Graph_Graph1->GetZaxis()->SetLabelSize(0.035);
   Graph_Graph1->GetZaxis()->SetTitleSize(0.035);
   Graph_Graph1->GetZaxis()->SetTitleFont(42);

   grtest->SetHistogram(Graph_Graph1);
   
   grtest->Draw("acp");
   
   // Predictions for pi-
   hadron = "pi-";
   double AUTpim[7];

   for (size_t i = 0; i < 7; i++)
   {
    AUTpim[i] = torino10.FUTSivers(target, hadron, S, x[i], z, Q2, PhT)/torino10.FUU(target, hadron, S, x[i], z, Q2, PhT);
    //cout << AUTpim[i] << endl;
   }

   TGraph* grpim = new TGraph(7,x,AUTpim);
   grpim->Draw("cp");
   grpim->SetLineColor(4);
   grpim->SetLineWidth(3);
   
  
   TLine *line = new TLine(0.009,0,0.381,0);
   line->Draw();
   
   TPaveText *pt = new TPaveText(0.78,0.755,0.98,0.995,"brNDC");
   pt->SetName("description");
   pt->SetBorderSize(2);
   pt->SetFillColor(0);
   pt->SetTextAlign(12);
   pt->AddText("HERMES");
   pt->AddText("lP^{#uparrow}#rightarrow l' #pi^{#pm} X");
   pt->AddText(" #sqrt{s} = 7.25374 (GeV) ");
   pt->Draw();
   
   // pi+ Hermes data
   Double_t Graph2_fx1001[7] = {
   0.036,
   0.056,
   0.076,
   0.098,
   0.133,
   0.186,
   0.275};
   Double_t Graph2_fy1001[7] = {
   0.0385,
   0.0264,
   0.0542,
   0.0516,
   0.0361,
   0.0531,
   0.0546};
   Double_t Graph2_fex1001[7] = {
   0,
   0,
   0,
   0,
   0,
   0,
   0};
   Double_t Graph2_fey1001[7] = {
   0.008856636,
   0.007300685,
   0.008160882,
   0.007738863,
   0.008044874,
   0.01072007,
   0.01502431};
   TGraphErrors* gre = new TGraphErrors(7,Graph2_fx1001,Graph2_fy1001,Graph2_fex1001,Graph2_fey1001);
   gre->SetName("Graph2");
   gre->SetTitle("");
   gre->SetFillColor(1);
   gre->SetLineColor(4);
   gre->SetMarkerColor(4);
   gre->SetMarkerStyle(8);
   gre->SetMarkerSize(1.2);
   gre->Draw("p");

   // pi- Hermes data
   TGraphErrors *grem = new TGraphErrors(7);
   grem->SetName("Graph");
   grem->SetTitle("Graph");
   grem->SetFillColor(4);
   grem->SetLineColor(4);
   grem->SetMarkerColor(4);
   grem->SetMarkerStyle(24);
   grem->SetMarkerSize(1.2);
   grem->SetPoint(0,0.035,0.0104);
   grem->SetPointError(0,0,0.0095);
   grem->SetPoint(1,0.055,0.0169);
   grem->SetPointError(1,0,0.0080);
   grem->SetPoint(2,0.076,-0.0157);
   grem->SetPointError(2,0,0.0089);
   grem->SetPoint(3,0.098,0.0103);
   grem->SetPointError(3,0,0.0082);
   grem->SetPoint(4,0.133,-0.0053);
   grem->SetPointError(4,0,0.0079);
   grem->SetPoint(5,0.186, 0.0204);
   grem->SetPointError(5,0,0.0108);
   grem->SetPoint(6,0.275,-0.0192);
   grem->SetPointError(6,0, 0.0154);
   grem->Draw("p");
 
   nemo->Modified();
   resultsivers->cd();
   resultsivers->Modified();
   resultsivers->cd();
   resultsivers->SetSelected(resultsivers);
}

void HermesCollins()
{
   TCanvas *resultcollins = new TCanvas("Collins Hermes", "Collins Hermes",10,45,700,500);
   resultcollins->SetHighLightColor(2);
   resultcollins->Range(0,0,1,1);
   resultcollins->SetFillColor(0);
   resultcollins->SetFillStyle(4000);
   resultcollins->SetBorderMode(0);
   resultcollins->SetBorderSize(2);
   resultcollins->SetFrameFillStyle(0);
   resultcollins->SetFrameBorderMode(0);
   TLatex *   tex = new TLatex(0.08,0.6,"A_{UT}^{sin (#phi_{S} + #phi_{h})}");
   tex->SetTextSize(0.08);
   tex->SetTextAngle(90);
   tex->SetLineWidth(2);
   tex->Draw();
      tex = new TLatex(0.85,0.05,"x");
   tex->SetLineWidth(2);
   tex->Draw();
  
// ------------>Primitives in pad: nemo
   TPad *nemo = new TPad("nemo", "This is nemo",0.1,0.1,1,1);
   nemo->Draw();
   nemo->cd();
   nemo->Range(-0.0375,-0.1875,0.4275,0.1875);
   nemo->SetFillColor(0);
   nemo->SetFillStyle(4000);
   nemo->SetBorderMode(0);
   nemo->SetBorderSize(0);
   nemo->SetFrameFillStyle(0);
   nemo->SetFrameBorderMode(0);
   nemo->SetFrameFillStyle(0);
   nemo->SetFrameBorderMode(0);

   // Predictions for pi+
   TMD torino10;
   string target = "proton";
   string hadron = "pi+";

   double S  = pow2(7.25374);
   double z  = 0.36;
   double Q2 = 2.41;
   double PhT = 0.36;

   double x[7] = {0.036,0.056,0.076,0.098,0.133,0.186,0.275};
   double AUTpip[7];

   for (size_t i = 0; i < 7; i++)
   {
    double y = Q2/(S*x[i]);
    double DNN = (1.-y)/(1.+pow2(1-y));
    AUTpip[i] = DNN*torino10.FUTCollins(target, hadron, S, x[i], z, Q2, PhT)/torino10.FUU(target, hadron, S, x[i], z, Q2, PhT);
    //AUTpip[i] = torino10.FUTCollins(target, hadron, S, x[i], z, Q2, PhT)/torino10.FUU(target, hadron, S, x[i], z, Q2, PhT);
    //cout << AUTpip[i] << endl;
   }

   TGraph* grtest = new TGraph(7,x,AUTpip);
   grtest->Draw("cp");
   grtest->SetLineColor(2);
   grtest->SetLineWidth(3);

   TH1F *Graph_Graph1 = new TH1F("Graph_Graph1","",100,0.009,0.381);
   Graph_Graph1->SetMinimum(-0.12);
   Graph_Graph1->SetMaximum(0.08);
   Graph_Graph1->SetDirectory(0);
   Graph_Graph1->SetStats(0);
   Graph_Graph1->GetXaxis()->SetNdivisions(8);
   Graph_Graph1->GetXaxis()->SetLabelFont(42);
   Graph_Graph1->GetXaxis()->SetLabelOffset(0.04);
   Graph_Graph1->GetXaxis()->SetLabelSize(0.05);
   Graph_Graph1->GetXaxis()->SetTitleSize(0.035);
   Graph_Graph1->GetXaxis()->SetTitleColor(0);
   Graph_Graph1->GetXaxis()->SetTitleFont(42);
   Graph_Graph1->GetYaxis()->SetNdivisions(8);
   Graph_Graph1->GetYaxis()->SetLabelFont(42);
   Graph_Graph1->GetYaxis()->SetLabelOffset(0.01);
   Graph_Graph1->GetYaxis()->SetLabelSize(0.05);
   Graph_Graph1->GetYaxis()->SetTitleSize(0.035);
   Graph_Graph1->GetYaxis()->SetTitleFont(42);
   Graph_Graph1->GetZaxis()->SetLabelFont(42);
   Graph_Graph1->GetZaxis()->SetLabelSize(0.035);
   Graph_Graph1->GetZaxis()->SetTitleSize(0.035);
   Graph_Graph1->GetZaxis()->SetTitleFont(42);

   grtest->SetHistogram(Graph_Graph1);
   
   grtest->Draw("acp");
   
   // Predictions for pi-
   hadron = "pi-";
   double AUTpim[7];

   for (size_t i = 0; i < 7; i++)
   {
    double y = Q2/(S*x[i]);
    double DNN = (1.-y)/(1.+pow2(1-y));
    AUTpim[i] = DNN*torino10.FUTCollins(target, hadron, S, x[i], z, Q2, PhT)/torino10.FUU(target, hadron, S, x[i], z, Q2, PhT);
    //AUTpim[i] = torino10.FUTCollins(target, hadron, S, x[i], z, Q2, PhT)/torino10.FUU(target, hadron, S, x[i], z, Q2, PhT);
    //cout << AUTpim[i] << endl;
   }

   TGraph* grpim = new TGraph(7,x,AUTpim);
   grpim->Draw("cp");
   grpim->SetLineColor(4);
   grpim->SetLineWidth(3);
   
  
   TLine *line = new TLine(0.009,0,0.381,0);
   line->Draw();
   
   TPaveText *pt = new TPaveText(0.78,0.755,0.98,0.995,"brNDC");
   pt->SetName("description");
   pt->SetBorderSize(2);
   pt->SetFillColor(0);
   pt->SetTextAlign(12);
   pt->AddText("HERMES");
   pt->AddText("lP^{#uparrow}#rightarrow l' #pi^{#pm} X");
   pt->AddText(" #sqrt{s} = 7.25374 (GeV) ");
   pt->Draw();
   
   // pi+ Hermes data
   //DATA pip x
   TGraphErrors *gre = new TGraphErrors(7);
   gre->SetName("Graph");
   gre->SetTitle("Graph");
   gre->SetFillColor(4);
   gre->SetLineColor(4);
   gre->SetMarkerColor(4);
   gre->SetMarkerStyle(20);
   gre->SetMarkerSize(1.2);
   gre->SetPoint(0,0.036,0.0038);
   gre->SetPointError(0,0,0.008683893);
   gre->SetPoint(1,0.056,0.0033);
   gre->SetPointError(1,0,0.009035486);
   gre->SetPoint(2,0.076,0.0036);
   gre->SetPointError(2,0,0.008772685);
   gre->SetPoint(3,0.098,0.0071);
   gre->SetPointError(3,0,0.008077747);
   gre->SetPoint(4,0.133,0.0279);
   gre->SetPointError(4,0,0.005923681);
   gre->SetPoint(5,0.186,0.0289);
   gre->SetPointError(5,0,0.007003571);
   gre->SetPoint(6,0.275,0.015);
   gre->SetPointError(6,0,0.01116288);
   gre->Draw("p");

   // pi- Hermes data
   TGraphErrors *grem = new TGraphErrors(7);
   grem->SetName("Graph");
   grem->SetTitle("Graph");
   grem->SetFillColor(4);
   grem->SetLineColor(4);
   grem->SetMarkerColor(4);
   grem->SetMarkerStyle(24);
   grem->SetMarkerSize(1.2);
   grem->SetPoint(0,0.035,0.0007);
   grem->SetPointError(0,0,0.0077);
   grem->SetPoint(1,0.055, -0.0234);
   grem->SetPointError(1,0,0.0074);
   grem->SetPoint(2,0.076,-0.0034);
   grem->SetPointError(2,0,0.0089);
   grem->SetPoint(3,0.098,-0.0232);
   grem->SetPointError(3,0,0.0085);
   grem->SetPoint(4,0.133,-0.0127);
   grem->SetPointError(4,0,0.0084);
   grem->SetPoint(5,0.186, -0.0531);
   grem->SetPointError(5,0,0.0115);
   grem->SetPoint(6,0.275,-0.0752);
   grem->SetPointError(6,0, 0.0160);
   grem->Draw("p");

   nemo->Modified();
   resultcollins->cd();
   resultcollins->Modified();
   resultcollins->cd();
   resultcollins->SetSelected(resultcollins);
}


//========================================================= main
int main(int argc, char **argv)
{

  PARAMETERS P;
  PARTONCONTENT part;

  double Q2 = 2.41;
 
// --- Root
  TStyle *plain  = new TStyle("Plain","Plain Style (no colors/fill areas)");

  plain->SetCanvasBorderMode(0);
  plain->SetPadBorderMode(0);
  plain->SetPadColor(0);
  plain->SetCanvasColor(0);
  plain->SetTitleColor(0);
  plain->SetStatColor(0);
  plain->SetPalette(1);
   

  gROOT->SetStyle("Plain");

  TApplication* theApp = new TApplication("App", &argc, argv);

  DrawSiversDistribution( part,  P, Q2);

  DrawTransversityDistribution( part,  P, Q2);

  DrawCollinsDistribution( part,  P, Q2);

  HermesSivers();

  HermesCollins();

  theApp->Run(true);

  return 0;
}










