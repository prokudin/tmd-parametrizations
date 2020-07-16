#include "sivers.h"

using namespace std;

//========================================================= Sivers effect...x dependence
double sivers_x_dependence(double x, double a, double b, double n)
{
  double sivers_x = n * pow(x,a) * pow((1.-x),b) * pow((a+b),(a+b)) /(pow(a,a) * pow(b,b));

  return sivers_x;
}

//========================================================= Sivers effect...
void sivers(struct PARTONCONTENT& partcontent, PARAMETERS& Params, double x ) {
  // Returns the parton content for Sivers effect...
 
  partcontent.up  *= 2.*
    sivers_x_dependence(x, Params.Sivers.parameters.a_up, 
			Params.Sivers.parameters.b_up, 
			Params.Sivers.parameters.n_up);

  partcontent.down  *= 2.*
    sivers_x_dependence(x, Params.Sivers.parameters.a_down, 
			Params.Sivers.parameters.b_down, 
			Params.Sivers.parameters.n_down); // 

  partcontent.anti_up *= 2.*
    sivers_x_dependence(x, Params.Sivers.parameters.a_anti_up, 
			Params.Sivers.parameters.b_anti_up, 
			Params.Sivers.parameters.n_anti_up); //

  partcontent.anti_down *= 2.* 
    sivers_x_dependence(x, Params.Sivers.parameters.a_anti_down, 
			Params.Sivers.parameters.b_anti_down, 
			Params.Sivers.parameters.n_anti_down); //  


  partcontent.strange *= 2.* 
    sivers_x_dependence(x, Params.Sivers.parameters.a_strange, 
			Params.Sivers.parameters.b_strange, 
			Params.Sivers.parameters.n_strange); // 


  partcontent.anti_strange *= 2.* 
    sivers_x_dependence(x, Params.Sivers.parameters.a_anti_strange, 
			Params.Sivers.parameters.b_anti_strange, 
			Params.Sivers.parameters.n_anti_strange); //   

  partcontent.charm     = 0.;
  partcontent.anti_charm  = 0.;
  partcontent.bottom     = 0.;
  partcontent.anti_bottom  = 0.;
  partcontent.top     = 0.;
  partcontent.anti_top  = 0.;

}


//========================================================= Sivers partcontent for the Proton
void SiversDistribution( struct PARTONCONTENT& partcontent,  PARAMETERS& Params, double x, double Q2)
  // Returns the parton content at x,Q2
{
  unpolarised(partcontent, x, Q2);
  sivers(partcontent, Params, x);
}


//========================================================= Sivers partcontent for the all targets
void SiversDistribution( std::string target, struct PARTONCONTENT& partcontent, PARAMETERS& Params, double x, double Q2)
  // Returns the parton content at x,Q2
{
  SiversDistribution( partcontent,  Params, x, Q2);
  
  switch( getTarget(target) ){
  case PROTON:
    break;
  case ANTIPROTON:
    Antiproton(partcontent);
    break;
  case DEUTERON:
    Deutron(partcontent);
    break;
  case NEUTRON:
    Neutron(partcontent);
    break;
  default:
    break;
  }
}


//========================================================= Sivers partcontent KT
void SiversDistributionKt( struct PARTONCONTENT& partcontent,  PARAMETERS& Params, double x, double kt, double Q2)
  // Returns the parton content at x,Q2
{
  unpolarised(partcontent, x, Q2);
  sivers(partcontent, Params, x);

  double kt2_average = Params.GetKt2Average();

  double unpolarised_kt = 1./(M_PI * kt2_average) * exp( -kt*kt / kt2_average);
 
  partcontent.up  *= sqrt( 2.*exp(1.)/ Params.Sivers.parameters.m2_up )  * kt * exp( -kt*kt/Params.Sivers.parameters.m2_up) * unpolarised_kt; // Eq (17.) hep-ph/0507181 

  partcontent.down  *= sqrt( 2.*exp(1.)/ Params.Sivers.parameters.m2_down )  * kt * exp( -kt*kt/Params.Sivers.parameters.m2_down) * unpolarised_kt; // 

  partcontent.anti_up *= sqrt( 2.*exp(1.)/ Params.Sivers.parameters.m2_anti_up )  * kt * exp( -kt*kt/Params.Sivers.parameters.m2_anti_up) * unpolarised_kt; //

  partcontent.anti_down *= sqrt( 2.*exp(1.)/ Params.Sivers.parameters.m2_anti_down )  * kt * exp( -kt*kt/Params.Sivers.parameters.m2_anti_down) * unpolarised_kt; //  


  partcontent.strange *= sqrt( 2.*exp(1.)/ Params.Sivers.parameters.m2_strange )  * kt * exp( -kt*kt/Params.Sivers.parameters.m2_strange) * unpolarised_kt; // 


  partcontent.anti_strange *= sqrt( 2.*exp(1.)/ Params.Sivers.parameters.m2_anti_strange )  * kt * exp( -kt*kt/Params.Sivers.parameters.m2_anti_strange) * unpolarised_kt; //   

  partcontent.charm     = 0.;
  partcontent.anti_charm  = 0.;
  partcontent.bottom     = 0.;
  partcontent.anti_bottom  = 0.;
  partcontent.top     = 0.;
  partcontent.anti_top  = 0.;
}


//========================================================= Sivers partcontent First moment
void SiversDistributionFirstMoment( struct PARTONCONTENT& partcontent,  PARAMETERS& Params, double x, double Q2)
// Returns the parton content at x,Q2
{
  unpolarised(partcontent, x, Q2);
  sivers(partcontent, Params, x);

  //double integrand = 1;
  double mpr = 0.93827203 ; // the proton mass
  double kt2_average = Params.GetKt2Average();

  partcontent.up  *= sqrt( 2.*exp(1.) )/( 4.*mpr ) * pow(Params.Sivers.parameters.m2_up, 3./2.)*kt2_average/
    pow(Params.Sivers.parameters.m2_up +kt2_average, 2); // Eq (22.) hep-ph/0507181 

  partcontent.down  *= sqrt( 2.*exp(1.) )/( 4.*mpr ) * pow(Params.Sivers.parameters.m2_down, 3./2.)*kt2_average/
    pow(Params.Sivers.parameters.m2_down +kt2_average, 2); // 

  partcontent.anti_up *= sqrt( 2.*exp(1.) )/( 4.*mpr ) * pow(Params.Sivers.parameters.m2_up, 3./2.)*kt2_average/
    pow(Params.Sivers.parameters.m2_up +kt2_average, 2); //

  partcontent.anti_down *= sqrt( 2.*exp(1.) )/( 4.*mpr ) * pow(Params.Sivers.parameters.m2_anti_down, 3./2.)*kt2_average/
    pow(Params.Sivers.parameters.m2_anti_down +kt2_average, 2); //  


  partcontent.strange *= sqrt( 2.*exp(1.) )/( 4.*mpr ) * pow(Params.Sivers.parameters.m2_strange, 3./2.)*kt2_average/
    pow(Params.Sivers.parameters.m2_strange +kt2_average, 2); // 


  partcontent.anti_strange *= sqrt( 2.*exp(1.) )/( 4.*mpr ) * pow(Params.Sivers.parameters.m2_anti_strange, 3./2.)*kt2_average/
    pow(Params.Sivers.parameters.m2_anti_strange +kt2_average, 2); //   

  partcontent.charm     = 0.;
  partcontent.anti_charm  = 0.;
  partcontent.bottom     = 0.;
  partcontent.anti_bottom  = 0.;
  partcontent.top     = 0.;
  partcontent.anti_top  = 0.;
}


//========================================================= Draw Sivers partcontent
TCanvas* DrawSiversDistribution( struct PARTONCONTENT& partcontent,  PARAMETERS& Params, double Q2){

  TCanvas* canvas1 = new TCanvas("sivers","sivers");
  canvas1->SetFillColor(0);
  canvas1->SetFrameFillStyle(0);
  canvas1->SetFillStyle(0);

  //  canvas1->Size(200,300);
  TLatex *t = new TLatex();
  t->SetTextAngle(90);
  t->SetTextSize(0.05);

  t->DrawLatex(0.04,0.74,"x#Delta^{N} f^{(1)}(x)");
  t->DrawLatex(0.10,0.92,"u");
  t->DrawLatex(0.10,0.8,"d");
  t->DrawLatex(0.10,0.62,"#bar{u}");
  t->DrawLatex(0.10,0.48,"#bar{d}");
  t->DrawLatex(0.10,0.34,"s");
  t->DrawLatex(0.10,0.2,"#bar{s}");
  t->DrawLatex(0.51,0.74,"x#Delta^{N} f(x, k_{#perp}  )");
  t->DrawLatex(0.560,0.92,"u");
  t->DrawLatex(0.560,0.8,"d");
  t->DrawLatex(0.560,0.62,"#bar{u}");
  t->DrawLatex(0.560,0.48,"#bar{d}");
  t->DrawLatex(0.560,0.34,"s");
  t->DrawLatex(0.560,0.2,"#bar{s}");

  t->SetTextAngle(0);
  t->SetTextSize(0.05);
  t->DrawLatex(0.45,0.05,"x");
  t->DrawLatex(0.8,0.05,"k_{#perp}   (GeV)");

  t->DrawLatex(0.1,0.01,"Eur.Phys.J.A39:89-100,2009");

  canvas1->Draw();

  TPad* tpad = new TPad("nemo"," ",0.1,0.1  ,1.,1.,0,0,0);
  tpad->Range(0,0,10,10);
  tpad->SetFillColor(0);
  tpad->SetFrameFillStyle(0);
  tpad->SetFillStyle(0);
  gStyle->SetPadBorderMode(0);
  gStyle->SetFrameLineWidth(0);  
  gStyle->SetLabelSize(0.14,"xy");
  gStyle->SetLabelOffset(0.06,"x");
  gStyle->SetLabelOffset(0.02,"y");
  gStyle->SetNdivisions(4,"y");
  gStyle->SetNdivisions(8,"x");
  Float_t small = 0;

  gPad->SetFillColor(0);

  tpad->Divide(2,6,small,small);
  TPad *subpad1 = (TPad*)tpad->GetPad(1);
  TPad *subpad2 = (TPad*)tpad->GetPad(2);
  TPad *subpad3 = (TPad*)tpad->GetPad(3);
  TPad *subpad4 = (TPad*)tpad->GetPad(4);
  TPad *subpad5 = (TPad*)tpad->GetPad(5);
  TPad *subpad6 = (TPad*)tpad->GetPad(6);
  TPad *subpad7 = (TPad*)tpad->GetPad(7);
  TPad *subpad8 = (TPad*)tpad->GetPad(8);
  TPad *subpad9 = (TPad*)tpad->GetPad(9);
  TPad *subpad10= (TPad*)tpad->GetPad(10);
  TPad *subpad11= (TPad*)tpad->GetPad(11);
  TPad *subpad12= (TPad*)tpad->GetPad(12);

     
  subpad1->SetFillStyle(0);
  subpad1->SetFrameFillStyle(0);
  subpad2->SetFillStyle(0);
  subpad2->SetFrameFillStyle(0);
  subpad3->SetFillStyle(0);
  subpad3->SetFrameFillStyle(0);
  subpad4->SetFillStyle(0);
  subpad4->SetFrameFillStyle(0);
  subpad5->SetFillStyle(0);
  subpad5->SetFrameFillStyle(0);
  subpad6->SetFillStyle(0);
  subpad6->SetFrameFillStyle(0);
  subpad7->SetFillStyle(0);
  subpad7->SetFrameFillStyle(0);
  subpad8->SetFillStyle(0);
  subpad8->SetFrameFillStyle(0);
  subpad9->SetFillStyle(0);
  subpad9->SetFrameFillStyle(0);
  subpad10->SetFillStyle(0);
  subpad10->SetFrameFillStyle(0);
  subpad11->SetFillStyle(0);
  subpad11->SetFrameFillStyle(0);
  subpad12->SetFillStyle(0);
  subpad12->SetFrameFillStyle(0);

  tpad->Draw();


  int    npoints = 30;

  double x_min = 0.001;
  double x_max = 0.95;
  double x_step  = (log( x_max ) - log( x_min ))/float(npoints - 1);

  double kt_min = 0.001;
  double kt_max = 1.3;
  double kt_step  = (kt_max  - kt_min )/float(npoints - 1);

  double xx[npoints];
  double ktt[npoints];
  double yy[6][npoints];
  double yy1[6][npoints];

  double x_draw = 0.1;
  
  for(int i = 0; i < npoints; ++i){
    if(i == 0){
	xx[i] = x_min;
	ktt[i] = kt_min;
    }
    else {
      xx[i] = exp( log( x_min ) + float(i) * x_step );
      ktt[i] = ktt[i-1] + kt_step;
    }
      
      SiversDistributionFirstMoment( partcontent,  Params, xx[i], Q2);
      yy[0][i] = xx[i] * partcontent.up;
      yy[1][i] = xx[i] * partcontent.down;
      yy[2][i] = xx[i] * partcontent.anti_up;
      yy[3][i] = xx[i] * partcontent.anti_down;
      yy[4][i] = xx[i] * partcontent.strange;
      yy[5][i] = xx[i] * partcontent.anti_strange;

      SiversDistributionKt( partcontent,  Params, x_draw, ktt[i], Q2);
      yy1[0][i] = x_draw * partcontent.up;
      yy1[1][i] = x_draw * partcontent.down;
      yy1[2][i] = x_draw * partcontent.anti_up;
      yy1[3][i] = x_draw * partcontent.anti_down;
      yy1[4][i] = x_draw * partcontent.strange;
      yy1[5][i] = x_draw * partcontent.anti_strange;


  }


  double right_mar = 0.17; // Normal aspect ration
  double left_mar = 0.12;
  double right_mar1 = 0.1; 

  TGraph* resultgraph;

  subpad1->cd();
  subpad1->SetLeftMargin(left_mar);
  subpad1->SetRightMargin(right_mar);
  subpad1->SetBottomMargin(small);

  //up
  resultgraph = new TGraph(npoints, xx, yy[0]);
  resultgraph->SetLineColor(2);
  resultgraph->SetLineWidth(2);
  resultgraph->SetTitle("");
  subpad1-> SetLogx();
  resultgraph->Draw("ACP");

  subpad2->cd();
  subpad2->SetLeftMargin(left_mar);
  subpad2->SetRightMargin(right_mar1);
  subpad2->SetBottomMargin(small);

  //up
  resultgraph = new TGraph(npoints, ktt, yy1[0]);
  resultgraph->SetLineColor(2);
  resultgraph->SetLineWidth(2);
  resultgraph->SetTitle("");
  resultgraph->Draw("ACP");

  subpad3->cd();
  subpad3->SetLeftMargin(left_mar);
  subpad3->SetRightMargin(right_mar);
  subpad3->SetTopMargin(small);
  subpad3->SetBottomMargin(small);

  //down
  resultgraph = new TGraph(npoints, xx, yy[1]);
  resultgraph->SetLineColor(2);
  resultgraph->SetLineWidth(2);
  resultgraph->SetTitle("");
  subpad3-> SetLogx();
  resultgraph->Draw("ACP");

  subpad4->cd();
  subpad4->SetLeftMargin(left_mar);
  subpad4->SetRightMargin(right_mar1);
  subpad4->SetTopMargin(small);
  subpad4->SetBottomMargin(small);

  //down
  resultgraph = new TGraph(npoints, ktt, yy1[1]);
  resultgraph->SetLineColor(2);
  resultgraph->SetLineWidth(2);
  resultgraph->SetTitle("");
  resultgraph->Draw("ACP");

  subpad5->cd();
  subpad5->SetLeftMargin(left_mar);
  subpad5->SetRightMargin(right_mar);
  subpad5->SetTopMargin(small);
  subpad5->SetBottomMargin(small);

  //anti up
  resultgraph = new TGraph(npoints, xx, yy[2]);
  resultgraph->SetLineColor(2);
  resultgraph->SetLineWidth(2);
  resultgraph->SetTitle("");
  resultgraph->SetMinimum(0.);
  resultgraph->SetMaximum(0.001);
  subpad5-> SetLogx();
  resultgraph->Draw("ACP");

  subpad6->cd();
  subpad6->SetLeftMargin(left_mar);
  subpad6->SetRightMargin(right_mar1);
  subpad6->SetTopMargin(small);
  subpad6->SetBottomMargin(small);

  //anti up
  resultgraph = new TGraph(npoints, ktt, yy1[2]);
  resultgraph->SetLineColor(2);
  resultgraph->SetLineWidth(2);
  resultgraph->SetTitle("");
  resultgraph->Draw("ACP");


  subpad7->cd();
  subpad7->SetLeftMargin(left_mar);
  subpad7->SetRightMargin(right_mar);
  subpad7->SetTopMargin(small);
  subpad7->SetBottomMargin(small);

  //anti down
  resultgraph = new TGraph(npoints, xx, yy[3]);
  resultgraph->SetLineColor(2);
  resultgraph->SetLineWidth(2);
  resultgraph->SetTitle("");
  subpad7-> SetLogx();
  resultgraph->Draw("ACP");

  subpad8->cd();
  subpad8->SetLeftMargin(left_mar);
  subpad8->SetRightMargin(right_mar1);
  subpad8->SetTopMargin(small);
  subpad8->SetBottomMargin(small);

  //anti down
  resultgraph = new TGraph(npoints, ktt, yy1[3]);
  resultgraph->SetLineColor(2);
  resultgraph->SetLineWidth(2);
  resultgraph->SetTitle("");
  resultgraph->Draw("ACP");

  subpad9->cd();
  subpad9->SetLeftMargin(left_mar);
  subpad9->SetRightMargin(right_mar);
  subpad9->SetTopMargin(small);
  subpad9->SetBottomMargin(small);

  //strange
  resultgraph = new TGraph(npoints, xx, yy[4]);
  resultgraph->SetLineColor(2);
  resultgraph->SetLineWidth(2);
  resultgraph->SetTitle("");
  subpad9-> SetLogx();
  resultgraph->Draw("ACP");


  subpad10->cd();
  subpad10->SetLeftMargin(left_mar);
  subpad10->SetRightMargin(right_mar1);
  subpad10->SetTopMargin(small);
  subpad10->SetBottomMargin(small);

  //strange
  resultgraph = new TGraph(npoints, ktt, yy1[4]);
  resultgraph->SetLineColor(2);
  resultgraph->SetLineWidth(2);
  resultgraph->SetTitle("");
  resultgraph->Draw("ACP");

  subpad11->cd();
  subpad11->SetLeftMargin(left_mar);
  subpad11->SetRightMargin(right_mar);
  subpad11->SetTopMargin(small);
  subpad11->SetBottomMargin(0.24);

  //anti strange
  resultgraph = new TGraph(npoints, xx, yy[5]);
  resultgraph->SetLineColor(2);
  resultgraph->SetLineWidth(2);
  resultgraph->SetTitle("");
  subpad11-> SetLogx();
  resultgraph->Draw("ACP");

  subpad12->cd();
  subpad12->SetLeftMargin(left_mar);
  subpad12->SetRightMargin(right_mar1);
  subpad12->SetTopMargin(small);
  subpad12->SetBottomMargin(0.24);

  //anti strange
  resultgraph = new TGraph(npoints, ktt, yy1[5]);
  resultgraph->SetLineColor(2);
  resultgraph->SetLineWidth(2);
  resultgraph->SetTitle("");
  resultgraph->Draw("ACP");

  return canvas1;
}