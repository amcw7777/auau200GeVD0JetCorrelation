#include "d0jet_corr_plotter.cc"
void runD0jetCorrPlotter()
{
  D0jetCorrPlotter *plotter = new D0jetCorrPlotter();
  plotter->init("root-files/noPx.root");
  pair<int,int> inputCentralityBin(1,10);
  int inputRebinFactor = 100;
  plotter->getCorrelation(inputCentralityBin,inputRebinFactor);
  plotter->plotCorrelation();
  TCanvas *c = new TCanvas();
  c->cd();
  TH1D *candFar = plotter->getCandCorrelationFar();
  TH1D *bkgFar = plotter->getBkgCorrelationFar();
  candFar->Draw();
  bkgFar->Draw("same");
  candFar->SetLineColor(2);
  TLegend *leg = new TLegend(0.1,0.5,0.5,0.9);
  leg->AddEntry(bkgFar,"Bkg correlation");
  leg->AddEntry(candFar,"candidate ecorrelation");
  leg->Draw("same");

  // plotter->getPxCut(0.1);
}

