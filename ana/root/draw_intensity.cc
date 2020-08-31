/**
 * @file draw_intensity.cc
 * @brief セルロースの偏極度依存性のシュミレーション
 * @author Daisuke Miura
 * @date 2020.7.29
 */

#include "VESTA.h"
#include "CANVAS.h"

void GraphOutput(string file, TGraph* g){
  ofstream ofs (file.c_str());
  for (int i = 0; i < g->GetN(); i++){
    double xdata, ydata;
    g->GetPoint(i, xdata, ydata);
    ofs << xdata << "\t" << ydata << endl;
  }
  ofs.close();
}

void draw_intensity(){
  auto * c = new TCanvas("c", "c", 1400, 1000);
  double mleft = 0.2;
  double mright = 0.05;
  double mtop = 0.1;
  double mbottom = 0.1;
  c -> SetMargin (mleft, mright, mbottom, mtop);

  TH1 *frame = MakeFrame("", "#it{Q} /#AA^{-1}", "Intensity /a.u.", 0.98, -1.0, 3.0, 25000);
  VESTA ves(0);
  for (int i = 0; i < ves.vmix_.size(); i++) {
    string strmix[7] = {"mix_100.dat", "mix_80.dat", "mix_60.dat", "mix_40.dat", "mix_30.dat", "mix_20.dat", "mix_010.dat"};
    string strself[7] = {"self_100.dat", "self_80.dat", "self_60.dat", "self_40.dat", "self_30.dat", "self_20.dat", "self_010.dat"};
    string strcross[7] = {"cross_100.dat", "cross_80.dat", "cross_60.dat", "cross_40.dat", "cross_30.dat", "cross_20.dat", "cross_010.dat"};
    GraphOutput(strmix[i].c_str(), ves.vmix_[i]);
    GraphOutput(strself[i].c_str(), ves.vself_[i]);
    GraphOutput(strcross[i].c_str(), ves.vcross_[i]);
  }
  SaveCanvas("result/intensity_pol_dep");
}
