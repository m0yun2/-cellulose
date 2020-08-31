/**
 * @file vesta_ana.h
 * @brief セルロースの偏極度依存性のシュミレーション
 * @author Daisuke Miura
 * @date 2020.7.29
 */

 /**
  * @brief VESTA で作ったcrytak structure factor のデータをグラフにする
  * @details それぞれのQがの偏極度依存性や、TAIKANでのI vs Q のグラフを作る
  */
class VESTA{
private:
  string path_ = "/Users/daisuke/J-PARC/cellulose/cif/";
  string mresfile_ = "/Users/daisuke/J-PARC/cellulose/TAIKAN/TAIKAN_MD_res.txt";
  string bresfile_ = "/Users/daisuke/J-PARC/cellulose/TAIKAN/TAIKAN_BW_res.txt";
  vector<string> vstr_ = {"neg100.txt", "neg080.txt", "neg060.txt", "neg040.txt", "neg030.txt", "neg020.txt", "neg010.txt", "pos000.txt", "pos010.txt", "pos020.txt", "pos030.txt", "pos040.txt", "pos060.txt", "pos080.txt", "pos100.txt"};
  vector<int> vcol_ = {kViolet + 10, kViolet + 8, kViolet + 6, kViolet + 2, kViolet -9 , kViolet -5, kViolet, kBlack,
  kYellow +4, kYellow +3, kOrange +1, kOrange +8, kPink-9, kPink-1, kRed};
  pair<double, double> pqrange_ = make_pair(0.9, 7.5);
  const double qborder_ = 1.62;
public:
  VESTA(int mode);
  TGraph *LoadResolution(string filename);
  TGraph *LoadStructure(string filename);
  TGraph *LoadGausGraph(string filename);
  TGraph *AnalysisVESTA(int polindex, int mode);
  TGraph *SetGaus(TGraph *sg);
  void GraphOutput(string file, TGraph* g);

  TGraph *gresM_;
  TGraph *gresB_;
  vector <TGraph*> vgaus_, vasy_, vself_, vcross_, vmix_;
};

VESTA::VESTA(int mode){
  //Load Resolution Data
  gresM_ = LoadResolution(mresfile_);
  gresB_ = LoadResolution(bresfile_);
  //Prepare output file
  vector<string> vostr;
  for (int i = 0; i < vstr_.size(); i++){
    string copy = vstr_[i];
    string ostr = path_ + copy.erase(6, 4) + "_gaus.txt";
    vostr.emplace_back(ostr);
  }
  //Load Structure Data
  switch (mode){
    case 1:
    for (int i = 0; i < vstr_.size(); i++) {
      vgaus_.emplace_back(SetGaus(LoadStructure(path_ + vstr_[i])));
      vgaus_[i]->SetLineColor(vcol_[i]);
      vgaus_[i]->SetMarkerColor(vcol_[i]);
      vgaus_[i]->SetMarkerStyle(24);
      vgaus_[i]->SetMarkerSize(0.7);
      cout << i+1 << "/" << vstr_.size() <<"\t Loading Data ... " << vstr_[i] << endl;
      string copy = vstr_[i];
      string ostr = path_ + copy.erase(6, 4) + "_gaus.txt";
      GraphOutput(ostr, vgaus_[i]);
    }
    break;
    default:
    for (int i = 0; i < vostr.size(); i++){
      vgaus_.emplace_back(LoadGausGraph(vostr[i]));
      vgaus_[i]->SetLineColor(vcol_[i]);
      vgaus_[i]->SetMarkerColor(vcol_[i]);
      vgaus_[i]->SetMarkerStyle(24);
      vgaus_[i]->SetMarkerSize(0.7);
    }
    break;
  }

  //i < 7 : Number of data from negative data index to unpolarized data
  for (int i = 0; i < 7; i++){
    vasy_.emplace_back(AnalysisVESTA(i, 0));
    vmix_.emplace_back(AnalysisVESTA(i, 1));
    vself_.emplace_back(AnalysisVESTA(i, 2));
    vcross_.emplace_back(AnalysisVESTA(i, 3));
  }
}

/**
     * @brief TAIKAN のQ Resolution を読み込む
     * @param fiilname : TAIKAN のQ resokution のデータ(高田さんの論文より)
     * @return Qresolution vs. Q のTGraph
     * @detail 引数で指定するファイルをTGraph に直す。Middle bank とBackward bank でResoluton が変わる。
     */
TGraph *VESTA::LoadResolution(string filename){
  ifstream ifs;
  ifs.open (filename.c_str());
  if (!ifs) cout << "\033[33mFile name is wrong in LoadResolution \033[m" << endl << "Please check the following file. (path, filename.etc...) " << endl <<  filename.c_str() << endl;
  double q, res;
  TGraph *g = new TGraph();
  while (ifs >> res >> q){
    g->SetPoint(g->GetN(), q, res);
  }
  return g;
}

/**
     * @brief VESTA で作ったStructurefactor を読み込む
     * @param fiilname : VESTA でexport した.txt
     * @return I vs. Q のTGraph
     * @detail 引数で指定するファイルをTGraph に直す。ローレンツファクターはq^{-2}, 同じd space に出るピークは足し合わせを行っている
     */
TGraph *VESTA::LoadStructure(string filename){
  ifstream ifs;
  ifs.open (filename.c_str());
  if (!ifs) cout << "\033[33m File name is wrong in Load Structure \033[m" << endl  << "Please check the following file. (path, filename.etc...) " << endl <<  filename.c_str() << endl;
  int h, k, l;
  double d, freal, fimg, fabs, theta, ml, inten;
  vector<double> vd, vfabs, vml, vfreal, vfimg;
  vector<double> vq, vinten;
  while (ifs >> h >> k >> l >> d >> freal >> fimg >> fabs >> theta >> inten >> ml){
    if (isnan(inten)) continue;
    vd.emplace_back(d);
    vfreal.emplace_back(freal);
    vfimg.emplace_back(fimg);
    vfabs.emplace_back(fabs);
    vml.emplace_back(ml);
  }
  for (int i = 1; i < vd.size()+1; i++){
    if(vd[i-1] == vd[i]){
      double q = 2.0*M_PI/vd[i-1];
      vq.emplace_back(q);
      vinten.emplace_back(pow(vfabs[i-1],2)*2*vml[i-1]/pow(q,2));
      i+= 1;
    } else{
      double q = 2.0*M_PI/vd[i-1];
      vq.emplace_back(q);
      vinten.emplace_back(pow(vfabs[i-1],2)*vml[i-1]/pow(q,2));
    }
  }
  TGraph *g = new TGraph(vinten.size(), &vq[0], &vinten[0]);
  return g;
}

/**
     * @brief TGraph をガウス化させる
     * @param fiilname : 基本的にはLoadStructure の棒グラフ
     * @return I(ガウス) vs. Q のTGraph
     * @detail 引数で指定したTGraph をガウス幅をもたせたTGraph に直す。ガウスの幅はTAIKANのQresolution で与える、つまりLoadResolution のTGraph を使う
     */
TGraph* VESTA::SetGaus(TGraph *sg){
  vector <double> vx, vy;
  vector <TF1*> vf;
  for (int i = 0; i < sg->GetN(); i++){
    double xdata, ydata;
    sg->GetPoint(i, xdata, ydata);
    double qmin = xdata*0.9;
    double qmax = xdata*1.1;
    TF1 *f1 = new TF1("f1", "gaus(0)", qmin, qmax);
    double qres = 0.0;
    if (xdata < qborder_) qres = gresM_->Eval(xdata)/2.35;
    else qres = gresB_->Eval(xdata)/2.35;
    f1->FixParameter(0, ydata);
    f1->FixParameter(1, xdata);
    f1->FixParameter(2, xdata*qres);
    vf.emplace_back(f1);
  }
  TGraph *g = new TGraph();
  double plot = pqrange_.first;
  while ( plot < pqrange_.second){
    double sumdata = 0.0;
    for (int i = 0; i < vf.size(); i++) sumdata += vf[i]->Eval(plot);
    g->SetPoint(g->GetN(), plot, sumdata);
    plot += 0.001;
  }
  return g;
}

/**
     * @brief VESTA の棒グラフをガウスにしたデータが保存されている.txt を読み込んでTGraph にする。
     * @param fiilname : TGraph にしたい _gaus.txt
     * @return I vs Q のTGraph
     * @detail コンパイルするたびにLoadStructure を呼ぶと非常に動作が遅い, 一度コンパイルした際にLoadStructure のデータを_gaus.txt に保存しておく。以降は_gaus.txt を呼び出すことで処理を軽くする。
     */
TGraph* VESTA::LoadGausGraph(string filename){
  ifstream ifs;
  ifs.open (filename.c_str());
  if (!ifs) cout << "\033[33m File name is wrong in LoadGausGraph \033[m" << endl  << "Please check the following file. (path, filename.etc...) " << endl <<  filename.c_str() << endl;
  double qdata, idata;
  TGraph *g = new TGraph();
  while (ifs >> qdata >> idata) g->SetPoint(g->GetN(), qdata, idata);
  return g;
}

/**
     * @brief Asymmetry, Self Term, Cross Term, Sel + Cross term を求める
     * @param npolindex : 解析したい偏極度に相当するvstr_ のIndex。Negative のものを選ぶ。
     * @param mode : 解析モード, 0 : Asymmetry, 1 : mix, 2 : Self, 3 : cross
     * @return I vs Q のTGraph
     * @detail
     */
TGraph* VESTA::AnalysisVESTA(int npolindex, int mode){
  string copy = vstr_[npolindex];
  copy.erase(0, 3);
  copy.erase(3, 4);
  double pol = stod(copy);
  pol /= 100.0;
  int ppolindex = vstr_.size() -1 -npolindex;
  int upolindex = (vstr_.size() -1)/2 ;
  vector <double> vni, vpi, vui, vq;
  TGraph *g = new TGraph();
  for (int i = 0; i < vgaus_[npolindex]->GetN(); i++){
    double q, ni, pi, ui;
    vgaus_[npolindex] -> GetPoint(i, q, ni);
    vgaus_[ppolindex] -> GetPoint(i, q, pi);
    vgaus_[upolindex] -> GetPoint(i, q, ui);
    double ana = 0.0;
    double k = 3.89;
    switch (mode){
      case 0: //Asymmetry
      ni += 1000;
      pi += 1000;
      ana = (ni - pi)/(ni + pi);
      break;
      case 1: //mix
      ana = (ni - pi)/2.0/pol/2.0/k;
      break;
      case 2: //self
      ana = (ni - pi - 2*ui)/2.0/pow(pol, 2)/pow(k,2);
      break;
      case 3: //cross
      ana = ((ni - pi)/2.0/pol/2.0/k) - ((ni - pi - 2*ui)/2.0/pow(pol, 2)/pow(k,2));
    }
    g->SetPoint(g->GetN(), q, ana);
  }
  return g;
}


/**
     * @brief TGraph のデータを.txt に
     * @param file: 保存するファイル名
     * @param g : .txt にしたいTGraph
     * @return Q "\t" I の.txt
     * @detail タブ区切りでQ I の.txt にする
     */
void VESTA::GraphOutput(string file, TGraph* g){
  ofstream ofs (file.c_str());
  for (int i = 0; i < g->GetN(); i++){
    double xdata, ydata;
    g->GetPoint(i, xdata, ydata);
    ofs << xdata << "\t" << ydata << endl;
  }
  ofs.close();
}
