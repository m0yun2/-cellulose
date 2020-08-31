/**
 * @file CAMVAS.h
 * @brief Modify of TCanvas of ROOT
 * @author Daisuke Miura
 * @date 2020.8.1
 */
 
TH1 *MakeFrame (const string title, const string xtitle, const string ytitle, const double xmin, const double ymin, const double xmax, const double ymax, int titleoption = 0, int gridoption = 0){
  TH1 *frame = gPad-> DrawFrame(xmin, ymin, xmax, ymax);
  gPad->SetTickx(1);
  gPad->SetTicky(1);
  TLatex *titlex = new TLatex();
  titlex->SetTextFont(42);
  titlex->SetTextAlign(22);
  TLatex *titley = new TLatex();
  titley->SetTextFont(42);
  titley->SetTextAlign(22);
  titley->SetTextAngle(90);

  titlex->SetTextSize( 0.05 );
  titley->SetTextSize( 0.05 );
  frame->GetXaxis()->SetLabelSize( 0.05 );
  frame->GetYaxis()->SetLabelSize( 0.05 );

  switch (titleoption) {
    default: //normal
      titlex->SetTextSize( 0.04 );
      titley->SetTextSize( 0.04 );
      titlex->DrawLatex( (xmax -xmin)/2.0 + xmin , ymin - (ymax-ymin)*0.08, xtitle.c_str());
      titley->DrawLatex( xmin - (xmax -xmin)*0.08 , (ymax - ymin)/2.0 + ymin , ytitle.c_str());
    break;
    case 1:
      titlex->SetTextSize( 0.08 );
      titlex->DrawLatex( (xmax -xmin)/2.0 + xmin , ymin- (ymax-ymin)*0.1, xtitle.c_str());
      break;
    case 2:
      titley->SetTextSize( 0.08 );
      titley->DrawLatex( xmin -0.15 , (ymax - ymin)/2.0 + ymin , ytitle.c_str());
      break;
    case 3:
      break;
    case 4: //Divided Canvas mode
      titlex->SetTextSize( 0.05 );
      titley->SetTextSize( 0.05 );
      titlex->DrawLatex( (xmax -xmin)/2.0 + xmin , ymin- (ymax-ymin)*0.1, xtitle.c_str());
      titley->DrawLatex( xmin -0.15 , (ymax - ymin)/2.0 + ymin , ytitle.c_str());
      frame->GetYaxis()->SetLabelSize(0.05);
      frame->GetXaxis()->SetLabelSize(0.05);
      break;
    case 5:
      titlex->SetTextSize( 0.05 );
      titley->SetTextSize( 0.05 );
      gPad->SetTicky(0);
      frame->GetYaxis()->SetLabelOffset(1000);
      frame->GetYaxis()->SetTickLength(0.0);
      titlex->DrawLatex( (xmax -xmin)/2.0 + xmin , ymin - (ymax-ymin)*0.1, xtitle.c_str());
      titley->DrawLatex( xmin - (xmax -xmin)*0.05 , (ymax - ymin)/2.0 + ymin , ytitle.c_str());
      break;
  }

  switch (gridoption) {
    case 0:
      break;
    case 1:
      gPad->SetGridx(1);
      break;
    case 2:
      gPad->SetGridx(1);
      gPad->SetGridy(1);
      break;
  }
  return frame;
}

int SaveCanvas(string savefile){
  string pngfile = savefile + ".png";
  string pdffile = savefile + ".pdf";
  string psfile = savefile + ".ps";
  gPad->Print(pngfile.c_str());
  string convert = "convert " + pngfile + " " + pdffile;
  string convertps = "convert " + pngfile + " " + psfile;
  system(convert.c_str());
  system(convertps.c_str());
  return 0;
}
