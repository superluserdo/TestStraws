// ---------------------------------------------------------------------
// Version of TestStraws for use with data taken during the QA tests
// of each module.
// ---------------------------------------------------------------------
#include "TROOT.h"
#include "TChain.h"
#include "TFile.h"
#include <TTree.h>

#include "TH1.h"
#include "TH2.h"
#include "TProfile.h"
#include "TStyle.h"
#include "TVector3.h"
#include "TCanvas.h"
#include "TText.h"
#include "TSystem.h"
#include "TPaveLabel.h"
#include "TLine.h"
#include "TLegend.h"
#include "TMath.h"
#include "TRandom3.h"
#include "TLeaf.h"
//
#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>
using namespace std;

#include <cmath>
#include <TF1.h>
#include <TF2.h>
#include <TGraph.h>
#include <TRandom3.h>
#include <TColor.h>
// --------------------------------------------------------------
// Event
// TriggerTime
// HitTime
// Station
// Module
// View
// Layer
// Wire
//
struct StrawData {
       unsigned int run;
       unsigned int event;
       double triggerTime;
       double hitTime;
       int station;
       int module;
       int iview;
       int layer;
       int wire;
};
struct ScintData {
       unsigned int run;
       unsigned int event;
       double triggerTime;
       double hitTime;
       int station;
       int module;
       int iview;
       int layer;
       int wire;
};
// --------------------------------------------------------------
// Now my stuff from TestStraws.h:
// --------------------------------------------------------------
      const int Nmax = 2; 
      TRandom3* Rnd = new TRandom3;
      int IcallRndm = 0;
// --------------------------------------------------------------
// Declare my global event variables:
// --------------------------------------------------------------
// First arrays to hold details of 100 sequential hits off the
// ProfTree file:
// --------------------------------------------------------------
      int Nhits;
      int Iruns[100];
      int Ievts[100];
      double TrigTimes[100];
      int Layers[100];
      int Istraw[100];
      int Module[100];
      double RawTimes[100];
      int Iviews[100];
// --------------------------------------------------------------
// Now some arrays to contain everything we need for a particular
// "event". From these we can make a decision whether to proceed
// with track reconstruction.
// --------------------------------------------------------------
      int NhitsInEvent;
      int ModEvent[100];
      int LayEvent[100];
      int IstEvent[100];
      float RawTevent[100];
      int Nhmods[3];
      int Nhlays[4];
// --------------------------------------------------------------
// Now all the arrays used in the actual track finding.
// --------------------------------------------------------------
      float Zpos;
      float Xzstraw[6][32][4][3];
      float Diam;
      float Ylength;
      float Size[4];
      float Xzstrau[6][32][4][3];
      float Xzstrav[6][32][4][3];
      float Xzaty[2][32][4][3];
      int   Iexist[32][4][3];
      float Pedest[32][4][3];
      float T0data;
      float Tofl;
      float Yped;
      float Hits[32][4][3];
      float Dtimes[32][4][3];
      float Tof[32][4][3];
      float Rtimes[32][4][3];
      int   Mask[32][4][3];
      float Yhits[32][4][3];
      float Zxhit[2][32][4][3];
      float Trkinf[60][3];
      float Trkin2[60][3];
      float Trkin3[60][3];
//
      float Thecol;
      float Sthpla;
      float Cthpla;
      int Nposte;
      float Xypost[2][4][3];
      int Nspose[4][3];
      int Nlpose[4][3];

// --------------------------------------------------------------
   struct StrHist_t {
          TH1F *fWire[Nmax];
          TH1F *fModule[Nmax];
          TH1F *fLayer[Nmax];
          TH1F *fdiff12[Nmax];
          TH1F *fdiff34[Nmax];
          TH1F *fDtimes[Nmax];
          TH1F *fDdists[Nmax];
          TH1F *fTrkMod[Nmax];
          TH1F *fYbest[Nmax];
          TH1F *fGrfit[Nmax];
          TH1F *fPhizxd[Nmax];
          TH2F *fZxte[Nmax];
          TH1F *fRest1[Nmax];
          TH1F *fRest2[Nmax];
          TH1F *fRest3[Nmax];
          TH1F *fRest4[Nmax];
          TH1F *fRest[Nmax];

   };
   StrHist_t          StrHist;
//
   void DefineColours();
   void BookHistograms();
   void BookStrHistograms();
//
   void Strini();
   void Bkzero(int imode);
   void WriteStraws();
   void StudyStraws(int iev);
   void PlotHistograms();
   void PlotStrHistograms(int Ipage);
   void WriteHistograms();
   void WriteStrHistograms(int Ipage);
   int  mod(int n1, int n2);
   float BkRndm();
//
   void Vfill(float* Array, int Len, float Value);
   void Vzero(float* Array, int Len);
   void Bkxy2u(float xx,float yy, float &xu, float &yu);
   void Bkxy2v(float xx,float yy, float &xv, float &yv);
   void Bku2xy(float xu,float yu, float &xx, float &yy);
   void Bkv2xy(float xv,float yv, float &xx, float &yy);
//
   void Strcht(float zz1,float xx1,float zz2,float xx2,int Nstraw,int Layer,int Nmod,float &Zpt, float &Xpt);
   void Strcor(float zz1,float xx1,float zz2,float xx2,int Nst,int Nlay,int Nmod,float &Angle, int &Iflag);
   void Strelm(int Nmod,int Nst1,int Layer1,int Nst2,int Layer2,int Icand,
               float &zz1,float &xx1,float &zz2,float &xx2,int &Iflag);
   void Strflr(int Nstraw, int Layer, int Nmod, float zpt, float xpt);
   void Strclp(float Gr,float c1,float Xpt,float Ypt,float &Xinter,float &Yinter);
   void Sttime1(float Time,float &Dist);
   void Sttime2(float Time,float &Dist);
   void Sttime3(float Time,float &Dist);
   void Linfit(float* Xpts,float* Ypts,int Npts,float &Grad,float &Cint);
//
   int  Bknint(float var);
   float Bkatan(float py , float px);
//
   void Strft4(int Nmod);
   void Strft3(int Nmod);
   void Strte2(int Nmod,int Ist1,int Ist2,int Ist3,int Ist4);
   void Strte3(int Nmod,int Ist1,int Ist2,int Ist3,int Ist4); 
   void Strana();
   void Strevt(int &Goodt0);
//
   double FitBack(double* x, double* par);
   double Lorentzian(double* x, double* par);
   double Gaussian1(double* x, double* par);
   double Gaussian2(double* x, double* par);
   double FitFunction(double* x, double* par);
   double FitFunction2(double* x, double* par);
// -----------------------------------------------------------------------------------
    int main() {
// -----------------------------------------------------------------------------------
// load data file
// --------------
//  TFile *inFile = new TFile("/hepstore/g2share/ProfTree/Lab3TreeDump_00294_00297_00298.root");
//  TFile *inFile = new TFile("/hepstore/btk/ProfTree/run01805_ArEth_1600V_300mV_Prof_Tree.root");
    TFile *inFile = new TFile("/hepstore/g2share/ProfTree/3ModuleData/Lab3TreeDump.root");
//
    TTree *strawTree = (TTree *)inFile->Get("professorTreeDumper/strawHits");
    TTree *scintTree = (TTree *)inFile->Get("professorTreeDumper/scintHits");
//
    StrawData strawHit;
    ScintData scintHit;
    strawTree->SetBranchAddress("strawHits",&strawHit.run);
    scintTree->SetBranchAddress("scintHits",&scintHit.run);
// -------------------------------------------------------------------
// Book all of our histograms:
// -------------------------------------------------------------------
    DefineColours();
    BookHistograms();
// -------------------------------------------------------------------
// set up Straw Geometry:
// -------------------------------------------------------------------
    int i1 = 1;
    Bkzero(i1);
    Strini();
// -------------------------------------------------------------------
    int nTotal = strawTree->GetEntries();
    cout<<" Ntotal: "<<nTotal<<endl;
// -------------------------------------------------------------------
// Loop round 10 times filling 100 hits into arrays.
// -------------------------------------------------------------------
    for(int iloop=1; iloop<=3000; iloop++) {
       int istart = 100*(iloop-1) + 1;
       int iend   = istart + 99;
       for(int ia=1; ia<=100; ia++) {
          Iruns[ia-1] = 0;
          Ievts[ia-1] = 0;
          TrigTimes[ia-1] = 0.0;
          Layers[ia-1] = 0;
          Istraw[ia-1] = 0;
          Module[ia-1] = 0;
          RawTimes[ia-1] = 0.0;
          Iviews[ia-1] = 0;
       }
//
       int iev = 0;
       for (int jentry=(istart-1); jentry <= (iend-1); jentry++) {
           strawTree->GetEntry(jentry);
           scintTree->GetEntry(jentry);
           int jstraw = strawHit.wire;
           int iok = 1;
//         if(jstraw == 1) {
//                         iok = 0;
//         }
           if(iok == 1) {
             iev = iev + 1;
//
             Iruns[iev-1]      = strawHit.run;
             Ievts[iev-1]      = strawHit.event;
             Istraw[iev-1]     = strawHit.wire;
             Layers[iev-1]     = strawHit.layer;
             Module[iev-1]     = strawHit.module;
             TrigTimes[iev-1]  = strawHit.triggerTime;
             RawTimes[iev-1]   = strawHit.hitTime; 
             Iviews[iev-1]     = strawHit.iview;
//
             int jrunscint     = scintHit.run;
             int jevscint      = scintHit.event;
             double Trigscint  = scintHit.triggerTime;
             double Hitscint   = scintHit.hitTime;
             int Statscint     = scintHit.station;
             int Modscint      = scintHit.module;
             int jviewscint    = scintHit.iview;
             int Layscint      = scintHit.layer;
             int jscint        = scintHit.wire;
//     
             if(Layers[iev-1] > -1) {
               if(Iviews[iev-1] == 1) {
                        Layers[iev-1] = Layers[iev-1] + 2;
               }
//
               Istraw[iev-1]     = Istraw[iev-1] + 1;
               Layers[iev-1]     = Layers[iev-1] + 1;
               Module[iev-1]     = Module[iev-1] + 1;
               StrHist.fWire[0]->Fill(Istraw[iev-1]-0.5);
               StrHist.fModule[0]->Fill(Module[iev-1]-0.5);
               StrHist.fLayer[0]->Fill(Layers[iev-1]-0.5);
//
//             cout <<" Run "<<Iruns[iev-1] 
//                  <<" Straw "<<Istraw[iev-1]<<" Layer: "<<Layers[iev-1]<< " Module: "<< Module[iev-1]
//                  <<" Event: "<<Ievts[iev-1]<<" Trig: "<<TrigTimes[iev-1]
//                  <<" HitTime: "<<RawTimes[iev-1]<<endl;
//
//                  <<" Iview    "<<Iviews[iev-1]<<endl;
//
//             cout <<" Scint "<<jscint<<" Layer: "<<Layscint<< " Run: "<< jrunscint
//                  <<" Event: "<<jevscint<<" Trig: "<<Trigscint
//                  <<" HitTime: "<<Hitscint
//                  <<" Iview    "<<jviewscint<<endl;
//
//                int Iblob;
//                cin>>Iblob;
             }
           }
       }
// -------------------------------------------------------------------
// We now have our local arrays filled with details for 100 hits.
// Try to make sense of them.
// Split the 100 hits according to the event number:
// -------------------------------------------------------------------
       int Jcon1   = 1;
       int Jcon2   = 1;
       int JevNo   = 0;
       int JevNext = 0;
       int istep1  = 0;
       int istep2  = 0;
       NhitsInEvent= 0;
       int Ngood   = 0;
//
       do {
          NhitsInEvent = 0;
          Ngood        = 0;
          istep1 = istep1 + 1;
          Jcon1 = 1;
          JevNo = Ievts[istep1 - 1];
          NhitsInEvent = NhitsInEvent + 1;
          Ngood = Ngood + 1;
          ModEvent[Ngood - 1] = Module[istep1-1];
          LayEvent[Ngood - 1] = Layers[istep1-1];
          IstEvent[Ngood - 1] = Istraw[istep1-1];
          RawTevent[Ngood - 1]= RawTimes[istep1-1];
          istep2 = istep1;
//        cout<<" In outer loop: istep1 is: "<<istep1<<endl;
//        cout<<" JevNo is: "<<JevNo<<endl;
          do {
             Jcon2  = 0;
             istep2 = istep2 + 1;
//           cout<<" In inner loop: istep2 is: "<<istep2<<endl;
             JevNext = Ievts[istep2-1];
//           cout<<" JevNext is: "<<JevNext<<endl;
             if(JevNext == JevNo) {
                                Jcon2 = 1;
                                NhitsInEvent = NhitsInEvent + 1;
                                int GoodHit = 1; 
                                for(int itest=1; itest<=Ngood; itest++) {
                                   int modtst = ModEvent[itest-1];
                                   int laytst = LayEvent[itest-1];
                                   int isttst = IstEvent[itest-1];
                                   if(Module[istep2-1] == modtst) {
                                     if(Layers[istep2-1] == laytst) {
                                       if(Istraw[istep2-1] == isttst) { 
                                                                      GoodHit = 0;
                                       }
                                     }
                                   }
                                }
                                if(GoodHit == 1) {
                                  Ngood = Ngood + 1;
                                  ModEvent[Ngood - 1] = Module[istep2-1];
                                  LayEvent[Ngood - 1] = Layers[istep2-1];
                                  IstEvent[Ngood - 1] = Istraw[istep2-1];
                                  RawTevent[Ngood - 1]= RawTimes[istep2-1];
                                }
             } 
             if(JevNext != JevNo) {
                                  Jcon2 = 0;
             }
             if(istep2 > 98) {
                             Jcon2 = 0;
             }
          } while (Jcon2 == 1);
// -------------------------------------------------------------------
          istep1 = istep2 - 1;
//        cout<<" Back out. istep1 is now: "<<istep1<<endl;
//        cout<<" Event "<<JevNo<<" has "<<NhitsInEvent<<" hits. "<<endl;
//        cout<<" Event "<<JevNo<<" has "<<Ngood<<" stored hits. "<<endl;
// -------------------------------------------------------------------
// For tracking to proceed demand that there is a module with one hit
// in each layer that look a reasonable bet to be a track:
// -------------------------------------------------------------------
          int DoTracking = 0;
          int i2 = 2;
          Bkzero(i2);
          Nhmods[0] = 0;
          Nhmods[1] = 0;
          Nhmods[2] = 0;
          for(int igood=1; igood<=Ngood; igood++) {
             int modd =  ModEvent[igood - 1];
             Nhmods[modd-1] = Nhmods[modd-1] + 1;
             int layy =  LayEvent[igood - 1];
             int istt =  IstEvent[igood - 1]; 
             Rtimes[istt-1][layy-1][modd-1] = RawTevent[igood - 1];
             float rrrr =  RawTevent[igood - 1];
//           cout<<" Unpacking hit in "<<istt<<" "<<layy<<" "<<modd<<" "<<rrrr<<endl;
          }
//        cout<<" Numbers of hits in each module: "<<endl;
//        cout<<" Module 1: "<<Nhmods[0]<<endl;
//        cout<<" Module 2: "<<Nhmods[1]<<endl;
//        cout<<" Module 3: "<<Nhmods[2]<<endl;
// --------------------------------------
// If a module has 4 hits check they are
// in the different layers:
// --------------------------------------
          for(int imod=1; imod<=3; imod++) {
             if(Nhmods[imod-1] == 4) {
// Module has 4 hits:
               Nhlays[0] = 0;     
               Nhlays[1] = 0;
               Nhlays[2] = 0;
               Nhlays[3] = 0;
               int is1   = 0;
               int is2   = 0;
               int is3   = 0;
               int is4   = 0;
               for(int ist=1; ist<=32; ist++) {
                  if(Rtimes[ist-1][0][imod-1] > 0.0) {
                    Nhlays[0] = Nhlays[0] + 1;
                    is1 = ist;
                  }
                  if(Rtimes[ist-1][1][imod-1] > 0.0) {
                    Nhlays[1] = Nhlays[1] + 1;
                    is2 = ist;
                  }
                  if(Rtimes[ist-1][2][imod-1] > 0.0) {
                    Nhlays[2] = Nhlays[2] + 1;
                    is3 = ist;
                  }
                  if(Rtimes[ist-1][3][imod-1] > 0.0) {
                    Nhlays[3] = Nhlays[3] + 1;
                    is4 = ist;
                  }
               }
               if(Nhlays[0] == 1 && Nhlays[1] == 1) {
                 if(Nhlays[2] == 1 && Nhlays[3] == 1) {
// Module has 4 hits in separate layers:  
                   if(abs(is1 - is2) < 2) {
                     if(abs(is1 - is3) < 3) {
                       if(abs(is3 - is4) < 2) {
                                              DoTracking = 1;
                       }
                     }
                   }
                 }
               }
             }
          }
//        if(JevNo == 9386) {
//                          DoTracking = 1;
//        }
          if(DoTracking == 1) {
//          cout<<" About to track in event: "<<JevNo<<endl;
            for(int imod=1;imod<=3; imod++) {
               for(int ilay=1; ilay<=4; ilay++) {
                  for(int istr=1; istr<=32; istr++) {
                     if(Rtimes[istr-1][ilay-1][imod-1] > 0.0) {
                       float rrr = Rtimes[istr-1][ilay-1][imod-1];
//                     cout<<" Hit in "<<istr<<" "<<ilay<<" "<<imod<<" "<<rrr<<endl;
                     }
                  }
               }
            }
//          int iblob2;
//          cin>>iblob2;
// -------------------------------------------------------------------
// The Goodt0 flag indicates whether a reliable t0 has been obtained:
// -------------------------------------------------------------------
            int Goodt0 = 0;
            Strevt(Goodt0);
            if(Goodt0 == 1) {
                           Strana();
            }
          }
// -------------------------------------------------------------------
          if(istep1 > 98) {
                          Jcon1 = 0;
          }
       } while (Jcon1 == 1);
// -------------------------------------------------------------------
// Get the next 100 hits:
// -------------------------------------------------------------------
    }
// -------------------------------------------------------------------
// And finally plot the histograms:
// -------------------------------------------------------------------
//  WriteStraws();
    PlotHistograms();
    WriteHistograms();
}
// ----------------------------------------------------------------------
void BookHistograms() {
// ----------------------------------------------------------------------
// Book all histograms for the analysis:
// ----------------------------------------------------------------------
   BookStrHistograms();
}
// --------------------------------------------------------------------
void DefineColours() {
// --------------------------------------------------------------------
       TColor *col229 = new TColor(229,0.01,0.11,0.18,"col229");
       TColor *col230 = new TColor(230,0.31,0.51,0.11,"col230");
       TColor *col231 = new TColor(231,0.11,0.51,0.71,"col231");
       TColor *col232 = new TColor(232,0.21,0.51,0.14,"col232");
       TColor *col233 = new TColor(233,0.43,0.21,0.45,"col233");
       TColor *col234 = new TColor(234,0.15,0.75,0.31,"col234");
       TColor *col235 = new TColor(235,0.22,0.14,0.16,"col235");
       TColor *col236 = new TColor(236,0.71,0.15,0.11,"col236");
       TColor *col237 = new TColor(237,0.61,0.25,0.21,"col237");
       TColor *col238 = new TColor(238,0.51,0.35,0.06,"col238");
       TColor *col239 = new TColor(239,0.41,0.45,0.16,"col239");
       TColor *col240 = new TColor(240,0.31,0.25,0.06,"col240");
       TColor *col241 = new TColor(241,0.99,0.61,0.61,"col241");
       TColor *col242 = new TColor(242,0.31,0.61,0.61,"col242");
       TColor *col243 = new TColor(243,0.81,0.61,0.41,"col243");
       TColor *col244 = new TColor(244,0.41,0.61,0.81,"col244");
       TColor *col245 = new TColor(245,0.21,0.91,0.41,"col245");
       TColor *col246 = new TColor(246,0.71,0.61,0.51,"col246");
       TColor *col247 = new TColor(247,0.11,0.41,0.81,"col247");
       TColor *col248 = new TColor(248,0.81,0.11,0.95,"col248");
       TColor *col249 = new TColor(249,0.61,0.11,0.41,"col249");
       TColor *col250 = new TColor(250,0.81,0.21,0.53,"col250");
       TColor *col251 = new TColor(251,0.31,0.41,0.61,"col251");
//
       TColor *sjm1   = new TColor(10001, 0.323, 0.189, 0.022, "mysludg1");
       TColor *sjm2   = new TColor(10002, 0.523, 0.189, 0.022, "mysludg2");
       TColor *sjm3   = new TColor(10003, 0.723, 0.189, 0.022, "mysludg3");
       TColor *sjm4   = new TColor(10004, 0.223, 0.389, 0.022, "mysludg4");
       TColor *sjm5   = new TColor(10005, 0.223, 0.589, 0.022, "mysludg5");
       TColor *sjm6   = new TColor(10006, 0.223, 0.789, 0.122, "mysludg6");
       TColor *sjm7   = new TColor(10007, 0.223, 0.089, 0.322, "mysludg7");
       TColor *sjm8   = new TColor(10008, 0.223, 0.089, 0.522, "mysludg8");
       TColor *sjm9   = new TColor(10009, 0.223, 0.089, 0.722, "mysludg9");
       TColor *sjm10  = new TColor(10010, 0.189, 0.589, 0.589, "mysludg10");
       TColor *sjm11  = new TColor(10011, 0.389, 0.589, 0.589, "mysludg11");
       TColor *sjm12  = new TColor(10012, 0.589, 0.589, 0.589, "mysludg12");
       TColor *sjm13  = new TColor(10013, 0.600, 0.600, 0.100, "mysludg13");
       TColor *sjm14  = new TColor(10014, 0.600, 0.600, 0.500, "mysludg14");
       TColor *sjm15  = new TColor(10015, 0.600, 0.600, 0.700, "mysludg15");
//
       TColor *col500 = new TColor(500,0.81,0.06,0.44,"col500");
       TColor *col499 = new TColor(499,0.97,1.00,0.16,"col499");
       TColor *col498 = new TColor(498,0.44,0.87,0.89,"col498");
       TColor *col497 = new TColor(497,0.89,0.47,0.10,"col497");
       TColor *col496 = new TColor(496,0.70,0.38,0.84,"col496");
       TColor *col495 = new TColor(495,0.27,0.86,0.66,"col495");
       TColor *col494 = new TColor(494,0.02,0.00,0.86,"col494");
       TColor *col493 = new TColor(493,0.00,0.86,0.02,"col493");
       TColor *col492 = new TColor(492,0.85,0.35,0.01,"col492");
       TColor *col491 = new TColor(491,1.00,0.53,0.57,"col491");
       TColor *col490 = new TColor(490,0.14,0.82,0.68,"col490");
       TColor *col489 = new TColor(489,0.45,0.11,0.29,"col489");
       TColor *col488 = new TColor(488,0.22,0.40,0.45,"col488");
       TColor *col487 = new TColor(487,0.22,0.50,0.16,"col487");
       TColor *col486 = new TColor(486,0.50,0.12,0.03,"col486");
       TColor *col485 = new TColor(485,0.98,0.50,0.14,"col485");
       TColor *col484 = new TColor(484,0.79,0.51,0.96,"col484");
       TColor *col483 = new TColor(483,0.96,0.75,0.63,"col483");
       TColor *col482 = new TColor(482,0.13,1.00,0.00,"col482");
       TColor *col481 = new TColor(481,0.98,1.00,0.11,"col481");
       TColor *col480 = new TColor(480,0.42,1.00,0.79,"col480");
       TColor *col479 = new TColor(479,0.57,0.31,1.00,"col479");
       TColor *col478 = new TColor(478,1.00,0.27,0.43,"col478");
       TColor *col477 = new TColor(477,0.11,0.86,0.05,"col477");
       TColor *col476 = new TColor(476,0.86,0.50,0.04,"col476");
       TColor *col475 = new TColor(475,0.09,0.75,0.05,"col475");
//
       TColor *col301 = new TColor(301,1.00,0.00,0.00,"col301");
       TColor *col302 = new TColor(302,0.00,1.00,0.00,"col302");
       TColor *col303 = new TColor(303,0.00,0.00,1.00,"col303");
       TColor *col304 = new TColor(304,0.00,1.00,1.00,"col304");
       TColor *col305 = new TColor(305,1.00,1.00,0.00,"col305");
       TColor *col306 = new TColor(306,1.00,0.00,1.00,"col306");
       TColor *col307 = new TColor(307,0.157,0.325,1.00,"col307");
       TColor *col308 = new TColor(308,1.00,0.50,0.00,"col308");
       TColor *col309 = new TColor(309,0.714,0.427,1.00,"col309");
       TColor *col310 = new TColor(310,0.50,0.40,1.00,"col310");
       TColor *col311 = new TColor(311,1.00,0.341,0.00,"col311");
       TColor *col312 = new TColor(312,0.255,0.827,1.00,"col312");
       TColor *col313 = new TColor(313,0.624,1.00,0.00,"col313");
       TColor *col314 = new TColor(314,0.620,1.00,0.827,"col314");
       TColor *col315 = new TColor(315,0.565,0.00,1.00,"col315");
       TColor *col316 = new TColor(316,0.5,0.627,1.00,"col316");
       TColor *col317 = new TColor(317,0.325,1.00,0.361,"col317");
       TColor *col318 = new TColor(318,0.00,0.00,0.00,"col318");
       TColor *col319 = new TColor(319,0.00,0.827,0.50,"col319");
       TColor *col320 = new TColor(320,0.50,0.5,0.50,"col320");
       TColor *col321 = new TColor(321,0.235,0.141,0.463,"col321");
       TColor *col322 = new TColor(322,0.180,0.055,0.482,"col322");
       TColor *col323 = new TColor(323,0.639,0.00,0.835,"col323");
       TColor *col324 = new TColor(324,0.216,0.078,0.289,"col324");
       TColor *col325 = new TColor(325,0.808,0.00,0.290,"col325");
       TColor *col326 = new TColor(326,0.910,0.718,0.00,"col326");
       TColor *col327 = new TColor(327,0.690,1.00,0.192,"col327");
       TColor *col328 = new TColor(328,0.051,0.922,0.847,"col328");
       TColor *col329 = new TColor(329,0.380,0.855,1.00,"col329");
//
       TColor *col526 = new TColor(526,1.0,0.30,0.30,"col526");
       TColor *col501 = new TColor(501,1.0,0.56,0.23,"col501");
       TColor *col502 = new TColor(502,1.0,0.90,0.25,"col502");
       TColor *col503 = new TColor(503,0.56,1.0,0.41,"col503");
       TColor *col504 = new TColor(504,0.55,0.89,0.91,"col504");
       TColor *col505 = new TColor(505,0.38,0.96,0.91,"col505");
       TColor *col506 = new TColor(506,0.02,0.76,0.91,"col506");
       TColor *col507 = new TColor(507,0.44,0.67,0.91,"col507");
       TColor *col508 = new TColor(508,0.51,0.29,0.91,"col508");
       TColor *col509 = new TColor(509,0.64,0.23,0.91,"col508");
       TColor *col510 = new TColor(510,0.85,0.48,0.91,"col510");
       TColor *col511 = new TColor(511,0.89,0.91,0.46,"col511");
       TColor *col512 = new TColor(512,0.55,0.91,0.53,"col512");
       TColor *col513 = new TColor(513,0.23,0.89,0.91,"col513");
       TColor *col514 = new TColor(514,0.42,0.53,0.91,"col514");
       TColor *col515 = new TColor(515,0.55,0.67,0.91,"col515");
       TColor *col516 = new TColor(516,0.72,0.59,0.91,"col516");
       TColor *col517 = new TColor(517,0.91,0.58,0.59,"col517");
       TColor *col518 = new TColor(518,0.91,0.72,0.89,"col518");
       TColor *col519 = new TColor(519,0.95,0.78,0.72,"col519");
       TColor *col520 = new TColor(520,0.95,0.95,0.70,"col520");
       TColor *col521 = new TColor(521,0.76,0.95,0.73,"col521");
       TColor *col522 = new TColor(522,0.69,0.95,0.92,"col522");
       TColor *col523 = new TColor(523,0.76,0.73,0.95,"col523");
       TColor *col524 = new TColor(524,0.78,1.00,0.36,"col524");
       TColor *col525 = new TColor(525,1.00,0.78,0.27,"col525");
}
// --------------------------------------------------------------------
void BookStrHistograms() {
// --------------------------------------------------------------------
// Book all the Histograms:
// --------------------------------------------------------------------
       for(int i=1; i<=Nmax; i++) {
          StrHist.fWire[i-1] = new TH1F(" Wire "," Wire ",33,0.0,33.0);
          StrHist.fWire[i-1]->SetFillColor(520-i);
          StrHist.fLayer[i-1] = new TH1F(" Layer "," Layer ",5,0.0,5.0);
          StrHist.fLayer[i-1]->SetFillColor(519-i);
          StrHist.fModule[i-1] = new TH1F(" Module "," Module ",4,0.0,4.0);
          StrHist.fModule[i-1]->SetFillColor(518-i);
          StrHist.fdiff12[i-1] = new TH1F(" diff12 "," diff12 ",120,-60.0,60.0);
          StrHist.fdiff12[i-1]->SetFillColor(517-i);
          StrHist.fdiff34[i-1] = new TH1F(" diff34 "," diff34 ",120,-60.0,60.0);
          StrHist.fdiff34[i-1]->SetFillColor(516-i);
          StrHist.fDtimes[i-1] = new TH1F(" Dtimes "," Dtimes ",100,0.0,100.0);
          StrHist.fDtimes[i-1]->SetFillColor(515-i);
          StrHist.fDdists[i-1] = new TH1F(" Ddists "," Ddists ",120,0.0,0.3);
          StrHist.fDdists[i-1]->SetFillColor(514-i);
          StrHist.fTrkMod[i-1] = new TH1F(" TrkMod "," TrkMod ",4,0.0,4.0);
          StrHist.fTrkMod[i-1]->SetFillColor(504-i);
          StrHist.fYbest[i-1] = new TH1F(" Ybest "," Ybest ",50,-5.0,5.0);
          StrHist.fYbest[i-1]->SetFillColor(513-i);
          StrHist.fGrfit[i-1] = new TH1F(" Grfit "," Grfit ",80,-2.0,2.0);
          StrHist.fGrfit[i-1]->SetFillColor(512-i);
          StrHist.fPhizxd[i-1] = new TH1F(" Phizxd "," Phizxd ",80,-40.0,40.0);
          StrHist.fPhizxd[i-1]->SetFillColor(511-i);
          StrHist.fZxte[i-1] = new TH2F(" Zxte "," Zxte ",100,0.0,50.0,90,0.0,30.0);
          StrHist.fRest1[i-1] = new TH1F(" Rest1 "," Rest1 ",100,-200.0,200.0);
          StrHist.fRest1[i-1]->SetFillColor(510-i);
          StrHist.fRest2[i-1] = new TH1F(" Rest2 "," Rest2 ",100,-200.0,200.0);
          StrHist.fRest2[i-1]->SetFillColor(509-i);
          StrHist.fRest3[i-1] = new TH1F(" Rest3 "," Rest3 ",80,-800.0,800.0);
          StrHist.fRest3[i-1]->SetFillColor(508-i);
          StrHist.fRest4[i-1] = new TH1F(" Rest4 "," Rest4 ",80,-800.0,800.0);
          StrHist.fRest4[i-1]->SetFillColor(507-i);
          StrHist.fRest[i-1] = new TH1F(" Rest "," Rest ",80,-800.0,800.0);
          StrHist.fRest[i-1]->SetFillColor(506-i);
       }
}
// --------------------------------------------------------------------
void PlotHistograms() {
// --------------------------------------------------------------------
     int Ipage = 0;
     PlotStrHistograms(Ipage);
}
// ----------------------------------------------------------------
void WriteHistograms() {
// ----------------------------------------------------------------
        int Ipage = 0;
        WriteStrHistograms(Ipage);
}
// ----------------------------------------------------------------
void WriteStrHistograms(int page) {
// ----------------------------------------------------------------
     if(page == 0) {
        TFile a("Straws.hist","recreate");
        for(int i=1; i<=Nmax; i++) {
          StrHist.fWire[i-1]->Write();
          StrHist.fLayer[i-1]->Write();
          StrHist.fModule[i-1]->Write();
          StrHist.fdiff12[i-1]->Write();
          StrHist.fdiff34[i-1]->Write();
          StrHist.fDtimes[i-1]->Write();
          StrHist.fDdists[i-1]->Write();
          StrHist.fTrkMod[i-1]->Write();
          StrHist.fYbest[i-1]->Write();
          StrHist.fGrfit[i-1]->Write();
          StrHist.fPhizxd[i-1]->Write();
          StrHist.fZxte[i-1]->Write();
          StrHist.fRest1[i-1]->Write();
          StrHist.fRest2[i-1]->Write();
          StrHist.fRest3[i-1]->Write();
          StrHist.fRest4[i-1]->Write();
          StrHist.fRest[i-1]->Write();
        }
     }
}
// ----------------------------------------------------------------
void PlotStrHistograms(int Ipage) {
// ----------------------------------------------------------------
// Plot out the Histograms:
// ----------------------------------------------------------------
   gStyle->SetMarkerStyle(8);
   gStyle->SetMarkerSize(0.5);
   gStyle->SetOptStat(111110);
   gStyle->SetStatH(0.25);
   gStyle->SetTitleH(0.085);
   gStyle->SetTitleW(0.6);
   gStyle->SetHistLineWidth(0.3);
   gStyle->SetPalette(1);
//
//   int Npar = 3;
//   double Binlow = -200.0;
//   double Binhi  =  200.0;
//   TF1 *FitFcn = new TF1("FitFcn",FitFunction,Binlow,Binhi,Npar);
//   FitFcn->SetNpx(500);
//   FitFcn->SetLineWidth(1.8);
//   FitFcn->SetLineColor(kSpring);
//   int Nparb = 3;
//   TF1 *Gau1Fcn = new TF1("Gau1Fcn",Gaussian1,Binlow,Binhi,Nparb);
//   Gau1Fcn->SetLineWidth(1.8);
//   Gau1Fcn->SetLineColor(kPink+10);
//   int Npars = 3;
//   TF1 *Gau2Fcn = new TF1("Gau2Fcn",Gaussian2,Binlow,Binhi,Npars);
//   Gau2Fcn->SetLineWidth(1.8);
//   Gau2Fcn->SetLineColor(kCyan+1);
//   Gau2Fcn->SetNpx(500);
// ----------------------------------------------------------------------------
   if(Ipage == 0) {
     TCanvas * c1 = new TCanvas("Straw Page 1","Straw Page 1",200,200,800,800);
     c1->Divide(4,2);
     for(int i=1; i<=8; i++) {
        c1->cd(i);
        if(i == 1) {
          StrHist.fWire[0]->Draw();
        }
        if(i == 2) {
          StrHist.fLayer[0]->Draw();
        }
        if(i == 3) {
          StrHist.fModule[0]->Draw();
        }
        if(i == 4) {
          StrHist.fdiff12[0]->Draw();
        }
        if(i == 5) {
          StrHist.fdiff34[0]->Draw();
        }
        if(i == 6) {
          StrHist.fDtimes[0]->Draw();
        }
        if(i == 7) {
          StrHist.fDdists[0]->Draw();
        }
        if(i == 8) {
          StrHist.fTrkMod[0]->Draw();
        }
     }
     c1->SaveAs("StrawsPage1.ps");
     c1->SaveAs("StrawsPage1.pdf");
     c1->SaveAs("StrawsPage1.gif");
   }
//
   if(Ipage == 0) {
     TCanvas * c1 = new TCanvas("Straw Page 2","Straw Page 2",200,200,800,800);
     c1->Divide(4,2);
     for(int i=1; i<=7; i++) {
        c1->cd(i);
        if(i == 1) {
          StrHist.fYbest[0]->Draw();
        }
        if(i == 2) {
          StrHist.fGrfit[0]->Draw();
        }
        if(i == 3) {
          StrHist.fPhizxd[0]->Draw();
        }
        if(i == 4) {
          StrHist.fZxte[0]->Draw();
        }
        if(i == 5) {
          StrHist.fRest1[0]->Draw();
        }
        if(i == 6) {
          StrHist.fRest2[0]->Draw();
        }
//      if(i == 7) {
//        StrHist.fRest3[0]->Draw();
//      }
        if(i == 7) {
          StrHist.fRest[0]->Draw();
        }
     }
     c1->SaveAs("StrawsPage2.ps");
     c1->SaveAs("StrawsPage2.pdf");
     c1->SaveAs("StrawsPage2.gif");
   }
//
}
// -----------------------------------------------------------------------------
// Track finding code:
// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void Vfill(float* Array, int Len , float Value) {
// -----------------------------------------------------------------------------
    for(int iloc=1; iloc<=Len; iloc++) {
                                       Array[iloc-1] = Value;
    }
}
// -----------------------------------------------------------------------------
void Vzero(float* Array, int Len) {
// -----------------------------------------------------------------------------
    for(int iloc=1; iloc<=Len; iloc++) {
                                       Array[iloc-1] = 0.0;
    }
}
// -----------------------------------------------------------------------------
void Strclp(float Gr,float c1,float Xpt,float Ypt,float &Xinter,float &Yinter) {
// -----------------------------------------------------------------------------
// Given a line specified by its gradient GR and intercept C1,the routine
// calculates the point (XINTER,YINTER) which is the point of closest
// approach of this line to a point (XPT,YPT)
// ----------------------------------------------------------------------
// The gradient of the perpendicular is -1./GR
// --------------------------------------------
// Get intercept of this perpendicular:
// intercept is the point on both lines:
// --------------------------------------------
      float Xxxbig = 9e9;
//
      Xinter = c1;
      Yinter = Ypt;
      if(Gr < Xxxbig) {
                      Xinter = ( (Ypt - c1) * Gr + Xpt)/(Gr*Gr + 1.0);
                      Yinter = Gr*Xinter + c1;
      }
}
// ------------------------------------------------------------------
void Sttime1(float Time,float &Dist) {
// ------------------------------------------------------------------
// Routine to use the time - distance relationship.
// ------------------------------------------------------------------
      float Dvel   = 4.7/1000.0;
      Dist = 0.0;
      if(Time > 0.0) {
                     Dist = Time*Dvel;
      }
      if(Dist > 0.3) {
                     Dist = 0.3;
      }
}
// ------------------------------------------------------------------
void Sttime2(float Time,float &Dist) {
// ------------------------------------------------------------------
// Routine to use the time - distance relationship.
// This version from Garfield: Ar/Eth 50:50  No magnetic field.
// Voltage = 1700V.
// ------------------------------------------------------------------
   float TimeArray[13];
   float DistArray[13];
//
   TimeArray[0] =  3.0;
   TimeArray[1] =  6.0;
   TimeArray[2] = 10.0;
   TimeArray[3] = 14.0;
   TimeArray[4] = 18.0;
   TimeArray[5] = 22.0;
   TimeArray[6] = 26.0;
   TimeArray[7] = 30.0;
   TimeArray[8] = 34.0;
   TimeArray[9] = 38.0;
   TimeArray[10]= 42.0;
   TimeArray[11]= 45.0;
   TimeArray[12]= 47.0;
//
   DistArray[0] = 0.01;
   DistArray[1] = 0.03;
   DistArray[2] = 0.05;
   DistArray[3] = 0.07;
   DistArray[4] = 0.09;
   DistArray[5] = 0.11;
   DistArray[6] = 0.13;
   DistArray[7] = 0.15;
   DistArray[8] = 0.17;
   DistArray[9] = 0.19;
   DistArray[10]= 0.21;
   DistArray[11]= 0.23;
   DistArray[12]= 0.25;
// ------------------------------------------------------------------
   if(Time < TimeArray[0]) {
                           Dist = 0.005;
                           return;
   }
   if(Time > TimeArray[12]) {
                           Dist = 0.25;
                           return;
   }
//
   int iloct = 0;
   for(int it=1; it<=13; it++) {
      if(Time >= TimeArray[it-1] && Time < TimeArray[it]) {
                                                       iloct = it; 
      }
   }
//
   float Widtht = TimeArray[iloct] - TimeArray[iloct-1];
   float Widthd = DistArray[iloct] - DistArray[iloct-1];
//
   float Extrat = Time - TimeArray[iloct-1];
   Dist = DistArray[iloct-1];
   Dist = Dist + (Widthd*Extrat/Widtht);
   if(Dist > 0.25) {
                   Dist = 0.25;
   }
}
// ------------------------------------------------------------------
void Sttime3(float Time,float &Dist) {
// ------------------------------------------------------------------
// Routine to use the time - distance relationship.
// This version from Garfield: Ar/Eth 50:50  No magnetic field.
// Voltage = 1700V.
// ------------------------------------------------------------------
   float TimeArray[13];
   float DistArray[13];
//
   float Fac = 55.0/47.0;
   TimeArray[0] =  3.0*Fac;
   TimeArray[1] =  6.0*Fac;
   TimeArray[2] = 10.0*Fac;
   TimeArray[3] = 14.0*Fac;
   TimeArray[4] = 18.0*Fac;
   TimeArray[5] = 22.0*Fac;
   TimeArray[6] = 26.0*Fac;
   TimeArray[7] = 30.0*Fac;
   TimeArray[8] = 34.0*Fac;
   TimeArray[9] = 38.0*Fac;
   TimeArray[10]= 42.0*Fac;
   TimeArray[11]= 45.0*Fac;
   TimeArray[12]= 47.0*Fac;
//
   DistArray[0] = 0.01;
   DistArray[1] = 0.03;
   DistArray[2] = 0.05;
   DistArray[3] = 0.07;
   DistArray[4] = 0.09;
   DistArray[5] = 0.11;
   DistArray[6] = 0.13;
   DistArray[7] = 0.15;
   DistArray[8] = 0.17;
   DistArray[9] = 0.19;
   DistArray[10]= 0.21;
   DistArray[11]= 0.23;
   DistArray[12]= 0.25;
// ------------------------------------------------------------------
   if(Time < TimeArray[0]) {
                           Dist = 0.005;
                           return;
   }
   if(Time > TimeArray[12]) {
                           Dist = 0.25;
                           return;
   }
//
   int iloct = 0;
   for(int it=1; it<=13; it++) {
      if(Time >= TimeArray[it-1] && Time < TimeArray[it]) {
                                                       iloct = it;
      }
   }
//
   float Widtht = TimeArray[iloct] - TimeArray[iloct-1];
   float Widthd = DistArray[iloct] - DistArray[iloct-1];
//
   float Extrat = Time - TimeArray[iloct-1];
   Dist = DistArray[iloct-1];
   Dist = Dist + (Widthd*Extrat/Widtht);
   if(Dist > 0.25) {
                   Dist = 0.25;
   }
}
// ------------------------------------------------------------------
void Strflr(int Nstraw, int Layer, int Nmod, float zpt, float xpt) {
// ------------------------------------------------------------------
// Fills array with real space points.
// Stores space points on track:
// ------------------------------------------------------------------
      Nposte = Nposte + 1;
      if(Nposte < 5) {
        Xypost[0][Nposte-1][Nmod-1] = zpt;
        Xypost[1][Nposte-1][Nmod-1] = xpt;
        Nspose[Nposte-1][Nmod-1]    = Nstraw;
        Nlpose[Nposte-1][Nmod-1]    = Layer;
//      cout<<" Strflr: "<<Nposte<<" "<<Nstraw<<" "<<Layer<<" "<<Nmod<<" "<<zpt<<" "<<xpt<<endl;
      }
}
// --------------------------------------------------------------------------------------------------------
void Strcht(float zz1,float xx1,float zz2,float xx2,int Nstraw,int Layer,int Nmod,float &Zpt, float &Xpt) {
// --------------------------------------------------------------------------------------------------------
// Compute the point(Zpt,Xpt) nearest to the tangent in cell Nstraw.
// The tangent is the line through ( zz1,xx1) , (zz2,xx2)
// ---------------------------------------------------------------------
      float Pibk  = 3.1415926535;
      float Twopi = 2.0*Pibk;
      float Piby2 = Pibk/2.0;
// ---------------------------------------------------------------------
//    if(Nstraw < 1 || Nstraw > 8) {
//        cout<<" In Strcht: "<<Nstraw<<endl;
//    }
      float Rr1 = sqrt(zz1*zz1 + xx1*xx1);
      float Rr2 = sqrt(zz2*zz2 + xx2*xx2);
      float  Za = zz1;
      float  Xa = xx1;
      float  Zb = zz2;
      float  Xb = xx2;
      if(Rr1 > Rr2) {
                    Za = zz2;
                    Xa = xx2;
                    Zb = zz1;
                    Xb = xx1;
      }
      int Iflag  = 0;
      float Thet = 0.0;
      Strcor(Za,Xa,Zb,Xb,Nstraw,Layer,Nmod,Thet,Iflag);
      float z3l = Hits[Nstraw-1][Layer-1][Nmod-1]*cos(Thet);
      float x3l = Hits[Nstraw-1][Layer-1][Nmod-1]*sin(Thet);
// ---------------------------------------------------------------
// Convert space point back to GM2 coordinates:
// ---------------------------------------------------------------
      float Sth = Sthpla;
      float Cth = Cthpla;
      if(Iflag == 1) {
                     float The = Thecol + (Piby2/4.0);
                     Sth = sin(The);
                     Cth = cos(The);
      }
      Zpt = z3l*Sth + x3l*Cth + Xzaty[0][Nstraw-1][Layer-1][Nmod-1];
      Xpt = x3l*Sth - z3l*Cth + Xzaty[1][Nstraw-1][Layer-1][Nmod-1];
}
// ------------------------------------------------------------------------------------------------------
void Strcor(float zz1,float xx1,float zz2,float xx2,int Nst,int Nlay,int Nmod,float &Angle, int &Iflag) {
// ------------------------------------------------------------------------------------------------------
// Routine to compute the angle of the track (zz1,xx1),(zz2,xx2) in the
// local frame based on the cell Nst,Nlay.
// --------------------------------------------------------------------
      float Pibk = 3.1415926535;
      float Piby2 = Pibk/2.0;
      float Pi3hlf = Piby2 * 3.0;
      float Twopi = 2.0*Pibk;
      float Dysmal = 1e-6;
//
      float Xlow = -10.0;
// ---------------------------------------------------------------------
// convert zz1,xx1,zz2,xx2 into local cell coordinates.
// ---------------------------------------------------------------------
      Iflag = 0;
      float Sth = Sthpla;
      float Cth = Cthpla;
//
      float z1l = (zz1 - Xzaty[0][Nst-1][Nlay-1][Nmod-1])*Sth - (xx1 - Xzaty[1][Nst-1][Nlay-1][Nmod-1])*Cth;
      float x1l = (xx1 - Xzaty[1][Nst-1][Nlay-1][Nmod-1])*Sth + (zz1 - Xzaty[0][Nst-1][Nlay-1][Nmod-1])*Cth;
      float Dz  = (zz2 - zz1)*Sth - (xx2 - xx1)*Cth;
      float Dx  = (xx2 - xx1)*Sth + (zz2 - zz1)*Cth;
// -----------------------------------------------------------------------------------------
      if(fabs(Dx) < Dysmal) {
            Iflag = 1;
            float The = Thecol + (Piby2/4.0);
            Sth = sin(The);
            Cth = cos(The);
            z1l = (zz1 - Xzaty[0][Nst-1][Nlay-1][Nmod-1])*Sth - (xx1 - Xzaty[1][Nst-1][Nlay-1][Nmod-1])*Cth;
            x1l = (xx1 - Xzaty[1][Nst-1][Nlay-1][Nmod-1])*Sth + (zz1 - Xzaty[0][Nst-1][Nlay-1][Nmod-1])*Cth;
            Dz  = (zz2 - zz1)*Sth - (xx2 - xx1)*Cth;
            Dx  = (xx2 - xx1)*Sth + (zz2 - zz1)*Cth;
      }
      float Alpha = atan2(Dx,Dz);
      Angle = Alpha + Pi3hlf;
      if(Alpha > atan2(-Xlow,(x1l-Xlow)*Dz/Dx - z1l)) {
                                                      Angle = Alpha + Piby2;
      }
      if(Angle < 0.0) {
                      Angle = Angle + Twopi;
      }
}
// --------------------------------------------------------------------
void Strelm(int Nmod,int Nst1,int Layer1,int Nst2,int Layer2,int Icand,
            float &zz1,float &xx1,float &zz2,float &xx2,int &Iflag) {
// --------------------------------------------------------------------
// Finds the common tangent to the drift circles in cells Nst1 and
// Nst2 , for the tangent requested by the value of ICAND.
// --------------------------------------------------------------------
      float Pibk = 3.1415926535;
      float Twopi = 2.0*Pibk;
      float Piby2 = Pibk/2.0;
// ----------------------------
// Pick up drift distances
// ----------------------------
      Iflag    = 0;
      float d1 = Hits[Nst1-1][Layer1-1][Nmod-1];
      float d2 = Hits[Nst2-1][Layer2-1][Nmod-1];
      float Xyphys1  = Size[Layer1-1];
      float Xyphys2  = Size[Layer2-1];
      if(d1 > Xyphys1) {
                       d1 = Xyphys1;
      }
      if(d2 > Xyphys2) {
                       d2 = Xyphys2;
      }
      float Zcell1 = Xzaty[0][Nst1-1][Layer1-1][Nmod-1];
      float Xcell1 = Xzaty[1][Nst1-1][Layer1-1][Nmod-1];
      float Zcell2 = Xzaty[0][Nst2-1][Layer2-1][Nmod-1];
      float Xcell2 = Xzaty[1][Nst2-1][Layer2-1][Nmod-1];
// ------------------------------------
// Distances between cell centres
// ------------------------------------
      float D = sqrt((Zcell1-Zcell2)*(Zcell1-Zcell2) + (Xcell1-Xcell2)*(Xcell1-Xcell2));
// --------------------------------------
// Angle of line joining cell centres
// --------------------------------------
      float Alpha = atan2((Xcell2-Xcell1) , (Zcell2-Zcell1));
//printf("hello %f\n", D);
      float Phi = 0.0;
      float Delta1 = 0.0;
      float Delta2 = 0.0;
      if(Icand == 3) {
// -----------------
// Crossed tangents
// -----------------
        if(fabs(d1+d2) > D) {
             Iflag = 1;
             return;
        }
        Phi    = acos( (d1+d2)/D);
        Delta1 = Alpha - Phi;
        Delta2 = Delta1 + Pibk;
      }
      if(Icand == 4) {
        if(fabs(d1+d2) > D) {
          Iflag = 1;
          return;
        }
        Phi    = acos((d1+d2)/D);
        Delta1 = Alpha + Phi;
        Delta2 = Delta1 + Pibk;
      }
      if(Icand == 2) {
// ------------------------
// 'outer edge' tangents
// ------------------------
        if(fabs(d2-d1) > D) {
           Iflag = 1;
           return;
        }
        Phi    = asin((d2-d1)/D);
        Delta1 = Alpha + Phi + Pibk/2.0;
        Delta2 = Delta1;
      }
      if(Icand == 1) {
        if(fabs(d2-d1) > D) {
          Iflag = 1;
          return;
        }
        Phi    = asin((d2-d1)/D);
        Delta1 = Alpha - Phi + 3.0*Pibk/2.0;
        Delta2 = Delta1;
      }
//
      zz1 = Zcell1 + d1*cos(Delta1);
      xx1 = Xcell1 + d1*sin(Delta1);
      zz2 = Zcell2 + d2*cos(Delta2);
      xx2 = Xcell2 + d2*sin(Delta2);
//printf("hello %f\n", Xzaty[1][Nst1-1][Layer1-1][Nmod-1]);
//printf("hello %f\n", zz1);
}
// --------------------------------------------------------------------
void Linfit(float* Xpts,float* Ypts,int Npts,float &Grad,float &Cint) {
// --------------------------------------------------------------------
// Simple straight line fitting routine.
// --------------------------------------------------------------------
      if(Npts < 3) {
                   Grad = 0.0;
                   Cint = 0.0;
                   return;
      }
//
      float Count = 0.0;
      float Sumx  = 0.0;
      float Sumy  = 0.0;
      float Sumxy = 0.0;
      float Sumxx = 0.0;
      float Sumyy = 0.0;
      for(int j=1; j<=Npts; j++) {
         if(Ypts[j-1] != 0.0) {
               Sumx  = Sumx + Xpts[j-1];
               Sumy  = Sumy + Ypts[j-1];
               Count = Count + 1.0;
         }
      }
//
      if(Count <= 1.0) {
                       Grad = 0.0;
                       Cint = 0.0;
                       return;
      }
      float Ymed = Sumy/Count;
      float Xmed = Sumx/Count;
//
      for(int j=1; j<=Npts; j++) {
         if(Ypts[j-1] != 0.0) {
            float Scartx = Xpts[j-1] - Xmed;
            float Scarty = Ypts[j-1] - Ymed;
            Sumxy  = Sumxy + Scartx*Scarty;
            Sumxx  = Sumxx + Scartx*Scartx;
            Sumyy  = Sumyy + Scarty*Scarty;
         }
      }
// -----------------------------------------------
// Fit Parameters:
// -----------------------------------------------
      if(Sumxx == 0.0) {
                       Grad = 0.0;
                       Cint = 0.0;
                       return;
      }
//
      float A = Sumxy/Sumxx;
      float B = Ymed - A*Xmed;
      float E = 0.0;
      if(Count >= 3.0) {
              E = (Sumyy - Sumxy*A)/(Count-2.0);
      }
      Grad = A;
      Cint = B;
}
// ----------------------------------------------------------
float Strphi()  {
// ----------------------------------------------------------
// Whilst trying to make a 3 hit TE in Module 1 get the
// angle of the nearest 4 hit TE.
// ----------------------------------------------------------
      float Pibk = 3.1415926535;
      float Twopi= 2.0*Pibk;
//
      float Angexp = 0.0;
      return Angexp;
}
// ------------------------------------------------------
float Acosbk(float Rawcos) {
// ------------------------------------------------------
      float Cosnew = Rawcos;
      if(Rawcos > 1.0) {
                       Cosnew = 0.99999;
      }
      if(Rawcos < -1.0) {
                        Cosnew = -0.99999;
      }
      float Result = acos(Cosnew);
      return Result;
}
// ------------------------------------------------------
float Bkatan(float py , float px) {
// ------------------------------------------------------
// Calculate phi between 0.0 and twopi.
// ------------------------------------------------------
      float tanphi;
      float phi;
      float phitemp;
      float pibk;
      pibk = 3.1415926535;
      phi = 0.0;
      if(px > 0.0 && py > 0.0) {
         tanphi = py / px ;
         phi    = atan(tanphi);
      }
      if(px < 0.0 && py > 0.0) {
         tanphi  = py / px ;
         phitemp = atan(tanphi);
         phi     = pibk + phitemp;
      }
      if(px < 0.0 && py < 0.0) {
         tanphi  = py / px ;
         phitemp = atan(tanphi);
         phi     = phitemp + pibk;
      }
      if(px > 0.0 && py < 0.0) {
         tanphi  = py / px ;
         phitemp = atan(tanphi);
         phi     = phitemp + (2.0*pibk);
      }
      if(px > 0.0 && py < 0.0) {
         tanphi  = py / px ;
         phitemp = atan(tanphi);
         phi     = phitemp + (2.0*pibk);
      }
      return phi;
}
// ------------------------------------------------------
int Bknint(float Var) {
// ------------------------------------------------------
// Return closest integer to real variable Var.
// ------------------------------------------------------
    int Intclo;
// Get integer below Var:
    int Ivar     = Var;
    float Vard   = Ivar;
    float Diff   = Var - Vard;
    if(Diff < 0.5) {
                   Intclo = Ivar;
    }
    if(Diff >= 0.5) {
                   Intclo = Ivar+1;
    }
    return Intclo;
}
// ---------------------------------------------------------------------
void Strcal() {
// ---------------------------------------------------------------------
// Routine to get the STRAWS calibration constants.
// ---------------------------------------------------------------------
//    Common /Strawg/Zpos,Xzstraw(6,32,4), Diam, Ylength,
//   +               Size(4),Xzstrau(6,32,4), Xzstrav(6,32,4)
//    Common /Strawp/Iexist(32,4),Pedest(32,4),T0data,Tofl,Yped
//
      float Value = 0.0;
      int Len = 32*4;
      for(int imod=1; imod<=3; imod++) {
         for(int ia=1; ia<=32; ia++) {
            for(int ib=1; ib<=4; ib++) {
               Pedest[ia-1][ib-1][imod-1] = Value;
            }
         }
      }
      T0data = 0.0;
      Yped   = 0.33;
//    Tofl   = Zpos / 30.0;
      Tofl   = 0.0;
}
// --------------------------------------------------------------------
void Strana() {
// --------------------------------------------------------------------
// Having obtained the drift times and hits now find the track.
// --------------------------------------------------------------------
// Start by assuming the wire positions for the middle of the cells:
// This gets changed in any case in Strft4 etc.
// ------------------------------------------------------------------
      for(int imod=1; imod<=3; imod++) {
         for(int ilay=1; ilay<=4; ilay++) {
            for(int istr=1; istr<=32; istr++) {
               Xzaty[0][istr-1][ilay-1][imod-1] = 0.0;
               Xzaty[1][istr-1][ilay-1][imod-1] = 0.0;
            }
         }
      }
      for(int imod=1; imod<=3; imod++) {
         for(int ilay=1; ilay<=4; ilay++) {
            for(int istr=1; istr<=32; istr++) {
// z coordinate:
               Xzaty[0][istr-1][ilay-1][imod-1] = Xzstraw[2][istr-1][ilay-1][imod-1];
// x coordinate half way up the straw:
               Xzaty[1][istr-1][ilay-1][imod-1] = Xzstraw[0][istr-1][ilay-1][imod-1] 
                                                + Xzstraw[3][istr-1][ilay-1][imod-1];
            }
         }
      }
// ------------------------------------------------------------------
// Now find 4-hit track elements in each of the Straw Modules:
// ------------------------------------------------------------------
      for(int imod=1; imod<=3; imod++) {
         Strft4(imod);
// ------------------------------------------------------------------
// Now look for 3-hit track elements.
// ------------------------------------------------------------------
//       Strft3(imod);
      }
}
// ------------------------------------------------------------------
void Strft3(int Nmod) {
// ------------------------------------------------------------------
// Finds a 3-hit TE in a Straw Module.
// It makes calls to Strte2 and Strte3 to deal with 3 hit cases.
// ------------------------------------------------------------------
// -----------------------------------------
// Reset the number of TE points:
// -----------------------------------------
      Nposte = 0;
//
      int Ist1  = 0;
      int Ist2  = 0;
      int Ist3  = 0;
      int Ist4  = 0;
      int Nhits = 0;
      for(int istr=1; istr<=32; istr++) {
         if(Hits[istr-1][0][Nmod-1] > 0.0) {
                                   Ist1  = istr;
                                   Nhits = Nhits + 1;
         }
         if(Hits[istr-1][1][Nmod-1] > 0.0) {
                                   Ist2  = istr;
                                   Nhits = Nhits + 1;
         }
         if(Hits[istr-1][2][Nmod-1] > 0.0) {
                                   Ist3  = istr;
                                   Nhits = Nhits + 1;
         }
         if(Hits[istr-1][3][Nmod-1] > 0.0) {
                                   Ist4  = istr;
                                   Nhits = Nhits + 1;
         }
      }
      cout<<" Strft3: "<<Nmod<<" "<<Ist1<<" "<<Ist2<<" "<<Ist3<<" "<<Ist4<<endl;
// ------------------------------------------------------------------
// See if we can do anything with 3 hit cases.
// ------------------------------------------------------------------
      if(Nhits != 3) {
                     return;
      }
      if(Nhits == 3) {
        if(Ist1 > 0 && Ist2 > 0) {
          if(Ist3 == 0 || Ist4 == 0) {
            Strte2(Nmod,Ist1,Ist2,Ist3,Ist4);
            return;
          }
        }
//
        if(Ist3 > 0 && Ist4 > 0) {
          if(Ist1 == 0 || Ist2 == 0) {
            Strte3(Nmod,Ist1,Ist2,Ist3,Ist4);
            return;
          }
        }
      }
}
// -------------------------------------------------------------------
void Bkxy2u(float xx , float yy , float &xu , float &yu) {
// -------------------------------------------------------------------
// Translate a point in the xy coordinate system into the u frame.
// -------------------------------------------------------------------
     float Pibk = 3.1415926535;
     float Piby2 = Pibk/2.0;
// --------------------------------------------------------------------
// Convert from global to local coordinates:
// --------------------------------------------------------------------
// Xwire Ywire Theta are the positions of the local frames origin and 
// orientation as expressed in the global frame.
// --------------------------------------------------------------------
     float Xwire  = 0.0;
     float Ywire  = 0.0;
     float Thetad = 82.5;
     float Theta  = Thetad * Pibk/180.0;
     float Cth    = cos(Theta-Piby2);
     float Sth    = sin(Theta-Piby2);
// --------------------------------------------------------------
// Translate into the local frame - in this case the u frame.
// --------------------------------------------------------------
     xu = ( xx - Xwire )*Cth + ( yy - Ywire )*Sth; 
     yu = ( yy - Ywire )*Cth - ( xx - Xwire )*Sth;
}
// ---------------------------------------------------------------
void Bkxy2v(float xx , float yy , float &xv , float &yv) {
// ---------------------------------------------------------------
// Translate a point in the xy coordinate system into the v frame.
// ---------------------------------------------------------------
     float Pibk  = 3.1415926535;
     float Piby2 = Pibk/2.0;
// --------------------------------------------------------------------
// Convert from global to local coordinates:
// --------------------------------------------------------------------
// Xwire Ywire Theta are the positions of the local frames origin and 
// orientation as expressed in the global frame.
// --------------------------------------------------------------------
     float Xwire  = 0.0;
     float Ywire  = 0.0;
     float Thetad = 97.5;
     float Theta  = Thetad * Pibk/180.0;
     float Cth    = cos(Theta-Piby2);
     float Sth    = sin(Theta-Piby2);
// --------------------------------------------------------------
// Translate into the local frame - in this case the v frame.
// --------------------------------------------------------------
     xv = ( xx - Xwire )*Cth + ( yy - Ywire )*Sth; 
     yv = ( yy - Ywire )*Cth - ( xx - Xwire )*Sth;
}
// --------------------------------------------------------------------
void Bku2xy(float xu , float yu , float &xx , float &yy) {
// --------------------------------------------------------------------
// Translate a point in the u coordinate system into the xy frame.
// --------------------------------------------------------------------
     float Pibk  = 3.1415926535;
     float Piby2 = Pibk/2.0;
// --------------------------------------------------------------------
// Xwire Ywire Theta are the positions of the local frames origin and 
// orientation as expressed in the global frame.
// --------------------------------------------------------------------
     float Xwire  = 0.0;
     float Ywire  = 0.0;
     float Thetad = 82.5;
     float Theta  = Thetad * Pibk/180.0;
     float Cth    = cos(Theta-Piby2);
     float Sth    = sin(Theta-Piby2);
// --------------------------------------------------------------
// Now translate (Xlocal,Ylocal) into the global frame.
// --------------------------------------------------------------
     xx = xu*Cth - yu*Sth + Xwire;
     yy = yu*Cth + xu*Sth + Ywire;
}
// -------------------------------------------------------------------
void Bkv2xy(float xv , float yv , float &xx , float &yy) {
// --------------------------------------------------------------------
// Translate a point in the u coordinate system into the xy frame.
// --------------------------------------------------------------------
     float Pibk = 3.1415926535;
     float Piby2 = Pibk/2.0;
// --------------------------------------------------------------------
// Xwire Ywire Theta are the positions of the local frames origin and 
// orientation as expressed in the global frame.
// --------------------------------------------------------------------
     float Xwire  = 0.0;
     float Ywire  = 0.0;
     float Thetad = 97.5;
     float Theta  = Thetad * Pibk/180.0;
     float Cth    = cos(Theta-Piby2);
     float Sth    = sin(Theta-Piby2);
// --------------------------------------------------------------
// Now translate (Xlocal,Ylocal) into the global frame.
// --------------------------------------------------------------
     xx = xv*Cth - yv*Sth + Xwire;
     yy = yv*Cth + xv*Sth + Ywire;
}
// ----------------------------------------------------------------
float BkRndm() {
// ----------------------------------------------------------------
      IcallRndm = IcallRndm + 1;
      double x1d = 0.0;
      x1d = Rnd->Rndm();
      float x1 = x1d;
      return x1;
}
// ----------------------------------------------------------------
double FitBack(double* x, double* par) {
// ----------------------------------------------------------------
//  double Ans = par[0] + par[1]*exp(-x[0]/par[2]);
    double Ans = par[0] + par[1]*x[0] + par[2]*x[0]*x[0];
    return Ans;
}
// ----------------------------------------------------------------
double Lorentzian(double* x, double* par) {
// ----------------------------------------------------------------
    double pibk = 3.1415926535;
    double Part1 = 0.5*par[0]*par[1]/pibk;
    double Part2 = (x[0]-par[2])*(x[0]-par[2]) + 0.25*par[1]*par[1];
    if(Part2 < 1.0e-10) {
      Part2 = 1.0e-10;
    }
    double Ans   = Part1/Part2;
    return Ans;
}
// ----------------------------------------------------------------
double Gaussian1(double* x, double* par) {
// ----------------------------------------------------------------
    double pibk = 3.1415926535;
    double Twopi = 2.0*pibk;
    double Ans1 = 0.0;
    double Blob = 0.0;
    double Ans2 = 0.0;
    double Ans  = 0.0;
    if(par[2] != 0.0) {
      Ans1 = 1.0/(sqrt(Twopi*par[2]*par[2]));
      Blob   = x[0] - par[1];
      Ans2   = exp(-(Blob*Blob)/(2.0*par[2]*par[2]));
      Ans    = par[0]*Ans1*Ans2;
    }
//  if( x[0] < 100.0 || x[0] > 140.0) {
//                                    Ans = 0.0;
//  }
    return Ans;
}
// ----------------------------------------------------------------
double Gaussian2(double* x, double* par) {
// ----------------------------------------------------------------
    double pibk = 3.1415926535;
    double Twopi = 2.0*pibk;
    double Ans1 = 0.0;
    double Blob = 0.0;
    double Ans2 = 0.0;
    double Ans  = 0.0;
    if(par[2] != 0.0) {
      Ans1 = 1.0/(sqrt(Twopi*par[2]*par[2]));   
      Blob   = x[0] - par[1];
      Ans2   = exp(-(Blob*Blob)/(2.0*par[2]*par[2]));
      Ans    = par[0]*Ans1*Ans2;
    }
//  if( x[0] < 100.0 || x[0] > 140.0) {   
//                                    Ans = 0.0;
//  }
    return Ans;
}

// --------------------------------------------------------------
double FitFunction(double* x, double* par) {
// --------------------------------------------------------------
//    double Back = FitBack(x,par);
      double Gau1 = Gaussian1(x,par);
//    double Gau2 = Gaussian2(x,&par[3]);
//    double Ans  = Gau1 + Gau2;
      double Ans  = Gau1;
      return Ans;
}
// --------------------------------------------------------------
double FitFunction2(double* x, double* par) {
// --------------------------------------------------------------
//    double Back = FitBack(x,par);
      double Gau1 = Gaussian1(x,par);
      double Gau2 = Gaussian2(x,&par[3]);   
      double Ans  = Gau1 + Gau2;
      return Ans;
}
// --------------------------------------------------------------------
void Bkzero(int Imode) {
// --------------------------------------------------------------------
// Zero counters at run start.
// --------------------------------------------------------------------
     if(Imode == 1) {
                    Zpos = 0.0;
                    for(int ia=1; ia<=6; ia++) {
                       for(int ib=1; ib<=32; ib++) {
                          for(int ic=1; ic<=4; ic++) {                             
                              for(int id=1; id<=3; id++) {
                                 Xzstraw[ia-1][ib-1][ic-1][id-1] = 0.0;
                                 Xzstrau[ia-1][ib-1][ic-1][id-1] = 0.0;
                                 Xzstrav[ia-1][ib-1][ic-1][id-1] = 0.0;
                              }
                          }
                       }
                    }
                    Diam    = 0.0;
                    Ylength = 0.0;
                    for(int ia=1; ia<=32; ia++) {
                       for(int ib=1; ib<=4; ib++) {
                          for(int ic=1; ic<=3; ic++) {
                             Iexist[ia-1][ib-1][ic-1] = 0;
                             Pedest[ia-1][ib-1][ic-1] = 0.0;
                          }                      
                       }
                    }
                    T0data  = 0.0;
                    Tofl    = 0.0;
                    Yped    = 0.0;
                    for(int ia=1; ia<=32; ia++) {
                       for(int ib=1; ib<=4; ib++) {
                          for(int ic=1; ic<=3; ic++) {
                             Tof[ia-1][ib-1][ic-1] = 0.0;
                          }
                       }
                    }
     }
// --------------------------------------------------------------------
// Zero counters for new event:
// --------------------------------------------------------------------
     for(int ia=1; ia<=32; ia++) {
         for(int ib=1; ib<=4; ib++) {
            for(int ic=1; ic<=3; ic++) {
               Rtimes[ia-1][ib-1][ic-1] = 0.0;
               Dtimes[ia-1][ib-1][ic-1] = 0.0;
               Hits[ia-1][ib-1][ic-1]   = 0.0;
               Yhits[ia-1][ib-1][ic-1]  = 0.0;
               Mask[ia-1][ib-1][ic-1]   = 0;
               Zxhit[0][ia-1][ib-1][ic-1] = 0.0;
               Zxhit[1][ia-1][ib-1][ic-1] = 0.0;
            }
         }
      }
//
      Nposte = 0;
      for(int ic=1; ic<=4; ic++) {
         for(int imod=1; imod<=3; imod++) {
            Xypost[0][ic-1][imod-1] = 0.0;
            Xypost[1][ic-1][imod-1] = 0.0;
            Nspose[ic-1][imod-1]    = 0;
            Nlpose[ic-1][imod-1]    = 0;
         }
      }
      for(int imod=1; imod<=3; imod++) {
         for(int id=1; id<=60; id++) {
            Trkinf[id-1][imod-1] = 0.0;
            Trkin2[id-1][imod-1] = 0.0;
            Trkin3[id-1][imod-1] = 0.0;
         }
      } 
}
// --------------------------------------------------------------------
void Strini() {
// --------------------------------------------------------------------
     float Pibk   = 3.1415926535;
     float Twopi  = 2.0*Pibk;
     float Piby2  = Pibk/2.0;
//
     Thecol = Piby2/2.0;
     Sthpla = sin(Thecol);
     Cthpla = cos(Thecol);
//
     Diam     = 0.5;
     Ylength  = 10.0;
     float Thet = 7.5*Pibk/180.0;
     float Tanthe   = tan(Thet);   
     float costhe   = cos(Thet);
     float Wdist    = 0.6/costhe;
     float Whalf    = Wdist/2.0;
     for(int Ilay=1; Ilay<=4; Ilay++) {
        Size[Ilay-1] = Diam / 2.0;
     }
// -------------------------------------------------------------------
// Now set the Z coordinate of each individual layer...
// goes from Z1 to Z4 for the 4 layers in the module.
// -------------------------------------------------------------------
     float Zpos  = 10.0;
     float Zoff1 = 13.0;
     float Zoff2 = 26.0;
     float Z1  = Zpos + 3.931;
     float Z2  = Zpos + 4.408;
     float Z3  = Zpos + 6.091;
     float Z4  = Zpos + 6.568;
     float Z5  = Zpos + Zoff1 + 3.931;
     float Z6  = Zpos + Zoff1 + 4.408;
     float Z7  = Zpos + Zoff1 + 6.091;
     float Z8  = Zpos + Zoff1 + 6.568;
     float Z9  = Zpos + Zoff2 + 3.931;
     float Z10 = Zpos + Zoff2 + 4.408;
     float Z11 = Zpos + Zoff2 + 6.091;
     float Z12 = Zpos + Zoff2 + 6.568;
// --------------------------------------------------------------------
// Allow for the module to be offset in x:
// --------------------------------------------------------------------
     float Xoff1 =  0.0;
// ----------------------------------------------------------------------
     for(int ia=1; ia<=6; ia++) {
        for(int ib=1; ib<=32; ib++) {
           for(int ic=1; ic<=4; ic++) {
              for(int im=1; im<=3; im++) {
                 Xzstraw[ia-1][ib-1][ic-1][im-1] = 0.0;
              }
           }
        }
     }  
     float Ybot   = -4.55;
     float Ytop   =  4.55;
// ----------------------------------------------------------------------
// First Module:
// ----------------------------------------------------------------------
     for(int istraw=1;istraw<=32; istraw++) {
        Xzstraw[2][istraw-1][0][0] = Z1;
        Xzstraw[5][istraw-1][0][0] = Z1;
        Xzstraw[2][istraw-1][1][0] = Z2;
        Xzstraw[5][istraw-1][1][0] = Z2;
        Xzstraw[2][istraw-1][2][0] = Z3;
        Xzstraw[5][istraw-1][2][0] = Z3;
        Xzstraw[2][istraw-1][3][0] = Z4;
        Xzstraw[5][istraw-1][3][0] = Z4;
//
        Xzstraw[2][istraw-1][0][1] = Z5;
        Xzstraw[5][istraw-1][0][1] = Z5;
        Xzstraw[2][istraw-1][1][1] = Z6;
        Xzstraw[5][istraw-1][1][1] = Z6;
        Xzstraw[2][istraw-1][2][1] = Z7;
        Xzstraw[5][istraw-1][2][1] = Z7;
        Xzstraw[2][istraw-1][3][1] = Z8;
        Xzstraw[5][istraw-1][3][1] = Z8;
//
        Xzstraw[2][istraw-1][0][2] = Z9;
        Xzstraw[5][istraw-1][0][2] = Z9;
        Xzstraw[2][istraw-1][1][2] = Z10;
        Xzstraw[5][istraw-1][1][2] = Z10;
        Xzstraw[2][istraw-1][2][2] = Z11;
        Xzstraw[5][istraw-1][2][2] = Z11;
        Xzstraw[2][istraw-1][3][2] = Z12;
        Xzstraw[5][istraw-1][3][2] = Z12;
// Layer 1 - Bottom.
        Xzstraw[0][istraw-1][0][0] = Xoff1 + 4.7984 + float(istraw-1)*Wdist;  
        Xzstraw[1][istraw-1][0][0] = Ybot;   
        Xzstraw[0][istraw-1][0][1] = Xoff1 + 4.7984 + float(istraw-1)*Wdist;
        Xzstraw[1][istraw-1][0][1] = Ybot;
        Xzstraw[0][istraw-1][0][2] = Xoff1 + 4.7984 + float(istraw-1)*Wdist;
        Xzstraw[1][istraw-1][0][2] = Ybot;
// Layer 1 - Top
        Xzstraw[3][istraw-1][0][0] = Xzstraw[0][istraw-1][0][0] + 9.1*Tanthe;
        Xzstraw[4][istraw-1][0][0] = Ytop;
        Xzstraw[3][istraw-1][0][1] = Xzstraw[0][istraw-1][0][1] + 9.1*Tanthe;
        Xzstraw[4][istraw-1][0][1] = Ytop;
        Xzstraw[3][istraw-1][0][2] = Xzstraw[0][istraw-1][0][2] + 9.1*Tanthe;
        Xzstraw[4][istraw-1][0][2] = Ytop;
// Layer 2 - Bottom
        Xzstraw[0][istraw-1][1][0] = Xoff1 + 4.7984 + Whalf + float(istraw-1)*Wdist;
        Xzstraw[1][istraw-1][1][0] = Ybot;
        Xzstraw[0][istraw-1][1][1] = Xoff1 + 4.7984 + Whalf + float(istraw-1)*Wdist;
        Xzstraw[1][istraw-1][1][1] = Ybot;
        Xzstraw[0][istraw-1][1][2] = Xoff1 + 4.7984 + Whalf + float(istraw-1)*Wdist;
        Xzstraw[1][istraw-1][1][2] = Ybot;
// Layer 2 - Top
        Xzstraw[3][istraw-1][1][0] = Xzstraw[0][istraw-1][1][0] + 9.1*Tanthe;
        Xzstraw[4][istraw-1][1][0] = Ytop;
        Xzstraw[3][istraw-1][1][1] = Xzstraw[0][istraw-1][1][1] + 9.1*Tanthe;
        Xzstraw[4][istraw-1][1][1] = Ytop;
        Xzstraw[3][istraw-1][1][2] = Xzstraw[0][istraw-1][1][2] + 9.1*Tanthe;
        Xzstraw[4][istraw-1][1][2] = Ytop;
// Layer 3 - Bottom
        Xzstraw[0][istraw-1][2][0] = Xoff1 + 6.378 + float(istraw-1)*Wdist;
        Xzstraw[1][istraw-1][2][0] = Ybot;
        Xzstraw[0][istraw-1][2][1] = Xoff1 + 6.378 + float(istraw-1)*Wdist;
        Xzstraw[1][istraw-1][2][1] = Ybot;
        Xzstraw[0][istraw-1][2][2] = Xoff1 + 6.378 + float(istraw-1)*Wdist;
        Xzstraw[1][istraw-1][2][2] = Ybot;
// Layer 3 - Top
        Xzstraw[3][istraw-1][2][0] = Xzstraw[0][istraw-1][2][0] - 9.1*Tanthe;
        Xzstraw[4][istraw-1][2][0] = Ytop;
        Xzstraw[3][istraw-1][2][1] = Xzstraw[0][istraw-1][2][1] - 9.1*Tanthe;
        Xzstraw[4][istraw-1][2][1] = Ytop;
        Xzstraw[3][istraw-1][2][2] = Xzstraw[0][istraw-1][2][2] - 9.1*Tanthe;
        Xzstraw[4][istraw-1][2][2] = Ytop;
// Layer 4 - Bottom
        Xzstraw[0][istraw-1][3][0] = Xoff1 + 6.378 - Whalf + float(istraw-1)*Wdist;
        Xzstraw[1][istraw-1][3][0] = Ybot;
        Xzstraw[0][istraw-1][3][1] = Xoff1 + 6.378 - Whalf + float(istraw-1)*Wdist;
        Xzstraw[1][istraw-1][3][1] = Ybot;
        Xzstraw[0][istraw-1][3][2] = Xoff1 + 6.378 - Whalf + float(istraw-1)*Wdist;
        Xzstraw[1][istraw-1][3][2] = Ybot;
// Layer 4 - Top
        Xzstraw[3][istraw-1][3][0] = Xzstraw[0][istraw-1][3][0] - 9.1*Tanthe;
        Xzstraw[4][istraw-1][3][0] = Ytop;
        Xzstraw[3][istraw-1][3][1] = Xzstraw[0][istraw-1][3][1] - 9.1*Tanthe;
        Xzstraw[4][istraw-1][3][1] = Ytop;
        Xzstraw[3][istraw-1][3][2] = Xzstraw[0][istraw-1][3][2] - 9.1*Tanthe;
        Xzstraw[4][istraw-1][3][2] = Ytop;
     }
// ---------------------------------------------------------------------
// Now calculate positions of the wire ends in the u and v frames.
// Z positions are the same in all frames:
// ----------------------------------------------------------------------
     for(int ilay=1; ilay<=4; ilay++) {
        for(int istr=1; istr<=32; istr++) {
           Xzstrau[2][istr-1][ilay-1][0] = Xzstraw[2][istr-1][ilay-1][0];
           Xzstrau[5][istr-1][ilay-1][0] = Xzstraw[5][istr-1][ilay-1][0];
           Xzstrav[2][istr-1][ilay-1][0] = Xzstraw[2][istr-1][ilay-1][0];
           Xzstrav[5][istr-1][ilay-1][0] = Xzstraw[5][istr-1][ilay-1][0];
           Xzstrau[2][istr-1][ilay-1][1] = Xzstraw[2][istr-1][ilay-1][1];
           Xzstrau[5][istr-1][ilay-1][1] = Xzstraw[5][istr-1][ilay-1][1];
           Xzstrav[2][istr-1][ilay-1][1] = Xzstraw[2][istr-1][ilay-1][1];
           Xzstrav[5][istr-1][ilay-1][1] = Xzstraw[5][istr-1][ilay-1][1];
           Xzstrau[2][istr-1][ilay-1][2] = Xzstraw[2][istr-1][ilay-1][2];
           Xzstrau[5][istr-1][ilay-1][2] = Xzstraw[5][istr-1][ilay-1][2];
           Xzstrav[2][istr-1][ilay-1][2] = Xzstraw[2][istr-1][ilay-1][2];
           Xzstrav[5][istr-1][ilay-1][2] = Xzstraw[5][istr-1][ilay-1][2];
        }
     }
// ----------------------------------------------------------------------
     float xx = 0.0;
     float yy = 0.0;
     float xu = 0.0;
     float yu = 0.0;
     float xv = 0.0;
     float yv = 0.0;
     for(int ilay=1; ilay<=4; ilay++) {
        for(int istr=1; istr<=32; istr++) {
// Straw bottom in u frame:  
           xx = Xzstraw[0][istr-1][ilay-1][0];
           yy = Xzstraw[1][istr-1][ilay-1][0];
           Bkxy2u(xx,yy,xu,yu);
           Xzstrau[0][istr-1][ilay-1][0] = xu;
           Xzstrau[1][istr-1][ilay-1][0] = yu;
           xx = Xzstraw[0][istr-1][ilay-1][1];
           yy = Xzstraw[1][istr-1][ilay-1][1];
           Bkxy2u(xx,yy,xu,yu);
           Xzstrau[0][istr-1][ilay-1][1] = xu;
           Xzstrau[1][istr-1][ilay-1][1] = yu;
           xx = Xzstraw[0][istr-1][ilay-1][2];
           yy = Xzstraw[1][istr-1][ilay-1][2];
           Bkxy2u(xx,yy,xu,yu);
           Xzstrau[0][istr-1][ilay-1][2] = xu;
           Xzstrau[1][istr-1][ilay-1][2] = yu;
// Straw top in u frame:
           xx = Xzstraw[3][istr-1][ilay-1][0];
           yy = Xzstraw[4][istr-1][ilay-1][0];
           Bkxy2u(xx,yy,xu,yu);
           Xzstrau[3][istr-1][ilay-1][0] = xu;
           Xzstrau[4][istr-1][ilay-1][0] = yu;
           xx = Xzstraw[3][istr-1][ilay-1][1];
           yy = Xzstraw[4][istr-1][ilay-1][1];
           Bkxy2u(xx,yy,xu,yu);
           Xzstrau[3][istr-1][ilay-1][1] = xu;
           Xzstrau[4][istr-1][ilay-1][1] = yu;
           xx = Xzstraw[3][istr-1][ilay-1][2];
           yy = Xzstraw[4][istr-1][ilay-1][2];
           Bkxy2u(xx,yy,xu,yu);
           Xzstrau[3][istr-1][ilay-1][2] = xu;
           Xzstrau[4][istr-1][ilay-1][2] = yu;
// Straw bottom in v frame:
           xx = Xzstraw[0][istr-1][ilay-1][0];
           yy = Xzstraw[1][istr-1][ilay-1][0];
           Bkxy2v(xx,yy,xv,yv);
           Xzstrav[0][istr-1][ilay-1][0] = xv;
           Xzstrav[1][istr-1][ilay-1][0] = yv;
           xx = Xzstraw[0][istr-1][ilay-1][1];
           yy = Xzstraw[1][istr-1][ilay-1][1];
           Bkxy2v(xx,yy,xv,yv);
           Xzstrav[0][istr-1][ilay-1][1] = xv;
           Xzstrav[1][istr-1][ilay-1][1] = yv;
           xx = Xzstraw[0][istr-1][ilay-1][2];
           yy = Xzstraw[1][istr-1][ilay-1][2];
           Bkxy2v(xx,yy,xv,yv);
           Xzstrav[0][istr-1][ilay-1][2] = xv;
           Xzstrav[1][istr-1][ilay-1][2] = yv;
// Straw top in v frame:
           xx = Xzstraw[3][istr-1][ilay-1][0];
           yy = Xzstraw[4][istr-1][ilay-1][0];
           Bkxy2v(xx,yy,xv,yv);
           Xzstrav[3][istr-1][ilay-1][0] = xv;
           Xzstrav[4][istr-1][ilay-1][0] = yv;
           xx = Xzstraw[3][istr-1][ilay-1][1];
           yy = Xzstraw[4][istr-1][ilay-1][1];
           Bkxy2v(xx,yy,xv,yv);
           Xzstrav[3][istr-1][ilay-1][1] = xv;
           Xzstrav[4][istr-1][ilay-1][1] = yv;
           xx = Xzstraw[3][istr-1][ilay-1][2];
           yy = Xzstraw[4][istr-1][ilay-1][2];
           Bkxy2v(xx,yy,xv,yv);
           Xzstrav[3][istr-1][ilay-1][2] = xv;
           Xzstrav[4][istr-1][ilay-1][2] = yv;
        }
     }
//for (int imod = 0; imod < 3; imod++){
//for (int ilay = 0; ilay < 4; ilay++){
//for (int istr = 0; istr < 32; istr++){
//for (int iindex = 0; iindex < 6;iindex++){
//printf("%f", Xzstraw[iindex][istr][ilay][imod]);
//}
//printf("hello\n");
//}
//}
//}
// ----------------------------------------------------------
// Now get all the calibration constants for the STRAWs.
// ----------------------------------------------------------
     Strcal();
}
// ----------------------------------------------------------
void Strevt(int &Goodt0) {
// ----------------------------------------------------------
// Starting from the Raw time values ... first subtract
// pedestals and time of flight.
// Then try to determine t0 for this event:
// Strevt assumes a hit in each layer has already been
// chosen.
// Goodt0 = 1 if a reliable t0 is obtained.
// ----------------------------------------------------------
      float t0est[6];
      float dt1[3];
      float dt2[3];
      float dt3[3];
      float dt4[3];
      float is1[3];
      float is2[3];
      float is3[3];
      float is4[3];
//
      for(int imod=1; imod<=3; imod++) {
         dt1[imod-1] = -999.0;
         dt2[imod-1] = -999.0;
         dt3[imod-1] = -999.0;
         dt4[imod-1] = -999.0;
         is1[imod-1] = -999;
         is2[imod-1] = -999;
         is3[imod-1] = -999;
         is4[imod-1] = -999;
      }
      float Dvel   = 4.7/1000.0;
// ----------------------------------------------------------
// Once this routine has been called there should be at most
// one hit in each layer - raw time values in Rtimes array.
// ----------------------------------------------------------
      for(int imod=1; imod<=3; imod++) {
         for(int ilay=1; ilay<=4; ilay++) {
            for(int istr=1; istr<=32; istr++) {
               Tof[istr-1][ilay-1][imod-1] = Tofl;
               if(Rtimes[istr-1][ilay-1][imod-1] > 0.0) {
                 Dtimes[istr-1][ilay-1][imod-1] = Rtimes[istr-1][ilay-1][imod-1] - Pedest[istr-1][ilay-1][imod-1] - Tofl;
                 if(Dtimes[istr-1][ilay-1][imod-1] <= 0.0) {
                                                Dtimes[istr-1][ilay-1][imod-1] = 1.0;
                 }
               }
            }
         }
      }
// ------------------------------------------------------------------
// Now try to determine the t0 for this event:
// ------------------------------------------------------------------
      T0data = 0.0;
//
      for(int imod=1; imod<=3; imod++) {
         for(int ist=1; ist<=32; ist++) {
            if(Dtimes[ist-1][0][imod-1] > 0.0) {
                 dt1[imod-1] = Dtimes[ist-1][0][imod-1];
                 is1[imod-1] = ist;
            }
            if(Dtimes[ist-1][1][imod-1] > 0.0) {
                 dt2[imod-1] = Dtimes[ist-1][1][imod-1];
                 is2[imod-1] = ist;
            }
            if(Dtimes[ist-1][2][imod-1] > 0.0) {
                 dt3[imod-1] = Dtimes[ist-1][2][imod-1];
                 is3[imod-1] = ist;
            }
            if(Dtimes[ist-1][3][imod-1] > 0.0) {
                 dt4[imod-1] = Dtimes[ist-1][3][imod-1];
                 is4[imod-1] = ist;
            }
         }
      }
// ------------------------------------------------------------------
// We now have the times in each of the 4 layers in each module.
// Use them to obtain 2 estimates of t0:
// Use the fact that we know the sum of the drift times in adjacent
// layers.
// ------------------------------------------------------------------
      int Nest   = 0;
      float t0aver = 0.0;
      for(int imod=1; imod<=3; imod++) {
         if(abs(is1[imod-1]-is2[imod-1]) < 2) {
           if(dt1[imod-1] > 0.0 && dt2[imod-1] > 0.0) {
             Nest = Nest + 1;
             t0est[Nest-1] = (dt1[imod-1] + dt2[imod-1] - 59.0)/2.0;
             t0aver = t0aver + t0est[Nest-1];
           }
         }
         if(abs(is3[imod-1]-is4[imod-1]) < 2) {
           if(dt3[imod-1] > 0.0 && dt4[imod-1] > 0.0) {
             Nest = Nest + 1;
             t0est[Nest-1] = (dt3[imod-1] + dt4[imod-1] - 59.0)/2.0;
             t0aver = t0aver + t0est[Nest-1];
           }
         }
      }
//
//    cout<<" Strevt: No of T0 estimates: "<<Nest<<endl;
//    for(int iest=1; iest<=Nest; iest++) {
//       cout<<iest<<" "<<t0est[iest-1]<<endl;
//    }
      float t0diff = t0est[0] - t0est[1];
//    cout<<" Difference in 1st two estimates: "<<t0diff<<endl;
//    int ibl4;
//    cin>>ibl4;
      Goodt0 = 0;
      if(Nest == 2 && fabs(t0diff) < 5.0) {
                                          Goodt0 = 1;
      }
      if(Nest > 0) {
                   t0aver = t0aver/float(Nest);
      }
      T0data = t0aver;
//    cout<<" Estimated T0 is: "<<T0data<<endl;
// ------------------------------------------------------------------
// Subtract the obtained t0 from the drift times:
// ------------------------------------------------------------------
      for(int imod=1; imod<=3; imod++) {
         for(int ilay=1; ilay<=4; ilay++) {
            for(int istr=1; istr<=32; istr++) {
               if(Dtimes[istr-1][ilay-1][imod-1] > 0.0) {
                 Dtimes[istr-1][ilay-1][imod-1] = Dtimes[istr-1][ilay-1][imod-1] - T0data;
                 Dtimes[istr-1][ilay-1][imod-1] = Dtimes[istr-1][ilay-1][imod-1] - 4.0;
                 if(Dtimes[istr-1][ilay-1][imod-1] <= 0.0) {
                         Dtimes[istr-1][ilay-1][imod-1] = 1.0;
                 }
                 float Time = Dtimes[istr-1][ilay-1][imod-1];
                 float Dist = 0.0;
                 Sttime3(Time,Dist);
//               float Dist = Time*Dvel;
                 if(Dist > 0.25) {
                              Dist = 0.25;
                 }
                 Hits[istr-1][ilay-1][imod-1] = Dist;
//               cout<<" Drift info: "<<istr<<" "<<ilay<<" "<<imod<<" "<<Time<<" "<<Dist<<endl;
               }
            }
         }
      }     
}
// --------------------------------------------------------------------
void Bkwrte() {
// --------------------------------------------------------------------
// Write each event out for further study.
// --------------------------------------------------------------------
   char StrFile[7];
   float pibk = 3.1415926535;
   float Twopi = 2.0*pibk;
//    
   int Jstr   = 2;
   sprintf(StrFile,"gm%i%s",Jstr,".var");
//    
   ofstream pdffile(StrFile);
   pdffile<<T0data<<endl;   
//
//   int ilocf = 0;
//   for(int ia=1; ia<=6; ia++) {
//      float a1   = Trkinf[ilocf];
//      float a2   = Trkinf[ilocf+1];
//      float a3   = Trkinf[ilocf+2];
//      float a4   = Trkinf[ilocf+3];
//      float a5   = Trkinf[ilocf+4];
//      float a6   = Trkinf[ilocf+5];
//      float a7   = Trkinf[ilocf+6];
//      float a8   = Trkinf[ilocf+7];
//      float a9   = Trkinf[ilocf+8];
//      float a10  = Trkinf[ilocf+9];
//      pdffile<<a1<<" "<<a2<<" "<<a3<<" "<<a4<<" "<<a5<<" "<<a6<<" "<<a7<<" "<<a8<<" "<<a9<<" "<<a10<<endl;
//      ilocf = ilocf + 10;
//   }
//
//   int iloc2 = 0;
//   for(int ib=1; ib<=6; ib++) {
//      float b1   = Trkin2[iloc2];            
//      float b2   = Trkin2[iloc2+1];          
//      float b3   = Trkin2[iloc2+2];          
//      float b4   = Trkin2[iloc2+3];          
//      float b5   = Trkin2[iloc2+4];          
//      float b6   = Trkin2[iloc2+5];          
//      float b7   = Trkin2[iloc2+6];          
//      float b8   = Trkin2[iloc2+7];          
//      float b9   = Trkin2[iloc2+8];
//      float b10  = Trkin2[iloc2+9];
//      pdffile<<b1<<" "<<b2<<" "<<b3<<" "<<b4<<" "<<b5<<" "<<b6<<" "<<b7<<" "<<b8<<" "<<b9<<" "<<b10<<endl;
//      iloc2 = iloc2 + 10;
//   }  
//
//   int iloc3 = 0;
//   for(int ic=1; ic<=6; ic++) {
//      float c1   = Trkin3[iloc3];
//      float c2   = Trkin3[iloc3+1];
//      float c3   = Trkin3[iloc3+2];
//      float c4   = Trkin3[iloc3+3];
//      float c5   = Trkin3[iloc3+4];
//      float c6   = Trkin3[iloc3+5];
//      float c7   = Trkin3[iloc3+6];
//      float c8   = Trkin3[iloc3+7];
//      float c9   = Trkin3[iloc3+8];
//      float c10  = Trkin3[iloc3+9];
//      pdffile<<c1<<" "<<c2<<" "<<c3<<" "<<c4<<" "<<c5<<" "<<c6<<" "<<c7<<" "<<c8<<" "<<c9<<" "<<c10<<endl;
//      iloc3 = iloc3 + 10;
//   }  
//
//      int Ntot = 0;
//      for(int istr=1; istr<=32; istr++) {
//         for(int ilay=1; ilay<=4; ilay++) {
//            if(Hits[istr-1][ilay-1] > 0.0) {
//                                           Ntot = Ntot + 1;
//            }
//         }
//      }
//      pdffile<<Ntot<<endl;
//
//      for(int istr=1; istr<=32; istr++) {
//         for(int ilay=1; ilay<=4; ilay++) {
//            if(Hits[istr-1][ilay-1] > 0.0) {
//              float aa1 = Hits[istr-1][ilay-1];
//              float aa2 = Dtimes[istr-1][ilay-1];
//              float aa3 = Rtimes[istr-1][ilay-1];
//              float aa4 = Yhits(Istraw,Ilay,Icol)
//              iaa = Mask(Istraw,Ilay,Icol)
//        aa4 = Zxhit(1,Istraw,Ilay,Icol)
//        aa5 = Zxhit(2,Istraw,Ilay,Icol)
//        aa6 = Xzaty(1,Istraw,Ilay,Icol)
//        aa7 = Xzaty(2,Istraw,Ilay,Icol)
//        bb1 = Hitsmc(Istraw,Ilay,Icol)
//        bb2 = Dtimmc(Istraw,Ilay,Icol)
//        bb3 = Yhitmc(Istraw,Ilay,Icol)
//        bb4 = Zxhmc(1,Istraw,Ilay,Icol)
//        bb5 = Zxhmc(2,Istraw,Ilay,Icol)
//        Write(26,1004) Istraw,Ilay,Icol,iaa,aa1,aa2,aa3,aa4,aa5, 
//     +                 aa6,aa7
//        Write(26,1005) bb1,bb2,bb3,bb4,bb5
//        If(Icall .eq. 1) Then
//          Write(61,*) ' Reco Hit in ',Istraw,Ilay,Icol
//          Write(61,*)
//     +    Dtimes(Istraw,Ilay,Icol),Hits(Istraw,Ilay,Icol)
//          Write(61,*)
//     +    Yhits(Istraw,Ilay,Icol),T0data
//          Write(61,*)
//     +    Zxhit(1,Istraw,Ilay,Icol),Zxhit(2,Istraw,Ilay,Icol)
//        Endif
//      Endif
//   22 Continue
//   21 Continue
}
// ------------------------------------------------------------------
void Strft4(int Nmod) {
// ------------------------------------------------------------------
// Finds a 4-hit TE in a Straw module.
// This routine deals with modules with 4 hits on a track.
// ------------------------------------------------------------------
      float Xptsfit[4];
      float Zptsfit[4];
      float Zw[4];
      float Xw[4];
      float Dr[4];
      float Resols[4000];
      int   Nused[4000];
      int   Isols[4000];
//
      float Pibk  = 3.1415926535;
      float Twopi = 2.0*Pibk;
      float Piby2 = Pibk/2.0;
// ------------------------------
// Size of vertical stepping bin:
// ------------------------------
      float Ysize = 0.01;
// ------------------------------
// Total number of vertical steps
// ------------------------------
      float Ystep = 10.0/Ysize;
// ------------------------------------------
// Size of array needed to store information:
// ------------------------------------------
      float Yloop = 4.0 * Ystep;
      int Kystep = Bknint(Ystep);
      int Kyloop = Bknint(Yloop);
//
      for(int i=1; i<=Kyloop; i++) {
         Resols[i-1] = 99999999.0;
         Nused[i-1]  = 0;
         Isols[i-1]  = 0;
      }
// -----------------------------------------
// Reset the number of TE points:
// -----------------------------------------
      Nposte = 0;
// ------------------------------------------------------------------
// Declare some variables for later use.
// ------------------------------------------------------------------
      float Xbot = 0.0;
      float Ybot = 0.0;
      float Zbot = 0.0;
      float Xtop = 0.0;
      float Ytop = 0.0;
      float Ztop = 0.0;
      float xx1  = 0.0;
      float xx2  = 0.0;
      float xx3  = 0.0;
      float xx4  = 0.0;
      float zz1  = 0.0;
      float zz2  = 0.0;
      float zz3  = 0.0;
      float zz4  = 0.0;
      float yy1  = 0.0;
      float yy2  = 0.0;
      float yy3  = 0.0;
      float yy4  = 0.0;
      float xu1  = 0.0;
      float xu2  = 0.0;
      float xu3  = 0.0;
      float xu4  = 0.0;
      float zu1  = 0.0;
      float zu2  = 0.0;
      float zu3  = 0.0;
      float zu4  = 0.0;
      float yu1  = 0.0;
      float yu2  = 0.0;
      float yu3  = 0.0;
      float yu4  = 0.0;
//
      float Gradu = 0.0;
      float Cintu = 0.0;
      float Gradv = 0.0;
      float Cintv = 0.0;
      float xv1  = 0.0;
      float xv2  = 0.0;
      float xv3  = 0.0;
      float xv4  = 0.0;
      float zv1  = 0.0;
      float zv2  = 0.0;
      float zv3  = 0.0;
      float zv4  = 0.0;
      float yv1  = 0.0;
      float yv2  = 0.0;
      float yv3  = 0.0;
      float yv4  = 0.0;
// ------------------------------------------------------------------
// Find the hits in each of the layers:
// ------------------------------------------------------------------
      int Layer1 = 1;
      int Layer2 = 2;
      int Layer3 = 3;
      int Layer4 = 4;
      int Ist1 = 0;
      int Ist2 = 0;
      int Ist3 = 0;
      int Ist4 = 0;
      int Nhits = 0;
      for(int istr=1; istr<=32; istr++) {
         if(Hits[istr-1][0][Nmod-1] > 0.0) {
                                   Ist1  = istr;
                                   Nhits = Nhits + 1;
         }
         if(Hits[istr-1][1][Nmod-1] > 0.0) {
                                   Ist2  = istr;
                                   Nhits = Nhits + 1;
         }
         if(Hits[istr-1][2][Nmod-1] > 0.0) {
                                   Ist3  = istr;
                                   Nhits = Nhits + 1;
         }
         if(Hits[istr-1][3][Nmod-1] > 0.0) {
                                   Ist4  = istr;
                                   Nhits = Nhits + 1;
         }
      }
// ------------------------------------------------------------------
// If we dont have hits in all four layers then return.
// ------------------------------------------------------------------
      if(Nhits != 4) {
                     return;
      }
      if(Ist1 == 0 || Ist2 == 0) {
                                 return;
      }
      if(Ist3 == 0 || Ist4 == 0) {
                                 return;
      }
// ------------------------------------------------------------------
// In a module with 4 trackable hits: Make some basic plots.
// ------------------------------------------------------------------
      float diff12 = Rtimes[Ist1-1][0][Nmod-1] - Rtimes[Ist2-1][1][Nmod-1];
      float diff34 = Rtimes[Ist3-1][2][Nmod-1] - Rtimes[Ist4-1][3][Nmod-1];
      StrHist.fdiff12[0]->Fill(diff12);
      StrHist.fdiff34[0]->Fill(diff34);
      StrHist.fDtimes[0]->Fill(Dtimes[Ist1-1][0][Nmod-1]);
      StrHist.fDtimes[0]->Fill(Dtimes[Ist2-1][1][Nmod-1]);
      StrHist.fDtimes[0]->Fill(Dtimes[Ist3-1][2][Nmod-1]);
      StrHist.fDtimes[0]->Fill(Dtimes[Ist4-1][3][Nmod-1]);
      StrHist.fDdists[0]->Fill(Hits[Ist1-1][0][Nmod-1]);
      StrHist.fDdists[0]->Fill(Hits[Ist2-1][1][Nmod-1]);
      StrHist.fDdists[0]->Fill(Hits[Ist3-1][2][Nmod-1]);
      StrHist.fDdists[0]->Fill(Hits[Ist4-1][3][Nmod-1]);
// ------------------------------------------------------------------
// Do track finding at various positions vertically and minimise the
// residuals:
// ------------------------------------------------------------------
      for (int Iybin=1; Iybin <= 910; Iybin++) {
          float Ytest = float(Iybin)/100.0;
          Ytest = Ytest - 4.55;
          for(int ilay=1; ilay<=4; ilay++) {
             for(int istr=1; istr<=32; istr++) {
                float Zendc = Xzstraw[2][istr-1][ilay-1][Nmod-1];
                float Xendc = Xzstraw[0][istr-1][ilay-1][Nmod-1];
                float Zenda = Xzstraw[5][istr-1][ilay-1][Nmod-1];
                float Xenda = Xzstraw[3][istr-1][ilay-1][Nmod-1];
// ------------------------------------------------------------
// Work out wire position at the vertical position (y) at which
// the track finding will take place:
// ------------------------------------------------------------
// Old calculation in the xyz frame:
                float Zaty  = Zendc;
                float Xaty  = Xendc + ((Xenda-Xendc)*(Ytest -(-4.55))/(9.1));
// ------------------------------------
// New calculation in the u or v frame:
// ------------------------------------
                if(ilay <= 2) {
                  Xbot = Xzstrau[0][istr-1][ilay-1][Nmod-1];
                  Ybot = Xzstrau[1][istr-1][ilay-1][Nmod-1];
                  Zbot = Xzstrau[2][istr-1][ilay-1][Nmod-1];
                  Xtop = Xzstrau[3][istr-1][ilay-1][Nmod-1];
                  Ytop = Xzstrau[4][istr-1][ilay-1][Nmod-1];
                  Ztop = Xzstrau[5][istr-1][ilay-1][Nmod-1];
                }
                if(ilay > 2) {
                  Xbot = Xzstrav[0][istr-1][ilay-1][Nmod-1];
                  Ybot = Xzstrav[1][istr-1][ilay-1][Nmod-1];
                  Zbot = Xzstrav[2][istr-1][ilay-1][Nmod-1];
                  Xtop = Xzstrav[3][istr-1][ilay-1][Nmod-1];
                  Ytop = Xzstrav[4][istr-1][ilay-1][Nmod-1];
                  Ztop = Xzstrav[5][istr-1][ilay-1][Nmod-1];
                }
// ------------------------------------------------------------------
// So Xzaty contains the z/x wire positions in the u frame for
// Layers 1 and 2. And in the v frame for Layers 3 and 4.
// ------------------------------------------------------------------
                Xzaty[0][istr-1][ilay-1][Nmod-1] = Ztop;
                Xzaty[1][istr-1][ilay-1][Nmod-1] = Xtop;
             }
          }
// ------------------------------------------------------------------
// Have found the hit straws in layers 1 and 2.
// Find the 4 possible tangents to these two drift circles.
// Extrapolate to layers 3 and 4 ...
// ------------------------------------------------------------------
// In u frame:
          Zw[0] = Xzaty[0][Ist1-1][0][Nmod-1];
          Xw[0] = Xzaty[1][Ist1-1][0][Nmod-1];
          Dr[0] = Hits[Ist1-1][0][Nmod-1];
          Ybot  = Xzstrau[1][Ist1-1][0][Nmod-1];
          Ytop  = Xzstrau[4][Ist1-1][0][Nmod-1];
          yu1   = Ybot + ((Ytop-Ybot)*(Ytest -(-4.55))/(9.1));
//
          Zw[1] = Xzaty[0][Ist2-1][1][Nmod-1];
          Xw[1] = Xzaty[1][Ist2-1][1][Nmod-1];
          Dr[1] = Hits[Ist2-1][1][Nmod-1];
          Ybot  = Xzstrau[1][Ist2-1][1][Nmod-1];
          Ytop  = Xzstrau[4][Ist2-1][1][Nmod-1];
          yu2   = Ybot + ((Ytop-Ybot)*(Ytest -(-4.55))/(9.1));
// In v frame:
          Zw[2] = Xzaty[0][Ist3-1][2][Nmod-1];
          Xw[2] = Xzaty[1][Ist3-1][2][Nmod-1];
          Dr[2] = Hits[Ist3-1][2][Nmod-1];
          Ybot  = Xzstrav[1][Ist3-1][2][Nmod-1];
          Ytop  = Xzstrav[4][Ist3-1][2][Nmod-1];
          yv3   = Ybot + ((Ytop-Ybot)*(Ytest -(-4.55))/(9.1));
// 
          Zw[3] = Xzaty[0][Ist4-1][3][Nmod-1];
          Xw[3] = Xzaty[1][Ist4-1][3][Nmod-1];
          Dr[3] = Hits[Ist4-1][3][Nmod-1];
          Ybot  = Xzstrav[1][Ist4-1][3][Nmod-1];
          Ytop  = Xzstrav[4][Ist4-1][3][Nmod-1];
          yv4   = Ybot + ((Ytop-Ybot)*(Ytest -(-4.55))/(9.1));
//
          for(int itan=1; itan<=4; itan++) {
// ------------------------------------------------------------------
// Routine Strelm returns the Z and X points of the common
// tangents to the drift circles for each of the 4 cases in turn.
// ------------------------------------------------------------------
//   Itan = 1 MEANS OUTER EDGE TANGENT , RHS OF DRIFT CIRCLES
//        = 2        ''    ''     ''     LHS  ''  ''     ''
//        = 3       CROSSED TANGENT , BOTTOM RHS TO TOP LHS
//        = 4           ''     ''       ''   LHS        RHS
// ------------------------------------------------------------------
             zu1 = 0.0;
             xu1 = 0.0;
             zu2 = 0.0;
             xu2 = 0.0;
             int iflag = 0;
             Strelm(Nmod,Ist1,Layer1,Ist2,Layer2,itan,zu1,xu1,zu2,xu2,iflag);
             int Irloc = (Iybin-1)*4 + itan;
             Gradu = 9999.9;
             int Icon = 1;
             if(fabs(zu2-zu1) > 1.0e-6) {
                                        Gradu =(xu2 - xu1)/(zu2 - zu1);
                                        Icon = 1;
             }
             if(fabs(zu2-zu1) <= 1.0e-6) {
                                         Resols[Irloc-1] = 99999.9;
                                         Icon = 0;
             }
             if(Icon == 1) {
               Cintu = xu1 - (Gradu * zu1);
// -------------------------------------------------------------------
// Sanity check: Residuals in layers used to define the line should
// be zero.
// -------------------------------------------------------------------
               float Zwire1u = Zw[0];
               float Xwire1u = Xw[0];
               float Zint1u  = 999.9;
               float Xint1u  = 999.9;
               Strclp(Gradu,Cintu,Zwire1u,Xwire1u,Zint1u,Xint1u);
               float Zwire2u = Zw[1];
               float Xwire2u = Xw[1];
               float Zint2u  = 999.9;
               float Xint2u  = 999.9;
               Strclp(Gradu,Cintu,Zwire2u,Xwire2u,Zint2u,Xint2u);
               float Zpt1u   = 999.9;
               float Xpt1u   = 999.9;  
               Strcht(zu1,xu1,zu2,xu2,Ist1,Layer1,Nmod,Zpt1u,Xpt1u);
               float Res1u = sqrt((Zpt1u - Zint1u)*(Zpt1u - Zint1u) + (Xpt1u - Xint1u)*(Xpt1u - Xint1u));
               Res1u = Res1u*10000.0;
               float Zpt2u   = 999.9;
               float Xpt2u   = 999.9;
               Strcht(zu1,xu1,zu2,xu2,Ist2,Layer2,Nmod,Zpt2u,Xpt2u);
               float Res2u = sqrt((Zpt2u - Zint2u)*(Zpt2u - Zint2u) + (Xpt2u - Xint2u)*(Xpt2u - Xint2u));
               Res2u = Res2u*10000.0;
// -------------------------------------------------------------------
// Find closest point to this line on the hit cells in layers 3 and 4
// We need to swap the extrapolated tangent into the v frame from the
// u frame: 
// -------------------------------------------------------------------
               zv1 = zu1;
               Bku2xy(xu1,yu1,xx1,yy1);
               Bkxy2v(xx1,yy1,xv1,yv1);
//
               zv2 = zu2;
               Bku2xy(xu2,yu2,xx2,yy2);
               Bkxy2v(xx2,yy2,xv2,yv2);
//
               Gradv =(xv2 - xv1)/(zv2 - zv1);
               Cintv = xv1 - (Gradv * zv1);
//
               float Zwire3v = Zw[2];
               float Xwire3v = Xw[2];
               float Zint3v  = 999.9;
               float Xint3v  = 999.9;
               Strclp(Gradv,Cintv,Zwire3v,Xwire3v,Zint3v,Xint3v);
               float Zwire4v = Zw[3];
               float Xwire4v = Xw[3];
               float Zint4v  = 999.9;
               float Xint4v  = 999.9;
               Strclp(Gradv,Cintv,Zwire4v,Xwire4v,Zint4v,Xint4v);
               float Dist3v = sqrt((Zwire3v-Zint3v)*(Zwire3v-Zint3v) + (Xwire3v-Xint3v)*(Xwire3v-Xint3v)); 
               int Icon3 = 1;
               if(Dist3v > 1.2*Size[2]) {
                             Resols[Irloc-1] = 99999.9;
                             Icon3 = 0;
               }
               float Dist4v = sqrt((Zwire4v-Zint4v)*(Zwire4v-Zint4v) + (Xwire4v-Xint4v)*(Xwire4v-Xint4v));
               int Icon4 = 1;
               if(Dist4v > 1.2*Size[3]) {
                             Resols[Irloc-1] = 99999.9;
                             Icon4 = 0;
               }
// -------------------------------------------------------------------
// Extrapolated line passes through the hit cells in layers 3 and 4:
// Calculate the residuals:
// -------------------------------------------------------------------
               if(Icon3 == 1 && Icon4 == 1) {
                 float Zpt3v = 999.9;
                 float Xpt3v = 999.9;
                 Strcht(zv1,xv1,zv2,xv2,Ist3,Layer3,Nmod,Zpt3v,Xpt3v);
                 float Res3v = sqrt((Zpt3v - Zint3v)*(Zpt3v - Zint3v) + (Xpt3v - Xint3v)*(Xpt3v - Xint3v));
                 Res3v = Res3v*10000.0;
                 float Zpt4v = 999.9;
                 float Xpt4v = 999.9;
                 Strcht(zv1,xv1,zv2,xv2,Ist4,Layer4,Nmod,Zpt4v,Xpt4v);
                 float Res4v = sqrt((Zpt4v - Zint4v)*(Zpt4v - Zint4v) + (Xpt4v - Xint4v)*(Xpt4v - Xint4v));
                 Res4v = Res4v*10000.0;
                 float Restot = Res3v + Res4v;
                 Resols[Irloc-1] = Restot;
               } 
             }
          }
      }
// ------------------------------------------------------------------
// We now have the residual information for extrapolation of each
// of the 4 tangents at each test Y position ... choose the best.
// ------------------------------------------------------------------
      int Locbst = 0;
      float Resbst = 99999.9;
      for(int isol=1; isol<=4000; isol++) {
         float Resid = Resols[isol-1];
         if(Resid < Resbst) {
                            Resbst = Resid;
                            Locbst = isol;
         }
      }
// ------------------------------------------------------------------
// Best chosen .... now refit it:
// ------------------------------------------------------------------
      int Iybest   = int((Locbst-1)/4) + 1;
      int Itnbst   = Locbst - ((Iybest-1)*4);
      float Ybest  = float(Iybest)/100.0;
      Ybest = Ybest - 4.55;
//
      for(int istr=1; istr<=32; istr++) {
         for(int ilay=1; ilay<=4; ilay++) {
            float Zendcb = Xzstraw[2][istr-1][ilay-1][Nmod-1];
            float Xendcb = Xzstraw[0][istr-1][ilay-1][Nmod-1];
            float Zendab = Xzstraw[5][istr-1][ilay-1][Nmod-1];
            float Xendab = Xzstraw[3][istr-1][ilay-1][Nmod-1];
// ------------------------------------------------------------
// Work out wire position at the vertical position (y) at which
// the best track residuals were found:
// ------------------------------------------------------------
            float Zatyb  = Zendcb;
            float Xatyb  = Xendcb + ((Xendab-Xendcb)*(Ybest -(-4.55))/(9.1));
// ------------------------------------
// New calculation in the u or v frame:
// ------------------------------------   
            if(ilay <= 2) {
              Xbot = Xzstrau[0][istr-1][ilay-1][Nmod-1];
              Ybot = Xzstrau[1][istr-1][ilay-1][Nmod-1];
              Zbot = Xzstrau[2][istr-1][ilay-1][Nmod-1];
              Xtop = Xzstrau[3][istr-1][ilay-1][Nmod-1];
              Ytop = Xzstrau[4][istr-1][ilay-1][Nmod-1];
              Ztop = Xzstrau[5][istr-1][ilay-1][Nmod-1];
            }
            if(ilay >= 3) {
              Xbot = Xzstrav[0][istr-1][ilay-1][Nmod-1];
              Ybot = Xzstrav[1][istr-1][ilay-1][Nmod-1];
              Zbot = Xzstrav[2][istr-1][ilay-1][Nmod-1];
              Xtop = Xzstrav[3][istr-1][ilay-1][Nmod-1];
              Ytop = Xzstrav[4][istr-1][ilay-1][Nmod-1];
              Ztop = Xzstrav[5][istr-1][ilay-1][Nmod-1];
            }
// ------------------------------------------------------------------
// So Xzaty contains the z/x wire positions in the u frame for
// Layers 1 and 2. And in the v frame for Layers 3 and 4.
// ------------------------------------------------------------------
            Xzaty[0][istr-1][ilay-1][Nmod-1] = Ztop;
            Xzaty[1][istr-1][ilay-1][Nmod-1] = Xtop;
         }
      }
//
      Zw[0] = Xzaty[0][Ist1-1][0][Nmod-1];
      Xw[0] = Xzaty[1][Ist1-1][0][Nmod-1];
      Dr[0] = Hits[Ist1-1][0][Nmod-1];
      Ybot  = Xzstrau[1][Ist1-1][0][Nmod-1];
      Ytop  = Xzstrau[4][Ist1-1][0][Nmod-1];
      yu1   = Ybot + ((Ytop-Ybot)*(Ybest -(-4.55))/(9.1));
//
      Zw[1] = Xzaty[0][Ist2-1][1][Nmod-1];
      Xw[1] = Xzaty[1][Ist2-1][1][Nmod-1];
      Dr[1] = Hits[Ist2-1][1][Nmod-1];
      Ybot  = Xzstrau[1][Ist2-1][1][Nmod-1];
      Ytop  = Xzstrau[4][Ist2-1][1][Nmod-1];
      yu2   = Ybot + ((Ytop-Ybot)*(Ybest -(-4.55))/(9.1));
//
      Zw[2] = Xzaty[0][Ist3-1][2][Nmod-1];
      Xw[2] = Xzaty[1][Ist3-1][2][Nmod-1];
      Dr[2] = Hits[Ist3-1][2][Nmod-1];
      Ybot  = Xzstrav[1][Ist3-1][2][Nmod-1];
      Ytop  = Xzstrav[4][Ist3-1][2][Nmod-1];
      yv3   = Ybot + ((Ytop-Ybot)*(Ybest -(-4.55))/(9.1));
//
      Zw[3] = Xzaty[0][Ist4-1][3][Nmod-1];
      Xw[3] = Xzaty[1][Ist4-1][3][Nmod-1];
      Dr[3] = Hits[Ist4-1][3][Nmod-1];
      Ybot  = Xzstrav[1][Ist4-1][3][Nmod-1];
      Ytop  = Xzstrav[4][Ist4-1][3][Nmod-1];
      yv4   = Ybot + ((Ytop-Ybot)*(Ybest -(-4.55))/(9.1));
// 
      int iflagb = 0;
      Strelm(Nmod,Ist1,Layer1,Ist2,Layer2,Itnbst,zu1,xu1,zu2,xu2,iflagb);
      Gradu = 999.9;
      if(fabs(zu2-zu1) > 1.0e-6) {
                                 Gradu = (xu2 - xu1) / (zu2 - zu1);
      }
      if(fabs(zu2-zu1) <= 1.0e-6) {
                                  return;
      }
      Cintu = xu1 - (Gradu * zu1);
//
      float Zwire1ub = Zw[0];
      float Xwire1ub = Xw[0];
      float Zint1ub  = 999.9;
      float Xint1ub  = 999.9;
      Strclp(Gradu,Cintu,Zwire1ub,Xwire1ub,Zint1ub,Xint1ub);
      float Dist1ub = sqrt((Zwire1ub-Zint1ub)*(Zwire1ub-Zint1ub) + (Xwire1ub-Xint1ub)*(Xwire1ub-Xint1ub));
      if(Dist1ub > 1.2*Size[0]) {
                                return;
      }
      float Zpt1ub = 999.9;
      float Xpt1ub = 999.9;
      Strcht(zu1,xu1,zu2,xu2,Ist1,Layer1,Nmod,Zpt1ub,Xpt1ub);
//
      float Zwire2ub = Zw[1];
      float Xwire2ub = Xw[1];
      float Zint2ub  = 999.9;
      float Xint2ub  = 999.9;
      Strclp(Gradu,Cintu,Zwire2ub,Xwire2ub,Zint2ub,Xint2ub);
      float Dist2ub = sqrt((Zwire2ub-Zint2ub)*(Zwire2ub-Zint2ub) + (Xwire2ub-Xint2ub)*(Xwire2ub-Xint2ub));
      if(Dist2ub > 1.2*Size[1]) {
                                return;
      }
      float Zpt2ub = 999.9;
      float Xpt2ub = 999.9;
      Strcht(zu1,xu1,zu2,xu2,Ist2,Layer2,Nmod,Zpt2ub,Xpt2ub);
// -------------------------------------------------------------------
// Find closest point to this line on the hit cells in layers 3 and 4
// We need to swap the extrapolated tangent into the v frame from the
// u frame:
// -------------------------------------------------------------------
      zv1 = zu1;
      Bku2xy(xu1,yu1,xx1,yy1);
      Bkxy2v(xx1,yy1,xv1,yv1);
//
      zv2 = zu2;
      Bku2xy(xu2,yu2,xx2,yy2);
      Bkxy2v(xx2,yy2,xv2,yv2);
//
      if(fabs(zv2-zv1) > 1.0e-6) {
                                 Gradv =(xv2 - xv1)/(zv2 - zv1);   
      }
      if(fabs(zv2-zv1) <= 1.0e-6) {
                                  return;
      }
      Cintv = xv1 - (Gradv * zv1); 
//
      float Zwire3vb = Zw[2];
      float Xwire3vb = Xw[2];
      float Zint3vb  = 999.9;
      float Xint3vb  = 999.9;
      Strclp(Gradv,Cintv,Zwire3vb,Xwire3vb,Zint3vb,Xint3vb);
      float Dist3vb = sqrt((Zwire3vb-Zint3vb)*(Zwire3vb-Zint3vb) + (Xwire3vb-Xint3vb)*(Xwire3vb-Xint3vb));
      if(Dist3vb > 1.2*Size[2]) {
                                return;
      }
      float Zpt3vb = 999.9;
      float Xpt3vb = 999.9;
      Strcht(zv1,xv1,zv2,xv2,Ist3,Layer3,Nmod,Zpt3vb,Xpt3vb);
      float Zwire4vb = Zw[3];
      float Xwire4vb = Xw[3];
      float Zint4vb  = 999.9;
      float Xint4vb  = 999.9;
      Strclp(Gradv,Cintv,Zwire4vb,Xwire4vb,Zint4vb,Xint4vb);
      float Dist4vb = sqrt((Zwire4vb-Zint4vb)*(Zwire4vb-Zint4vb) + (Xwire4vb-Xint4vb)*(Xwire4vb-Xint4vb));
      if(Dist4vb > 1.2*Size[3]) {
                                return;
      }
      float Zpt4vb = 999.9;
      float Xpt4vb = 999.9;
      Strcht(zv1,xv1,zv2,xv2,Ist4,Layer4,Nmod,Zpt4vb,Xpt4vb);
      float Res1ub = sqrt((Zpt1ub - Zint1ub)*(Zpt1ub - Zint1ub) + (Xpt1ub - Xint1ub)*(Xpt1ub - Xint1ub));
      Res1ub = Res1ub*10000.0;
      float Res2ub = sqrt((Zpt2ub - Zint2ub)*(Zpt2ub - Zint2ub) + (Xpt2ub - Xint2ub)*(Xpt2ub - Xint2ub));
      Res2ub = Res2ub*10000.0;
      float Res3vb = sqrt((Zpt3vb - Zint3vb)*(Zpt3vb - Zint3vb) + (Xpt3vb - Xint3vb)*(Xpt3vb - Xint3vb));
      Res3vb = Res3vb*10000.0;
      float Res4vb = sqrt((Zpt4vb - Zint4vb)*(Zpt4vb - Zint4vb) + (Xpt4vb - Xint4vb)*(Xpt4vb - Xint4vb));
      Res4vb = Res4vb*10000.0;
// -----------------------------------------------------------
// Convert all points back into the x-y frame for fitting.
// -----------------------------------------------------------
      float Xpt1b = 0.0;
      Bku2xy(Xpt1ub,yu1,Xpt1b,yy1);
      float Zpt1b = Zpt1ub;
//
      float Xpt2b = 0.0;
      Bku2xy(Xpt2ub,yu2,Xpt2b,yy2);
      float Zpt2b = Zpt2ub;
//
      float Xpt3b = 0.0;
      Bkv2xy(Xpt3vb,yv3,Xpt3b,yy3);
      float Zpt3b = Zpt3vb;
//
      float Xpt4b = 0.0;    
      Bkv2xy(Xpt4vb,yv4,Xpt4b,yy4);
      float Zpt4b = Zpt4vb;
//
      Strflr(Ist1,Layer1,Nmod,Zpt1b,Xpt1b);
      Strflr(Ist2,Layer2,Nmod,Zpt2b,Xpt2b);
      Strflr(Ist3,Layer3,Nmod,Zpt3b,Xpt3b);
      Strflr(Ist4,Layer4,Nmod,Zpt4b,Xpt4b);
//
      Yhits[Ist1-1][Layer1-1][Nmod-1] = Ybest;
      Yhits[Ist2-1][Layer2-1][Nmod-1] = Ybest;
      Yhits[Ist3-1][Layer3-1][Nmod-1] = Ybest;
      Yhits[Ist4-1][Layer4-1][Nmod-1] = Ybest;
      Mask[Ist1-1][Layer1-1][Nmod-1] = 1;
      Mask[Ist2-1][Layer2-1][Nmod-1] = 1;
      Mask[Ist3-1][Layer3-1][Nmod-1] = 1;
      Mask[Ist4-1][Layer4-1][Nmod-1] = 1;
      Trkinf[0][Nmod-1]  = float(Ist1);
      Trkinf[1][Nmod-1]  = float(Ist2);
      Trkinf[2][Nmod-1]  = float(Ist3);
      Trkinf[3][Nmod-1]  = float(Ist4);
//
      Trkinf[4][Nmod-1]  = Hits[Ist1-1][Layer1-1][Nmod-1];
      Trkinf[5][Nmod-1]  = Zpt1b;
      Trkinf[6][Nmod-1]  = Xpt1b;
      Trkinf[7][Nmod-1]  = Ybest;
      StrHist.fYbest[0]->Fill(Ybest);
      float Zwire1 = Zwire1ub;
      float Xwire1 = 0.0;
      Bku2xy(Xwire1ub,yu1,Xwire1,yy1);
      Trkinf[8][Nmod-1]  = Zwire1;
      Trkinf[9][Nmod-1]  = Xwire1;
      float Zint1 = Zint1ub;
      float Xint1 = 0.0;
      Bku2xy(Xint1ub,yu1,Xint1,yy1);
      Trkinf[10][Nmod-1] = Zint1;
      Trkinf[11][Nmod-1] = Xint1;
//
      Trkinf[12][Nmod-1] = Hits[Ist2-1][Layer2-1][Nmod-1];
      Trkinf[13][Nmod-1] = Zpt2b;
      Trkinf[14][Nmod-1] = Xpt2b;
      Trkinf[15][Nmod-1] = Ybest;
      float Zwire2 = Zwire2ub;
      float Xwire2 = 0.0;
      Bku2xy(Xwire2ub,yu2,Xwire2,yy2);
      float Zint2 = Zint2ub;
      float Xint2 = 0.0;
      Bku2xy(Xint2ub,yu2,Xint2,yy2);
      Trkinf[16][Nmod-1] = Zwire2;
      Trkinf[17][Nmod-1] = Xwire2;
      Trkinf[18][Nmod-1] = Zint2;
      Trkinf[19][Nmod-1] = Xint2;
//
      Trkinf[20][Nmod-1] = Hits[Ist3-1][Layer3-1][Nmod-1];
      Trkinf[21][Nmod-1] = Zpt3b;
      Trkinf[22][Nmod-1] = Xpt3b;
      Trkinf[23][Nmod-1] = Ybest;
      float Zwire3 = Zwire3vb;
      float Xwire3 = 0.0;
      Bkv2xy(Xwire3vb,yv3,Xwire3,yy3);
      float Zint3 = Zint3vb;
      float Xint3 = 0.0;
      Bkv2xy(Xint3vb,yv3,Xint3,yy3);
      Trkinf[24][Nmod-1] = Zwire3;
      Trkinf[25][Nmod-1] = Xwire3;
      Trkinf[26][Nmod-1] = Zint3;
      Trkinf[27][Nmod-1] = Xint3;
//
      Trkinf[28][Nmod-1] = Hits[Ist4-1][Layer4-1][Nmod-1];
      Trkinf[29][Nmod-1] = Zpt4b;
      Trkinf[30][Nmod-1] = Xpt4b;
      Trkinf[31][Nmod-1] = Ybest;
      float Zwire4 = Zwire4vb;
      float Xwire4 = 0.0;
      Bkv2xy(Xwire4vb,yv4,Xwire4,yy4);
      float Zint4 = Zint4vb;
      float Xint4 = 0.0;
      Bkv2xy(Xint4vb,yv4,Xint4,yy4);
      Trkinf[32][Nmod-1] = Zwire4;
      Trkinf[33][Nmod-1] = Xwire4;
      Trkinf[34][Nmod-1] = Zint4;
      Trkinf[35][Nmod-1] = Xint4;
// ------------------------------------------------------------------
// Finally make the fitted TE from the obtained space points.
// Make Xzaty contain wire coordinates in xy frame:
// ------------------------------------------------------------------
      float Zzendc = Xzstraw[2][Ist1-1][0][Nmod-1];
      float Xxendc = Xzstraw[0][Ist1-1][0][Nmod-1];
      float Zzenda = Xzstraw[5][Ist1-1][0][Nmod-1];
      float Xxenda = Xzstraw[3][Ist1-1][0][Nmod-1];
      float Zzaty  = Zzendc;
      float Xxaty  = Xxendc + ((Xxenda-Xxendc)*(Ybest -(-4.55))/(9.1));
      Xzaty[0][Ist1-1][0][Nmod-1] = Zzaty;
      Xzaty[1][Ist1-1][0][Nmod-1] = Xxaty;
//    cout<<" Ist1: "<<Ist1<<" "<<Nmod<<" "<<Xzaty[0][Ist1-1][0][Nmod-1]<<" "<<Xzaty[1][Ist1-1][0][Nmod-1]<<endl;
//    cout<<Zzendc<<" "<<Xxendc<<" "<<Zzenda<<" "<<Xxenda<<endl;
//
      Zzendc = Xzstraw[2][Ist2-1][1][Nmod-1];
      Xxendc = Xzstraw[0][Ist2-1][1][Nmod-1];
      Zzenda = Xzstraw[5][Ist2-1][1][Nmod-1];
      Xxenda = Xzstraw[3][Ist2-1][1][Nmod-1];
      Zzaty  = Zzendc;
      Xxaty  = Xxendc + ((Xxenda-Xxendc)*(Ybest -(-4.55))/(9.1));
      Xzaty[0][Ist2-1][1][Nmod-1] = Zzaty;
      Xzaty[1][Ist2-1][1][Nmod-1] = Xxaty;
//    cout<<" Ist2: "<<Ist2<<" "<<Nmod<<" "<<Xzaty[0][Ist2-1][0][Nmod-1]<<" "<<Xzaty[1][Ist2-1][0][Nmod-1]<<endl;
//    cout<<Zzendc<<" "<<Xxendc<<" "<<Zzenda<<" "<<Xxenda<<endl;
//
      Zzendc = Xzstraw[2][Ist3-1][2][Nmod-1];
      Xxendc = Xzstraw[0][Ist3-1][2][Nmod-1];
      Zzenda = Xzstraw[5][Ist3-1][2][Nmod-1];
      Xxenda = Xzstraw[3][Ist3-1][2][Nmod-1];
      Zzaty  = Zzendc;
      Xxaty  = Xxendc + ((Xxenda-Xxendc)*(Ybest -(-4.55))/(9.1));
      Xzaty[0][Ist3-1][2][Nmod-1] = Zzaty;
      Xzaty[1][Ist3-1][2][Nmod-1] = Xxaty;
//    cout<<" Ist3: "<<Ist3<<" "<<Nmod<<" "<<Xzaty[0][Ist3-1][0][Nmod-1]<<" "<<Xzaty[1][Ist3-1][0][Nmod-1]<<endl;
//    cout<<Zzendc<<" "<<Xxendc<<" "<<Zzenda<<" "<<Xxenda<<endl;
//
      Zzendc = Xzstraw[2][Ist4-1][3][Nmod-1];
      Xxendc = Xzstraw[0][Ist4-1][3][Nmod-1];
      Zzenda = Xzstraw[5][Ist4-1][3][Nmod-1];
      Xxenda = Xzstraw[3][Ist4-1][3][Nmod-1];
      Zzaty  = Zzendc;
      Xxaty  = Xxendc + ((Xxenda-Xxendc)*(Ybest -(-4.55))/(9.1));
      Xzaty[0][Ist4-1][3][Nmod-1] = Zzaty;
      Xzaty[1][Ist4-1][3][Nmod-1] = Xxaty;
//    cout<<" Ist4: "<<Ist4<<" "<<Nmod<<" "<<Xzaty[0][Ist4-1][0][Nmod-1]<<" "<<Xzaty[1][Ist4-1][0][Nmod-1]<<endl;
//    cout<<Zzendc<<" "<<Xxendc<<" "<<Zzenda<<" "<<Xxenda<<endl;
//
      int Ntofit = Nposte;
      for(int ifit=1; ifit<=Nposte; ifit++) {
         Xptsfit[ifit-1] = Xypost[1][ifit-1][Nmod-1];
         Zptsfit[ifit-1] = Xypost[0][ifit-1][Nmod-1];
      }
      float Grfit = 0.0;
      float Cfit  = 0.0;
      Linfit(Zptsfit,Xptsfit,Ntofit,Grfit,Cfit);
      StrHist.fGrfit[0]->Fill(Grfit);
      float zte1 = 0.0;
      float xte1 = 0.0;
      float zte2 = 0.0;
      float xte2 = 0.0;
      Strclp(Grfit,Cfit,Zpt1b,Xpt1b,zte1,xte1);
      Strclp(Grfit,Cfit,Zpt4b,Xpt4b,zte2,xte2);
//
//    cout<<" zte1 xte1 zte2 xte2: "<<zte1<<" "<<xte1<<" "<<zte2<<" "<<xte2<<endl;
      float xdiff  = xte2 - xte1;
      float zdiff  = zte2 - zte1;
      float Phizx  = Bkatan(xdiff,zdiff);
      float Phizxd = Phizx*180.0/Pibk;
      if(Phizxd > 180.0) {
                         Phizxd = Phizxd - 360.0;
      }
      StrHist.fPhizxd[0]->Fill(Phizxd);
// --------------------------------------------
// Turn into milliradians:
// --------------------------------------------
      float Phizxmr = Phizxd*1000.0*Pibk/180.0;
//
      Trkinf[36][Nmod-1] = zte1;
      Trkinf[37][Nmod-1] = xte1;
      StrHist.fZxte[0]->Fill(zte1,xte1);
      Trkinf[38][Nmod-1] = zte2;
      Trkinf[39][Nmod-1] = xte2;
      float Grrr = (xte2 - xte1) / (zte2 - zte1);
      float Ctrk = xte1 - (Grrr * zte1);
      Trkinf[40][Nmod-1] = Grrr;
      Trkinf[41][Nmod-1] = Ctrk;
// ------------------------------------------------------------------
// Finally store residuals to fitted track.
// ------------------------------------------------------------------
      float Zintt1 = 0.0;
      float Xintt1 = 0.0; 
      Strclp(Grrr,Ctrk,Zwire1,Xwire1,Zintt1,Xintt1);
      float Zptt1 = 0.0;
      float Xptt1 = 0.0;
      Strcht(zte1,xte1,zte2,xte2,Ist1,Layer1,Nmod,Zptt1,Xptt1);
      float Zintt2 = 0.0;
      float Xintt2 = 0.0;
      Strclp(Grrr,Ctrk,Zwire2,Xwire2,Zintt2,Xintt2);
      float Zptt2 = 0.0;
      float Xptt2 = 0.0;
      Strcht(zte1,xte1,zte2,xte2,Ist2,Layer2,Nmod,Zptt2,Xptt2);
      float Zintt3 = 0.0;
      float Xintt3 = 0.0;
      Strclp(Grrr,Ctrk,Zwire3,Xwire3,Zintt3,Xintt3);
      float Zptt3 = 0.0;
      float Xptt3 = 0.0;
      Strcht(zte1,xte1,zte2,xte2,Ist3,Layer3,Nmod,Zptt3,Xptt3);
      float Zintt4 = 0.0;
      float Xintt4 = 0.0;
      Strclp(Grrr,Ctrk,Zwire4,Xwire4,Zintt4,Xintt4);
      float Zptt4 = 0.0;
      float Xptt4 = 0.0;
      Strcht(zte1,xte1,zte2,xte2,Ist4,Layer4,Nmod,Zptt4,Xptt4);
// ------------------------------------------------------------------
// (Zpt1,Xpt1) etc are the space points on the drift circles that 
// have been fitted to make the track.
// (Zptt1,Xptt1) are the closest points on the drift circles to the
// fitted track...
// ------------------------------------------------------------------
      cout<< " Layer1: "<<Ist1<<" "<<Zpt1b<<" "<<Xpt1b<<" "<<Zptt1<<" "<<Xptt1<<endl;
      cout<< " Layer2: "<<Ist2<<" "<<Zpt2b<<" "<<Xpt2b<<" "<<Zptt2<<" "<<Xptt2<<endl;
      cout<< " Layer3: "<<Ist3<<" "<<Zpt3b<<" "<<Xpt3b<<" "<<Zptt3<<" "<<Xptt3<<endl;
      cout<< " Layer4: "<<Ist4<<" "<<Zpt4b<<" "<<Xpt4b<<" "<<Zptt4<<" "<<Xptt4<<endl;
//
      float Rest1 = sqrt((Zptt1 - Zintt1)*(Zptt1 - Zintt1) + (Xptt1 - Xintt1)*(Xptt1 - Xintt1));
      Rest1 = Rest1*10000.0;
      float rrr = 0.0;
      rrr = BkRndm();
      if(rrr > 0.5) {
                    Rest1 = -Rest1;
      }
      float Rest2 = sqrt((Zptt2 - Zintt2)*(Zptt2 - Zintt2) + (Xptt2 - Xintt2)*(Xptt2 - Xintt2));
      Rest2 = Rest2*10000.0;
      rrr = BkRndm();
      if(rrr > 0.5) {
                    Rest2 = -Rest2;
      }

      float Rest3 = sqrt((Zptt3 - Zintt3)*(Zptt3 - Zintt3) + (Xptt3 - Xintt3)*(Xptt3 - Xintt3));
      Rest3 = Rest3*10000.0;
      rrr = BkRndm();
      if(rrr > 0.5) {
                    Rest3 = -Rest3;
      }
      float Rest4 = sqrt((Zptt4 - Zintt4)*(Zptt4 - Zintt4) + (Xptt4 - Xintt4)*(Xptt4 - Xintt4));
      Rest4 = Rest4*10000.0;
      rrr = BkRndm();
      if(rrr > 0.5) {
                    Rest4 = -Rest4;
      }
      StrHist.fRest1[0]->Fill(Rest1);
      StrHist.fRest2[0]->Fill(Rest2);
      StrHist.fRest3[0]->Fill(Rest3);
      StrHist.fRest4[0]->Fill(Rest4);
      StrHist.fRest[0]->Fill(Rest1);
      StrHist.fRest[0]->Fill(Rest2);
      StrHist.fRest[0]->Fill(Rest3);
      StrHist.fRest[0]->Fill(Rest4);
//
      Trkinf[42][Nmod-1] = Rest1;
      Trkinf[43][Nmod-1] = Rest2;
      Trkinf[44][Nmod-1] = Rest3;
      Trkinf[45][Nmod-1] = Rest4;
      Trkinf[46][Nmod-1] = Zptt1;
      Trkinf[47][Nmod-1] = Xptt1;
      Trkinf[48][Nmod-1] = Zptt2;
      Trkinf[49][Nmod-1] = Xptt2;
      Trkinf[50][Nmod-1] = Zptt3;
      Trkinf[51][Nmod-1] = Xptt3;
      Trkinf[52][Nmod-1] = Zptt4;
      Trkinf[53][Nmod-1] = Xptt4;
//
      Zxhit[0][Ist1-1][Layer1-1][Nmod-1] = Zptt1;
      Zxhit[1][Ist1-1][Layer1-1][Nmod-1] = Xptt1;
      Zxhit[0][Ist2-1][Layer2-1][Nmod-1] = Zptt2;
      Zxhit[1][Ist2-1][Layer2-1][Nmod-1] = Xptt2;
      Zxhit[0][Ist3-1][Layer3-1][Nmod-1] = Zptt3;
      Zxhit[1][Ist3-1][Layer3-1][Nmod-1] = Xptt3;
      Zxhit[0][Ist4-1][Layer4-1][Nmod-1] = Zptt4;
      Zxhit[1][Ist4-1][Layer4-1][Nmod-1] = Xptt4;
      cout<<" Made a 4 hit Track Element in Module "<<Nmod<<endl;
      StrHist.fTrkMod[0]->Fill(Nmod-0.5);
}
// -----------------------------------------------------------------
void Strte2(int Nmod,int Ist1,int Ist2,int Ist3,int Ist4) {
// ------------------------------------------------------------------
// Finds a TE in a Straw Module where hits exist in Layers 1
// and 2 and either layer 3 Or layer 4. 
// The TE is created if the position of the track in the layer 
// without a hit is compatible with having gone through the cell 
// wall.
// ------------------------------------------------------------------
      float Xptsfit[4];
      float Zptsfit[4];
//
      float Zw[3];
      float Xw[3];
      float Dr[3];
      float Resols[400];
      float Angloc[400];
      int   IsWall[400];
//
      float Pibk  = 3.1415926535;
      float Twopi = 2.0*Pibk;
      float Piby2 = Pibk/2.0;
      float Ysize = 0.1;
      float Ystep = 10.0/Ysize;
      float Yloop = 4.0 * Ystep;
      int   Kystep = Bknint(Ystep);
      int   Kyloop = Bknint(Yloop);
//
      for(int i=1; i<=Kyloop; i++) {
         Resols[i-1] = 99999999.0;
         Angloc[i-1] = 99999999.0;
         IsWall[i-1] = 0;
      }
// ------------------------------------------------------------------
// Declare some variables for later use.
// ------------------------------------------------------------------
      float Xbot = 0.0;
      float Ybot = 0.0;
      float Zbot = 0.0;
      float Xtop = 0.0;
      float Ytop = 0.0;
      float Ztop = 0.0;
      float xx1  = 0.0;
      float xx2  = 0.0;
      float xx3  = 0.0;
      float xx4  = 0.0;
      float zz1  = 0.0;
      float zz2  = 0.0;
      float zz3  = 0.0;
      float zz4  = 0.0;
      float yy1  = 0.0;
      float yy2  = 0.0;
      float yy3  = 0.0;
      float yy4  = 0.0;
      float xu1  = 0.0;
      float xu2  = 0.0;
      float xu3  = 0.0;
      float xu4  = 0.0;
      float zu1  = 0.0;
      float zu2  = 0.0;
      float zu3  = 0.0;
      float zu4  = 0.0;
      float yu1  = 0.0;
      float yu2  = 0.0;
      float yu3  = 0.0;
      float yu4  = 0.0;
//
      float Gradu = 0.0;
      float Cintu = 0.0;
      float Gradv = 0.0;
      float Cintv = 0.0;
      float xv1  = 0.0;
      float xv2  = 0.0;
      float xv3  = 0.0;
      float xv4  = 0.0;
      float zv1  = 0.0;
      float zv2  = 0.0;
      float zv3  = 0.0;
      float zv4  = 0.0;
      float yv1  = 0.0;
      float yv2  = 0.0;
      float yv3  = 0.0;
      float yv4  = 0.0;
//
      int Layer1 = 1;
      int Layer2 = 2;
      int Layer3 = 3;
      int Layer4 = 4;
      float Angexp = 0.0;
// ------------------------------------------------------------------
// Do track finding at various positions vertically and minimise the
// residuals:
// ------------------------------------------------------------------
      for(int Iybin=1; Iybin<=91; Iybin++) {
         float Ytest = float(Iybin)/10.0;
         Ytest = Ytest - 4.55;
         for(int ilay=1; ilay<=4; ilay++) {
            for(int istr=1; istr<=32; istr++) {
               float Zendc = Xzstraw[2][istr-1][ilay-1][Nmod-1];
               float Xendc = Xzstraw[0][istr-1][ilay-1][Nmod-1];
               float Zenda = Xzstraw[5][istr-1][ilay-1][Nmod-1];
               float Xenda = Xzstraw[3][istr-1][ilay-1][Nmod-1];
// ------------------------------------------------------------
// Old calculation in the xyz frame:
// ------------------------------------------------------------
               float Zaty  = Zendc;
               float Xaty  = Xendc + ((Xenda-Xendc)*(Ytest -(-4.55))/(9.1));
// ------------------------------------
// New calculation in the u or v frame:
// ------------------------------------   
               if(ilay <= 2) {
                 Xbot = Xzstrau[0][istr-1][ilay-1][Nmod-1];
                 Ybot = Xzstrau[1][istr-1][ilay-1][Nmod-1];
                 Zbot = Xzstrau[2][istr-1][ilay-1][Nmod-1];
                 Xtop = Xzstrau[3][istr-1][ilay-1][Nmod-1];
                 Ytop = Xzstrau[4][istr-1][ilay-1][Nmod-1];
                 Ztop = Xzstrau[5][istr-1][ilay-1][Nmod-1];
               }
               if(ilay > 2) {
                 Xbot = Xzstrav[0][istr-1][ilay-1][Nmod-1];
                 Ybot = Xzstrav[1][istr-1][ilay-1][Nmod-1];
                 Zbot = Xzstrav[2][istr-1][ilay-1][Nmod-1];
                 Xtop = Xzstrav[3][istr-1][ilay-1][Nmod-1];
                 Ytop = Xzstrav[4][istr-1][ilay-1][Nmod-1];
                 Ztop = Xzstrav[5][istr-1][ilay-1][Nmod-1];
               }
               Xzaty[0][istr-1][ilay-1][Nmod-1] = Ztop;
               Xzaty[1][istr-1][ilay-1][Nmod-1] = Xtop;
            }
         }
// ------------------------------------------------------------------
// Have found the hit straws in layers 1 and 2.
// Find the 4 possible tangents to these two drift circles.
// Extrapolate to layers 3 and 4 ...
// ------------------------------------------------------------------
         Zw[0] = Xzaty[0][Ist1-1][0][Nmod-1];
         Xw[0] = Xzaty[1][Ist1-1][0][Nmod-1];
         Dr[0] = Hits[Ist1-1][0][Nmod-1];
         Ybot  = Xzstrau[1][Ist1-1][0][Nmod-1];
         Ytop  = Xzstrau[4][Ist1-1][0][Nmod-1];
         yu1   = Ybot + ((Ytop-Ybot)*(Ytest -(-4.55))/(9.1));
//
         Zw[1] = Xzaty[0][Ist2-1][1][Nmod-1];
         Xw[1] = Xzaty[1][Ist2-1][1][Nmod-1];
         Dr[1] = Hits[Ist2-1][1][Nmod-1];
         Ybot  = Xzstrau[1][Ist2-1][1][Nmod-1];
         Ytop  = Xzstrau[4][Ist2-1][1][Nmod-1];
         yu2   = Ybot + ((Ytop-Ybot)*(Ytest -(-4.55))/(9.1));
         if(Ist3 > 0) {
           Zw[2] = Xzaty[0][Ist3-1][2][Nmod-1];
           Xw[2] = Xzaty[1][Ist3-1][2][Nmod-1];
           Dr[2] = Hits[Ist3-1][2][Nmod-1];
           Ybot  = Xzstrav[1][Ist3-1][2][Nmod-1];
           Ytop  = Xzstrav[4][Ist3-1][2][Nmod-1];
           yv3   = Ybot + ((Ytop-Ybot)*(Ytest -(-4.55))/(9.1));
         }
         if(Ist4 > 0) {
           Zw[2] = Xzaty[0][Ist4-1][3][Nmod-1];
           Xw[2] = Xzaty[1][Ist4-1][3][Nmod-1];
           Dr[2] = Hits[Ist4-1][3][Nmod-1];
           Ybot  = Xzstrav[1][Ist4-1][3][Nmod-1];
           Ytop  = Xzstrav[4][Ist4-1][3][Nmod-1];
           yv3   = Ybot + ((Ytop-Ybot)*(Ytest -(-4.55))/(9.1));
         }
//
         for(int itan=1; itan<=4; itan++) {
// ------------------------------------------------------------------
// Routine Strelm returns the Z and X points of the common
// tangents to the drift circles for each of the 4 cases in turn.
// ------------------------------------------------------------------
//   Itan = 1 MEANS OUTER EDGE TANGENT , RHS OF DRIFT CIRCLES
//       = 2        ''    ''     ''     LHS  ''  ''     ''
//       = 3       CROSSED TANGENT , BOTTOM RHS TO TOP LHS
//       = 4           ''     ''       ''   LHS        RHS
// ------------------------------------------------------------------
            int iflag = 0;
            zu1 = 0.0;
            xu1 = 0.0;
            zu2 = 0.0;
            xu2 = 0.0;
            Strelm(Nmod,Ist1,Layer1,Ist2,Layer2,itan,zu1,xu1,zu2,xu2,iflag);
            int Irloc = (Iybin-1)*4 + itan;
            int Icon = 1;
            if(fabs(zu2-zu1) > 1.0e-6) {
                    Gradu = (xu2 - xu1) / (zu2 - zu1);
                    Bku2xy(xu1,yu1,xx1,yy1);
                    Bku2xy(xu2,yu2,xx2,yy2);
                    zz2 = zu2;
                    zz1 = zu1;
                    Angloc[Irloc-1] = Bkatan((xx2-xx1),(zz2-zz1));
                    if(Angloc[Irloc-1] < 0.0) {
                      Angloc[Irloc-1] = Angloc[Irloc-1] + Twopi;
                    }
                    Angloc[Irloc-1] = Angloc[Irloc-1]*180.0/Pibk;
                    Icon = 1;
            }
            if(fabs(zu2-zu1) <= 1.0e-6) {
                    Resols[Irloc-1] = 99999.9;
                    Icon = 0;
            }
            if(Icon == 1) {
              Cintu = xu1 - (Gradu * zu1);
// -------------------------------------------------------------------
// Sanity check: Residuals in layers used to define the line should
// be zero.
// -------------------------------------------------------------------
              float Zwire1u = Zw[0];
              float Xwire1u = Xw[0];
              float Zint1u  = 0.0;
              float Xint1u  = 0.0;   
              Strclp(Gradu,Cintu,Zwire1u,Xwire1u,Zint1u,Xint1u);
              float Zwire2u = Zw[1];
              float Xwire2u = Xw[1];
              float Zint2u  = 0.0;
              float Xint2u  = 0.0;
              Strclp(Gradu,Cintu,Zwire2u,Xwire2u,Zint2u,Xint2u);
              float Zpt1u   = 0.0;
              float Xpt1u   = 0.0;
              Strcht(zu1,xu1,zu2,xu2,Ist1,Layer1,Nmod,Zpt1u,Xpt1u);
              float Res1u = sqrt((Zpt1u - Zint1u)*(Zpt1u - Zint1u) + (Xpt1u - Xint1u)*(Xpt1u - Xint1u));
              Res1u = Res1u*10000.0;
              float Zpt2u = 0.0;
              float Xpt2u = 0.0;
              Strcht(zu1,xu1,zu2,xu2,Ist2,Layer2,Nmod,Zpt2u,Xpt2u);
              float Res2u = sqrt((Zpt2u - Zint2u)*(Zpt2u - Zint2u) + (Xpt2u - Xint2u)*(Xpt2u - Xint2u));
              Res2u = Res2u*10000.0;
// -------------------------------------------------------------------
// Find closest point to this line on the hit cell in layer 3 or 4 
// We need to swap the extrapolated tangent into the v frame from the
// u frame:
// -------------------------------------------------------------------
              zv1 = zu1;
              Bku2xy(xu1,yu1,xx1,yy1);
              Bkxy2v(xx1,yy1,xv1,yv1);
//
              zv2 = zu2;
              Bku2xy(xu2,yu2,xx2,yy2);
              Bkxy2v(xx2,yy2,xv2,yv2);
//
              Gradv =(xv2 - xv1)/(zv2 - zv1);
              Cintv = xv1 - (Gradv * zv1); 
              float Zwire3v = Zw[2];
              float Xwire3v = Xw[2];
              float Zint3v  = 0.0;
              float Xint3v  = 0.0;
              Strclp(Gradv,Cintv,Zwire3v,Xwire3v,Zint3v,Xint3v);
              float Dist3v = sqrt((Zwire3v-Zint3v)*(Zwire3v-Zint3v) + (Xwire3v-Xint3v)*(Xwire3v-Xint3v));
              int Icon2 = 1;
              if(Dist3v > 1.2*Size[2]) {
                                  Resols[Irloc-1] = 99999.9;
                                  Icon2 = 0;
              }
// -------------------------------------------------------------------
// Extrapolated line passes through the hit cell in layer 3 or 4:
// Calculate the residuals:
// -------------------------------------------------------------------
              if(Icon2 == 1) {
                float Zpt3v = 0.0;
                float Xpt3v = 0.0;
                if(Ist3 > 0) {
                  Strcht(zv1,xv1,zv2,xv2,Ist3,Layer3,Nmod,Zpt3v,Xpt3v);
                }
                if(Ist4 > 0) {
                  Strcht(zv1,xv1,zv2,xv2,Ist4,Layer4,Nmod,Zpt3v,Xpt3v);
                }
                float Res3v = sqrt((Zpt3v - Zint3v)*(Zpt3v - Zint3v) + (Xpt3v - Xint3v)*(Xpt3v - Xint3v));
                Res3v = Res3v*10000.0;
                Resols[Irloc-1] = Res3v;
// -------------------------------------------------------------------
// We now need to decide if the layer that didnt have a hit is
// compatible with the track having gone through a cell wall.
// -------------------------------------------------------------------
                IsWall[Irloc-1] = 1;
                int Istart = 0;
                int Iend   = 0;
                if(Ist3 > 0) {
// ------------------
// No hit in Layer 4:
// ------------------
                  Istart = Ist3 - 3;  
                  if(Istart < 1) {
                               Istart = 1;
                  }
                  Iend = Ist3 + 3;
                  if(Iend > 32) {
                              Iend = 32;
                  }
                  for(int Is=Istart; Is<=Iend; Is++) {
                     float Zwirev = Xzaty[0][Is-1][3][Nmod-1];
                     float Xwirev = Xzaty[1][Is-1][3][Nmod-1];
                     float Zintv = 0.0;
                     float Xintv = 0.0;
                     Strclp(Gradv,Cintv,Zwirev,Xwirev,Zintv,Xintv);
                     float Distv = sqrt((Zwirev-Zintv)*(Zwirev-Zintv) + (Xwirev-Xintv)*(Xwirev-Xintv));
                     if(Distv < 0.25) {
                       IsWall[Irloc-1] = 0;
                     }
                  }
                }
//
                if(Ist4 > 0) {
// ------------------
// No hit in Layer 3:
// ------------------     
                  Istart = Ist4 - 3;
                  if(Istart < 1) {
                               Istart = 1;
                  }
                  Iend = Ist4 + 3;
                  if(Iend > 32) {
                              Iend = 32;
                  }
                  for(int Is=Istart; Is<=Iend; Is++) {
                     float Zwirevv = Xzaty[0][Is-1][2][Nmod-1];
                     float Xwirevv = Xzaty[1][Is-1][2][Nmod-1];
                     float Zintvv  = 0.0;
                     float Xintvv  = 0.0;
                     Strclp(Gradv,Cintv,Zwirevv,Xwirevv,Zintvv,Xintvv);
                     float Distvv = sqrt((Zwirevv-Zintvv)*(Zwirevv-Zintvv) + (Xwirevv-Xintvv)*(Xwirevv-Xintvv));
                     if(Distvv < 0.25) {
                       IsWall[Irloc-1] = 0;
                     }
                  }
                }
              }
            }
         }
      }
// ------------------------------------------------------------------
// We now have the residual information for extrapolation of each
// of the 4 tangents at each test Y position ... make a list of the
// acceptable possibilities.
// ------------------------------------------------------------------
      int Locbst   = 0;
      float Resbst = 99999.9;
// ------------------------------------------------------------------
// Keep solution with smallest Residual that is most 
// consistent with track going straight through.
// ------------------------------------------------------------------
      for(int isol=1; isol<=400; isol++) {
         float Resid = Resols[isol-1];
         float Angte = Angloc[isol-1];
         if(IsWall[isol-1] == 1) {
           float Angdif = fabs(Angte - Angexp);
           if(Angdif > 180.0) {
             Angdif = 360.0 - Angdif;
           }
           if(Angdif < 5.0) {
             if(Resid < Resbst) {
                        Resbst = Resid;
                        Locbst  = isol;
             }
           }
         } 
      }
//
      if(Resbst > 1000.0) {
                          return;
      }
// ------------------------------------------------------------------
// Refit the chosen solution:
// ------------------------------------------------------------------
      int Iybest  = int((Locbst - 1)/4) + 1;
      int Itnbst  = Locbst - ((Iybest-1)*4);
      float Ybest  = float(Iybest)/10.0;
      Ybest = Ybest - 4.55;
//
      for(int istr=1; istr<=32; istr++) {
         for(int ilay=1; ilay<=4; ilay++) {
            float Zendcb = Xzstraw[2][istr-1][ilay-1][Nmod-1];
            float Xendcb = Xzstraw[0][istr-1][ilay-1][Nmod-1];
            float Zendab = Xzstraw[5][istr-1][ilay-1][Nmod-1];
            float Xendab = Xzstraw[3][istr-1][ilay-1][Nmod-1];
            float Zatyb  = Zendcb;
            float Xatyb  = Xendcb + ((Xendab-Xendcb)*(Ybest -(-4.55))/(9.1));
            if(ilay <= 2) {
              Xbot = Xzstrau[0][istr-1][ilay-1][Nmod-1];
              Ybot = Xzstrau[1][istr-1][ilay-1][Nmod-1];
              Zbot = Xzstrau[2][istr-1][ilay-1][Nmod-1];
              Xtop = Xzstrau[3][istr-1][ilay-1][Nmod-1];
              Ytop = Xzstrau[4][istr-1][ilay-1][Nmod-1];
              Ztop = Xzstrau[5][istr-1][ilay-1][Nmod-1];
            }
            if(ilay > 2) {
              Xbot = Xzstrav[0][istr-1][ilay-1][Nmod-1];
              Ybot = Xzstrav[1][istr-1][ilay-1][Nmod-1];
              Zbot = Xzstrav[2][istr-1][ilay-1][Nmod-1];
              Xtop = Xzstrav[3][istr-1][ilay-1][Nmod-1];
              Ytop = Xzstrav[4][istr-1][ilay-1][Nmod-1];
              Ztop = Xzstrav[5][istr-1][ilay-1][Nmod-1];
            }
            Xzaty[0][istr-1][ilay-1][Nmod-1] = Ztop;
            Xzaty[1][istr-1][ilay-1][Nmod-1] = Xtop;
         }
      }
//
      Zw[0] = Xzaty[0][Ist1-1][0][Nmod-1];
      Xw[0] = Xzaty[1][Ist1-1][0][Nmod-1];
      Dr[0] = Hits[Ist1-1][0][Nmod-1];
      Ybot  = Xzstrau[1][Ist1-1][0][Nmod-1];
      Ytop  = Xzstrau[4][Ist1-1][0][Nmod-1];
      yu1   = Ybot + ((Ytop-Ybot)*(Ybest -(-4.55))/(9.1));
//
      Zw[1] = Xzaty[0][Ist2-1][1][Nmod-1];
      Xw[1] = Xzaty[1][Ist2-1][1][Nmod-1];
      Dr[1] = Hits[Ist2-1][1][Nmod-1];
      Ybot  = Xzstrau[1][Ist2-1][1][Nmod-1];
      Ytop  = Xzstrau[4][Ist2-1][1][Nmod-1];
      yu2   = Ybot + ((Ytop-Ybot)*(Ybest -(-4.55))/(9.1));
//
      if(Ist3 > 0) {
        Zw[2] = Xzaty[0][Ist3-1][2][Nmod-1];
        Xw[2] = Xzaty[1][Ist3-1][2][Nmod-1];
        Dr[2] = Hits[Ist3-1][2][Nmod-1];
        Ybot  = Xzstrav[1][Ist3-1][2][Nmod-1];
        Ytop  = Xzstrav[4][Ist3-1][2][Nmod-1];
        yv3   = Ybot + ((Ytop-Ybot)*(Ybest -(-4.55))/(9.1));
      }
      if(Ist4 > 0) {
        Zw[2] = Xzaty[0][Ist4-1][3][Nmod-1];
        Xw[2] = Xzaty[1][Ist4-1][3][Nmod-1];
        Dr[2] = Hits[Ist4-1][3][Nmod-1];
        Ybot  = Xzstrav[1][Ist4-1][3][Nmod-1];
        Ytop  = Xzstrav[4][Ist4-1][3][Nmod-1];
        yv3   = Ybot + ((Ytop-Ybot)*(Ybest -(-4.55))/(9.1));
      }
//
      int iflagb = 0;
      Strelm(Nmod,Ist1,Layer1,Ist2,Layer2,Itnbst,zu1,xu1,zu2,xu2,iflagb);
      if(fabs(zu2-zu1) > 1.0e-6) {
                                 Gradu = (xu2 - xu1)/(zu2 - zu1);
      }
      if(fabs(zu2-zu1) <= 1.0e-6) {
                                  return;
      }
      Cintu = xu1 - (Gradu * zu1);
      float Zwire1ub = Zw[0];
      float Xwire1ub = Xw[0];
      float Zint1ub  = 0.0;
      float Xint1ub  = 0.0;
      Strclp(Gradu,Cintu,Zwire1ub,Xwire1ub,Zint1ub,Xint1ub);
      float Dist1ub = sqrt((Zwire1ub-Zint1ub)*(Zwire1ub-Zint1ub) + (Xwire1ub-Xint1ub)*(Xwire1ub-Xint1ub));
      if(Dist1ub > 1.2*Size[0]) {
                                return;
      }
      float Zpt1ub = 0.0;
      float Xpt1ub = 0.0;
      Strcht(zu1,xu1,zu2,xu2,Ist1,Layer1,Nmod,Zpt1ub,Xpt1ub);
//
      float Zwire2ub = Zw[1];
      float Xwire2ub = Xw[1];
      float Zint2ub  = 0.0;
      float Xint2ub  = 0.0;
      Strclp(Gradu,Cintu,Zwire2ub,Xwire2ub,Zint2ub,Xint2ub);
      float Dist2ub = sqrt((Zwire2ub-Zint2ub)*(Zwire2ub-Zint2ub) + (Xwire2ub-Xint2ub)*(Xwire2ub-Xint2ub));
      if(Dist2ub > 1.2*Size[1]) {
                                return;
      }
      float Zpt2ub = 0.0;
      float Xpt2ub = 0.0;
      Strcht(zu1,xu1,zu2,xu2,Ist2,Layer2,Nmod,Zpt2ub,Xpt2ub);
//
      zv1 = zu1;
      Bku2xy(xu1,yu1,xx1,yy1);
      Bkxy2v(xx1,yy1,xv1,yv1);
//
      zv2 = zu2;
      Bku2xy(xu2,yu2,xx2,yy2);
      Bkxy2v(xx2,yy2,xv2,yv2);
//
      if(fabs(zv2-zv1) > 1.0e-6) {
                                 Gradv =(xv2 - xv1)/(zv2 - zv1);
      }
      if(fabs(zv2-zv1) <= 1.0e-6) {
                                  return;
      }
      Cintv = xv1 - (Gradv * zv1);
//
      float Zwire3vb = Zw[2];
      float Xwire3vb = Xw[2];
      float Zint3vb = 0.0;
      float Xint3vb = 0.0;
      Strclp(Gradv,Cintv,Zwire3vb,Xwire3vb,Zint3vb,Xint3vb);
      float Dist3vb = sqrt((Zwire3vb-Zint3vb)*(Zwire3vb-Zint3vb) + (Xwire3vb-Xint3vb)*(Xwire3vb-Xint3vb));
      if(Dist3vb > 1.2*Size[2]) {
                                return;
      }
      float Zpt3vb = 0.0;
      float Xpt3vb = 0.0;
      if(Ist3 > 0) {
        Strcht(zv1,xv1,zv2,xv2,Ist3,Layer3,Nmod,Zpt3vb,Xpt3vb);
      }
      if(Ist4 > 0) {
        Strcht(zv1,xv1,zv2,xv2,Ist4,Layer4,Nmod,Zpt3vb,Xpt3vb);
      }
      float Res1ub = sqrt((Zpt1ub - Zint1ub)*(Zpt1ub - Zint1ub) + (Xpt1ub - Xint1ub)*(Xpt1ub - Xint1ub));
      Res1ub = Res1ub*10000.0;
      float Res2ub = sqrt((Zpt2ub - Zint2ub)*(Zpt2ub - Zint2ub) + (Xpt2ub - Xint2ub)*(Xpt2ub - Xint2ub));
      Res2ub = Res2ub*10000.0;
      float Res3vb = sqrt((Zpt3vb - Zint3vb)*(Zpt3vb - Zint3vb) + (Xpt3vb - Xint3vb)*(Xpt3vb - Xint3vb));
      Res3vb = Res3vb*10000.0;
// ---------------------------------------------------------------------
// Convert all points back into xy frame for fitting:
// ---------------------------------------------------------------------
      float Zpt1b = Zpt1ub;
      float Xpt1b = 0.0;
      Bku2xy(Xpt1ub,yu1,Xpt1b,yy1);
//
      float Zpt2b = Zpt2ub;
      float Xpt2b = 0.0;
      Bku2xy(Xpt2ub,yu2,Xpt2b,yy2);
//
      float Zpt3b = Zpt3vb;
      float Xpt3b = 0.0;
      Bkv2xy(Xpt3vb,yv3,Xpt3b,yy3);
//
      Strflr(Ist1,Layer1,Nmod,Zpt1b,Xpt1b);
      Strflr(Ist2,Layer2,Nmod,Zpt2b,Xpt2b);
      if(Ist3 > 0) {
        Strflr(Ist3,Layer3,Nmod,Zpt3b,Xpt3b);
      }
      if(Ist4 > 0) {
        Strflr(Ist4,Layer4,Nmod,Zpt3b,Xpt3b);
      }
      Yhits[Ist1-1][Layer1-1][Nmod-1] = Ybest;
      Yhits[Ist2-1][Layer2-1][Nmod-1] = Ybest;
      if(Ist3 > 0) {
        Yhits[Ist3-1][Layer3-1][Nmod-1] = Ybest;
      }
      if(Ist4 > 0) {
        Yhits[Ist4-1][Layer4-1][Nmod-1] = Ybest;
      }
      Mask[Ist1-1][Layer1-1][Nmod-1] = 1;
      Mask[Ist2-1][Layer2-1][Nmod-1] = 1;
      if(Ist3 > 0) {
        Mask[Ist3-1][Layer3-1][Nmod-1] = 1;
      }
      if(Ist4 > 0) {
        Mask[Ist4-1][Layer4-1][Nmod-1] = 1;
      }
      Trkin2[0][Nmod-1] = float(Ist1);
      Trkin2[1][Nmod-1] = float(Ist2);
      Trkin2[2][Nmod-1] = float(Ist3);
      Trkin2[3][Nmod-1] = float(Ist4);
//
      Trkin2[4][Nmod-1] = Hits[Ist1-1][Layer1-1][Nmod-1];
      Trkin2[5][Nmod-1] = Zpt1b;
      Trkin2[6][Nmod-1] = Xpt1b;
      Trkin2[7][Nmod-1] = Ybest;
      float Zwire1b = Zwire1ub;
      float Xwire1b = 0.0;
      Bku2xy(Xwire1ub,yu1,Xwire1b,yy1);
      Trkin2[8][Nmod-1]  = Zwire1b;
      Trkin2[9][Nmod-1]  = Xwire1b;
      float Zint1b = Zint1ub;
      float Xint1b = 0.0;
      Bku2xy(Xint1ub,yu1,Xint1b,yy1);
      Trkin2[10][Nmod-1] = Zint1b;
      Trkin2[11][Nmod-1] = Xint1b;
//
      Trkin2[12][Nmod-1] = Hits[Ist2-1][Layer2-1][Nmod-1];
      Trkin2[13][Nmod-1] = Zpt2b;
      Trkin2[14][Nmod-1] = Xpt2b;
      Trkin2[15][Nmod-1] = Ybest;
      float Zwire2b = Zwire2ub;
      float Xwire2b = 0.0;
      Bku2xy(Xwire2ub,yu2,Xwire2b,yy2);
      float Zint2b = Zint2ub;
      float Xint2b = 0.0;
      Bku2xy(Xint2ub,yu2,Xint2b,yy2);
      Trkin2[16][Nmod-1] = Zwire2b;
      Trkin2[17][Nmod-1] = Xwire2b;
      Trkin2[18][Nmod-1] = Zint2b;
      Trkin2[19][Nmod-1] = Xint2b;
      float Zwire3b = 0.0;
      float Xwire3b = 0.0;
      float Zint3b  = 0.0;
      float Xint3b  = 0.0;
      if(Ist3 > 0) {
        Trkin2[20][Nmod-1] = Hits[Ist3-1][Layer3-1][Nmod-1];
        Trkin2[21][Nmod-1] = Zpt3b;
        Trkin2[22][Nmod-1] = Xpt3b;
        Trkin2[23][Nmod-1] = Ybest;
        Zwire3b = Zwire3vb;
        Xwire3b = 0.0;
        Bkv2xy(Xwire3vb,yv3,Xwire3b,yy3);
        Zint3b = Zint3vb;
        Xint3b = 0.0;
        Bkv2xy(Xint3vb,yv3,Xint3b,yy3);
        Trkin2[24][Nmod-1] = Zwire3b;
        Trkin2[25][Nmod-1] = Xwire3b;
        Trkin2[26][Nmod-1] = Zint3b;
        Trkin2[27][Nmod-1] = Xint3b;
      }
      if(Ist4 > 0) {
        Trkin2[28][Nmod-1] = Hits[Ist4-1][Layer4-1][Nmod-1];
        Trkin2[29][Nmod-1] = Zpt3b;
        Trkin2[30][Nmod-1] = Xpt3b;
        Trkin2[31][Nmod-1] = Ybest;
        Zwire3b = Zwire3vb;
        Xwire3b = 0.0;
        Bkv2xy(Xwire3vb,yv3,Xwire3b,yy3);
        Zint3b = Zint3vb;
        Xint3b = 0.0;
        Bkv2xy(Xint3vb,yv3,Xint3b,yy3);
        Trkin2[32][Nmod-1] = Zwire3b;
        Trkin2[33][Nmod-1] = Xwire3b;
        Trkin2[34][Nmod-1] = Zint3b;
        Trkin2[35][Nmod-1] = Xint3b;
      }
// ------------------------------------------------------------------
// Finally make the fitted TE from the obtained space points.
// Must set up xzaty in xy frame for this:
// ------------------------------------------------------------------
      float Zzendc = Xzstraw[2][Ist1-1][0][Nmod-1];
      float Xxendc = Xzstraw[0][Ist1-1][0][Nmod-1];
      float Zzenda = Xzstraw[5][Ist1-1][0][Nmod-1];
      float Xxenda = Xzstraw[3][Ist1-1][0][Nmod-1];
      float Zzaty  = Zzendc;
      float Xxaty  = Xxendc + ((Xxenda-Xxendc)*(Ybest -(-4.55))/(9.1));
      Xzaty[0][Ist1-1][0][Nmod-1] = Zzaty;   
      Xzaty[1][Ist1-1][0][Nmod-1] = Xxaty; 
//
      Zzendc = Xzstraw[2][Ist2-1][1][Nmod-1];
      Xxendc = Xzstraw[0][Ist2-1][1][Nmod-1];
      Zzenda = Xzstraw[5][Ist2-1][1][Nmod-1];
      Xxenda = Xzstraw[3][Ist2-1][1][Nmod-1];
      Zzaty  = Zzendc;
      Xxaty  = Xxendc + ((Xxenda-Xxendc)*(Ybest -(-4.55))/(9.1));
      Xzaty[0][Ist2-1][1][Nmod-1] = Zzaty;   
      Xzaty[1][Ist2-1][1][Nmod-1] = Xxaty;  
      if(Ist3 > 0) {
        Zzendc = Xzstraw[2][Ist3-1][2][Nmod-1];
        Xxendc = Xzstraw[0][Ist3-1][2][Nmod-1];
        Zzenda = Xzstraw[5][Ist3-1][2][Nmod-1];
        Xxenda = Xzstraw[3][Ist3-1][2][Nmod-1];
        Zzaty  = Zzendc;
        Xxaty  = Xxendc + ((Xxenda-Xxendc)*(Ybest -(-4.55))/(9.1));
        Xzaty[0][Ist3-1][2][Nmod-1] = Zzaty;   
        Xzaty[1][Ist3-1][2][Nmod-1] = Xxaty;   
      }
      if(Ist4 > 0) {
        Zzendc = Xzstraw[2][Ist4-1][3][Nmod-1];
        Xxendc = Xzstraw[0][Ist4-1][3][Nmod-1];
        Zzenda = Xzstraw[5][Ist4-1][3][Nmod-1];
        Xxenda = Xzstraw[3][Ist4-1][3][Nmod-1];
        Zzaty  = Zzendc;
        Xxaty  = Xxendc + ((Xxenda-Xxendc)*(Ybest -(-4.55))/(9.1));
        Xzaty[0][Ist4-1][3][Nmod-1] = Zzaty;
        Xzaty[1][Ist4-1][3][Nmod-1] = Xxaty;
      }
//
      int Ntofit = Nposte;
      for(int ifit=1; ifit<=Nposte; ifit++) {
         Xptsfit[ifit-1] = Xypost[1][ifit-1][Nmod-1];
         Zptsfit[ifit-1] = Xypost[0][ifit-1][Nmod-1];
      }
      float Grfit = 0.0;
      float Cfit  = 0.0;
      Linfit(Zptsfit,Xptsfit,Ntofit,Grfit,Cfit);
      float zte1 = 0.0;
      float xte1 = 0.0;
      float zte2 = 0.0;
      float xte2 = 0.0;
      Strclp(Grfit,Cfit,Zpt1b,Xpt1b,zte1,xte1);
      Strclp(Grfit,Cfit,Zpt3b,Xpt3b,zte2,xte2);
//
      float xdiff  = xte2 - xte1;
      float zdiff  = zte2 - zte1;
      float Phizx  = Bkatan(xdiff,zdiff);
      float Phizxd = Phizx*180.0/Pibk;
      if(Phizxd > 180.0) {
                         Phizxd = Phizxd - 360.0;
      }
// --------------------------------------------
// Turn into milliradians:
// --------------------------------------------
      float Phizxmr = Phizxd*1000.0*Pibk/180.0;
//
      Trkin2[36][Nmod-1] = zte1;
      Trkin2[37][Nmod-1] = xte1;
      Trkin2[38][Nmod-1] = zte2;
      Trkin2[39][Nmod-1] = xte2;
      float Grrr = (xte2 - xte1) / (zte2 - zte1);
      float Ctrk = xte1 - (Grrr * zte1);
      Trkin2[40][Nmod-1] = Grrr;
      Trkin2[41][Nmod-1] = Ctrk;
// ------------------------------------------------------------------
// Finally store residuals to fitted track.
// ------------------------------------------------------------------
      float Zintt1 = 0.0;
      float Xintt1 = 0.0;
      Strclp(Grrr,Ctrk,Zwire1b,Xwire1b,Zintt1,Xintt1);
      float Zptt1  = 0.0;
      float Xptt1  = 0.0;
      Strcht(zte1,xte1,zte2,xte2,Ist1,Layer1,Nmod,Zptt1,Xptt1);
      float Zintt2 = 0.0;
      float Xintt2 = 0.0;
      Strclp(Grrr,Ctrk,Zwire2b,Xwire2b,Zintt2,Xintt2);
      float Zptt2  = 0.0;
      float Xptt2  = 0.0;
      Strcht(zte1,xte1,zte2,xte2,Ist2,Layer2,Nmod,Zptt2,Xptt2);
      float Zintt3 = 0.0;
      float Xintt3 = 0.0;
      Strclp(Grrr,Ctrk,Zwire3b,Xwire3b,Zintt3,Xintt3);
      float Zptt3  = 0.0;
      float Xptt3  = 0.0;
      if(Ist3 > 0) {
        Strcht(zte1,xte1,zte2,xte2,Ist3,Layer3,Nmod,Zptt3,Xptt3);
      }
      if(Ist4 > 0) {
        Strcht(zte1,xte1,zte2,xte2,Ist4,Layer4,Nmod,Zptt3,Xptt3);
      }
// ------------------------------------------------------------------
// (Zpt1,Xpt1) etc are the space points on the drift circles that 
// have been fitted to make the track.
// (Zptt1,Xptt1) are the closest points on the drift circles to the
// fitted track...
// ------------------------------------------------------------------
      float Rest1 = sqrt((Zptt1 - Zintt1)*(Zptt1 - Zintt1) + (Xptt1 - Xintt1)*(Xptt1 - Xintt1));
      Rest1 = Rest1*10000.0;
      float Rest2 = sqrt((Zptt2 - Zintt2)*(Zptt2 - Zintt2) + (Xptt2 - Xintt2)*(Xptt2 - Xintt2));
      Rest2 = Rest2*10000.0;
      float Rest3 = sqrt((Zptt3 - Zintt3)*(Zptt3 - Zintt3) + (Xptt3 - Xintt3)*(Xptt3 - Xintt3));
      Rest3 = Rest3*10000.0;
      Trkin2[42][Nmod-1] = Rest1;
      Trkin2[43][Nmod-1] = Rest2;
      if(Ist3 > 0) {
        Trkin2[44][Nmod-1] = Rest3;
      }
      if(Ist4 > 0) {
        Trkin2[45][Nmod-1] = Rest3;
      }
      Trkin2[46][Nmod-1] = Zptt1;
      Trkin2[47][Nmod-1] = Xptt1;
      Trkin2[48][Nmod-1] = Zptt2;
      Trkin2[49][Nmod-1] = Xptt2;
      if(Ist3 > 0) {
        Trkin2[50][Nmod-1] = Zptt3;
        Trkin2[51][Nmod-1] = Xptt3;
      }
      if(Ist4 > 0) {
        Trkin2[52][Nmod-1] = Zptt3;
        Trkin2[53][Nmod-1] = Xptt3;
      }
      Zxhit[0][Ist1-1][Layer1-1][Nmod-1] = Zptt1;
      Zxhit[1][Ist1-1][Layer1-1][Nmod-1] = Xptt1;
      Zxhit[0][Ist2-1][Layer2-1][Nmod-1] = Zptt2;
      Zxhit[1][Ist2-1][Layer2-1][Nmod-1] = Xptt2;
      if(Ist3 > 0) {
        Zxhit[0][Ist3-1][Layer3-1][Nmod-1] = Zptt3;
        Zxhit[1][Ist3-1][Layer3-1][Nmod-1] = Xptt3;
      }
      if(Ist4 > 0) {
        Zxhit[0][Ist4-1][Layer4-1][Nmod-1] = Zptt3;
        Zxhit[1][Ist4-1][Layer4-1][Nmod-1] = Xptt3;
      }
}
// --------------------------------------------------------------
void Strte3(int Nmod,int Ist1,int Ist2,int Ist3,int Ist4) {
// ------------------------------------------------------------------
// Finds a TE in a Straw Collection where hits exist in Layers 3
// and 4 and either layer 1 Or layer 2. 
// ------------------------------------------------------------------
      float Xptsfit[4];
      float Zptsfit[4];
//
      float Zw[3];
      float Xw[3];
      float Dr[3];
      float Resols[400];
      float Angloc[400];
      int   IsWall[400];
//
      float Pibk  = 3.1415926535;
      float Twopi = 2.0*Pibk;
      float Piby2 = Pibk/2.0;
      float Ysize = 0.1;
      float Ystep = 10.0/Ysize;
      float Yloop = 4.0 * Ystep;
      int Kystep = Bknint(Ystep);
      int Kyloop = Bknint(Yloop);
//
      for(int i=1; i<=Kyloop; i++) {
         Resols[i-1] = 9999999.0;
         Angloc[i-1] = 9999999.9;
         IsWall[i-1] = 0;
      }
// ------------------------------------------------------------------
// Declare some variables for later use.
// ------------------------------------------------------------------
      float Xbot = 0.0;
      float Ybot = 0.0;
      float Zbot = 0.0;
      float Xtop = 0.0;
      float Ytop = 0.0;
      float Ztop = 0.0;
      float xx1  = 0.0;
      float xx2  = 0.0;
      float xx3  = 0.0;
      float xx4  = 0.0;
      float zz1  = 0.0;
      float zz2  = 0.0;
      float zz3  = 0.0;
      float zz4  = 0.0;
      float yy1  = 0.0;
      float yy2  = 0.0;
      float yy3  = 0.0;
      float yy4  = 0.0;
      float xu1  = 0.0;
      float xu2  = 0.0;
      float xu3  = 0.0;
      float xu4  = 0.0;
      float zu1  = 0.0;
      float zu2  = 0.0;
      float zu3  = 0.0;
      float zu4  = 0.0;
      float yu1  = 0.0;
      float yu2  = 0.0;
      float yu3  = 0.0;
      float yu4  = 0.0;
//
      float Gradu = 0.0;
      float Cintu = 0.0;
      float Gradv = 0.0;
      float Cintv = 0.0;
      float xv1  = 0.0;
      float xv2  = 0.0;
      float xv3  = 0.0;
      float xv4  = 0.0;
      float zv1  = 0.0;
      float zv2  = 0.0;
      float zv3  = 0.0;
      float zv4  = 0.0;
      float yv1  = 0.0;
      float yv2  = 0.0;
      float yv3  = 0.0;
      float yv4  = 0.0;
//
      int Layer1 = 1;
      int Layer2 = 2;
      int Layer3 = 3;
      int Layer4 = 4;
      float Angexp = 0.0;
// ------------------------------------------------------------------
// Do track finding at various positions vertically and minimise the
// residuals:
// ------------------------------------------------------------------
      for(int Iybin=1; Iybin<=91; Iybin++) {
         float Ytest = float(Iybin)/10.0;
         Ytest = Ytest - 4.55;
         for(int ilay=1; ilay<=4; ilay++) {
            for(int istr=1; istr<=32; istr++) {
               float Zendc = Xzstraw[2][istr-1][ilay-1][Nmod-1];
               float Xendc = Xzstraw[0][istr-1][ilay-1][Nmod-1];
               float Zenda = Xzstraw[5][istr-1][ilay-1][Nmod-1];
               float Xenda = Xzstraw[3][istr-1][ilay-1][Nmod-1];
// ---------------------------------
// Old calculation in the xyz frame:
// ---------------------------------
               float Zaty  = Zendc;
               float Xaty  = Xendc + ((Xenda-Xendc)*(Ytest -(-4.55))/(9.1));
// ------------------------------------
// New calculation in the u or v frame:
// ------------------------------------   
               if(ilay <= 2) {
                 Xbot = Xzstrau[0][istr-1][ilay-1][Nmod-1];
                 Ybot = Xzstrau[1][istr-1][ilay-1][Nmod-1];
                 Zbot = Xzstrau[2][istr-1][ilay-1][Nmod-1];
                 Xtop = Xzstrau[3][istr-1][ilay-1][Nmod-1];
                 Ytop = Xzstrau[4][istr-1][ilay-1][Nmod-1];
                 Ztop = Xzstrau[5][istr-1][ilay-1][Nmod-1];
               }
               if(ilay > 2) {
                 Xbot = Xzstrav[0][istr-1][ilay-1][Nmod-1];
                 Ybot = Xzstrav[1][istr-1][ilay-1][Nmod-1];
                 Zbot = Xzstrav[2][istr-1][ilay-1][Nmod-1];
                 Xtop = Xzstrav[3][istr-1][ilay-1][Nmod-1];
                 Ytop = Xzstrav[4][istr-1][ilay-1][Nmod-1];
                 Ztop = Xzstrav[5][istr-1][ilay-1][Nmod-1];
               }
               Xzaty[0][istr-1][ilay-1][Nmod-1] = Ztop;
               Xzaty[1][istr-1][ilay-1][Nmod-1] = Xtop;
            }
         }
// ------------------------------------------------------------------
// Have found the hit straws in layers 3 and 4.
// Find the 4 possible tangents to these two drift circles.
// Extrapolate back to layers 2 and 1 ...
// ------------------------------------------------------------------
         Zw[0] = Xzaty[0][Ist3-1][2][Nmod-1];
         Xw[0] = Xzaty[1][Ist3-1][2][Nmod-1];
         Dr[0] = Hits[Ist3-1][2][Nmod-1];
         Ybot  = Xzstrav[1][Ist3-1][2][Nmod-1];
         Ytop  = Xzstrav[4][Ist3-1][2][Nmod-1];
         yv1   = Ybot + ((Ytop-Ybot)*(Ytest -(-4.55))/(9.1));
//
         Zw[1] = Xzaty[0][Ist4-1][3][Nmod-1];
         Xw[1] = Xzaty[1][Ist4-1][3][Nmod-1];
         Dr[1] = Hits[Ist4-1][3][Nmod-1];
         Ybot  = Xzstrav[1][Ist4-1][3][Nmod-1];
         Ytop  = Xzstrav[4][Ist4-1][3][Nmod-1];
         yv2   = Ybot + ((Ytop-Ybot)*(Ytest -(-4.55))/(9.1));
//
         if(Ist1 > 0) {
           Zw[2] = Xzaty[0][Ist1-1][0][Nmod-1];
           Xw[2] = Xzaty[1][Ist1-1][0][Nmod-1];
           Dr[2] = Hits[Ist1-1][0][Nmod-1];
           Ybot  = Xzstrau[1][Ist1-1][0][Nmod-1];
           Ytop  = Xzstrau[4][Ist1-1][0][Nmod-1];
           yu3   = Ybot + ((Ytop-Ybot)*(Ytest -(-4.55))/(9.1));
         }
//
         if(Ist2 > 0) {
           Zw[2] = Xzaty[0][Ist2-1][1][Nmod-1];
           Xw[2] = Xzaty[1][Ist2-1][1][Nmod-1];
           Dr[2] = Hits[Ist2-1][1][Nmod-1];
           Ybot  = Xzstrau[1][Ist2-1][1][Nmod-1];
           Ytop  = Xzstrau[4][Ist2-1][1][Nmod-1];
           yu3   = Ybot + ((Ytop-Ybot)*(Ytest -(-4.55))/(9.1));
         }
//
         for(int itan=1; itan<=4; itan++) {
// ------------------------------------------------------------------
// Routine Strelm returns the Z and X points of the common
// tangents to the drift circles for each of the 4 cases in turn.
// ------------------------------------------------------------------
//   Itan = 1 MEANS OUTER EDGE TANGENT , RHS OF DRIFT CIRCLES
//        = 2        ''    ''     ''     LHS  ''  ''     ''
//        = 3       CROSSED TANGENT , BOTTOM RHS TO TOP LHS
//        = 4           ''     ''       ''   LHS        RHS
// ------------------------------------------------------------------
            int iflag = 0;
            Strelm(Nmod,Ist3,Layer3,Ist4,Layer4,itan,zv1,xv1,zv2,xv2,iflag);
            int Irloc = (Iybin-1)*4 + itan;
            int Icon = 1;
            if(fabs(zv2-zv1) > 1.0e-6) {
              Gradv = (xv2 - xv1) / (zv2 - zv1);
              Bkv2xy(xv1,yv1,xx1,yy1);
              Bkv2xy(xv2,yv2,xx2,yy2);
              zz1 = zv1;
              zz2 = zv2;
              Angloc[Irloc-1] = Bkatan((xx2-xx1),(zz2-zz1));
              if(Angloc[Irloc-1] < 0.0) {
                Angloc[Irloc-1] = Angloc[Irloc-1] + Twopi;
              }
              Angloc[Irloc-1] = Angloc[Irloc-1]*180.0/Pibk;
              Icon = 1;
            } 
            if(fabs(zv2-zv1) <= 1.0e-6) {
              Resols[Irloc-1] = 99999.9;
              Icon = 0;
            }
            if(Icon == 1) {
              Cintv = xv1 - (Gradv * zv1);
// -------------------------------------------------------------------
// Sanity check: Residuals in layers used to define the line should
// be zero.
// -------------------------------------------------------------------
              float Zwire3v = Zw[0];
              float Xwire3v = Xw[0];
              float Zint3v = 0.0;
              float Xint3v = 0.0;
              Strclp(Gradv,Cintv,Zwire3v,Xwire3v,Zint3v,Xint3v);
              float Zwire4v = Zw[1];
              float Xwire4v = Xw[1];
              float Zint4v = 0.0;
              float Xint4v = 0.0;
              Strclp(Gradv,Cintv,Zwire4v,Xwire4v,Zint4v,Xint4v);
              float Zpt3v = 0.0;
              float Xpt3v = 0.0;
              Strcht(zv1,xv1,zv2,xv2,Ist3,Layer3,Nmod,Zpt3v,Xpt3v);
              float Res3v = sqrt((Zpt3v - Zint3v)*(Zpt3v - Zint3v) + (Xpt3v - Xint3v)*(Xpt3v - Xint3v));
              Res3v = Res3v*10000.0;
              float Zpt4v = 0.0;
              float Xpt4v = 0.0;
              Strcht(zv1,xv1,zv2,xv2,Ist4,Layer4,Nmod,Zpt4v,Xpt4v);
              float Res4v = sqrt((Zpt4v - Zint4v)*(Zpt4v - Zint4v) + (Xpt4v - Xint4v)*(Xpt4v - Xint4v));
              Res4v = Res4v*10000.0;
// -------------------------------------------------------------------
// Find closest point to this line on the hit cell in layer 2 or 1 
// We need to swap the extrapolated tangent into the u frame from the
// v frame:
// -------------------------------------------------------------------
              zu1 = zv1;
              Bkv2xy(xv1,yv1,xx1,yy1);
              Bkxy2u(xx1,yy1,xu1,yu1);
//
              zu2 = zv2;
              Bkv2xy(xv2,yv2,xx2,yy2);
              Bkxy2u(xx2,yy2,xu2,yu2);
//
              Gradu =(xu2 - xu1)/(zu2 - zu1);
              Cintu = xu1 - (Gradu * zu1); 
//
              float Zwire21u = Zw[2];
              float Xwire21u = Xw[2];
              float Zint21u  = 0.0;
              float Xint21u  = 0.0;
              Strclp(Gradu,Cintu,Zwire21u,Xwire21u,Zint21u,Xint21u);
              float Dist21u = sqrt((Zwire21u-Zint21u)*(Zwire21u-Zint21u) + (Xwire21u-Xint21u)*(Xwire21u-Xint21u));
              int Icon2 = 1;
              if(Dist21u > 1.2*Size[1]) {
                             Resols[Irloc-1] = 99999.9;
                             Icon2 = 0;
              }
// -------------------------------------------------------------------
// Extrapolated line passes through the hit cell in layer 2 or 1:
// Calculate the residuals:
// -------------------------------------------------------------------
              if(Icon2 == 1) {
                float Zpt21u = 0.0;
                float Xpt21u = 0.0;
                if(Ist2 > 0) {
                  Strcht(zu1,xu1,zu2,xu2,Ist2,Layer2,Nmod,Zpt21u,Xpt21u);
                }
                if(Ist1 > 0) {
                  Strcht(zu1,xu1,zu2,xu2,Ist1,Layer1,Nmod,Zpt21u,Xpt21u);
                }
                float Res21u = sqrt((Zpt21u - Zint21u)*(Zpt21u - Zint21u) + (Xpt21u - Xint21u)*(Xpt21u - Xint21u));
                Res21u = Res21u*10000.0;
                Resols[Irloc-1] = Res21u;
// -------------------------------------------------------------------
// We now need to decide if the layer that didnt have a hit is
// compatible with the track having gone through a cell wall.
// -------------------------------------------------------------------
                IsWall[Irloc-1] = 1;
                int Istart = 0;
                int Iend   = 0;
                if(Ist2 > 0) {
// ------------------
// No hit in Layer 1:
// ------------------
                  Istart = Ist2 - 3;  
                  if(Istart < 1) {
                           Istart = 1;
                  }
                  Iend = Ist2 + 3;
                  if(Iend > 32) {
                           Iend = 32;
                  }
                  for(int is=Istart; is<=Iend; is++) {
                     float Zwireu = Xzaty[0][is-1][0][Nmod-1];
                     float Xwireu = Xzaty[1][is-1][0][Nmod-1];
                     float Zintu = 0.0;
                     float Xintu = 0.0;
                     Strclp(Gradu,Cintu,Zwireu,Xwireu,Zintu,Xintu);
                     float Distu = sqrt((Zwireu-Zintu)*(Zwireu-Zintu) + (Xwireu-Xintu)*(Xwireu-Xintu));
                     if(Distu < 0.25) {
                       IsWall[Irloc-1] = 0;
                     }
                  }
                }
//
                if(Ist1 > 0) {
// ------------------
// No hit in Layer 2:
// ------------------     
                  Istart = Ist1 - 3;
                  if(Istart < 1) {
                                 Istart = 1;
                  }
                  Iend = Ist1 + 3;
                  if(Iend > 32) {
                                Iend = 32;
                  }
                  for(int is=Istart; is<=Iend; is++) {
                     float Zwireuu = Xzaty[0][is-1][1][Nmod-1];
                     float Xwireuu = Xzaty[1][is-1][1][Nmod-1];
                     float Zintuu = 0.0;
                     float Xintuu = 0.0;
                     Strclp(Gradu,Cintu,Zwireuu,Xwireuu,Zintuu,Xintuu);
                     float Distuu = sqrt((Zwireuu-Zintuu)*(Zwireuu-Zintuu) + (Xwireuu-Xintuu)*(Xwireuu-Xintuu));
                     if(Distuu < 0.25) {
                       IsWall[Irloc-1] = 0;
                     }
                  }
                }
              }
            }
         }
      }
// ------------------------------------------------------------------
// We now have the residual information for extrapolation of each
// of the 4 tangents at each test Y position ... choose the best.
// ------------------------------------------------------------------
      int Locbst = 0;
      float Resbst = 99999.9;
      for(int isol=1; isol<=400; isol++) {
         float Resid = Resols[isol-1];
         float Angte = Angloc[isol-1];
         if(IsWall[isol-1] == 1) {
           float Angdif = fabs(Angte - Angexp);
           if(Angdif > 180.0) {
             Angdif = 360.0 - Angdif;
           }
           if(Angdif < 5.0) {
             if(Resid < Resbst) {
                        Resbst = Resid;
                        Locbst  = isol;
             }
           } 
         }
      }   
//
      if(Resbst > 1000.0) {
                          return;
      }
// ------------------------------------------------------------------
// Best chosen .... now refit it:
// ------------------------------------------------------------------
      int Iybest  = int((Locbst-1)/4) + 1;
      int Itnbst  = Locbst - ((Iybest-1)*4);
      float Ybest = float(Iybest)/10.0;
      Ybest = Ybest - 4.55;
//
      for(int istr=1; istr<=32; istr++) {
         for(int ilay=1; ilay<=4; ilay++) {
            float Zendcb = Xzstraw[2][istr-1][ilay-1][Nmod-1];
            float Xendcb = Xzstraw[0][istr-1][ilay-1][Nmod-1];
            float Zendab = Xzstraw[5][istr-1][ilay-1][Nmod-1];
            float Xendab = Xzstraw[3][istr-1][ilay-1][Nmod-1];
            float Zatyb  = Zendcb;   
            float Xatyb  = Xendcb + ((Xendab-Xendcb)*(Ybest -(-4.55))/(9.1));
// ------------------------------------
// New calculation in the u or v frame:
// ------------------------------------   
            if(ilay <= 2) {
              Xbot = Xzstrau[0][istr-1][ilay-1][Nmod-1];
              Ybot = Xzstrau[1][istr-1][ilay-1][Nmod-1];
              Zbot = Xzstrau[2][istr-1][ilay-1][Nmod-1];
              Xtop = Xzstrau[3][istr-1][ilay-1][Nmod-1];
              Ytop = Xzstrau[4][istr-1][ilay-1][Nmod-1];
              Ztop = Xzstrau[5][istr-1][ilay-1][Nmod-1];
            }
            if(ilay > 2) {
              Xbot = Xzstrav[0][istr-1][ilay-1][Nmod-1];
              Ybot = Xzstrav[1][istr-1][ilay-1][Nmod-1];
              Zbot = Xzstrav[2][istr-1][ilay-1][Nmod-1];
              Xtop = Xzstrav[3][istr-1][ilay-1][Nmod-1];
              Ytop = Xzstrav[4][istr-1][ilay-1][Nmod-1];
              Ztop = Xzstrav[5][istr-1][ilay-1][Nmod-1];
            }
            Xzaty[0][istr-1][ilay-1][Nmod-1] = Ztop;
            Xzaty[1][istr-1][ilay-1][Nmod-1] = Xtop;
         }
      }
//
      Zw[0] = Xzaty[0][Ist3-1][2][Nmod-1];
      Xw[0] = Xzaty[1][Ist3-1][2][Nmod-1];
      Dr[0] = Hits[Ist3-1][2][Nmod-1];
      Ybot  = Xzstrav[1][Ist3-1][2][Nmod-1];
      Ytop  = Xzstrav[4][Ist3-1][2][Nmod-1];
      yv1   = Ybot + ((Ytop-Ybot)*(Ybest -(-4.55))/(9.1));
//
      Zw[1] = Xzaty[0][Ist4-1][3][Nmod-1];
      Xw[1] = Xzaty[1][Ist4-1][3][Nmod-1];
      Dr[1] = Hits[Ist4-1][3][Nmod-1];
      Ybot  = Xzstrav[1][Ist4-1][3][Nmod-1];
      Ytop  = Xzstrav[4][Ist4-1][3][Nmod-1];
      yv2   = Ybot + ((Ytop-Ybot)*(Ybest -(-4.55))/(9.1));
//
      if(Ist1 > 0) {
        Zw[2] = Xzaty[0][Ist1-1][0][Nmod-1];
        Xw[2] = Xzaty[1][Ist1-1][0][Nmod-1];
        Dr[2] = Hits[Ist1-1][0][Nmod-1];
        Ybot  = Xzstrau[1][Ist1-1][0][Nmod-1];
        Ytop  = Xzstrau[4][Ist1-1][0][Nmod-1];
        yu3   = Ybot + ((Ytop-Ybot)*(Ybest -(-4.55))/(9.1));
      }
      if(Ist2 > 0) {
        Zw[2] = Xzaty[0][Ist2-1][1][Nmod-1];
        Xw[2] = Xzaty[1][Ist2-1][1][Nmod-1];
        Dr[2] = Hits[Ist2-1][1][Nmod-1];
        Ybot  = Xzstrau[1][Ist2-1][1][Nmod-1];
        Ytop  = Xzstrau[4][Ist2-1][1][Nmod-1];
        yu3   = Ybot + ((Ytop-Ybot)*(Ybest -(-4.55))/(9.1));
      }
//
      int iflagb = 0;
      Strelm(Nmod,Ist3,Layer3,Ist4,Layer4,Itnbst,zv1,xv1,zv2,xv2,iflagb);
      if(fabs(zv2-zv1) > 1.0e-6) {
                                 Gradv = (xv2 - xv1) / (zv2 - zv1);
      }
      if(fabs(zv2-zv1) <= 1.0e-6) {
                                  return;
      }
      Cintv = xv1 - (Gradv * zv1);
//
      float Zwire3vb = Zw[0];
      float Xwire3vb = Xw[0];
      float Zint3vb  = 0.0;
      float Xint3vb  = 0.0;
      Strclp(Gradv,Cintv,Zwire3vb,Xwire3vb,Zint3vb,Xint3vb);
      float Dist3vb = sqrt((Zwire3vb-Zint3vb)*(Zwire3vb-Zint3vb) + (Xwire3vb-Xint3vb)*(Xwire3vb-Xint3vb));
      if(Dist3vb > 1.2*Size[0]) {
                                return;
      }
      float Zpt3vb = 0.0;
      float Xpt3vb = 0.0;
      Strcht(zv1,xv1,zv2,xv2,Ist3,Layer3,Nmod,Zpt3vb,Xpt3vb);
      float Zwire4vb = Zw[1];
      float Xwire4vb = Xw[1];
      float Zint4vb  = 0.0;
      float Xint4vb  = 0.0;
      Strclp(Gradv,Cintv,Zwire4vb,Xwire4vb,Zint4vb,Xint4vb);
      float Dist4vb = sqrt((Zwire4vb-Zint4vb)*(Zwire4vb-Zint4vb) + (Xwire4vb-Xint4vb)*(Xwire4vb-Xint4vb));
      if(Dist4vb > 1.2*Size[1]) {
                                return;
      }
      float Zpt4vb = 0.0;
      float Xpt4vb = 0.0;
      Strcht(zv1,xv1,zv2,xv2,Ist4,Layer4,Nmod,Zpt4vb,Xpt4vb);
//
      zu1 = zv1;
      Bkv2xy(xv1,yv1,xx1,yy1);
      Bkxy2u(xx1,yy1,xu1,yu1);
//
      zu2 = zv2;
      Bkv2xy(xv2,yv2,xx2,yy2);
      Bkxy2u(xx2,yy2,xu2,yu2);
//
      if(fabs(zu2-zu1) > 1.0e-6) {
                                 Gradu =(xu2 - xu1)/(zu2 - zu1);
      }
      if(fabs(zu2-zu1) <= 1.0e-6) {
                                  return;
      }
      Cintu = xu1 - (Gradu * zu1);
//
      float Zwire12ub = Zw[2];
      float Xwire12ub = Xw[2];
      float Zint12ub  = 0.0;
      float Xint12ub  = 0.0;
      Strclp(Gradu,Cintu,Zwire12ub,Xwire12ub,Zint12ub,Xint12ub);
      float Dist12ub = sqrt((Zwire12ub-Zint12ub)*(Zwire12ub-Zint12ub) + (Xwire12ub-Xint12ub)*(Xwire12ub-Xint12ub));
      if(Dist12ub > 1.2*Size[2]) {
                                 return;
      }
      float Zpt12ub = 0.0;
      float Xpt12ub = 0.0;
      if(Ist2 > 0) {
        Strcht(zu1,xu1,zu2,xu2,Ist2,Layer2,Nmod,Zpt12ub,Xpt12ub);
      }
      if(Ist1 > 0) {
        Strcht(zu1,xu1,zu2,xu2,Ist1,Layer1,Nmod,Zpt12ub,Xpt12ub);
      }
      float Res12ub = sqrt((Zpt12ub - Zint12ub)*(Zpt12ub - Zint12ub) + (Xpt12ub - Xint12ub)*(Xpt12ub - Xint12ub));
      Res12ub = Res12ub*10000.0;
      float Res4vb = sqrt((Zpt4vb - Zint4vb)*(Zpt4vb - Zint4vb) + (Xpt4vb - Xint4vb)*(Xpt4vb - Xint4vb));
      Res4vb = Res4vb*10000.0;
      float Res3vb = sqrt((Zpt3vb - Zint3vb)*(Zpt3vb - Zint3vb) + (Xpt3vb - Xint3vb)*(Xpt3vb - Xint3vb));
      Res3vb = Res3vb*10000.0;
// ---------------------------------------------------------------------
// Convert all points back into xy frame for fitting:
// ---------------------------------------------------------------------
      float Zpt3b = Zpt3vb;
      float Xpt3b = 0.0;
      Bkv2xy(Xpt3vb,yv1,Xpt3b,yy1); 
//
      float Zpt4b = Zpt4vb;
      float Xpt4b = 0.0;
      Bkv2xy(Xpt4vb,yv2,Xpt4b,yy2);
//
      float Zpt12b = Zpt12ub;
      float Xpt12b = 0.0;
      Bku2xy(Xpt12ub,yu3,Xpt12b,yy3);
//
      Strflr(Ist3,Layer3,Nmod,Zpt3b,Xpt3b);
      Strflr(Ist4,Layer4,Nmod,Zpt4b,Xpt4b);
      if(Ist2 > 0) {
        Strflr(Ist2,Layer2,Nmod,Zpt12b,Xpt12b);
      }
      if(Ist1 > 0) {
        Strflr(Ist1,Layer1,Nmod,Zpt12b,Xpt12b);
      }
      Yhits[Ist3-1][Layer3-1][Nmod-1] = Ybest;
      Yhits[Ist4-1][Layer4-1][Nmod-1] = Ybest;
      if(Ist2 > 0) {
        Yhits[Ist2-1][Layer2-1][Nmod-1] = Ybest;
      }
      if(Ist1 > 0) {
        Yhits[Ist1-1][Layer1-1][Nmod-1] = Ybest;
      }
      Mask[Ist3-1][Layer3-1][Nmod-1] = 1;
      Mask[Ist4-1][Layer4-1][Nmod-1] = 1;
      if(Ist2 > 0) {
        Mask[Ist2-1][Layer2-1][Nmod-1] = 1;
      }
      if(Ist1 > 0) {
        Mask[Ist1-1][Layer1-1][Nmod-1] = 1;
      }
      Trkin3[0][Nmod-1]  = float(Ist1);
      Trkin3[1][Nmod-1]  = float(Ist2);
      Trkin3[2][Nmod-1]  = float(Ist3);
      Trkin3[3][Nmod-1]  = float(Ist4);
      float Zwire12b = 0.0;
      float Xwire12b = 0.0;
      float Zint12b  = 0.0;
      float Xint12b  = 0.0;
      if(Ist1 > 0) {
        Trkin3[4][Nmod-1]  = Hits[Ist1-1][Layer1-1][Nmod-1];
        Trkin3[5][Nmod-1]  = Zpt12b;
        Trkin3[6][Nmod-1]  = Xpt12b;
        Trkin3[7][Nmod-1]  = Ybest;
        Zwire12b = Zwire12ub;
        Xwire12b = 0.0;
        Bku2xy(Xwire12ub,yu3,Xwire12b,yy3);
        Zint12b = Zint12ub;
        Xint12b = 0.0;
        Bku2xy(Xint12ub,yu3,Xint12b,yy3);
        Trkin3[8][Nmod-1]  = Zwire12b;
        Trkin3[9][Nmod-1]  = Xwire12b;
        Trkin3[10][Nmod-1] = Zint12b;
        Trkin3[11][Nmod-1] = Xint12b;
      }
      if(Ist2 > 0) {
        Trkin3[12][Nmod-1] = Hits[Ist2-1][Layer2-1][Nmod-1];
        Trkin3[13][Nmod-1] = Zpt12b;
        Trkin3[14][Nmod-1] = Xpt12b;
        Trkin3[15][Nmod-1] = Ybest;
        Zwire12b = Zwire12ub;
        Xwire12b = 0.0;
        Bku2xy(Xwire12ub,yu3,Xwire12b,yy3);
        Zint12b  = Zint12ub;
        Xint12b  = 0.0;
        Bku2xy(Xint12ub,yu3,Xint12b,yy3);
        Trkin3[16][Nmod-1] = Zwire12b;
        Trkin3[17][Nmod-1] = Xwire12b;
        Trkin3[18][Nmod-1] = Zint12b;
        Trkin3[19][Nmod-1] = Xint12b;
      }
      Trkin3[20][Nmod-1] = Hits[Ist3-1][Layer3-1][Nmod-1];
      Trkin3[21][Nmod-1] = Zpt3b;
      Trkin3[22][Nmod-1] = Xpt3b;
      Trkin3[23][Nmod-1] = Ybest;
//
      float Zwire3b = Zwire3vb;
      float Xwire3b = 0.0; 
      Bkv2xy(Xwire3vb,yv1,Xwire3b,yy1);
      float Zint3b = Zint3vb;
      float Xint3b = 0.0;
      Bkv2xy(Xint3vb,yv1,Xint3b,yy1);
      Trkin3[24][Nmod-1] = Zwire3b;
      Trkin3[25][Nmod-1] = Xwire3b;
      Trkin3[26][Nmod-1] = Zint3b;
      Trkin3[27][Nmod-1] = Xint3b;
//
      Trkin3[28][Nmod-1] = Hits[Ist4-1][Layer4-1][Nmod-1];
      Trkin3[29][Nmod-1] = Zpt4b;
      Trkin3[30][Nmod-1] = Xpt4b;
      Trkin3[31][Nmod-1] = Ybest;
      float Zwire4b = Zwire4vb;
      float Xwire4b = 0.0; 
      Bkv2xy(Xwire4vb,yv2,Xwire4b,yy2);
      float Zint4b = Zint4vb;
      float Xint4b = 0.0; 
      Bkv2xy(Xint4vb,yv2,Xint4b,yy2);
      Trkin3[32][Nmod-1] = Zwire4b;
      Trkin3[33][Nmod-1] = Xwire4b;
      Trkin3[34][Nmod-1] = Zint4b;
      Trkin3[35][Nmod-1] = Xint4b;
// ------------------------------------------------------------------
// Finally make the fitted TE from the obtained space points.
// Must set up xzaty in xy frame for this:
// ------------------------------------------------------------------
      float Zzendc = Xzstraw[2][Ist3-1][2][Nmod-1];
      float Xxendc = Xzstraw[0][Ist3-1][2][Nmod-1];
      float Zzenda = Xzstraw[5][Ist3-1][2][Nmod-1];
      float Xxenda = Xzstraw[3][Ist3-1][2][Nmod-1];
      float Zzaty  = Zzendc;
      float Xxaty  = Xxendc + ((Xxenda-Xxendc)*(Ybest -(-4.55))/(9.1));
      Xzaty[0][Ist3-1][2][Nmod-1] = Zzaty;
      Xzaty[1][Ist3-1][2][Nmod-1] = Xxaty;
//
      Zzendc = Xzstraw[2][Ist4-1][3][Nmod-1];
      Xxendc = Xzstraw[0][Ist4-1][3][Nmod-1];
      Zzenda = Xzstraw[5][Ist4-1][3][Nmod-1];
      Xxenda = Xzstraw[3][Ist4-1][3][Nmod-1];
      Zzaty  = Zzendc;
      Xxaty  = Xxendc + ((Xxenda-Xxendc)*(Ybest -(-4.55))/(9.1));
      Xzaty[0][Ist4-1][3][Nmod-1] = Zzaty;   
      Xzaty[1][Ist4-1][3][Nmod-1] = Xxaty;  
      if(Ist1 > 0) {
        Zzendc = Xzstraw[2][Ist1-1][0][Nmod-1];
        Xxendc = Xzstraw[0][Ist1-1][0][Nmod-1];
        Zzenda = Xzstraw[5][Ist1-1][0][Nmod-1];
        Xxenda = Xzstraw[3][Ist1-1][0][Nmod-1];
        Zzaty  = Zzendc;
        Xxaty  = Xxendc + ((Xxenda-Xxendc)*(Ybest -(-4.55))/(9.1));
        Xzaty[0][Ist1-1][0][Nmod-1] = Zzaty;   
        Xzaty[1][Ist1-1][0][Nmod-1] = Xxaty;   
      }
      if(Ist2 > 0) {
        Zzendc = Xzstraw[2][Ist2-1][1][Nmod-1];
        Xxendc = Xzstraw[0][Ist2-1][1][Nmod-1];
        Zzenda = Xzstraw[5][Ist2-1][1][Nmod-1];
        Xxenda = Xzstraw[3][Ist2-1][1][Nmod-1];
        Zzaty  = Zzendc;
        Xxaty  = Xxendc + ((Xxenda-Xxendc)*(Ybest -(-4.55))/(9.1));
        Xzaty[0][Ist2-1][1][Nmod-1] = Zzaty;
        Xzaty[1][Ist2-1][1][Nmod-1] = Xxaty;
      }
//
      int Ntofit = Nposte;
      for(int ifit=1; ifit<=Nposte; ifit++) {
         Xptsfit[ifit-1] = Xypost[1][ifit-1][Nmod-1];
         Zptsfit[ifit-1] = Xypost[0][ifit-1][Nmod-1];
      }
      float Grfit = 0.0;
      float Cfit  = 0.0;
      Linfit(Zptsfit,Xptsfit,Ntofit,Grfit,Cfit);
      float zte1 = 0.0;
      float xte1 = 0.0;
      float zte2 = 0.0;
      float xte2 = 0.0;
      Strclp(Grfit,Cfit,Zpt12b,Xpt12b,zte1,xte1);
      Strclp(Grfit,Cfit,Zpt4b,Xpt4b,zte2,xte2);
//
      float xdiff  = xte2 - xte1;
      float zdiff  = zte2 - zte1;
      float Phizx  = Bkatan(xdiff,zdiff);
      float Phizxd = Phizx*180.0/Pibk;
      if(Phizxd > 180.0) {
                         Phizxd = Phizxd - 360.0;
      }
// --------------------------------------------
// Turn into milliradians:
// --------------------------------------------
      float Phizxmr = Phizxd*1000.0*Pibk/180.0;
//
      Trkin3[36][Nmod-1] = zte1;
      Trkin3[37][Nmod-1] = xte1;
      Trkin3[38][Nmod-1] = zte2;
      Trkin3[39][Nmod-1] = xte2;
      float Grrr = (xte2 - xte1) / (zte2 - zte1);
      float Ctrk = xte1 - (Grrr * zte1);
      Trkin3[40][Nmod-1] = Grrr;
      Trkin3[41][Nmod-1] = Ctrk;
// ------------------------------------------------------------------
// Finally store residuals to fitted track.
// ------------------------------------------------------------------
      float Zintt3 = 0.0;
      float Xintt3 = 0.0;
      Strclp(Grrr,Ctrk,Zwire3b,Xwire3b,Zintt3,Xintt3);
      float Zptt3  = 0.0;
      float Xptt3  = 0.0;
      Strcht(zte1,xte1,zte2,xte2,Ist3,Layer3,Nmod,Zptt3,Xptt3);
      float Zintt4 = 0.0;
      float Xintt4 = 0.0;
      Strclp(Grrr,Ctrk,Zwire4b,Xwire4b,Zintt4,Xintt4);
      float Zptt4  = 0.0;
      float Xptt4  = 0.0;
      Strcht(zte1,xte1,zte2,xte2,Ist4,Layer4,Nmod,Zptt4,Xptt4);
//
      float Zintt12 = 0.0;
      float Xintt12 = 0.0;
      Strclp(Grrr,Ctrk,Zwire12b,Xwire12b,Zintt12,Xintt12);
      float Zptt12 = 0.0;
      float Xptt12 = 0.0;
      if(Ist2 > 0) {
        Strcht(zte1,xte1,zte2,xte2,Ist2,Layer2,Nmod,Zptt12,Xptt12);
      }
      if(Ist1 > 0) {
        Strcht(zte1,xte1,zte2,xte2,Ist1,Layer1,Nmod,Zptt12,Xptt12);
      }
// ------------------------------------------------------------------
// (Zpt1,Xpt1) etc are the space points on the drift circles that 
// have been fitted to make the track.
// (Zptt1,Xptt1) are the closest points on the drift circles to the
// fitted track...
// ------------------------------------------------------------------
      float Rest3 = sqrt((Zptt3 - Zintt3)*(Zptt3 - Zintt3) + (Xptt3 - Xintt3)*(Xptt3 - Xintt3));
      Rest3 = Rest3*10000.0;
      float Rest4 = sqrt((Zptt4 - Zintt4)*(Zptt4 - Zintt4) + (Xptt4 - Xintt4)*(Xptt4 - Xintt4));
      Rest4 = Rest4*10000.0;
      float Rest12 = sqrt((Zptt12 - Zintt12)*(Zptt12 - Zintt12) + (Xptt12 - Xintt12)*(Xptt12 - Xintt12));
      Rest12 = Rest12*10000.0;
//
      if(Ist1 > 0) {
        Trkin3[42][Nmod-1] = Rest12;
      }
      if(Ist2 > 0) {
        Trkin3[43][Nmod-1] = Rest12;
      }
      Trkin3[44][Nmod-1] = Rest3;
      Trkin3[45][Nmod-1] = Rest4;
      if(Ist1 > 0) {
        Trkin3[46][Nmod-1] = Zptt12;
        Trkin3[47][Nmod-1] = Xptt12;
      }
      if(Ist2 > 0) {
        Trkin3[48][Nmod-1] = Zptt12;
        Trkin3[49][Nmod-1] = Xptt12;
      }
      Trkin3[50][Nmod-1] = Zptt3;
      Trkin3[51][Nmod-1] = Xptt3;
      Trkin3[52][Nmod-1] = Zptt4;
      Trkin3[53][Nmod-1] = Xptt4;
//
      Zxhit[0][Ist3-1][Layer3-1][Nmod-1] = Zptt3;
      Zxhit[1][Ist3-1][Layer3-1][Nmod-1] = Xptt3;
      Zxhit[0][Ist4-1][Layer4-1][Nmod-1] = Zptt4;
      Zxhit[1][Ist4-1][Layer4-1][Nmod-1] = Xptt4;
      if(Ist2 > 0) {
        Zxhit[0][Ist2-1][Layer2-1][Nmod-1] = Zptt12;
        Zxhit[1][Ist2-1][Layer2-1][Nmod-1] = Xptt12;
      }
      if(Ist1 > 0) {
        Zxhit[0][Ist1-1][Layer1-1][Nmod-1] = Zptt12;
        Zxhit[1][Ist1-1][Layer1-1][Nmod-1] = Xptt12;
      }
}
