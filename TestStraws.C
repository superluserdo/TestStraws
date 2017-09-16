 // ---------------------------------------------------------------------
// Version of TestStraws for use with data taken during the QA tests
// of each module. Altered by Tom (changes finished 16/09/2017).
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

#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <algorithm>
#include <getopt.h>
using namespace std;

#include <cmath>
#include <TF1.h>
#include <TF2.h>
#include <TGraph.h>
#include <TRandom3.h>
#include <TColor.h>
#include "TestStraws.h"

int n_HitsOut = 0;
int allow_3Hit = 0;
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

struct singleStraw straw;
struct strawEndsStruct strawEnds;
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
     float Xzstraw[6][STRAWS][LAYERS][MODULES];
      float strawGeom[MODULES][LAYERS][STRAWS][2][3];
      float strawGeom_u[MODULES][LAYERS][STRAWS][2][3];
      float strawGeom_v[MODULES][LAYERS][STRAWS][2][3];
	float (*strawGeoms[3])[MODULES][LAYERS][STRAWS][2][3] = {
	&strawGeom_u, &strawGeom_v, &strawGeom};
      float Diam;
      float Ylength;
      float Size[4];
//      float Xzstrau[6][STRAWS][LAYERS][MODULES];
//      float Xzstrav[6][STRAWS][LAYERS][MODULES];
      float Xzaty[2][STRAWS][LAYERS][MODULES];
      int   Iexist[STRAWS][LAYERS][MODULES];
      float Pedest[STRAWS][LAYERS][MODULES];
      float T0data;
      float Tofl;
      float Yped;
      float Hits[STRAWS][LAYERS][MODULES];
      float Dtimes[STRAWS][LAYERS][MODULES];
      float Tof[STRAWS][LAYERS][MODULES];
      float Rtimes[STRAWS][LAYERS][MODULES];
      int   Mask[STRAWS][LAYERS][MODULES];
      float Yhits[STRAWS][LAYERS][MODULES];
      float Zxhit[2][STRAWS][LAYERS][MODULES];
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

    struct options_struct options = {
	options.help = 0,
	options.file = 0,
	options.verbose = 0,
	options.output = 0
    };

  int create_gif = 1;
  int create_pdf = 1;
  int create_ps = 1;
  int create_hist = 1;

  int first = 1;
// -----------------------------------------------------------------------------------
    int main(int argc, char **argv) {
// -----------------------------------------------------------------------------------
      ofstream trackeventsfile;
      trackeventsfile.open ("trackevents.txt");
      trackeventsfile << "layout = {\n\n";
      	
      trackeventsfile << "	    modules = " << MODULES << "\n";
      trackeventsfile << "	    layers = " << LAYERS << "\n";
      trackeventsfile << "	    straws = " << STRAWS << "\n";
      trackeventsfile << "}\n\n";
      trackeventsfile << "trackevents = (\n\n";
      trackeventsfile << fixed;
      trackeventsfile << setprecision(5);

    TFile *inFile;
    inFile = new TFile("/hepstore/g2share/ProfTree/3ModuleData/Lab3TreeDump.root");
// Arguments
// ---------

  char *createopts;
  int c;
  int index;

  while (1)
    {
      static struct option long_options[] =
        {
          /* These options set a flag. */
          /* These options don’t set a flag.
             We distinguish them by their indices. */
          {"help",    no_argument,       0, 'h'},
          {"append",  no_argument,       0, 'b'},
          {"verbose", no_argument,       0, 'v'},
          {"output",  required_argument, 0, 'o'},
          {"file",    required_argument, 0, 'f'},
          {"create",  required_argument, 0, 'c'},
          {0, 0, 0, 0}
        };
      /* getopt_long stores the option index here. */
      int option_index = 0;

      c = getopt_long (argc, argv, "hbvf:o:c:",
                       long_options, &option_index);

      /* Detect the end of the options. */
      if (c == -1)
        break;

      switch (c)
        {
        case 0:
          /* If this option set a flag, do nothing else now. */
          if (long_options[option_index].flag != 0)
            break;
          printf ("option %s", long_options[option_index].name);
          if (optarg)
            printf (" with arg %s", optarg);
          printf ("\n");
          break;

        case 'h':
          printf(
"Teststraws options:\n"
"	-h	--help			print this help\n"
"	-f	--file	[path]		specify data file to use\n"
"	-v	--verbose		print output to stdout\n"
"	-o	--output [file]		output to file\n"
"	-c	--create [formats]	create file formats:\n"
"						g (gif)\n"
"						p (pdf)\n"
"						s (ps)\n"
"						h (hist)\n"
"					eg. `--create gph'\n");

          return 0;

        case 'v':
	  options.verbose = 1;
          puts ("`verbose' flag is set. Output from Teststraws will be limited.");
          break;

        case 'o':
          printf ("Sending output to file `%s'\n", optarg);
          break;

        case 'f':
          printf ("Using file `%s'\n", optarg);
          inFile = new TFile(optarg);
          if (inFile->IsZombie()) {
		printf("Cannot open file '%s' as a valid root file. Quitting...\n", optarg);
		exit(-1);
          }
          break;

        case 'c':
          create_gif = 0;
          create_pdf = 0;
          create_ps = 0;
          create_hist = 0;
          if (optarg) {
            createopts = optarg;
            for (index = 0; index < strlen(createopts); index++) {
  		switch (createopts[index]) {
                    case '0':
                      create_gif = 0;
                      create_pdf = 0;
                      create_ps = 0;
                      create_hist = 0;
                      printf("Not creating any files\n");
                      break;
                    case 'g':
                      create_gif = 1;
                      break;
                    case 'p':
                      create_pdf = 1;
                      break;
                    case 's':
                      create_ps = 1;
                      break;
                    case 'h':
                      create_hist = 1;
                      break;
                    default:
                      printf("Unrecognised format `%c'\n", createopts[index]);
                  }
            }
          }
          break;

        case '?':
          /* getopt_long already printed an error message. */
          break;

        default:
          abort ();
        }
    }

  /* Print any remaining command line arguments (not options). */
  if (optind < argc)
    {
      printf ("non-option ARGV-elements: ");
      while (optind < argc)
        printf ("%s ", argv[optind++]);
      putchar ('\n');
    }

  if (create_gif || create_pdf || create_ps || create_hist) {
    printf("Will create histograms in the following formats:\n");
    if (create_gif)
      printf("	.gif\n");
    if (create_pdf)
      printf("	.pdf\n");
    if (create_ps)
      printf("	.ps\n");
    if (create_hist)
      printf("	.hist\n");
  }
  //exit (0);

    if (inFile->IsZombie()) {
       printf("Something wrong with input file\n");
  	exit(-1);
    }

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
    zero_Init(i1);
    geom_Init();
// -------------------------------------------------------------------
    int nTotal = strawTree->GetEntries();
    if (options.verbose) {
        cout<<" Ntotal: "<<nTotal<<endl;
    }
// -------------------------------------------------------------------
// Loop round 10 times filling 100 hits into arrays.
// -------------------------------------------------------------------
    for(int iloop = 0; iloop < 3000; iloop++) {
       int istart = 100*(iloop);
       int iend   = istart + 100;
       for(int ia = 0; ia < 100; ia++) {
          Iruns[ia] = 0;
          Ievts[ia] = 0;
          TrigTimes[ia] = 0.0;
          Layers[ia] = 0;
          Istraw[ia] = 0;
          Module[ia] = 0;
          RawTimes[ia] = 0.0;
          Iviews[ia] = 0;
       }
//
       int iev = 0;
       for (int jentry = istart; jentry < iend; jentry++) {
           strawTree->GetEntry(jentry);
           scintTree->GetEntry(jentry);
           int jstraw = strawHit.wire;
           int iok = 1;
//         if(jstraw == 1) {
//                         iok = 0;
//         }
           if(iok == 1) {
             Iruns[iev]      = strawHit.run;
             Ievts[iev]      = strawHit.event;
             Istraw[iev]     = strawHit.wire;
             Layers[iev]     = strawHit.layer;
             Module[iev]     = strawHit.module;
             TrigTimes[iev]  = strawHit.triggerTime;
             RawTimes[iev]   = strawHit.hitTime; 
             Iviews[iev]     = strawHit.iview;
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
             if(Layers[iev] > -1) {
               if(Iviews[iev] == 1) {
                        Layers[iev] = Layers[iev] + 2;
               }
//
               Istraw[iev]++;
               Layers[iev]++;
               Module[iev]++;
               StrHist.fWire[0]->Fill(Istraw[iev]-0.5);
               StrHist.fModule[0]->Fill(Module[iev]-0.5);
               StrHist.fLayer[0]->Fill(Layers[iev]-0.5);
//
//             cout <<" Run "<<Iruns[iev] 
//                  <<" Straw "<<Istraw[iev]<<" Layer: "<<Layers[iev]<< " Module: "<< Module[iev]
//                  <<" Event: "<<Ievts[iev]<<" Trig: "<<TrigTimes[iev]
//                  <<" HitTime: "<<RawTimes[iev]<<endl;
//
//                  <<" Iview    "<<Iviews[iev]<<endl;
//
//             cout <<" Scint "<<jscint<<" Layer: "<<Layscint<< " Run: "<< jrunscint
//                  <<" Event: "<<jevscint<<" Trig: "<<Trigscint
//                  <<" HitTime: "<<Hitscint
//                  <<" Iview    "<<jviewscint<<endl;
//
//                int Iblob;
//                cin>>Iblob;
             }
             iev++;
//
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
          istep1++;
          Jcon1 = 1;
          JevNo = Ievts[istep1 - 1];
          NhitsInEvent++;
          Ngood++;
          ModEvent[Ngood - 1] = Module[istep1-1];
          LayEvent[Ngood - 1] = Layers[istep1-1];
          IstEvent[Ngood - 1] = Istraw[istep1-1];
          RawTevent[Ngood - 1]= RawTimes[istep1-1];
          istep2 = istep1;
//        cout<<" In outer loop: istep1 is: "<<istep1<<endl;
//        cout<<" JevNo is: "<<JevNo<<endl;
          do {
             Jcon2  = 0;
             istep2++;
//           cout<<" In inner loop: istep2 is: "<<istep2<<endl;
             JevNext = Ievts[istep2-1];
//           cout<<" JevNext is: "<<JevNext<<endl;
             if(JevNext == JevNo) {
                                Jcon2 = 1;
                                NhitsInEvent++;
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
                                  Ngood++;
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
          zero_Init(i2);
          Nhmods[0] = 0;
          Nhmods[1] = 0;
          Nhmods[2] = 0;
          for(int igood = 0; igood < Ngood; igood++) {
             int modd =  ModEvent[igood];
             Nhmods[modd-1]++;
             int layy =  LayEvent[igood];
             int istt =  IstEvent[igood]; 
             Rtimes[istt-1][layy-1][modd-1] = RawTevent[igood];
             float rrrr =  RawTevent[igood];
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
          for(int imod = 0; imod < MODULES; imod++) {
             if(Nhmods[imod] == 4) {
// Module has 4 hits:
               Nhlays[0] = 0;     
               Nhlays[1] = 0;
               Nhlays[2] = 0;
               Nhlays[3] = 0;
               int is1   = 0;
               int is2   = 0;
               int is3   = 0;
               int is4   = 0;
               for(int ist = 0; ist < 32; ist++) {
                  if(Rtimes[ist][0][imod] > 0.0) {
                    Nhlays[0]++;
                    is1 = ist;
                  }
                  if(Rtimes[ist][1][imod] > 0.0) {
                    Nhlays[1]++;
                    is2 = ist;
                  }
                  if(Rtimes[ist][2][imod] > 0.0) {
                    Nhlays[2]++;
                    is3 = ist;
                  }
                  if(Rtimes[ist][3][imod] > 0.0) {
                    Nhlays[3]++;
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
            for(int imod = 0;imod < MODULES; imod++) {
               for(int ilay = 0; ilay < LAYERS; ilay++) {
                  for(int istr = 0; istr < STRAWS; istr++) {
                     if(Rtimes[istr][ilay][imod] > 0.0) {
                       float rrr = Rtimes[istr][ilay][imod];
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
            straw_Evt(Goodt0);
            if(Goodt0 == 1) {
                           straw_Analyse(trackeventsfile);
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

// -------------------------------------------------------------------
// Close file
// -------------------------------------------------------------------
    trackeventsfile << "	}\n";
    trackeventsfile << ");\n";
    trackeventsfile.close();

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
       for(int i = 0; i < Nmax; i++) {
          StrHist.fWire[i] = new TH1F(Form(" Wire %d", i)," Wire ",33,0.0,33.0);
          StrHist.fWire[i]->SetFillColor(520-i+1);
          StrHist.fLayer[i] = new TH1F(Form(" Layer %d", i)," Layer ",5,0.0,5.0);
          StrHist.fLayer[i]->SetFillColor(519-i+1);
          StrHist.fModule[i] = new TH1F(Form(" Module %d", i)," Module ",4,0.0,4.0);
          StrHist.fModule[i]->SetFillColor(518-i+1);
          StrHist.fdiff12[i] = new TH1F(Form(" diff12 %d", i)," diff12 ",120,-60.0,60.0);
          StrHist.fdiff12[i]->SetFillColor(517-i+1);
          StrHist.fdiff34[i] = new TH1F(Form(" diff34 %d", i)," diff34 ",120,-60.0,60.0);
          StrHist.fdiff34[i]->SetFillColor(516-i+1);
          StrHist.fDtimes[i] = new TH1F(Form(" Dtimes %d", i)," Dtimes ",100,0.0,100.0);
          StrHist.fDtimes[i]->SetFillColor(515-i+1);
          StrHist.fDdists[i] = new TH1F(Form(" Ddists %d", i)," Ddists ",120,0.0,0.3);
          StrHist.fDdists[i]->SetFillColor(514-i+1);
          StrHist.fTrkMod[i] = new TH1F(Form(" TrkMod %d", i)," TrkMod ",4,0.0,4.0);
          StrHist.fTrkMod[i]->SetFillColor(504-i+1);
          StrHist.fYbest[i] = new TH1F(Form(" Ybest %d", i)," Ybest ",50,-5.0,5.0);
          StrHist.fYbest[i]->SetFillColor(513-i+1);
          StrHist.fGrfit[i] = new TH1F(Form(" Grfit %d", i)," Grfit ",80,-2.0,2.0);
          StrHist.fGrfit[i]->SetFillColor(512-i+1);
          StrHist.fPhizxd[i] = new TH1F(Form(" Phizxd %d", i)," Phizxd ",80,-40.0,40.0);
          StrHist.fPhizxd[i]->SetFillColor(511-i+1);
          StrHist.fZxte[i] = new TH2F(Form(" Zxte %d", i)," Zxte ",100,0.0,50.0,90,0.0,30.0);
          StrHist.fRest1[i] = new TH1F(Form(" Rest1 %d", i)," Rest1 ",100,-200.0,200.0);
          StrHist.fRest1[i]->SetFillColor(510-i+1);
          StrHist.fRest2[i] = new TH1F(Form(" Rest2 %d", i)," Rest2 ",100,-200.0,200.0);
          StrHist.fRest2[i]->SetFillColor(509-i+1);
          StrHist.fRest3[i] = new TH1F(Form(" Rest3 %d", i)," Rest3 ",80,-800.0,800.0);
          StrHist.fRest3[i]->SetFillColor(508-i+1);
          StrHist.fRest4[i] = new TH1F(Form(" Rest4 %d", i)," Rest4 ",80,-800.0,800.0);
          StrHist.fRest4[i]->SetFillColor(507-i+1);
          StrHist.fRest[i] = new TH1F(Form(" Rest %d", i)," Rest ",80,-800.0,800.0);
          StrHist.fRest[i]->SetFillColor(506-i+1);
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
        for(int i = 0; i < Nmax; i++) {
          StrHist.fWire[i]->Write();
          StrHist.fLayer[i]->Write();
          StrHist.fModule[i]->Write();
          StrHist.fdiff12[i]->Write();
          StrHist.fdiff34[i]->Write();
          StrHist.fDtimes[i]->Write();
          StrHist.fDdists[i]->Write();
          StrHist.fTrkMod[i]->Write();
          StrHist.fYbest[i]->Write();
          StrHist.fGrfit[i]->Write();
          StrHist.fPhizxd[i]->Write();
          StrHist.fZxte[i]->Write();
          StrHist.fRest1[i]->Write();
          StrHist.fRest2[i]->Write();
          StrHist.fRest3[i]->Write();
          StrHist.fRest4[i]->Write();
          StrHist.fRest[i]->Write();
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
     if (create_ps)
       c1->SaveAs("StrawsPage1.ps");
     if (create_pdf)
       c1->SaveAs("StrawsPage1.pdf");
     if (create_gif)
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
     if (create_ps)
       c1->SaveAs("StrawsPage2.ps");
     if (create_pdf)
       c1->SaveAs("StrawsPage2.pdf");
     if (create_gif)
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
    for(int iloc = 0; iloc < Len; iloc++) {
                                       Array[iloc] = Value;
    }
}
// -----------------------------------------------------------------------------
void Vzero(float* Array, int Len) {
// -----------------------------------------------------------------------------
    for(int iloc = 0; iloc < Len; iloc++) {
                                       Array[iloc] = 0.0;
    }
}
// -----------------------------------------------------------------------------
void straw_ClosestPoint(float Gr,float c1,float Xpt,float Ypt,float &Xinter,float &Yinter) {
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
void straw_TimeDist_1(float Time,float &Dist) {
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
void straw_TimeDist_2(float Time,float &Dist) {
// ------------------------------------------------------------------
// Routine to use the time - distance relationship.
// This version from Garfield: Ar/Eth 50:50  No magnetic field.
// Voltage = 1700V.
// ------------------------------------------------------------------
   float TimeArray[13];
   float DistArray[13];
//
   fill_TimeDist_arrays(TimeArray, DistArray, 1);
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
   for(int it = 0; it < 13; it++) {
      if(Time >= TimeArray[it] && Time < TimeArray[it+1]) {
                                                       iloct = it+1; 
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
void straw_TimeDist_3(float Time,float &Dist) {
// ------------------------------------------------------------------
// Routine to use the time - distance relationship.
// This version from Garfield: Ar/Eth 50:50  No magnetic field.
// Voltage = 1700V.
// ------------------------------------------------------------------
   float TimeArray[13];
   float DistArray[13];
//
   float Fac = 55.0/47.0;

   fill_TimeDist_arrays(TimeArray, DistArray, Fac);
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
   for(int it = 0; it < 13; it++) {
      if(Time >= TimeArray[it] && Time < TimeArray[it+1]) {
                                                       iloct = it+1;
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
void straw_FillArray(int Nstraw, int Layer, int Nmod, float zpt, float xpt) {
// ------------------------------------------------------------------
// Fills array with real space points.
// Stores space points on track:
// ------------------------------------------------------------------
      Nposte++;
      if(Nposte < 5) {
        Xypost[0][Nposte-1][Nmod] = zpt;
        Xypost[1][Nposte-1][Nmod] = xpt;
        Nspose[Nposte-1][Nmod]    = Nstraw;
        Nlpose[Nposte-1][Nmod]    = Layer;
//      cout<<" straw_FillArray: "<<Nposte<<" "<<Nstraw<<" "<<Layer<<" "<<Nmod<<" "<<zpt<<" "<<xpt<<endl;
      }
}
// --------------------------------------------------------------------------------------------------------
void straw_NearestTangent(float zz1,float xx1,float zz2,float xx2,int Nstraw,int Layer,int Nmod,float &Zpt, float &Xpt) {
// --------------------------------------------------------------------------------------------------------
// Compute the point(Zpt,Xpt) nearest to the tangent in cell Nstraw.
// The tangent is the line through ( zz1,xx1) , (zz2,xx2)
// ---------------------------------------------------------------------
      float Pibk  = 3.1415926535;
      float Twopi = 2.0*Pibk;
      float Piby2 = Pibk/2.0;
// ---------------------------------------------------------------------
//    if(Nstraw < 1 || Nstraw > 8) {
//        cout<<" In straw_NearestTangent: "<<Nstraw<<endl;
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
      straw_LocalAngle(Za,Xa,Zb,Xb,Nstraw,Layer,Nmod,Thet,Iflag);
      float z3l = Hits[Nstraw][Layer-1][Nmod]*cos(Thet);
      float x3l = Hits[Nstraw][Layer-1][Nmod]*sin(Thet);
// ---------------------------------------------------------------
// Convert space point back to GM2 coordinates:
// ---------------------------------------------------------------
      float Sth = sin(PI/4.0);//Sthpla;
      float Cth = cos(PI/4.0);//Cthpla;
      if(Iflag == 1) {
                     float The = 3 * PI / 8.0;//Thecol + (Piby2/4.0);
                     Sth = sin(The);
                     Cth = cos(The);
      }
      Zpt = z3l*Sth + x3l*Cth + Xzaty[0][Nstraw][Layer-1][Nmod];
      Xpt = x3l*Sth - z3l*Cth + Xzaty[1][Nstraw][Layer-1][Nmod];
}
// ------------------------------------------------------------------------------------------------------
void straw_LocalAngle(float zz1,float xx1,float zz2,float xx2,int Nst,int Nlay,int Nmod,float &Angle, int &Iflag) {
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
      float Sth = sin(PI/4.0);//Sthpla;
      float Cth = cos(PI/4.0);//Cthpla;
//
      float z1l = (zz1 - Xzaty[0][Nst][Nlay-1][Nmod])*Sth - (xx1 - Xzaty[1][Nst][Nlay-1][Nmod])*Cth;
      float x1l = (xx1 - Xzaty[1][Nst][Nlay-1][Nmod])*Sth + (zz1 - Xzaty[0][Nst][Nlay-1][Nmod])*Cth;
      float Dz  = (zz2 - zz1)*Sth - (xx2 - xx1)*Cth;
      float Dx  = (xx2 - xx1)*Sth + (zz2 - zz1)*Cth;
// -----------------------------------------------------------------------------------------
      if(fabs(Dx) < Dysmal) {
            Iflag = 1;
            float The = 3 * PI / 8.0;//Thecol + (Piby2/4.0);
            Sth = sin(The);
            Cth = cos(The);
            z1l = (zz1 - Xzaty[0][Nst][Nlay-1][Nmod])*Sth - (xx1 - Xzaty[1][Nst][Nlay-1][Nmod])*Cth;
            x1l = (xx1 - Xzaty[1][Nst][Nlay-1][Nmod])*Sth + (zz1 - Xzaty[0][Nst][Nlay-1][Nmod])*Cth;
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
void straw_CommonTangent(int Nmod,int Nst1,int Layer1,int Nst2,int Layer2,int Icand,
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
      float d1 = Hits[Nst1][Layer1][Nmod];
      float d2 = Hits[Nst2][Layer2][Nmod];
      float Xyphys1  = Size[Layer1];
      float Xyphys2  = Size[Layer2];
      if(d1 > Xyphys1) {
                       d1 = Xyphys1;
      }
      if(d2 > Xyphys2) {
                       d2 = Xyphys2;
      }
      float Zcell1 = Xzaty[0][Nst1][Layer1][Nmod];
      float Xcell1 = Xzaty[1][Nst1][Layer1][Nmod];
      float Zcell2 = Xzaty[0][Nst2][Layer2][Nmod];
      float Xcell2 = Xzaty[1][Nst2][Layer2][Nmod];
// ------------------------------------
// Distances between cell centres
// ------------------------------------
      float D = sqrt((Zcell1-Zcell2)*(Zcell1-Zcell2) + (Xcell1-Xcell2)*(Xcell1-Xcell2));
// --------------------------------------
// Angle of line joining cell centres
// --------------------------------------
      float Alpha = atan2((Xcell2-Xcell1) , (Zcell2-Zcell1));
      float Phi = 0.0;
      float Delta1 = 0.0;
      float Delta2 = 0.0;
      if(Icand == 2) {
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
      if(Icand == 3) {
        if(fabs(d1+d2) > D) {
          Iflag = 1;
          return;
        }
        Phi    = acos((d1+d2)/D);
        Delta1 = Alpha + Phi;
        Delta2 = Delta1 + Pibk;
      }
      if(Icand == 1) {
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
      if(Icand == 0) {
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
      for(int j = 0; j < Npts; j++) {
         if(Ypts[j] != 0.0) {
               Sumx  = Sumx + Xpts[j];
               Sumy  = Sumy + Ypts[j];
               Count++;
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
      for(int j = 0; j < Npts; j++) {
         if(Ypts[j] != 0.0) {
            float Scartx = Xpts[j] - Xmed;
            float Scarty = Ypts[j] - Ymed;
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
float phi_arctan(float py , float px) {
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
int round(float Var) {
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
      int Len = STRAWS*LAYERS;
      for(int imod = 0; imod < MODULES; imod++) {
         for(int ilay = 0; ilay < LAYERS; ilay++) {
            for(int istr = 0; istr < STRAWS; istr++) {
               Pedest[istr][ilay][imod] = Value;
            }
         }
      }
      T0data = 0.0;
      Yped   = 0.33;
//    Tofl   = Zpos / 30.0;
      Tofl   = 0.0;
}
// --------------------------------------------------------------------
void straw_Analyse(ofstream &trackeventsfile) {
// --------------------------------------------------------------------
// Having obtained the drift times and hits now find the track.
// --------------------------------------------------------------------
// Start by assuming the wire positions for the middle of the cells:
// This gets changed in any case in straw_FindTrack_4 etc.
// ------------------------------------------------------------------
      for(int imod = 0; imod < MODULES; imod++) {
         for(int ilay = 0; ilay < LAYERS; ilay++) {
            for(int istr = 0; istr < STRAWS; istr++) {
// z coordinate:
               Xzaty[0][istr][ilay][imod] = strawGeom[imod][ilay][istr][0][2];
// x coordinate half way up the straw:
               Xzaty[1][istr][ilay][imod] = strawGeom[imod][ilay][istr][0][0]
						+ strawGeom[imod][ilay][istr][1][0];
            }
         }
      }
// ------------------------------------------------------------------
// Now find 4-hit track elements in each of the Straw Modules:
// ------------------------------------------------------------------
      int iStraws[LAYERS] = {0};
      for(int imod = 0; imod < MODULES; imod++) {
         straw_FindEvent(imod, trackeventsfile);
// ------------------------------------------------------------------
// Now look for 3-hit track elements.
// ------------------------------------------------------------------
//       straw_FindTrack_3(imod);
      }
}
// -------------------------------------------------------------------
void coordConvert(int uvType, int direction, float xIn , float yIn , float &xOut , float &yOut) {
// --------------------------------------------------------------------
// Translate a point in the u OR v coordinate system into the xy frame,
// or vice versa.
// uvType = 0 (u) or 1 (v)
// direction = 0 (u/v to xy) or 1 (xy to u/v)
// --------------------------------------------------------------------
     float Pibk = 3.1415926535;
     float Piby2 = Pibk/2.0;
// --------------------------------------------------------------------
// Xwire Ywire Theta are the positions of the local frames origin and 
// orientation as expressed in the global frame.
// --------------------------------------------------------------------
     float Xwire  = 0.0;
     float Ywire  = 0.0;
     float Thetad;
     if (uvType == 0) { // u
         Thetad = 82.5;
     }
     else { // v
         Thetad = 97.5;
     }
     float Theta  = Thetad * Pibk/180.0;
     float Cth    = cos(Theta-Piby2);
     float Sth    = sin(Theta-Piby2);
     if (direction == 0) { // xIn, yIn are in u or v
                           // xOut, yOut are in xy
// --------------------------------------------------------------
// Now translate (Xlocal,Ylocal) into the global frame.
// --------------------------------------------------------------
         xOut = xIn*Cth - yIn*Sth + Xwire;
         yOut = yIn*Cth + xIn*Sth + Ywire;
     }
// --------------------------------------------------------------
// Translate into the local (u or v) frame.
// --------------------------------------------------------------
     else { // xIn, yIn are in xy
            // xOut, yOut are u or v
         xOut = ( xIn - Xwire )*Cth + ( yIn - Ywire )*Sth; 
         yOut = ( yIn - Ywire )*Cth - ( xIn - Xwire )*Sth;
     }
}
float randomFloat() {
// ----------------------------------------------------------------
      IcallRndm++;
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
void zero_Init(int Imode) {
// --------------------------------------------------------------------
// Zero counters at run start.
// --------------------------------------------------------------------
     if(Imode == 1) {
                    Zpos = 0.0;
		for(int module = 0; module < MODULES; module++) {
			for(int layer = 0; layer < LAYERS; layer++) {
				for(int straw = 0; straw < STRAWS; straw++) {
					for(int ia=0; ia<2; ia++) {
						for(int ib=0; ib<3; ib++) {
							strawGeom[module][layer][straw][ia][ib] = 0;
							strawGeom_u[module][layer][straw][ia][ib] = 0;
							strawGeom_v[module][layer][straw][ia][ib] = 0;
						}
					}
				}
			}
		}
						
                    Diam    = 0.0;
                    Ylength = 0.0;
                    for(int ic = 0; ic < MODULES; ic++) {
                        for(int ib = 0; ib < LAYERS; ib++) {
                            for(int ia = 0; ia < STRAWS; ia++) {
                             Iexist[ia][ib][ic] = 0;
                             Pedest[ia][ib][ic] = 0.0;
                          }                      
                       }
                    }
                    T0data  = 0.0;
                    Tofl    = 0.0;
                    Yped    = 0.0;
                    for(int ic = 0; ic < MODULES; ic++) {
                        for(int ib = 0; ib < LAYERS; ib++) {
                            for(int ia = 0; ia < STRAWS; ia++) {
                             Tof[ia][ib][ic] = 0.0;
                          }
                       }
                    }
     }
// --------------------------------------------------------------------
// Zero counters for new event:
// --------------------------------------------------------------------
   for(int imod = 0; imod < MODULES; imod++) {
       for(int ilay = 0; ilay < LAYERS; ilay++) {
           for(int istr = 0; istr < STRAWS; istr++) {
               Rtimes[istr][ilay][imod] = 0.0;
               Dtimes[istr][ilay][imod] = 0.0;
               Hits[istr][ilay][imod]   = 0.0;
               Yhits[istr][ilay][imod]  = 0.0;
               Mask[istr][ilay][imod]   = 0;
               Zxhit[0][istr][ilay][imod] = 0.0;
               Zxhit[1][istr][ilay][imod] = 0.0;
            }
            Xypost[0][imod][imod] = 0.0;
            Xypost[1][imod][imod] = 0.0;
            Nspose[imod][imod]    = 0;
            Nlpose[imod][imod]    = 0;
         }
      }
//
      Nposte = 0;
      for(int imod = 0; imod < MODULES; imod++) {
         for(int id = 0; id < 60; id++) {
            Trkinf[id][imod] = 0.0;
            Trkin2[id][imod] = 0.0;
            Trkin3[id][imod] = 0.0;
         }
      } 
}
// --------------------------------------------------------------------
void geom_Init() {
// --------------------------------------------------------------------
     float xx = 0.0;
     float yy = 0.0;
     float xu = 0.0;
     float yu = 0.0;
     float xv = 0.0;
     float yv = 0.0;


	float theta = 7.5*PI/180;
	float wdist = 0.6/(cos(theta));
	float Xoff1 = 0.0; /* Allow for the module to be offset in X */
	straw.Ybot = -4.55;
	straw.Ytop = 4.55;
	float Zpos = 10.0;
	float Zmoddiff = 13.0;
	float layeroffsets_Z[LAYERS] = {3.931, 4.408, 6.091, 6.568};
	float layeroffsets_X[LAYERS] = {4.7984, 4.7984 + wdist/2, 6.378, 6.378 - wdist/2};
	float stereodiff;
	stereodiff = 9.1 * tan(theta);
     Diam     = 0.5;
     for(int ilay = 0; ilay < LAYERS; ilay++) {
        Size[ilay] = Diam / 2.0;
     }

	/* NOTE: in GM2 coordinates, "Z" is along beam line,
	 * which is drawn HORIZONTALLY here!
	 * Likewise, "X" is orthogonal, and is drawn vertically
	 * on the screen!
	 *
	 * I am using capital notation (X, Y, Z) for beamline
	 * coordinates and lowercase (x, y) for program
	 * coordinates.
	 */

	for (int imod = 0; imod < MODULES; imod ++) {
		for (int ilay = 0; ilay < LAYERS; ilay ++) {
			for (int istr = 0; istr < STRAWS; istr++) {

				/*	X	*/

				strawGeom[imod][ilay][istr][0][0] = Xoff1  + layeroffsets_X[ilay] + istr * wdist;
				if (ilay == 0 || ilay == 1)
					strawGeom[imod][ilay][istr][1][0] = Xoff1  + layeroffsets_X[ilay] + istr * wdist + stereodiff;
				else
					strawGeom[imod][ilay][istr][1][0] = Xoff1  + layeroffsets_X[ilay] + istr * wdist - stereodiff;

				/*	Y	*/
				
				strawGeom[imod][ilay][istr][0][1] = straw.Ybot;
				strawGeom[imod][ilay][istr][1][1] = straw.Ytop;
				
				/*	Z	*/
				
				strawGeom[imod][ilay][istr][0][2] = Zpos + Zmoddiff * imod + layeroffsets_Z[ilay];
				strawGeom[imod][ilay][istr][1][2] = Zpos + Zmoddiff * imod + layeroffsets_Z[ilay];

// ---------------------------------------------------------------------
// Now calculate positions of the wire strawEnds.in the u and v frames.
// Z positions are the same in all frames:
// ----------------------------------------------------------------------
               strawGeom_u[imod][ilay][istr][0][2] = strawGeom[imod][ilay][istr][0][2];
               strawGeom_u[imod][ilay][istr][1][2] = strawGeom[imod][ilay][istr][1][2];
               strawGeom_v[imod][ilay][istr][0][2] = strawGeom[imod][ilay][istr][0][2];
               strawGeom_v[imod][ilay][istr][1][2] = strawGeom[imod][ilay][istr][1][2];
// ----------------------------------------------------------------------
                for(int bottop = 0; bottop < 2; bottop++) {	
           // Straw in u frame:  
                      xx = strawGeom[imod][ilay][istr][bottop][0];
                      yy = strawGeom[imod][ilay][istr][bottop][1];
                      coordConvert(0, TO_UV, xx,yy,xu,yu);
                      strawGeom_u[imod][ilay][istr][bottop][0] = xu;
                      strawGeom_u[imod][ilay][istr][bottop][1] = yu;
           // Straw in v frame:  
                      xx = strawGeom[imod][ilay][istr][bottop][0];
                      yy = strawGeom[imod][ilay][istr][bottop][1];
                      coordConvert(1, TO_UV, xx,yy,xv,yv);
                      strawGeom_v[imod][ilay][istr][bottop][0] = xv;
                      strawGeom_v[imod][ilay][istr][bottop][1] = yv;
                }
            }
         }
     }
 
// ----------------------------------------------------------
// Now get all the calibration constants for the STRAWs.
// ----------------------------------------------------------
     Strcal();
}			

// --------------------------------------------------------------------
void straw_Evt(int &Goodt0) {
// ----------------------------------------------------------
// Starting from the Raw time values ... first subtract
// pedestals and time of flight.
// Then try to determine t0 for this event:
// straw_Evt assumes a hit in each layer has already been
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
      int Nest   = 0;
      float t0aver = 0.0;
//
      for(int imod = 0; imod < MODULES; imod++) {
         dt1[imod] = -999.0;
         dt2[imod] = -999.0;
         dt3[imod] = -999.0;
         dt4[imod] = -999.0;
         is1[imod] = -999;
         is2[imod] = -999;
         is3[imod] = -999;
         is4[imod] = -999;
      float Dvel   = 4.7/1000.0;
// ----------------------------------------------------------
// Once this routine has been called there should be at most
// one hit in each layer - raw time values in Rtimes array.
// ----------------------------------------------------------
    	for (int ilay = 0; ilay < LAYERS; ilay ++) {
		for (int istr = 0; istr < STRAWS; istr++) {
               Tof[istr][ilay][imod] = Tofl;
               if(Rtimes[istr][ilay][imod] > 0.0) {
                 Dtimes[istr][ilay][imod] = Rtimes[istr][ilay][imod] - Pedest[istr][ilay][imod] - Tofl;
                 if(Dtimes[istr][ilay][imod] <= 0.0) {
                                                Dtimes[istr][ilay][imod] = 1.0;
                 }
               }
            }
         }
// ------------------------------------------------------------------
// Now try to determine the t0 for this event:
// ------------------------------------------------------------------
      T0data = 0.0;
//
	 for (int istr = 0; istr < STRAWS; istr++) {
            if(Dtimes[istr][0][imod] > 0.0) {
                 dt1[imod] = Dtimes[istr][0][imod];
                 is1[imod] = istr;
            }
            if(Dtimes[istr][1][imod] > 0.0) {
                 dt2[imod] = Dtimes[istr][1][imod];
                 is2[imod] = istr;
            }
            if(Dtimes[istr][2][imod] > 0.0) {
                 dt3[imod] = Dtimes[istr][2][imod];
                 is3[imod] = istr;
            }
            if(Dtimes[istr][3][imod] > 0.0) {
                 dt4[imod] = Dtimes[istr][3][imod];
                 is4[imod] = istr;
            }
         }
// ------------------------------------------------------------------
// We now have the times in each of the 4 layers in each module.
// Use them to obtain 2 estimates of t0:
// Use the fact that we know the sum of the drift times in adjacent
// layers.
// ------------------------------------------------------------------
         if(abs(is1[imod]-is2[imod]) < 2) {
           if(dt1[imod] > 0.0 && dt2[imod] > 0.0) {
             Nest++;
             t0est[Nest-1] = (dt1[imod] + dt2[imod] - 59.0)/2.0;
             t0aver = t0aver + t0est[Nest-1];
           }
         }
         if(abs(is3[imod]-is4[imod]) < 2) {
           if(dt3[imod] > 0.0 && dt4[imod] > 0.0) {
             Nest++;
             t0est[Nest-1] = (dt3[imod] + dt4[imod] - 59.0)/2.0;
             t0aver = t0aver + t0est[Nest-1];
           }
         }
      }
      float t0diff = t0est[0] - t0est[1];
      Goodt0 = 0;
      if(Nest == 2 && fabs(t0diff) < 5.0) {
                                          Goodt0 = 1;
      }
      if(Nest > 0) {
                   t0aver = t0aver/float(Nest);
      }
      T0data = t0aver;
// ------------------------------------------------------------------
// Subtract the obtained t0 from the drift times:
// ------------------------------------------------------------------
     for(int imod = 0; imod < MODULES; imod++) {
    	for (int ilay = 0; ilay < LAYERS; ilay ++) {
		for (int istr = 0; istr < STRAWS; istr++) {
               if(Dtimes[istr][ilay][imod] > 0.0) {
                 Dtimes[istr][ilay][imod] = Dtimes[istr][ilay][imod] - T0data;
                 Dtimes[istr][ilay][imod] = Dtimes[istr][ilay][imod] - 4.0;
                 if(Dtimes[istr][ilay][imod] <= 0.0) {
                         Dtimes[istr][ilay][imod] = 1.0;
                 }
                 float Time = Dtimes[istr][ilay][imod];
                 float Dist = 0.0;
                 straw_TimeDist_3(Time,Dist);
                 if(Dist > 0.25) {
                              Dist = 0.25;
                 }
                 Hits[istr][ilay][imod] = Dist;
               }
            }
         }
      }     
}
// --------------------------------------------------------------------
void writeEvents() {
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
}
// -----------------------------------------------------------------
void straw_FindEvent(int Nmod, ofstream &trackeventsfile) {
// ------------------------------------------------------------------
// Finds a TE in a Straw Module where hits exist either in all
// four layers, or all but one (if allow_3Hit is enabled). 
// The TE is created if the position of the track in the layer 
// without a hit is compatible with having gone through the cell 
// wall.
// ------------------------------------------------------------------
      int exclude = 0;
      int inc[LAYERS] = {0};
      int side = 0;

      float Xptsfit[4];
      float Zptsfit[4];
//
      float Zw[4];
      float Xw[4];
      float Dr[4];
//
      float Pibk  = 3.1415926535;
      float Twopi = 2.0*Pibk;
      float Piby2 = Pibk/2.0;
      float Ysize = 0.01;
      float Ystep = 10.0/Ysize;
      float Yloop = 4.0 * Ystep;
      int   Kystep = round(Ystep);
      int   Kyloop = round(Yloop);
      float Resols[Kyloop];
      float Angloc[Kyloop];
      int   IsWall[Kyloop];
//
      for(int i=0; i < Kyloop; i++) {
         Resols[i] = 99999999.0;
         Angloc[i] = 99999999.0;
         IsWall[i] = 0;
      }
// -----------------------------------------
// Reset the number of TE points:
// -----------------------------------------
//       Nposte = 0;
//
// ------------------------------------------------------------------
// Declare some variables for later use.
// ------------------------------------------------------------------
      straw.Xbot = 0.0;
      straw.Ybot = 0.0;
      straw.Zbot = 0.0;
      straw.Xtop = 0.0;
      straw.Ytop = 0.0;
      straw.Ztop = 0.0;

	float coords[3][3][LAYERS] = {0};
//		     ^  ^  ^
//	x/y/z--------|  |  |----Layer (1/2/3/4)
//			|---Reference frame (u/v/Normal)
//
//	eg. what was yu3 is now coords[<y>][<u>][<2>] = coords[1][0][3]

      float Grads[2]; // = {Gradu. Gradv}
      float Cints[2]; // = {Cintu. Cintv}
      float Angexp = 0.0;

// ------------------------------------------------------------------
// Find the hits in each of the layers:
// ------------------------------------------------------------------

      int Layer1 = 1;
      int Layer2 = 2;
      int Layer3 = 3;
      int Layer4 = 4;
      int iStraws[LAYERS] = {0};
      int Nhits = 0;
      for(int ilay=0; ilay< LAYERS; ilay++) {
          for(int istr = 0; istr < STRAWS; istr++) {
             if(Hits[istr][ilay][Nmod] > 0.0) {
                                       iStraws[ilay]  = istr;
                                       Nhits++;
                                       inc[ilay] = 1;
             }
          }
      }
// ------------------------------------------------------------------
// If we dont have hits in all four layers then return.
// ------------------------------------------------------------------
 
      if(Nhits != 4) {
          if (Nhits == 3) {
              exclude = 1;
              int excluded;
              for (int ilay = 0; ilay < LAYERS; ilay++) {
                  if (inc[ilay] = 0) {
                     excluded = ilay;
                  }
              }
              side = excluded/2;
              if (!allow_3Hit) {
                  return;
              }
          }
          else {
              return;
          }
      }
// ------------------------------------------------------------------
// In a module with 4 trackable hits: Make some basic plots.
// ------------------------------------------------------------------
 
      float diff12 = Rtimes[iStraws[0]][0][Nmod] - Rtimes[iStraws[1]][1][Nmod];
      float diff34 = Rtimes[iStraws[2]][2][Nmod] - Rtimes[iStraws[3]][3][Nmod];
      StrHist.fdiff12[0]->Fill(diff12);
      StrHist.fdiff34[0]->Fill(diff34);
      StrHist.fDtimes[0]->Fill(Dtimes[iStraws[0]][0][Nmod]);
      StrHist.fDtimes[0]->Fill(Dtimes[iStraws[1]][1][Nmod]);
      StrHist.fDtimes[0]->Fill(Dtimes[iStraws[2]][2][Nmod]);
      StrHist.fDtimes[0]->Fill(Dtimes[iStraws[3]][3][Nmod]);
      StrHist.fDdists[0]->Fill(Hits[iStraws[0]][0][Nmod]);
      StrHist.fDdists[0]->Fill(Hits[iStraws[1]][1][Nmod]);
      StrHist.fDdists[0]->Fill(Hits[iStraws[2]][2][Nmod]);
      StrHist.fDdists[0]->Fill(Hits[iStraws[3]][3][Nmod]);

// ------------------------------------------------------------------
// Do track finding at various positions vertically and minimise the
// residuals:
// ------------------------------------------------------------------
      int Iybin_Max = 910; //Was 91 for 3 hit case
      float Iybin_Div = 100.0; // Was 10.0 for 3 hit case
      for(int Iybin=1; Iybin<=Iybin_Max; Iybin++) {
         float Ytest = float(Iybin)/Iybin_Div;
         Ytest = Ytest - 4.55;
	getEndsLoop(Nmod, Ytest, NULL, 1);
// ------------------------------------------------------------------
// Have found the hit straws in layers 1 and 2.
// Find the 4 possible tangents to these two drift circles.
// Extrapolate to layers 3 and 4 ...
// ------------------------------------------------------------------
 
        for (int i = 0; i < LAYERS; i++) {
            if (inc[i]) {
                Zw[i] = Xzaty[0][iStraws[i]][i][Nmod];
                Xw[i] = Xzaty[1][iStraws[i]][i][Nmod];
                Dr[i] = Hits[iStraws[i]][i][Nmod];
                straw.Ybot = (*strawGeoms[i/2])[Nmod][i][iStraws[i]][0][1];
                straw.Ytop = (*strawGeoms[i/2])[Nmod][i][iStraws[i]][1][1];
                coords[1][i/2][i] = straw.Ybot + ((straw.Ytop-straw.Ybot)*(Ytest -(-4.55))/(9.1));
            }
        }
         for(int itan = 0; itan < LAYERS; itan++) {
// ------------------------------------------------------------------
// Routine straw_CommonTangent returns the Z and X points of the common
// tangents to the drift circles for each of the 4 cases in turn.
// ------------------------------------------------------------------
//   Itan = 1 MEANS OUTER EDGE TANGENT , RHS OF DRIFT CIRCLES
//       = 2        ''    ''     ''     LHS  ''  ''     ''
//       = 3       CROSSED TANGENT , BOTTOM RHS TO TOP LHS
//       = 4           ''     ''       ''   LHS        RHS
// ------------------------------------------------------------------
            int iflag = 0;
            straw_CommonTangent(Nmod,iStraws[0],side,iStraws[1],side+1,itan,
		coords[2][side][0],
		coords[0][side][0],
		coords[2][side][1],
		coords[0][side][1],
		iflag);
            int Irloc = (Iybin-1)*4 + itan;
            Grads[1-side] = 9999.9;
            int Icon[4];
            Icon[0] = 1;
            if(fabs(coords[2][side][1] - coords[2][side][0]) > 1.0e-6) {
                    Grads[side] = 
                            (coords[0][side][1] - coords[0][side][0])/
                            (coords[2][side][1] - coords[2][side][0]);
                    for (int i = 0; i < 2; i++) {
                        coordConvert(side, TO_XY,
                             coords[0][side][i],
                             coords[1][side][i],
                             coords[0][2][i],
                             coords[1][2][i]);
                    }
                    coords[2][2][1] = coords[2][side][1];
                    coords[2][2][0] = coords[2][side][0];
                    Angloc[Irloc] = phi_arctan((coords[0][2][1]-coords[0][2][0]),
                                      (coords[2][2][1]-coords[2][2][0]));
                    if(Angloc[Irloc] < 0.0) {
                      Angloc[Irloc] = Angloc[Irloc] + Twopi;
                    }
                    Angloc[Irloc] = Angloc[Irloc]*180.0/Pibk;
                    Icon[0] = 1;
            }
            else {
                    Resols[Irloc] = 99999.9;
                    Icon[0] = 0;
            }

            if(Icon[0] == 1) {
              Cints[side] = coords[0][side][0] - (Grads[side] * coords[2][side][0]);
// -------------------------------------------------------------------
// Sanity check: Residuals in layers used to define the line should
// be zero.
// -------------------------------------------------------------------
              float wire[2][4]; // [X/Z][1u/2u/3v/4v]
              float inter[2][4]; // [X/Z][1u/2u/3v/4v]
              float pt[2][4] = {0}; // [X/Z][1u/2u/3v/4v]

              float res[4]; // [1u,2u,3v,4v]

              for (int j = 0; j < 2; j++) {
                  wire[1][side*2+j] = Zw[j];
                  wire[0][side*2+j] = Xw[j];
                  inter[1][side*2+j]  = 0.0;
                  inter[0][side*2+j]  = 0.0;
                  straw_ClosestPoint(Grads[side],Cints[side],
                                     wire[1][side*2+j],
                                     wire[0][side*2+j],
                                     inter[1][side*2+j],
                                     inter[0][side*2+j]);
              }
              for (int j = 0; j < 2; j++) {
                  straw_NearestTangent(coords[2][side][0],
                                       coords[0][side][0],
                                       coords[2][side][1],
                                       coords[0][side][1],
                                       iStraws[side*2+j], 2*side+1+j, Nmod,
                                       pt[1][side*2+j],
                                       pt[0][side*2+j]);
    
                  res[side*2+j] = sqrt( pow(pt[1][2*side+j]-inter[1][2*side+j], 2)
                                    + pow(pt[0][2*side+j]-inter[0][2*side+j], 2));
                  res[side*2+j] *= 10000.0;
              }
// -------------------------------------------------------------------
// Find closest point to this line on the hit cell in layer 3 or 4 
// We need to swap the extrapolated tangent into the v frame from the
// u frame:
// -------------------------------------------------------------------
 
              for (int i = 0; i < 2; i++) {
                  coords[2][1-side][i] = coords[2][side][i];
                  coordConvert(side, TO_XY,
                               coords[0][side][i],
                               coords[1][side][i],
                               coords[0][2][i],
                               coords[1][2][i]);
                  coordConvert(1-side, TO_UV,
                               coords[0][2][i],
                               coords[1][2][i],
                               coords[0][1-side][i],
                               coords[1][1-side][i]);
              }

              Grads[1-side] =(coords[0][1-side][1] - coords[0][1-side][0])/
                             (coords[2][1-side][1] - coords[2][1-side][0]);
              Cints[1-side] = coords[0][1-side][0] - (Grads[1-side] * coords[2][1-side][0]);

              for (int i = 0; i < 2; i++) {
                  wire[1][2-2*side+i] = Zw[2-2*side+i];
                  wire[0][2-2*side+i] = Xw[2-2*side+i];
                  inter[1][2-2*side+i] = 999.9;
                  inter[0][2-2*side+i] = 999.9;
                  straw_ClosestPoint(Grads[1-side],Cints[1-side],
                                     wire[1][2-2*side+i],wire[0][2-2*side+i],
                                     inter[1][2-2*side+i],inter[0][2-2*side+i]);
                  float Dist[4];
                  Dist[2-2*side+i] = sqrt(pow(wire[1][2-2*side+i] - inter[1][2-2*side+i], 2)
                                        + pow(wire[0][2-2*side+i] - inter[0][2-2*side+i], 2));
                  Icon[2-2*side+i] = 1;
                  if(Dist[2-2*side+i] > 1.2*Size[2-2*side+i]) {
                                      Resols[Irloc] = 99999.9;
                                      Icon[2-2*side+i] = 0;
                  }
               }

// -------------------------------------------------------------------
// Extrapolated line passes through the hit cell in layer 3 or 4:
// Calculate the residuals:
// -------------------------------------------------------------------
              int condition;
              float Restot = 0.0;
              if (!exclude) {
                  condition = (Icon[2] == 1 && Icon[3] == 1);
              }
              else {
                  condition = (Icon[2-2*side] == 1);
              } 
              if(condition) {
                for (int j = 0; j < 2; j++) {
                    if(inc[2-2*side + j]) {
                      straw_NearestTangent(coords[2][1-side][0],
                                           coords[0][1-side][0],
                                           coords[2][1-side][1],
                                           coords[0][1-side][1],
                                           iStraws[2-2*side+j], 2-2*side+1+j, Nmod,
                                           pt[1][2-2*side+j],
                                           pt[0][2-2*side+j]);
                      res[2-2*side+j] = sqrt(pow(pt[1][2-2*side+j] - inter[1][2-2*side+j], 2)
                                           + pow(pt[0][2-2*side+j] - inter[0][2-2*side+j], 2));
                      res[2-2*side+j] *= 10000.0;
                      Restot += res[2-2*side+j];
                    }
                }
                Resols[Irloc] = Restot;
                if (exclude) {
// -------------------------------------------------------------------
// We now need to decide if the layer that didnt have a hit is
// compatible with the track having gone through a cell wall.
// -------------------------------------------------------------------
                    IsWall[Irloc] = 1;
                    int Istart = 0;
                    int Iend   = 0;
                    for (int k = 0; k < 2; k++) {
                        if(iStraws[2-2*side + k] + 1 > 0) {
// ------------------
// No hit in Layer 4:
// ------------------
                        Istart = iStraws[2-2*side + k] - 1;  
                        if(Istart < 0) {
                                     Istart = 0;
                        }
                        Iend = iStraws[2-2*side + k] + 3;
                        if(Iend > 32) {
                                    Iend = 32;
                        }
                            for(int Is = Istart; Is < Iend; Is++) {
                               float Zwire = Xzaty[0][Is][3-2*side-k][Nmod];
                               float Xwire = Xzaty[1][Is][3-2*side-k][Nmod];
                               float Zint = 0.0;
                               float Xint = 0.0;
                               straw_ClosestPoint(Grads[1-side],Cints[1-side],Zwire,Xwire,Zint,Xint);
                               float Dist = sqrt((Zwire-Zint)*(Zwire-Zint) + (Xwire-Xint)*(Xwire-Xint));
                               if(Dist < 0.25) {
                                 IsWall[Irloc] = 0;
                                }
                            }
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
      float Resid;
      float Angte;
      float Angdif;
// ------------------------------------------------------------------
// Keep solution with smallest Residual that is most 
// consistent with track going straight through.
// ------------------------------------------------------------------
      if (exclude) {
          for(int isol = 0; isol < Kyloop; isol++) {
             Resid = Resols[isol];
             Angte = Angloc[isol];
             if(IsWall[isol] == 1) {
               Angdif = fabs(Angte - Angexp);
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
          if(Resbst > 1000.0) {
                              return;
          }
      }
      else {
          Locbst = 0;
          Resbst = 99999.9;
          for(int isol = 0; isol < Kyloop; isol++) {
             Resid = Resols[isol];
             if(Resid < Resbst) {
                                Resbst = Resid;
                                Locbst = isol;
             }
          }
      }
// ------------------------------------------------------------------
// Refit the chosen solution:
// ------------------------------------------------------------------
      int Iybest  = (Locbst)/4;
      int Itnbst  = Locbst - Iybest * 4;
      float Ybest  = float(Iybest + 1)/Iybin_Div;
      Ybest -= 4.55;
//
	getEndsLoop(Nmod, Ybest, NULL, 2);
//
      for (int i = 0; i < LAYERS; i++) {
//
//        int l = (k + 2*side)%4;
//        int m = (k < 2) ? k : 2;
//        int n = l/2;

          if (inc[i]) {
          Zw[i] = Xzaty[0][iStraws[i]][i][Nmod];
          Xw[i] = Xzaty[1][iStraws[i]][i][Nmod];
          Dr[i] = Hits[iStraws[i]][i][Nmod];
          straw.Ybot = (*strawGeoms[i/2])[Nmod][i][iStraws[i]][0][1];
          straw.Ytop = (*strawGeoms[i/2])[Nmod][i][iStraws[i]][1][1];
          coords[1][i/2][i] = straw.Ybot + ((straw.Ytop-straw.Ybot)*(Ybest -(-4.55))/(9.1));
          }
      }
//
//
      int iflagb = 0;
      straw_CommonTangent(Nmod,iStraws[2*side],2*side,iStraws[2*side+1],2*side+1,Itnbst,
                          coords[2][side][0],
                          coords[0][side][0],
                          coords[2][side][1],
                          coords[0][side][1],
                          iflagb);
      if(fabs(coords[2][side][1] - coords[2][side][0]) > 1.0e-6) {
                                 Grads[side] = (coords[0][side][1] - coords[0][side][0])/
                                               (coords[2][side][1] - coords[2][side][0]);
      }
      else {
                                  return;
      }
      Cints[side] = coords[0][side][0] - (Grads[side] * coords[2][side][0]);
      float wireb[2][4] = {0.0}; // [X,Z][1u,2u,3v,4v]
      float interb[2][4] = {0.0}; // [X,Z][1u,2u,3v,4v]
      float ptb[2][4] = {0.0}; // [X,Z][1u,2u,3v,4v]
      float Distb[4] = {0.0}; // [1u,2u,3v,4v]
      float resb[4] = {0.0}; // [1u,2u,3v,4v]

      for (int i = 0; i < 2; i++) {
          int j, k;
          if (side == 1)
              j = (i+2)%4;
          else
              j = i;
          k = j/2;
          if (inc[j]) {
              wireb[1][j] = Zw[j];
              wireb[0][j] = Xw[j];
              straw_ClosestPoint(Grads[k],Cints[k],
                                 wireb[1][j],wireb[0][j],
                                 interb[1][j],interb[0][j]);
              Distb[j] = sqrt(pow(wireb[1][j] - interb[1][j], 2)
                                 + pow(wireb[0][j] - interb[0][j], 2));
              if(Distb[j] > 1.2*Size[j]) {
                                        return;
              }
              straw_NearestTangent(
                                   coords[2][k][0], 
                                   coords[0][k][0], 
                                   coords[2][k][1], 
                                   coords[0][k][1], 
                                   iStraws[j], j+1,Nmod,
                                   ptb[1][j],
                                   ptb[0][j]);
          }
      }
// -------------------------------------------------------------------
// Find closest point to this line on the hit cells in layers 3 and 4
// We need to swap the extrapolated tangent into the v frame from the
// u frame:
// -------------------------------------------------------------------
//
      for (int i = 0; i < 2; i++) {
          coords[2][1-side][i] = coords[2][side][i];
          coordConvert(side, TO_XY,
                       coords[0][side][i],
                       coords[1][side][i],
                       coords[0][2][i],
                       coords[1][2][i]);
          coordConvert(1-side, TO_UV,
                       coords[0][2][i],
                       coords[1][2][i],
                       coords[0][1-side][i],
                       coords[1][1-side][i]);
      }
      if(fabs(coords[2][1-side][1] - coords[2][1-side][0]) > 1.0e-6) {
          Grads[1-side] = (coords[0][1-side][1] - coords[0][1-side][0])/
                          (coords[2][1-side][1] - coords[2][1-side][0]);
      }
      else {
          return;
      }
      Cints[1-side] = coords[0][1-side][0] - (Grads[1-side] * coords[2][1-side][0]);

      for (int i = 2; i < 4; i++) {
          int j, k;
          if (side == 1)
              j = (i+2)%4;
          else
              j = i;
          k = j/2;
          if (inc[j]) {
              wireb[1][j] = Zw[j];
              wireb[0][j] = Xw[j];
              straw_ClosestPoint(Grads[k],Cints[k],
                                 wireb[1][j],wireb[0][j],
                                 interb[1][j],interb[0][j]);
              Distb[j] = sqrt(pow(wireb[1][j] - interb[1][j], 2)
                                 + pow(wireb[0][j] - interb[0][j], 2));
              if(Distb[j] > 1.2*Size[j]) {
                                        return;
              }
              straw_NearestTangent(
                                   coords[2][k][0], 
                                   coords[0][k][0], 
                                   coords[2][k][1], 
                                   coords[0][k][1], 
                                   iStraws[j], j+1,Nmod,
                                   ptb[1][j],
                                   ptb[0][j]);
          }
      }
      for (int i = 2; i < 4; i++) {
          if (inc[i]) {
              resb[i] = sqrt( pow(ptb[1][i]-interb[1][i], 2) // Layer 1 (or 1-or-2)
                            + pow(ptb[0][i]-interb[0][i], 2));
              resb[i] *= 10000.0;
          }
      }
// ---------------------------------------------------------------------
// Convert all points back into xy frame for fitting:
// ---------------------------------------------------------------------
      float ptb_xy[2][4] = {0}; // [X,Z][1,2,3,4]
      float wireb_xy[2][4] = {0}; // [X,Z][1,2,3,4]
      float interb_xy[2][4] = {0}; // [X,Z][1,2,3,4]
      for (int i = 0; i < LAYERS; i++) {
          int j, k;
          if (side == 1)
              j = (i+2)%4;
          else
              j = i;
          k = j/2;
          if (inc[j]) {
              ptb_xy[1][j] = ptb[1][j];
              coordConvert(k, TO_XY, ptb[0][j],coords[1][k][j],ptb_xy[0][j],coords[1][2][j]);
              straw_FillArray(iStraws[j]+1,j+1,Nmod,ptb_xy[1][j],ptb_xy[0][j]);
              Yhits[iStraws[j]][j][Nmod] = Ybest;
              Mask[iStraws[j]][j][Nmod] = 1;

              wireb_xy[1][j] = wireb[1][j];
              coordConvert(k, TO_XY, wireb[0][j], coords[1][k][j], wireb_xy[0][j], coords[1][2][j]);
              interb_xy[1][j] = interb[1][j];
              coordConvert(k, TO_XY, interb[0][j], coords[1][k][j], interb_xy[0][j], coords[1][2][j]);
          }
      }
//
      StrHist.fYbest[0]->Fill(Ybest);

// ------------------------------------------------------------------
// Finally make the fitted TE from the obtained space points.
// Must set up xzaty in xy frame for this:
// ------------------------------------------------------------------
	getEndsLoop(Nmod, Ybest, iStraws, 3);
//
      int Ntofit = Nposte;
      for(int ifit = 0; ifit < Nposte; ifit++) {
         Xptsfit[ifit] = Xypost[1][ifit][Nmod];
         Zptsfit[ifit] = Xypost[0][ifit][Nmod];
      }
      float Grfit = 0.0;
      float Cfit  = 0.0;
      Linfit(Zptsfit,Xptsfit,Ntofit,Grfit,Cfit);
      float zte1 = 0.0;
      float xte1 = 0.0;
      float zte2 = 0.0;
      float xte2 = 0.0;
      straw_ClosestPoint(Grfit,Cfit,ptb_xy[1][0],ptb_xy[0][0],zte1,xte1);
      if (inc[3]) {
          straw_ClosestPoint(Grfit,Cfit,ptb_xy[1][3],ptb_xy[0][3],zte2,xte2);
      }
      else {
          straw_ClosestPoint(Grfit,Cfit,ptb_xy[1][2],ptb_xy[0][2],zte2,xte2);
      }
//
      float xdiff  = xte2 - xte1;
      float zdiff  = zte2 - zte1;
      float Phizx  = phi_arctan(xdiff,zdiff);
      float Phizxd = Phizx*180.0/Pibk;
      if(Phizxd > 180.0) {
                         Phizxd = Phizxd - 360.0;
      }
      StrHist.fPhizxd[0]->Fill(Phizxd);

// --------------------------------------------
// Turn into milliradians:
// --------------------------------------------
      float Phizxmr = Phizxd*1000.0*Pibk/180.0;
      StrHist.fZxte[0]->Fill(zte1,xte1);

//
      float Grrr = (xte2 - xte1) / (zte2 - zte1);
      float Ctrk = xte1 - (Grrr * zte1);
// ------------------------------------------------------------------
// Finally store residuals to fitted track.
// ------------------------------------------------------------------
      float intt[2][4] = {0};
      float ptt[2][4] = {0};
      for (int i = 0; i < LAYERS; i++) {
          int j, k;
          if (side == 1)
              j = (i+2)%4;
          else
              j = i;
          k = j/2;
          if (inc[j]) {
              straw_ClosestPoint(Grrr,Ctrk,wireb_xy[1][j],wireb_xy[0][j],intt[1][j],intt[0][j]);
              straw_NearestTangent(zte1,xte1,zte2,xte2,iStraws[j],j+1,Nmod,ptt[1][j],ptt[0][j]);
              if (options.verbose) {
                  cout<< " Layer" << j+1 << ": "<<iStraws[0]+1<<" "<<ptb[1][j]<<" "<<ptb[0][j]<<" "<<ptt[1][j]<<" "<<ptb[0][j]<<endl;
              }
          }
      }

// ------------------------------------------------------------------
// (Zpt1,Xpt1) etc are the space points on the drift circles that 
// have been fitted to make the track.
// (Zptt1,Xptt1) are the closest points on the drift circles to the
// fitted track...
// ------------------------------------------------------------------

      float rest[4]; // [1u,2u,3v,4v]
      for (int i = 0; i < LAYERS; i++) {
          int j, k;
          if (side == 1)
              j = (i+2)%4;
          else
              j = i;
          k = j/2;
          if (inc[j]) {
              rest[j] = sqrt( pow(ptt[1][j]-intt[1][j], 2)
                            + pow(ptt[0][j]-intt[0][j], 2));
              rest[j] *= 10000.0;
              float rrr = randomFloat();
              if(rrr > 0.5) {
                  rest[j] *= -1;
              }

              Zxhit[0][iStraws[i]][i][Nmod] = ptt[1][i];
              Zxhit[1][iStraws[i]][i][Nmod] = ptt[0][i];

          }
      }

//------------------------------------------------
//	Write to file for visualisation
//------------------------------------------------


      if (!first) {
      trackeventsfile << "	},\n";
      }
      else
	first = 0;
      trackeventsfile << "	{\n";
      trackeventsfile << "	module = "<<Nmod+1<<"\n";
      trackeventsfile << "	hitnum = 4\n";
      trackeventsfile << "	hits = (\n";

      for (int i = 0; i < LAYERS; i++) {
          trackeventsfile << "		{\n";
          trackeventsfile << "		X = "<<ptt[0][i]<<"\n";
          trackeventsfile << "		Z = "<<ptt[1][i]<<"\n";
          trackeventsfile << "		layer = "<<i+1<<"\n";
          trackeventsfile << "		straw = "<<iStraws[i]+1<<"\n";
	if (i == 3)
          trackeventsfile << "		}\n";
	else
          trackeventsfile << "		},\n";
      }

      trackeventsfile << "	);\n";

      trackeventsfile << "	line = {\n";

      trackeventsfile << "		Z1 = "<<zte1<<"\n";
      trackeventsfile << "		X1 = "<<xte1<<"\n";
      trackeventsfile << "		Z2 = "<<zte2<<"\n";
      trackeventsfile << "		X2 = "<<xte2<<"\n";
      trackeventsfile << "		}\n";

      trackeventsfile << "	Ybest = "<<Ybest<<"\n";

      if (options.verbose) {
          cout<<" Made a 4 hit Track Element in Module "<<Nmod+1<<endl;
      }
	n_HitsOut++;
      StrHist.fTrkMod[0]->Fill(Nmod+0.5);

}
float getSingleStraw(int Nmod, int ilay, int istr) {

                if(ilay < 2) {
                  straw.Xbot = strawGeom_u[Nmod][ilay][istr][0][0];
                  straw.Ybot = strawGeom_u[Nmod][ilay][istr][0][1];
                  straw.Zbot = strawGeom_u[Nmod][ilay][istr][0][2];
                  straw.Xtop = strawGeom_u[Nmod][ilay][istr][1][0];
                  straw.Ytop = strawGeom_u[Nmod][ilay][istr][1][1];
                  straw.Ztop = strawGeom_u[Nmod][ilay][istr][1][2];

                }
                if(ilay >= 2) {
                  straw.Xbot = strawGeom_v[Nmod][ilay][istr][0][0];
                  straw.Ybot = strawGeom_v[Nmod][ilay][istr][0][1];
                  straw.Zbot = strawGeom_v[Nmod][ilay][istr][0][2];
                  straw.Xtop = strawGeom_v[Nmod][ilay][istr][1][0];
                  straw.Ytop = strawGeom_v[Nmod][ilay][istr][1][1];
                  straw.Ztop = strawGeom_v[Nmod][ilay][istr][1][2];
                }
}
float getEnds(int Nmod, int ilay, int istr, int arg) {
	if (arg == 1) {
	
                strawEnds.Zendc = strawGeom[Nmod][ilay][istr][0][2];
                strawEnds.Xendc = strawGeom[Nmod][ilay][istr][0][0];
                strawEnds.Zenda = strawGeom[Nmod][ilay][istr][1][2];
                strawEnds.Xenda = strawGeom[Nmod][ilay][istr][1][0];
	}
	else if (arg == 2) {
		strawEnds.Zendcb = strawGeom[Nmod][ilay][istr][0][2];
		strawEnds.Xendcb = strawGeom[Nmod][ilay][istr][0][0];
		strawEnds.Zendab = strawGeom[Nmod][ilay][istr][1][2];
		strawEnds.Xendab = strawGeom[Nmod][ilay][istr][1][0];
	}
	else if (arg == 3) {
		strawEnds.Zzendc = strawGeom[Nmod][ilay][istr][0][2];
		strawEnds.Xxendc = strawGeom[Nmod][ilay][istr][0][0];
		strawEnds.Zzenda = strawGeom[Nmod][ilay][istr][1][2];
		strawEnds.Xxenda = strawGeom[Nmod][ilay][istr][1][0];
	}
}
int getEndsLoop(int Nmod, float Y, int iStraws[LAYERS], int arg) {

	if (arg == 1) {
          for(int ilay = 0; ilay < LAYERS; ilay++) {
             for(int istr = 0; istr < STRAWS; istr++) {
                getEnds(Nmod, ilay, istr, 1);
// ------------------------------------------------------------
// Work out wire position at the vertical position (y) at which
// the track finding will take place:
// ------------------------------------------------------------
// Old calculation in the xyz frame:
                float Zaty  = strawEnds.Zendc;
                float Xaty  = strawEnds.Xendc + ((strawEnds.Xenda-strawEnds.Xendc)*(Y-(-4.55))/(9.1));
// ------------------------------------
// New calculation in the u or v frame:
// ------------------------------------
 		getSingleStraw(Nmod, ilay, istr);
// ------------------------------------------------------------------
// So Xzaty contains the z/x wire positions in the u frame for
// Layers 1 and 2. And in the v frame for Layers 3 and 4.
// ------------------------------------------------------------------
                Xzaty[0][istr][ilay][Nmod] = straw.Ztop;
                Xzaty[1][istr][ilay][Nmod] = straw.Xtop;
             }
          }
	}

	if (arg == 2) {
          for(int ilay = 0; ilay < LAYERS; ilay++) {
             for(int istr = 0; istr < STRAWS; istr++) {
		getEnds(Nmod, ilay, istr, 2);
// ------------------------------------------------------------
// Work out wire position at the vertical position (y) at which
// the best track residuals were found:
// ------------------------------------------------------------
            float Zatyb  = strawEnds.Zendcb;
            float Xatyb  = strawEnds.Xendcb + ((strawEnds.Xendab-strawEnds.Xendcb)*(Y-(-4.55))/(9.1));
// ------------------------------------
// New calculation in the u or v frame:
// ------------------------------------   
 		getSingleStraw(Nmod, ilay, istr);
// ------------------------------------------------------------------
// So Xzaty contains the z/x wire positions in the u frame for
// Layers 1 and 2. And in the v frame for Layers 3 and 4.
// ------------------------------------------------------------------
            Xzaty[0][istr][ilay][Nmod] = straw.Ztop;
            Xzaty[1][istr][ilay][Nmod] = straw.Xtop;
         }
      }
		}

	if (arg == 3) {

	float Ybest = Y;
	for (int ilay = 0; ilay < LAYERS; ilay++) {
		getEnds(Nmod, ilay, iStraws[ilay], 3);
      float Zzaty  = strawEnds.Zzendc;
      float Xxaty  = strawEnds.Xxendc + ((strawEnds.Xxenda-strawEnds.Xxendc)*(Ybest -(-4.55))/(9.1));
      Xzaty[0][iStraws[ilay]][ilay][Nmod] = Zzaty;
      Xzaty[1][iStraws[ilay]][ilay][Nmod] = Xxaty;
	}

	}
}

int fill_TimeDist_arrays(float *TimeArray, float *DistArray, float Fac) {
   
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
}
