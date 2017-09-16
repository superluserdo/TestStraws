#define PI 3.1415926535 
#define MODULES 3
#define LAYERS 4
#define STRAWS 32

// Structures
struct num_struct {
	int modules;
	int layers;
	int straws;
};

 struct options_struct {
	int help;
	int file;
	int verbose;
	int output;
};

extern const struct num_struct num;

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
struct singleStraw {
	float Xbot;
	float Xtop;
	float Ybot;
	float Ytop;
	float Zbot;
	float Ztop;
};
struct strawEndsStruct {
	float Zendc;
	float Xendc;
	float Zenda;
	float Xenda;
	float Zendcb;
	float Xendcb;
	float Zendab;
	float Xendab;
	float Zzendc;
	float Xxendc;
	float Zzenda;
	float Xxenda;
};
// Histograms
void DefineColours();
void BookHistograms();
void BookStrHistograms();
void PlotHistograms();
void PlotStrHistograms(int Ipage);
void WriteHistograms();
void WriteStrHistograms(int Ipage);
void WriteStraws();

// Initialise
void geom_Init();
void zero_Init(int imode);
float randomFloat();

void Vfill(float* Array, int Len, float Value);
void Vzero(float* Array, int Len);
int fill_TimeDist_arrays(float *TimeArray, float *DistArray, float Fac);

// Functions on straws
void straw_NearestTangent(float zz1,float xx1,float zz2,float xx2,int Nstraw,int Layer,int Nmod,float &Zpt, float &Xpt);
void straw_LocalAngle(float zz1,float xx1,float zz2,float xx2,int Nst,int Nlay,int Nmod,float &Angle, int &Iflag);
void straw_CommonTangent(int Nmod,int Nst1,int Layer1,int Nst2,int Layer2,int Icand,
float &zz1,float &xx1,float &zz2,float &xx2,int &Iflag);
void straw_FillArray(int Nstraw, int Layer, int Nmod, float zpt, float xpt);
void straw_ClosestPoint(float Gr,float c1,float Xpt,float Ypt,float &Xinter,float &Yinter);
void straw_TimeDist_1(float Time,float &Dist);
void straw_TimeDist_2(float Time,float &Dist);
void straw_TimeDist_3(float Time,float &Dist);

// Track finding
void straw_Evt(int &Goodt0);
void straw_Analyse(ofstream &trackeventsfile);
void straw_FindEvent(int Nmod, ofstream &trackeventsfilee);

// Maths & fitting
void Linfit(float* Xpts,float* Ypts,int Npts,float &Grad,float &Cint);

int  round(float var);
float phi_arctan(float py , float px);

double FitBack(double* x, double* par);
double Lorentzian(double* x, double* par);
double Gaussian1(double* x, double* par);
double Gaussian2(double* x, double* par);
double FitFunction(double* x, double* par);
double FitFunction2(double* x, double* par);

// Convenience
float getSingleStraw(int Nmod, int ilay, int istr);
float getEnds(int Nmod, int ilay, int istr, int arg);
int getEndsLoop(int Nmod, float Y, int iStraws[LAYERS], int arg);

// Coordinate conversion
#define TYPE_U 0
#define TYPE_V 1
#define TO_XY 0
#define TO_UV 1
void coordConvert(int uvType, int direction, float xIn , float yIn , float &xOut , float &yOut);
