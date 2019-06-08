#ifndef XRAYLINE_HH
#define XRAYLINE_HH 1

#include "TString.h"

#include <vector>
#include <stdlib.h>

#include "mathfunc.h"

struct xrayline
{
public:
  TString element;
  TString linetype;
  std::vector<double> energies;
  std::vector<double> energies_err;
  std::vector<double> lorentzian_widths;
  std::vector<double> lorentzian_widths_err;
  std::vector<double> normalized_intensities;
  
  Double_t pdf(Double_t x) {
    Double_t result=0;
    for(Int_t i=0; i != int(energies.size()); ++i) {
      result += lorentzian(&x,normalized_intensities[i],energies[i],lorentzian_widths[i]);
    }
    return result;
  }
  
  // Tatsuno's LE tail function
  Double_t operator() (Double_t *x, Double_t *par) {
    Double_t e2c=par[0];// eV/ch
    Double_t ampl=par[1];// amplitude
    Double_t mean=par[2];// energy
    Double_t gaussian_sigma=par[3];// energy resolution
    Double_t letail_ratio=par[4];// tail ratio (area)
    Double_t letail_beta=par[5]/e2c;// tail beta [eV]
    Double_t result=0;
    for(Int_t i=0; i != int(energies.size()); ++i) {
      Double_t mi = (energies[i]-energies[0])/e2c+mean;
      Double_t wi = lorentzian_widths[i]/e2c;
      Double_t ii = ampl*normalized_intensities[i]/normalized_intensities[0];
      result += voigt(x,ii*(1.-letail_ratio),mi,gaussian_sigma,wi);
      result += letail(x,letail_ratio*ii,mi,gaussian_sigma,letail_beta);
    }
    return result;
  }
  
  void SetElement(TString elm) {element=elm;}
  void SetLineType(TString ltyp) {linetype=ltyp;}
  TString GetElement() {return element;}
  TString GetLineType() {return linetype;}
  
  void Print()
  {
    cout << "-------------------------------------------------------" << endl;
    cout << Form("element:  %s",element.Data()) << endl;
    cout << Form("lineType: %s",linetype.Data()) << endl;
    for(Int_t i=0; i != int(energies.size()); ++i) {
      cout << Form("%d: energy %9.3f+-%5.3f (eV), width %6.3f+-%5.3f (eV), norm %5.3f ",i,energies[i],energies_err[i],lorentzian_widths[i],lorentzian_widths_err[i],normalized_intensities[i]) << endl;
    }
    cout << "-------------------------------------------------------" << endl;
  }
};

// all objects are created when you read.
// Ti: Ref. PRA 73 (2006) 012508 & PRL 91 (2003) 240801
// Ni: Ref. PRA 56 (1997) 4554
xrayline TiKa={"Ti",
               "KAlpha",
               {4510.926,4509.467,4507.735,4513.848,4504.914,4502.611},
               {0.014,0.141,0.217,0.109,0.020,0.566},
               {1.32,1.54,2.77,1.75,1.73,3.30},
               {0.11,0.47,0.93,0.28,0.16,1.06},
               {0.512,0.073,0.067,0.032,0.292,0.024}
};
xrayline TiKb={"Ti",
               "KBeta",
               {4925.37,4930.096,4931.967,4935.59},
               {0.50,0.075,0.016,0.16},
               {16.3,4.25,0.42,0.47},
               {1.0,0.19,0.22,0.44},
               {0.199,0.455,0.326,0.0192}
};
xrayline VKa={"V",
              "KAlpha",
              {4952.237,4950.656,4948.266,4955.269,4944.672,4943.014},
              {0.012,0.184,0.261,0.141,0.021,0.303},
              {1.45,2.00,1.81,1.76,2.94,3.09},
              {0.02,0.03,0.70,0.30,0.04,0.26},
              {0.546,0.114,0.032,0.020,0.274,0.013}
};
xrayline VKb={"V",
              "KBeta",
              {5418.19, 5424.50, 5426.992},
              {0.1,0.1,0.1},
              {18.86, 5.48, 2.499},
              {0.1,0.1,0.1},
              {0.25774226, 0.23576424, 0.50649351}
};
xrayline CrKa={"Cr",
               "KAlpha",
               {5414.874,5414.099,5412.745,5410.583,5418.304,5405.551,5403.986},
               {0.002,0.006,0.016,0.030,0.038,0.003,0.003},
               {1.457,1.760,3.138,5.149,1.988,2.224,4.740},
               {0.002,0.007,0.020,0.051,0.058,0.004,0.050},
               {0.378,0.132,0.084,0.073,0.009,0.271,0.054}
};
xrayline CrKb={"Cr",
               "KBeta",
               {5947.00,5935.31,5946.24,5942.04,5944.93},
               {0.01,0.09,0.01,0.05,0.03},
               {1.70,15.98,1.90,6.69,3.37},
               {0.01,0.21,0.02,0.08,0.04},
               {0.307,0.236,0.172,0.148,0.137}
};
xrayline MnKa={"Mn",
               "KAlpha",
               {5898.853,5897.867,5894.829,5896.532,5899.417,5887.743,5886.495},
               {0.002,0.006,0.022,0.015,0.017,0.003,0.013},
               {1.715,2.043,4.499,2.663,0.969,2.361,4.216},
               {0.002,0.007,0.033,0.020,0.020,0.005,0.021},
               {0.353,0.141,0.079,0.066,0.005,0.229,0.110}
};
xrayline MnKb={"Mn",
               "KBeta",
               {6490.89,6486.31,6477.73,6490.06,6488.83},
               {0.01,0.07,0.08,0.02,0.03},
               {1.83,9.40,13.22,1.81,2.81},
               {0.01,0.09,0.15,0.02,0.04},
               {0.2547,0.2345,0.233,0.1645,0.1132}
};
xrayline FeKa={"Fe",
               "KAlpha",
               {6404.148,6403.295,6400.653,6402.077,6391.190,6389.106,6390.275},
               {0.002,0.005,0.019,0.012,0.004,0.022,0.033},
               {1.613,1.965,4.833,2.803,2.487,2.339,4.433},
               {0.003,0.005,0.026,0.015,0.006,0.019,0.012},
               {0.278,0.182,0.106,0.094,0.207,0.066,0.065}
};
xrayline FeKb={"Fe",
               "KBeta",
               {7046.90,7057.21,7058.36,7054.75},
               {0.09,0.02,0.01,0.06},
               {14.17,3.12,1.97,6.38},
               {0.17,0.02,0.02,0.08},
               {0.301,0.279,0.241,0.179}
};
xrayline CoKa={"Co",
               "KAlpha",
               {6930.425,6929.388,6927.676,6930.941,6915.713,6914.659,6913.078},
               {0.002,0.006,0.012,0.027,0.003,0.007,0.023},
               {1.795,2.695,4.555,0.808,2.406,2.773,4.463},
               {0.002,0.008,0.015,0.027,0.004,0.010,0.035},
               {0.37858,0.144,0.127,0.0088,0.197,0.095,0.050}
};
xrayline CoKb={"Co",
               "KBeta",
               {7649.60,7647.83,7639.87,7645.49,7636.21,7654.13},
               {0.08,0.03,0.08,0.06,0.20,0.11},
               {3.05,3.58,9.78,4.89,13.59,3.79},
               {0.01,0.03,0.12,0.08,0.37,0.16},
               {0.449,0.189,0.153,0.103,0.082,0.025}
};
xrayline NiKa={"Ni",
               "KAlpha",
               {7478.281,7476.529,7461.131,7459.874,7458.029},
               {0.002,0.013,0.004,0.015,0.048},
               {2.013,4.711,2.674,3.039,4.476},
               {0.002,0.015,0.005,0.021,0.078},
               {0.487,0.171,0.250,0.064,0.028}
};
xrayline NiKb={"Ni",
               "KBeta",
               {8265.01,8263.01,8256.67,8268.70},
               {0.01,0.02,0.10,0.07},
               {3.76,4.34,13.70,5.18},
               {0.02,0.03,0.16,0.09},
               {0.450,0.258,0.203,0.089}
};
xrayline CuKa={"Cu",
               "KAlpha",
               {8047.837,8045.367,8027.993,8026.504},
               {0.002,0.022,0.005,0.014},
               {2.285,3.358,2.666,3.571},
               {0.003,0.027,0.007,0.014},
               {0.579,0.080,0.236,0.105}
};
xrayline CuKb={"Cu",
               "KBeta",
               {8905.532,8903.109,8908.462,8897.387,8911.393},
               {0.002,0.010,0.020,0.050,0.057},
               {3.52,3.52,3.55,8.08,5.31},
               {0.01,0.01,0.03,0.08,0.08},
               {0.485,0.248,0.110,0.100,0.055}
};




xrayline KHe4La={"KHe4",
                 "LAlpha",
                 {6463.46},// from T.Koike
                 {0.15},// kokodeha tekito
                 {1.2},// chiral & pheno
                 {0.1},// tekito
                 {1.0}
};

xrayline KHe3La={"KHe3",
                 "LAlpha",
                 {6224.57},// from Yamagata-Hirenzaki
                 {0.15},// kokodeha tekito
                 {2.0},// pheno
                 {0.1},// tekito
                 {1.0}
};


inline xrayline get_xl(TString line_name)
{
  // this is not smart...
  Bool_t b_ka = (line_name.Contains("KAlpha")) || (line_name.Contains("Ka")) || (line_name.Contains("ka"));
  Bool_t b_kb = (line_name.Contains("KBeta")) || (line_name.Contains("Kb")) || (line_name.Contains("kb"));
  Bool_t b_la = (line_name.Contains("LAlpha")) || (line_name.Contains("La")) || (line_name.Contains("la"));
  if ((b_ka==false) && (b_kb==false) && (b_la==false)){
    cout << "Error: cannot find linetype for " << line_name << endl;
    exit (EXIT_FAILURE);
    xrayline hoge{};
    return hoge;
  }
  if ((line_name.Contains("Ti")) || (line_name.Contains("ti"))) {
    if (b_ka)      {return TiKa;}
    else if (b_kb) {return TiKb;}
  }else if ((line_name.Contains("V")) || (line_name.Contains("v"))) {
    if (b_ka)      {return VKa;}
    else if (b_kb) {return VKb;}
  }else if ((line_name.Contains("Cr")) || (line_name.Contains("cr"))) {
    if (b_ka)      {return CrKa;}
    else if (b_kb) {return CrKb;}
  }else if ((line_name.Contains("Mn")) || (line_name.Contains("mn"))) {
    if (b_ka)      {return MnKa;}
    else if (b_kb) {return MnKb;}
  }else if ((line_name.Contains("Fe")) || (line_name.Contains("fe"))) {
    if (b_ka)      {return FeKa;}
    else if (b_kb) {return FeKb;}
  }else if ((line_name.Contains("Co")) || (line_name.Contains("co"))) {
    if (b_ka)      {return CoKa;}
    else if (b_kb) {return CoKb;}
  }else if ((line_name.Contains("Ni")) || (line_name.Contains("ni"))) {
    if (b_ka)      {return NiKa;}
    else if (b_kb) {return NiKb;}
  }else if ((line_name.Contains("Cu")) || (line_name.Contains("cu"))) {
    if (b_ka)      {return CuKa;}
    else if (b_kb) {return CuKb;}
  }else if ((line_name.Contains("KHe4")) || (line_name.Contains("khe4"))) {
    if (b_la)      {return KHe4La;}
  }else if ((line_name.Contains("KHe3")) || (line_name.Contains("khe3"))) {
    if (b_la)      {return KHe3La;}
  }
  cout << "Error: cannot find xrayline for " << line_name << endl;
  exit (EXIT_FAILURE);
  xrayline hoge{};
  return hoge;
}


#endif
