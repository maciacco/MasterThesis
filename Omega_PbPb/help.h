struct MiniLambda {
  Double32_t pt;
  Double32_t eta;
  Double32_t mass;
  Double32_t ct;
  Double32_t radius;      //[0,101.6,8]
  Double32_t dcaV0PV;     //[0,10.16,16]
  Double32_t dcaPiPV;     //[0,20.32,8]
  Double32_t dcaPrPV;     //[0,20.32,8]
  Double32_t dcaV0tracks; //[0,2.54,8]
  Double32_t cosPA;       //[0.95,1,16]
  Double32_t tpcNsigmaPi; //[-5,5,8]
  Double32_t tpcNsigmaPr; //[-5,5,8]
  unsigned char tpcClV0Pi;
  unsigned char tpcClV0Pr;
  unsigned char centrality;
  bool matter;
  bool hasTOFhit;
  bool hasITSrefit;
};

struct MiniLambdaMC : public MiniLambda {
  float ptMC;
  float etaMC;
  float ctMC;
  float yMC;
  float ptMotherMC;
  int pdg;
  bool isPrimary;
  bool isReconstructed;
  unsigned char flag;
};

struct MiniCascade {
  Double32_t pt;
  Double32_t eta;
  Double32_t mass;
  Double32_t ct;
  Double32_t radius; //[0,25.4,8]
  Double32_t radiusV0; //[0,25.4,8]
  Double32_t dcaBachPV; //[0,2.54,8]
  Double32_t dcaV0PV; //[0,2.54,8]
  Double32_t dcaV0piPV; //[0,2.54,8]
  Double32_t dcaV0prPV; //[0,2.54,8]
  Double32_t dcaV0tracks; //[0,2.54,8]
  Double32_t dcaBachV0; //[0,2.54,8]
  Double32_t cosPA; //[0.95,1,16]
  Double32_t cosPAV0; //[0.95,1,16]
  Double32_t V0invMassDelta; //[-0.01,0.01,8]
  Double32_t tpcNsigmaBach; //[-5,5,8]
  Double32_t tpcNsigmaV0Pr; //[-5,5,8]
  Double32_t tpcNsigmaV0Pi; //[-5,5,8]
  Double32_t competingMass; //[0,0.254,8]
  Double32_t bachBarCosPA;  //[0.9999, 1., 8]
  unsigned char tpcClBach;
  unsigned char tpcClV0Pr;
  unsigned char tpcClV0Pi;
  unsigned char centrality;
  bool matter;
  bool hasTOFhit;
  bool hasITSrefit;
  bool isOmega;
};

struct MiniCascadeMC : public MiniCascade {
  float ptMC;
  float etaMC;
  float ctMC;
  float yMC;
  int pdg;
  bool isReconstructed;
  unsigned char flag;
};

struct MiniLambdaBDTOut {
  Double32_t pt;
  Double32_t ct;
  Double32_t mass;
  Double32_t bdtOutputPrompt;
  Double32_t bdtOutputNonPrompt;
  Double32_t bdtOutputBackground;
  bool matter;
  bool pileUpCheck;
};