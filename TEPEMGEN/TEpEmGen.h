#ifndef ROOT_TEpEmGen
#define ROOT_TEpEmGen
/* Copyright(c) 1998-2002, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */
//------------------------------------------------------------------------
// TEpEmGen is an interface class to fortran event generator of
// single e+e- pair production in ultraperipheral PbPb collisions
// at 5.5 GeV/c
//%
// Yuri.Kharlov@cern.ch
// 9 October 2002
//------------------------------------------------------------------------

#include "TGenerator.h"

#include <memory>

namespace o2::aegis::tepemgen { class EpEmSampler; }

// c++ interface to the f77 program - event generator of
// e+e- pair production in ultraperipheral ion collisions
// Author: Yuri Kharlov, 20 September 2002
//
// Revised on September 2018 for ALICEo2: Roberto Preghenella (preghenella@bo.infn.it)

class TEpEmGen : public TGenerator {

 public:
  /// Which implementation of the generator to use.
  ///
  /// kFortran is the historical epemgen.f/diffcross.f path; kCpp is the port
  /// in ../TEPEMGENCPP, which evaluates the cross section in long double and
  /// escalates to __float128 where cancellation demands it (see O2-6340 and
  /// TEPEMGENCPP/doc/02-findings.md). The two produce identical kinematics
  /// from the same seeded gRandom; the weights differ by the precision
  /// correction, at the 1e-3..1e-4 level on the total cross section.
  ///
  /// Both are kept while the port is being validated. The default can be set
  /// per process without touching code via TEPEMGEN_BACKEND=cpp|fortran.
  enum EBackend { kFortran = 0, kCpp = 1 };

  /// Default for objects constructed afterwards. Reads TEPEMGEN_BACKEND on
  /// first use.
  static EBackend DefaultBackend();
  static void SetDefaultBackend(EBackend b);

  void SetBackend(EBackend b) { fBackend = b; }
  EBackend GetBackend() const { return fBackend; }

  TEpEmGen();
  virtual ~TEpEmGen();
  
  void Initialize(Double_t ymin, Double_t ymax, Double_t ptmin, Double_t ptmaxm, Double_t cm_energy = 5160., Double_t Z = 82.);
  virtual void GenerateEvent() {TGenerator::GenerateEvent();};
  Int_t ImportParticles(TClonesArray *particles, Option_t *option);

 protected:
  void GenerateEvent (Double_t ymin, Double_t ymax, Double_t ptmin, Double_t ptmax,
	       	      Double_t &yElectron, Double_t &yPositron,
		      Double_t &xElectron, Double_t &xPositron,
		      Double_t &phi12,     Double_t &weight);
  Double_t GetXsection();
  Double_t GetDsection();

  EBackend fBackend;   //! selected implementation

  // Transient: not streamable, and meaningless to persist.
  std::unique_ptr<o2::aegis::tepemgen::EpEmSampler> fSampler;  //!

  ClassDef(TEpEmGen,1);  //Interface to EpEmGen Event Generator
};

#endif
