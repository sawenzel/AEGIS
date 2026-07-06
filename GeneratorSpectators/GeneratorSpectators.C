#include "FairGenerator.h"
#include "GeneratorSpectators.h"
#include "Generators/GeneratorTGenerator.h"
#include <fairlogger/Logger.h>

class GeneratorSpectatorsO2 : public o2::eventgen::GeneratorTGenerator {

public:
  GeneratorSpectatorsO2() : GeneratorTGenerator("GeneratorSpectators", "GeneratorSpectators")
  {
    spec = new GeneratorSpectators();
    setTGenerator(spec);
  };

  ~GeneratorSpectatorsO2() { delete spec; };

  Bool_t Init() override
  {
    spec->Init();
    return true;
  }

  void updateHeader(o2::dataformats::MCEventHeader *eventHeader) final
  {
    // not updating the header, just printing particles for the time being
    LOG(info) << Form("--- Number of particles: %zu ---", mParticles.size());
    for (int i = 0; i < mParticles.size(); ++i) {
      auto &p = mParticles[i];
      LOG(info) << Form(
          "PdgCode: %d, FirstMother: %d, FirstDaughter: %d, LastDaughter: %d",
          p.GetPdgCode(), p.GetFirstMother(), p.GetFirstDaughter(),
          p.GetLastDaughter());
    }
  }

  GeneratorSpectators* getGenerator() const { return spec; }

private:
  GeneratorSpectators *spec = nullptr;
};


FairGenerator *GeneratorSingleNeutron(float mom = 2680.,
                                      float beamDiv = 0.000032,
                                      float beamDivMin = -999.,
                                      float beamDivMax = -999.,
                                      float beamCrossAngle = 0.,
                                      float beamCrossAngleMin = -999.,
                                      float beamCrossAngleMax = -999.,
                                      int beamCrossPlane = 2)
{
  auto wrap = new GeneratorSpectatorsO2();
  auto spec = wrap->getGenerator();
  spec->SetNpart(1);
  spec->SetParticle(2112);
  spec->SetMomentum(mom);
  spec->SetDirection(0, 0., 0., 1.);
  spec->SetDivergence(beamDiv);
  spec->SetSampleDivergence(beamDivMin, beamDivMax);
  spec->SetCrossing(beamCrossAngle, beamCrossPlane);
  spec->SetSampleCrossing(beamCrossAngleMin, beamCrossAngleMax);
  return wrap;
}

FairGenerator *GeneratorSingleProton(float mom = 2680.,
                                     float beamDiv = 0.000032,
                                     float beamDivMin = -999.,
                                     float beamDivMax = -999.,
                                     float beamCrossAngle = 0.,
                                     float beamCrossAngleMin = -999.,
                                     float beamCrossAngleMax = -999.,
                                     int beamCrossPlane = 2)
{
  auto wrap = new GeneratorSpectatorsO2();
  auto spec = wrap->getGenerator();
  spec->SetNpart(1);
  spec->SetParticle(2212);
  spec->SetMomentum(mom);
  spec->SetDirection(0, 0., 0., 1.);
  spec->SetDivergence(beamDiv);
  spec->SetSampleDivergence(beamDivMin, beamDivMax);
  spec->SetCrossing(beamCrossAngle, beamCrossPlane);
  spec->SetSampleCrossing(beamCrossAngleMin, beamCrossAngleMax);
  return wrap;
}

FairGenerator *GeneratorNeutrons(int nNeutrons = -1,
                                 float b = -1.,
                                 bool useFluctuation = false,
                                 float mom = 2680.,
                                 float beamDiv = 0.000032,
                                 float beamDivMin = -999.,
                                 float beamDivMax = -999.,
                                 float beamCrossAngle = 0.,
                                 float beamCrossAngleMin = -999.,
                                 float beamCrossAngleMax = -999.,
                                 int beamCrossPlane = 2)
{
  auto wrap = new GeneratorSpectatorsO2();
  auto spec = wrap->getGenerator();
  spec->SetParticle(2112);
  spec->SetImpactParameter(b);
  spec->SetNpartFluctuation(useFluctuation);
  spec->SetNpart(nNeutrons);
  spec->SetMomentum(mom);
  spec->SetDirection(0, 0., 0., 1.);
  spec->SetDivergence(beamDiv);
  spec->SetSampleDivergence(beamDivMin, beamDivMax);
  spec->SetCrossing(beamCrossAngle, beamCrossPlane);
  spec->SetSampleCrossing(beamCrossAngleMin, beamCrossAngleMax);
  return wrap;
}

FairGenerator *GeneratorProtons(int nProtons = -1,
                                float b = -1.,
                                bool useFluctuation = false,
                                float mom = 2680.,
                                float beamDiv = 0.000032,
                                float beamDivMin = -999.,
                                float beamDivMax = -999.,
                                float beamCrossAngle = 0.,
                                float beamCrossAngleMin = -999.,
                                float beamCrossAngleMax = -999.,
                                int beamCrossPlane = 2)
{
  auto wrap = new GeneratorSpectatorsO2();
  auto spec = wrap->getGenerator();
  spec->SetParticle(2212);
  spec->SetImpactParameter(b);
  // TODO: implement parameterizations for the number of protons produced
  // At the moment, those for the neutrons are used
  spec->SetNpartFluctuation(useFluctuation);
  spec->SetNpart(nProtons);
  spec->SetMomentum(mom);
  spec->SetDirection(0, 0., 0., 1.);
  spec->SetDivergence(beamDiv);
  spec->SetSampleDivergence(beamDivMin, beamDivMax);
  spec->SetCrossing(beamCrossAngle, beamCrossPlane);
  spec->SetSampleCrossing(beamCrossAngleMin, beamCrossAngleMax);
  return wrap;
}

