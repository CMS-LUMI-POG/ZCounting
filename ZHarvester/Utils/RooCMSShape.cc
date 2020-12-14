/*****************************************************************************
 * Project: CMS detector at the CERN
 *
 * Package: PhysicsTools/TagAndProbe/RooCMSShape
 *
 *
 * Authors:
 *   Nadia Adam, Princeton - neadam@princeton.edu
 *   Adam Hunt, Princeton  - ahunt@princeton.edu
 *   Kalanand Mishra, Fermilab - kalanand@fnal.gov
 *
 * Description:
 *   Defines a probability density function which has exponential decay
 *   distribution at high mass beyond the pole position (say, Z peak)
 *   but turns over (i.e., error function) at low mass due to threshold
 *   effect. We use this to model the background shape in Z->ll invariant
 *   mass.
 * History:
 *
 *
 *****************************************************************************/

#include "RooCMSShape.h"

ClassImp(RooCMSShape);

RooCMSShape::RooCMSShape(const char* name,
                         const char* title,
                         RooAbsReal& _x,
                         RooAbsReal& _alpha,
                         RooAbsReal& _beta,
                         RooAbsReal& _gamma,
                         RooAbsReal& _peak)
    : RooAbsPdf(name, title),
      x("x", "x", this, _x),
      alpha("alpha", "alpha", this, _alpha),
      beta("beta", "beta", this, _beta),
      gamma("gamma", "gamma", this, _gamma),
      peak("peak", "peak", this, _peak) {
          std::cout<<"make RooCMSShape 1"<<std::endl;

      }

RooCMSShape::RooCMSShape(const RooCMSShape& other, const char* name)
    : RooAbsPdf(other, name),
      x("x", this, other.x),
      alpha("alpha", this, other.alpha),
      beta("beta", this, other.beta),
      gamma("gamma", this, other.gamma),
      peak("peak", this, other.peak) {
          std::cout<<"make RooCMSShape 2"<<std::endl;

      }

Double_t RooCMSShape::evaluate() const {
  // ENTER EXPRESSION IN TERMS OF VARIABLE ARGUMENTS HERE
  std::cout<<"evaluate"<<std::endl;
  Double_t v = (alpha - x) * beta;
  // Double_t erf = RooMath::erfc((alpha - x) * beta);
  // Double_t erf = 1 - TMath::Erf((alpha - x) * beta);
  Double_t u = (x - peak) * gamma;

  if (u < -70)
    u = 1e20;
  else if (u > 70)
    u = 0;
  else
    u = exp(-u);  //exponential decay

  if(v > 5)
    v = 2.;
  else if(v < -5.)
    v = 0.;
  else
    v = 1 - TMath::Erf(v);

  return v * u;
}