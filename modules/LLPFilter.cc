/*
 *  Delphes: a framework for fast simulation of a generic collider experiment
 *  Copyright (C) 2012-2014  Universite catholique de Louvain (UCL), Belgium
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

/** \class LLPFilter
 *
 *  Removes particles with specific PDG codes
 *
 *  \author M. Selvaggi
 *
 */

#include "modules/LLPFilter.h"

#include "classes/DelphesClasses.h"
#include "classes/DelphesFactory.h"
#include "classes/DelphesFormula.h"

#include "ExRootAnalysis/ExRootClassifier.h"
#include "ExRootAnalysis/ExRootFilter.h"
#include "ExRootAnalysis/ExRootResult.h"

#include "TDatabasePDG.h"
#include "TFormula.h"
#include "TLorentzVector.h"
#include "TMath.h"
#include "TObjArray.h"
#include "TRandom3.h"
#include "TString.h"

#include <algorithm>
#include <iostream>
#include <sstream>
#include <stdexcept>

using namespace std;

//------------------------------------------------------------------------------

LLPFilter::LLPFilter() :
  fItInputArray(0)
{
}

//------------------------------------------------------------------------------

LLPFilter::~LLPFilter()
{
}

//------------------------------------------------------------------------------

void LLPFilter::Init()
{

  ExRootConfParam param;
  Size_t i, size;

  // PT threshold
  fPTMin = GetDouble("PTMin", 0.0);

  fInvert = GetBool("Invert", false);
  fCscDecay = GetBool("CscDecay", true);

  // no pileup
  fRequireNotPileup = GetBool("RequireNotPileup", false);

  fRequireStatus = GetBool("RequireStatus", false);
  fStatus = GetInt("Status", 1);

  fRequireCharge = GetBool("RequireCharge", false);
  fCharge = GetInt("Charge", 1);

  // import input array
  fInputArray = ImportArray(GetString("InputArray", "Delphes/allParticles"));
  fItInputArray = fInputArray->MakeIterator();

  param = GetParam("PdgCode");
  size = param.GetSize();

  // read PdgCodes to be filtered out from the data card

  fPdgCodes.clear();
  for(i = 0; i < size; ++i)
  {
    fPdgCodes.push_back(param[i].GetInt());
  }

  // create output array
  fOutputArray = ExportArray(GetString("OutputArray", "filteredParticles"));
}

//------------------------------------------------------------------------------

void LLPFilter::Finish()
{
  if(fItInputArray) delete fItInputArray;
}

//------------------------------------------------------------------------------

void LLPFilter::Process()
{
  Candidate *candidate;
  Int_t pdgCode;
  Bool_t pass;
  Double_t pt, eta;

  fItInputArray->Reset();
  while((candidate = static_cast<Candidate *>(fItInputArray->Next())))
  {
    //all distance units are in mm
    pdgCode = candidate->PID;
    const TLorentzVector &candidateMomentum = candidate->Momentum;
    const TLorentzVector &candidateProdPosition = candidate->Position;
    const TLorentzVector &candidateDecayPosition = candidate->DecayPosition;
    pt = candidateMomentum.Pt();
    eta = candidateMomentum.Eta();
    if(pt < fPTMin) continue;
    if (candidate->D2-candidate->D1 != 2) continue;//require 3 daughters
    if(find(fPdgCodes.begin(), fPdgCodes.end(), pdgCode) == fPdgCodes.end()) continue; //require pdgID is one of the LLP id
    if (pdgCode != 9900012)cout<<"something is wrong"<<endl;
    if (fCscDecay)
    {
      if (abs(eta) < 2.4 && abs(eta) > 0.9
         && abs(candidateDecayPosition.Z())<11000 && abs(candidateDecayPosition.Z())>5680
         && sqrt(pow(candidateDecayPosition.X(),2)+pow(candidateDecayPosition.Y(),2)) < 6955)
      {
        fOutputArray->Add(candidate);
      }

    }
    else{
      fOutputArray->Add(candidate);
    }
  }
}
