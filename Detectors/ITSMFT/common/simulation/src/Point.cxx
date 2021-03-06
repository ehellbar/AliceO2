// Copyright CERN and copyright holders of ALICE O2. This software is
// distributed under the terms of the GNU General Public License v3 (GPL
// Version 3), copied verbatim in the file "COPYING".
//
// See https://alice-o2.web.cern.ch/ for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

/// \file Point.cxx
/// \brief Implementation of the Point class

#include "ITSMFTSimulation/Point.h"

#include <iostream>

ClassImp(o2::ITSMFT::Point)

using std::cout;
using std::endl;
using namespace o2::ITSMFT;
using namespace o2; //::Base;


void Point::Print(const Option_t *opt) const
{
  printf("Det: %5d Track: %6d E.loss: %.3e P: %+.3e %+.3e %+.3e\n"
	 "PosIn: %+.3e %+.3e %+.3e PosOut: %+.3e %+.3e %+.3e\n",
	 GetDetectorID(),GetTrackID(),GetEnergyLoss(),GetPx(),GetPy(),GetPz(),
	 GetStartX(),GetStartY(),GetStartZ(),GetX(),GetY(),GetZ() );
}


