// Copyright CERN and copyright holders of ALICE O2. This software is
// distributed under the terms of the GNU General Public License v3 (GPL
// Version 3), copied verbatim in the file "COPYING".
//
// See https://alice-o2.web.cern.ch/ for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

#ifndef AliceO2_TPC_DigitPos_H
#define AliceO2_TPC_DigitPos_H

#include "TPCBase/Defs.h"
#include "TPCBase/CRU.h"
#include "TPCBase/PadPos.h"
#include "TPCBase/PadSecPos.h"

namespace o2 {
namespace TPC {

class DigitPos {
public:
  DigitPos() {}
  DigitPos(CRU c, PadPos pad) : mCRU(c), mPadPos(pad) {}
  const CRU&  getCRU() const { return mCRU; }
  CRU& cru()       { return mCRU; }
  PadPos    getPadPos()       const { return mPadPos; }
  PadPos    getGlobalPadPos() const;
  PadSecPos getPadSecPos()    const;

  PadPos& padPos()       { return mPadPos; }

  bool isValid() const { return mPadPos.isValid(); }

  bool    operator==(const DigitPos& other)  const { return (mCRU==other.mCRU) && (mPadPos==other.mPadPos); }
  bool    operator!=(const DigitPos& other)  const { return (mCRU!=other.mCRU) || (mPadPos!=other.mPadPos); }
  bool    operator< (const DigitPos& other)  const { return (mCRU <other.mCRU) && (mPadPos <other.mPadPos); }

private:
  CRU mCRU{};
  PadPos mPadPos{};          /// Pad position in the local partition coordinates: row starts from 0 for each partition
};

}
}
#endif
