// Copyright CERN and copyright holders of ALICE O2. This software is
// distributed under the terms of the GNU General Public License v3 (GPL
// Version 3), copied verbatim in the file "COPYING".
//
// See https://alice-o2.web.cern.ch/ for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

#ifndef ALICEO2_TPC_CALIBRAWBASE_H_
#define ALICEO2_TPC_CALIBRAWBASE_H_

/// \file   CalibRawBase.h
/// \author Jens Wiechula, Jens.Wiechula@ikf.uni-frankfurt.de

#include <vector>
#include <memory>

#include "TString.h"
#include "Rtypes.h"

#include "TPCBase/Defs.h"
#include "TPCBase/Mapper.h"

#include "TPCReconstruction/GBTFrameContainer.h"


namespace o2
{
namespace TPC
{


/// \brief Base class for raw data calibrations
///
/// This class is the base class for raw data calibrations
/// It implements base raw reader functionality and calls 
/// an 'Update' function for each digit
///
/// origin: TPC
/// \author Jens Wiechula, Jens.Wiechula@ikf.uni-frankfurt.de

class CalibRawBase
{
  public:
    CalibRawBase(PadSubset padSubset = PadSubset::ROC) : mMapper(Mapper::instance()), mNevents(0), mTimeBinsPerCall(500), mPadSubset(padSubset) {;}

    /// Update function called once per digit
    ///
    /// \param sector
    virtual Int_t Update(const Int_t roc, const Int_t row, const Int_t pad,
                         const Int_t timeBin, const Float_t signal) = 0;

    /// add GBT frame container to process
    void addGBTFrameContainer(GBTFrameContainer *cont) { mGBTFrameContainers.push_back(std::unique_ptr<GBTFrameContainer>(cont)); }

    /// set number of time bins to process in one call to ProcessEvent
    void setTimeBinsPerCall(Int_t nTimeBins) { mTimeBinsPerCall = nTimeBins; }

    /// return the number of time bins processed in one call to ProcessEvent
    Int_t getTimeBinsPerCall() const { return mTimeBinsPerCall; }

    /// return pad subset type used
    PadSubset getPadSubset() const { return mPadSubset; }

    void setupContainers(TString fileInfo);
    
    /// Process one event with mTimeBinsPerCall length
    Bool_t ProcessEvent(); 
    
    /// Rewind the events
    void RewindEvents();

    /// Dump the relevant data to file
    virtual void dumpToFile(TString filename) {}
    
  protected:
    const Mapper&  mMapper;    //!< TPC mapper

  private:
    size_t    mNevents;                //!< number of processed events
    Int_t     mTimeBinsPerCall;        //!< numver of time bins to process in ProcessEvent
    PadSubset mPadSubset;              //!< pad subset type used
    std::vector<std::unique_ptr<GBTFrameContainer>> mGBTFrameContainers; //! raw reader pointer

    virtual void ResetEvent() = 0;
    virtual void EndEvent() {++mNevents; }
};

//----------------------------------------------------------------
// Inline Functions
//----------------------------------------------------------------
inline Bool_t CalibRawBase::ProcessEvent()
{
  if (!mGBTFrameContainers.size()) return kFALSE;
  ResetEvent();
  
  // loop over raw readers, fill digits for 500 time bins and process
  // digits
  
  const int nRowIROC = mMapper.getNumberOfRowsROC(0);

  std::vector<DigitData> digits(80);
  for (auto& reader_ptr : mGBTFrameContainers) {
    auto reader = reader_ptr.get();
    for (int i=0; i<mTimeBinsPerCall; ++i) {
      if (reader->getData(digits)) {
        for (auto& digi : digits) {
          CRU cru(digi.getCRU());
          const int sector = cru.sector().getSector();
          // TODO: OROC case needs subtraction of number of pad rows in IROC
          const PadRegionInfo& regionInfo = mMapper.getPadRegionInfo(cru.region());
          const PartitionInfo& partInfo = mMapper.getPartitionInfo(cru.partition());

          // row is local in region (CRU)
          const int row    = digi.getRow();
          const int pad    = digi.getPad();
          if (row==255 || pad==255) continue;

          int rowOffset = 0;
          switch (mPadSubset) {
            case PadSubset::ROC: {
                rowOffset = regionInfo.getGlobalRowOffset();
                rowOffset -= (cru.rocType()==RocType::OROC)*nRowIROC;
              break;
            }
            case PadSubset::Region: {
              break;
            }
            case PadSubset::Partition: {
              rowOffset = regionInfo.getGlobalRowOffset();
              rowOffset -= partInfo.getGlobalRowOffset();
              break;
            }
          }

          // modify row depending on the calibration type used
          const int timeBin= i; //digi.getTimeStamp();
          const float signal = digi.getChargeFloat();
          //const FECInfo& fecInfo = mTPCmapper.getFECInfo(PadSecPos(sector, row, pad));
          //printf("Call update: %d, %d, %d, %d (%d), %.3f -- reg: %02d -- FEC: %02d, Chip: %02d, Chn: %02d\n", sector, row, pad, timeBin, i, signal, cru.region(), fecInfo.getIndex(), fecInfo.getSampaChip(), fecInfo.getSampaChannel());
          Update(sector, row+rowOffset, pad, timeBin, signal );
        }
      }
      else {
        return kFALSE;
      }
      digits.clear();
    }
  }

  EndEvent();
  return kTRUE;
}

} // namespace TPC

} // namespace o2
#endif
