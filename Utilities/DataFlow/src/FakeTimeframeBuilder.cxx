// Copyright CERN and copyright holders of ALICE O2. This software is
// distributed under the terms of the GNU General Public License v3 (GPL
// Version 3), copied verbatim in the file "COPYING".
//
// See https://alice-o2.web.cern.ch/ for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

#include "DataFlow/FakeTimeframeBuilder.h"
#include "TimeFrame/TimeFrame.h"
#include "Headers/DataHeader.h"
#include <string>
#include <vector>
#include <functional>
#include <cstring>

using DataHeader = o2::Header::DataHeader;
using DataDescription = o2::Header::DataDescription;
using DataOrigin = o2::Header::DataOrigin;
using IndexElement = o2::DataFormat::IndexElement;

namespace {
  o2::Header::DataDescription lookupDataDescription(const char *key) {
    if (strcmp(key, "RAWDATA") == 0)
      return o2::Header::gDataDescriptionRawData;
    else if (strcmp(key, "CLUSTERS") == 0)
      return o2::Header::gDataDescriptionClusters;
    else if (strcmp(key, "TRACKS") == 0)
      return o2::Header::gDataDescriptionTracks;
    else if (strcmp(key, "CONFIG") == 0)
      return o2::Header::gDataDescriptionConfig;
    else if (strcmp(key, "INFO") == 0)
      return o2::Header::gDataDescriptionInfo;
    return o2::Header::gDataDescriptionInvalid;
  }

  o2::Header::DataOrigin lookupDataOrigin(const char *key) {
    if (strcmp(key, "TPC") == 0)
      return o2::Header::gDataOriginTPC;
    if (strcmp(key, "TRD") == 0)
      return o2::Header::gDataOriginTRD;
    if (strcmp(key, "TOF") == 0)
      return o2::Header::gDataOriginTOF;
    if (strcmp(key, "ITS") == 0)
      return o2::Header::gDataOriginITS;
    return o2::Header::gDataOriginInvalid;
  }


}

namespace o2 { namespace DataFlow {

std::unique_ptr<char[]> fakeTimeframeGenerator(std::vector<FakeTimeframeSpec> &specs, std::size_t &totalSize) {
  // Calculate the total size of your timeframe. This is
  // given by:
  // - N*The size of the data header (this should actually depend on the
  //   kind of data as different dataDescriptions will probably have
  //   different headers).
  // - Sum_N(The size of the buffer_i)
  // - The size of the index header
  // - N*sizeof(dataheader)
  // Assuming all the data header
  size_t sizeOfHeaders = specs.size()*sizeof(DataHeader);
  size_t sizeOfBuffers = 0;
  for (auto && spec : specs) {
    sizeOfBuffers += spec.bufferSize;
  }
  size_t sizeOfIndexHeader = sizeof(DataHeader);
  size_t sizeOfIndex = sizeof(IndexElement)*specs.size();
  totalSize = sizeOfHeaders + sizeOfBuffers + sizeOfIndexHeader + sizeOfIndex;

  // Add the actual - data
  auto buffer = std::make_unique<char[]>(totalSize);
  char *bi = buffer.get();
  std::vector<IndexElement> headers;
  int count = 0;
  for (auto &&spec : specs) {
    IndexElement el;
    el.first.dataDescription = lookupDataDescription(spec.dataDescription);
    el.first.dataOrigin = lookupDataOrigin(spec.origin);
    el.first.payloadSize = spec.bufferSize;
    el.first.headerSize = sizeof(el.first);
    el.second = count++;
    // Let's zero at least the header...
    memset(bi, 0, sizeof(el.first));
    memcpy(bi, &el, sizeof(el.first));
    headers.push_back(el);
    bi += sizeof(el.first);
    spec.bufferFiller(bi, spec.bufferSize);
    bi += spec.bufferSize;
  }

  // Add the index
  DataHeader index;
  index.dataDescription = DataDescription("TIMEFRAMEINDEX");
  index.dataOrigin = DataOrigin("FKE");
  index.headerSize = sizeOfIndexHeader;
  index.payloadSize = sizeOfIndex;
  memcpy(bi, &index, sizeof(index));
  memcpy(bi+sizeof(index), headers.data(), headers.size() * sizeof(IndexElement));
  return std::move(buffer);
}

}} // o2::Headers
