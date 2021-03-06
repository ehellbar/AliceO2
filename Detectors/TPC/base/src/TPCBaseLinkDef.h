// Copyright CERN and copyright holders of ALICE O2. This software is
// distributed under the terms of the GNU General Public License v3 (GPL
// Version 3), copied verbatim in the file "COPYING".
//
// See https://alice-o2.web.cern.ch/ for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

#ifdef __CLING__

#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;

#pragma link C++ class o2::TPC::GlobalPosition2D;
#pragma link C++ class o2::TPC::GlobalPosition3D;
#pragma link C++ class o2::TPC::LocalPosition2D;
#pragma link C++ class o2::TPC::LocalPosition3D;

#pragma link C++ class o2::TPC::CalArray<float>+;
#pragma link C++ class o2::TPC::CalArray<double>+;
#pragma link C++ class o2::TPC::CalArray<int>+;
#pragma link C++ class o2::TPC::CalArray<unsigned>+;
#pragma link C++ class o2::TPC::CalArray<short>+;
#pragma link C++ class o2::TPC::CalArray<bool>+;
#pragma link C++ class o2::TPC::CalDet<float>+;
#pragma link C++ class o2::TPC::CalDet<double>+;
#pragma link C++ class o2::TPC::CalDet<int>+;
#pragma link C++ class o2::TPC::CalDet<unsigned>+;
#pragma link C++ class o2::TPC::CalDet<short>+;
#pragma link C++ class o2::TPC::CalDet<bool>+;
#pragma link C++ class o2::TPC::ContainerFactory;
#pragma link C++ class o2::TPC::CRU;
#pragma link C++ class o2::TPC::Digit+;
#pragma link C++ class o2::TPC::DigitPos;
#pragma link C++ class o2::TPC::FECInfo;
#pragma link C++ class o2::TPC::Mapper;
#pragma link C++ class o2::TPC::PadCentre;
#pragma link C++ class o2::TPC::PadInfo;
#pragma link C++ class o2::TPC::PadPos;
#pragma link C++ class o2::TPC::PadRegionInfo;
#pragma link C++ class o2::TPC::PadROCPos;
#pragma link C++ class o2::TPC::PadSecPos;
#pragma link C++ class o2::TPC::PartitionInfo;
//#pragma link C++ class o2::TPC::RandomRing;
#pragma link C++ class o2::TPC::ROC;
#pragma link C++ class o2::TPC::Sector;

#pragma link C++ namespace o2::TPC::Painter;
#pragma link C++ function o2::TPC::Painter::Draw(CalArray<float>);
#pragma link C++ function o2::TPC::Painter::Draw(CalArray<double>);
#pragma link C++ function o2::TPC::Painter::Draw(CalArray<int>);
#pragma link C++ function o2::TPC::Painter::Draw(CalArray<short>);
#pragma link C++ function o2::TPC::Painter::Draw(CalArray<bool>);
#pragma link C++ function o2::TPC::Painter::Draw(CalDet<float>);
#pragma link C++ function o2::TPC::Painter::Draw(CalDet<double>);
#pragma link C++ function o2::TPC::Painter::Draw(CalDet<int>);
#pragma link C++ function o2::TPC::Painter::Draw(CalDet<short>);
#pragma link C++ function o2::TPC::Painter::Draw(CalDet<bool>);

#endif
