// Copyright 2019-2020 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

#include "../src/ControlServiceHelpers.h"
#include <catch_amalgamated.hpp>
#include <iostream>
#include <regex>
#include <string_view>

TEST_CASE("TestControlServiceHelper")
{
  using namespace o2::framework;
  std::match_results<std::string_view::const_iterator> match;
  REQUIRE(ControlServiceHelpers::parseControl("foo", match) == false);
  std::match_results<std::string_view::const_iterator> match2;
  std::string sv = "CONTROL_ACTION: READY_TO_QUIT_ME";
  REQUIRE(ControlServiceHelpers::parseControl(sv, match2) == true);
  REQUIRE(match2[1].str() == "QUIT");
  REQUIRE(match2[2].str() == "ME");
  std::string sv2 = "   dsca  CONTROL_ACTION: READY_TO_QUIT_ME";
  REQUIRE(ControlServiceHelpers::parseControl(sv2, match2) == true);
  REQUIRE(match2[1].str() == "QUIT");
  REQUIRE(match2[2].str() == "ME");
  const static std::regex controlRE2("^(NOTIFY_STREAMING_STATE) (IDLE|STREAMING|EOS)", std::regex::optimize);
  const static std::regex controlRE3("^(NOTIFY_DEVICE_STATE) ([A-Z ]*)", std::regex::optimize);
  std::string sv3 = "   dsca  CONTROL_ACTION: NOTIFY_STREAMING_STATE IDLE";
  REQUIRE(ControlServiceHelpers::parseControl(sv3, match2) == true);
  REQUIRE(match2[1].str() == "NOTIFY_STREAMING_STATE");
  REQUIRE(match2[2].str() == "IDLE");
  std::string_view sv4 = "   dsca  CONTROL_ACTION: NOTIFY_STREAMING_STATE IDLE";
  REQUIRE(ControlServiceHelpers::parseControl(sv4, match2) == true);
  REQUIRE(match2[1].str() == "NOTIFY_STREAMING_STATE");
  REQUIRE(match2[2].str() == "IDLE");
  std::string_view sv5 = "   asdsca  CONTROL_ACTION: NOTIFY_DEVICE_STATE STOP";
  REQUIRE(ControlServiceHelpers::parseControl(sv5, match2) == true);
  REQUIRE(match2[1].str() == "NOTIFY_DEVICE_STATE");
  REQUIRE(match2[2].str() == "STOP");
  std::string_view sv6 = "   asdsca  CONTROL_ACTION: ";
  REQUIRE(ControlServiceHelpers::parseControl(sv6, match2) == false);
}
