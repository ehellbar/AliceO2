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
#ifndef O2_FRAMEWORK_OBJECTCACHE_H_
#define O2_FRAMEWORK_OBJECTCACHE_H_

#include "Framework/DataRef.h"
#include <unordered_map>
#include <map>
#include <string>

namespace o2::framework
{

/// A cache for CCDB objects or objects in general
/// which have more than one timeframe of lifetime.
///
/// The cache is keyed *per path* rather than by a global id-derived hash.
/// Earlier versions stored a `matcherToId` (path -> id) map alongside an
/// `idToObject` (id -> deserialised object) map keyed by the SHM payload
/// pointer of the incoming message. Because SHM addresses are recycled by
/// the FairMQ allocator once a chunk is freed, two distinct CCDB paths could
/// transiently share the same id at different points in time. Within a single
/// timeframe an earlier path's deserialisation could then overwrite the
/// `idToObject` slot a later path's `matcherToId` was still pointing at,
/// turning the next `delete reinterpret_cast<T*>(idToObject[oldId])` into a
/// destructor call on the wrong object type. Storing the (id, object) pair
/// per path closes that hole: every path looks up its own slot only.
struct ObjectCache {
  struct Id {
    int64_t value;
    static Id fromRef(DataRef& ref)
    {
      return {reinterpret_cast<int64_t>(ref.payload)};
    }
    bool operator==(const Id& other) const
    {
      return value == other.value;
    }

    struct hash_fn {
      std::size_t operator()(const Id& id) const
      {
        return id.value;
      }
    };
  };

  /// Per-path cache entry for a deserialised CCDB object.
  /// `id` is the version marker — compared against the incoming message's id
  /// to detect "did this object change?". `obj` is the heap-owned, type-
  /// erased pointer to the deserialised value; the path that put it here is
  /// the only one allowed to delete it.
  struct Entry {
    Id id{0};
    void* obj{nullptr};
  };

  /// Per-path cache entry for the CCDB metadata map.
  struct MetadataEntry {
    Id id{0};
    std::map<std::string, std::string> metadata;
  };

  /// Path -> (id, object). Replaces the former matcherToId / idToObject pair.
  std::unordered_map<std::string, Entry> matcherToEntry;

  /// Path -> (id, metadata). Replaces the former matcherToMetadataId /
  /// idToMetadata pair.
  std::unordered_map<std::string, MetadataEntry> matcherToMetadata;
};

} // namespace o2::framework

#endif // O2_FRAMEWORK_OBJECTCACHE_H_
