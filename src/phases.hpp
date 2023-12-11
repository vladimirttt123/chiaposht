// Copyright 2018 Chia Network Inc

// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at

//    http://www.apache.org/licenses/LICENSE-2.0

// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.

#ifndef SRC_CPP_PHASES_HPP_
#define SRC_CPP_PHASES_HPP_

#pragma once

#include <cstdint>

enum phase_flags : uint8_t {
    ENABLE_BITFIELD = 1 << 0,
    SHOW_PROGRESS = 1 << 1,
		NO_COMPACTION = 1 << 2,
		TABLE_7_FULL_SCAN = 1 << 3,
		BUFFER_AS_CACHE = 1 << 4,
		PARALLEL_READ = 1 << 5,
		PHASE_1_ONLY = 1 << 6,
		SKIP_PHASE_1 = 1 << 7
};

#endif  // SRC_CPP_PHASES_HPP
