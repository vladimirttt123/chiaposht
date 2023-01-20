//
// Created by Mariano Sorgente on 2020/09/28.
//
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

#ifndef SRC_CPP_ENTRY_SIZES_HPP_
#define SRC_CPP_ENTRY_SIZES_HPP_
#include <math.h>
#include <stdio.h>
#define NOMINMAX

#include "calculate_bucket.hpp"
#include "pos_constants.hpp"
#include "util.hpp"

class EntrySizes {
public:
    static uint32_t GetMaxEntrySize(uint8_t k, uint8_t table_index, bool phase_1_size)
    {
        // This represents the largest entry size that each table will have, throughout the
        // entire plotting process. This is useful because it allows us to rewrite tables
        // on top of themselves without running out of space.
        switch (table_index) {
            case 1:
                // Represents f1, x
                if (phase_1_size) {
                    return Util::ByteAlign(k + kExtraBits + k) / 8;
                } else {
                    // After computing matches, table 1 is rewritten without the f1, which
                    // is useless after phase1.
                    return Util::ByteAlign(k) / 8;
                }
            case 2:
            case 3:
            case 4:
            case 5:
            case 6:
                if (phase_1_size)
                    // If we are in phase 1, use the max size, with metadata.
                    // Represents f, pos, offset, and metadata
                    return Util::ByteAlign(
                               k + kExtraBits + (k) + kOffsetSize +
                               k * kVectorLens[table_index + 1]) /
                           8;
                else
                    // If we are past phase 1, we can use a smaller size, the smaller between
                    // phases 2 and 3. Represents either:
                    //    a:  sort_key, pos, offset        or
                    //    b:  line_point, sort_key
                    return Util::ByteAlign(
                               std::max(static_cast<uint32_t>(2 * k + kOffsetSize),
                                   static_cast<uint32_t>(3 * k - 1))) /
                           8;
            case 7:
            default:
                // Represents line_point, f7
                return Util::ByteAlign(3 * k - 1) / 8;
        }
    }

    // Get size of entries containing (sort_key, pos, offset). Such entries are
    // written to table 7 in phase 1 and to tables 2-7 in phase 2.
    static uint32_t GetKeyPosOffsetSize(uint8_t k)
    {
        return cdiv(2 * k + kOffsetSize, 8);
    }

    // Calculates the size of one C3 park. This will store bits for each f7 between
    // two C1 checkpoints, depending on how many times that f7 is present. For low
    // values of k, we need extra space to account for the additional variability.
    static uint32_t CalculateC3Size(uint8_t k)
    {
        if (k < 20) {
            return Util::ByteAlign(8 * kCheckpoint1Interval) / 8;
        } else {
            return Util::ByteAlign(kC3BitsPerEntry * kCheckpoint1Interval) / 8;
        }
    }

    static uint32_t CalculateLinePointSize(uint8_t k) { return Util::ByteAlign(2 * k) / 8; }

    // This is the full size of the deltas section in a park. However, it will not be fully filled
    static uint32_t CalculateMaxDeltasSize(uint8_t k, uint8_t table_index)
    {
        if (table_index == 1) {
            return Util::ByteAlign((kEntriesPerPark - 1) * kMaxAverageDeltaTable1) / 8;
        }
        return Util::ByteAlign((kEntriesPerPark - 1) * kMaxAverageDelta) / 8;
    }

    static uint32_t CalculateStubsSize(uint32_t k)
    {
        return Util::ByteAlign((kEntriesPerPark - 1) * (k - kStubMinusBits)) / 8;
    }

    static uint32_t CalculateParkSize(uint8_t k, uint8_t table_index)
    {
        return CalculateLinePointSize(k) + CalculateStubsSize(k) +
               CalculateMaxDeltasSize(k, table_index);
    }

		static uint64_t EvaluateTableFileSize( uint8_t k, uint8_t table_index, uint8_t phase ){
			return 0;
		}

		static uint64_t EvalueateBucketFileSize( uint8_t k, uint8_t table_index, uint8_t phase, uint32_t buckets_count, uint32_t bucket_num ){
			return EvaluateTableFileSize( k, table_index, phase ) / (uint64_t)buckets_count;
			/* k=25, -u 16,
					Phase 1:
						t1: 113318K - bucket_sizes: [14330K, 14353K] ~ 1/45 final
						t2: 141622K - bucket_sizes: [30682K, 30756K] ~ 1/21 final
						t3: 141567K - bucket_sizes: [42924K, 43033K]
						t4: 141468K - bucket_sizes: [42909K, 43006K]
						t5: 141264K - bucket_sizes: [36748K, 36810K]
						t6: 140873K - bucket_sizes: [30559K, 30615K]
						t7: 256944K - bucket_sizes: [0, 0]
					Phase 2: t1 & t7 does not touched
						t6:	0		buckets_sizes: 0-12=[17175K, 17198K], 13=1988K, 14-15=0
						t5:	0		buckets_sizes: 0-11=[16663K, 16691K], 12=14795K, 13-15=0
						t4:	0		buckets_sizes: 0-11=[16480K, 16506K], 12=13219K, 13-15=0
						t3:	0		buckets_sizes: 0-11=[16415K, 16443K], 12=12623K, 13-15=0
						t2:	0		buckets_sizes: 0-11=[16380K, 16409K], 12=12384K, 13-15=0
					Phase 3:
						t1: 0
						t2:     buckets_sezes: [16202K, 16833K]

						t3p1:		buckets_sizes: 0-9=[25707K, 25753K], 10=4874K, 11-15=0
						t3p2:		buckets_sizes: 0-11=[14336K], 12=11264K, 13-15=0 (p3s)
						t4p1:		buckets_sizes: 0-9=[25751K, 25798K], 10=6213K 11-15=0
						t4p2:		buckets_sizes: 0-11=[14336K],	12=12288K, 13-15=0
						t5p1:		buckets_sizes: 0-9=[25843K, 25905K], 10=9818K, 11-15=0
						t5p2:		buckets_sizes: 0-12=[14336K],	13=1024K, 14-15=0
						t6p1:		buckets_sizes: 0-9=[26184K, 26212K], 10=19704K, 11-15=0
						t6p2:		buckets_sizes: 0-12=[14336K],	13=10240K, 14-15=0
						t7p1:		buckets_sizes: 0-10=[27454K], 11=22767K,	12-15=0
						t7p2:		buckets_sizes: 0-15=[13312K]

					final: 652687K


				K=32 -u 128
					Phase1:
						t1: bs=0.281GiB
						t2: bs=0.563GiB
						t3: bs=0.813GiB
						t4: bs=0.813GiB
						t5: bs=0.688GiB
						t6: bs=0.563GiB
					Phase2:
						t1: bs=0-101[0.313GiB], 0.025GiB
						t2: bs=0-101[0.313GiB], 0.065GiB
						t3: 0-101[0.315Gib],  0.165GiB
					Phase3:
						t1p1: [0.299, 0.330]
						t1p2: 0-101[0.250GiB], 0.052GiB
						t1p2: 0-80[0.471GiB], 0.431GiB
			*/
		}
};

#endif  // CHIAPOS_ENTRY_SIZES_HPP
