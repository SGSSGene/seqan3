// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <gtest/gtest.h>

#include <algorithm>
#include <type_traits>

#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <seqan3/alphabet/quality/phred42.hpp>
#include <seqan3/alphabet/quality/qualified.hpp>
#include <seqan3/range/views/persist.hpp>
#include <seqan3/search/configuration/hit.hpp>
#include <seqan3/search/configuration/max_error.hpp>
#include <seqan3/search/configuration/on_result.hpp>
#include <seqan3/search/search_ng.hpp>
#include <seqan3/test/expect_range_eq.hpp>
#include "helper.hpp"

using seqan3::operator""_dna4;
using seqan3::operator""_phred42;

TEST(search_ng_test, demo) {
    using index_t = seqan3::bi_fm_index<seqan3::dna4, seqan3::text_layout::single>;

    seqan3::dna4_vector text{"ACGTACGTACGT"_dna4};
    index_t index{text};

    // Search for "ACGT" with 0 errors, should find 3 positions
    {
        auto result = search_ng("ACGT"_dna4, index, 0);
        std::sort(begin(result), end(result)); // This step only needed to make the result reliable
        EXPECT_RANGE_EQ(result, (std::vector<std::size_t>{0, 4, 8}));
    }
    // Search for "ACGG" with 0 errors, should find 0 positions
    {
        auto result = search_ng("ACGG"_dna4, index, 0);
        std::sort(begin(result), end(result)); // This step only needed to make the result reliable
        EXPECT_RANGE_EQ(result, (std::vector<std::size_t>{}));
    }
    // Search for "ACGG" with 1 errors, should find 3 unique positions
    // it will find 9 positions, since multiple alignments are possible.
    // Substitue or insertion
    // ACGT  or ACG  or AC G
    // |||S     |||I    ||I|
    // ACGG     ACGG    ACGG
    {
        auto result = search_ng("ACGG"_dna4, index, 1);
        std::sort(begin(result), end(result)); // This step only needed to make the result reliable
        EXPECT_RANGE_EQ(result, (std::vector<std::size_t>{0, 0, 0, 4, 4, 4, 8, 8, 8}));
    }

}
