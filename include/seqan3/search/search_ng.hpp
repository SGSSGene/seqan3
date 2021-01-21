// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \author Christopher Pockrandt <christopher.pockrandt AT fu-berlin.de>
 * \brief Provides the public interface for search algorithms.
 */

#pragma once

#include <seqan3/core/configuration/configuration.hpp>
#include <seqan3/range/views/persist.hpp>
#include <seqan3/search/search.hpp>
#include <seqan3/search/detail/search_traits.hpp>
#include <seqan3/search/fm_index/all.hpp>
#include <seqan3/std/algorithm>
#include <seqan3/std/ranges>
#include <array>
#include <seqan3/search/fm_index/bi_fm_index_cursor_ng.hpp>


#include <sdsl/suffix_trees.hpp>

#include <seqan3/alphabet/all.hpp>
#include <seqan3/core/range/type_traits.hpp>
#include <seqan3/range/views/join.hpp>
#include <seqan3/range/views/slice.hpp>
#include <seqan3/search/fm_index/bi_fm_index.hpp>
#include <seqan3/std/ranges>

namespace seqan3
{

template <typename query_t,
          fm_index_specialisation index_t,
          typename delegate_t>

class Search_ng_impl {
    using cursor_t = bi_fm_index_cursor_ng<index_t>;

    query_t const & query;
    std::size_t     maxAllowedErrors;
    delegate_t      report;


public:
    Search_ng_impl(query_t const & query,
                   index_t const & index,
                   std::size_t     maxAllowedErrors,
                   delegate_t      report)
        : query            {query}
        , maxAllowedErrors {maxAllowedErrors}
        , report           {report}
    {
        auto cur = cursor_t{index};
        search(cur, 0, 0);
    }

    void search(cursor_t cur,
                std::size_t e,
                std::size_t pos) {

        using index_alphabet_type = typename cursor_t::alphabet_type;
        using index_rank_type     = std::decay_t<decltype(std::declval<index_alphabet_type>().to_rank())>;

        constexpr static std::size_t index_sigma = alphabet_size<index_alphabet_type>;


        // cursor is empty, no need to continue search
        if (cur.count() == 0) {
            return;
        }

        // check if it breaks upper limit
        if (e > maxAllowedErrors) {
            return;
        }

        // check if reaches end of search sequence
        if (pos == query.size()) {
            report(cur);
            return;
        }


        // search matches
        // expected next character converted to rank
        auto rank = cur.convert(query[pos]);
        search(cur.extend_right(rank), e,  pos+1);

        // search substitution
        if (e+1 <= maxAllowedErrors) {
            for (index_rank_type i{1}; i <= index_sigma; ++i) {
                if (i == rank) continue; // skip character that is a match
                search(cur.extend_right(i), e+1, pos+1);
            }
        }

        //search deletions
        if (e+1 <= maxAllowedErrors) {
            for (index_rank_type i{1}; i <= index_sigma; ++i) {
                search(cur.extend_right(i), e+1, pos);
            }
        }

        // search insertion
        if (e+1 <= maxAllowedErrors) {
            search(cur, e+1, pos+1);
        }
    }

};

template <fm_index_specialisation index_t,
          typename query_t>
auto search_ng(query_t const & query,
               index_t const & index,
               std::size_t maxAllowedErrors) {

    std::vector<size_t> results;
    // report function that translates cursors into positions
    auto report = [&](auto const& cursor) {
        // translates positions of the cursor into positions in the sequences (since we only have one sequence "ACGTACGTACGT" we know it is always zero
        cursor.locate([&](auto sequenceIdx, auto position) {
            assert(sequenceIdx == 0);
            results.push_back(position);
        });
    };
    // performance the actual backtracking search
    Search_ng_impl(query, index, maxAllowedErrors, report);
    return results;
}


} // namespace seqan3
