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

#include <seqan3/core/algorithm/configuration.hpp>
#include <seqan3/range/views/persist.hpp>
#include <seqan3/search/search.hpp>
#include <seqan3/search/detail/search_traits.hpp>
#include <seqan3/search/fm_index/all.hpp>
#include <seqan3/std/algorithm>
#include <seqan3/std/ranges>
#include <array>

#include <sdsl/suffix_trees.hpp>

#include <seqan3/alphabet/all.hpp>
#include <seqan3/core/type_traits/range.hpp>
#include <seqan3/range/views/join.hpp>
#include <seqan3/range/views/slice.hpp>
#include <seqan3/search/fm_index/bi_fm_index.hpp>
#include <seqan3/std/ranges>

#include "search_ng2.hpp"

namespace seqan3
{

enum class Dir_ng7 : uint8_t { Left, Right };
template <typename T>
struct Block_ng7 {
	T       rank;
	uint8_t l;
	uint8_t u;
	Dir_ng7 dir;
};


template <typename index_t, typename search_scheme_t, typename delegate_t>
struct Search_ng7 {
	using cursor_t = bi_fm_index_cursor_ng2<index_t>;

	index_t const& index;
	search_scheme_t const& search;
	size_t qidx;
	delegate_t const& delegate;

	using BlockIter = typename search_scheme_t::const_iterator;

	mutable std::vector<char> actions;

	using index_alphabet_type = typename cursor_t::alphabet_type;
    using index_rank_type     = std::decay_t<decltype(std::declval<index_alphabet_type>().to_rank())>;

	constexpr static size_t index_sigma = alphabet_size<index_alphabet_type>;


	Search_ng7(index_t const& _index, search_scheme_t const& _search, size_t _qidx, delegate_t const& _delegate)
		: index     {_index}
		, search    {_search}
		, qidx      {_qidx}
		, delegate  {_delegate}
	{
		search_next(bi_fm_index_cursor_ng2<index_t>{index}, 0, search.begin());
	}

	auto extend(cursor_t const& cur, Dir_ng7 dir, index_rank_type rank) const noexcept {
		if (dir == Dir_ng7::Right) {
			return cur.extend_right(rank);
		}
		return cur.extend_left(rank);
	}

	void search_next(cursor_t const& cur, int e, BlockIter blockIter) const noexcept {
		if (cur.count() == 0) {
			return;
		}

		if (blockIter == end(search)) {
			delegate(qidx, cur, actions);
			return;
		}

		auto rank = blockIter->rank;

		// match
		if (blockIter->l <= e and e <= blockIter->u) {
			auto newCur = extend(cur, blockIter->dir, rank);
			actions.push_back('M');
			search_next(newCur, e, blockIter+1);
			actions.pop_back();
		}

		/*if (pos > 0 and (actions.back() == 'M' or actions.back() == 'S') and search.l[pos-1] <= e+1 and e+1 <= search.u[pos-1]) {
			for (index_rank_type i{1}; i <= index_sigma; ++i) {
				auto newCur = extend(cur, pos-1, i);
				actions.push_back('D');
				search_next(newCur, e+1, pos-1); // deletion occurred in query
				actions.pop_back();
			}

			actions.push_back('I');
			search_next(cur, e+1, pos+1); // insertion occurred in query
			actions.pop_back();
		}*/


		if (blockIter->l <= e+1 and e+1 <= blockIter->u) {
			for (index_rank_type i{1}; i < rank; ++i) {
				auto newCur = extend(cur, blockIter->dir, i);
				actions.push_back('D');
				search_next(newCur, e+1, blockIter); // deletion occurred in query
				actions.pop_back();
				actions.push_back('S');
				search_next(newCur, e+1, blockIter+1); // as substitution
				actions.pop_back();
			}

			/*{
				index_rank_type i = rank;
				auto newCur = extend(cur, blockIter->dir, i);
				actions.push_back('D');
				search_next(newCur, e+1, blockIter); // deletion occurred in query
				actions.pop_back();
			}*/

			for (index_rank_type i(rank+1); i <= index_sigma; ++i) {
				auto newCur = extend(cur, blockIter->dir, i);
				actions.push_back('D');
				search_next(newCur, e+1, blockIter); // deletion occurred in query
				actions.pop_back();
				actions.push_back('S');
				search_next(newCur, e+1, blockIter+1); // as substitution
				actions.pop_back();
			}


			actions.push_back('I');
			search_next(cur, e+1, blockIter+1); // insertion occurred in query
			actions.pop_back();
		}
	}
};


template <typename index_t, typename queries_t, typename search_schemes_t, typename delegate_t>
void search_ng7(index_t const & index, queries_t && queries, uint8_t _max_error, search_schemes_t const & search_scheme, delegate_t && delegate)
{
	auto len = queries[0].size();
    auto internal_delegate = [&delegate, len] (size_t qidx, auto const & it, auto const& actions)
    {
        it.locate([&](auto p1, auto p2) {
            delegate(qidx, p1, p2/*, actions*/);
        }, len);
    };

    std::vector<std::vector<Block_ng7<size_t>>> search_scheme2;
    for (auto s : search_scheme) {
        std::vector<Block_ng7<size_t>> search2;
        for (size_t i{0}; i < s.pi.size(); ++i) {
            auto dir = [&]() {
                if (i == 0) {
                    return s.pi[i] < s.pi[i+1]?Dir_ng7::Right:Dir_ng7::Left;
                } else {
                    return s.pi[i-1] < s.pi[i]?Dir_ng7::Right:Dir_ng7::Left;
                }
            }();
            search2.emplace_back(Block_ng7<size_t>{{}, s.l[i], s.u[i], dir});
        }
        search_scheme2.emplace_back(move(search2));
    }

	using query_alphabet_t = range_innermost_value_t<queries_t>;
    for (size_t i{0}; i < queries.size(); ++i) {
        auto const& query = queries[i];
        for (size_t j{0}; j < search_scheme.size(); ++j) {
            auto& search = search_scheme2[j];
            for (size_t k {0}; k < search.size(); ++k) {
                search[k].rank = index.convert(query[search_scheme[j].pi[k]]);
            }
            Search_ng7{index, search, i, internal_delegate};
        }
    }

}

//!\}

} // namespace seqan3
