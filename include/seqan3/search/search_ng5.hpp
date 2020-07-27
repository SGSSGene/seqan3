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

template <typename cursor_t, typename query_t, typename search_scheme_t, typename delegate_t>
struct Search_ng5 {
	query_t const& query;
	std::vector<int> const& dir;
	search_scheme_t const& search;
	size_t qidx;
	delegate_t const& delegate;

	using index_alphabet_type = typename cursor_t::alphabet_type;
	using query_alphabet_type = range_innermost_value_t<query_t>;
    using index_rank_type     = std::decay_t<decltype(std::declval<index_alphabet_type>().to_rank())>;

	constexpr static size_t index_sigma = alphabet_size<index_alphabet_type>;

	static constexpr bool sparse_search = false;


	Search_ng5(cursor_t const& _cursor, query_t const& _query, std::vector<int> const& _dir, search_scheme_t const& _search, size_t _qidx, delegate_t const& _delegate)
		: query     {_query}
		, dir       {_dir}
		, search    {_search}
		, qidx      {_qidx}
		, delegate  {_delegate}
	{
		if constexpr (sparse_search) {
			search_next<false, false, true>(_cursor, 0, 0);
		} else {
			search_next<true, true, true>(_cursor, 0, 0);
		}
	}

	auto extend(cursor_t const& cur, int pos, index_rank_type rank) const noexcept {
		assert(pos >= 0 and pos < dir.size());
		if (dir[pos] > 0) {
			return cur.extend_right(rank);
		}
		return cur.extend_left(rank);
	}

	template <bool substitution, bool insertion, bool deletion>
	void search_next(cursor_t const& cur, int e, size_t pos) const noexcept {
		if (not cur.valid()) {
			return;
		}
		if (pos == query.size()) {
			//if constexpr(substitution) {
				delegate(qidx, cur);
			//}
			return;
		}


		// match
		if (search.l[pos] <= e and e <= search.u[pos]) {
			auto rank = cur.convert(query[search.pi[pos]]);
			auto newCur = extend(cur, pos, rank);
			search_next<true, true, true>(newCur, e, pos+1);
		}

		if (search.l[pos] <= e+1 and e+1 <= search.u[pos]) {
			if constexpr (substitution or insertion) {
				auto rank = cur.convert(query[search.pi[pos]]);
				#if 1
				for (index_rank_type i{1}; i < rank; ++i) {
					auto newCur = extend(cur, pos, i);
					if constexpr (substitution) {
						search_next<true, true, true>(newCur, e+1, pos+1); // as substitution
					}
					if constexpr (insertion) {
						if constexpr (sparse_search) {
							search_next<false, true, false>(newCur, e+1, pos);   // as insertion
						} else {
							search_next<true, true, true>(newCur, e+1, pos);   // as insertion
						}
					}
				}

				if constexpr (insertion) {
					auto newCur = extend(cur, pos, rank);
					//search_next<false, true, false>(newCur, e+1, pos);   // as insertion
					//search_next<true, true, true>(newCur, e+1, pos);   // as insertion
				}
				for (index_rank_type i(rank+1); i <= index_sigma; ++i) {
					auto newCur = extend(cur, pos, i);
					if constexpr (substitution) {
						search_next<true, true, true>(newCur, e+1, pos+1); // as substitution
					}
					if constexpr (insertion) {
						if constexpr (sparse_search) {
							search_next<false, true, false>(newCur, e+1, pos);   // as insertion
						} else {
							search_next<true, true, true>(newCur, e+1, pos);   // as insertion
						}
					}
				}

				/*for (index_rank_type i{1}; i <= index_sigma; ++i) {
					if (i == rank) continue;

					auto newCur = extend(cur, pos, i);
					if constexpr (substitution) {
						search_next<true, true, true>(newCur, e+1, pos+1); // as substitution
					}
					if constexpr (insertion) {
						search_next<false, true, false>(newCur, e+1, pos);   // as insertion
						//search_next<true, true, true>(newCur, e+1, pos);   // as insertion
					}
				}*/
			}
			// deletion
			if constexpr (deletion) {
				//search_next<false, false, true>(cur, e+1, pos+1);
				search_next<true, true, true>(cur, e+1, pos+1);
			}
			#endif
		}
	}
};

template <typename index_t, typename queries_t, typename search_schemes_t, typename delegate_t>
void search_ng5(index_t const & index, queries_t && queries, uint8_t _max_error, search_schemes_t const & search_scheme, delegate_t && delegate)
{
    if (search_scheme.empty()) return;

    auto len = queries[0].size();
    auto internal_delegate = [&delegate, len] (size_t qidx, auto const & it)
    {
        it.locate([&](auto p1, auto p2) {
            delegate(qidx, p1, p2);
        }, len);
    };

    std::vector<std::vector<int>> dirs;
    for (auto const& search : search_scheme) {
        dirs.push_back(std::vector<int>{1});
        for (size_t i{1}; i < search.pi.size(); ++i) {
            dirs.back().push_back(search.pi[i] > search.pi[i-1]?1:-1);
        }
    }

	using query_alphabet_t = range_innermost_value_t<queries_t>;
    auto rootCursor = bi_fm_index_cursor_ng2<index_t>{index};
    for (size_t i{0}; i < queries.size(); ++i) {
        auto const& query = queries[i];
        for (size_t j{0}; j < search_scheme.size(); ++j) {
            auto const& search = search_scheme[j];
            auto const& dir   = dirs[j];
            Search_ng5<std::decay_t<decltype(rootCursor)>, decltype(query), decltype(search), decltype(internal_delegate)>{rootCursor, query, dir, search, i, internal_delegate};
        }
    }

}

//!\}

} // namespace seqan3
