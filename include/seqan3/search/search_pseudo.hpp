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

template <bool EditDistance, typename index_t, typename search_scheme_t, typename query_t, typename delegate_t>
struct Search_pseudo {
	using cursor_t = bi_fm_index_cursor_ng2<index_t>;

	index_t const& index;
	delegate_t const& delegate;

	decltype(search_scheme_t::pi) const& pi;
	decltype(search_scheme_t::l) const& l;
	decltype(search_scheme_t::u) const& u;
	query_t const& query;

	using index_alphabet_type = typename cursor_t::alphabet_type;
    using index_rank_type     = std::decay_t<decltype(std::declval<index_alphabet_type>().to_rank())>;

	constexpr static size_t index_sigma = alphabet_size<index_alphabet_type>;


	Search_pseudo(index_t const& _index, search_scheme_t const& _search, query_t const& _query, delegate_t const& _delegate) noexcept
		: index     {_index}
		, pi{_search.pi}
		, l{_search.l}
		, u{_search.u}
		, query{_query}
		, delegate  {_delegate}
	{
		auto cur       = bi_fm_index_cursor_ng2<index_t>{index};

		if constexpr (EditDistance) {
			search_distance(cur, 0, 0);
		} else {
			search_hm(cur, 0, 0);
		}
	}

	auto extend(cursor_t const& cur, index_rank_type rank, std::size_t pos) const noexcept {
		if (pos == 0 or pi[pos-1] < pi[pos]) {
			return cur.extend_right(rank);
		} else {
			return cur.extend_left(rank);
		}
	}

	void search_hm(cursor_t const& cur, int e, std::size_t pos) const noexcept {

		if (cur.count() == 0) {
			return;
		}

		if (pos == query.size()) {
			if (l[pos-1] <= e and e <= u[pos-1]) {
				delegate(cur);
			}
			return;
		}

		if (e > u[pos]) {
			return;
		}


		// expected next character
		auto rank = index.convert(query[pi[pos]]);

		// cursors extended by one character
		auto cursors = std::array<cursor_t, index_sigma+1>{};
		if (e+1 <= u[pos]) {
			for (auto i{1}; i <= index_sigma; ++i) {
				cursors[i] = extend(cur, i, pos);
			}
		} else {
			cursors[rank] = extend(cur, rank, pos);
		}

		// search matches
		if (l[pos] <= e) {
			search_hm(cursors[rank], e,  pos+1);
		}

		// search substitution
		if (l[pos] <= e+1 and e+1 <= u[pos]) {
			for (index_rank_type i{1}; i <= index_sigma; ++i) {
				if (i == rank) continue;
				search_hm(cursors[i], e+1, pos+1);
			}
		}
	}

	void search_distance(cursor_t const& cur, int e, std::size_t pos) const noexcept {

		if (cur.count() == 0) {
			return;
		}

		if (pos == query.size()) {
			if (l[pos-1] <= e and e <= u[pos-1]) {
				delegate(cur);
			}
			return;
		}

		if (e > u[pos]) {
			return;
		}

		// expected next character
		auto rank = index.convert(query[pi[pos]]);

		// cursors extended by one character
		auto cursors = std::array<cursor_t, index_sigma+1>{};
		if (e+1 <= u[pos]) {
			for (auto i{1}; i <= index_sigma; ++i) {
				cursors[i] = extend(cur, i, pos);
			}
		} else {
			cursors[rank] = extend(cur, rank, pos);
		}

		// search matches
		if (l[pos] <= e) {
			search_distance(cursors[rank], e,  pos+1);
		}

		// search substitution
		if (l[pos] <= e+1 and e+1 <= u[pos]) {
			for (index_rank_type i{1}; i <= index_sigma; ++i) {
				if (i == rank) continue;
				search_distance(cursors[i], e+1, pos+1);
			}
		}

		//search deletions
		if (e+1 <= u[pos]) {
			for (index_rank_type i{1}; i <= index_sigma; ++i) {
				search_distance(cursors[i], e+1, pos);
			}
		}

		// search insertion
		if (l[pos] <= e+1 and e+1 <= u[pos]) {
			search_distance(cur, e+1, pos+1);
		}
	}
};


template <bool EditDistance, typename index_t, typename queries_t, typename search_schemes_t, typename delegate_t, typename sm_t>
void search_pseudo(index_t const & index, queries_t && queries, uint8_t _max_error, search_schemes_t const & search_scheme, delegate_t && delegate, sm_t const& sm)
{
    std::size_t qidx;
    auto internal_delegate = [&qidx, &delegate] (auto const & it)
    {
        delegate(qidx, it);
    };

    for (qidx = {0}; qidx < queries.size(); ++qidx) {
        for (size_t j{0}; j < search_scheme.size(); ++j) {
            Search_pseudo<EditDistance, std::decay_t<decltype(index)>, std::decay_t<decltype(search_scheme[j])>, std::decay_t<decltype(queries[qidx])>, std::decay_t<decltype(internal_delegate)>> {index, search_scheme[j], queries[qidx], internal_delegate};

        }
    }

}

//!\}

} // namespace seqan3
