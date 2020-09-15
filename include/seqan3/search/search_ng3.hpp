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


template <typename iter_t, typename Pred>
auto move_to_end_ng3(iter_t begin, iter_t end, Pred pred) -> iter_t {
	auto iter = begin;
	auto tmp_end = end;

	if (iter == tmp_end) {
		return tmp_end;
	}

	while (iter+1 != tmp_end) {
		if (pred(*iter)) {
			--tmp_end;
			std::swap(*iter, *tmp_end);
		} else {
			++iter;
		}
	}

	if (pred(*iter)) {
		return iter;
	}
	return tmp_end;
}

namespace seqan3
{

template <typename cursor_t, typename queries_t, typename search_scheme_t, typename delegate_t>
struct Search_ng3 {
	queries_t const& queries;
	using query_t = typename queries_t::value_type;
	std::vector<std::tuple<size_t, int>>& errors;
	std::vector<int> const& dir;
	search_scheme_t const& search;
	delegate_t const& delegate;
	using iter_t = typename std::decay_t<decltype(errors)>::iterator;

	using index_alphabet_type = typename cursor_t::alphabet_type;
	using query_alphabet_type = range_innermost_value_t<query_t>;
    using index_rank_type     = std::decay_t<decltype(std::declval<index_alphabet_type>().to_rank())>;

	constexpr static size_t index_sigma = alphabet_size<index_alphabet_type>;


	Search_ng3(cursor_t const& _cursor, queries_t const& _queries,  std::vector<std::tuple<size_t, int>>& _errors, std::vector<int> const& _dir, search_scheme_t const& _search, delegate_t const& _delegate)
		: queries  {_queries}
		, errors   {_errors}
		, dir      {_dir}
		, search   {_search}
		, delegate {_delegate}
	{
		search_next(_cursor, 0, end(errors));
	}

	auto extend(cursor_t const& cur, int pos, index_rank_type rank) const noexcept {
		assert(pos >= 0 and pos < dir.size());
		if (dir[pos] > 0) {
			return cur.extend_right(rank);
		}
		return cur.extend_left(rank);
	}

	void search_next(cursor_t const cur, size_t pos, iter_t _end) noexcept {
		if (not cur.valid() or begin(errors) == _end) {
			return;
		}
		if (pos == queries[0].size()) {
			for (auto iter = begin(errors); iter != _end; ++iter) {
				auto& [id, e] = *iter;
				delegate(id, cur);
			}
			return;
		}

		auto idx = search.pi[pos];
		for (size_t i{1}; i < index_sigma; ++i) {
			auto newEnd = move_to_end_ng3(begin(errors), _end, [&](auto& p) {
				auto& [id, e] = p;
				auto rank = cur.convert(queries[id][idx]);
				e += (i != rank);
				if (search.l[pos] > e or e > search.u[pos]) {
					e -= (i != rank);
					return true;
				}
				return false;
			});
			if (begin(errors) == newEnd) {
				continue;
			}
			search_next(extend(cur, pos, i), pos+1, newEnd);
			for (auto iter = begin(errors); iter != newEnd; ++iter) {
				auto& [id, e] = *iter;
				auto rank = cur.convert(queries[id][idx]);
				e -= (i != rank);
			}
		}
	}
};

template <typename index_t, typename queries_t, typename search_schemes_t, typename delegate_t>
void search_ng3(index_t const & index, queries_t && queries, uint8_t _max_error, search_schemes_t const & search_scheme, delegate_t && delegate)
{
    if (search_scheme.empty()) return;

    auto length = queries[0].size();
    auto internal_delegate = [&delegate, length] (size_t qidx, auto const & it)
    {
        it.locate([&](auto p1, auto p2) {
            delegate(qidx, p1, p2);
        }, length);
    };
    std::vector<std::vector<int>> dirs;
    for (auto const& search : search_scheme) {
        dirs.push_back(std::vector<int>{1});
        for (size_t i{1}; i < search.pi.size(); ++i) {
            dirs.back().push_back(search.pi[i] > search.pi[i-1]?1:-1);
        }
    }

    auto rootCursor = bi_fm_index_cursor_ng2{index};

    std::vector<std::tuple<size_t, int>> errors;
    for (size_t idx{0}; idx < queries.size(); ++idx) {
        errors.emplace_back(idx, 0);
    }

    for (size_t j{0}; j < search_scheme.size(); ++j) {
        auto const& search = search_scheme[j];
        auto const& dir   = dirs[j];
        Search_ng3{rootCursor, queries, errors, dir, search, internal_delegate};
    }
}

//!\}

} // namespace seqan3
