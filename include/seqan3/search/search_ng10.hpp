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

enum class Dir_ng10 : uint8_t { Left, Right };
template <typename T>
struct Block_ng10 {
	T       rank;
	uint8_t l;
	uint8_t u;
	Dir_ng10 dir;
};

template <bool _S, bool _D, bool _I, char _LastAction>
struct Info_ng10 {
	static constexpr bool S = _S;
	static constexpr bool D = _D;
	static constexpr bool I = _I;
	static constexpr char lastAction = _LastAction;
};

template <typename index_t, typename search_scheme_t, typename delegate_t>
struct Search_ng10 {
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

	struct Stack {
		cursor_t  cursor;
		int       e;
		BlockIter blockIter;

		using SearchF = void(Search_ng10::*)(cursor_t const&, int, BlockIter);
		SearchF searchF;
	};
	std::vector<Stack> stack;

	Search_ng10(index_t const& _index, search_scheme_t const& _search, size_t _qidx, delegate_t const& _delegate)
		: index     {_index}
		, search    {_search}
		, qidx      {_qidx}
		, delegate  {_delegate}
	{
		auto cur       = bi_fm_index_cursor_ng2<index_t>{index};
		auto blockIter = search.begin();

		using _Info = Info_ng10<true, true, true, 'M'>;

		stack.push_back({cur, 0, blockIter, &Search_ng10::search_next<_Info, _Info>});
	}

	bool nextStep() {
		if (stack.empty()) {
			return false;
		}
		auto [cur, e, blockIter, f] = stack.back();
		stack.pop_back();
		(this->*f)(cur, e, blockIter);
		return true;
	};

	template <bool Right>
	static auto extend(cursor_t const& cur, index_rank_type rank) {
		if constexpr (Right) {
			return cur.extend_right(rank);
		} else {
			return cur.extend_left(rank);
		}
	}

	template <typename LInfo, typename RInfo>
	void search_next(cursor_t const& cur, int e, BlockIter blockIter) {
		if (cur.count() == 0) {
			return;
		}

		if (blockIter == end(search)) {
			if constexpr ((LInfo::lastAction == 'M' or LInfo::lastAction == 'I') and (RInfo::lastAction == 'M' or RInfo::lastAction == 'I')) {
				delegate(qidx, cur, actions);
			}
			return;
		}
		if (blockIter->dir == Dir_ng10::Right) {
			search_next_dir<LInfo, RInfo, true>(cur, e, blockIter);
		} else {
			search_next_dir<LInfo, RInfo, false>(cur, e, blockIter);
		}
	}

	template <typename LInfo, typename RInfo, bool Right>
	void search_next_dir(cursor_t const& cur, int e, BlockIter blockIter) {
		constexpr bool Substitution = (Right?RInfo::S:LInfo::S);
		constexpr bool Deletion     = (Right?RInfo::D:LInfo::D);
		constexpr bool Insertion    = (Right?RInfo::I:LInfo::I);


		auto rank = blockIter->rank;

		using OnMatchL      = std::conditional_t<Right, LInfo, Info_ng10<true, true, true, 'M'>>;
		using OnMatchR      = std::conditional_t<Right, Info_ng10<true, true, true, 'M'>, RInfo>;
		using OnSubstituteL = std::conditional_t<Right, LInfo, Info_ng10<true, false, false, 'S'>>;
		using OnSubstituteR = std::conditional_t<Right, Info_ng10<true, false, false, 'S'>, RInfo>;
		using OnDeletionL   = std::conditional_t<Right, LInfo, Info_ng10<true, true, false, 'D'>>;
		using OnDeletionR   = std::conditional_t<Right, Info_ng10<true, true, false, 'D'>, RInfo>;
		using OnInsertionL  = std::conditional_t<Right, LInfo, Info_ng10<true, false, true, 'I'>>;
		using OnInsertionR  = std::conditional_t<Right, Info_ng10<true, false, true, 'I'>, RInfo>;

		// match
		if (blockIter->l <= e and e <= blockIter->u) {
			auto newCur = extend<Right>(cur, rank);
			stack.push_back({newCur, e, blockIter+1, &Search_ng10::search_next<OnMatchL, OnMatchR>});
		}

		if (blockIter->l <= e+1 and e+1 <= blockIter->u) {
			for (index_rank_type i{1}; i < rank; ++i) {
				auto newCur = extend<Right>(cur, i);
				if constexpr (Deletion) {
					// deletion occurred in query
					stack.push_back({newCur, e+1, blockIter, &Search_ng10::search_next<OnDeletionL, OnDeletionR>});
				}
				if constexpr (Substitution) {
					stack.push_back({newCur, e+1, blockIter+1, &Search_ng10::search_next<OnSubstituteL, OnSubstituteR>});

				}
			}

			for (index_rank_type i(rank+1); i <= index_sigma; ++i) {
				auto newCur = extend<Right>(cur, i);

				if constexpr (Deletion) {
					// deletion occurred in query
					stack.push_back({newCur, e+1, blockIter, &Search_ng10::search_next<OnDeletionL, OnDeletionR>});
				}
				if constexpr (Substitution) {
					stack.push_back({newCur, e+1, blockIter+1, &Search_ng10::search_next<OnSubstituteL, OnSubstituteR>});
				}
			}


			if constexpr (Insertion) {
				// insertion occurred in query
				stack.push_back({cur, e+1, blockIter+1, &Search_ng10::search_next<OnInsertionL, OnInsertionR>});
			}
		}
	}
};


template <typename index_t, typename queries_t, typename search_schemes_t, typename delegate_t>
void search_ng10(index_t const & index, queries_t && queries, uint8_t _max_error, search_schemes_t const & search_scheme, delegate_t && delegate)
{
    if (search_scheme.empty()) return;

    auto len = queries[0].size();
    auto internal_delegate = [&delegate, len] (size_t qidx, auto const & it, auto const& actions)
    {
        it.locate([&](auto p1, auto p2) {
            delegate(qidx, p1, p2/*, actions*/);
        }, len);
    };

    std::vector<std::vector<Block_ng10<size_t>>> search_scheme2;
    for (auto s : search_scheme) {
        std::vector<Block_ng10<size_t>> search2;
        for (size_t i{0}; i < s.pi.size(); ++i) {
            auto dir = [&]() {
                if (i == 0) {
                    return s.pi[i] < s.pi[i+1]?Dir_ng10::Right:Dir_ng10::Left;
                } else {
                    return s.pi[i-1] < s.pi[i]?Dir_ng10::Right:Dir_ng10::Left;
                }
            }();
            search2.emplace_back(Block_ng10<size_t>{{}, s.l[i], s.u[i], dir});
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
            Search_ng10 searchClass{index, search, i, internal_delegate};
            while(searchClass.nextStep()) {};
        }
    }

}

//!\}

} // namespace seqan3
