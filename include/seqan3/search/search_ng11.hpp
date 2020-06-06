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
#include <seqan3/search/algorithm/detail/search.hpp>
#include <seqan3/search/algorithm/detail/search_traits.hpp>
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

enum class Dir_ng11 : uint8_t { Left, Right };
template <typename T>
struct Block_ng11 {
	T       rank;
	uint8_t l;
	uint8_t u;
	Dir_ng11 dir;
};

template <typename index_t, typename search_scheme_t, typename delegate_t>
struct Search_ng11 {
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


	Search_ng11(index_t const& _index, search_scheme_t const& _search, size_t _qidx, delegate_t const& _delegate)
		: index     {_index}
		, search    {_search}
		, qidx      {_qidx}
		, delegate  {_delegate}
	{
		auto cur       = bi_fm_index_cursor_ng2<index_t>{index};
		auto blockIter = search.begin();

		search_next<'M', 'M'>(cur, 0, blockIter);

	}

	template <bool Right>
	static auto extend(cursor_t const& cur, index_rank_type rank) noexcept {
		if constexpr (Right) {
			return cur.extend_right(rank);
		} else {
			return cur.extend_left(rank);
		}
	}

	template <bool Right, typename CB>
	static auto extend_cb(cursor_t const& cur, CB const& cb) noexcept {
		if constexpr (Right) {
			return cur.extend_right_cb(cb);
		} else {
			return cur.extend_left_cb(cb);
		}
	}


	template <char LInfo, char RInfo>
	void search_next(cursor_t const& cur, int e, BlockIter blockIter) const noexcept {
		if (cur.count() == 0) {
			return;
		}

		if (blockIter == end(search)) {
			if constexpr ((LInfo == 'M' or LInfo == 'I') and (RInfo == 'M' or RInfo == 'I')) {
				delegate(qidx, cur, actions);
			}
			return;
		}
		if (blockIter->dir == Dir_ng11::Right) {
			search_next_dir<LInfo, RInfo, true>(cur, e, blockIter);
		} else {
			search_next_dir<LInfo, RInfo, false>(cur, e, blockIter);
		}
	}

	template <char LInfo, char RInfo, bool Right>
	void search_next_dir(cursor_t const& cur, int e, BlockIter blockIter) const noexcept {
		static constexpr char TInfo = Right ? RInfo : LInfo;

		constexpr bool Deletion     = TInfo == 'M' or TInfo == 'D';
		constexpr bool Insertion    = TInfo == 'M' or TInfo == 'I';

		constexpr char OnMatchL      = Right ? LInfo : 'M';
		constexpr char OnMatchR      = Right ? 'M'   : RInfo;
		constexpr char OnSubstituteL = Right ? LInfo : 'S';
		constexpr char OnSubstituteR = Right ? 'S'   : RInfo;
		constexpr char OnDeletionL   = Right ? LInfo : 'D';
		constexpr char OnDeletionR   = Right ? 'D'   : RInfo;
		constexpr char OnInsertionL  = Right ? LInfo : 'I';
		constexpr char OnInsertionR  = Right ? 'I'   : RInfo;

		auto rank = blockIter->rank;

		bool matchAllowed    = blockIter->l <= e and e <= blockIter->u;
		bool mismatchAllowed = blockIter->l <= e+1 and e+1 <= blockIter->u;

		if (mismatchAllowed) {
			auto cursors = std::array<cursor_t, index_sigma+1>{};
			extend_cb<Right>(cur, [&](auto c, auto newCur) {
				if (c > index_sigma) return;
				cursors[c] = newCur;
			});

			if (matchAllowed) {
				auto newCur = cursors[rank];
				search_next<OnMatchL, OnMatchR>(newCur, e, blockIter+1);
			}

			for (index_rank_type i{1}; i < rank; ++i) {
				auto newCur = cursors[i];

				if constexpr (Deletion) {
					search_next<OnDeletionL, OnDeletionR>(newCur, e+1, blockIter); // deletion occurred in query
				}
				search_next<OnSubstituteL, OnSubstituteR>(newCur, e+1, blockIter+1); // as substitution
			}

			for (index_rank_type i(rank+1); i <= index_sigma; ++i) {
				auto newCur = cursors[i];

				if constexpr (Deletion) {
					search_next<OnDeletionL, OnDeletionR>(newCur, e+1, blockIter); // deletion occurred in query
				}
				search_next<OnSubstituteL, OnSubstituteR>(newCur, e+1, blockIter+1); // as substitution
			}


			if constexpr (Insertion) {
				search_next_ins<OnInsertionL, OnInsertionR, Right>(cur, e+1, blockIter+1, cursors); // insertion occurred in query
			}


		} else if (matchAllowed) {
			auto newCur = extend<Right>(cur, rank);
			search_next<OnMatchL, OnMatchR>(newCur, e, blockIter+1);
		}
	}

	template <char LInfo, char RInfo, bool lastDir>
	void search_next_ins(cursor_t const& cur, int e, BlockIter blockIter, std::array<cursor_t, index_sigma+1> const& cursors) const noexcept {
		if (cur.count() == 0) {
			return;
		}

		if (blockIter == end(search)) {
			if constexpr ((LInfo == 'M' or LInfo == 'I') and (RInfo == 'M' or RInfo == 'I')) {
				delegate(qidx, cur, actions);
			}
			return;
		}
		if (lastDir == (blockIter->dir == Dir_ng11::Right)) {
			if (blockIter->dir == Dir_ng11::Right) {
				search_next_dir_ins<LInfo, RInfo, true>(cur, e, blockIter, cursors);
			} else {
				search_next_dir_ins<LInfo, RInfo, false>(cur, e, blockIter, cursors);
			}
		} else {
			if (blockIter->dir == Dir_ng11::Right) {
				search_next_dir<LInfo, RInfo, true>(cur, e, blockIter);
			} else {
				search_next_dir<LInfo, RInfo, false>(cur, e, blockIter);
			}

		}
	}

	template <char LInfo, char RInfo, bool Right>
	void search_next_dir_ins(cursor_t const& cur, int e, BlockIter blockIter, std::array<cursor_t, index_sigma+1> const& cursors) const noexcept {
		static constexpr char TInfo = Right ? RInfo : LInfo;

		constexpr bool Deletion     = TInfo == 'M' or TInfo == 'D';
		constexpr bool Insertion    = TInfo == 'M' or TInfo == 'I';

		constexpr char OnMatchL      = Right ? LInfo : 'M';
		constexpr char OnMatchR      = Right ? 'M'   : RInfo;
		constexpr char OnSubstituteL = Right ? LInfo : 'S';
		constexpr char OnSubstituteR = Right ? 'S'   : RInfo;
		constexpr char OnDeletionL   = Right ? LInfo : 'D';
		constexpr char OnDeletionR   = Right ? 'D'   : RInfo;
		constexpr char OnInsertionL  = Right ? LInfo : 'I';
		constexpr char OnInsertionR  = Right ? 'I'   : RInfo;

		auto rank = blockIter->rank;

		bool matchAllowed    = blockIter->l <= e and e <= blockIter->u;
		bool mismatchAllowed = blockIter->l <= e+1 and e+1 <= blockIter->u;

		if (matchAllowed) {
			auto newCur = cursors[rank];
			search_next<OnMatchL, OnMatchR>(newCur, e, blockIter+1);
		}

		if (mismatchAllowed) {
			for (index_rank_type i{1}; i < rank; ++i) {
				auto newCur = cursors[i];

				if constexpr (Deletion) {
					search_next<OnDeletionL, OnDeletionR>(newCur, e+1, blockIter); // deletion occurred in query
				}
				search_next<OnSubstituteL, OnSubstituteR>(newCur, e+1, blockIter+1); // as substitution
			}

			for (index_rank_type i(rank+1); i <= index_sigma; ++i) {
				auto newCur = cursors[i];

				if constexpr (Deletion) {
					search_next<OnDeletionL, OnDeletionR>(newCur, e+1, blockIter); // deletion occurred in query
				}
				search_next<OnSubstituteL, OnSubstituteR>(newCur, e+1, blockIter+1); // as substitution
			}


			if constexpr (Insertion) {
				search_next_ins<OnInsertionL, OnInsertionR, Right>(cur, e+1, blockIter+1, cursors); // insertion occurred in query
			}
		}
	}

};


template <typename index_t, typename queries_t, typename search_schemes_t, typename delegate_t>
void search_ng11(index_t const & index, queries_t && queries, uint8_t _max_error, search_schemes_t const & search_scheme, delegate_t && delegate)
{
    auto internal_delegate = [&delegate] (size_t qidx, auto const & it, auto const& actions)
    {
        it.locate([&](auto p1, auto p2) {
            delegate(qidx, p1, p2/*, actions*/);
        });
    };

    std::vector<std::vector<Block_ng11<size_t>>> search_scheme2;
    for (auto s : search_scheme) {
        std::vector<Block_ng11<size_t>> search2;
        for (size_t i{0}; i < s.pi.size(); ++i) {
            auto dir = [&]() {
                if (i == 0) {
                    return s.pi[i] < s.pi[i+1]?Dir_ng11::Right:Dir_ng11::Left;
                } else {
                    return s.pi[i-1] < s.pi[i]?Dir_ng11::Right:Dir_ng11::Left;
                }
            }();
            search2.emplace_back(Block_ng11<size_t>{{}, s.l[i], s.u[i], dir});
        }
        search_scheme2.emplace_back(move(search2));
    }

	using query_alphabet_t = innermost_value_type_t<queries_t>;
    for (size_t i{0}; i < queries.size(); ++i) {
        auto const& query = queries[i];
        for (size_t j{0}; j < search_scheme.size(); ++j) {
            auto& search = search_scheme2[j];
            for (size_t k {0}; k < search.size(); ++k) {
                search[k].rank = index.convert(query[search_scheme[j].pi[k]]);
            }
            Search_ng11{index, search, i, internal_delegate};
        }
    }

}

//!\}

} // namespace seqan3
