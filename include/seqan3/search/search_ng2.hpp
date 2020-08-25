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
#include <fmt/format.h>

namespace seqan3 {

/*!\addtogroup submodule_fm_index
 * \{
 */

/*!\brief The SeqAn Bidirectional FM Index Cursor.
 * \implements seqan3::bi_fm_index_cursor_specialisation
 * \tparam index_t The type of the underlying index; must model seqan3::bi_fm_index_specialisation.
 * \details
 *
 * The cursor's interface provides searching a string both from left to right as well as from right to left in the
 * indexed text. It extends the interface of the unidirectional seqan3::fm_index_cursor.
 * All methods modifying the cursor (e.g. extending by a character with extend_right()) return a `bool` value whether
 * the operation was successful or not. In case of an unsuccessful operation the cursor remains unmodified, i.e. an
 * cursor can never be in an invalid state except default constructed cursors that are always invalid.
 *
 * \if DEV
 *     The behaviour is equivalent to a prefix and suffix tree with the space and time efficiency of the underlying pure
 *     FM indices. The cursor traverses the implicit prefix and suffix trees beginning at the root node. The implicit
 *     prefix and suffix trees are not compacted, i.e. going down an edge using extend_right(char) will increase the
 *     query by only one character.
 * \endif
 *
 * The asymptotic running times for using the cursor depend on the SDSL index configuration. To determine the exact
 * running times, you have to additionally look up the running times of the used traits (configuration).
 */
template <typename index_t>
class bi_fm_index_cursor_ng2
{
public:
    using size_type     = typename index_t::size_type;
    using alphabet_type = typename index_t::alphabet_type;
    using rank_type     = std::decay_t<decltype(std::declval<alphabet_type>().to_rank())>;

private:
public: //!TODO

    index_t const* index {};

    size_type fwd_lb {};
    size_type rev_lb {};
    size_type range  {};

public:
    /*!\name Constructors, destructor and assignment
     * \{
     */
    //!\brief Default constructor. Accessing member functions on a default constructed object is undefined behavior.
    //        Default construction is necessary to make this class semi-regular and e.g., to allow construction of
    //        std::array of cursors.
    bi_fm_index_cursor_ng2() noexcept                                           = default;
    bi_fm_index_cursor_ng2(bi_fm_index_cursor_ng2 const &) noexcept             = default;
    bi_fm_index_cursor_ng2 & operator=(bi_fm_index_cursor_ng2 const &) noexcept = default;
    bi_fm_index_cursor_ng2(bi_fm_index_cursor_ng2 &&) noexcept                  = default;
    bi_fm_index_cursor_ng2 & operator=(bi_fm_index_cursor_ng2 &&) noexcept      = default;
    ~bi_fm_index_cursor_ng2()                                                   = default;

    //! \brief Construct from given index.
    bi_fm_index_cursor_ng2(index_t const & _index) noexcept
      : index {&_index}
      , range {index->size()}
    {
#ifndef NDEBUG
		for (rank_type r{1}; r <= alphabet_size<alphabet_type>; ++r) {
			assert(index->fwd_fm.index.char2comp[r] == index->rev_fm.index.char2comp[r]);
		}
#endif
	}

protected:
public://!TODO
    bi_fm_index_cursor_ng2(index_t const & _index, size_type _fwd_lb, size_type _rev_lb, size_type _range) noexcept
      : index  {&_index}
      , fwd_lb {_fwd_lb}
      , rev_lb {_rev_lb}
      , range {_range}
    {}

public:

	bool overlap(bi_fm_index_cursor_ng2 const& _other) const noexcept {
		if (rev_lb + range < _other.rev_lb) {
			return false;
		}
		if (_other.rev_lb + _other.range < rev_lb) {
			return false;
		}
		return true;
	}
	auto join(bi_fm_index_cursor_ng2 const& _other) const noexcept {
		assert(overlap(_other));
		auto new_rev_lb = std::min(rev_lb, _other.rev_lb);
		auto new_rev_rb = std::max(rev_lb + range, _other.rev_lb + _other.range);
		auto new_range = new_rev_rb - new_rev_lb;
		//!TODO fwd_lb is wrong, this should be a uni cursor
		return bi_fm_index_cursor_ng2{*index, fwd_lb, new_rev_lb, new_range};
	}

	bool operator<(bi_fm_index_cursor_ng2 const& _other) const noexcept {
		//!TODO
		return rev_lb < _other.rev_lb;
	}


    //\}

    /*!\brief Tries to extend the query by the character `c` to the right.
     * \tparam char_t Type of the character; needs to be convertible to the character type `char_type` of the index.
     * \param[in] c Character to extend the query with to the right.
     * \returns `true` if the cursor could extend the query successfully.
     *
     * ### Complexity
     *
     * \f$O(T_{BACKWARD\_SEARCH})\f$
     *
     * ### Exceptions
     *
     * No-throw guarantee.
     */
    auto extend_right(rank_type const c) const noexcept
    {
        assert(index != nullptr);
        auto&      csa            = index->fwd_fm.index;
        auto const c_begin        = csa.C[c];
        auto const [rank_l, s, b] = csa.wavelet_tree.lex_count(fwd_lb, fwd_lb + range, c);
        return bi_fm_index_cursor_ng2{*index, c_begin + rank_l, rev_lb + s, range -b -s};
    }

    /*!\brief Tries to extend the query by the character `c` to the left.
     * \tparam char_t Type of the character needs to be convertible to the character type `char_type` of the index.
     * \param[in] c Character to extend the query with to the left.
     * \returns `true` if the cursor could extend the query successfully.
     *
     * ### Complexity
     *
     * \f$O(T_{BACKWARD\_SEARCH})\f$
     *
     * ### Exceptions
     *
     * No-throw guarantee.
     */
    auto extend_left(rank_type const c) const noexcept
    {
        assert(index != nullptr);
        auto&      csa            = index->rev_fm.index;
        auto const c_begin        = csa.C[c];
        auto const [rank_l, s, b] = csa.wavelet_tree.lex_count(rev_lb, rev_lb + range, c);
        return bi_fm_index_cursor_ng2{*index, fwd_lb + s, c_begin + rank_l, range -b -s};
    }

    template <typename CB>
    void extend_right_cb(CB const& cb) const noexcept {
        assert(index != nullptr);
//    	fmt::print("coming from right\n");
		auto& csa = index->fwd_fm.index;
		csa.wavelet_tree.lex_count_fast_cb(fwd_lb, fwd_lb + range, [&](size_t rank_l, size_t s, size_t b, size_t, size_type c) noexcept {
	        size_type const c_begin = csa.C[c];
	        //fmt::print("report_r_cb {} {} ({} {} {}) ({} {} {})\n", c, c_begin, rank_l, s, b, fwd_lb, rev_lb, range);
	        cb(c, bi_fm_index_cursor_ng2{*index, c_begin + rank_l, rev_lb + s, range -b -s});
		});
    }

    template <typename CB>
    void extend_left_cb(CB const& cb) const noexcept
    {
        assert(index != nullptr);
//    	fmt::print("coming from left\n");
        auto& csa = index->rev_fm.index;
		csa.wavelet_tree.lex_count_fast_cb(rev_lb, rev_lb +range, [&](size_t rank_l, size_t s, size_t b, size_t, size_type c) noexcept {
	        size_type const c_begin = csa.C[c];
	        //fmt::print("report_l_cb {} {} ({} {} {}) ({} {} {})\n", c, c_begin, rank_l, s, b, fwd_lb, rev_lb, range);
	        cb(c, bi_fm_index_cursor_ng2{*index, fwd_lb + s, c_begin + rank_l, range -b -s});
		});
    }

    bool valid() const noexcept {
        return range > 0;
    }





    /*!\brief Counts the number of occurrences of the searched query in the text.
     * \returns Number of occurrences of the searched query in the text.
     *
     * ### Complexity
     *
     * Constant.
     *
     * ### Exceptions
     *
     * No-throw guarantee.
     */
    size_type count() const noexcept
    {
        return range;
    }


    template <typename delegate_t>
    void locate(delegate_t delegate) const
    //!\cond
        requires index_t::text_layout_mode == text_layout::collection
    //!\endcond
    {
        assert(index != nullptr);
        for (size_type i = 0; i < count(); ++i)
        {
            size_type loc               = index->rev_fm.index[rev_lb + i];
            size_type sequence_rank     = index->rev_fm.text_begin_rs.rank(loc + 1);
            size_type sequence_position = loc - index->rev_fm.text_begin_ss.select(sequence_rank);
            delegate(sequence_rank - 1, sequence_position);
        }

    }

    template <typename query_alphabet_type>
    auto convert(query_alphabet_type const query_c) const -> size_type {
        static_assert(std::convertible_to<query_alphabet_type, alphabet_type>,
                     "The character must be convertible to the alphabet of the index.");

        auto index_c = to_rank(static_cast<alphabet_type>(query_c)) + 1;
        if constexpr(!std::same_as<alphabet_type, sdsl::plain_byte_alphabet>)
        {
            assert(index->fwd_fm.index.char2comp[index_c] == index->rev_fm.index.char2comp[index_c]);
            assert(not (index_c == 0 && index->fwd_fm.index.char2comp[index_c] > 0));
            return index->fwd_fm.index.char2comp[index_c];
        } else {
            return index_c;
        }
    }

    /*template <typename char_t>
    auto convert(std::vector<char_t> const& query) -> std::vector<size_t> {
        using fwd_alphabet_t = typename decltype(index->fwd_fm.index)::alphabet_type;
        using rev_alphabet_t = typename decltype(index->rev_fm.index)::alphabet_type;

        static_assert(std::convertible_to<char_t, index_alphabet_type>,
                     "The character must be convertible to the alphabet of the index.");
        static_assert(std::is_same_v<fwd_alphabet_t, rev_alphabet_t>);
//        static_assert(std::is_same_v<fwd_alphabet_t, index_alphabet_type>);

        auto realQuery = std::vector<size_t>(query.size());
        for (size_t k{0}; k < query.size(); ++k) {
            auto c = to_rank(static_cast<index_alphabet_type>(query[k])) + 1;
            size_t cc = c;
            if constexpr(!std::same_as<fwd_alphabet_t, sdsl::plain_byte_alphabet>)
            {
                assert(index->fwd_fm.char2comp[c] == index->rev_fm.char2comp[c]);
                cc = index->fwd_fm.char2comp[c];
                assert(not (cc == 0 && c > 0));
            }
            realQuery[k] = cc;
        }
        return realQuery;
    }*/

};

//!\}

} // namespace seqan3


namespace seqan3
{

template <bool fast_lex, bool use_scoring_matrix, int cacheLevel, typename cursor_t, typename query_t, typename search_scheme_t, typename delegate_t, typename scoring_matrix_t>
struct Search_ng2 {
	query_t const& query;
	std::vector<int> const& dir;
	search_scheme_t const& search;
	size_t qidx;
	delegate_t const& delegate;
	scoring_matrix_t const& scoring_matrix;

	using index_alphabet_type = typename cursor_t::alphabet_type;
	using query_alphabet_type = range_innermost_value_t<query_t>;
    using index_rank_type     = std::decay_t<decltype(std::declval<index_alphabet_type>().to_rank())>;

	constexpr static size_t index_sigma = alphabet_size<index_alphabet_type>;

	static constexpr size_t fac(size_t f) {
		if (f > 0) {
			return 6 * fac(f-1);
		}
		return 1;
	}
	//inline static std::array<std::optional<cursor_t>, fac(cacheLevel)> cache2;
	using Cache2 = std::array<std::optional<cursor_t>, fac(cacheLevel)>;
	inline static std::unique_ptr<Cache2> cache2;

//	inline static std::unordered_map<size_t, cursor_t> cache;


	auto get_cursor(int id, int pos, query_t const& _query) -> cursor_t {
		if (not cache2->at(id).has_value()) {
			auto cursor = get_cursor(id / 6, pos-1, _query);
			if (cursor.valid()) {
				auto rank_c = cursor.convert((query[search.pi[pos]]));
				cache2->at(id) = extend(cursor, pos, rank_c);
			} else {
				cache2->at(id) = cursor;
			}
		}
		return cache2->at(id).value();
	}

	Search_ng2(cursor_t const& _cursor, query_t const& _query, std::vector<int> const& _dir, search_scheme_t const& _search, size_t _qidx, size_t _search_idx, delegate_t const& _delegate, scoring_matrix_t const& _scoring_matrix)
		: query    {_query}
		, dir      {_dir}
		, search   {_search}
		, qidx     {_qidx}
		, delegate {_delegate}
		, scoring_matrix {_scoring_matrix}
	{
		if constexpr (cacheLevel > 0) {
			auto cursor = cache2->at(0).value();
			std::size_t id{};
			for (auto i{0}; i < cacheLevel; ++i) {
				auto rank_c = cursor.convert((query[search.pi[i]]));
				id = id * 6 + rank_c;
			}
			cursor = get_cursor(id, cacheLevel-1, _query);

			search_next(cursor, 0, cacheLevel);

/*			auto iter = cache.find(id);
			if (iter == end(cache)) {
				auto cursor = _cursor;
				for (int i{0}; i < cacheLevel and cursor.valid(); ++i) {
					auto rank_c = cursor.convert(query[search.pi[i]]);
					cursor = extend(cursor, i, rank_c);
				}
				cache[id] = cursor;
				search_next(cursor, 0, cacheLevel);
			} else {
				search_next(iter->second, 0, cacheLevel);
			}*/
		} else {
			search_next(_cursor, 0, 0);
		}
	}

	auto extend(cursor_t const& cur, int pos, index_rank_type rank) const noexcept {
		assert(pos >= 0 and pos < dir.size());
		if (dir[pos] > 0) {
			return cur.extend_right(rank);
		}
		return cur.extend_left(rank);
	}

	template <typename CB>
	void extend_cb(cursor_t const& cur, int pos, CB const& cb) const noexcept {
		assert(pos >= 0 and pos < dir.size());
		if (dir[pos] > 0) {
			cur.extend_right_cb(cb);
		} else {
			cur.extend_left_cb(cb);
		}
	}


	void search_next(cursor_t const& cur, int e, size_t pos) const noexcept {
		if (not cur.valid()) {
			return;
		}

		if (pos == query.size()) {
			delegate(qidx, cur);
			return;
		}
		assert(pos >= 0);
		assert(pos < search.l.size());
		assert(pos < search.u.size());
		assert(pos < search.pi.size());
		assert(search.pi[pos] >= 0);
		assert(search.pi[pos] < query.size());

		if constexpr (use_scoring_matrix and fast_lex) {
			auto rank_c = cur.convert(query[search.pi[pos]]);

			extend_cb(cur, pos, [&](auto c, auto newCur) {
				if (c > index_sigma) return;
				auto score = scoring_matrix[c-1][rank_c-1];
				if (search.l[pos] <= e+score and e+score <= search.u[pos]) {
					search_next(newCur, e+score, pos+1);
				}
			});


		} else if constexpr (use_scoring_matrix) {
			auto rank_c = cur.convert(query[search.pi[pos]]);

			for (index_rank_type i{1}; i <= index_sigma; ++i) {
				auto score = scoring_matrix[i-1][rank_c-1];
				if (search.l[pos] <= e+score and e+score <= search.u[pos]) {
					auto newCur = extend(cur, pos, i);
					search_next(newCur, e+score, pos+1);
				}
			}
		} else if constexpr (fast_lex) {
			auto rank_c = cur.convert(query[search.pi[pos]]);

			if (search.l[pos] <= e and e+1 <= search.u[pos]) {
				extend_cb(cur, pos, [&](auto c, auto newCur) {
					if (c > index_sigma) return;
					int error = (c != rank_c);
					search_next(newCur, e+error, pos+1);
				});

			} else if (search.l[pos] <= e and e <= search.u[pos]) {
				auto newCur = extend(cur, pos, rank_c);
				search_next(newCur, e, pos+1);

			} else if (search.l[pos] <= e+1 and e+1 <= search.u[pos]) {
				extend_cb(cur, pos, [&](auto c, auto newCur) {
					if (c > index_sigma) return;
					if (c == rank_c) return;
					search_next(newCur, e+1, pos+1);
				});
			}
		} else {
			auto rank_c = cur.convert(query[search.pi[pos]]);

			if (search.l[pos] <= e and e+1 <= search.u[pos]) {
				for (index_rank_type i{1}; i <= index_sigma; ++i) {
					int error = (i != rank_c);
					auto newCursor = extend(cur, pos, i);
					search_next(newCursor, e+error, pos+1);
				}
			} else if (search.l[pos] <= e and e <= search.u[pos]) {
				auto newCur = extend(cur, pos, rank_c);
				search_next(newCur, e, pos+1);
			} else if (search.l[pos] <= e+1 and e+1 <= search.u[pos]) {
				for (index_rank_type i{1}; i <= index_sigma; ++i) {
					if (i == rank_c) continue;
					auto newCur = extend(cur, pos, i);
					search_next(newCur, e+1, pos+1);
				}
			}
		}
	}
};

template <bool fast_lex = false, bool use_scoring_matrix = false, int cache_level = 0, typename index_t, typename queries_t, typename search_schemes_t, typename delegate_t, typename sm_t>
void search_ng2(index_t const & index, queries_t && queries, uint8_t _max_error, search_schemes_t const & search_scheme, delegate_t && delegate, sm_t const& sm)
{
    if (search_scheme.empty()) return;

    auto len = queries[0].size();
    auto internal_delegate = [&delegate, len] (size_t qidx, auto const & it)
    {
        delegate(qidx, it);
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

	using S = Search_ng2<fast_lex, use_scoring_matrix, cache_level, std::decay_t<decltype(rootCursor)>, decltype(queries[0]), decltype(search_scheme[0]), decltype(internal_delegate), sm_t>;

    S::cache2 = std::make_unique<typename S::Cache2>();

    for (size_t j{0}; j < search_scheme.size(); ++j) {
        auto const& search = search_scheme[j];
        auto const& dir   = dirs[j];

		*S::cache2 = {};
		S::cache2->at(0) = rootCursor;
        for (size_t i{0}; i < queries.size(); ++i) {
        auto const& query = queries[i];
        //auto realQuery = rootCursor.convert(query);
            Search_ng2<fast_lex, use_scoring_matrix, cache_level, std::decay_t<decltype(rootCursor)>, decltype(query), decltype(search), decltype(internal_delegate), sm_t>{rootCursor, query, dir, search, i, j, internal_delegate, sm};
        }
    }
}

//!\}

} // namespace seqan3
