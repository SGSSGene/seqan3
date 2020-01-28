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

namespace seqan3
{

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
    using size_type  = typename index_t::size_type;
private:
    using index_alphabet_type = typename index_t::alphabet_type;

    index_t const * index;

    size_type fwd_lb;
    size_type rev_lb;
    size_type length;

public:
    /*!\name Constructors, destructor and assignment
     * \{
     */
    //!\brief Default constructor. Accessing member functions on a default constructed object is undefined behavior.
    //        Default construction is necessary to make this class semi-regular and e.g., to allow construction of
    //        std::array of cursors.
    bi_fm_index_cursor_ng2() noexcept = default;                                       //!< Defaulted.
    bi_fm_index_cursor_ng2(bi_fm_index_cursor_ng2 const &) noexcept = default;             //!< Defaulted.
    bi_fm_index_cursor_ng2 & operator=(bi_fm_index_cursor_ng2 const &) noexcept = default; //!< Defaulted.
    bi_fm_index_cursor_ng2(bi_fm_index_cursor_ng2 &&) noexcept = default;                  //!< Defaulted.
    bi_fm_index_cursor_ng2 & operator=(bi_fm_index_cursor_ng2 &&) noexcept = default;      //!< Defaulted.
    ~bi_fm_index_cursor_ng2() = default;                                               //!< Defaulted.

    //! \brief Construct from given index.
    bi_fm_index_cursor_ng2(index_t const & _index) noexcept : index(&_index),
                                                          fwd_lb(0), rev_lb(0), length(index->size())
    {}
    bi_fm_index_cursor_ng2(index_t const & _index, size_t _fwd_lb, size_t _rev_lb, size_t _length) noexcept : index(&_index),
                                                          fwd_lb(_fwd_lb), rev_lb(_rev_lb), length(_length)
    {}


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
    auto extend_right(size_type const c) const noexcept
    {
        auto& csa = index->fwd_fm.index;
        size_type const c_begin = csa.C[c];
        auto      const [rank_l, s, b] = csa.wavelet_tree.lex_count(fwd_lb, fwd_lb+length, c);
        return bi_fm_index_cursor_ng2{*index, c_begin + rank_l, rev_lb + s, length -b -s};
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
    auto extend_left(size_type c) const noexcept
    {
        auto& csa = index->rev_fm.index;
        size_type const c_begin = csa.C[c];
        auto      const [rank_l, s, b] = csa.wavelet_tree.lex_count(rev_lb, rev_lb+length, c);
        return bi_fm_index_cursor_ng2{*index, fwd_lb + s, c_begin + rank_l, length -b -s};

    }
    bool valid() const noexcept {
        return length > 0;
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
        assert(index != nullptr);
        return length;
    }


    template <typename delegate_t>
    void locate(delegate_t delegate, size_t length) const
    //!\cond
        requires index_t::text_layout_mode == text_layout::collection
    //!\endcond
    {
        assert(index != nullptr);
        auto os = index->size() - length - 1;

        for (size_type i = 0; i < count(); ++i)
        {
            size_type loc               = os - index->fwd_fm.index[fwd_lb + i];
            size_type sequence_rank     = index->fwd_fm.text_begin_rs.rank(loc + 1);
            size_type sequence_position = loc - index->fwd_fm.text_begin_ss.select(sequence_rank);
            delegate(sequence_rank - 1, sequence_position);
        }
    }
    template <typename char_t>
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
    }

};

//!\}

} // namespace seqan3


namespace seqan3
{

template <typename query_t, typename search_scheme_t, typename delegate_t>
struct Search {
	query_t const& query;
	std::vector<int> const& dir;
	search_scheme_t const& search;
	size_t qidx;
	delegate_t const& delegate;

	template <typename cursor_t>
	Search(cursor_t const& _cursor, query_t const& _query, std::vector<int> const& _dir, search_scheme_t const& _search, size_t _qidx, delegate_t const& _delegate)
		: query    {_query}
		, dir      {_dir}
		, search   {_search}
		, qidx     {_qidx}
		, delegate {_delegate}
	{
		search_next(_cursor, 0, 0);
	}

	template <typename cursor_t, typename char_t>
	auto extend(cursor_t&& cur, int pos, char_t c) const noexcept {
		if (dir[pos] > 0) {
			return cur.extend_right(c);
		}
		return cur.extend_left(c);
	}

	template <typename cursor_t>
	void search_next(cursor_t&& cur, int e, size_t pos) const noexcept {
		if (not cur.valid()) {
			return;
		}
		if (pos == query.size()) {
			delegate(qidx, cur);
			return;
		}
		if (search.l[pos] <= e and e <= search.u[pos]) {
			search_next(extend(cur, pos, query[search.pi[pos]]), e, pos+1);
		}

		if (search.l[pos] <= e+1 and e+1 <= search.u[pos]) {
			auto rank = query[search.pi[pos]];
			for (size_t i{1}; i < rank; ++i) {
				search_next(extend(cur, pos, i), e+1, pos+1);
			}
			for (size_t i{rank+1ul}; i < 5ul; ++i) {
				search_next(extend(cur, pos, i), e+1, pos+1);
			}
		}
	}
};

template <typename index_t, typename queries_t, typename search_schemes_t, typename delegate_t>
void search_ng2(index_t const & index, queries_t && queries, uint8_t _max_error, search_schemes_t const & search_scheme, delegate_t && delegate)
{
    auto length = queries[0].size();
    auto internal_delegate = [&delegate, length] (size_t qidx, auto const & it)
    {
        it.locate([&](auto p1, auto p2) {
            delegate(qidx, std::make_pair(p1, p2));
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
    for (size_t i{0}; i < queries.size(); ++i) {
        auto const& query = queries[i];
        auto realQuery = rootCursor.convert(query);
        for (size_t j{0}; j < search_scheme.size(); ++j) {
            auto const& search = search_scheme[j];
            auto const& dir   = dirs[j];
            Search{rootCursor, realQuery, dir, search, i, internal_delegate};
        }
    }
}

//!\}

} // namespace seqan3
