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

#include <sdsl/suffix_trees.hpp>

#include <seqan3/alphabet/all.hpp>
#include <seqan3/core/range/type_traits.hpp>
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
 *
 * =================TODO================
 * This documentation might by outdated (by SimonGG)
 * =====================================
 *
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
class bi_fm_index_cursor_ng
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
    bi_fm_index_cursor_ng() noexcept                                          = default;
    bi_fm_index_cursor_ng(bi_fm_index_cursor_ng const &) noexcept             = default;
    bi_fm_index_cursor_ng & operator=(bi_fm_index_cursor_ng const &) noexcept = default;
    bi_fm_index_cursor_ng(bi_fm_index_cursor_ng &&) noexcept                  = default;
    bi_fm_index_cursor_ng & operator=(bi_fm_index_cursor_ng &&) noexcept      = default;
    ~bi_fm_index_cursor_ng()                                                  = default;

    //! \brief Construct from given index.
    bi_fm_index_cursor_ng(index_t const & _index) noexcept
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
    bi_fm_index_cursor_ng(index_t const & _index, size_type _fwd_lb, size_type _rev_lb, size_type _range) noexcept
      : index  {&_index}
      , fwd_lb {_fwd_lb}
      , rev_lb {_rev_lb}
      , range {_range}
    {}

public:



    //\}

    /*!\brief Tries to extend the query by the character `c` to the left.
     * \param[in] c Character to extend the query with to the left.
     * \returns returns a new cursor covering the extended query
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
        return bi_fm_index_cursor_ng{*index, c_begin + rank_l, rev_lb + s, range -b -s};
    }

    /*!\brief Tries to extend the query by the character `c` to the left.
     * \param[in] c Character to extend the query with to the left.
     * \returns returns a new cursor covering the extended query
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
        return bi_fm_index_cursor_ng{*index, fwd_lb + s, c_begin + rank_l, range -b -s};
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
    void locate(delegate_t && delegate) const
    //!\cond
        requires (index_t::text_layout_mode == text_layout::single)
    //!\endcond
    {
        assert(index != nullptr);

        for (size_type i = 0; i < count(); ++i)
        {
            size_type sequence_position = index->rev_fm.index[rev_lb + i];
            delegate(0, sequence_position);
        }
    }



    template <typename delegate_t>
    void locate(delegate_t && delegate) const
    //!\cond
        requires (index_t::text_layout_mode == text_layout::collection)
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

    /*!TODO
     * This function converts a query_c to a fitting rank_type
     * This function should not be part of the iterator, this belongs into the index.
     */
    template <typename query_alphabet_type>
    auto convert(query_alphabet_type const query_c) const noexcept -> rank_type {
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
};

//!\}

} // namespace seqan3
