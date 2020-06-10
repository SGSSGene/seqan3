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


namespace seqan3
{

// forward declaration
template <typename cursor_t, typename query_t, typename search_t, typename blocks_length_t, typename delegate_t>
void search_ss_ng(cursor_t cur, query_t const& query,
                  typename cursor_t::size_type const lb, typename cursor_t::size_type const rb,
                  uint8_t const errors_spent, uint8_t const block_id, bool const go_right, search_t const & search,
                  blocks_length_t const & blocks_length, detail::search_param const error_left, delegate_t && delegate);


template <typename cursor_t, typename query_t, typename search_t, typename blocks_length_t,
          typename delegate_t>
void search_ss_exact_ng(cursor_t cur, query_t const & query,
                            typename cursor_t::size_type const lb, typename cursor_t::size_type const rb,
                            uint8_t const errors_spent, uint8_t const block_id, bool const go_right,
                            search_t const & search, blocks_length_t const & blocks_length,
                            detail::search_param const error_left, delegate_t && delegate)
{
    using size_type = typename cursor_t::size_type;

    uint8_t const block_id2 = std::min<uint8_t>(block_id + 1, search.blocks() - 1);
    bool const go_right2 = (block_id < search.blocks() - 1) && (search.pi[block_id + 1] > search.pi[block_id]);

    if (go_right)
    {
        size_type const infix_lb = rb - 1; // inclusive
        size_type const infix_rb = lb + blocks_length[block_id] - 1; // exclusive

        if (!cur.extend_right(query | views::slice(infix_lb, infix_rb + 1)))
            return;

        search_ss_ng(cur, query, lb, infix_rb + 2, errors_spent, block_id2, go_right2, search, blocks_length, error_left, delegate);
    }
    else
    {
        size_type const infix_lb = rb - blocks_length[block_id] - 1; // inclusive
        size_type const infix_rb = lb - 1; // inclusive

        if (!cur.extend_left(query | views::slice(infix_lb, infix_rb + 1)))
            return;

        search_ss_ng(cur, query, infix_lb, rb, errors_spent, block_id2, go_right2, search, blocks_length, error_left, delegate);
    }
}

/*!\brief Searches a query sequence in a bidirectional index using a single search of a search schemes.
 *        Sub-function for deletions at the end of a block.
 *
 * \copydetails search_ss_exact
 */
template <typename cursor_t, typename query_t, typename search_t, typename blocks_length_t,
          typename delegate_t>
void search_ss_deletion_ng(cursor_t cur, query_t const & query,
                           typename cursor_t::size_type const lb, typename cursor_t::size_type const rb,
                           uint8_t const errors_spent, uint8_t const block_id, bool const go_right,
                           search_t const & search, blocks_length_t const & blocks_length,
                           detail::search_param const error_left, delegate_t && delegate)
{
    uint8_t const max_error_left_in_block = search.u[block_id] - errors_spent;
    uint8_t const min_error_left_in_block = std::max(search.l[block_id] - errors_spent, 0);

    // Switch to the next block when the min number of errors is reached
    if (min_error_left_in_block == 0)
    {
        uint8_t const block_id2 = std::min<uint8_t>(block_id + 1, search.blocks() - 1);
        bool const go_right2 = search.pi[block_id2] > search.pi[block_id2 - 1];

        search_ss_ng(cur, query, lb, rb, errors_spent, block_id2, go_right2, search, blocks_length, error_left, delegate);
    }

    // Insert deletions into the current block as long as possible
    // Do not allow deletions at the beginning of the leftmost block
    // Do not allow deletions at the end of the rightmost block
    if (!(search.pi[block_id] == 1 && !go_right) &&
        !(search.pi[block_id] == search.blocks() && go_right) &&
        max_error_left_in_block > 0 && error_left.total > 0 && error_left.deletion > 0 &&
        ((go_right && cur.extend_right()) || (!go_right && cur.extend_left())))
    {
        detail::search_param error_left2{error_left};
        error_left2.total--;
        error_left2.deletion--;
        do
        {
            search_ss_deletion_ng(cur, query, lb, rb, errors_spent + 1, block_id, go_right, search, blocks_length, error_left2, delegate);
        } while ((go_right && cur.cycle_back()) || (!go_right && cur.cycle_front()));
    }
}

/*!\brief Searches a query sequence in a bidirectional index using a single search of a search schemes.
 *        Sub-function for approximate search step (iterating over all children in a conceptual suffix tree).
 *
 * \copydetails search_ss_exact
 *
 * \param[in] min_error_left_in_block Number of remaining errors that need to be spent in the current block.
 */
template <typename cursor_t, typename query_t, typename search_t, typename blocks_length_t, typename delegate_t>
void search_ss_children_ng(cursor_t cur, query_t const & query,
                               typename cursor_t::size_type const lb, typename cursor_t::size_type const rb,
                               uint8_t const errors_spent, uint8_t const block_id, bool const go_right,
                               uint8_t const min_error_left_in_block, search_t const & search,
                               blocks_length_t const & blocks_length, detail::search_param const error_left,
                               delegate_t && delegate)
{
    using size_type = typename cursor_t::size_type;
    if ((go_right && cur.extend_right()) || (!go_right && cur.extend_left()))
    {
        size_type const chars_left = blocks_length[block_id] - (rb - lb - 1);

        size_type lb2 = lb - !go_right;
        size_type rb2 = rb + go_right;

        do
        {
            bool const delta = cur.last_rank() != to_rank(query[(go_right ? rb : lb) - 1]);

            // skip if there are more min errors left in the current block than characters in the block
            // i.e. chars_left - 1 < min_error_left_in_block - delta
            // TODO: move that outside the if / do-while struct
            // TODO: incorporate error_left.deletion into formula
            if (error_left.deletion == 0 && chars_left + delta < min_error_left_in_block + 1u)
                continue;

            if (!delta || error_left.substitution > 0)
            {
                detail::search_param error_left2{error_left};
                error_left2.total -= delta;
                error_left2.substitution -= delta;

                // At the end of the current block
                if (rb - lb == blocks_length[block_id])
                {
                    // Leave the possibility for one or multiple deletions at the end of a block.
                    // Thus do not change the direction (go_right) yet.
                    if (error_left.deletion > 0)
                    {
                        //search_ss_deletion_ng(cur, query, lb2, rb2, errors_spent + delta, block_id, go_right, search, blocks_length, error_left2, delegate);
                    }
                    else
                    {
                        uint8_t const block_id2 = std::min<uint8_t>(block_id + 1, search.blocks() - 1);
                        bool const go_right2 = search.pi[block_id2] > search.pi[block_id2 - 1];

                        search_ss_ng(cur, query, lb2, rb2, errors_spent + delta, block_id2, go_right2, search, blocks_length, error_left2, delegate);
                    }
                }
                else
                {
                    search_ss_ng(cur, query, lb2, rb2, errors_spent + delta, block_id, go_right, search, blocks_length, error_left2, delegate);
                }
            }

            // Deletion
            // TODO: check whether the conditions for deletions at the beginning/end of the query are really necessary
            // No deletion at the beginning of the leftmost block.
            // No deletion at the end of the rightmost block.
            if (error_left.deletion > 0 &&
                !(go_right && (rb == 1 || rb == std::ranges::size(query) + 1)) &&
                !(!go_right && (lb == 0 || lb == std::ranges::size(query))))
            {
                detail::search_param error_left3{error_left};
                error_left3.total--;
                error_left3.deletion--;
                search_ss_ng(cur, query, lb, rb, errors_spent + 1, block_id, go_right, search, blocks_length, error_left3, delegate);
            }
        } while ((go_right && cur.cycle_back()) || (!go_right && cur.cycle_front()));
    }
}


template <typename cursor_t, typename query_t, typename search_t,
          typename blocks_length_t, typename delegate_t>
void search_ss_ng(cursor_t cur, query_t const & query,
                      typename cursor_t::size_type const lb, typename cursor_t::size_type const rb,
                      uint8_t const errors_spent, uint8_t const block_id, bool const go_right, search_t const & search,
                      blocks_length_t const & blocks_length, detail::search_param const error_left, delegate_t && delegate)
{
    uint8_t const max_error_left_in_block = search.u[block_id] - errors_spent;
    uint8_t const min_error_left_in_block = std::max(search.l[block_id] - errors_spent, 0); // NOTE: changed

    // Done.
    if (min_error_left_in_block == 0 && lb == 0 && rb == std::ranges::size(query) + 1)
    {
        delegate(cur);
    }
    // Exact search in current block.
    else if (((max_error_left_in_block == 0) && (rb - lb - 1 != blocks_length[block_id])) ||
             (error_left.total == 0 && min_error_left_in_block == 0))
    {
        search_ss_exact_ng(cur, query, lb, rb, errors_spent, block_id, go_right, search, blocks_length, error_left, delegate);
    }
    // Approximate search in current block.
    // i.e. blocks_length[block_id] - (rb - lb - (lb != rb)) >= min_error_left_in_block
    else if (error_left.total > 0)
    {
        // Insertion
        /*if (error_left.insertion > 0)
        {
            using size_type = typename cursor_t::size_type;

            size_type const lb2 = lb - !go_right;
            size_type const rb2 = rb + go_right;

            detail::search_param error_left2{error_left};
            error_left2.total--;
            error_left2.insertion--;
            // At the end of the current block
            if (rb - lb == blocks_length[block_id])
            {
                // Leave the possibility for one or multiple deletions at the end of a block.
                // Thus do not change the direction (go_right) yet.
                // TODO: benchmark the improvement on preventing insertions followed by a deletion and vice versa. Does
                // it pay off the additional complexity and documentation for the user? (Note that the user might only
                // allow for insertions and deletion and not for mismatches).
                //search_ss_deletion_ng(cur, query, lb2, rb2, errors_spent + 1, block_id, go_right, search, blocks_length, error_left2, delegate);
            }
            else
            {
                search_ss_ng(cur, query, lb2, rb2, errors_spent + 1, block_id, go_right, search, blocks_length, error_left2, delegate);
            }
        }*/
        search_ss_children_ng(cur, query, lb, rb, errors_spent, block_id, go_right, min_error_left_in_block, search, blocks_length, error_left, delegate);
    }
}


template <typename index_t, typename query_t, typename block_info_t, typename search_scheme_t, typename delegate_t>
void search_ss_ng(index_t const & index, query_t const & query, block_info_t const& block_info, detail::search_param const error_left,
                      search_scheme_t const & search_scheme, delegate_t && delegate)
{
    for (uint8_t search_id = 0; search_id < search_scheme.size(); ++search_id)
    {
        auto const & search = search_scheme[search_id];
        auto const & [blocks_length, start_pos] = block_info[search_id];

        search_ss_ng(
            index.cursor(),           // cursor on the index
            query,                    // query to be searched
            start_pos, start_pos + 1, // infix range already searched (open interval)
                                      // the first character of `query` has the index 1 (not 0)
            0,                        // errors spent
            0,                        // current block id in search scheme
            true,                     // search the first block from left to right
            search, blocks_length,    // search scheme information
            error_left,               // errors left (broken down by error types)
            delegate                  // delegate function called on hit
        );
    }
}


template <typename index_t, typename query_t, typename block_info_t, typename search_schemes_t, typename delegate_t>
void search_single_ng(index_t const & index, query_t const & query, block_info_t const& block_info, uint8_t _max_error, search_schemes_t const & search_schemes, delegate_t && delegate)
{
    //using search_traits_t = search_traits<configuration_t>;

    // retrieve error numbers / rates
    auto max_error = detail::search_param{_max_error, _max_error, 0, 0};

    // TODO: if total not set: max_error.total = max_error.deletion + max_error.substitution + max_error.insertion;
    // TODO: throw exception when any error number or rate is higher than the total error number/rate
    // throw std::invalid_argument("The total number of errors is set to zero while there is a positive number"
    //                             " of errors for a specific error type.");

    // construct internal delegate for collecting hits for later filtering (if necessary)
    auto internal_delegate = [&delegate] (auto const & it)
    {
        for (auto const & text_pos : it.locate()) {
            delegate(text_pos);
        }
    };

    search_ss_ng(index, query, block_info, max_error, search_schemes, internal_delegate);
}


template <typename index_t, typename queries_t, typename search_schemes_t, typename delegate_t>
void search_ng(index_t const & index, queries_t && queries, uint8_t max_error, search_schemes_t const & search_schemes, delegate_t && delegate)
{
    // retrieve cumulative block lengths and starting position
    auto const block_info = detail::search_scheme_block_info(search_schemes, std::ranges::size(queries[0]));

    for (size_t i{0}; i < queries.size(); ++i) {
       search_single_ng(index, queries[i], block_info, max_error, search_schemes, [i, delegate](auto const& it){
           auto [gidx, pos] = it;
           delegate(i, gidx, pos);
       });
    }
}

//!\}

} // namespace seqan3
