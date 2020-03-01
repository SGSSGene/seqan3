#pragma once

namespace seqan3 {

struct Iter {
	size_t idx;
	size_t start;
	size_t end;

	bool valid() const {
		return start <= end and end < std::numeric_limits<size_t>::max()-100;
	}
};

template <bool reversed, typename alphabet_t>
struct Index {
	size_t sequenceLength;
	size_t sequenceCount;
	std::vector<size_t> db;
	std::vector<size_t> jump;
	std::vector<size_t> totalJump;
	static constexpr size_t alphabet_size = alphabet_t::alphabet_size;

	Index(std::vector<std::vector<alphabet_t>> const& sequences) {
		sequenceLength = sequences[0].size();
		sequenceCount = sequences.size();

		db        = std::vector<size_t>((sequenceLength+1) * sequenceCount, 0);
		jump      = std::vector<size_t>(db.size() * alphabet_size, 0);
		totalJump = std::vector<size_t>(alphabet_size*(sequenceLength+1), 0);

		auto currentBegin = db.size();
		// sort past last column (just use initial order)
		{
			auto i = sequenceLength;
			currentBegin -= sequenceCount;
			//fmt::print("has char amounts: {}\n", fmt::join(amount, ", "));
			for (size_t j{0}; j < sequenceCount; ++j) {
				db[currentBegin+j] = j;
			}
		}
		// sort all other columns
		for (auto i{sequenceLength}; i > 0; --i) {
			currentBegin -= sequences.size();
			auto amount = std::vector<int>(alphabet_size+1, 0);
			for (auto const& s : sequences) {
				if constexpr (not reversed) {
					amount[s[i-1].to_rank()+1] += 1;
				} else {
					amount[s[sequenceLength - i].to_rank()+1] += 1;
				}

			}
			for (size_t j{1}; j < amount.size(); ++j) {
				amount[j] = amount[j-1] + amount[j];
			}
			//fmt::print("has char amounts: {}\n", fmt::join(amount, ", "));
			for (size_t j{0}; j < sequences.size(); ++j) {
				auto lastIdx    = currentBegin + sequences.size() + j;
				auto const& sid = db[lastIdx];
				auto const& s   = sequences[sid];
				auto const c    = [&]() {
					if constexpr (not reversed) {
						return s[i-1].to_rank();
					} else {
						return s[sequenceLength - i].to_rank();
					}
				}();
				auto& a = amount[c];
				auto iter = db.begin() + currentBegin + a;
				db[currentBegin+a] = sid;

				++a;

				//fmt::print("j:{}, lastIdx:{}, sid:{}, a:({})\n", j, lastIdx, sid, fmt::join(amount, ", "));

				totalJump[(i-1)*alphabet_size+c] += 1;

				// !TODO Where do all the other systems jump to????
				for (size_t k{0}; k < alphabet_size; ++k) {
					jump[lastIdx*alphabet_size+k] = totalJump[(i-1)*alphabet_size+k];
				}
			}
		}
		// set missing jumps into the first column
		{
			auto i = sequenceLength;
			// !TODO Where do all the other systems jump to????
			for (size_t k{0}; k < alphabet_size; ++k) {
				jump[k] = totalJump[i*alphabet_size+k];
			}

		}
		totalJump.insert(totalJump.begin(), 0);
		// integrate over totalJump
		for (size_t i{1}; i < totalJump.size(); ++i) {
			totalJump[i] = totalJump[i-1] + totalJump[i];
		}
	}

	auto begin(size_t idx) const -> Iter {
		auto seq_count = sequenceCount;
		return Iter{idx, seq_count * idx, (seq_count * idx) + seq_count-1};
	}
	auto next(Iter const& iter, size_t c) const -> Iter {
		if (not iter.valid()) {
			return iter;
		}
		auto total = totalJump[(iter.idx-1)*alphabet_size + c];

		auto newIdx = total + jump[(iter.start-1)*alphabet_size + c];
		if ((iter.start % sequenceCount) == 0) {
			newIdx = total + 0;
		}

		auto endIdx = total + jump[iter.end  *alphabet_size + c]-1;

		return Iter{iter.idx-1, newIdx, endIdx};
	}
};

struct BiIter {
	Iter fwdIter;
	Iter revIter;

	bool valid() const noexcept {
		return fwdIter.valid();
	}
	size_t stepsLeft() const noexcept {
		return fwdIter.idx + revIter.idx;
	}
	size_t size() const noexcept {
		return fwdIter.end - fwdIter.start +1;
	}
};


template <typename alphabet_t>
struct BiIndex {
	using alphabet_type = alphabet_t;

	Index<false, alphabet_t> fwdIndex;
	Index<true, alphabet_t>  revIndex;
	BiIndex(std::vector<std::vector<alphabet_t>> const& sequences)
		: fwdIndex(sequences)
		, revIndex(sequences)
	{}

	auto begin(size_t idx) const noexcept -> BiIter {
		return BiIter{fwdIndex.begin(idx), revIndex.begin(revIndex.sequenceLength - idx)};
	}
	auto extend_left(BiIter const& iter, size_t c) const noexcept -> BiIter {
		if (not iter.valid()) {
			return iter;
		}
		auto newFwd = fwdIndex.next(iter.fwdIter, c);
		size_t smaller = 0;
		for (size_t i{0}; i < c; ++i) {
			auto tIter = fwdIndex.next(iter.fwdIter, i);
			smaller += tIter.end - tIter.start + 1;
		}
		auto newRev = iter.revIter;
		newRev.start += smaller;
		newRev.end = newRev.start + newFwd.end - newFwd.start;
		return BiIter{newFwd, newRev};
	}
	auto extend_right(BiIter const& iter, size_t c) const noexcept -> BiIter {
		if (not iter.valid()) {
			return iter;
		}
		auto newRev = revIndex.next(iter.revIter, c);
		size_t smaller = 0;
		for (size_t i{0}; i < c; ++i) {
			auto tIter = revIndex.next(iter.revIter, i);
			smaller += tIter.end - tIter.start + 1;
		}
		auto newFwd = iter.fwdIter;
		newFwd.start += smaller;
		newFwd.end = newFwd.start + newRev.end - newRev.start;
		return BiIter{newFwd, newRev};
	}


};

}
