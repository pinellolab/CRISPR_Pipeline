import pathlib
import random
import sys
import time
from collections import Counter

REPO_ROOT = pathlib.Path(__file__).resolve().parents[1]
BIN_DIR = REPO_ROOT / "bin"
if str(BIN_DIR) not in sys.path:
    sys.path.insert(0, str(BIN_DIR))

import seqSpecCheck as ssc


def _naive_analyze_guides_in_reads(reads, guide_list):
    """Reference oracle: the original O(reads * guides) implementation,
    kept here only to differentially test the optimized version against."""
    positions = []
    upstream_map = {}
    guide_hits = Counter({g: 0 for g in guide_list})

    if not reads or not guide_list:
        return positions, upstream_map, guide_hits

    for seq in reads:
        for guide in guide_list:
            idx = seq.find(guide)
            if idx != -1:
                positions.append(idx)
                guide_hits[guide] += 1
                start_cut = max(0, idx - 12)
                upstream_fragment = seq[start_cut:idx]
                if len(upstream_fragment) < 12:
                    padding = '-' * (12 - len(upstream_fragment))
                    upstream_fragment = padding + upstream_fragment
                if idx not in upstream_map:
                    upstream_map[idx] = []
                upstream_map[idx].append(upstream_fragment)
                break
    return positions, upstream_map, guide_hits


def _random_seq(rng, length, alphabet="ACGT"):
    return "".join(rng.choice(alphabet) for _ in range(length))


def _read_contains_at_most_one_guide(seq, guide_list):
    hits = [g for g in guide_list if g in seq]
    return len(hits) <= 1


def test_single_guide_match_matches_naive_reference():
    guides = ["ACGTACGTACGTACGTACGT"]  # 20bp
    reads = [
        "NNNNNNNNNN" + guides[0] + "TTTT",  # guide with room upstream
        "AC" + guides[0],  # guide near the very start (needs '-' padding)
        "NNNNNNNNNNNNNNNNNNNNNNNNNN",  # no match
    ]
    expected = _naive_analyze_guides_in_reads(reads, guides)
    observed = ssc.analyze_guides_in_reads(reads, guides)
    assert observed == expected


def test_upstream_padding_when_guide_near_read_start():
    guide = "A" * 20
    read = "CC" + guide  # only 2bp upstream, needs 10 '-' padding chars
    positions, upstream_map, guide_hits = ssc.analyze_guides_in_reads([read], [guide])
    assert positions == [2]
    assert upstream_map[2] == ["----------CC"]
    assert guide_hits[guide] == 1


def test_no_match_returns_empty_results():
    positions, upstream_map, guide_hits = ssc.analyze_guides_in_reads(
        ["NNNNNNNNNNNNNNNNNNNN"], ["A" * 20]
    )
    assert positions == []
    assert upstream_map == {}
    assert guide_hits["A" * 20] == 0


def test_empty_reads_or_guides_short_circuit():
    assert ssc.analyze_guides_in_reads([], ["A" * 20]) == ([], {}, Counter())
    guide_hits = ssc.analyze_guides_in_reads(["ACGT"], [])[2]
    assert guide_hits == Counter()


def test_mixed_length_guides_are_all_found():
    guide19 = "A" * 19
    guide21 = "C" * 21
    reads = [
        "GG" + guide19 + "TT",
        "GG" + guide21 + "TT",
    ]
    positions, upstream_map, guide_hits = ssc.analyze_guides_in_reads(
        reads, [guide19, guide21]
    )
    assert positions == [2, 2]
    assert guide_hits[guide19] == 1
    assert guide_hits[guide21] == 1


def test_finds_leftmost_position_when_read_start_has_earlier_guide():
    # Two distinct 20bp guides; guide_b appears earlier in the read than
    # guide_a, but guide_a is listed first. The optimized implementation
    # finds the leftmost match (guide_b); the naive implementation finds
    # whichever guide (by list order) has any match at all (guide_a), since
    # it never compares positions across different guides. This is the one
    # documented, intentional behavior difference from the original
    # algorithm -- leftmost-position is the more meaningful semantic for
    # position-distribution analysis, which is the actual purpose of
    # `positions`. Real guide libraries are designed to be sequence-distinct,
    # so in practice a read essentially never contains two different full
    # guide sequences -- this scenario doesn't arise on real data.
    guide_a = "A" * 20
    guide_b = "C" * 20
    read = "GG" + guide_b + "TT" + guide_a
    naive_positions, _, _ = _naive_analyze_guides_in_reads([read], [guide_a, guide_b])
    optimized_positions, _, _ = ssc.analyze_guides_in_reads([read], [guide_a, guide_b])
    assert naive_positions == [2 + 20 + 2]  # finds guide_a, listed first
    assert optimized_positions == [2]  # finds guide_b, leftmost


def test_differential_random_reads_with_realistic_unique_guides():
    """On reads containing at most one guide (the realistic case for
    sequence-distinct guide libraries), the optimized implementation must
    match the original exactly."""
    rng = random.Random(0)
    guide_length = 20
    guides = [_random_seq(rng, guide_length) for _ in range(200)]
    guide_set = set(guides)

    reads = []
    for _ in range(300):
        read_len = rng.randint(30, 150)
        if rng.random() < 0.5:
            # Embed exactly one real guide at a random position.
            guide = rng.choice(guides)
            pos = rng.randint(0, read_len - guide_length)
            prefix = _random_seq(rng, pos)
            suffix = _random_seq(rng, read_len - pos - guide_length)
            read = prefix + guide + suffix
        else:
            read = _random_seq(rng, read_len)
        # Regenerate any read that happens to contain more than one guide by
        # chance, so this test stays within the "realistic" case being
        # verified (see test_finds_leftmost_position_... for the documented
        # divergent case).
        while not _read_contains_at_most_one_guide(read, guide_set):
            read = _random_seq(rng, read_len)
        reads.append(read)

    expected = _naive_analyze_guides_in_reads(reads, guides)
    observed = ssc.analyze_guides_in_reads(reads, guides)
    assert observed == expected


def test_optimized_implementation_scales_with_reads_not_guides():
    """Coarse regression guard: with a large guide library, the optimized
    path should stay fast (no guides-count multiplier), while the naive
    reference slows down substantially. Uses a generous margin to avoid
    CI flakiness -- this is a algorithmic-shape check, not a strict
    benchmark."""
    rng = random.Random(1)
    guide_length = 20
    guides = [_random_seq(rng, guide_length) for _ in range(4000)]
    reads = [_random_seq(rng, 150) for _ in range(2000)]

    start = time.perf_counter()
    ssc.analyze_guides_in_reads(reads, guides)
    optimized_elapsed = time.perf_counter() - start

    start = time.perf_counter()
    _naive_analyze_guides_in_reads(reads[:200], guides)  # 1/10th the reads
    naive_partial_elapsed = time.perf_counter() - start

    # The optimized run over ALL reads should still be faster than the naive
    # run over just 1/10th of them -- if this fails, the guides-count
    # multiplier has crept back in.
    assert optimized_elapsed < naive_partial_elapsed
