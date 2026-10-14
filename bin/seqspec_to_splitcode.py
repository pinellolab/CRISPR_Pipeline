#!/usr/bin/env python

"""Write a splitcode config that cuts a guide read down to its spacer.

The guide read carries a leader, the spacer and the scaffold. kallisto sees a
20 nt feature as two 19-mers, and in a 90 bp read its jumping heuristic can
stride over both, so a third of the reads that do carry a library spacer never
pseudoalign. Handing kb a read that is only the spacer removes the problem.

Everything here comes from the seqspec, so a library with a different leader,
a different scaffold or a different spacer length needs no code change.

The one thing the seqspec cannot say is whether a declared flank is usable as
an anchor. A leader ending in a homopolymer stutters in the read, so its 3'
boundary moves; anchoring on it slides into the spacer. `--max-homopolymer`
drops such an anchor instead of trusting the declaration.
"""

import argparse
import sys

import yaml

# seqspec region types. The spacer is the sgRNA target; a flank is only an
# anchor if its own sequence is fixed.
SPACER_TYPES = ("sgrna_target", "grna_target", "protospacer")
FIXED = ("fixed", "joined")


def _leaves(node):
    """Yield the leaf regions under node in order: the library as sequenced."""
    if not isinstance(node, dict):
        return
    kids = [r for r in (node.get("regions") or []) if isinstance(r, dict)]
    if not kids:
        yield node
        return
    for r in kids:
        yield from _leaves(r)


class _IgnoreTags(yaml.SafeLoader):
    """Read a seqspec without the seqspec package.

    A seqspec is tagged YAML (!Assay, !Region, !Read). The tags carry no
    information this script needs, so map every one of them onto the plain
    mapping or sequence underneath.
    """


_IgnoreTags.add_multi_constructor(
    "!", lambda loader, suffix, node:
        loader.construct_mapping(node, deep=True)
        if isinstance(node, yaml.MappingNode)
        else loader.construct_sequence(node, deep=True)
        if isinstance(node, yaml.SequenceNode)
        else loader.construct_scalar(node))


def load_modality(path, modality):
    """Return the ordered leaf regions of one modality of a seqspec."""
    with open(path) as fh:
        spec = yaml.load(fh, Loader=_IgnoreTags)
    for top in spec.get("library_spec") or []:
        if top.get("region_id") == modality or top.get("region_type") == modality:
            return list(_leaves(top))
    raise SystemExit(f"no library_spec region for modality {modality!r} in {path}")


def homopolymer_run(seq, at_end):
    """Length of the homopolymer run at one end of seq."""
    if not seq:
        return 0
    s = seq[::-1] if at_end else seq
    n = 1
    while n < len(s) and s[n] == s[0]:
        n += 1
    return n


def build(regions, anchor_len, max_homopolymer):
    """Pick the spacer region and whichever flanking anchors are trustworthy."""
    idx = next((i for i, r in enumerate(regions)
                if (r.get("region_type") or "").lower() in SPACER_TYPES), None)
    if idx is None:
        raise SystemExit("no sgrna_target region in the seqspec: nothing to extract")
    spacer = regions[idx]

    def anchor(region, side):
        """A prefix or suffix of a fixed flank, or None with the reason why not."""
        if region is None:
            return None, "no flanking region (the read starts or ends at the spacer)"
        if (region.get("sequence_type") or "").lower() not in FIXED:
            return None, f"{region['region_id']} is not a fixed sequence"
        seq = (region.get("sequence") or "").upper()
        if not seq:
            return None, f"{region['region_id']} declares no sequence"
        # The bases that touch the spacer are the ones that locate the boundary.
        run = homopolymer_run(seq, at_end=(side == "left"))
        if run > max_homopolymer:
            return None, (f"{region['region_id']} ends in a {run}-base "
                          f"homopolymer, so its boundary with the spacer is not fixed")
        return (seq[-anchor_len:] if side == "left" else seq[:anchor_len]), None

    left, why_left = anchor(regions[idx - 1] if idx > 0 else None, "left")
    right, why_right = anchor(regions[idx + 1] if idx + 1 < len(regions) else None, "right")
    return spacer, left, right, why_left, why_right


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--seqspec", required=True)
    p.add_argument("--modality", default="crispr")
    p.add_argument("--output", required=True, help="splitcode config to write")
    p.add_argument("--anchor-len", type=int, default=20,
                   help="bases of each flank used as the anchor (default 20). The "
                        "whole scaffold is a worse anchor than its start, because a "
                        "late spacer leaves it truncated at the read end.")
    p.add_argument("--distance", type=int, default=2,
                   help="substitutions tolerated in an anchor (default 2)")
    p.add_argument("--max-homopolymer", type=int, default=1,
                   help="longest homopolymer allowed at a flank's spacer boundary "
                        "before the anchor is dropped as unreliable (default 1)")
    p.add_argument("--read-index", type=int, default=1,
                   help="index of the file holding the guide read (default 1)")
    args = p.parse_args()

    regions = load_modality(args.seqspec, args.modality)
    spacer, left, right, why_left, why_right = build(
        regions, args.anchor_len, args.max_homopolymer)

    lo = int(spacer.get("min_len") or 0)
    hi = int(spacer.get("max_len") or lo)
    if not lo:
        raise SystemExit(f"spacer region {spacer['region_id']} declares no length")

    for side, why in (("5'", why_left), ("3'", why_right)):
        if why:
            print(f"note: no {side} anchor -- {why}", file=sys.stderr)

    if left is None and right is None:
        # The read is the spacer already; splitcode has nothing to do.
        print("note: the guide read carries no flanks, so no extraction is needed",
              file=sys.stderr)
        open(args.output, "w").close()
        return

    rows = []
    if left is not None:
        rows.append(("left_anchor", left, f"{args.read_index}:0:0"))
    if right is not None:
        # Constrain the scaffold to start no earlier than the shortest spacer,
        # so a chance early match cannot ask for a negative-length extract.
        rows.append(("right_anchor", right, f"{args.read_index}:{lo}:0"))

    if left is not None and right is not None:
        extract = "{left_anchor}<spacer[%d-%d]>{right_anchor}" % (lo, hi)
    elif right is not None:
        # Without a left boundary the spacer's length cannot be read off the
        # sequence, so take the longest it can be: kallisto matches k-mers, and
        # a shorter spacer carrying extra 5' bases still contains all of its own.
        extract = "<spacer[%d]>{right_anchor}" % hi
    else:
        extract = "{left_anchor}<spacer[%d-%d]>" % (lo, hi)

    with open(args.output, "w") as fh:
        fh.write(f"@extract {extract}\n")
        fh.write("ids\ttags\tdistances\tlocations\n")
        for name, seq, loc in rows:
            fh.write(f"{name}\t{seq}\t{args.distance}\t{loc}\n")

    print(f"spacer {spacer['region_id']} [{lo}-{hi}] nt; "
          f"anchors: 5'={left or '-'} 3'={right or '-'}", file=sys.stderr)


if __name__ == "__main__":
    main()
