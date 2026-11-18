"""Pure coordinate helpers for one-turn circular replicon intervals."""


def interval_segments(start, end, length):
    """Return inclusive native segments for an unwrapped interval.

    Coordinates are 1-based and ``end`` may exceed ``length`` only by one
    turn.  Returned segments are ordered along the unwrapped interval.
    """
    start, end, length = int(start), int(end), int(length)
    if length < 1 or start < 1 or start > length or end < start or end > start + length - 1:
        raise ValueError("invalid single-turn circular interval")
    if end <= length:
        return [(start, end)]
    return [(start, length), (1, end - length)]


def circular_slice(dna, start, end):
    """Extract an unoriented one-turn interval from a native DNA string."""
    segments = interval_segments(start, end, len(dna))
    return "".join(dna[left - 1:right] for left, right in segments)


def valid_cds_interval(gene, allow_unknown_length=True):
    """Check a called CDS against its declared native coordinate system."""
    start, end = int(gene["start"]), int(gene["end"])
    length = int(gene.get("contig_length") or 0)
    if start < 1 or end < start or length < 0:
        return False
    if not length:
        return allow_unknown_length
    if end <= length:
        return True
    return (str(gene.get("topology", "")).lower() == "circular" and
            start <= length and end <= start + length - 1)


def circular_gaps(intervals, length):
    """Return unoccupied gaps between the union of native intervals."""
    segments = []
    for start, end in intervals:
        segments.extend(interval_segments(start, end, length))
    if not segments:
        return []
    merged = []
    for start, end in sorted(segments):
        if merged and start <= merged[-1][1] + 1:
            merged[-1] = (merged[-1][0], max(merged[-1][1], end))
        else:
            merged.append((start, end))
    gaps = [right[0] - left[1] - 1 for left, right in zip(merged, merged[1:])]
    gaps.append(length - merged[-1][1] + merged[0][0] - 1)
    return gaps


def circular_span(intervals, length):
    """Return the shortest occupied circular arc covering intervals."""
    gaps = circular_gaps(intervals, length)
    return length - max(gaps) if gaps else 0


def circular_overlap(left, right, length):
    """Whether two one-turn intervals share any native base."""
    return any(a <= d and c <= b
               for a, b in interval_segments(*left, length)
               for c, d in interval_segments(*right, length))
