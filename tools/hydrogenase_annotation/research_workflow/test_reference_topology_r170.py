import unittest

from reference_reactions import _coordinate_distance, _local_core_pair


class R170ReferenceTopologyTests(unittest.TestCase):
    thresholds = {"local_max_span_bp": 200, "local_max_intergenic_gap_bp": 100}

    def gene(self, start, end, *, circular=False):
        return {"id": str(start), "contig": "ctg_circular-yes" if circular else "ctg",
                "start": start, "end": end, "strand": 1,
                "contig_length": 1000 if circular else 0}

    def test_documented_circular_origin_pair_uses_shortest_arc(self):
        anchor = self.gene(950, 980, circular=True)
        partner = self.gene(10, 40, circular=True)
        self.assertTrue(_local_core_pair(anchor, partner, self.thresholds))
        self.assertEqual(_coordinate_distance(anchor, partner), 60)

    def test_linear_and_incomplete_topology_keep_original_distance_gate(self):
        anchor = self.gene(950, 980)
        partner = self.gene(10, 40)
        self.assertFalse(_local_core_pair(anchor, partner, self.thresholds))
        self.assertEqual(_coordinate_distance(anchor, partner), 940)
        circular = self.gene(950, 980, circular=True)
        invalid = dict(circular, contig_length=0, start=10, end=40)
        self.assertFalse(_local_core_pair(circular, invalid, self.thresholds))
        self.assertEqual(_coordinate_distance(circular, invalid), float('inf'))

    def test_pair_prefilter_keeps_candidates_with_room_for_a_bridge(self):
        anchor = self.gene(925, 940, circular=True)
        partner = self.gene(55, 80, circular=True)
        self.assertTrue(_local_core_pair(anchor, partner, self.thresholds))

    def test_equal_circular_distances_have_same_evidence(self):
        anchor = self.gene(500, 530, circular=True)
        a = self.gene(600, 630, circular=True)
        b = self.gene(400, 430, circular=True)
        self.assertEqual(_coordinate_distance(anchor, a), _coordinate_distance(anchor, b))


if __name__ == "__main__":
    unittest.main()
