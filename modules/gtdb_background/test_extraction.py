"""Coordinate and annotation edge cases that can corrupt contextual datasets."""
import unittest
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature,FeatureLocation,BeforePosition
from extract_windows import neighborhood_indices,ordered_cds,decode_feature

class ExtractionTests(unittest.TestCase):
    def test_reverse_center_preserves_all_CDS_positions(self):
        x=neighborhood_indices(30,15,-1)
        self.assertEqual([i for _,i,_ in x],list(range(22,7,-1)))
        self.assertEqual(x[7],(0,15,False))
    def test_linear_end_is_missing_and_circular_wraps_once(self):
        linear=neighborhood_indices(20,1,1)
        self.assertEqual(sum(i is None for _,i,_ in linear),6)
        circular=neighborhood_indices(20,1,1,circular=True)
        self.assertEqual(len({i for _,i,_ in circular}),15)
        self.assertEqual(circular[0],(-7,14,True))
    def test_small_circular_record_does_not_repeat_neighbors(self):
        x=neighborhood_indices(10,0,1,circular=True)
        self.assertTrue(any(i is None for _,i,_ in x))
    def test_invalid_CDS_does_not_collapse_index(self):
        r=SeqRecord(Seq('A'*400))
        r.features=[SeqFeature(FeatureLocation(200,230,strand=1),type='CDS',qualifiers={'translation':['M'+'A'*9]}),
                    SeqFeature(FeatureLocation(50,80,strand=-1),type='CDS',qualifiers={'pseudo':['']}),
                    SeqFeature(FeatureLocation(100,130,strand=1),type='CDS',qualifiers={'translation':['M'+'A'*9]})]
        fs=ordered_cds(r);self.assertEqual(len(fs),3)
        self.assertEqual(decode_feature(fs[0],0,set())['sequence_status'],'pseudo')
        self.assertEqual(int(fs[1].location.start),100)
    def test_source_unknown_and_partial_are_not_negative_labels(self):
        f=SeqFeature(FeatureLocation(BeforePosition(0),30,strand=1),type='CDS',qualifiers={'translation':['M'+'A'*9],'product':['hypothetical protein']})
        g=decode_feature(f,0,set());self.assertTrue(g['source_unknown']);self.assertEqual(g['sequence_status'],'partial_CDS')

if __name__=='__main__':unittest.main()
