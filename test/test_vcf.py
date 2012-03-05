import unittest
import doctest
import os
import commands
from StringIO import StringIO
import sys

import vcf
from vcf import utils

suite = doctest.DocTestSuite(vcf.parser)


def fh(fname):
    return file(os.path.join(os.path.dirname(__file__), fname))


class TestVcfSpecs(unittest.TestCase):

    def test_vcf_4_0(self):
        reader = vcf.Reader(fh('example-4.0.vcf'))
        assert reader.metadata['fileformat'] == 'VCFv4.0'

        # test we can walk the file at least
        for r in reader:

            if r.POS == 1230237:
                assert r.is_monomorphic
            else:
                assert not r.is_monomorphic

            if 'AF' in r.INFO:
                self.assertEqual(type(r.INFO['AF']),  type([]))

            for c in r:
                assert c

                # issue 19, in the example ref the GQ is length 1
                if c.called:
                    self.assertEqual(type(c.data['GQ']),  type(1))
                    if 'HQ' in c.data and c.data['HQ'] is not None:
                        self.assertEqual(type(c.data['HQ']),  type([]))



    def test_vcf_4_1(self):
        return
        reader = vcf.Reader(fh('example-4.1.vcf'))
        self.assertEqual(reader.metadata['fileformat'],  'VCFv4.1')

        # contigs were added in vcf4.1
        # probably need to add a reader.contigs attribute
        assert 'contig' in reader.metadata

        # test we can walk the file at least
        for r in reader:
            for c in r:
                assert c

        # asserting False while I work out what to check
        assert False

    def test_vcf_4_1_sv(self):
        return
       
        reader = vcf.Reader(fh('example-4.1-sv.vcf'))

        assert 'SVLEN' in reader.infos

        # test we can walk the file at least
        for r in reader:
            print r
            for c in r:
                print c
                assert c

        # asserting False while I work out what to check
        assert False

        
class TestGatkOutput(unittest.TestCase):

    filename = 'gatk.vcf'

    samples = ['BLANK', 'NA12878', 'NA12891', 'NA12892',
            'NA19238', 'NA19239', 'NA19240']
    formats = ['AD', 'DP', 'GQ', 'GT', 'PL']
    infos = ['AC', 'AF', 'AN', 'BaseQRankSum', 'DB', 'DP', 'DS',
            'Dels', 'FS', 'HRun', 'HaplotypeScore', 'InbreedingCoeff',
            'MQ', 'MQ0', 'MQRankSum', 'QD', 'ReadPosRankSum']

    n_calls = 37

    def setUp(self):
        self.reader = vcf.Reader(fh(self.filename))

    def testSamples(self):
        self.assertEqual(self.reader.samples, self.samples)

    def testFormats(self):
        self.assertEqual(set(self.reader.formats), set(self.formats))

    def testInfos(self):
        self.assertEqual(set(self.reader.infos), set(self.infos))


    def testCalls(self):
        n = 0

        for site in self.reader:
            n += 1
            self.assertEqual(len(site.samples), len(self.samples))


            # check sample name lookup
            for s in self.samples:
                assert site.genotype(s)

            # check ordered access
            self.assertEqual([x.sample for x in site.samples], self.samples)

        self.assertEqual(n,  self.n_calls)


class TestFreebayesOutput(TestGatkOutput):

    filename = 'freebayes.vcf'
    formats = ['AO', 'DP', 'GL', 'GLE', 'GQ', 'GT', 'QA', 'QR', 'RO']
    infos = ['AB', 'ABP', 'AC', 'AF', 'AN', 'AO', 'BVAR', 'CIGAR',
            'DB', 'DP', 'DPRA', 'EPP', 'EPPR', 'HWE', 'LEN', 'MEANALT',
            'NUMALT', 'RPP', 'MQMR', 'ODDS', 'MQM', 'PAIREDR', 'PAIRED',
            'SAP', 'XRM', 'RO', 'REPEAT', 'XRI', 'XAS', 'XAI', 'SRP',
            'XAM', 'XRS', 'RPPR', 'NS', 'RUN', 'CpG', 'TYPE']
    n_calls = 104


    def testParse(self):
        reader = vcf.Reader(fh('freebayes.vcf'))




        print reader.samples
        self.assertEqual(len(reader.samples), 7)
        n = 0
        for r in reader:
            n+=1
            for x in r:
                assert x
        assert n == self.n_calls

class TestSamtoolsOutput(unittest.TestCase):

    def testParse(self):
        reader = vcf.Reader(fh('samtools.vcf'))

        self.assertEqual(len(reader.samples), 1)
        self.assertEqual(sum(1 for _ in reader), 11)


class Test1kg(unittest.TestCase):

    def testParse(self):
        reader = vcf.Reader(fh('1kg.vcf.gz'))

        self.assertEqual(len(reader.samples), 629)
        for _ in reader:
            pass

class TestReader(unittest.TestCase):

    def setUp(self):
        self.reader = vcf.Reader(fh('issue-order-4.1.vcf'))
        # TODO contig should be a list like filters, formats, etc.
        # and then this test will change
        self.ordered_metadata = [('fileformat', 'VCFv4.1'), ('fileDate', '20090805'), \
            ('source', 'myImputationProgramV3.1'), \
            ('reference', 'file:///seq/references/1000GenomesPilot-NCBI36.fasta'), \
            ('contig', '<ID=20,length=62435964,assembly=B36,md5=f126cdf8a6e0c7f379d618ff66beb2da,species="Homo sapiens",taxonomy=x>'), \
            ('phasing', 'partial'), ("test_metadata", "testing")]
        self.ordered_infos = ["NS", "DP", "AF", "AA", "DB", "H2"] 
        self.ordered_filters = ["s50", "q10"]
        self.ordered_formats = ["GT", "GQ", "DP", "HQ"]
        self.ordered_samples = ["NA00002", "NA00001", "NA00003"]

    def test_metadata_order(self):
        """Tests order and content"""
        self.assertIsNotNone(self.reader) 
        for expected_pair, actual_key in zip(self.ordered_metadata, self.reader.metadata.keys()):
            actual_value = self.reader.metadata[actual_key]
            self.assertEquals(expected_pair[0], actual_key)
            self.assertEquals(expected_pair[1], actual_value)

    def test_infos_order(self):
        """Tests order only TODO content"""
        self.assertIsNotNone(self.reader) 
        self.assertEquals(self.ordered_infos, self.reader.infos.keys())

    def test_filters_order(self):
        """Tests order only TODO content"""
        self.assertIsNotNone(self.reader) 
        self.assertEquals(self.ordered_filters, self.reader.filters.keys())

    def test_formats_order(self):
        """Tests order only TODO content"""
        self.assertIsNotNone(self.reader) 
        self.assertEquals(self.ordered_formats, self.reader.formats.keys())

    def test_samples_order(self):
        """Tests order only TODO content"""
        self.assertIsNotNone(self.reader) 
        self.assertEquals(self.ordered_samples, self.reader.samples)

    def test_pass_none(self):
        """Tests that filter of PASS is None.  If 
            this changes, the Writer needs to change."""
        self.assertIsNotNone(self.reader) 
        filters = []
        for record in self.reader:
            filters.append(record.FILTER)
        # PASS is turned into None
        self.assertEqual(4, filters.count(None))
    
    def test_info_number(self):
        """Tests INFO.Number is parsed and stored correctly."""   
        old_stderr = sys.stderr
        sys.stderr = StringIO()
        reader = vcf.Reader(fh('issue-info-4.1.vcf')) 
        self.assertTrue("INFO[INVALID].Number = -1 which is invalid" in sys.stderr.getvalue())
        sys.stderr.close()
        sys.stderr = old_stderr
        self.assertEquals(1, reader.infos["STR_1"].num) 
        self.assertEquals(2, reader.infos["STR_2"].num) 
        self.assertEquals("A", reader.infos["STR_A"].num) 
        self.assertEquals("G", reader.infos["STR_G"].num) 
        self.assertEquals(".", reader.infos["STR_VARIES"].num) 
        self.assertEquals(0, reader.infos["DB"].num) 
        self.assertEquals(-1, reader.infos["INVALID"].num) 

    def test_format_number(self):
        """Tests FORMAT.Number is parsed and stored correctly."""
        old_stderr = sys.stderr
        sys.stderr = StringIO()
        reader = vcf.Reader(fh('issue-info-4.1.vcf'))
        self.assertTrue("FORMAT[INVALID_ZERO].Number = 0 which is invalid" in sys.stderr.getvalue())
        self.assertTrue("FORMAT[INVALID_NEG].Number = -1 which is invalid" in sys.stderr.getvalue())
        sys.stderr.close()
        sys.stderr = old_stderr
        self.assertEquals(1, reader.formats["GT"].num)
        self.assertEquals(2, reader.formats["HQ"].num)
        self.assertEquals("A", reader.formats["A"].num)
        self.assertEquals("G", reader.formats["G"].num)
        self.assertEquals(".", reader.formats["U"].num)
        self.assertEquals(0, reader.formats["INVALID_ZERO"].num)
        self.assertEquals(-1, reader.formats["INVALID_NEG"].num)

class TestWriter(unittest.TestCase):

    def testWrite(self):

        reader = vcf.Reader(fh('gatk.vcf'))
        out = StringIO()
        writer = vcf.Writer(out, reader)

        records = list(reader)

        map(writer.write_record, records)
        out.seek(0)
        reader2 = vcf.Reader(out)

        self.assertEquals(reader.samples, reader2.samples)
        self.assertEquals(reader.formats, reader2.formats)
        self.assertEquals(reader.infos, reader2.infos)

        for l, r in zip(records, reader2):
            self.assertEquals(l.samples, r.samples)

    def test_header(self):
        in_file = fh('issue-order-4.1.vcf')
        reader = vcf.Reader(in_file)
        out = StringIO()
        writer = vcf.Writer(out, reader)
        # write records as well 
        map(writer.write_record, list(reader))
        out_str = out.getvalue()
        # header lines does not include:
        # #CHROM POS ID REF ALT QUAL FILTER INFO FORMAT"
        self.assertTrue(out_str.startswith("".join(reader._header_lines)))

    def test_header_info(self):
        # suppress errors 
        old_stderr = sys.stderr
        sys.stderr = StringIO()

        reader_in = vcf.Reader(fh('issue-info-4.1.vcf'))
        out = StringIO()
        writer = vcf.Writer(out, reader_in)
        out.seek(0)
        reader_out = vcf.Reader(out)

        set_numbers = set()
        for in_key, out_key in zip(sorted(reader_in.infos.keys()), sorted(reader_out.infos.keys())):
            self.assertEquals(in_key, out_key)
            self.assertEquals(reader_in.infos[in_key].num, reader_out.infos[out_key].num)
            set_numbers.add(reader_out.infos[out_key].num)

        self.assertTrue(set_numbers >= set(["A", "G", 1, 2, ".", 0, -1]))

        # we suppressed errors above
        sys.stderr.close()
        sys.stderr = old_stderr

    def test_record_info(self):
        """Make sure whatever was in the input file 
            is preserved in the output file."""
        old_stderr = sys.stderr
        sys.stderr = StringIO()
        # this file has errors issued we don't want 
        # runners of the unit tests to be concerned with
        reader_in = vcf.Reader(fh('issue-info-4.1.vcf'))
        out = StringIO()

        writer = vcf.Writer(out, reader_in)
        
        in_records = list(reader_in)

        map(writer.write_record, in_records)
        out.seek(0)
        reader_out = vcf.Reader(out)
        # now reset stderr
        sys.stderr.close()
        sys.stderr = old_stderr

        lines = out.getvalue().split("\n")
        self.assertTrue(lines >= 29)

        count_found = 0
        for line in lines:
            if line.startswith("20\t14370"):
                count_found += 1
                self.assertTrue("NS=3;DP=14;AF=0.5;DB;H2;STR_1=Test1;STR_VARIES=0/0,0/1;STR_2=FIRST,SECOND;STR_A=True;STR_G=Something" in line)
            elif line.startswith("20\t17330"):
                count_found += 1
                self.assertTrue("NS=3;DP=11;AF=0.017;STR_A=False;STR_G=Something;STR_1=Test2;STR_2=Aa,Bb;STR_VARIES=0/0,0/1,1/1" in line)
            elif line.startswith("20\t1110696"):
                count_found += 1
                self.assertTrue("NS=2;DP=10;AF=0.333,0.667;AA=T;DB;STR_1=Test3;STR_2=Bb,Cc;STR_A=False,True;STR_G=Something_else" in line)
            elif line.startswith("20\t1230237"):
                count_found += 1
                self.assertTrue("NS=3;DP=13;AA=T;STR_1=Test4;STR_VARIES=0/0;STR_2=FIRST,SECOND;STR_G=Something" in line) 
            elif line.startswith("20\t1234567"):
                count_found += 1
                self.assertTrue("NS=3;DP=9;AA=G;STR_A=False,True;STR_G=Something;STR_1=Test5;STR_2=Aa,Bb;STR_VARIES=0/0,0/1,1/1" in line)

        self.assertEqual(5, count_found)
        out.close() 

    def test_record_filter(self):
        """Make sure whatever is read in is also written out.""" 
        reader_in = vcf.Reader(fh('issue-filter-4.1.vcf')) 
        out = StringIO()

        writer = vcf.Writer(out, reader_in)     
        in_records = list(reader_in)

        map(writer.write_record, in_records)        
        out.seek(0)
        reader_out = vcf.Reader(out)

        out_records = list(reader_out)

        # check records
        in_filters = []
        out_filters = []
        for record in in_records:
            in_filters.append(record.FILTER)

        for record in out_records:
            out_filters.append(record.FILTER)

        self.assertEquals(in_filters, out_filters)

        # make sure that the user gets ['q10', 's50']
        self.assertIn(['q10', 's50'], in_filters)
        self.assertIn(['q10', 's50'], out_filters)

        # check output

        lines = out.getvalue().split("\n")        
        self.assertTrue(lines >= 29) 

        count_found = 0 
        for line in lines:            
            if line.startswith("20\t14370"):
                count_found += 1
                self.assertTrue("29\tPASS\tNS=3" in line)
            elif line.startswith("20\t17330"):
                count_found += 1
                self.assertTrue("3\tq10\tNS=3" in line)
            elif line.startswith("20\t1110696"):
                count_found += 1
                self.assertTrue("67\tPASS\tNS=2" in line)
            elif line.startswith("20\t1230237"):
                count_found += 1
                self.assertTrue("47\tq10;s50\tNS=3" in line) 
            elif line.startswith("20\t1234567"):
                count_found += 1
                self.assertTrue("50\t.\tNS=" in line)

        self.assertEqual(5, count_found)
        out.close()

    def test_header_format(self):
        # suppress errors 
        old_stderr = sys.stderr
        sys.stderr = StringIO()

        reader_in = vcf.Reader(fh('issue-info-4.1.vcf'))
        out = StringIO()
        writer = vcf.Writer(out, reader_in)
        out.seek(0)
        reader_out = vcf.Reader(out)

        set_numbers = set()
        for in_key, out_key in zip(reader_in.formats.keys(), reader_out.formats.keys()):
            self.assertEquals(in_key, out_key)
            self.assertEquals(reader_in.formats[in_key].num, reader_out.formats[out_key].num)
            set_numbers.add(reader_out.formats[out_key].num)

        self.assertTrue(set_numbers >= set(["A", "G", 1, 2, ".", 0, -1]))

        # we suppressed errors above
        sys.stderr.close()
        sys.stderr = old_stderr

    def test_new_metadata(self):
        """Tests that we can add new metadata, info,
        format, and filter fields, that they are included
        in the output header and the order is preserved."""
        reader = vcf.Reader(fh('gatk.vcf'))

        unmodified_number_of_filters = len(reader.filters)
        unmodified_out = StringIO()
        unmodified_writer = vcf.Writer(unmodified_out, reader)
        unmodified_out_str = unmodified_out.getvalue()

        print "********* unmod out str *************"
        print unmodified_out_str
        print "**********************"

        # modify the information our reader has 
        reader.filters["NEW_FILTER"] = vcf.parser._Filter("NEW_FILTER", "Testing")
        print "**********************"
        print reader.filters
        print "**********************"

        reader.formats["NEW_FORMAT"] = vcf.parser._Format("NEW_FORMAT", 1, "String", "Testing")
        reader.infos["NEW_INFO"] = vcf.parser._Info("NEW_INFO", 1, "String", "Testing")
        reader.metadata["NEW_META"] = "Testing"

        modified_number_of_filters = len(reader.filters)

        # note this is a good test because the original file had no filters,
        # and now we are adding one -- first time around it didn't get written
        # make sure this test doesn't get changed by anyone
        self.assertEquals(0, unmodified_number_of_filters)
        self.assertEquals(1, modified_number_of_filters)

        expected_output = "replace with real thing, this is here to make sure we read it"
        with open(os.path.join(os.path.dirname(__file__), "gatk_modified.vcf"), 'r') as f:
            lines = f.readlines()
            self.assertTrue(len(lines) >= 160)
            expected_output = "".join(lines)

        modified_out = StringIO()
        modified_writer = vcf.Writer(modified_out, reader)

        modified_out_str = modified_out.getvalue()

        print "******** mod out str **************"
        print modified_out_str
        print "**********************"

        # header lines does not include:
        # #CHROM POS ID REF ALT QUAL FILTER INFO FORMAT"
        self.assertTrue(unmodified_out_str.startswith("".join(reader._header_lines)))
        self.assertNotEquals(unmodified_out_str, modified_out_str)


        # expected output has data too
        self.assertTrue(expected_output.startswith(modified_out_str))

        # TODO also test data
        # TODO this test is far too long!!!

class TestRecord(unittest.TestCase):

    def test_num_calls(self):
        reader = vcf.Reader(fh('example-4.0.vcf'))
        for var in reader:
            num_calls = (var.num_hom_ref + var.num_hom_alt + \
                         var.num_het + var.num_unknown)
            self.assertEqual(len(var.samples), num_calls)

    def test_call_rate(self):
        reader = vcf.Reader(fh('example-4.0.vcf'))
        for var in reader:
            call_rate = var.call_rate
            if var.POS == 14370:
                self.assertEqual(3.0/3.0, call_rate)
            if var.POS == 17330:
                self.assertEqual(3.0/3.0, call_rate)
            if var.POS == 1110696:
                self.assertEqual(3.0/3.0, call_rate)
            if var.POS == 1230237:
                self.assertEqual(3.0/3.0, call_rate)
            elif var.POS == 1234567:
                self.assertEqual(2.0/3.0, call_rate)

    def test_aaf(self):
        reader = vcf.Reader(fh('example-4.0.vcf'))
        for var in reader:
            aaf = var.aaf
            if var.POS == 14370:
                self.assertEqual(3.0/6.0, aaf)
            if var.POS == 17330:
                self.assertEqual(1.0/6.0, aaf)
            if var.POS == 1110696:
                self.assertEqual(None, aaf)
            if var.POS == 1230237:
                self.assertEqual(0.0/6.0, aaf)
            elif var.POS == 1234567:
                self.assertEqual(None, aaf)

    def test_pi(self):
        reader = vcf.Reader(fh('example-4.0.vcf'))
        for var in reader:
            pi = var.nucl_diversity
            if var.POS == 14370:
                self.assertEqual(6.0/10.0, pi)
            if var.POS == 17330:
                self.assertEqual(1.0/3.0, pi)
            if var.POS == 1110696:
                self.assertEqual(None, pi)
            if var.POS == 1230237:
                self.assertEqual(0.0/6.0, pi)
            elif var.POS == 1234567:
                self.assertEqual(None, pi)

    def test_info(self):
        old_stderr = sys.stderr
        sys.stderr = StringIO()
        # this file has errors issued we don't want 
        # runners of the unit tests to be concerned with
        reader = vcf.Reader(fh('issue-info-4.1.vcf'))
        sys.stderr.close()
        sys.stderr = old_stderr
        for record in reader:
            info = record.INFO
            if record.POS == 14370:
                self.assertEqual("Test1", info["STR_1"])
                self.assertEqual(["FIRST", "SECOND"], info["STR_2"])
                self.assertEqual(["True"], info["STR_A"])
                self.assertEqual(["Something"], info["STR_G"])
                self.assertEqual(["0/0","0/1"], info["STR_VARIES"])
                self.assertEqual(3, info["NS"])
                self.assertEqual(14, info["DP"])
                self.assertEqual([0.5], info["AF"])
                self.assertNotIn("AA", info)
                self.assertTrue(info["DB"])
                self.assertTrue(info["H2"])
            elif record.POS == 17330:
                self.assertEqual("Test2", info["STR_1"])
                self.assertEqual(["Aa", "Bb"], info["STR_2"])
                self.assertEqual(["A"], record.ALT)
                self.assertEqual(["False"], info["STR_A"])
                self.assertEqual(["Something"], info["STR_G"])
                self.assertEqual(["0/0","0/1", "1/1"], info["STR_VARIES"])
                self.assertEqual(3, info["NS"])
                self.assertEqual(11, info["DP"])
                self.assertEqual([0.017], info["AF"])
                self.assertNotIn("AA", info)
                self.assertNotIn("DB", info)
                self.assertNotIn("H2", info)
            elif record.POS == 1110696:
                self.assertEqual("Test3", info["STR_1"])
                self.assertEqual(["Bb", "Cc"], info["STR_2"])
                self.assertTrue(isinstance(record.ALT, list))
                self.assertEqual(2, len(record.ALT))
                self.assertEqual(["False", "True"], info["STR_A"])
                self.assertEqual(["Something_else"], info["STR_G"])
                self.assertNotIn("STR_VARIES", info)
                self.assertEqual(2, info["NS"])
                self.assertEqual(10, info["DP"])
                self.assertEqual([0.333, 0.667], info["AF"])
                self.assertEqual("T", info["AA"])
                self.assertTrue(info["DB"])
                self.assertNotIn("H2", info)
            elif record.POS == 1230237:
                self.assertEqual("Test4", info["STR_1"])
                self.assertEqual(["FIRST", "SECOND"], info["STR_2"])
                self.assertEqual([None], record.ALT)
                self.assertNotIn("STR_A", info)
                self.assertEqual(["Something"], info["STR_G"])
                self.assertEqual(["0/0"], info["STR_VARIES"])
                self.assertEqual(3, info["NS"])
                self.assertEqual(13, info["DP"])
                self.assertNotIn("AF", info)
                self.assertEqual("T", info["AA"])
                self.assertNotIn("DB", info)
                self.assertNotIn("H2", info)
            elif record.POS == 1234567:
                self.assertEqual("Test5", info["STR_1"])
                self.assertEqual(["Aa", "Bb"], info["STR_2"])
                self.assertEqual(["False", "True"], info["STR_A"])
                self.assertEqual(["Something"], info["STR_G"])
                self.assertEqual(["0/0", "0/1", "1/1"], info["STR_VARIES"])
                self.assertEqual(3, info["NS"])
                self.assertEqual(9, info["DP"])
                self.assertNotIn("AF", info)
                self.assertEqual("G", info["AA"])
                self.assertNotIn("DB", info)
                self.assertNotIn("H2", info)
            else:
                self.assertFalse("Unexpected position: " + str(record.POS))

class TestCall(unittest.TestCase):

    def test_phased(self):
        reader = vcf.Reader(fh('example-4.0.vcf'))
        for var in reader:
            phases = [s.phased for s in var.samples]
            if var.POS == 14370:
                self.assertEqual([True, True, False], phases)
            if var.POS == 17330:
                self.assertEqual([True, True, False], phases)
            if var.POS == 1110696:
                self.assertEqual([True, True, False], phases)
            if var.POS == 1230237:
                self.assertEqual([True, True, False], phases)
            elif var.POS == 1234567:
                self.assertEqual([False, False, False], phases)

    def test_gt_bases(self):
        reader = vcf.Reader(fh('example-4.0.vcf'))
        for var in reader:
            gt_bases = [s.gt_bases for s in var.samples]
            if var.POS == 14370:
                self.assertEqual(['G|G', 'A|G', 'A/A'], gt_bases)
            elif var.POS == 17330:
                self.assertEqual(['T|T', 'T|A', 'T/T'], gt_bases)
            elif var.POS == 1110696:
                self.assertEqual(['G|T', 'T|G', 'T/T'], gt_bases)
            elif var.POS == 1230237:
                self.assertEqual(['T|T', 'T|T', 'T/T'], gt_bases)
            elif var.POS == 1234567:
                self.assertEqual([None, 'GTCT/GTACT', 'G/G'], gt_bases)

    def test_gt_types(self):
        reader = vcf.Reader(fh('example-4.0.vcf'))
        for var in reader:
            for s in var:
                print s.data
            gt_types = [s.gt_type for s in var.samples]
            if var.POS == 14370:
                self.assertEqual([0,1,2], gt_types)
            elif var.POS == 17330:
                self.assertEqual([0,1,0], gt_types)
            elif var.POS == 1110696:
                self.assertEqual([1,1,2], gt_types)
            elif var.POS == 1230237:
                self.assertEqual([0,0,0], gt_types)
            elif var.POS == 1234567:
                self.assertEqual([None,1,2], gt_types)

class TestTabix(unittest.TestCase):

    def setUp(self):
        self.reader = vcf.Reader(fh('tb.vcf.gz'))

        self.run = vcf.parser.pysam is not None


    def testFetchRange(self):
        if not self.run:
            return
        lines = list(self.reader.fetch('20', 14370, 14370))
        self.assertEquals(len(lines), 1)
        self.assertEqual(lines[0].POS, 14370)

        lines = list(self.reader.fetch('20', 14370, 17330))
        self.assertEquals(len(lines), 2)
        self.assertEqual(lines[0].POS, 14370)
        self.assertEqual(lines[1].POS, 17330)


        lines = list(self.reader.fetch('20', 1110695, 1234567))
        self.assertEquals(len(lines), 3)

    def testFetchSite(self):
        if not self.run:
            return
        site = self.reader.fetch('20', 14370)
        assert site.POS == 14370

        site = self.reader.fetch('20', 14369)
        assert site is None




class TestOpenMethods(unittest.TestCase):

    samples = 'NA00001 NA00002 NA00003'.split()

    def testOpenFilehandle(self):
        r = vcf.Reader(fh('example-4.0.vcf'))
        self.assertEqual(self.samples, r.samples)
        self.assertEqual('example-4.0.vcf', os.path.split(r.filename)[1])

    def testOpenFilename(self):
        r = vcf.Reader(filename='test/example-4.0.vcf')
        self.assertEqual(self.samples, r.samples)

    def testOpenFilehandleGzipped(self):
        r = vcf.Reader(fh('tb.vcf.gz'))
        self.assertEqual(self.samples, r.samples)

    def testOpenFilenameGzipped(self):
        r = vcf.Reader(filename='test/tb.vcf.gz')
        self.assertEqual(self.samples, r.samples)


class TestFilter(unittest.TestCase):


    def testApplyFilter(self):
        s, out = commands.getstatusoutput('python scripts/vcf_filter.py --site-quality 30 test/example-4.0.vcf sq')
        #print out
        assert s == 0
        buf = StringIO()
        buf.write(out)
        buf.seek(0)

        print buf.getvalue()
        reader = vcf.Reader(buf)

        # check filter got into output file
        assert 'sq30' in reader.filters

        print reader.filters

        # check sites were filtered
        n = 0
        for r in reader:
            if r.QUAL < 30:
                assert 'sq30' in r.FILTER
                n += 1
            else:
                # "PASS" is changed to None
                assert not r.FILTER or 'sq30' not in r.FILTER
        assert n == 2


    def testApplyMultipleFilters(self):
        s, out = commands.getstatusoutput('python scripts/vcf_filter.py --site-quality 30 '
        '--genotype-quality 50 test/example-4.0.vcf sq mgq')
        assert s == 0
        #print out
        buf = StringIO()
        buf.write(out)
        buf.seek(0)
        reader = vcf.Reader(buf)

        print reader.filters

        assert 'mgq50' in reader.filters
        assert 'sq30' in reader.filters

        
class TestRegression(unittest.TestCase):

    def test_issue_16(self):
        reader = vcf.Reader(fh('issue-16.vcf'))
        assert reader.next().QUAL == None

    def test_null_mono(self):
        # null qualities were written as blank, causing subsequent parse to fail
        print os.path.abspath(os.path.join(os.path.dirname(__file__),  'null_genotype_mono.vcf') )
        p = vcf.Reader(fh('null_genotype_mono.vcf'))
        assert p.samples
        out = StringIO()
        writer = vcf.Writer(out, p)
        map(writer.write_record, p)
        out.seek(0)
        print out.getvalue()
        p2 = vcf.Reader(out)
        rec = p2.next()
        assert rec.samples


class TestUtils(unittest.TestCase):

    def test_walk(self):
        # easy case: all same sites
        reader1 = vcf.Reader(fh('example-4.0.vcf'))
        reader2 = vcf.Reader(fh('example-4.0.vcf'))
        reader3 = vcf.Reader(fh('example-4.0.vcf'))

        n = 0
        for x in utils.walk_together(reader1, reader2, reader3):
            assert len(x) == 3
            assert (x[0] == x[1]) and (x[1] == x[2])
            n+= 1
        assert n == 5

        # artificial case 2 from the left, 2 from the right, 2 together, 1 from the right, 1 from the left

        expected = 'llrrttrl'
        reader1 = vcf.Reader(fh('walk_left.vcf'))
        reader2 = vcf.Reader(fh('example-4.0.vcf'))

        for ex, recs in zip(expected, utils.walk_together(reader1, reader2)):

            if ex == 'l':
                assert recs[0] is not None
                assert recs[1] is None
            if ex == 'r':
                assert recs[1] is not None
                assert recs[0] is None
            if ex == 't':
                assert recs[0] is not None
                assert recs[1] is not None





suite.addTests(unittest.TestLoader().loadTestsFromTestCase(TestGatkOutput))
suite.addTests(unittest.TestLoader().loadTestsFromTestCase(TestFreebayesOutput))
suite.addTests(unittest.TestLoader().loadTestsFromTestCase(TestSamtoolsOutput))
suite.addTests(unittest.TestLoader().loadTestsFromTestCase(TestReader))
suite.addTests(unittest.TestLoader().loadTestsFromTestCase(TestWriter))
suite.addTests(unittest.TestLoader().loadTestsFromTestCase(TestTabix))
suite.addTests(unittest.TestLoader().loadTestsFromTestCase(TestOpenMethods))
suite.addTests(unittest.TestLoader().loadTestsFromTestCase(TestFilter))
suite.addTests(unittest.TestLoader().loadTestsFromTestCase(Test1kg))
suite.addTests(unittest.TestLoader().loadTestsFromTestCase(TestRecord))
suite.addTests(unittest.TestLoader().loadTestsFromTestCase(TestCall))
suite.addTests(unittest.TestLoader().loadTestsFromTestCase(TestRegression))
