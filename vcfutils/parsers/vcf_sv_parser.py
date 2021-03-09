import vcfutils.utils.helpers as helpers
import pandas as pd
import vcf
import os
import itertools


class SV_vcf():
    def __init__(self, filepath, caller):
        self.filepath = filepath
        self.reader = helpers._get_reader(filepath)
        self.caller = caller


    def _parse_vcf(self):

        for record in self.reader:

            data = {
                'chrom': record.CHROM,
                'pos': record.POS,
                'ref': record.REF,
                'alt': record.ALT,
                'qual': record.QUAL,
                'id': record.ID,
                'filter': record.FILTER,
            }

            info = record.INFO

            for k, v in info.items():
                if isinstance(v, list):
                    v = ';'.join(map(str, v))
                data[k] = v
            
            assert len(record.samples)==2
            
            sample_type = "normal"
            for sample in record.samples:
                sample_name = sample.sample
                sample_data = sample.data
                for k, v in sample_data._asdict().items():
                    if isinstance(v, list):
                        v = ';'.join([str(val) for val in v])
                    k = '{}_{}'.format(sample_type, k)
                    data[k] = v
                sample_type="tumor"
            yield data


    @staticmethod
    def _group_bnds(calls):
        bnds = {}

        for record in calls:
            if record['SVTYPE'] == 'BND':
                if 'MATEID' not in record:
                    continue

                if record['MATEID'] in bnds:
                    
                    yield (record, bnds[record['MATEID']])
                    bnds.pop(record['MATEID'])
                else:
                    bnds[record['id']] = record
            else:
                yield record,

        assert len(bnds) == 0

    def _get_mates(self, records):
        ends_with_val = {
            'lumpy': '_1',
            'svaba': ':1',
            'gridss': 'h'
        }

        if records[0]['id'].endswith(ends_with_val[self.caller]):
            mate1, mate2 = records
        else:
            mate2, mate1 = records

        return mate1, mate2

    @staticmethod
    def _get_strand_from_alt(alt):
        """
        If the nucleotide comes first, then it is a "+" facing break at
        that site (see element "W" of the VCF4.2 specs, pg 13), otherwise
        it is "-" at that site. Then, if the bracket is ], then the
        partner breakpoint is "+", otherwise it is left facing
        :param alt:
        :type alt:
        :return:
        :rtype:
        """
        return alt[0].orientation

    def _get_strands(self, mate1, mate2):
        if self.caller == 'lumpy':
            strands_1 = mate1['STRANDS'].split(':')[0]
            strands_2 = mate2['STRANDS'].split(':')[0]
            assert strands_1 == strands_2[::-1]
            return strands_1[0], strands_2[0]
        else:
            strand_1 = self._get_strand_from_alt(mate1['alt'])
            strand_2 = self._get_strand_from_alt(mate2['alt'])

        return strand_1, strand_2

    @staticmethod
    def _process_lumpy_unmatched_record(record):
        record = record[0]

        strands = record['STRANDS'].split(':')[0]
        assert len(strands) == 2
        outdata = {
            'chromosome_1': record['chrom'],
            'position_1': record['pos'],
            'chromosome_2': record['chrom'],
            'position_2': record['END'],
            'strand_1': strands[0],
            'strand_2': strands[1],
            'type': record['SVTYPE']
        }

        return outdata

    def _process_bnd_call(self, record):

        assert len(record) == 2

        mate1, mate2 = self._get_mates(record)
        strand_1, strand_2 = self._get_strands(mate1, mate2)

        assert mate1['SVTYPE'] == mate2['SVTYPE']

        outdata = {
            'chromosome_1': mate1['chrom'],
            'position_1': mate1['pos'],
            'strand_1': strand_1,
            'chromosome_2': mate2['chrom'],
            'position_2': mate2['pos'],
            'strand_2': strand_2,
            'type': mate1['SVTYPE'],
            'filter': mate1["filter"],
            "quality": mate1["qual"]
        }

        return outdata

    def _filter_low_qual_calls(self, calls):

        for call in calls:
            if len(call) == 1 and self.caller == 'lumpy':

                if call[0]['filter'] and 'LOW_QUAL' in call[0]['filter']:
                    continue
            else:
                assert len(call) == 2

                if call[0]['filter'] and 'LOW_QUAL' in call[0]['filter'] and 'LOW_QUAL' in call[1]['filter']:
                    continue

            yield call

    def gather_records(self):
        records = self._parse_vcf()
        records = self._group_bnds(records)
        records = self._filter_low_qual_calls(records)
        return records

    def _classify_calls(self, data):

        data['type'] = None
        data.loc[(data['position_1'] < data['position_2']) & (data['strand_1'] == '+') & (data['strand_2'] == '-'), 'rearrangement_type'] = 'deletion'
        data.loc[(data['position_2'] < data['position_1']) & (data['strand_2'] == '+') & (data['strand_1'] == '-'), 'rearrangement_type'] = 'deletion'
        data.loc[(data['position_1'] < data['position_2']) & (data['strand_1'] == '-') & (data['strand_2'] == '+'), 'rearrangement_type'] = 'duplication'       
        data.loc[(data['position_2'] < data['position_1']) & (data['strand_2'] == '-') & (data['strand_1'] == '+'), 'rearrangement_type'] = 'duplication'       
        data.loc[(data['strand_1'] == data['strand_2']), 'rearrangement_type'] = 'inversion'
        data.loc[(data['chromosome_1'] != data['chromosome_2']), 'rearrangement_type'] = 'translocation'
        # then for size...
        data['length'] = (data['position_1'] - data['position_2']).abs()
        data['size_class'] = pd.cut(data['length'], [0, 1e4, 1e6, 1e10], labels=['S', 'M', 'L'])

        data.loc[(data['rearrangement_type'] == 'translocation'), 'length'] = float('inf')
        data.loc[(data['rearrangement_type'] == 'translocation'), 'size_class'] = 'L'

        return data


    def as_data_frame(self):
        data = pd.DataFrame(self.record_data)
        data['caller'] = self.caller
        data['breakpoint_id'] = data.index
        data['breakpoint_id'] = data['breakpoint_id'].astype(str) + '_' + data['caller']
        return data

    def to_csv(self, output):
        df = self.as_data_frame()
        df = self._classify_calls(df)
        df.to_csv(output, sep="\t", index=False)


class Lumpy_vcf(SV_vcf):
    def __init__(self, filepath):
        super(Lumpy_vcf, self).__init__(filepath, "lumpy")
        self.record_data = []
        
    def parse(self):
        records = self.gather_records()
        self.record_data = [self.process_record(record) for record in records]

    def process_record(self, record):
        mate1 = record[0]

        if len(record) == 1:
            processed = self._process_lumpy_unmatched_record(record)
        else:
            processed = self._process_bnd_call(record)
            
        processed["tumor_depth"] = mate1["tumor_SU"]
        processed["tumor_split_reads"] = mate1["tumor_SR"]
        processed["tumor_discordant_reads"] =  mate1["tumor_PE"]
        processed["normal_depth"] = mate1["normal_SU"]
        processed["normal_split_reads"] = mate1["normal_SR"]
        processed["normal_discordant_reads"] = mate1["normal_PE"]

        return processed


class Gridss_vcf(SV_vcf):
    def __init__(self, filepath):
        super(Gridss_vcf, self).__init__(filepath, "gridss")
        self.record_data = []
        
    def parse(self):
        records = self.gather_records()
        self.record_data =  [self.process_record(record) for record in records]

    def process_record(self, record):
        mate1 = record[0]
        processed = self._process_bnd_call(record)
        processed["tumor_depth"] = mate1["tumor_VF"]
        processed["tumor_split_reads"] = mate1["tumor_SR"]
        processed["tumor_discordant_reads"] =  mate1["tumor_RP"]
        processed["normal_depth"] = mate1["normal_VF"]
        processed["normal_split_reads"] = mate1["normal_SR"]
        processed["normal_discordant_reads"] = mate1["normal_RP"]

        return processed


class Svaba_vcf(SV_vcf):
    def __init__(self, filepath):
        super(Svaba_vcf, self).__init__(filepath, "svaba")
        self.record_data = []
        
    def parse(self):
        records = self.gather_records()
        self.record_data =  [self.process_record(record) for record in records]

    def process_record(self, record):
        mate1 = record[0]
        processed = self._process_bnd_call(record)
        processed["tumor_depth"] = mate1["tumor_AD"]
        processed["tumor_split_reads"] = mate1["tumor_SR"]
        processed["tumor_discordant_reads"] =  mate1["tumor_DR"]
        processed["normal_depth"] = mate1["normal_AD"]
        processed["normal_split_reads"] = mate1["normal_SR"]
        processed["normal_discordant_reads"] = mate1["normal_DR"]

        return processed