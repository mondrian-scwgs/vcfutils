import vcfutils.helpers as helpers
# import helpers
import pandas as pd
import vcf
import os
import itertools

# from pandas.api.types import CategoricalDtype

ANNOTATIONS = ['ann', 'ma', 'dbsnp', 'cosmic', 'lof', 'nmd', '1000gen', 'low_mappability']

########################
# file-->file function #
# no parsing ###########
########################
def concatenate_vcf(infiles, outfile):
    '''
    Concatenate VCF files
    :param infiles: list of input VCF files to be concatenated
    :param outfile: output VCF file
    '''

    with open(outfile, 'w') as ofile:
        header = None

        for ifile in infiles:

            if os.path.getsize(ifile) == 0:

                warnings.warn('input file {} is empty'.format(ifile))
                continue

            with open(ifile) as f:

                if not header:
                    header = helpers._get_header(f)

                    for line in header:
                        ofile.write(line)
                else:
                    if not helpers._get_header(f) == header:
                        warnings.warn('merging vcf files with mismatching headers')

                for l in f:
                    ofile.write(l)


def sort_vcf(infile, outfile, docker_image=None):
    cmd = ['cat', infile, '|', 'vcf-sort', '>', outfile]

    pypeliner.commandline.execute(*cmd, docker_image=docker_image)


def index_vcf(vcf_file, docker_image=None):
    """ Create a tabix index for a VCF file
    :param vcf_file: Path of VCF to create index for. Should compressed by bgzip.
    :param index_file: Path of index file.
    This is meant to be used from pypeliner so it does some name mangling to add .tmp to the index file.
    """

    pypeliner.commandline.execute('tabix', '-f', '-p', 'vcf', vcf_file, docker_image=docker_image)


def compress_vcf(in_file, out_file, docker_image=None):
    """ Compress a VCF file using bgzip.
    :param in_file: Path of uncompressed VCF file.
    :param out_file: Path were compressed VCF file will be written.
    """
    pypeliner.commandline.execute('bgzip', '-c', in_file, '>', out_file, docker_image=docker_image)


def split_vcf(in_file, out_dir, lines_per_file):
    """ Split a VCF file into smaller files.
    :param in_file: Path of VCF file to split.
    :param out_dir: output directory for split vcf files
    :param lines_per_file: Maximum number of lines to be written per file.
     """
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

    def line_group(line, line_idx=itertools.count()):
        return int(next(line_idx) / lines_per_file)

    reader = vcf.Reader(filename=in_file)

    for file_idx, records in itertools.groupby(reader, key=line_group):
        file_name = os.path.join(out_dir, str(file_idx) + "_split.vcf")

        with open(file_name, 'w') as out_fh:
            writer = vcf.Writer(out_fh, reader)

            for record in records:
                writer.write_record(record)

            writer.close()


class SV_vcf():
    def __init__(self, filepath, caller):
        self.filepath = filepath
        self.reader = helpers._get_reader(filepath)
        self.caller = caller


    @staticmethod

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

            for sample in record.samples:
                sample_name = sample.sample
                sample_data = sample.data
                for k, v in sample_data._asdict().items():
                    if isinstance(v, list):
                        v = ';'.join([str(val) for val in v])
                    k = '{}_{}'.format(sample_name, k)
                    data[k] = v

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
            'type': mate1['SVTYPE']
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

    def fetch(self):
        records = self._parse_vcf()
        records = self._group_bnds(records)
        records = self._filter_low_qual_calls(records)

        for record in records:
            if len(record) == 1 and self.caller == 'lumpy':
                yield self._process_lumpy_unmatched_record(record)
            else:
                yield self._process_bnd_call(record)

    def as_data_frame(self):
        data = [record for record in self.fetch()]

        data = pd.DataFrame(data)
        data['caller'] = self.caller

        data['breakpoint_id'] = data.index

        data['breakpoint_id'] = data['breakpoint_id'].astype(str) + '_' + data['caller']

        return data

    def to_csv(self, output):
        df = self.as_data_frame()
        df.to_csv(output, sep="\t", index=False)

class Lumpy_vcf(SV_vcf):
    def __init__(self, filepath):
        super(Lumpy_vcf, self).__init__(filepath, "lumpy")

class Gridss_vcf(SV_vcf):
    def __init__(self, filepath):
        super(Gridss_vcf, self).__init__(filepath, "gridss")

class Svaba_vcf(SV_vcf):
    def __init__(self, filepath):
        super(Svaba_vcf, self).__init__(filepath, "svaba")


class SNV_vcf():
    def __init__(self, filepath, caller):
        self.filepath = filepath
        self.reader = helpers._get_reader(filepath)
        self.caller = caller


    def parse_main_cols(self, record):
        data = {
            'chrom': record.CHROM,
            'pos': record.POS,
            'ref': record.REF,
            'alt': record.ALT,
            'qual': record.QUAL,
        }

        data['alt'] = ';'.join(map(str, data['alt']))

        info = record.INFO

        for k, v in info.items():
            if k.lower() in ANNOTATIONS:
                continue
            if isinstance(v, list):
                v = ';'.join(v)
            data[k] = v

        for sample in record.samples:

            sample_type = sample.sample
            sample_data = sample.data

            for k, v in sample_data._asdict().items():
                k += "_" + sample_type
                if isinstance(v, list):
                    v = ';'.join([str(val) for val in v])
                data[k] = v

        return data


    def parse_record(self, record):
        data = self.parse_main_cols(record)
        info = record.INFO

        # snpeff_annotations = self.parse_snpeff(self.snpeff_cols, info['ANN'], record.CHROM, record.POS)
        # ma_annotations = self.parse_mutation_assessor(self.ma_cols, info['MA'], record.CHROM, record.POS)

        # id_annotations = self.parse_list_annotations(info['DBSNP'], record.CHROM, record.POS, 'dbsnp')
        # id_annotations += self.parse_list_annotations(info['Cosmic'], record.CHROM, record.POS, 'cosmic')
        # ## TODO: I expected to see a 1000gen: False in the record. but the key is missing.
        # id_annotations += self.parse_flag_annotation(info.get('1000Gen'), record.CHROM, record.POS, '1000Gen')
        # id_annotations += self.parse_flag_annotation(info.get('LOW_MAPPABILITY'), record.CHROM, record.POS,
        #                                              'LOW_MAPPABILITY')
        return data


    def parse_vcf(self):
        for record in self.reader:

            data = self.parse_record(record)

            # if self.filter_records(data, id_annotations):
            #     continue

            yield data

    def as_data_frame(self):
        return pd.DataFrame(list(self.parse_vcf()))

    def parse_snpeff(self, snpeff_cols, snpeff_entries, chrom, pos):
        records = []
        for record in snpeff_entries:
            record = record.strip().split('|')

            record = dict(zip(snpeff_cols, record))

            record['chrom'] = chrom
            record['pos'] = pos

            records.append(record)

        return records

    def parse_mutation_assessor(self, cols, record, chrom, pos):
        if not record:
            return
        record = record.strip().split('|')

        record = dict(zip(cols, record))

        record['chrom'] = chrom
        record['pos'] = pos

        return record

    def parse_list_annotations(self, records, chrom, pos, label):
        if records == [None]:
            return []

        parsed_records = []
        for entry in records:
            parsed_records.append(
                {'chrom': chrom, 'pos': pos, 'value': entry, 'type': label}
            )

        return parsed_records

    def parse_flag_annotation(self, record, chrom, pos, label):
        if not record:
            return []

        return [{'chrom': chrom, 'pos': pos, 'value': record, 'type': label}]



class Mutect_vcf(SNV_vcf):
    def __init__(self, filepath):
        super(Mutect_vcf, self).__init__(filepath, "mutect")

class Samtools_vcf(SNV_vcf):
    def __init__(self, filepath):
        super(Samtools_vcf, self).__init__(filepath, "samtools")

class Freebayes_vcf(SNV_vcf):
    def __init__(self, filepath):
        super(Freebayes_vcf, self).__init__(filepath, "freebayes")

class Bcftools_vcf(SNV_vcf):
    def __init__(self, filepath):
        super(Bcftools_vcf, self).__init__(filepath, "bcftools")

class Rtg_vcf(SNV_vcf):
    def __init__(self, filepath):
        super(Rtg_vcf, self).__init__(filepath, "rtg")

class Strelka_vcf(SNV_vcf):
    def __init__(self, filepath):
        super(Strelka_vcf, self).__init__(filepath, "strelka")

class Museq_vcf(SNV_vcf):
    def __init__(self, filepath):
        super(Museq_vcf, self).__init__(filepath, "museq")

# class SNV_vcf():


mu = Museq_vcf("/juno/work/shah/tantalus/SC-3281/results/variant_calling/sample_SA1202LA/strelka_indel.vcf.gz")
print(mu.as_data_frame())

# split_vcf("/juno/work/shah/users/grewald/SVBENCH2/trimming_primary_only/wgs/output/breakpoints/SA1256PP/SA1256PP_lumpy.vcf", 500)