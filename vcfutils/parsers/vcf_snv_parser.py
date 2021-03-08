import vcfutils.utils.helpers as helpers
import pandas as pd
import vcf
import os
import itertools
from numpy import NaN
import gzip

class VcfUndeclaredSampleError(Exception):
    pass

class SNV_vcf():
    def __init__(self, filepath, filters):
        self.filepath = filepath
        self.reader = helpers._get_reader(filepath)
        self.filters = filters
        self.normal, self.tumor = self._read_sample_names(filepath)


    def _read_sample_names(self, filepath):
        header = helpers._get_header(gzip.open(filepath, "rt"))
        tumor_sample_line = list(filter(lambda line: 'tumor_sample' in line, header))
        normal_sample_line = list(filter(lambda line: 'normal_sample' in line, header))
        assert len(normal_sample_line) == 1
        normal = normal_sample_line[0].split("=")[-1]
        normal = normal.strip().strip("\n")
        tumor = None

        if len(tumor_sample_line) != 0:
            assert len(tumor_sample_line) == 1
            tumor = tumor_sample_line[0].split("=")[-1]
            tumor = tumor.strip().strip("\n")

        return normal, tumor


    def parse_main_cols(self, record):
        data = {"main_cols":{}, "samples":{}, "info":{}}
        data["main_cols"]= {
            'chrom': record.CHROM,
            'pos': record.POS,
            'ref': record.REF,
            'alt': record.ALT[:],
            'qual': record.QUAL,
            'filter': record.FILTER
        }
        data["main_cols"]['alt'] = ';'.join(map(str, data["main_cols"]['alt']))

        info = record.INFO
        data["info"] = info

        for sample in record.samples:

            sample_type = sample.sample
            sample_data = sample.data

            for k, v in sample_data._asdict().items():
                if sample_type == self.normal:
                    k += "_" + "NORMAL"
                elif sample_type == self.tumor:
                    k += "_" + "TUMOR"
                else:
                    raise VcfUndeclaredSampleError("unrecognized sample in vcf undefined in header")
                
                if isinstance(v, list):
                    v = ';'.join([str(val) for val in v])
                data["samples"][k] = v

        return data
    # def parse_record(self, record):
    #     data = self.parse_main_cols(record)
    #     # info = record.INFO

    #     # snpeff_annotations = self.parse_snpeff(self.snpeff_cols, info['ANN'], record.CHROM, record.POS)
    #     # ma_annotations = self.parse_mutation_assessor(self.ma_cols, info['MA'], record.CHROM, record.POS)

    #     # id_annotations = self.parse_list_annotations(info['DBSNP'], record.CHROM, record.POS, 'dbsnp')
    #     # id_annotations += self.parse_list_annotations(info['Cosmic'], record.CHROM, record.POS, 'cosmic')
    #     # ## TODO: I expected to see a 1000gen: False in the record. but the key is missing.
    #     # id_annotations += self.parse_flag_annotation(info.get('1000Gen'), record.CHROM, record.POS, '1000Gen')
    #     # id_annotations += self.parse_flag_annotation(info.get('LOW_MAPPABILITY'), record.CHROM, record.POS,
    #     #                                              'LOW_MAPPABILITY')
    #     return data

    def filter_records(self, record):

        for vcf_filter in self.filters:
            filter_name, relationship, value = vcf_filter

            if filter_name in record:
                if self.eval_expr(record[filter_name], relationship, value):
                    return True

            # for annotation in annotations:
            #     if filter_name == annotation['type']:
            #         if self.eval_expr(record[filter_name], relationship, value):
            #             return True

    def gather_records(self):
        for record in self.reader:
            data = self.parse_main_cols(record)

            if self.filters:
                if self.filter_records(data["main_cols"]):
                    continue
            yield data

    def as_data_frame(self):
        return pd.DataFrame(list(self.record_data))



class Mutect_vcf(SNV_vcf):
    def __init__(self, filepath, filters=None):
        super(Mutect_vcf, self).__init__(filepath, filters)
        self.record_data = []
        
    def parse(self):
        data = self.gather_records()
        self.record_data = [self.process_record(record) for record in data]

    def process_record(self, record):
        main_cols = record["main_cols"]
        info = record["info"]
        samples = record["samples"]
        
        main_cols["tumor_depth"] = samples["DP_TUMOR"]
        main_cols["tumor_alt_depth"] = samples["AD_TUMOR"].split(";")[1]
        main_cols["tumor_ref_depth"] = samples["AD_TUMOR"].split(";")[0]
        main_cols["normal_depth"] = samples["DP_NORMAL"]
        main_cols["normal_alt_depth"] = samples["AD_NORMAL"].split(";")[1]
        main_cols["normal_ref_depth"] = samples["AD_NORMAL"].split(";")[0]
        return main_cols

class Samtools_vcf(SNV_vcf):
    def __init__(self, filepath, filters=None):
        super(Samtools_vcf, self).__init__(filepath, filters=None)
        self.record_data = []
        
    def parse(self):
        data = self.gather_records()
        self.record_data = [self.process_record(record) for record in data]

    def process_record(self, record):
        main_cols = record["main_cols"]
        info = record["info"]

        main_cols["tumor_depth"] = NaN
        main_cols["tumor_alt_depth"] = NaN
        main_cols["tumor_ref_depth"] = NaN
        main_cols["normal_depth"] = info["DP"]
        main_cols["normal_alt_depth"] = NaN
        main_cols["normal_ref_depth"] = NaN
        return main_cols


class Freebayes_vcf(SNV_vcf):
    def __init__(self, filepath, filters=None):
        super(Freebayes_vcf, self).__init__(filepath, filters)
        self.record_data = []
        
    def parse(self):
        data = self.gather_records()
        self.record_data = [self.process_record(record) for record in data]

    def process_record(self, record):
        main_cols = record["main_cols"]
        samples = record["samples"]

        main_cols["tumor_depth"] = NaN
        main_cols["tumor_alt_depth"] = NaN
        main_cols["tumor_ref_depth"] = NaN
        main_cols["normal_depth"] = samples["DP_NORMAL"]
        main_cols["normal_alt_depth"] = samples["AO_NORMAL"]
        main_cols["normal_ref_depth"] = samples["RO_NORMAL"]
        return main_cols

class Bcftools_vcf(SNV_vcf):
    def __init__(self, filepath, filters=None):
        super(Bcftools_vcf, self).__init__(filepath, filters)

class Rtg_vcf(SNV_vcf):
    def __init__(self, filepath, filters=None):
        super(Rtg_vcf, self).__init__(filepath, filters)
        self.record_data = []
        
    def parse(self):
        data = self.gather_records()
        self.record_data = [self.process_record(record) for record in data]

    def process_record(self, record):
        main_cols = record["main_cols"]
        samples = record["samples"]

        main_cols["tumor_depth"] = NaN
        main_cols["tumor_alt_depth"] = NaN
        main_cols["tumor_ref_depth"] = NaN
        main_cols["normal_depth"] = samples["DP_NORMAL"]
        if isinstance(samples["AD_NORMAL"],int):
            main_cols["normal_alt_depth"] = NaN
            main_cols["normal_ref_depth"] = samples["AD_NORMAL"]
        else:
            main_cols["normal_alt_depth"] = samples["AD_NORMAL"].split(";")[1]
            main_cols["normal_ref_depth"] = samples["AD_NORMAL"].split(";")[0]

        return main_cols


class Strelka_vcf(SNV_vcf):
    def __init__(self, filepath, filters=None):
        super(Strelka_vcf, self).__init__(filepath, filters)
        self.record_data = []
        
    def parse(self):
        data = self.gather_records()
        self.record_data = [self.process_record(record) for record in data]

    def process_record(self, record):
        main_cols = record["main_cols"]
        samples = record["samples"]
        main_cols["tumor_depth"] = samples["DP_TUMOR"]
        main_cols["tumor_alt_depth"] = NaN
        main_cols["tumor_ref_depth"] = NaN
        main_cols["normal_depth"] = samples["DP_NORMAL"]
        main_cols["normal_alt_depth"] = NaN
        main_cols["normal_ref_depth"] = NaN

        return main_cols


class Museq_vcf(SNV_vcf):
    def __init__(self, filepath, filters):
        super(Museq_vcf, self).__init__(filepath, filters)
        self.record_data = []

    def parse(self):
        data = self.gather_records()
        self.record_data = [self.process_record(record) for record in data]

    def process_record(self, record):
        main_cols = record["main_cols"]
        samples = record["samples"]
        if self.tumor!=None:
            main_cols["tumor_depth"] = samples["DP_TUMOR"]
            main_cols["tumor_alt_depth"] = samples["AC_TUMOR"]
            main_cols["tumor_ref_depth"] = samples["RC_TUMOR"]
        else:
            main_cols["tumor_depth"] = NaN
            main_cols["tumor_alt_depth"] = NaN
            main_cols["tumor_ref_depth"] = NaN
        main_cols["normal_depth"] = samples["DP_NORMAL"]
        main_cols["normal_alt_depth"] = samples["AC_NORMAL"]
        main_cols["normal_ref_depth"] = samples["RC_NORMAL"]

        return main_cols
