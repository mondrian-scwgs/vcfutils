import vcfutils.utils.helpers as helpers
import pandas as pd
import vcf
import os
import itertools
from numpy import NaN
import gzip
from csverve import csverve 
import yaml
from typing import List, Dict, Tuple, Union, Iterable



class VcfUndeclaredSampleError(Exception):
    pass


class SNV_vcf():
    def __init__(self, filepath: str) -> None:
        """
        Sv vcf files (lumpy, svaba, gridss)
        @param filepath: Path of vcf.
        """
        self.filepath: str = filepath
        self.reader: vcf.Reader = helpers._get_reader(filepath)
        self.normal: str
        self.tumor: str
        self.normal, self.tumor = self._read_sample_names(filepath)

    def _read_sample_names(self, filepath: str) -> Tuple[str, str]:
        """
        read sample names from vcf
        @param filepath: Path of vcf.
        @return tumor/normal sample names
        """
        header: list[str] = helpers._get_header(gzip.open(filepath, "rt"))
        tumor_sample_line: str = list(filter(lambda line: 'tumor_sample' in line, header))
        normal_sample_line: str = list(filter(lambda line: 'normal_sample' in line, header))
        assert len(normal_sample_line) == 1
        normal: str = normal_sample_line[0].split("=")[-1]
        normal: str = normal.strip().strip("\n")
        tumor: str = None
        if len(tumor_sample_line) != 0:
            assert len(tumor_sample_line) == 1
            tumor: str = tumor_sample_line[0].split("=")[-1]
            tumor: str = tumor.strip().strip("\n")
        return normal, tumor

    def parse_main_cols(self,
        record: vcf.model._Record
    ) -> Iterable[Dict[str, Dict[str, Union[str, int, float]]]]:
        """
        parse commmon snv vcf columns from vcf record
        @param record: vcf record
        @return parsed record
        """
        data: Dict[str, Dict[str, Union[str, int, float]]]
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

        info: Dict[str, Union[str, int, float]] = record.INFO
        data["info"] = info

        for sample in record.samples:

            sample_type: str = sample.sample
            sample_data: Dict[str, Union[str, int, float]] = sample.data

            for k, v in sample_data._asdict().items():
                if sample_type == self.normal:
                    k += "_" + "NORMAL"
                elif sample_type == self.tumor:
                    k += "_" + "TUMOR"
                else:
                    raise VcfUndeclaredSampleError("unrecognized sample in vcf undefined in header")
                
                if isinstance(v, list):
                    v: str = ';'.join([str(val) for val in v])
                data["samples"][k] = v

        return data


    def gather_records(
            self
        ) ->  Iterable[Dict[str, Dict[str, Union[str, int, float]]]]:
        """
        gather parsed records from vcf
        @return parsed record
        """
        for record in self.reader:
            data:  Dict[str, Dict[str, Union[str, int, float]]]
            data = self.parse_main_cols(record)
            yield data


    def to_csv(self, output: str) -> None:
        """
        write parsed vcf data to a csv
        @param output: output file
        """
        if not output.endswith(".gz"):
            output += ".gz"
        print(itertools.islice(self.record_data, 100000-1))
        dataframes: Mapping[pd.DataFrame, Dict[str, Union[str, int, float]]]
        dataframes = map(pd.DataFrame, helpers._group_iterator(self.record_data))
        write_header = True
        for dataframe in dataframes:
            print("here")
            cc
            csverve.write_dataframe_to_csv_and_yaml(dataframe, output, 
                dataframe.dtypes, write_header=write_header
            )
            write_header=False
        yaml_file: str = output + ".yaml"
        metadata: Dict[str, str] = yaml.load(open(yaml_file))
        metadata["samples"] = {"tumor": self.tumor, "normal": self.normal}
        with open(yaml_file, 'wt') as f:
            yaml.safe_dump(metadata, f, default_flow_style=False)    
    

class Mutect_vcf(SNV_vcf):
    def __init__(self, filepath: str) -> None:
        """
        Mutect Snv vcf.
        @param filapth: filepath of vcf
        """
        super(Mutect_vcf, self).__init__(filepath)
        self.record_data: Iterable[Dict[str, Union[str, int, float]]]  = []
        
    def parse(self) -> None:
        """
        parse records from vcf.
        """
        records: Iterable[Dict[str, Union[str, int, float]]]
        records = iter(self.gather_records())
        self.record_data: Iterable[Dict[str, Union[str, int, float]]]
        self.record_data = map(self.process_record, records)

    def process_record(self, 
        record: Dict[str, Dict[str, Union[str, int, float]]]
    ) -> Dict[str, Union[str, int, float]]:
        """
        add tumor read data to vcf records.
        @param: record: vcf record
        @return: record with additional data
        """
        main_cols: Dict[str, Union[str, int, float]]
        info: Dict[str, Union[str, int, float]]
        samples: Dict[str, Union[str, int, float]]
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
    def __init__(self, filepath: str) -> None:
        """
        Samtools SNV vcf
        @param filepath: filepath to vcf
        """
        super(Samtools_vcf, self).__init__(filepath)
        self.record_data: Iterable[Dict[str, Union[str, int, float]]] = []
        
    def parse(self) -> None:
        """
        parse records from samtools vcf
        """
        records: Iterable[Dict[str, Union[str, int, float]]]
        records = iter(self.gather_records())
        self.record_data: Iterable[Dict[str, Union[str, int, float]]]
        self.record_data = map(self.process_record, records)

    def process_record(self, 
        record:Dict[str, Dict[str, Union[str, int, float]]] 
    ) -> Dict[str, Union[str, int, float]]:
        """
        add tumor read data to vcf records.
        @param: record: vcf record
        @return: record with additional data
        """
        main_cols: Dict[str, Union[str, int, float]]
        info: Dict[str, Union[str, int, float]]
        samples: Dict[str, Union[str, int, float]]
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
    def __init__(self, filepath: str) -> None:
        """
        Freebayes snv vcf.
        @param: filepath: vcf filepath
        """
        super(Freebayes_vcf, self).__init__(filepath)
        self.record_data: Iterable[Dict[str, Union[str, int, float]]]  = []
        
    def parse(self) -> None:
        """
        parse records from samtools vcf
        """
        records: Iterable[Dict[str, Union[str, int, float]]]
        records = iter(self.gather_records())
        self.record_data: Iterable[Dict[str, Union[str, int, float]]]
        self.record_data = map(self.process_record, records)


    def process_record(self, record):
        """
        add tumor read data to vcf records.
        @param: record: vcf record
        @return: record with additional data
        """
        main_cols: Dict[str, Union[str, int, float]]
        info: Dict[str, Union[str, int, float]]
        samples: Dict[str, Union[str, int, float]]
        main_cols = record["main_cols"]
        samples = record["samples"]

        main_cols["tumor_depth"] = NaN
        main_cols["tumor_alt_depth"] = NaN
        main_cols["tumor_ref_depth"] = NaN
        main_cols["normal_depth"] = samples["DP_NORMAL"]
        main_cols["normal_alt_depth"] = samples["AO_NORMAL"]
        main_cols["normal_ref_depth"] = samples["RO_NORMAL"]
        return main_cols


class Rtg_vcf(SNV_vcf):
    def __init__(self, filepath):
        """
        Rtg SNV vcf.
        """
        super(Rtg_vcf, self).__init__(filepath)
        self.record_data: Iterable[Dict[str, Union[str, int, float]]] = []
        
    def parse(self) -> None:
        """
        parse records from samtools vcf
        """
        records: Iterable[Dict[str, Union[str, int, float]]]
        records = iter(self.gather_records())
        self.record_data: Iterable[Dict[str, Union[str, int, float]]]
        self.record_data = map(self.process_record, records)


    def process_record(self, 
        record:Dict[str, Dict[str, Union[str, int, float]]] 
    ) -> Dict[str, Union[str, int, float]]:
        """
        add tumor read data to vcf records.
        @param: record: vcf record
        @return: record with additional data
        """
        main_cols: Dict[str, Union[str, int, float]]
        info: Dict[str, Union[str, int, float]]
        samples: Dict[str, Union[str, int, float]]
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
    def __init__(self, filepath):
        """
        Strelka snv vcf.
        """
        super(Strelka_vcf, self).__init__(filepath)
        self.record_data: Iterable[Dict[str, Union[str, int, float]]] = []
        
    def parse(self) -> None:
        """
        parse records from samtools vcf
        """
        records: Iterable[Dict[str, Union[str, int, float]]]
        records = iter(self.gather_records())
        self.record_data: Iterable[Dict[str, Union[str, int, float]]]
        self.record_data = map(self.process_record, records)


    def process_record(self, 
        record:Dict[str, Dict[str, Union[str, int, float]]] 
    ) -> Dict[str, Union[str, int, float]]:
        """
        add tumor read data to vcf records.
        @param: record: vcf record
        @return: record with additional data
        """
        main_cols: Dict[str, Union[str, int, float]]
        info: Dict[str, Union[str, int, float]]
        samples: Dict[str, Union[str, int, float]]
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
    def __init__(self, filepath):
        """
        Museq SNV vcf.
        """
        super(Museq_vcf, self).__init__(filepath)
        self.record_data: Iterable[Dict[str, Union[str, int, float]]] = []

    def parse(self) -> None:
        """
        parse records from samtools vcf
        """
        records: Iterable[Dict[str, Union[str, int, float]]]
        records = iter(self.gather_records())
        self.record_data: Iterable[Dict[str, Union[str, int, float]]]
        self.record_data = map(self.process_record, records)

    def process_record(self, 
        record:Dict[str, Dict[str, Union[str, int, float]]] 
    ) -> Dict[str, Union[str, int, float]]:
        """
        add tumor read data to vcf records.
        @param: record: vcf record
        @return: record with additional data
        """
        main_cols: Dict[str, Union[str, int, float]]
        info: Dict[str, Union[str, int, float]]
        samples: Dict[str, Union[str, int, float]]
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
