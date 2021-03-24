from vcfutils.parsers import vcf_sv_parser
import pytest
import os
import yaml
import pandas as pd
from numpy import isnan
from itertools import tee
import vcf

def _raises_correct_error(function, *args,expected_error,
                        **kwargs):
    raised = False
    try:
        function(*args, **kwargs)
    except Exception as e:
        if type(e) == expected_error:
            raised = True
        else:
            print("raised wrong error: raised: {}, expected: {}"
                   .format(type(e), expected_error))
    finally:
        return raised

@pytest.fixture
def lumpy():
    file = os.path.join(os.path.dirname(os.path.realpath(__file__)), "testdata", "lumpy.vcf.gz")
    assert os.path.exists(file)
    return file

@pytest.fixture
def gridss():
    file = os.path.join(os.path.dirname(os.path.realpath(__file__)), "testdata", "gridss.vcf.gz")
    assert os.path.exists(file)
    return file

@pytest.fixture
def svaba():
    file = os.path.join(os.path.dirname(os.path.realpath(__file__)), "testdata", "svaba.vcf.gz")
    assert os.path.exists(file)
    return file

@pytest.fixture
def bnd():
    mateid = "_1"
    id = "_2"
    return [{"SVTYPE": "BND", "MATEID": mateid, "id": id, "STRANDS": "+-:5"}, 
        {"SVTYPE": "BND", "MATEID": id, "id": mateid, "STRANDS": "-+:5"}]


@pytest.fixture
def test_output():
    o = "/juno/work/shah/abramsd/CODE/vcfutils/testoutput.csv.gz"
    if os.path.exists(o):
        os.path.remove(o)
    return o

def _check_records(a, b):
    assert a["chrom"] == b.CHROM
    assert a["pos"] == b.POS
    assert a["ref"] == b.REF
    assert a["alt"] == b.ALT
    assert a["qual"] == b.QUAL
    assert a["filter"] == b.FILTER


def _check_yaml_formatting(y_file, parser):
    yaml_data = yaml.load(open(y_file))
    print(yaml_data)


def test_parse_lumpy(lumpy):
    p = vcf_sv_parser.Lumpy_vcf(lumpy)
    record = next(vcf.Reader(filename=lumpy))
    parsed_record = next(p._parse_vcf())
    _check_records(parsed_record, record)

def test_parse_gridss(gridss):
    p = vcf_sv_parser.Gridss_vcf(gridss)
    record = next(vcf.Reader(filename=gridss))
    parsed_record = next(p._parse_vcf())
    _check_records(parsed_record, record)

def test_parse_svaba(svaba):
    p = vcf_sv_parser.Svaba_vcf(svaba)
    record = next(vcf.Reader(filename=svaba))
    parsed_record = next(p._parse_vcf())
    _check_records(parsed_record, record)


def test_bnd_matched(lumpy, bnd):
    matched = next(vcf_sv_parser.SV_vcf(lumpy, caller="lumpy")._group_bnds(bnd)) 
    assert isinstance(matched, tuple)
    
    
# def test_bnd_unmatched(lumpy, bnd):
#     bnd[1]["MATEID"] = -1
#     matched = next(vcf_sv_parser.SV_vcf(lumpy, caller="lumpy")._group_bnds(bnd) )
#     print(type(matched), matched)
#     cc
#     assert isinstance(matched, list)


def test_get_mates_lumpy(lumpy, bnd):
    bnd[0]["id"] = "_1"
    matched = vcf_sv_parser.SV_vcf(lumpy, caller="lumpy")._get_mates(bnd)
    assert bnd[0] == matched[0]
    assert bnd[1] == matched[1] 

def test_get_mates_reversed_lumpy(lumpy, bnd):
    bnd[0]["id"] = "_2"
    matched = vcf_sv_parser.SV_vcf(lumpy, caller="lumpy")._get_mates(bnd)
    assert bnd[0] == matched[1]
    assert bnd[1] == matched[0] 

def test_get_mates_gridss(gridss, bnd):
    bnd[0]["id"] = "_h"
    matched = vcf_sv_parser.SV_vcf(gridss, caller="gridss")._get_mates(bnd)
    assert bnd[0] == matched[0]
    assert bnd[1] == matched[1] 

def test_get_mates_reversed_gridss(gridss, bnd):
    bnd[0]["id"] = "_2"
    matched = vcf_sv_parser.SV_vcf(gridss, caller="gridss")._get_mates(bnd)
    assert bnd[0] == matched[1]
    assert bnd[1] == matched[0] 

def test_get_mates_svaba(svaba, bnd):
    bnd[0]["id"] = ":1"
    matched = vcf_sv_parser.SV_vcf(svaba, caller="svaba")._get_mates(bnd)
    assert bnd[0] == matched[0]
    assert bnd[1] == matched[1] 

def test_get_mates_reversed_svaba(svaba, bnd):
    bnd[0]["id"] = ":2"
    matched = vcf_sv_parser.SV_vcf(svaba, caller="svaba")._get_mates(bnd)
    assert bnd[0] == matched[1]
    assert bnd[1] == matched[0] 

def test_get_strands_lumpy(lumpy, bnd):
    bnd[0]["id"] = "_1"
    p = vcf_sv_parser.SV_vcf(lumpy, caller="lumpy")
    mate1, mate2 = p._get_mates(bnd)
    strand1, strand2 = p._get_strands(mate1, mate2)
    assert strand1 == "+"
    assert strand2 == "-"


def test_get_strands_lumpy_unmatched(lumpy, bnd):
    bnd[0]["id"] = "_1"
    bnd[0]["STRANDS"] = bnd[1]["STRANDS"]
    p = vcf_sv_parser.SV_vcf(lumpy, caller="lumpy")
    mate1, mate2 = p._get_mates(bnd)
    assert _raises_correct_error(p._get_strands, mate1, mate2,
                              expected_error=AssertionError)


def test_filter_low_qual_calls_lumpy_filtered(lumpy):
    p = vcf_sv_parser.SV_vcf(lumpy, caller="lumpy")
    records = p._parse_vcf()
    record = [list(next(p._group_bnds(records)))]
    record[0][0]["filter"] = "LOW_QUAL"
    filt_records = list(p._filter_low_qual_calls(record))
    assert filt_records == []

def test_filter_low_qual_calls_lumpy_unfiltered(lumpy):
    p = vcf_sv_parser.SV_vcf(lumpy, caller="lumpy")
    records = p._parse_vcf()
    record = [list(next(p._group_bnds(records)))]
    filt_records = list(p._filter_low_qual_calls(record))
    assert record[0][0] == filt_records[0][0]

def test_filter_low_qual_calls_non_svaba_filtered(svaba):
    p = vcf_sv_parser.SV_vcf(svaba, caller="svaba")
    records = p._parse_vcf()
    record = [list(next(p._group_bnds(records)))]
    record[0][0]["filter"] = "LOW_QUAL"
    record[0][1]["filter"] = "LOW_QUAL"
    filt_records = list(p._filter_low_qual_calls(record))
    assert filt_records == []

def test_filter_low_qual_calls_svaba_unfiltered(svaba):
    p = vcf_sv_parser.SV_vcf(svaba, caller="svaba")
    records = p._parse_vcf()
    record = [list(next(p._group_bnds(records)))]
    filt_records = list(p._filter_low_qual_calls(record))
    assert record[0][0] == filt_records[0][0]

def test_as_data_frame(lumpy):
    p = vcf_sv_parser.Lumpy_vcf(lumpy)
    records = p.gather_records()
    p.record_data = [p.process_record(record) for record in records]
    df = p.as_data_frame()
    assert df.position_1.tolist() == [r["position_1"] for r  in p.record_data]
 
def test_lumpy_parse(lumpy):
    p = vcf_sv_parser.Lumpy_vcf(lumpy)
    records, dup = tee(p.gather_records())
    processed_records = [p.process_record(record) for record in records]
    for record, processed in zip(dup, processed_records):
        record = record[0]
        processed["tumor_depth"] == record["tumor_SU"]
        processed["tumor_split_reads"] == record["tumor_SR"]
        processed["tumor_discordant_reads"] ==  record["tumor_PE"]
        processed["normal_depth"] == record["normal_SU"]
        processed["normal_split_reads"] == record["normal_SR"]
        processed["normal_discordant_reads"] == record["normal_PE"]
 
def test_svaba_parse(svaba):
    p = vcf_sv_parser.Svaba_vcf(svaba)
    records, dup = tee(p.gather_records())
    processed_records = [p.process_record(record) for record in records]
    for record, processed in zip(dup, processed_records):
        record = record[0]
        assert processed["tumor_depth"] == record["tumor_AD"]
        assert processed["tumor_split_reads"] == record["tumor_SR"]
        assert processed["tumor_discordant_reads"] ==  record["tumor_DR"]
        assert processed["normal_depth"] == record["normal_AD"]
        assert processed["normal_split_reads"] == record["normal_SR"]
        assert processed["normal_discordant_reads"] == record["normal_DR"]

def test_gridss_parse(gridss):
    p = vcf_sv_parser.Gridss_vcf(gridss)
    records, dup = tee(p.gather_records())
    processed_records = [p.process_record(record) for record in records]
    for record, processed in zip(dup, processed_records):
        record = record[0]
        assert processed["tumor_depth"] == record["tumor_VF"]
        assert processed["tumor_split_reads"] == record["tumor_SR"]
        assert processed["tumor_discordant_reads"] ==  record["tumor_RP"]
        assert processed["normal_depth"] == record["normal_VF"]
        assert processed["normal_split_reads"] == record["normal_SR"]
        assert processed["normal_discordant_reads"] == record["normal_RP"]

    # print(matched)
# def test_read_sample_names_freebayes(freebayes):
#     p = vcf_snv_parser.SNV_vcf(freebayes)
#     n, t = p._read_sample_names(freebayes)
#     assert t == None
#     assert n == "SA219N"

# def test_read_sample_names_rtg(rtg):
#     p = vcf_snv_parser.SNV_vcf(rtg)
#     n, t = p._read_sample_names(rtg)
#     assert t == None
#     assert n == "SA219N"

# def test_read_sample_names_samtools(samtools):
#     p = vcf_snv_parser.SNV_vcf(samtools)
#     n, t = p._read_sample_names(samtools)
#     assert t == None
#     assert n == "SA219N"

# def test_read_sample_names_mutect(mutect):
#     p = vcf_snv_parser.SNV_vcf(mutect)
#     n, t = p._read_sample_names(mutect)
#     assert t == "SA225"
#     assert n == "SA225N"
    
# def test_read_sample_names_strelka(strelka):
#     p = vcf_snv_parser.SNV_vcf(strelka)
#     n, t = p._read_sample_names(strelka)
#     assert t == "SA224"
#     assert n == "SA224N"

# def test_read_sample_names_museq_germline(museq_germline):
#     p = vcf_snv_parser.SNV_vcf(museq_germline)
#     n, t = p._read_sample_names(museq_germline)
#     assert t == None
#     assert n == "SA220N"


# def test_to_process_record_museq_somatic(museq_somatic):
#     p = vcf_snv_parser.Museq_vcf(museq_somatic)
#     records = p.gather_records()
#     record = next(records)
#     processed = p.process_record(record)
#     samples = record["samples"]
#     assert processed["tumor_depth"] == samples["DP_TUMOR"]
#     assert processed["tumor_alt_depth"] == samples["AC_TUMOR"]
#     assert processed["tumor_ref_depth"] == samples["RC_TUMOR"]
#     assert processed["normal_depth"] == samples["DP_NORMAL"]
#     assert processed["normal_alt_depth"] == samples["AC_NORMAL"]
#     assert processed["normal_ref_depth"] == samples["RC_NORMAL"]

# def test_to_process_record_freebayes(freebayes):
#     p = vcf_snv_parser.Freebayes_vcf(freebayes)
#     records = p.gather_records()
#     record = next(records)
#     processed = p.process_record(record)
#     samples = record["samples"]
#     assert isnan(processed["tumor_depth"])
#     assert isnan(processed["tumor_alt_depth"])
#     assert isnan(processed["tumor_ref_depth"])
#     assert processed["normal_depth"] == samples["DP_NORMAL"]
#     assert processed["normal_alt_depth"] == samples["AO_NORMAL"]
#     assert processed["normal_ref_depth"] == samples["RO_NORMAL"]

# def test_to_process_record_samtools(samtools):
#     p = vcf_snv_parser.Samtools_vcf(samtools)
#     records = p.gather_records()
#     record = next(records)
#     processed = p.process_record(record)
#     info = record["info"]
#     assert isnan(processed["tumor_depth"])
#     assert isnan(processed["tumor_alt_depth"])
#     assert isnan(processed["tumor_ref_depth"])
#     assert processed["normal_depth"] == info["DP"]
#     assert isnan(processed["normal_alt_depth"])
#     assert isnan(processed["normal_ref_depth"])

# def test_to_process_record_museq_germline(museq_germline):
#     p = vcf_snv_parser.Museq_vcf(museq_germline)
#     records = p.gather_records()
#     record = next(records)
#     processed = p.process_record(record)
#     samples = record["samples"]
#     assert isnan(processed["tumor_depth"])
#     assert isnan(processed["tumor_alt_depth"])
#     assert isnan(processed["tumor_ref_depth"])
#     assert processed["normal_depth"] == samples["DP_NORMAL"]
#     assert processed["normal_alt_depth"] == samples["AC_NORMAL"]
#     assert processed["normal_ref_depth"] == samples["RC_NORMAL"]

# def test_to_process_record_mutect(mutect):
#     p = vcf_snv_parser.Mutect_vcf(mutect)
#     records = p.gather_records()
#     record = next(records)
#     processed = p.process_record(record)
#     samples = record["samples"]
#     assert processed["tumor_depth"] == samples["DP_TUMOR"]
#     assert processed["tumor_alt_depth"] == samples["AD_TUMOR"].split(";")[1]
#     assert processed["tumor_ref_depth"] == samples["AD_TUMOR"].split(";")[0]
#     assert processed["normal_depth"] == samples["DP_NORMAL"]
#     assert processed["normal_alt_depth"] == samples["AD_NORMAL"].split(";")[1]
#     assert processed["normal_ref_depth"] == samples["AD_NORMAL"].split(";")[0]
# # def test_to_csv_museq_somatic(museq_somatic, test_output):
# #     p = vcf_snv_parser.Museq_vcf(museq_somatic)
# #     p.parse()
# #     data = pd.DataFrame(p.record_data)
# #     p.to_csv(test_output)
# #     yaml = "testoutput.csv.gz"

# def test_parse_main_cols_museq_somatic(museq_somatic):
#     p = vcf_snv_parser.SNV_vcf(museq_somatic)
#     record = next(p.reader)
#     main_cols = p.parse_main_cols(record)["main_cols"]
#     assert main_cols["chrom"] == record.CHROM
#     assert main_cols["pos"] == record.POS
#     assert main_cols["ref"] == record.REF
#     assert main_cols["alt"] == record.ALT[0]
#     assert main_cols["qual"] == record.QUAL
#     assert main_cols["filter"] == record.FILTER

# def test_parse_main_cols_freebayes(freebayes):
#     p = vcf_snv_parser.SNV_vcf(freebayes)
#     record = next(p.reader)
#     main_cols = p.parse_main_cols(record)["main_cols"]
#     assert main_cols["chrom"] == record.CHROM
#     assert main_cols["pos"] == record.POS
#     assert main_cols["ref"] == record.REF
#     assert main_cols["alt"] == record.ALT[0]
#     assert main_cols["qual"] == record.QUAL
#     assert main_cols["filter"] == record.FILTER

# def test_parse_main_cols_rtg(rtg):
#     p = vcf_snv_parser.SNV_vcf(rtg)
#     record = next(p.reader)
#     main_cols = p.parse_main_cols(record)["main_cols"]
#     assert main_cols["chrom"] == record.CHROM
#     assert main_cols["pos"] == record.POS
#     assert main_cols["ref"] == record.REF
#     assert main_cols["alt"] == record.ALT[0]
#     assert main_cols["qual"] == record.QUAL
#     assert main_cols["filter"] == record.FILTER

# def test_parse_main_cols_samtools(samtools):
#     p = vcf_snv_parser.SNV_vcf(samtools)
#     record = next(p.reader)
#     main_cols = p.parse_main_cols(record)["main_cols"]
#     assert main_cols["chrom"] == record.CHROM
#     assert main_cols["pos"] == record.POS
#     assert main_cols["ref"] == record.REF
#     assert main_cols["alt"] == record.ALT[0]
#     assert main_cols["qual"] == record.QUAL
#     assert main_cols["filter"] == record.FILTER

# def test_parse_main_cols_mutect(mutect):
#     p = vcf_snv_parser.SNV_vcf(mutect)
#     record = next(p.reader)
#     main_cols = p.parse_main_cols(record)["main_cols"]
#     assert main_cols["chrom"] == record.CHROM
#     assert main_cols["pos"] == record.POS
#     assert main_cols["ref"] == record.REF
#     assert main_cols["alt"] == record.ALT[0]
#     assert main_cols["qual"] == record.QUAL
#     assert main_cols["filter"] == record.FILTER
    
# def test_parse_main_cols_strelka(strelka):
#     p = vcf_snv_parser.SNV_vcf(strelka)
#     record = next(p.reader)
#     main_cols = p.parse_main_cols(record)["main_cols"]
#     assert main_cols["chrom"] == record.CHROM
#     assert main_cols["pos"] == record.POS
#     assert main_cols["ref"] == record.REF
#     assert main_cols["alt"] == record.ALT[0]
#     assert main_cols["qual"] == record.QUAL
#     assert main_cols["filter"] == record.FILTER

# def test_parse_main_cols_museq_germline(museq_germline):
#     p = vcf_snv_parser.SNV_vcf(museq_germline)
#     record = next(p.reader)
#     main_cols = p.parse_main_cols(record)["main_cols"]
#     assert main_cols["chrom"] == record.CHROM
#     assert main_cols["pos"] == record.POS
#     assert main_cols["ref"] == record.REF
#     assert main_cols["alt"] == record.ALT[0]
#     assert main_cols["qual"] == record.QUAL
#     assert main_cols["filter"] == record.FILTER

