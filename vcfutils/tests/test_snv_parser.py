from vcfutils.parsers import vcf_snv_parser
import pytest
import os
import yaml
import pandas as pd
from numpy import isnan


@pytest.fixture
def freebayes():
    file = os.path.join(os.path.dirname(os.path.realpath(__file__)), "testdata", "freebayes.vcf.gz")
    assert os.path.exists(file)
    return file

@pytest.fixture
def rtg():
    file = os.path.join(os.path.dirname(os.path.realpath(__file__)), "testdata", "rtg.vcf.gz")
    assert os.path.exists(file)
    return file

@pytest.fixture
def samtools():
    file = os.path.join(os.path.dirname(os.path.realpath(__file__)), "testdata", "samtools.vcf.gz")
    assert os.path.exists(file)
    return file

@pytest.fixture
def museq_germline():
    file = os.path.join(os.path.dirname(os.path.realpath(__file__)), "testdata", "museq_germline.vcf.gz")
    assert os.path.exists(file)
    return file

@pytest.fixture
def test_output():
    o = "/juno/work/shah/abramsd/CODE/vcfutils/testoutput.csv.gz"
    if os.path.exists(o):
        os.path.remove(o)
    return o

@pytest.fixture
def museq_somatic():
    file = os.path.join(os.path.dirname(os.path.realpath(__file__)), "testdata", "museq_somatic.vcf.gz")
    assert os.path.exists(file)
    return file

@pytest.fixture
def strelka():
    file = os.path.join(os.path.dirname(os.path.realpath(__file__)), "testdata", "strelka.vcf.gz")
    assert os.path.exists(file)
    return file


@pytest.fixture
def mutect():
    file = os.path.join(os.path.dirname(os.path.realpath(__file__)), "testdata", "mutect.vcf.gz")
    assert os.path.exists(file)
    return file

def _check_yaml_formatting(y_file, parser):
    yaml_data = yaml.load(open(y_file))
    print(yaml_data)


def test_read_sample_names_museq_somatic(museq_somatic):
    p = vcf_snv_parser.SNV_vcf(museq_somatic)
    n, t = p._read_sample_names(museq_somatic)
    assert t == "SA224"
    assert n == "SA224N"

def test_read_sample_names_freebayes(freebayes):
    p = vcf_snv_parser.SNV_vcf(freebayes)
    n, t = p._read_sample_names(freebayes)
    assert t == None
    assert n == "SA219N"

def test_read_sample_names_rtg(rtg):
    p = vcf_snv_parser.SNV_vcf(rtg)
    n, t = p._read_sample_names(rtg)
    assert t == None
    assert n == "SA219N"

def test_read_sample_names_samtools(samtools):
    p = vcf_snv_parser.SNV_vcf(samtools)
    n, t = p._read_sample_names(samtools)
    assert t == None
    assert n == "SA219N"

def test_read_sample_names_mutect(mutect):
    p = vcf_snv_parser.SNV_vcf(mutect)
    n, t = p._read_sample_names(mutect)
    assert t == "SA225"
    assert n == "SA225N"
    
def test_read_sample_names_strelka(strelka):
    p = vcf_snv_parser.SNV_vcf(strelka)
    n, t = p._read_sample_names(strelka)
    assert t == "SA224"
    assert n == "SA224N"

def test_read_sample_names_museq_germline(museq_germline):
    p = vcf_snv_parser.SNV_vcf(museq_germline)
    n, t = p._read_sample_names(museq_germline)
    assert t == None
    assert n == "SA220N"


def test_to_process_record_museq_somatic(museq_somatic):
    p = vcf_snv_parser.Museq_vcf(museq_somatic)
    records = p.gather_records()
    record = next(records)
    processed = p.process_record(record)
    samples = record["samples"]
    assert processed["tumor_depth"] == samples["DP_TUMOR"]
    assert processed["tumor_alt_depth"] == samples["AC_TUMOR"]
    assert processed["tumor_ref_depth"] == samples["RC_TUMOR"]
    assert processed["normal_depth"] == samples["DP_NORMAL"]
    assert processed["normal_alt_depth"] == samples["AC_NORMAL"]
    assert processed["normal_ref_depth"] == samples["RC_NORMAL"]

def test_to_process_record_freebayes(freebayes):
    p = vcf_snv_parser.Freebayes_vcf(freebayes)
    records = p.gather_records()
    record = next(records)
    processed = p.process_record(record)
    samples = record["samples"]
    assert isnan(processed["tumor_depth"])
    assert isnan(processed["tumor_alt_depth"])
    assert isnan(processed["tumor_ref_depth"])
    assert processed["normal_depth"] == samples["DP_NORMAL"]
    assert processed["normal_alt_depth"] == samples["AO_NORMAL"]
    assert processed["normal_ref_depth"] == samples["RO_NORMAL"]

def test_to_process_record_samtools(samtools):
    p = vcf_snv_parser.Samtools_vcf(samtools)
    records = p.gather_records()
    record = next(records)
    processed = p.process_record(record)
    info = record["info"]
    assert isnan(processed["tumor_depth"])
    assert isnan(processed["tumor_alt_depth"])
    assert isnan(processed["tumor_ref_depth"])
    assert processed["normal_depth"] == info["DP"]
    assert isnan(processed["normal_alt_depth"])
    assert isnan(processed["normal_ref_depth"])

def test_to_process_record_museq_germline(museq_germline):
    p = vcf_snv_parser.Museq_vcf(museq_germline)
    records = p.gather_records()
    record = next(records)
    processed = p.process_record(record)
    samples = record["samples"]
    assert isnan(processed["tumor_depth"])
    assert isnan(processed["tumor_alt_depth"])
    assert isnan(processed["tumor_ref_depth"])
    assert processed["normal_depth"] == samples["DP_NORMAL"]
    assert processed["normal_alt_depth"] == samples["AC_NORMAL"]
    assert processed["normal_ref_depth"] == samples["RC_NORMAL"]

def test_to_process_record_mutect(mutect):
    p = vcf_snv_parser.Mutect_vcf(mutect)
    records = p.gather_records()
    record = next(records)
    processed = p.process_record(record)
    samples = record["samples"]
    assert processed["tumor_depth"] == samples["DP_TUMOR"]
    assert processed["tumor_alt_depth"] == samples["AD_TUMOR"].split(";")[1]
    assert processed["tumor_ref_depth"] == samples["AD_TUMOR"].split(";")[0]
    assert processed["normal_depth"] == samples["DP_NORMAL"]
    assert processed["normal_alt_depth"] == samples["AD_NORMAL"].split(";")[1]
    assert processed["normal_ref_depth"] == samples["AD_NORMAL"].split(";")[0]
# def test_to_csv_museq_somatic(museq_somatic, test_output):
#     p = vcf_snv_parser.Museq_vcf(museq_somatic)
#     p.parse()
#     data = pd.DataFrame(p.record_data)
#     p.to_csv(test_output)
#     yaml = "testoutput.csv.gz"

def test_parse_main_cols_museq_somatic(museq_somatic):
    p = vcf_snv_parser.SNV_vcf(museq_somatic)
    record = next(p.reader)
    main_cols = p.parse_main_cols(record)["main_cols"]
    assert main_cols["chrom"] == record.CHROM
    assert main_cols["pos"] == record.POS
    assert main_cols["ref"] == record.REF
    assert main_cols["alt"] == record.ALT[0]
    assert main_cols["qual"] == record.QUAL
    assert main_cols["filter"] == record.FILTER

def test_parse_main_cols_freebayes(freebayes):
    p = vcf_snv_parser.SNV_vcf(freebayes)
    record = next(p.reader)
    main_cols = p.parse_main_cols(record)["main_cols"]
    assert main_cols["chrom"] == record.CHROM
    assert main_cols["pos"] == record.POS
    assert main_cols["ref"] == record.REF
    assert main_cols["alt"] == record.ALT[0]
    assert main_cols["qual"] == record.QUAL
    assert main_cols["filter"] == record.FILTER

def test_parse_main_cols_rtg(rtg):
    p = vcf_snv_parser.SNV_vcf(rtg)
    record = next(p.reader)
    main_cols = p.parse_main_cols(record)["main_cols"]
    assert main_cols["chrom"] == record.CHROM
    assert main_cols["pos"] == record.POS
    assert main_cols["ref"] == record.REF
    assert main_cols["alt"] == record.ALT[0]
    assert main_cols["qual"] == record.QUAL
    assert main_cols["filter"] == record.FILTER

def test_parse_main_cols_samtools(samtools):
    p = vcf_snv_parser.SNV_vcf(samtools)
    record = next(p.reader)
    main_cols = p.parse_main_cols(record)["main_cols"]
    assert main_cols["chrom"] == record.CHROM
    assert main_cols["pos"] == record.POS
    assert main_cols["ref"] == record.REF
    assert main_cols["alt"] == record.ALT[0]
    assert main_cols["qual"] == record.QUAL
    assert main_cols["filter"] == record.FILTER

def test_parse_main_cols_mutect(mutect):
    p = vcf_snv_parser.SNV_vcf(mutect)
    record = next(p.reader)
    main_cols = p.parse_main_cols(record)["main_cols"]
    assert main_cols["chrom"] == record.CHROM
    assert main_cols["pos"] == record.POS
    assert main_cols["ref"] == record.REF
    assert main_cols["alt"] == record.ALT[0]
    assert main_cols["qual"] == record.QUAL
    assert main_cols["filter"] == record.FILTER
    
def test_parse_main_cols_strelka(strelka):
    p = vcf_snv_parser.SNV_vcf(strelka)
    record = next(p.reader)
    main_cols = p.parse_main_cols(record)["main_cols"]
    assert main_cols["chrom"] == record.CHROM
    assert main_cols["pos"] == record.POS
    assert main_cols["ref"] == record.REF
    assert main_cols["alt"] == record.ALT[0]
    assert main_cols["qual"] == record.QUAL
    assert main_cols["filter"] == record.FILTER

def test_parse_main_cols_museq_germline(museq_germline):
    p = vcf_snv_parser.SNV_vcf(museq_germline)
    record = next(p.reader)
    main_cols = p.parse_main_cols(record)["main_cols"]
    assert main_cols["chrom"] == record.CHROM
    assert main_cols["pos"] == record.POS
    assert main_cols["ref"] == record.REF
    assert main_cols["alt"] == record.ALT[0]
    assert main_cols["qual"] == record.QUAL
    assert main_cols["filter"] == record.FILTER




# import os
# ls
# test_data_dir = "/juno/work/shah/abramsd/CODE/testdata"
# svaba = "/juno/work/shah/abramsd/CODE/testdata/out.svaba.germline.sv.vcf.gz"
# lumpy = "/juno/work/shah/abramsd/CODE/testdata/SA1256PP_lumpy.vcf"
# gridss = "/juno/work/shah/abramsd/CODE/testdata/gridss.vcf.gz"

# freebayes = "/juno/work/shah/isabl_data_lake/analyses/74/87/7487/results/germline/SA219N/SA219N_freebayes_germline.vcf.gz"
# rtg = "/juno/work/shah/isabl_data_lake/analyses/74/87/7487/results/germline/SA219N/SA219N_rtg_germline.vcf.gz"
# samtools = "/juno/work/shah/isabl_data_lake/analyses/74/87/7487/results/germline/SA219N/SA219N_samtools_germline.vcf.gz"
# museq = "/juno/work/shah/isabl_data_lake/analyses/74/86/7486/results/germline/SA220N/SA220N_museq_single_annotated.vcf.gz"

# mutect = "/juno/work/shah/isabl_data_lake/analyses/77/97/7797/results/somatic/SA225/SA225_mutect.vcf.gz"
# strelka = "/juno/work/shah/isabl_data_lake/analyses/77/96/7796/results/somatic/SA224/SA224_strelka_snv_annotated.vcf.gz"
# museq_somatic = "/juno/work/shah/isabl_data_lake/analyses/77/96/7796/results/somatic/SA224/SA224_museq_paired_annotated.vcf.gz"
# # testoutput_dir = "/juno/work/shah/abramsd/CODE/vcfutils/vcfutils/tests/test_outs"
# # parser = vcf_sv_parser.Lumpy_vcf(lumpy)
# # parser.parse()
# # parser.to_
# 
# (os.path.join(testoutput_dir, "lumpy_csv.csv"))

# # parser = vcf_sv_parser.Gridss_vcf(gridss)
# # parser.parse()
# # parser.to_csv(os.path.join(testoutput_dir, "gridss_csv.csv"))

# # parser = vcf_sv_parser.Svaba_vcf(svaba)
# # parser.parse()
# # parser.to_csv(os.path.join(testoutput_dir, "svaba_csv.csv"))

# filter_out = []
# filter_out.append(('LOW_MAPPABILITY', 'eq', True))
# chromosomes = list((map(str,list(range(1, 23))))) + ["X"]
# filter_out.append(('CHROM', 'notin', chromosomes))
# filter_out.append(('PR', 'lt', 0))

# # p = vcf_snv_parser.Freebayes_vcf(freebayes, filter_out)
# # p.parse()
# # p.to_csv("/juno/work/shah/abramsd/CODE/test_outs/freebayes.csv.gz")

# # p = vcf_snv_parser.Rtg_vcf(rtg, filter_out)
# # p.parse()
# # p.to_csv("/juno/work/shah/abramsd/CODE/test_outs/rtg.csv.gz")

# # p = vcf_snv_parser.Samtools_vcf(samtools, filter_out)
# # p.parse()
# # p.to_csv("/juno/work/shah/abramsd/CODE/test_outs/samtools.csv.gz")

# # p = vcf_snv_parser.Museq_vcf(museq, filter_out)
# # p.parse()
# # p.to_csv("/juno/work/shah/abramsd/CODE/test_outs/museq_germline.csv.gz")

# # p = vcf_snv_parser.Mutect_vcf(mutect, filter_out)
# # p.parse()
# # p.to_csv("/juno/work/shah/abramsd/CODE/test_outs/mutect.csv.gz")

# # p = vcf_snv_parser.Strelka_vcf(strelka, filter_out)
# # p.parse()
# # p.to_csv("/juno/work/shah/abramsd/CODE/test_outs/strelka.csv.gz")

# p = vcf_snv_parser.Museq_vcf(museq_somatic)
# p.parse()
# print(p.to_csv("/juno/work/shah/abramsd/CODE/test_outs/museq_csv.csv.gz"))