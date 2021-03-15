from vcfutils.parsers import vcf_snv_parser
from vcfutils.parsers import vcf_sv_parser







import os

test_data_dir = "/juno/work/shah/abramsd/CODE/testdata"
svaba = "/juno/work/shah/abramsd/CODE/testdata/out.svaba.germline.sv.vcf.gz"
lumpy = "/juno/work/shah/abramsd/CODE/testdata/SA1256PP_lumpy.vcf"
gridss = "/juno/work/shah/abramsd/CODE/testdata/gridss.vcf.gz"

freebayes = "/juno/work/shah/isabl_data_lake/analyses/74/87/7487/results/germline/SA219N/SA219N_freebayes_germline.vcf.gz"
rtg = "/juno/work/shah/isabl_data_lake/analyses/74/87/7487/results/germline/SA219N/SA219N_rtg_germline.vcf.gz"
samtools = "/juno/work/shah/isabl_data_lake/analyses/74/87/7487/results/germline/SA219N/SA219N_samtools_germline.vcf.gz"
museq = "/juno/work/shah/isabl_data_lake/analyses/74/86/7486/results/germline/SA220N/SA220N_museq_single_annotated.vcf.gz"

mutect = "/juno/work/shah/isabl_data_lake/analyses/77/97/7797/results/somatic/SA225/SA225_mutect.vcf.gz"
strelka = "/juno/work/shah/isabl_data_lake/analyses/77/96/7796/results/somatic/SA224/SA224_strelka_snv_annotated.vcf.gz"
museq_somatic = "/juno/work/shah/isabl_data_lake/analyses/77/96/7796/results/somatic/SA224/SA224_museq_paired_annotated.vcf.gz"
# testoutput_dir = "/juno/work/shah/abramsd/CODE/vcfutils/vcfutils/tests/test_outs"
# parser = vcf_sv_parser.Lumpy_vcf(lumpy)
# parser.parse()
# parser.to_csv(os.path.join(testoutput_dir, "lumpy_csv.csv"))

# parser = vcf_sv_parser.Gridss_vcf(gridss)
# parser.parse()
# parser.to_csv(os.path.join(testoutput_dir, "gridss_csv.csv"))

# parser = vcf_sv_parser.Svaba_vcf(svaba)
# parser.parse()
# parser.to_csv(os.path.join(testoutput_dir, "svaba_csv.csv"))

filter_out = []
filter_out.append(('LOW_MAPPABILITY', 'eq', True))
chromosomes = list((map(str,list(range(1, 23))))) + ["X"]
filter_out.append(('CHROM', 'notin', chromosomes))
filter_out.append(('PR', 'lt', 0))

# p = vcf_snv_parser.Freebayes_vcf(freebayes, filter_out)
# p.parse()
# p.to_csv("/juno/work/shah/abramsd/CODE/test_outs/freebayes.csv.gz")

# p = vcf_snv_parser.Rtg_vcf(rtg, filter_out)
# p.parse()
# p.to_csv("/juno/work/shah/abramsd/CODE/test_outs/rtg.csv.gz")

# p = vcf_snv_parser.Samtools_vcf(samtools, filter_out)
# p.parse()
# p.to_csv("/juno/work/shah/abramsd/CODE/test_outs/samtools.csv.gz")

# p = vcf_snv_parser.Museq_vcf(museq, filter_out)
# p.parse()
# p.to_csv("/juno/work/shah/abramsd/CODE/test_outs/museq_germline.csv.gz")

# p = vcf_snv_parser.Mutect_vcf(mutect, filter_out)
# p.parse()
# p.to_csv("/juno/work/shah/abramsd/CODE/test_outs/mutect.csv.gz")

# p = vcf_snv_parser.Strelka_vcf(strelka, filter_out)
# p.parse()
# p.to_csv("/juno/work/shah/abramsd/CODE/test_outs/strelka.csv.gz")

p = vcf_snv_parser.Museq_vcf(museq_somatic)
p.parse()
print(p.to_csv("/juno/work/shah/abramsd/CODE/test_outs/museq_csv.csv.gz"))