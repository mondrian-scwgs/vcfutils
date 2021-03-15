from vcfutils.parsers import vcf_snv_parser, vcf_sv_parser
import click
import warnings

COMMANDS = ("vcf2csv") 

@click.command()
@click.argument('cmd')
@click.argument('output')

@click.option('--input_vcf', help='input vcfs for concat',  multiple=True, default=None)
@click.option('--vcf_type', help='type of vcf', default=None)

def main(cmd, output, input_vcf=None, vcf_type=None):
    if input_vcf==():
        warnings.warn("you must provide at least 1 vcf")
        return 

    if cmd not in COMMANDS:
        warnings.warn("invalid command, choose from (concat, vcf2csv)")
        return


    if cmd == "vcf2csv":
        if len(input_vcf) > 1:
            warnings.warn("vcf2csv can only take 1 vcf")

        input_vcf = input_vcf[0]
        
        if vcf_type == "lumpy":
            parser = vcf_sv_parser.Lumpy_vcf(input_vcf)

        elif vcf_type == "svaba":
            parser = vcf_sv_parser.Svaba_vcf(input_vcf)

        elif vcf_type == "gridss":
            parser = vcf_sv_parser.Gridss_vcf(input_vcf)

        elif vcf_type == "strelka":
            parser = vcf_snv_parser.Strelka_vcf(input_vcf)

        elif vcf_type == "freebayes":
            parser = vcf_snv_parser.Freebayes_vcf(input_vcf)

        elif vcf_type == "museq":
            parser = vcf_snv_parser.Museq_vcf(input_vcf)
        
        elif vcf_type == "rtg":
            parser = vcf_snv_parser.Rtg_vcf(input_vcf)

        elif vcf_type == "samtools":
            parser = vcf_snv_parser.Samtools_vcf(input_vcf)

        elif vcf_type == "mutect":
            parser = vcf_snv_parser.Mutect_vcf(input_vcf)

        else:
            warnings.warn("invalid vcf type, choose from (lumpy, destruct, svaba, gridss, strelka, museq)")
            return 
        parser.parse()
        parser.to_csv(output)
