


import click
from vcfutils import vcfutils
import warnings

COMMANDS = ("concat", "vcf2csv") 

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

    if cmd == "concat":
        vcfutils.concatenate_vcf(list(input_vcf), output)

    if cmd == "vcf2csv":
        if len(input_vcf) > 1:
            warnings.warn("vcf2csv can only take 1 vcf")

        input_vcf = input_vcf[0]
        
        if vcf_type == "lumpy":
            vcfutils.Lumpy_vcf(input_vcf).to_csv(output)

        elif vcf_type == "destruct":
            pass
            # vcfutils.Lumpy_vcf(input_vcf).to_csv(output)

        elif vcf_type == "svaba":
            vcfutils.Svaba_vcf(input_vcf).to_csv(output)

        elif vcf_type == "gridss":
            vcfutils.Gridss_vcf(input_vcf).to_csv(output)

        elif vcf_type == "strelka":
            pass
        elif vcf_type == "museq":
            pass

        else:
            warnings.warn("invalid vcf type, choose from (lumpy, destruct, svaba, gridss, strelka, museq)")
            return 
    pass
