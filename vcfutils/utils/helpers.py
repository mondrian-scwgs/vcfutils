import warnings
import vcf
def _get_header(infile):
    '''
    Extract header from the VCF file
    :param infile: input VCF file vcf.Reader
    :return: header
    '''

    header = []
    for line in infile:
        if line.startswith('##'):
            header.append(line)
        elif line.startswith('#'):
            header.append(line)
            return header
        else:
            raise Exception('invalid header: missing #CHROM line')

    warnings.warn("One of the input files is empty")

    return []

def _get_reader(vcf_file):
    if vcf_file.endswith(".gz"):
        return vcf.Reader(filename=vcf_file)
    else:
        return vcf.Reader(open(vcf_file, "r"))