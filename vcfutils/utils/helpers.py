import warnings
import vcf
import itertools


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

def _group_iterator(iterable, chunk=100000):    
    while True:
        yield itertools.chain((next(iterable),), itertools.islice(iterable, chunk-1))

# def _group_iterator(iterable,n=100):
#     it = iter(iterable)
#     while True:
#         chunk = tuple(itertools.islice(it, n))
#         if not chunk:
#             return
#         yield chunk
