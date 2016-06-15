"""Convert vcf to a long format usable for 
'normal' database queries with columns:

  chrom
  pos
  individual
  alele position (diploids!)
  genotype
"""
import sys
import collections

def parse_sample(fields, sample):

    # lookup list for genotypes
    # reference and a list of alternatives
    genotypes = [fields[3]] + fields[4].split(',')

    # parse the top level given FORMAT
    fmt = fields[8]
    sfields = dict(zip(fmt.split(':'), sample.split(':')))

    # decode the genotype
    gt = sfields['GT']
    gtl = [genotypes[int(gt[0])], genotypes[int(gt[2])]]
    sfields['GT'] = gtl

    return sfields

def variants(it):
    """parse vcf
    """
    for line in it:
        # skip field definitions
        if line[:2] == "##":
            continue

        # read column names
        if line[0] == "#":
            colnames = line[1:].strip().split('\t')
            stdcols = colnames[:8]
            # column 9 contains sample format string
            samples = colnames[9:]
            Var = collections.namedtuple('var', stdcols + ['samples'])
            continue

        # parse standard rows
        fields = line.strip().split("\t")
        samples = {sample:parse_sample(fields, f) for sample, f in zip(samples, fields[9:])}
        yield Var(*(fields[:8] + [samples]))

def main():
    for var in variants(sys.stdin):
        for sample, genotypes in var.samples.iteritems():
            for i, gt in enumerate(genotypes['GT']):
                print "\t".join([var.CHROM, var.POS, sample, i, gt])

if __name__ == '__main__':
    main()
