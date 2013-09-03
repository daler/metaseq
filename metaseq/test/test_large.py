"""
module for testing the larger files (x.bam, x.bed.gz, etc)
"""
import multiprocessing
import metaseq
import pybedtools

CPUS = multiprocessing.cpu_count()

gs = {}
for kind in ['bam', 'bigwig', 'bed', 'bigbed']:
    if kind == 'bed':
        ext = 'bed.gz'
    else:
        ext = kind
    gs[kind] = metaseq.genomic_signal(
        metaseq.example_filename('x.%s' % ext), kind)

# generate the test features
features = pybedtools.BedTool()\
        .window_maker(
            b=pybedtools.BedTool('chr2L 0 500000',
                                 from_string=True).fn,
            w=1000)\
        .shuffle(seed=1,
                 genome={'chr2L': (0, 5000000)})

args = (features,)
kwargs = dict(processes=CPUS, bins=100)
bam_array = gs['bam'].array(*args, **kwargs)
bed_array = gs['bed'].array(*args, **kwargs)
bw_array = gs['bigwig'].array(*args, method='get_as_array', **kwargs)
bb_array = gs['bigbed'].array(*args, **kwargs)


assert (bam_array == bed_array).all()
assert (bb_array == bed_array).all()
assert (bw_array == bed_array).all()

assert (bb_array == bam_array).all()
assert (bw_array == bam_array).all()

assert (bw_array == bb_array).all()
