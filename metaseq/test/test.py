import metaseq
gs = {}
for kind in ['bed', 'bam', 'bigbed', 'bigwig']:
    gs[kind] = metaseq.genomic_signal(metaseq.example_filename('gdc.%s' % kind), kind)

def check_local_count(kind, coord, expected, stranded):
    result = gs[kind].local_count(coord, stranded=stranded)
    print gs[kind].fn
    assert result == expected, (kind, coord, result)

def test_local_count():
    for kind in ['bam', 'bigbed', 'bed']:
        for coord, expected, stranded in (
            ('chr2L:1-80', 3, False),       #  easy case
            ('chr2L:1000-3000', 0, False),  #  above upper boundary
            ('chr2L:1-9', 0, False),        #  below lower boundary
            ('chr2L:71-73[-]', 2, False),   #  unstranded = 2
            ('chr2L:71-73[-]', 1, True),    #  stranded = 1
            ('chr2L:70-71', 2, False),      #  pathological corner case
            ('chr2L:75-76', 0, False),      #  pathological corner case
        ):
            yield check_local_count, kind, coord, expected, stranded
