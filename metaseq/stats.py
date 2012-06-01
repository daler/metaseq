import fisher


def hypergeom(m, n, n1, n2):
    """
    From Fury et al., www.nslij-genetics.org/wli/pub/ieee-embs06.pdf

    :param m: overlapping genes
    :param n: total genes that could be sampled
    :param n1: number of genes in set 1
    "param n2: number of genes in set 2
    """
    return fisher.pvalue(*_fury_table(m, n, n1, n2)).right_tail


def _fury_table(m, n, n1, n2):
    """
    Convert to 2x2 table described in Fury et al, ready for use by `fisher`
    """
    return (m, n1-m, n2-m, n-n1-n2+m)


if __name__ == "__main__":
    from rpy2.robjects import r
    def _test_hypergeom(m, n, n1, n2):
        R_pval = r.phyper(min(n1, n2), n1, n - n1, n2)[0] \
                - r.phyper(m - 1, n1, n - n1, n2)[0]
        f_pval = fisher.pvalue(*_fury_table(m, n, n1, n2)).right_tail

        # at least to 10 sig figs
        R_str = ('%.10f' % R_pval)
        f_str = ('%.10f' % f_pval)

        print 'R:', R_str, 'Fisher:', f_str
        assert R_str == f_str

    to_try = [
            dict(m=40, n=1000, n1=100, n2=300),
            dict(m=10, n=1000, n1=100, n2=300),
            dict(m=40, n=10000, n1=1000, n2=3000),
            dict(m=40, n=10000, n1=100, n2=3000),
            ]
    for i in to_try:
        _test_hypergeom(**i)


