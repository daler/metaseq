import metaseq


class Chipseq(object):
    def __init__(self, ip_bam, control_bam):
        self.ip = metaseq.genomic_signal(ip_bam, kind='bam')
        self.control = metaseq.genomic_signal(control_bam, kind='bam')
        self.ip_array = None
        self.control_array = None


    def diffed_array(self, features, force=False, func=None, **kwargs):
        """
        `func` is a function to apply to the diffed arrays, by default
        metaseq.plotutils.nice_log; another option might be `lambda x: x`, 
        or `lambda x: 1e6*x`
        """
        if not self.ip_array or force:
            self.ip_array = self.ip.array(features, **kwargs)
            self.ip_array /= self.ip.million_mapped_reads()
        if not self.control_array or force:
            self.control_array = self.control.array(features, **kwargs)
            self.controL_array /= self.control.million_mapped_reads()

        if func is None:
            metaseq.plotutils.nice_log

        return func(self.ip_array - self.control_array)
