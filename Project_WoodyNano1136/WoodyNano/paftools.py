import os
from WoodyNano import utils

'''
Define class alignment_results()
'''
'''
--for-only	Only map to the forward strand of the reference sequences. For paired-end reads in the forward-reverse orientation, the first read is mapped to forward strand of the reference and the second read to the reverse stand.
--rev-only	Only map to the reverse complement strand of the reference sequences.
-a	Generate CIGAR and output alignments in the SAM format. Minimap2 outputs in PAF by default.
-x splice	Long-read spliced alignment (-k15 -w5 --splice -g2000 -G200k -A1 -B2 -O2,32 -E1,0 
            -C9 -z200 -ub --junc-bonus=9 --splice-flank=yes). In the splice mode, 
            1) long deletions are taken as introns and represented as the ‘N’ CIGAR operator; 
            2) long insertions are disabled; 
            3) deletion and insertion gap costs are different during chaining; 
            4) the computation of the ‘ms’ tag ignores introns to demote hits to pseudogenes.
-C INT	Cost for a non-canonical GT-AG splicing (effective with --splice) [0]
-u CHAR	How to find canonical splicing sites GT-AG 
            - f: transcript strand; 
            - b: both strands; 
            - n: no attempt to match GT-AG [n]
'''



class PAF:
    def __init__(self):
        self.alignment = []
    
    class Alignment:
        def __init__(self):
            self.query_name = str
            self.query_len = int
            self.query_start = int
            self.query_end = int
            self.strand = chr
            self.target_name = str
            self.target_len = int
            self.target_start = int
            self.target_end = int
            self.residual_matched = int
            self.align_block_len = int
            self.mapq = int 
            self.opt = None

        def __repr__(self):
            return f"qname: {self.query_name}\
                    \tquery_length: {self.query_len}\
                    \tquery_start: {self.query_start}\
                    \tquery_end: {self.query_end}\
                    \nstrand: {self.strand}\
                    \ntarget_name: {self.target_name}\
                    \ttarget_len: {self.target_len}\
                    \ttarget_start: {self.target_start}\
                    \ttarget_end: {self.target_end}\
                    \nresidual_matched: {self.residual_matched}\
                    \talign_block_len: {self.align_block_len}\
                    \tmap_qual: {self.mapq}\
                    \nopt: {self.opt}\n"

        def effective_read_length(self):
            return int(self.query_end - self.query_start)

        def clipped_base(self):
            return int(self.query_len-self.effective_read_length())

        @staticmethod
        def asAlignment(lst):
            o = PAF.Alignment()
            o.query_name = lst[0]
            o.query_len = int(lst[1])
            o.query_start = int(lst[2])
            o.query_end = int(lst[3])
            o.strand = lst[4]
            o.target_name = lst[5]
            o.target_len = int(lst[6])
            o.target_start = int(lst[7])
            o.target_end = int(lst[8])
            o.residual_matched = int(lst[9])
            o.align_block_len = int(lst[10])
            o.mapq = int(lst[11])
            try:
                o.opt = {i.split(':')[0]: i.split(':')[-1] for i in lst[12:]}
            except IndexError:
                o.opt = None
            return o

        def same_loci(self, paf_align):
            _same = False
            if self.target_name != paf_align.target_name:
                pass
            else:
                if (self.target_start == paf_align.target_start) or (self.target_end, paf_align.target_end):
                    _same = True
                
                else: 
                    #check overlap
                    start_max = max(self.target_start,paf_align.target_start)
                    end_min = min(self.target_end, paf_align.target_end)
                    if start_max <= end_min:
                        _same = True
            return _same
        
    def get_primary_paf(self):
        new = PAF()
        for align in self.alignment:
            if align.opt:
                if align.opt['tp'] == 'P':
                    new.alignment.append(align)
        return new
    
    @staticmethod
    def Import(fname=None):
        if not fname:
            print('No file name provided')
            pass
        output = PAF()
        with open(fname) as o:
            alignment_section = o.readlines()
        for line in alignment_section:
            line = line.strip().split(sep='\t')
            output.alignment.append(PAF.Alignment.asAlignment(line))

        return output
    

    

    



    
