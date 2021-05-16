import os
from WoodyNano import utils
from parsy import regex
import re
import numpy as np

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


class Minimap2:
    '''
    Class to manage functions that call Minimap2 via os.system()
    '''
    @staticmethod
    def aligner(in_dir, ref, output_dir, suffix, parameters='-ax splice -uf -k 14 -t 48 --for-only'):
        fname = in_dir.split('/')[-1]
        out_sam = utils.change_suffix(fname, suffix)
        out_log = utils.change_suffix(fname, '.log')
        cmd = "nohup minimap2 " + parameters + " " + ref + " "+in_dir +\
            " 1> " + output_dir+out_sam+"  2> "+output_dir+out_log+' &'
        os.system(cmd)
        return


class SAM:
    def __init__(self):
        self.__header = []
        self.alignment = []

    @staticmethod
    def flag_tags():
        '''
        Infromation from: https://github.com/samtools/hts-specs
        '''
        return {
            '1': 'template having multiple segments in sequencing',
            '2': 'each segment properly aligned according to the aligner',
            '4': 'segment unmapped',
            '8': 'next segment in the template unmapped',
            '16': 'SEQ being reverse complemented',
            '32': 'SEQ of the next segment in the template being reverse complemented',
            '64': 'the first segment in the template',
            '128': 'the last segment in the template',
            '256': 'secondary alignment',
            '512': 'not passing filters, such as platform/vendor quality controls',
            '1024': 'PCR or optical duplicate',
            '2048': 'supplementary alignment'
        }

    @staticmethod
    def cigar_tags():
        return {
            'M': 'alignment match(can be a sequence match or mismatch)',
            'I': 'insertion to the reference',
            'D': 'deletion from the reference',
            'N': 'skipped region from the reference',
            'S': 'soft clipping(clipped sequences present in SEQ)',
            'H': 'hard clipping(clipped sequences NOT present in SEQ)',
            'P': 'padding(silent deletion from padded reference)',
            '=': 'sequence match',
            'X': 'sequence mismatch'
        }


    @staticmethod
    def Import(fname=None, get_primary=True):
        if not fname:
            print('No file name provided')
            pass
        output = SAM()
        with open(fname) as o:
            f = o.readline()
            while f[0] == '@':
                output.__header.append(f.strip())
                f = o.readline()
            alignment_section = [f] + o.readlines()

        for line in alignment_section:
            line = line.strip().split(sep='\t')
            new = SAM.Alignment.asAlignment(line)

            if get_primary:
                if 'tp:A:P' in new.opt:
                    output.alignment.append(new)
            else:
                output.alignment.append(new)

        return output


    def export(self, fname):
        with open(fname, 'w') as o:
            """
            write header
            """
            for lines in self.get_header(get_value=True):
                o.write(f'{lines}\n')

            for a in self.alignment:
                o.write(f'{a.to_line().strip()}\n')
        
        return 0

    
    def get_header(self, get_value=False):
        if get_value:
            return self.__header
        else:
            print(*self.__header, sep='\n')
            pass


    def flag_stats(self):
        flags_all = [_align.flag for _align in self.alignment]
        o = dict()
        for tag in set(flags_all):
            o[tag] = flags_all.count(tag)
        return o


    def extract_by_flag(self, flag_str):
        o = list()
        for align in self.alignment:
            if align.flag == flag_str:
                o.append(align)
        return o

    class Alignment:
        def __init__(self):
            self.qname = str()
            self.flag = int()
            self.rname = str()
            self.pos = int()
            self.mapq = int()
            self.cigar = str()
            self.rnext = str()
            self.pnext = int()
            self.tlen = int()
            self.seq = str()
            self.qual = str()
            self.opt = None
            self.cigar_summary = None
            self.mapped_read_length = None
            self.matching_nucleotides = None
            self.clipping = None
            self.read_length = None
            self.reference_span = None

        def __repr__(self):
            return f"qname: {self.qname}\
                    \tflag: {self.flag}\
                    \trname: {self.rname}\
                    \tpos: {self.pos}\
                    \tmapq: {self.mapq}\
                    \ncigar: {self.cigar}\
                    \nrnext: {self.rnext}\
                    \tpnext: {self.pnext}\
                    \ttlen: {self.tlen}\
                    \nseq: {self.seq}\
                    \nqual: {self.qual}\
                    \nopt:\n{self.opt}\
                    \nread_length:{self.read_length}\
                    \nmapped_read_length:{self.mapped_read_length}\
                    \nclipping:{self.clipping}\
                    \nmatching_nucleotides:{self.matching_nucleotides}\
                    \nreference_span:{self.reference_span}\
                    \ncigar_summary:{self.cigar_summary}\n"

        @staticmethod
        def asAlignment(lst):
            o = SAM.Alignment()
            o.qname = lst[0]
            o.flag = int(lst[1])
            o.rname = lst[2]
            o.pos = int(lst[3])
            o.mapq = int(lst[4])
            o.cigar = lst[5]
            o.rnext = lst[6]
            o.pnext = int(lst[7])
            o.tlen = int(lst[8])
            o.seq = lst[9]
            o.qual = lst[10]
            try:
                o.opt = lst[11:]
            except IndexError:
                return o
            return o
        

        def to_line(self):
            
            opt_line = ''
            
            for i in self.opt:
                opt_line += f'{i}\t'
            
            opt_line[:-3]

            return f"{self.qname}\t{self.flag}\t{self.rname}\t{self.pos}\t{self.mapq}\t{self.cigar}\t{self.rnext}\t{self.pnext}\t{self.tlen}\t{self.seq}\t{self.qual}\t{opt_line}\n"

        

        def summarize_cigar(self, exon_cigar, keep=False, out=False):
            
            def odd(length):
                return [(i+1)*2-1 for i in range(0, length)]

            def even(length):
                return [(i)*2 for i in range(0, length)]

            cigar = self.cigar
            char = np.array(re.split('[0-9]+', cigar)[1:])
            num = np.array([int(i) for i in re.split('[A-Z]|=', cigar)[:-1]])

            # record clipping
            clipping = [0, 0]
            for pos in [0,-1]:
                if char[pos] == 'S':
                    clipping[pos] = num[pos]
                    np.delete(char, pos)
                    np.delete(num, pos)

            # record exon-intron
            # record intron
            is_intron = char == 'N'
            num_intron = sum(is_intron)
            exon_intron_length = [0]*(2*num_intron+1)
            
            for idx, value in zip(odd(num_intron), num[is_intron]):
                exon_intron_length[idx] = value

            # record exon
            tmp = 0
            nth_exon = 0

            for idx in range(0, len(num)):
                if not is_intron[idx]:
                    if char[idx] in exon_cigar:
                        tmp += num[idx]
                else:
                    exon_intron_length[nth_exon*2] = tmp
                    tmp = 0
                    nth_exon += 1

            exon_intron_length[nth_exon*2] = tmp
            
            if keep:
                self.cigar_summary = {
                    'clipping': clipping,
                    'extron': exon_intron_length,
                    'exon_cigar':exon_cigar
                }
            
            if out:
                return {
                    'clipping': clipping,
                    'extron': exon_intron_length,
                    'exon_cigar':exon_cigar
                }

            pass

        def cal_ref_span(self, keep=False, out=False):
            spanned = None
            if not self.cigar_summary:
                spanned = sum(self.summarize_cigar(exon_cigar='=XD', keep=False, out=True)['extron']) - 1
            else:
                spanned = sum(self.cigar_summary['extron']) - 1
            
            if keep:
                    self.reference_span = self.pos, self.pos + spanned
            if out:
                return self.pos, self.pos + spanned

            return

        def cal_read_length(self, keep=False, out=False):
            if keep:
                    self.read_length = len(self.seq)
            if out:
                return len(self.seq)

            return

        def cal_mapped_read_length(self, keep=False, out=False):
            if self.clipping:
                clipped_bases = sum(self.clipping)
            else:
                clipped_bases = sum(self.summarize_cigar(exon_cigar='=XD', keep=False, out=True)['clipping'])

            if self.read_length:
                rlen = self.read_length
            else:
                rlen = len(self.seq)

            mapped_read_length = rlen - clipped_bases

            if keep:
                    self.mapped_read_length = mapped_read_length
            if out:
                return mapped_read_length
                
            return


        def cal_clipping(self, keep=False, out=False):
            if self.cigar_summary:
                clipping = self.cigar_summary['clipping']
            else:
                clipping = self.summarize_cigar(exon_cigar='=XD', keep=False, out=True)['clipping']

            if keep:
                    self.clipping = clipping[0], clipping[1]
            if out:
                return clipping[0], clipping[1]
            
            return

        def cal_matching_nucleotides(self, keep=False, out=False):
            cigar = self.cigar
            char = np.array(re.split('[0-9]+', cigar)[1:])
            num = np.array([int(i) for i in re.split('[A-Z]|=', cigar)[:-1]])

            sum(num[char == '='])
            matching = sum(num[char == '='])

            if keep:
                    self.matching_nucleotides = matching
            if out:
                return matching
            
            return 

        def cal_all_stats(self, keep=True):
            self.summarize_cigar(exon_cigar='=XD', keep=True, out=False)
            self.cal_read_length(keep=keep, out=False)
            self.cal_clipping(keep=keep, out=False)
            self.cal_mapped_read_length(keep=keep, out=False)
            self.cal_ref_span(keep=keep, out=False)
            self.cal_matching_nucleotides(keep=keep, out=False)
            return 