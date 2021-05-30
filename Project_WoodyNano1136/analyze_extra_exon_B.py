#%%
import pandas as pd
import sys
from WoodyNano import samtools


def gene_picker(gene_file, chromosome, ref_span):
    name_correct = gene_file['seqname'] == chromosome
    left_in = gene_file['start'] <= ref_span[0]+200
    right_in = gene_file['end'] >= ref_span[1]-200
    tmp = gene_file[name_correct & left_in & right_in]
    # if not sum(name_correct & left_in & right_in) == 1:
    #     print(f'Not one {chromosome}')
    #     return
    # else:
    #     return [int(gene_file[name_correct & left_in & right_in][i]) for i in [3,4]]
    return tmp


def gene_picker2(gene_file, chromosome, sticky_site, sites):
    chr_slice = gene_file[gene_file['seqname'] == chromosome]
    if sites == 'L':
        tmp = abs(chr_slice['start'] - sticky_site)

    if sites == 'R':
        tmp = abs(chr_slice['end'] - sticky_site)

    gene_id = list(chr_slice[tmp == min(tmp)]['geneid'])

    try:
        out_slice = chr_slice[chr_slice['geneid'] == str(gene_id[0])]
        return out_slice

    except IndexError:
        return None


def gene_picker3(gene_file, chromosome, sticky_site, strand):
    name_correct = gene_file['seqname'] == chromosome
    left_in = gene_file['start'] <= sticky_site + 200
    right_in = gene_file['end'] >= sticky_site - 200
    strand_match = gene_file['strand'] == strand
    tmp = gene_file[name_correct & left_in & right_in]
    # if set(tmp['geneid']).__len__() == 1:
    #     return tmp
    # else:
    #     print(f'Weird loci {chromosome}:{sticky_site}')
    return tmp


def nth_extron(n, feature):
    """
    The nth intron is the (2n)th extron
    The nth exon is the (2n-1)th extron
    """
    if feature == 'exon':
        return 2*n - 1
    if feature == 'intron':
        return 2*n


def nth_extron_position(extron_list, ref_start, n):
    if n <= 0 or n > extron_list.__len__():
        raise Exception(
            'Extrons are counted from 1, 2,..., N; N is the numbers of extron')

    if n == 1:
        return ref_start, ref_start + extron_list[n-1] - 1

    elif n <= extron_list.__len__():
        last_end = nth_extron_position(extron_list, ref_start, n-1)[1] + 1
        return last_end, last_end + extron_list[n-1] - 1


def fragment_position(read, n, feature):
    extron_list = read.cigar_summary['extron']
    ref_start = read.reference_span[0]
    n = nth_extron(n, feature)
    return nth_extron_position(extron_list, ref_start, n)


def aligned_toward(sam_read1, sam_read2):
    start1, end1 = sam_read1.reference_span
    start2, end2 = sam_read2.reference_span
    diff_start, diff_end = abs(start1-start2), abs(end1-end2)

    if diff_start == diff_end:
        # raise Exception(f'Error: aligned_toward()\nRead: {sam_read1.qname} diff_start == diff_end\n')
        print(f'Read: {sam_read1.qname} diff_start == diff_end\n')
        # print(f'WoodyNano\n{sam_read1}\nPychopper\n{sam_read2}\n')
        return None, None
    else:
        return ('L', diff_start) if diff_start < diff_end else ('R', diff_end)


def position_extra_intron(max_distance, **kwargs):
    read_w = kwargs['woodynano']
    read_p = kwargs['pychopper']

    aligned_to, diff = aligned_toward(read_w, read_p)
    if aligned_to == None:
        return None, None

    else:
        if diff > max_distance:
            # raise Exception(f'Error: position_extra_intron()\n\
            #         Distance between {read_w.qname} and {read_p.qname} is to long: {diff}bp\n')
            print(f'Error: position_extra_intron()\n\
                Distance between {read_w.qname} and {read_p.qname} is to long: {diff}bp')
            # print(f'WoodyNano:\n{read_w}\nPychopper:\n{read_p}\n')

            return None, None

        else:
            n_extra_intron = (
                len(read_w.cigar_summary['extron']) - len(read_p.cigar_summary['extron'])) / 2

            if n_extra_intron > 0:
                software = 'woodynano'
                n_extra_intron = int(n_extra_intron)
            elif n_extra_intron < 0:
                software = 'pychopper'
                n_extra_intron = abs(int(n_extra_intron))

            # extract extra intron position from kwargs[software]
            read = kwargs[software]

            if aligned_to == 'R':
                first_extra = 0
            elif aligned_to == 'L':
                ttl_intron = (len(read.cigar_summary['extron'])-1)/2
                first_extra = ttl_intron - n_extra_intron

            positions = []
            for n in range(int(first_extra), int(first_extra + n_extra_intron)):
                positions.append(fragment_position(read, n + 1, 'intron'))

            return software, positions


def get_gene_locus(gff, read1, read2):
    aligned_to, diff = aligned_toward(read1, read2)
    if aligned_to == 'R':
        sticky_site = min(read1.reference_span[1], read2.reference_span[1])
    elif aligned_to == 'L':
        sticky_site = max(read1.reference_span[0], read2.reference_span[0])

    strand = '+' if not ('16' in read1.split_flag()) else '-'
    transcript_slice = gene_picker3(
        gene_file=gff,
        chromosome=read1.rname,
        sticky_site=sticky_site,
        strand=strand
    )

    if transcript_slice.shape[0] == 0:
        # print('Novel gene expression locus')
        return (None, None)

    else:
        annot = min(transcript_slice['start']), max(transcript_slice['end'])
        return annot


def classify_extra_intron(annot, positions):
    n_extra = len(positions)
    no_ref = 0
    n_in = 0
    n_out = 0

    if annot[0]:
        for pos in positions:
            if (pos[0] > annot[0]) and pos[1] < annot[1]:
                n_in += 1
            else:
                n_out += 1
    else:
        no_ref = n_extra

    return n_extra, no_ref, n_in, n_out


def main(path_align_summary, path_gtf, path_out, path_sam={'woodynano': '', 'pychopper': ''}):
    def last_name(alignment):
        return alignment.qname.split('|')[-1]

    # read align_summary file
    df = pd.read_csv(path_align_summary)
    c = (df['read_length_woodynano'] - df['read_length_pychopper'])
    cond1 = (c >= min(c)) & (c < 0)  # D
    cond2 = c == 0  # C
    cond3 = (c > 0) & (c < -min(c))  # B
    cond4 = c >= -min(c)  # A
    cond0 = df['numbers_of_exon_woodynano'] != df['numbers_of_exon_pychopper']
    df = df[cond0 & cond3]

    # read gene annotation file
    gtf = pd.read_csv(path_gtf, sep='\t', header=None)
    gtf.columns = [
        'seqname', 'source', 'feature',
        'start', 'end', 'score', 'strand',
        'frame', 'attribute']

    gtf['geneid'] = [_.split("\"")[-2] for _ in gtf['attribute']]
    transcript_gtf = gtf[gtf['feature'] == 'transcript']

    # read samfiles
    sam_woodynano = samtools.SAM.Import(
        fname=path_sam['woodynano'],
        get_primary=True
    )
    sam_pychopper = samtools.SAM.Import(
        fname=path_sam['pychopper'],
        get_primary=True
    )

    dict_woodynano = {}
    dict_pychopper = {}

    for i in sam_woodynano.alignment:
        if not last_name(i) in dict_woodynano.keys():
            dict_woodynano[last_name(i)] = i
    for i in sam_pychopper.alignment:
        if not last_name(i) in dict_pychopper.keys():
            dict_pychopper[last_name(i)] = i

    data = {
        'woodynano': {
            'n_in': 0,
            'n_out': 0,
            'no_ref': 0,
            'total': 0
        },
        'pychopper': {
            'n_in': 0,
            'n_out': 0,
            'no_ref': 0,
            'total': 0
        }
    }

    for readname in df['readname']:
        read_w = dict_woodynano[readname]
        read_p = dict_pychopper[readname]
        read_w.cal_all_stats()
        read_p.cal_all_stats()

        software, positions = position_extra_intron(
            max_distance=276,
            **{'woodynano': read_w, 'pychopper': read_p}
        )
        # print(type(software), type(positions))

        if software == None:
            print(f'Skip read: {read_w.qname}\n\n')
            continue

        annot = get_gene_locus(transcript_gtf, read_w, read_p)

        n_extra, no_ref, n_in, n_out = classify_extra_intron(annot, positions)

        data[software]['total'] += n_extra
        data[software]['no_ref'] += no_ref
        data[software]['n_in'] += n_in
        data[software]['n_out'] += n_out

    pd.DataFrame(data).to_csv(path_out, index=True)
    print('DONE!')
    return


# %%
path_woodynano = sys.argv[1]
path_pychopper = sys.argv[2]
path_align_summary = sys.argv[3]
path_gtf = sys.argv[4]
path_out = sys.argv[5]

# path_woodynano = '/Users/zhujiachen/Desktop/WoodyNano_Revision/remap/SAM/Ptr_bio1_full_length_woodynano.sam'
# path_pychopper = '/Users/zhujiachen/Desktop/WoodyNano_Revision/remap/SAM/Ptr_bio1_full_length_pychopper.sam'
# path_align_summary = '/Users/zhujiachen/Desktop/WoodyNano_Revision/alignment_summary.csv'
# path_gtf = "/Users/zhujiachen/Desktop/WoodyNano_Revision/GFF/Ptrichocarpa_533_v4.1.gene_exons.gtf"
# path_out = '/Users/zhujiachen/Desktop/WoodyNano_Revision/Ptr_bio1_types_of_extra_intron.csv'

main(
    path_align_summary=path_align_summary,
    path_gtf=path_gtf,
    path_out=path_out,
    path_sam={
        'woodynano': path_woodynano,
        'pychopper': path_pychopper
    }
)

# %%
