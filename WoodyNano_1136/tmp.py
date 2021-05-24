
def position_extra_intron(max_distance, **kwargs):
    read_w = kwargs['woodynano']
    read_p = kwargs['pychopper']
    
    aligned_to, diff = aligned_toward(read_w, read_p)
    if diff > max_distance:
        raise Exception(f'Distance to long')

    n_extra_intron = (len(read_w.cigar_summary['extron']) -\
                    len(read_p.cigar_summary['extron'])) / 2

    if n_extra_intron > 0:
        software = 'woodynano'
        n_extra = int(n_extra_intron)
    elif n_extra_intron < 0:
        software = 'pychopper'
        n_extra = abs(int(n_extra_intron))

    # extract extra intron position from kwargs[software]
    read = kwargs[software]
    
    if aligned_to == 'R':
        first_extra = 0
    elif aligned_to == 'L':
        ttl_intron = (len(read.cigar_summary['extron'])-1)/2
        first_extra = ttl_intron - n_extra
    
    positions = []
    for n in range(first_extra, first_extra + n_extra_intron):
        positions.append(fragment_position(read, n + 1, 'intron'))
    
    return software, positions
    

def get_gene_locus(gff, read1, read2):
    aligned_to, diff = aligned_toward(read1, read2)
    if aligned_to == 'R':
        sticky_site = min(read1.reference_span[1], read2).reference_span[1]
    elif aligned_to == 'L':
        sticky_site = max(read1.reference_span[0], read2).reference_span[0]
    
    transcript_slice = gene_picker2(
        gene_file=gff,
        chromosome=read1.rname,
        sticky_site=sticky_site, 
        sites=aligned_to
        )
    
    if type(transcript_slice) == 'NoneType':
        # print('Novel gene expression locus')
        return (None,None)
        
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