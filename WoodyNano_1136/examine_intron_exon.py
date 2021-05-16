# %%
from WoodyNano import samtools, paftools

# %%
def last_name(alignment):
    return alignment.qname.split('|')[-1]

tmp = samtools.SAM.Import(
    fname='/Users/zhujiachen/Desktop/remap/Ptr_bio1_full_length_woodynano.sam',
    get_primary=True
)
for i in tmp.alignment:
    i.cal_all_stats(keep=True)

samfile_w = {last_name(i):i for i in tmp.alignment}

#%%
paffile_W = paftools.PAF.Import(
    fname='/Users/zhujiachen/Desktop/remap/Ptr_bio1_full_length_woodynano.paf'
)

# tmp = samtools.SAM.Import(
#     fname='/Users/zhujiachen/Desktop/SAM/Ptr_bio1_full_length_pychopper.sam',
#     get_primary=True
# )

# for i in tmp.alignment:
#     i.cal_all_stats(keep=True)

# samfile_p = {last_name(i):i for i in tmp.alignment}

# both_mapped = set(samfile_p.keys()) & set(samfile_w.keys())
# print(len(both_mapped))

# for name in list(both_mapped):
#     # check if overlap
#     align_w = samfile_w[name]
#     align_p = samfile_p[name]

#     if not align_w.rname == align_p.rname:
#         both_mapped.remove(name)
#     else:
#         span_w = align_w.reference_span
#         span_p = align_p.reference_span
#         # check overlap
#         if not min(span_w[-1],span_p[-1]) > max(span_w[0],span_p[0]):
#             both_mapped.remove(name)



# # %%
# samfile_w = {i: samfile_w[i] for i in both_mapped}
# samfile_p = {i: samfile_p[i] for i in both_mapped}

# # %%
# %%












#%%


