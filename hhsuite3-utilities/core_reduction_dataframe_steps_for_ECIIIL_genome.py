# MEANT TO BE RUN WITH `Checking for reduction in size...` or `Showing reductuion in size...` notebooks
# with `%run -i core_reduction_dataframe_steps_for_ECIIIL_genome.py`
#------------------------------------------------------------------------------------------------------
# merge and drop the `hit_num` column from each df. Also remove the columns of 
# superfluous values and alignment info since we want to focus on size change.
# (The alignment info is the last few columns and so can use index to specify those.)
# We'll keep some of the assessements to highlight some that don't have convincing 
# orthologs, yet -- in other words, may need futher sorting as 'scRPR2' and 'scRPR2' did.
import pandas as pd
df = pd.concat(size_change_summary_df_lines_list, ignore_index=True)
columns_to_drop = (['hit_num', 'E-value','Identities','Similarity', 'Aligned_cols',
                    'Sum_probs','Template_Neff'] + list(df.columns[-6:]))
df = df.drop(columns=columns_to_drop)
# Add a new column which is the percent reduction in size
def calculate_percent_reduction(items):
    return ("{:.1%}".format(1 - (items[1]/float(items[0]))))
df['%reduction'] = df[['query_length','hit_length']].apply(calculate_percent_reduction, axis=1)
# Add a new column which is the estimated lost size
def estimate_mass_reduction(items):
    return ("{:.1f}".format(items.size_diff * 110/1000))
df['lost kDa est.'] = df.apply(estimate_mass_reduction, axis=1)
df = df.sort_values('qid', ascending=True)
df = df.reset_index(drop=True)
df = df.rename(columns={'qid':'id'})
df = df.rename(columns={'hid':'hit_id'})
# Add a new column which is the hit locus derived from the description, a.k.a. htitle
def extract_locus(items):
    chr = items.htitle.split("ECIIIL:CH",1)[1].split("_ECIII-L:",1)[0]
    start_n_end_n_frame_data = items.htitle.split("_ECIII-L:",1)[1]
    start,end,frameNmore,*_= start_n_end_n_frame_data.split(":")
    frame = frameNmore.split()[0]
    return ("chr {}:{}-{} ({})".format(chr,start,end,frame))
df.insert(6, 'hit_locus (strand)', df.apply(extract_locus, axis=1))
df = df.drop('htitle', axis=1)
df
