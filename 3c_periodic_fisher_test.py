
import pandas as pd
import scipy.stats as st
df_cog = pd.read_csv('cog_fishertest.csv',header=0)
df_cog.head(2)
globalsum = sum(df_cog['global'])
mutantsum = sum(df_cog['mutant'])

pvals = {}
for i,r in df_cog.iterrows():
    category = [r['mutant'],r['global'] ]
    noncat = [mutantsum - r['mutant'],globalsum - r['global'] ]
    oddsratio, pvalue = st.fisher_exact([category, noncat])

    pvals[r['cog']] = pvalue

for x in pvals:
    if pvals[x]<0.05:
        print x, ':', pvals[x]