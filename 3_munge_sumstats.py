#!/usr/bin/env python
# coding: utf-8

# In[1]:


import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import os
import glob
import time
import subprocess as sp
import multiprocessing as mp
import psutil

from basic_tools import *


debug=False


# # functions used for loading gwas result files

# In[2]:


def save_sumstats(gwas_result,pheno_code):
    dest=premunge_path.format(pheno_code)
    if os.path.exists(dest) or os.path.exists(sumstats_path.format(pheno_code)+'.sumstats.gz'):
        print("Congratulations!. premunged sumstat of",pheno_code,"exists. passed.")
        return
    elif os.path.exists(sumstats_path.format(pheno_code)+'.sumstats.gz'):
        print("Congratulations!. munged sumstat of",pheno_code,"exists. passed.")
        return       
    print("Saving sumstats:",dest)
    tosave=gwas_result
    tosave.to_csv(dest,sep=' ',columns=['SNP','A1','A2','N','BETA','P','MAF','INFO'],index=False)
    print("Finished saving sumstats:",dest)

def munge_sumstats(pheno_code):
    if os.path.exists(sumstats_path.format(pheno_code)+'.sumstats.gz'):
        print("Congratulations!. munged sumstats of",pheno_code,"exists. passed.")
        return
    script=['munge_sumstats.py','--sumstats',premunge_path.format(pheno_code),'--merge-alleles',snplist_path,'--out',sumstats_path.format(pheno_code)]
    print('Started:',' '.join(script))
    sp.call(script) 
    print('Finished:',' '.join(script))
    #sp.call(['munge_sumstats.py','--sumstats',premunge_path.format(pheno_code),'--merge-alleles',snplist_path,'--out',sumstats_path.format(pheno_code)])


# In[3]:


if debug:
    gwas_result=read_gwas2('23115')
    gwas_result.to_pickle('23115.pickle')


# In[4]:


def for_mp(pheno_code):
    gwas_result=read_gwas2(pheno_code)
    save_sumstats(gwas_result,pheno_code)
    munge_sumstats(pheno_code)


# In[6]:


#pheno_code_list=ukbb_table_filtered['Phenotype Code'].str.replace('_irnt','')


# In[4]:


start=int(sys.argv[1])
end=int(sys.argv[2])


# In[23]:


#phenotypes_filtered


# ```
# #!/bin/sh
# chrN=$1
# 
# echo $chrN "parallel started"
# 
# qctool -g "/data/UKbiobank/ukb_genetic_info/imp/ukb_imp_chr${chrN}_v3.bgen" -s "/data/UKbiobank/ukb_genetic_info/imp/ukb46263_imp_chr${chrN}_v3_s487320.sample" -og "/data01/ch6845/ldsc_data/ukb_imp_plink/ukb_imp_chr${chrN}_v3.bed" -incl-samples "/data01/ch6845/ldsc_data/ukb_sub.sample"
# ```

# ```
# jupyter nbconvert 3_munge_sumstats.ipynb --to script && mkdir ~/premunge &&python 3_munge_sumstats.py 0 50
# ```

# In[6]:


#start,end=383,384


# In[9]:


#phenotypes_filtered['phenotype']


# In[21]:


phenotypes_filtered.shape


# In[ ]:


for idx,row in phenotypes_filtered.iloc[start:end].iterrows():
    pheno_code=idx
    if os.path.exists(sumstats_path.format(pheno_code)+'.sumstats.gz'):
        print("Congratulations!. passed",pheno_code,"exists. passed.")
        continue
    print(pheno_code)
    gwas_result=read_gwas2(pheno_code)
    save_sumstats(gwas_result,pheno_code)
    munge_sumstats(pheno_code)
    
    run_command('rm {}'.format(premunge_path.format(pheno_code)))
    
#    print(pheno_code)


# In[23]:


idx_list=[]
for idx,row in phenotypes_filtered.iterrows():
    pheno_code=idx
    if os.path.exists(sumstats_path.format(pheno_code)+'.sumstats.gz'):
        #print("Congratulations!. passed",row['phenotype'],"exists. passed.")
        continue
    idx_list.append(idx)
    print(idx)
    continue
    gwas_result=read_gwas2(pheno_code)
    save_sumstats(gwas_result,pheno_code)
    munge_sumstats(pheno_code)
    
    run_command('rm {}'.format(premunge_path.format(pheno_code)))
    
#    print(pheno_code)


# In[18]:


for i in range(0, len(idx_list), int(len(idx_list)/25)):
    print(idx_list[i:i + int(len(idx_list)/25)])


# In[5]:


"""
cnt=0
for idx,pheno_code in ukbb_table_filtered['Phenotype Code'].str.replace('_irnt','').iteritems():
    if os.path.exists(sumstats_path.format(pheno_code)+'.sumstats.gz'):
        continue
    print(idx,pheno_code)
    cnt+=1
"""


# In[ ]:


"""
a=pd.read_csv('/home/ch6845/premunge/20426.premunge',sep=' ')
a.shape,a[np.logical_not(a['BETA']=='Infinity')].shape,a[np.logical_not(a['BETA']=='-Infinity')].shape
a[np.logical_not(a['BETA']=='Infinity')].to_csv('/home/ch6845/premunge/20426.premunge',sep=' ',index=None)
"""

