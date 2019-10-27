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

from path_configure import *

debug=False


# In[14]:


def save_annot_chr_cm(pheno_code,chrN,scale=1,mode='cm'):
    dest=annot_path.format(pheno_code,chrN)
    if os.path.exists(dest):
        print("Congratulations!. annot of",pheno_code,chrN,"exists. passed.")
        return
    tosave=pd.read_csv(bim_path.format(chrN),sep='\t',names=['CHR','SNP','CM','BP','A1','A2'])

    if mode=='cm':
        cm_series=tosave['CM']
    elif mode=='bp':
        cm_series=tosave['BP']/1000000
    cm_max=int(np.ceil(np.max(cm_series)/scale)*scale)
    cm_min=int(np.floor(np.min(cm_series)/scale)*scale)
    print(pheno_code,"chr:",chrN,"scale:",scale,"{} min,max".format(mode),cm_min,cm_max)
    category=np.zeros(len(cm_series),dtype=np.int32)
    category_name=[]
    category_count=0
    for i in range(cm_min,cm_max,scale):
        cm_index=np.logical_and(cm_series>=(i),cm_series<(i+scale))
        print(pheno_code,"chr:",chrN,"SNPs in {} range".format(mode),i,i+scale,"are",sum(cm_index),end='')  
        if sum(cm_index)>100:
            category_count=category_count+1
            category[cm_index]=category_count
            category_name.append(str(chrN)+'.'+str(i))
            print('-> included')
        else:
            print('-> ignored')
            #print(pheno_code,"chr:",chrN,"SNPs in {} range".format(mode),i,i+scale,"are ignored")  

    category_onehot=pd.DataFrame(np.vstack([np.zeros(category_count),np.eye(category_count)])[category],columns=category_name,dtype='int32')

    tosave.reset_index(drop=False,inplace=True)
    tosave=pd.concat([tosave,category_onehot],axis=1, join='inner')
    
    tosave.to_csv(dest,sep='\t',columns=['CHR','BP','SNP','CM']+category_name,index=False)
    print("Saved",dest)  


# In[3]:


#Process level parallelism for shell commands
def estimate_ldscore_chrN(annot_name,chrN,mode='annot'):
    """Defines the work unit on an input file"""
    #ldsc.py --l2 --bfile /scratch/ch6845/ukb_ld/ukb_imp_plink_gwas/ukb_imp_chr1_v3 --ld-wind-cm 1 --annot /scratch/ch6845/ukb_ld/annot/23115.chr1.annot --out --print-snps /scratch/ch6845/ukb_ld/hapmap3_snps/hm.1.snp
    #ldsc.py --l2 --bfile /data01/ch6845/ukb_imp_plink_gwas/ukb_imp_chr9_v3 --ld-wind-cm 1 --annot /scratch/ch6845/ukb_ld/annot/23115.9.annot --out test2 --print-snps /scratch/ch6845/ukb_ld/hapmap3_snps/hm.9.snp
    if mode=='annot':
        if os.path.exists(ld_path.format(annot_name,chrN)+'.l2.ldscore.gz'):
            print("Congratulations!. ld score file of",annot_name,chrN,"exists. passed.")
            return 0
        scripts=['ldsc.py','--l2','--bfile',bfile_path.format(chrN), '--ld-wind-cm', '1','--annot',annot_path.format(annot_name,chrN),'--out',ld_path.format(annot_name,chrN),'--print-snps',print_snps_path.format(chrN)]
    elif mode=='univariate':
        annot_name='uni'
        if os.path.exists(ld_path.format(annot_name,chrN)+'.l2.ldscore.gz'):
            print("Congratulations!. ld score file of",chrN,"exists. passed.")
            return 0
        scripts=['ldsc.py','--l2','--bfile',bfile_path.format(chrN), '--ld-wind-cm', '1','--out',ld_path.format(annot_name,chrN),'--print-snps',print_snps_path.format(chrN)]
    print(' '.join(scripts))
    sp.call(scripts)    
    print("estimating ld score of ",annot_name,chrN,"finished")
    return 0


def estimate_ldscore(annot_name):
    for chrN in range(1,22+1):
        estimate_ldscore_chrN(chrN)
    print("estimating LD score of",annot_name,"for all chr finished")


def estimate_ldscore_parallel(annot_name,mode='annot'):
    pool = mp.Pool(processes=22)
    if mode=='annot':
        pool.starmap(estimate_ldscore_chrN,[[annot_name,chrN] for chrN in range(1,22+1)])
    elif mode=='univariate':
        pool.starmap(estimate_ldscore_chrN,[[annot_name,chrN,'univariate'] for chrN in range(1,22+1)])
    print("estimating LD score of",annot_name,"finished")


# In[16]:


mode=sys.argv[1]
scale=int(sys.argv[2])


# ```
# jupyter nbconvert 4_make_annot_ldscore.ipynb --to script && python 4_make_annot_ldscore.py 64 bp
# python 4_make_annot_ldscore.py lb 200
# ```

# In[178]:


if scale!=0:
    for chrN in range(1,22+1):
        save_annot_chr_cm(mode+str(scale),chrN,scale=scale,mode=mode)


# In[7]:


if scale==0:
    estimate_ldscore_parallel('uni','univariate')
else:
    estimate_ldscore_parallel(mode+str(scale),'annot')


# ```
# for i in {1..22}; do gzip -d "./lb64.${i}.l2.ldscore.gz"; sed -i -e "s/CHR\tSNP\tBP\tL2/CHR\tSNP\tBP\t${i}.0L2/g" "./lb64.${i}.l2.ldscore"; gzip "./lb64.${i}.l2.ldscore" ; done
# ```

# ```
# for i in {9..9}; do gzip -d "./cm128.${i}.l2.ldscore.gz"; sed -i -e "s/CHR\tSNP\tBP\tL2/CHR\tSNP\tBP\t${i}.0L2/g" "./cm128.${i}.l2.ldscore"; gzip "./cm128.${i}.l2.ldscore" ; done
# ```
