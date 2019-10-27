#!/usr/bin/env python
# coding: utf-8

# In[1]:


import sys
import os
import glob
import subprocess as sp
import multiprocessing as mp
import pandas as pd
import numpy as np

from basic_tools import *

debug=False


# In[2]:


def run_ldsc(pheno_code,ld,output,mode='original',samp_prev=np.nan,pop_prev=np.nan):
    if os.path.exists(ldsc_path.format(pheno_code)+'.log'):
        print("Congratulations!. ldsc result of",pheno_code,"exists. passed.")
        return
    if mode=='original':
        script=['ldsc.py','--h2',sumstats_path.format(pheno_code)+'.sumstats.gz', 
                 '--ref-ld-chr',ld_path.format(ld,''),
                 '--w-ld-chr',wld_path,
                 '--out',ldsc_path.format(output)]
    elif mode=='my':
        script=['ldsc_my.py','--h2',sumstats_path.format(pheno_code)+'.sumstats.gz', 
         '--ref-ld-chr',ld_path.format(ld,''),
         '--w-ld-chr',wld_path,
         '--out',ldsc_path.format(output)] 
    else:
        print("run_ldsc mode Error!!!!!!!")
    
    if np.isnan(samp_prev)==False and np.isnan(pop_prev)==False:
        script+=['--samp-prev',str(samp_prev),'--pop-prev',str(pop_prev)]
    
    print('Started:',' '.join(script))
    sp.call(script) 
    print('Finished:',' '.join(script))


# In[3]:


def run_ldsc_wrapper(prefix,scale,pheno_code,samp_prev=np.nan,pop_prev=np.nan):
    run_ldsc(pheno_code,prefix,'{}.{}'.format(prefix,pheno_code),mode='original' if mode=='uni' else 'my',samp_prev=samp_prev,pop_prev=pop_prev)
        


# In[3]:


sys.argv#uni 0 20 x x


# In[4]:


mode=sys.argv[1]
scale=int(sys.argv[2])
cores=int(sys.argv[3])
start=int(sys.argv[4])
end=int(sys.argv[5])

if mode=='uni':
    prefix=mode
else:
    prefix=mode+str(scale)


# In[25]:


#start,end,prefix=0,1000,'bp300'


# In[3]:


phenotypes_uni_filtered['prevalence']=phenotypes_uni_filtered['n_cases']/phenotypes_uni_filtered['n_non_missing']


# In[4]:


phenotypes_uni_filtered.shape


# In[5]:


pheno_code_list_todo=[]

for idx,row in phenotypes_uni_filtered.iloc[start:end].iterrows():
    if os.path.exists(ldsc_path.format('{}.{}'.format(prefix,idx))+'.log'):
        #print(ldsc_path.format('{}.{}'.format(prefix,idx))+'.log','exists')
        continue
    print(idx,end=' ')
    pheno_code_list_todo.append((idx,row['prevalence']))


# In[11]:


"""
phenotypes_filtered['prevalence']=phenotypes_filtered['n_cases']/phenotypes_filtered['n_non_missing']

phenotypes_filtered.shape

pheno_code_list_todo=[]
for idx,row in phenotypes_filtered.iloc[start:end].iterrows():
    if os.path.exists(ldsc_path.format('{}.{}'.format(prefix,idx))+'.log'):
        continue
    print(idx,end=' ')
    pheno_code_list_todo.append((idx,row['prevalence']))
"""


# ```
# jupyter nbconvert 5_run_ldsc.ipynb --to script
# 
# export SCREENDIR=$HOME/.screen
# 
# start=0;end=600;mode=uni
# python 5_run_ldsc.py $mode 0 10 $start $end 
# 
# start=0;end=600;mode=bp
# python 5_run_ldsc.py $mode 300 10 $start $end && python 5_run_ldsc.py $mode 128 10 $start $end && python 5_run_ldsc.py $mode 64 5 $start $end && python 5_run_ldsc.py $mode 32 5 $start $end && python 5_run_ldsc.py $mode 16 5 $start $end && python 5_run_ldsc.py $mode 8 2 $start $end
# 
# 
# python 5_run_ldsc.py bp 128 10 $start $end && python 5_run_ldsc.py bp 64 5 $start $end && python 5_run_ldsc.py bp 32 5 $start $end && python 5_run_ldsc.py bp 16 5 $start $end && python 5_run_ldsc.py bp 8 2 $start $end
# ```

# In[ ]:


#pool = mp.Pool(processes=15)
#pool.starmap(run_ldsc_wrapper,[(mode,scale,pheno_code,prevelence,prevelence) for (pheno_code,prevelence) in pheno_code_list_todo])


# In[ ]:


pool = mp.Pool(processes=cores)
#pool.starmap(run_ldsc_wrapper,[(mode,scale,pheno_code) for pheno_code in pheno_code_list_todo])
pool.starmap(run_ldsc_wrapper,[(prefix,scale,pheno_code,prevelence,prevelence) for (pheno_code,prevelence) in pheno_code_list_todo])

