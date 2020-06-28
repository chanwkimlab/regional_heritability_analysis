#!/usr/bin/env python
# coding: utf-8

# ```
# jupyter nbconvert 1_download_gwas_gwasatlas.ipynb --to script
# 
# python 1_download_gwas_gwasatlas.py 0 200
# 
# http://geneatlas.roslin.ed.ac.uk/downloads/?traits=0
# ftp://ftp.igmm.ed.ac.uk/pub/GeneATLAS/
# wget -r ftp://ftp.igmm.ed.ac.uk/pub/GeneATLAS/
# ```

# In[ ]:


#html = urlopen('http://geneatlas.roslin.ed.ac.uk/downloads/?traits=0')
#html_read=html.read()
#soup = BeautifulSoup(html_read)
#table_list=pd.read_html(html_read)


# In[49]:


import sys
import os
from urllib.request import urlopen

#from bs4 import BeautifulSoup
#import pandas as pd


# In[21]:


from ftplib import FTP


# In[45]:


start=int(sys.argv[1])
end=int(sys.argv[2])


# In[27]:


ftp = FTP('ftp.igmm.ed.ac.uk')
ftp.login()
ftp.cwd('pub/GeneATLAS')


# In[42]:


file_list=ftp.nlst()


# In[44]:


file_list_clip=file_list[start:end]


# In[40]:


for file in file_list_clip:
    path_temp='/data01/ch6845/GeneATLAS/data/' + file+'.tmp'
    path='/data01/ch6845/GeneATLAS/data/' + file
    
    if os.path.exists(path):
        print('{} already exists in {}'.format(file,path))
        continue
    else:
        print('{} is being fetched to {}'.format(file,path))
    with open(path_temp, 'wb') as localfile:
        ftp.retrbinary('RETR ' + file, localfile.write);
    os.rename(path_temp,path)
    print('Downloading {} finished'.format(file));

