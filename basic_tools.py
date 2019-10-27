import subprocess
import pandas as pd
import numpy as np
import glob
import log_parser

from path_configure import *

def run_command(command,quiet=False):
    print("------{}-----".format("RUN"))
    print(command)
    process = subprocess.Popen(command.split(), stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    stdout, stderr = process.communicate()
    if quiet==False:
        print("------{}-----".format("ERROR"))
        print(stderr.decode())
        print("------{}-----".format("OUTPUT"))
        print(stdout.decode())
    return stdout,stderr


def code_to_filename(phenotype_code):
    return phenotypes_filtered.loc[phenotype_code]['File']
#[phenotypes_filtered['phenotype']==phenotype_code].iloc[0]['File']
    """
    ukbb_table_select=ukbb_table_filtered[ukbb_table_filtered['Phenotype Code'].str.replace('_irnt','')==phenotype_code]
    if ukbb_table_select.shape[0]==0:
        print("phenotype code does not exist in ukbb_table")
    elif ukbb_table_select.shape[0]==2:
        print("duplicated phenotype code")
    elif ukbb_table_select.shape[0]==1:
        return ukbb_table_select['File'].values[0]
    else:
        print("code to path conversion error.")
    """
    
def code_to_description(phenotype_code):
    return phenotypes_filtered.loc[phenotype_code]['description']
    """
    ukbb_table_select=ukbb_table_filtered[ukbb_table_filtered['Phenotype Code'].str.replace('_irnt','')==phenotype_code]
    return ukbb_table_select.iloc[0]['Phenotype Description']
    #phenotypes_select=phenotypes[phenotypes['phenotype'].str.contains(phenotype_code)]
    #return phenotypes_select.iloc[0]['description']
    """
def read_gwas(phenotype_code):
    filename=code_to_filename(phenotype_code)
    gwas_result=pd.read_csv(gwas_path+filename,sep='\t',compression='gzip')
    print("Loading gwas result(not munged) finished",gwas_path+filename)
    return gwas_result

def read_gwas2(phenotype_code):
    load_variables()
    print("read_gwas2({}) was called (munge version of read_gwas). It returns subset of gwas result".format(phenotype_code))
    gwas_result=read_gwas(phenotype_code)#.loc[variants.index,:]
    gwas_result['chr']=variants['chr']
    gwas_result=gwas_result[gwas_result['chr']!='X']
    gwas_result['chr'] = gwas_result['chr'].astype(int)
    
    gwas_result['pos']=variants['pos']
    gwas_result['ref']=variants['ref']
    gwas_result['alt']=variants['alt']
    gwas_result['rsid']=variants['rsid']
    gwas_result['varid']=variants['varid']
    #gwas_result['cm']=variants['cm']
    gwas_result['info']=variants['info']
    gwas_result['minor_allele']=variants['minor_allele']
    gwas_result['minor_AF']=variants['minor_AF']
    #gwas_result=gwas_result.loc[index_filter,:]
    print("Munging loaded gwas result finished.",phenotype_code)
    cname_dict={'chr':'CHR','pos':'BP','rsid':'SNP','cm':'CM','ref':'A1','alt':'A2','pval':'P','beta':'BETA','minor_AF':'MAF','info':'INFO','n_complete_samples':'N'}
    gwas_result.rename(columns=cname_dict,inplace=True)
    return gwas_result


def index_format(index,mode='minimal',compensate=0):
    index_new=[]    
    for col in index:
        col_dot=col.split('.')
        col_dot_L2=col_dot[3].split("L2")

        sub_cat="{:03d}".format(int(col_dot_L2[0])+compensate)
        col_dot[3]='L2'.join([sub_cat]+col_dot_L2[1:])
        if mode=='minimal':
            a='.'.join(["{:02d}".format(int(col_dot[2]))]+[sub_cat])
        elif mode=='chr':
            a="{:02d}".format(int(col_dot[2]))
        else:
            a='.'.join(col_dot[:2]+["{:02d}".format(int(col_dot[2]))]+col_dot[3:])
        index_new.append(a)
    return index_new

"""
def minimal_to_chr(index):
    index_new=[]
    for col in index:
        col_dot=col.split('.')
        a="{:02d}".format(int(col_dot[0]))
        index_new.append(a)
    return index_new
"""
"""
def cm_minimal_to_float(index,compensate=0.0):
    index_new=[]
    max_cm={1:294,10:184,11:162,12:176,13:130,14:117,15:151,16:132,17:129,18:121,19:107,2:275,20:111,21:65,22:76,3:228,4:221,5:209,6:199,7:191,8:179,9:181}
    max_cm=[max_cm[i] for i in range(1,22+1)]
    for col in index:
        #col_dot=col.split('.')
        #a="{:02d}".format(int(col_dot[0]))
        col=float(col)
        col_adjusted=int(col)+1000*(col-int(col)+compensate)/(max_cm[int(col)-1])
        index_new.append(col_adjusted)
        #print(col,col_adjusted)
    return index_new

def float_to_chrcm(list):
    array=np.array(list)
    chr_array=np.floor(array).astype(int)
    scaled_cm_array=array-chr_array
    print(scaled_cm_array)
    
    max_cm={1:294,10:184,11:162,12:176,13:130,14:117,15:151,16:132,17:129,18:121,19:107,2:275,20:111,21:65,22:76,3:228,4:221,5:209,6:199,7:191,8:179,9:181}
    max_cm=[max_cm[i] for i in range(1,22+1)]
    
    max_cm_array=scaled_cm_array*pd.Series(max_cm).loc[chr_array-1]
    cm_array=max_cm_array
    df=pd.DataFrame([chr_array,cm_array],index=['CHR','CM']).T
    #sumstats['CMcoord']=chr_array+sumstats['CM']/(pd.Series(max_cm).loc[chr_array-1]).values
    return df

#float_to_chrcm([1.5,20.5])
"""

def search_col(columns,queries,exact=False):
    cnames=[]
    
    for col in columns:
        #print(col,"a",query in col)
        check=False
        for query in queries:
            if query in col and (col.find(query)==0 if exact else True):
                check=True
        if check:
            cnames.append(col)
    return cnames

def search_df(df,column_name,queries):
    series=df[column_name]
    for item in series:
        print(item)
    #print(series)
    pass

def nearest_points(x,y,refx,refy,rank,index=None):
    x=np.array(x)
    y=np.array(y)
    relative_dist=np.sqrt(np.power((x-refx)/np.nanmax(np.abs(x)),2)+np.power((y-refy)/np.nanmax(np.abs(y)),2))
    df=pd.DataFrame([x,y,relative_dist],index=['x','y','relative_dist'],columns=index).T
    return df.sort_values('relative_dist').head(rank)


def read_sumstats(pheno_code,mhc=False,quiet=False):
    load_variables2()
    sumstats=pd.read_csv(sumstats_path.format(pheno_code)+'.sumstats.gz',sep='\t',compression='gzip')
    sumstats=sumstats[np.logical_not(np.isnan(sumstats['Z']))]
    if quiet==False:
        print("SNP in sumstats before merging",sumstats.shape[0])
    sumstats_merged=pd.merge(sumstats,refld_stack,how='inner',on='SNP')
    if quiet==False:
        print("SNP in sumstats after merging with refld",sumstats_merged.shape[0])
    if mhc:
        sumstats_merged=pd.merge(sumstats_merged,wldmhc_stack,how='inner',on='SNP')
    else:
        sumstats_merged=pd.merge(sumstats_merged,wld_stack,how='inner',on='SNP')
    if quiet==False:
        print("SNP in sumstats after merging with wld",sumstats_merged.shape[0])
    sumstats_merged=pd.merge(sumstats_merged,cm_map,how='inner',on='SNP')
    if quiet==False:
        print("SNP in sumstats after merging with cm_map",sumstats_merged.shape[0])
    final=sumstats_merged[['SNP','Z','refld.L2','CM','CHR','BP','N']]
    final['chisq']=final['Z']**2
    del sumstats
    del sumstats_merged
    return final#.sort_values(['CHR','BP'])

def sumstats_filter(sumstats,option={}):
    chr_filter=option.get('chrN',[i for i in range(1,22+1)])
    
    cm_min=option.get('cm_min',None)
    cm_max=option.get('cm_max',None)
    if len(chr_filter)<22:
        sumstats=sumstats[np.logical_or.reduce([sumstats.CHR==chrN for chrN in chr_filter])]
    if cm_min!=None:
        sumstats=sumstats[sumstats.CM>=cm_min]
    if cm_max!=None:
        sumstats=sumstats[sumstats.CM<cm_max]
        
    return sumstats
#sumstats=read    


def split_se(parsed,cname):
    parsed[[cname,cname+'_se']]=parsed[cname].str.split("(",expand=True)
    parsed[[cname+'_se','temp']]=parsed[cname+'_se'].str.split(")",n=1,expand=True)
    parsed[cname]=parsed[cname].str.strip().astype(float)
    parsed[cname+'_se']=parsed[cname+'_se'].str.strip().astype(float)
    parsed.drop(columns='temp',inplace=True)
    cname_index=np.flatnonzero(parsed.columns==cname)[0]
    cols = parsed.columns.tolist()
    cols = cols[:cname_index+1]+cols[-1:] +cols[cname_index+1:-1] 
    parsed=parsed[cols]
    return parsed



def read_ldsc(suffix="cm1",verbose=True):

    
    log_files=glob.glob(ldsc_path.format('{}.{}'.format(suffix,'*'))+'.log')
    if verbose:
        parse_list=[\
        ('Total Observed scale h2',0),('Total Liability scale h2',0),('Lambda GC',0),('Mean Chi^2',0),('Intercept',0),('Ratio',0),\
        ('Categories',1),('Observed scale h2',1),('Observed scale h2 SE',1),('Liability scale h2',1),('Liability scale h2 SE',1),\
        ('Proportion of SNPs',1),('Proportion of h2g',1),('Enrichment',1),\
        ('Coefficients',1),('Coefficient SE',1),\
        ]
    else:
        parse_list=[\
        ('Total Observed scale h2',0),('Total Liability scale h2',0),('Lambda GC',0),('Mean Chi^2',0),('Intercept',0),('Ratio',0),\
        ]

    parsed=log_parser.file_todf(log_files,parse_list,dtype=float)
    parsed['phenotype']=pd.DataFrame(pd.DataFrame(parsed['filename'].str.split('/',expand=False).tolist()).iloc[:,-1].str.split('.').to_list())[1]

    parsed=split_se(parsed,'Intercept')
    parsed=split_se(parsed,'Total Observed scale h2')
    parsed=split_se(parsed,'Total Liability scale h2')
    parsed=parsed.set_index('phenotype')
    h2_cm_sorted=parsed.sort_values('Total Observed scale h2',ascending=False)
    h2_cm_sorted.columns = h2_cm_sorted.columns.map(lambda x: suffix+'.'+str(x))
    return h2_cm_sorted

def parse_uni_regression_result(h2_total,pheno_code):
    h2=h2_total['uni.Total Observed scale h2'][pheno_code]
    h2_se=h2_total['uni.Total Observed scale h2_se'][pheno_code]
    intercept=h2_total['uni.Intercept'][pheno_code]
    
    m_5_50=pd.concat([pd.read_csv(m_5_50_chr,sep='\t',header=None) for m_5_50_chr in glob.glob(ld_path.format('uni','*')+".l2.M_5_50")])
    
    regression_result_uni={"h2":h2,"h2_se":h2_se,"intercept":intercept}
    regression_result_uni['category']='uni'
    regression_result_uni['M']=m_5_50.sum().values[0]
    
    print("finished loading","uni")
    return regression_result_uni

def parse_par_regression_result(h2_total,pheno_code,suffix,quiet=False):    
    categories=search_col(h2_total.columns,['{}.Categories'.format(suffix)])
    categories_strip=[category[category.index('Categories')+11:] for category in categories]
    categories_strip=[category.replace('L2_0','') for category in categories_strip]

    regression_result_par_list=[]

    for category_strip in categories_strip:
        coef=h2_total['{}.Coefficients.{}L2_0'.format(suffix,category_strip)][pheno_code]
        coef_se=h2_total['{}.Coefficient SE.{}L2_0'.format(suffix,category_strip)][pheno_code]
        intercept=h2_total['{}.Intercept'.format(suffix)][pheno_code]
        h2=h2_total['{}.Observed scale h2.{}L2_0'.format(suffix,category_strip)][pheno_code]
        h2_se=h2_total['{}.Observed scale h2 SE.{}L2_0'.format(suffix,category_strip)][pheno_code]

        regression_result_par={"coef":coef,"coef_se":coef_se,"h2":h2,"h2_se":h2_se,"intercept":intercept}
        regression_result_par['category']=suffix+'.'+category_strip
        regression_result_par['M']=h2/coef
            
        regression_result_par_list.append(regression_result_par)
        if not quiet:
            print("finished loading",category_strip)

    regression_result_par_df=pd.DataFrame(regression_result_par_list)
    regression_result_par_df['category']=pd.Index(regression_result_par_df['category'])
    
    return regression_result_par_df

def make_regression_result_list(h2_total,pheno_code,suffix_list=['cm300','cm128','cm64','cm32','cm16','cm8','cm4']):
    regression_result_list=[]
    for suffix in suffix_list:
        regression_result_par_df=parse_par_regression_result(h2_total,pheno_code,suffix,quiet=True)
        regression_result_list.append(regression_result_par_df)
    return pheno_code,regression_result_list


def description_to_short(desc,mode='pub',suffix="basic"):

    description_dict=description_dict_merge
    
    
    if not mode in ['ori','pub','abbr','abbr_pub','abbr_pub_ori']:
        raise    
    if desc in description_dict.index:
        if mode=='ori':
            return desc
        elif mode=='pub':
            #print(correlation_description_dict['description_pub'].loc['Tobacco smoking: Ex-smoker'],description_dict['description_pub'].loc['Tobacco smoking: Ex-smoker'])
            return description_dict['description_pub'+"_"+suffix].loc[desc]
        elif mode=='abbr':
            return description_dict['description_abbr'+"_"+suffix].loc[desc]
        elif mode=='abbr_pub':
            abbr=description_to_short(desc,mode='abbr',suffix=suffix)
            return description_to_short(desc,mode='pub',suffix=suffix) if type(abbr)==float or abbr=='' else abbr
        elif mode=='abbr_pub_ori':
            abbr_pub=description_to_short(desc,mode='abbr_pub',suffix=suffix)
            #print(abbr_pub)
            return desc if type(abbr_pub)==float or abbr_pub=='' else abbr_pub
    else:
        return np.nan

"""
def description_to_short(desc,mode='pub',data=0):
    if data==0:
        description_dict=phenotypes_par_filtered_description_dict
    elif data==1:
        description_dict=pleiotropic_loci_description_dict
    elif data==2:
        description_dict=correlation_description_dict
        
    if not mode in ['ori','pub','abbr','abbr_pub','abbr_pub_ori']:
        raise    
    if desc in description_dict.index:
        if mode=='ori':
            return desc
        elif mode=='pub':
            #print(correlation_description_dict['description_pub'].loc['Tobacco smoking: Ex-smoker'],description_dict['description_pub'].loc['Tobacco smoking: Ex-smoker'])
            return description_dict['description_pub'].loc[desc]
        elif mode=='abbr':
            return description_dict['description_abbr'].loc[desc]
        elif mode=='abbr_pub':
            abbr=description_to_short(desc,mode='abbr',data=data)
            return description_to_short(desc,mode='pub',data=data) if type(abbr)==float or abbr=='' else abbr
        elif mode=='abbr_pub_ori':
            abbr_pub=description_to_short(desc,mode='abbr_pub',data=data)
            #print(abbr_pub)
            return desc if type(abbr_pub)==float or abbr_pub=='' else abbr_pub
    else:
        return np.nan
"""

def category_to_format(category,mode):
    scale,chrN,start=category.split('.')
    chrN=int(chrN);start=int(start)
    end=start+int(scale.replace('bp',''))
    if type(mode)==int:
        scale_int=int(scale.replace('bp',''))
        return '{}.{}.{}'.format(scale,chrN,start-scale_int*((start//scale_int)%2**mode))
    elif mode=='padding':
        return '{}.{:02d}.{:03d}'.format(scale,chrN,start)    
    elif mode=='scale':
        scale_int=int(scale.replace('bp',''))
        return scale_int
    elif mode=='chr':
        return chrN
    elif mode=='bp_range':
        load_variables3()
        end=chr_bp_max[chrN-1]/1000000 if end>chr_bp_max[chrN-1]/1000000 else end
        return end-start    
    elif mode=='chr_bp':
        load_variables3()
        end=chr_bp_max[chrN-1]/1000000 if end>chr_bp_max[chrN-1]/1000000 else end
        return 'chr{}: {}~{:d}Mb'.format(chrN,start,int(end))
    else:
        raise
        
try:
    phenotypes_filtered
    phenotypes_uni_filtered
    phenotypes_par_filtered
    """
    phenotypes_par_filtered_description_dict
    pleiotropic_loci_description_dict
    corelation_description_dict
    """
    description_dict_merge
except:
    phenotypes_filtered=pd.read_csv(phenotypes_filtered_path,index_col='phenotype')
    phenotypes_uni_filtered=pd.read_csv(phenotypes_uni_filtered_path,index_col='phenotype')
    phenotypes_par_filtered=pd.read_csv(phenotypes_par_filtered_path,index_col='phenotype')
    """
    phenotypes_par_filtered_description_dict=pd.read_csv(phenotypes_par_filtered_description_dict_path,sep='\t',index_col='description')
    pleiotropic_loci_description_dict=pd.read_csv(pleiotropic_loci_description_dict_path,sep='\t',index_col='description',converters={'description_pub':lambda x:x.replace('\\n','\n')})
    correlation_description_dict=pd.read_csv(correlation_description_dict_path,sep='\t',index_col='description',converters={'description_pub':lambda x:x.replace('\\n','\n')})
    """
    description_dict_merge=pd.read_csv(description_dict_merge_path,sep='\t',index_col='description')
"""
try:
	ukbb_table_filtered
	phenotypes
except:
	ukbb_table_filtered=pd.read_csv(ukbb_table_filtered_path)
	ukbb_table_filtered['phenotype']=ukbb_table_filtered['Phenotype Code']

	phenotypes_both_sexes=pd.read_csv(phenotypes_both_sexes_file_path,sep='\t',compression='gzip');phenotypes_both_sexes['Sex']='both_sexes'
	phenotypes_male=pd.read_csv(phenotypes_male_file_path,sep='\t',compression='gzip');phenotypes_male['Sex']='male'
	phenotypes_female=pd.read_csv(phenotypes_female_file_path,sep='\t',compression='gzip');phenotypes_female['Sex']='female'
    
	biomarkers_both_sexes=pd.read_csv(biomarkers_both_sexes_file_path,sep='\t',compression='gzip');biomarkers_both_sexes['Sex']='both_sexes'
	biomarkers_male=pd.read_csv(biomarkers_male_file_path,sep='\t',compression='gzip');biomarkers_male['Sex']='male'
	biomarkers_female=pd.read_csv(biomarkers_female_file_path,sep='\t',compression='gzip');biomarkers_female['Sex']='female'
	phenotypes=pd.concat([phenotypes_both_sexes,phenotypes_male,phenotypes_female,biomarkers_both_sexes,biomarkers_male,biomarkers_female],sort='False')

	phenotypes_both_sexes_v2=pd.read_csv(phenotypes_both_sexes_v2_file_path,sep='\t',compression='gzip')
	phenotypes_male_v2=pd.read_csv(phenotypes_male_v2_file_path,sep='\t',compression='gzip')
	phenotypes_female_v2=pd.read_csv(phenotypes_female_v2_file_path,sep='\t',compression='gzip')
	phenotypes_both_sexes_v2['Sex']='both_sexes';phenotypes_male_v2['Sex']='male';phenotypes_female_v2['Sex']='female';
	phenotypes=pd.concat([phenotypes_both_sexes_v2,phenotypes_male_v2,phenotypes_female_v2],sort='False')
	#phenotypes=pd.concat([phenotypes,phenotypes_v2],sort='False')
	
	ukbb_table_filtered=ukbb_table_filtered.merge(phenotypes,how='left',on=['phenotype','Sex'])
"""        


        
def load_variables():
    global variants
    global index_filter
    try:
        variants
        #index_filter
    except:
        print("found that variables 'variants' was not loaded. trying to load the variables")
        variants=pd.read_csv(variants_file_path,sep='\t',compression='gzip')
        #cm_map_=pd.read_csv(variants_w_cm_path,sep='\t',header=None,names=['chr','rsid','cm','pos','ref','alt'])
        #cm_map=pd.read_csv(variants_w_cm_path,sep='\t',header=None,names=['CHR','SNP','CM','BP','map.ref','map.alt'])
        
        # pre-calculated 'filter.index' must be made already one time.
        #index_filter=pd.read_csv(filter_index_path,header=None).iloc[:,0].values

        #variants_['cm']=cm_map_['cm']
        #variants=variants_.loc[index_filter,:]
        #variants['chr'] = variants['chr'].astype(int)

        #del variants_
        #del cm_map_

def load_variables2():
    global cm_map
    global wld_stack
    global refld_stack
    try:
        pass
        #cm_map
        #wld_stack
        #refld_stack
    except:
        print("found that variables 'cm_map', 'w_ld_stack','refld_stack' were not loaded. trying to load the variables")
        cm_map=pd.read_csv(variants_w_cm_path,sep='\t',header=None,names=['CHR','SNP','CM','BP','map.ref','map.alt'])
        print("SNPs in cm_map",cm_map.shape[0])
    
        wld_stack=pd.concat([pd.read_csv(wld,sep='\t',compression='gzip') for wld in glob.glob(wld_path+"*")])
        wld_stack.columns = wld_stack.columns.map(lambda x: 'wld.'+str(x))
        wld_stack.rename(columns={'wld.SNP':'SNP'},inplace=True)
        print("SNPs in wld",wld_stack.shape[0])
        """
        wldmhc_stack=pd.concat([pd.read_csv(wld,sep='\t',compression='gzip') for wld in glob.glob(wldmhc_path+"*ldscore.gz")])
        wldmhc_stack.columns = wldmhc_stack.columns.map(lambda x: 'wld.'+str(x))
        wldmhc_stack.rename(columns={'wld.SNP':'SNP'},inplace=True)
        print("SNPs in wldmhc",wldmhc_stack.shape[0])
        """
        refld_stack=pd.concat([pd.read_csv(refld,sep='\t',compression='gzip') for refld in glob.glob(ld_path.format('uni','')+"*"+"ldscore.gz")])
        refld_stack.columns = refld_stack.columns.map(lambda x: 'refld.'+str(x))
        refld_stack.rename(columns={'refld.SNP':'SNP'},inplace=True)
        print("SNPs in refld",refld_stack.shape[0])

        
def load_variables3():
    global gwas_result_sample
    global chr_bp_max
    try:
        gwas_result_sample
        chr_bp_max
    except:
        print("found that variables 'gwas_result_sample','chr_bp_max' were not loaded. trying to load the variables")
        gwas_result_sample=pd.read_pickle('23115.pickle')
        chr_bp_max=[gwas_result_sample[gwas_result_sample.CHR==chrN].BP.iloc[-1] for chrN in range(1,22+1)]
        print("gwas_result_sample.shape:",gwas_result_sample.shape)
        