import math
import sys
import os
import gzip
import scipy

def process(file):
    sepN='\t'
    work_dir="/home/website/gwashug-server/job/TEST/"
    data_dir="/home/website/gwashug-server/job/TEST/trait1.txt.gz"
    rsid="rsID"
    beta="BETA"
    p="P"
    size="20000"
    pop="EUR"
    with open(file,'r')as fr:
        for line in fr.readlines():
            inf=line.strip('\n').split()
            if len(inf)<2:
                continue
            if inf[0]=='Separator':
                sepN=inf[1]
                if sepN=='tab':
                    sep='\t'
                elif sepN=='space':
                    sep=' '
                else:
                    sep=','
            elif inf[0]=='Work_dir':
                work_dir=inf[1]
            elif inf[0]=='Data_file':
                data_file=inf[1]
            elif inf[0]=='Rsid_alias':
                rsid=inf[1]
            elif inf[0]=='Beta_alias':
                beta=inf[1]
            elif inf[0]=='PValue_alias':
                p=inf[1]
            elif inf[0]=='Sample_Size_alias':
                size=inf[1]
            elif inf[0]=='Ancestry':
                pop=inf[1]
    snpinfo = {}
    with gzip.open(data_file,'r') as fr:
        name=fr.readline().decode().strip('\n').split(sep)
        rsid_col=name.index(rsid)
        beta_col=name.index(beta)
        p_col=name.index(p)
        for line in fr.readlines():
            snp=line.decode().strip('\n').split(sep)
            id=snp[rsid_col]
            snpinfo.update({id:[snp[beta_col],snp[p_col]]})
    if os.path.exists('/home/preprocess_gwas/'+pop+'.map'):
        fr=open('/home/preprocess_gwas/'+pop+'.map','r')
    else:
        fr=open('/home/preprocess_gwas/EUR.map','r')
    fw=open(work_dir+'/trait1.txt','w')
    fw.write('chr rs ps n_mis n_obs allele1 allele0 af beta se p_wald\n')
    for line in fr.readlines():
        data1=line.strip('\n')
        data=data1.split('\t')
        id=data[2]
        if data[2] in snpinfo:
            try:
                p_val=float(snpinfo[id][1])
                beta_val=float(snpinfo[id][0])
                if(p_val>=0.0 and p_val<1.0):
                    mid=scipy.stats.norm.ppf(p_val/2.0)
                    #z=abs(mid)
                    #if beta_val<0:
                    #    z=z*(-1)
                    se=abs(beta_val/mid)
                    if not math.isclose(se,0,rel_tol=0,abs_tol=1E-6):
                        fw.write(data[0]+' '+data[2]+' '+data[1]+' 0 '+size+' '+data[4]+' '+data[5]+' '+data[3]+' '+str(beta_val)+' '+str(se)+' '+str(p_val)+'\n')
            except Exception as e:
                continue
    fr.close()
    fw.close()
    #os.system('gzip -cf '+work_dir+'/trait1.txt > '+work_dir+'/trait1.txt.gz')
def processX(file):
    sepN='\t'
    work_dir="/home/website/gwashug-server/job/TEST/"
    data_dir="/home/website/gwashug-server/job/TEST/trait2.txt.gz"
    rsid="rsID"
    beta="BETA"
    p="P"
    size="20000"
    pop="EUR"
    with open(file,'r')as fr:
        for line in fr.readlines():
            inf=line.strip('\n').split()
            if len(inf)<2:
                continue
            if inf[0]=='Separator':
                sepN=inf[1]
                if sepN=='tab':
                    sep='\t'
                elif sepN=='space':
                    sep=' '
                else:
                    sep=','
            elif inf[0]=='Work_dir':
                work_dir=inf[1]
            elif inf[0]=='Data_file':
                data_file=inf[1]
            elif inf[0]=='Rsid_alias':
                rsid=inf[1]
            elif inf[0]=='Beta_alias':
                beta=inf[1]
            elif inf[0]=='PValue_alias':
                p=inf[1]
            elif inf[0]=='Sample_Size_alias':
                size=inf[1]
            elif inf[0]=='Ancestry':
                pop=inf[1]
    snpinfo = {}
    with gzip.open(data_file,'r') as fr:
        name=fr.readline().decode().strip('\n').split(sep)
        rsid_col=name.index(rsid)
        beta_col=name.index(beta)
        p_col=name.index(p)
        for line in fr.readlines():
            snp=line.decode().strip('\n').split(sep)
            id=snp[rsid_col]
            snpinfo.update({id:[snp[beta_col],snp[p_col]]})
    if os.path.exists('/home/preprocess_gwas/'+pop+'.map'):
        fr=open('/home/preprocess_gwas/'+pop+'.map','r')
    else:
        fr=open('/home/preprocess_gwas/EUR.map','r')
    fw=open(work_dir+'/trait2.txt','w')
    fw.write('chr rs ps n_mis n_obs allele1 allele0 af beta se p_wald\n')
    for line in fr.readlines():
        data1=line.strip('\n')
        data=data1.split('\t')
        id=data[2]
        if data[2] in snpinfo:
            try:
                p_val=float(snpinfo[id][1])
                beta_val=float(snpinfo[id][0])
                if(p_val>=0.0 and p_val<1.0):
                    mid=scipy.stats.norm.ppf(p_val/2.0)
                    #z=abs(mid)
                    #if beta_val<0:
                    #    z=z*(-1)
                    se=abs(beta_val/mid)
                    if not math.isclose(se,0,rel_tol=0,abs_tol=1E-6):
                        fw.write(data[0]+' '+data[2]+' '+data[1]+' 0 '+size+' '+data[4]+' '+data[5]+' '+data[3]+' '+str(beta_val)+' '+str(se)+' '+str(p_val)+'\n')
            except Exception as e:
                continue
    fr.close()
    fw.close()
    #os.system('gzip -cf '+work_dir+'/trait2.txt > '+work_dir+'/trait2.txt.gz')
if __name__=='__main__':
    path=sys.argv[1]
    files=os.listdir(path)
    flag=False
    flagx=False
    for file in files:
        if '.param' in file:
            paramPath=path+'/'+file
            print(paramPath)
            flag=True
        if '.xparam' in file:
            xparamPath=path+'/'+file
            print(xparamPath)
            flagx=True
        if flag and flagx:
            break
    if flag:
        process(paramPath)
    else:
        print('no param file!!!')
    if flagx:
        processX(xparamPath)

