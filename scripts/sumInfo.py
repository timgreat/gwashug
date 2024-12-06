import os
import sys
import argparse
import scipy.stats

#global paras
gwas1='/gwashug/jobs/trait1/trait1.txt'
gwas2='/gwashug/jobs/trait2/trait2.txt'
name1='trait1'
name2='trait2'
trait1='trait1'
trait2='trait2'
pop1='EUR'
pop2='EUR'



#
def LDSCFormat(gwas,ldsc):
    with open(gwas,'r')as fr,open(ldsc,'w')as fw:
        fr.readline()
        fw.write('SNP\tA1\tA2\tZ\tN\tP\n')
        for line in fr.readlines():
            data=line.strip('\n').split(' ')
            try:
                p_val=float(data[10])
                beta_val=float(data[8])
                mid=scipy.stats.norm.ppf(p_val/2.0)
                z=abs(mid)
                if beta_val<0:
                    z=z*(-1)
                fw.write(data[1]+'\t'+data[5]+'\t'+data[6]+'\t'+str(z)+'\t'+data[4]+'\t'+data[10]+'\n')
            except Exception as e:
                continue


def LDSC(ldsc1,ldsc2):
    path='/gwashug/data/cogwas/'+name1+'_'+name2+'/LDSC'
    if not os.path.exists(path):
        os.system('mkdir '+path)
        os.system('mkdir '+path+'/output')
     
#    ldsc1=path+'/'+name1+'.ldsc'
#    ldsc2=path+'/'+name2+'.ldsc'
#    LDSCFormat(gwas1,ldsc1)
#    LDSCFormat(gwas2,ldsc2)
#    os.system('/anaconda3/envs/ldsc/bin/python /arc/2/PGS/software/ldsc/munge_sumstats.py --sumstats '+ldsc1+' --out '+path+'/'+name1+' --merge-alleles /arc/2/PGS/data/ref/EUR_w_ld_chr/w_hm3.snplist')
#    os.system('/anaconda3/envs/ldsc/bin/python /arc/2/PGS/software/ldsc/munge_sumstats.py --sumstats '+ldsc2+' --out '+path+'/'+name2+' --merge-alleles /arc/2/PGS/data/ref/EUR_w_ld_chr/w_hm3.snplist')
    os.system('/anaconda3/envs/ldsc/bin/python /arc/2/PGS/software/ldsc/ldsc.py --rg '+ldsc1+','+ldsc2+' --ref-ld-chr /arc/2/PGS/data/ref/'+pop1+'_w_ld_chr/ --w-ld-chr /arc/2/PGS/data/ref/'+pop1+'_w_ld_chr/ --out '+path+'/output/'+name1+'_'+name2)
#    os.system('rm -rf '+ldsc1)
#    os.system('rm -rf '+ldsc2)
 #   os.system('rm -rf '+path+'/'+name1+'.log')
 #   os.system('rm -rf '+path+'/'+name2+'.log')
#    os.system('rm -rf '+path+'/'+name1+'.sumstats.gz')
#    os.system('rm -rf '+path+'/'+name2+'.sumstats.gz')

def MR(mrpvalue):
    if not os.path.exists('/gwashug/data/cogwas/'+name1+'_'+name2+'/MR'):
        os.system('mkdir /gwashug/data/cogwas/'+name1+'_'+name2+'/MR')
    with open('/gwashug/data/cogwas/'+name1+'_'+name2+'/MR/summary_info.txt','w')as fw:
        fw.write("GWAS1:\t"+gwas1+"\n")
        fw.write("Pop1:\t"+pop1+"\n")
        fw.write("Trait1:\t"+trait1+"\n")

        fw.write("GWAS2:\t"+gwas2+"\n")
        fw.write("Pop2:\t"+pop2+"\n")
        fw.write("Trait2:\t"+trait2+"\n")

        fw.write("MRPvalue:\t"+mrpvalue+'\n')
    os.system('/anaconda3/envs/gwashug/bin/Rscript /gwashug/scripts/cogwas/trait_level/MR.R --info /gwashug/data/cogwas/'+name1+'_'+name2+'/MR/summary_info.txt')

def REMR(mrpvalue):
    if not os.path.exists('/gwashug/data/cogwas/'+name1+'_'+name2+'/reverse_MR'):
        os.system('mkdir /gwashug/data/cogwas/'+name1+'_'+name2+'/reverse_MR')
    with open('/gwashug/data/cogwas/'+name1+'_'+name2+'/reverse_MR/summary_info.txt','w')as fw:
        fw.write("GWAS1:\t"+gwas2+"\n")
        fw.write("Pop1:\t"+pop2+"\n")
        fw.write("Trait1:\t"+trait2+"\n")

        fw.write("GWAS2:\t"+gwas1+"\n")
        fw.write("Pop2:\t"+pop1+"\n")
        fw.write("Trait2:\t"+trait1+"\n")

        fw.write("MRPvalue:\t"+mrpvalue+'\n')
    os.system('/anaconda3/envs/gwashug/bin/Rscript /gwashug/scripts/cogwas/trait_level/MR.R --info /gwashug/data/cogwas/'+name1+'_'+name2+'/reverse_MR/summary_info.txt')

def GECKO(samestudy):
    if not os.path.exists('/gwashug/data/cogwas/'+name1+'_'+name2+'/GECKO'):
        os.system('mkdir /gwashug/data/cogwas/'+name1+'_'+name2+'/GECKO')
    with open('/gwashug/data/cogwas/'+name1+'_'+name2+'/GECKO/summary_info.txt','w')as fw:
        fw.write("GWAS1:\t"+gwas1+"\n")
        fw.write("Pop1:\t"+pop1+"\n")

        fw.write("GWAS2:\t"+gwas2+"\n")
        fw.write("Pop2:\t"+pop2+"\n")
        
        fw.write("SameStudy:\t"+samestudy+"\n")
    os.system('/anaconda3/envs/gwashug/bin/Rscript /gwashug/scripts/cogwas/trait_level/GECKO.R --info /gwashug/data/cogwas/'+name1+'_'+name2+'/GECKO/summary_info.txt')

def LAVA(type1,type2,case1,case2):
    if not os.path.exists('/gwashug/data/cogwas/'+name1+'_'+name2+'/LAVA'):
        os.system('mkdir /gwashug/data/cogwas/'+name1+'_'+name2+'/LAVA')
    with open('/gwashug/data/cogwas/'+name1+'_'+name2+'/LAVA/summary_info.txt','w')as fw:
        fw.write("GWAS1:\t"+gwas1+"\n")
        fw.write("Pop1:\t"+pop1+"\n")
        fw.write("Trait1:\t"+trait1+"\n")
        fw.write("Type1:\t"+type1+"\n")

        fw.write("GWAS2:\t"+gwas2+"\n")
        fw.write("Pop2:\t"+pop2+"\n")
        fw.write("Trait2:\t"+trait2+"\n")
        fw.write("Type2:\t"+type2+"\n")
        if not (case1=='non'):
            fw.write("Case1:\t"+case1+"\n")
        if not (case2=='non'):
            fw.write("Case2:\t"+case2+"\n")
    os.system('/anaconda3/envs/gwashug/bin/Rscript /gwashug/scripts/cogwas/trait_level/LAVA.R --info /gwashug/data/cogwas/'+name1+'_'+name2+'/LAVA/summary_info.txt')
    
def COLOC(type1,type2,pph4):
    if not os.path.exists('/gwashug/data/cogwas/'+name1+'_'+name2+'/COLOC'):
        os.system('mkdir /gwashug/data/cogwas/'+name1+'_'+name2+'/COLOC')
    with open('/gwashug/data/cogwas/'+name1+'_'+name2+'/COLOC/summary_info.txt','w')as fw:
        fw.write("GWAS1:\t"+gwas1+"\n")
        fw.write("Pop1:\t"+pop1+"\n")
        fw.write("Trait1:\t"+trait1+"\n")
        fw.write("Type1:\t"+type1+"\n")

        fw.write("GWAS2:\t"+gwas2+"\n")
        fw.write("Pop2:\t"+pop2+"\n")
        fw.write("Trait2:\t"+trait2+"\n")
        fw.write("Type2:\t"+type2+"\n")

        fw.write("PPH4:\t"+pph4+"\n")
    os.system('/anaconda3/envs/gwashug/bin/Rscript /gwashug/scripts/cogwas/trait_level/coloc.R --info /gwashug/data/cogwas/'+name1+'_'+name2+'/COLOC/summary_info.txt')

def COLOCEQTL(type1,type2,usemodel1,usemodel2,pph4):
    if not os.path.exists('/gwashug/data/cogwas/'+name1+'_'+name2+'/COLOCEQTL'):
        os.system('mkdir /gwashug/data/cogwas/'+name1+'_'+name2+'/COLOCEQTL')
    with open('/gwashug/data/cogwas/'+name1+'_'+name2+'/COLOCEQTL/summary_info.txt','w')as fw:
        fw.write("GWAS1:\t"+gwas1+"\n")
        fw.write("Type1:\t"+type1+'\n')

        fw.write("GWAS2:\t"+gwas2+'\n')
        fw.write("Type2:\t"+type2+'\n')
        
        fw.write("Use_model1:\t"+usemodel1+'\n')
        fw.write("Use_model2:\t"+usemodel2+'\n')

        fw.write("PPH4:\t"+pph4+'\n')
    os.system('/anaconda3/envs/gwashug/bin/Rscript /gwashug/scripts/cogwas/molecular_level/coloc_eqtl1.R --info /gwashug/data/cogwas/'+name1+'_'+name2+'/COLOCEQTL/summary_info.txt')
    
def SMR(usemodel1,usemodel2):
    if not os.path.exists('/gwashug/data/cogwas/'+name1+'_'+name2+'/SMR'):
        os.system('mkdir /gwashug/data/cogwas/'+name1+'_'+name2+'/SMR')
    with open('/gwashug/data/cogwas/'+name1+'_'+name2+'/SMR/summary_info.txt','w')as fw:
        fw.write("GWAS1:\t"+gwas1+"\n")
        fw.write("GWAS2:\t"+gwas2+'\n')

        fw.write("Use_model1:\t"+usemodel1+'\n')
        fw.write("Use_model2:\t"+usemodel2+'\n')

    os.system('/anaconda3/envs/gwashug/bin/Rscript /gwashug/scripts/cogwas/molecular_level/SMR1.R --info /gwashug/data/cogwas/'+name1+'_'+name2+'/SMR/summary_info.txt')

def INTACT(usemodel1,usemodel2):
    if not os.path.exists('/gwashug/data/cogwas/'+name1+'_'+name2+'/INTACT'):
        os.system('mkdir /gwashug/data/cogwas/'+name1+'_'+name2+'/INTACT')
    if not os.path.exists('/gwashug/data/cogwas/'+name1+'_'+name2+'/INTACT/'+usemodel1+'_'+usemodel2):
        os.system('mkdir /gwashug/data/cogwas/'+name1+'_'+name2+'/INTACT/'+usemodel1+'_'+usemodel2) 
    with open('/gwashug/data/cogwas/'+name1+'_'+name2+'/INTACT/'+usemodel1+'_'+usemodel2+'/summary_info.txt','w')as fw:
        fw.write("GWAS1:\t"+gwas1+"\n")
        fw.write("GWAS2:\t"+gwas2+'\n')

        fw.write("Use_model1:\t"+usemodel1+'\n')
        fw.write("Use_model2:\t"+usemodel2+'\n')

    os.system('/anaconda3/envs/gwashug/bin/Rscript /gwashug/scripts/cogwas/molecular_level/INTACT1.R --info /gwashug/data/cogwas/'+name1+'_'+name2+'/INTACT/'+usemodel1+'_'+usemodel2+'/summary_info.txt')

def checkMRPvalue(ff,mrpvalue):
    pva=float(mrpvalue)
    minPva=1.0
    with open(ff,'r')as fr:
        fr.readline()
        for line in fr.readlines():
            now=float(line.strip('\n').split(' ')[10])
            if  minPva > now:
                minPva=now
    if pva > minPva:
        return mrpvalue
    else:
       return str(minPva)

if __name__=='__main__':
    opt=argparse.ArgumentParser()
    opt.add_argument('--method', type=str,default='ALL',required=True,help='method')
    opt.add_argument('--g1', type=str,required=True,help='gwas1')
    opt.add_argument('--g2', type=str,required=True,help='gwas2')
    opt.add_argument('--p1', type=str,default='EUR',help='pop1')
    opt.add_argument('--p2', type=str,default='EUR',help='pop2')
    opt.add_argument('--r1', type=str,default='trait1',help='trait1')
    opt.add_argument('--r2', type=str,default='trait2',help='trait2')
    opt.add_argument('--m', type=str,default='1E-5',help='MRPvalue')
    opt.add_argument('--s', type=str,default='T',help='SameStudy')
    opt.add_argument('--t1', type=str,default='quant',help='type1')
    opt.add_argument('--t2', type=str,default='quant',help='type2')
    opt.add_argument('--c1', type=str,default='non',help='case1')
    opt.add_argument('--c2', type=str,default='non',help='case2')
    opt.add_argument('--u1', type=str,default='GTEx_Whole_Blood',help='use mode1')
    opt.add_argument('--u2', type=str,default='GTEx_Whole_Blood',help='use model2')
    opt.add_argument('--pp', type=str,default='0.5',help='pph4')
    
    args = opt.parse_args() 
    
    gwas1=args.g1
    gwas2=args.g2
    name1=args.g1.split('/')[-1].replace('.gemma','')
    name2=args.g2.split('/')[-1].replace('.gemma','')
    if not os.path.exists('/gwashug/data/cogwas/'+name1+'_'+name2):
        os.system('mkdir /gwashug/data/cogwas/'+name1+'_'+name2)
    pop1=args.p1
    pop2=args.p2
    trait1=args.r1
    trait2=args.r2
    
    if args.method=='LDSC':
        ldsc1=gwas1.replace('.gemma','.sumstats.gz')
        ldsc2=gwas2.replace('.gemma','.sumstats.gz')
        LDSC(ldsc1,ldsc2)
    elif args.method=='MR':
        mrpvalue=checkMRPvalue(gwas1,args.m)
        MR(mrpvalue)
    elif args.method=='REMR':
        mrpvalue=checkMRPvalue(gwas1,args.m)
        REMR(mrpvalue)
    elif args.method=='GECKO':
        GECKO(args.s)
    elif args.method=='LAVA':
        LAVA(args.t1,args.t2,args.c1,args.c2)
    elif args.method=='COLOC':
        COLOC(args.t1,args.t2,args.pp)    
    elif args.method=='COLOCEQTL':
        COLOCEQTL(args.t1,args.t2,args.u1,args.u2,args.pp)
    elif args.method=='SMR':
        SMR(args.u1,args.u2)
    elif args.method=='INTACT':
        INTACT(args.u1,args.u2)
    else:
        try:
            ldsc1=gwas1.replace('gemma','sumstats.gz')
            ldsc2=gwas2.replace('gemma','sumstats.gz')
            LDSC(ldsc1,ldsc2)
        except Exception as e:
            pass
        try:
            mrpvalue=checkMRPvalue(gwas1,args.m)
            MR(mrpvalue)
        except Exception as e:
            pass
        try:
            GECKO(args.s)
        except Exception as e:
            pass
        try:
            LAVA(args.t1,args.t2,args.c1,args.c2)
        except Exception as e:
            pass
        try:
            COLOC(args.t1,args.t2,args.pp)
        except Exception as e:
            pass
        try:
            COLOCEQTL(args.t1,args.t2,args.u1,args.u2,args.pp)
        except Exception as e:
            pass
        try:
            SMR(args.u1,args.u2)
        except Exception as e:
            pass
        try:
            INTACT(args.u1,args.u2)
        except Exception as e:
            pass
