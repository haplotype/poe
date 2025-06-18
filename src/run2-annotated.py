import os
import subprocess

for chr in list(range(1,23))[::-1]:
    fi = "wgs4rnaseq.chr"+str(chr)+".2snps.vcf.gz"
    fo = "filter-pass.chr"+str(chr)+".vcf"
    cmd = "bcftools view -f 'PASS' " + fi + " > " + fo  
    print(cmd) 
#    os.system(cmd)            
##the PASS is all you need. 

    f2 = "filter-pass-1477.chr"+str(chr)
    cmd = "/home/guany4/work/poe/poe  -f fa-mo-son.1477 -i " + fo + " -o " + f2  + " " 
    print(cmd) 
#    os.system(cmd)
##here poe is a suite of tools developed to aid the parent-of-origin analysis using idul. 
##poe will be available with documentation at github.com/haplotype/poe 
##in this command -f takes a trio file (ID only) and vcf to extract child genotype in the order specified in the trio file. 
##the child genotypes are partially phased by Mendelian inheritance, and the triple heterozygous sites are marked as missing. 
##output are in bimbam format. 



    gt=f2+".bimbam.gt"
    pos=f2+".bimbam.pos" 
    out=f2+'.mt'
    cmd = "/home/guany4/work/elai/elai -g " + gt + " -p 1 -exclude-maf 0.01 -C 2 -c 10 -s 20 -w 0 -wmg -weg -pos " + pos + " -o " + out + " -e 5 -nthreads 64" 
    print(cmd) 
#    os.system(cmd)
##this runs multi-threading elai, which can be found at github.com/haplotype/elai. 
##this special elai assume genotypes are partially phased, and impute the missing into paternal and maternal haplotype background. 
##the parental background with a high imptued dosage is assigned reference allele. 
##the output is in bimbam mgt format. each SNP occupy three rows: genotype, paternal, maternal. 


    kin = "~/data/rnaseq/geno/kinship.1477.rkm.gz"
    gt=f2+".snpab-bgt.txt"
    cov="~/data/rnaseq/pheno/age.sex.bmi.blood.cells.1477"

    phchr == chr

#    for phchr in list(range(1, 23))[::-1]:
#       if phchr == chr:
#            continue
##uncomment these three lines to do trans test for phenotype on other chromosomes. 

       ph="~/data/rnaseq/pheno/rnaseq.1477.16824.chr"+str(phchr)+".qqnorm"
       out = f2+".16824.chr"+str(phchr)+".qqnorm"
       cmd = "~/work/idul/idul -b -g " + gt + " -c " + cov + " -p " + ph  + " -k " + kin + " -f 0.001 -t 100 -o " +  out
       print(cmd)
#       os.system(cmd)
##this does single marker test for genotype, paternal and maternal alleles for each SNP. 
##-b output Bayes factors; without -b output p-values. 

#       cmd = "~/work/idul/idul -g " + gt + " -c " + cov + " -p " + ph  + " -k " + kin + " -j 2 -b -a -f 0.001 -t 100 -o " +  out 
#       print(cmd) 
#       os.system(cmd)
##this does joint analysis, that is, ignoring genotypes, but do y=Wa+pb1+mb0+u+e; 
##-b output Bayes factors (recommended for 2 d.f. tests). 


##residuals
#       out = "residual.chr"+str(phchr) 
#       cmd = "~/work/idul/idul -r -c " + cov + " -p " + ph  + " -k " + kin + " -o " +  out 
#       print(cmd) 
#       os.system(cmd)

#        cmd = "~/work/idul/idul -g " + gt + " -c " + cov + " -p " + ph  + " -k " + kin + " -0 -j 1 -b -f 0.001 -t 100 -o " +  out 
#        print(cmd) 
#        os.system(cmd)

#        cmd = "~/work/idul/idul -g " + gt + " -c " + cov + " -p " + ph  + " -k " + kin + " -1 -j 1 -b -f 0.001 -t 100 -o " +  out 
#        print(cmd) 
#        os.system(cmd)
