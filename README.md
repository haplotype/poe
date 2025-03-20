# poe
Software and pipeline for parent-of-origin analysis for trios

The following is the "pipeline" that produces POE eQTL results in the manuscript [Abundant Parent-of-origin Effect eQTL in Humans: The Framingham Heart Study](https://www.biorxiv.org/content/10.1101/2024.06.05.597677v1.abstract). The software [poe](https://www.github.com/haplotype/poe) perform data preparation, partial phasing based Mendelian inheritance, and masking triple heterozygous SNPs. Then [elai-mt](https://www.github.com/haplotype/elai) performs imputation of masked heterozygous into paternal and maternal haplotype background to obtain complete phasing of the child. Finally [idul](https://www.github.com/haplotype/idul) performs association analysis by fitting linear mixed models. [Kindred](https://www.github.com/haplotype/kindred) can be used to infer genetic relatedness matrix (GRM) used in idul (other methods to infer GRM can also work).   

```
import os
import subprocess

for chr in list(range(1,23))[::-1]:
    fi = "wgs4rnaseq.chr"+str(chr)+".2snps.vcf.gz"
    fo = "filter-pass.chr"+str(chr)+".vcf"
    cmd = "bcftools view -f 'PASS' " + fi + " > " + fo  
    print(cmd) 
#    os.system(cmd)            
##the PASS is all you need. The vcf file contain genotypes of all samples in multiple trios.  

    f2 = "filter-pass-1477.chr"+str(chr)
    cmd = "~/work/poe/poe  -f fa-mo-son.1477 -i " + fo + " -o " + f2  + " " 
    print(cmd) 
#    os.system(cmd)
##here the program poe is in poe-linux.tar, which contains a suite of tools developed to aid the parent-of-origin analysis using idul. 
##in command -f takes a trio file (ID only) and vcf to extract child genotype in the order specified in the trio file.
##one line in the trio file looks like this: FatherID MotherID ChildID
##the child genotypes are partially phased by Mendelian inheritance, and the triple heterozygous sites are marked as missing. 
##output genotypes of the children are in phased bimbam format that will feed to elai-mt.  

    gt=f2+".bimbam.gt"
    pos=f2+".bimbam.pos" 
    out=f2+'.mt'
    cmd = "~/work/elai/elai-mt -g " + gt + " -p 1 -exclude-maf 0.01 -C 2 -c 10 -s 20 -w 0 -wmg -weg -pos " + pos + " -o " + out + " -e 10 -nthreads 64" 
    print(cmd) 
#    os.system(cmd)
##this runs multi-threading elai, which can be found at github.com/haplotype/elai. 
##this special edition of elai takes phased genotype, fit a LD model,
##and impute the missing genotypes (those triple heterozygous) into paternal and maternal haplotype background. 
##the parental background with a high imptued dosage is assigned reference allele, the other alterantive allele.  
##this step will produce two files *.bgt.txt and *.snpab.txt. 

    fmgt="output/"+out+".bgt.txt"
    fsab="output/"+out+".snpab.txt"
    cmd="~/work/poe/poe -s "+fsab + " -g " + fmgt + " -o "+f2
    print(cmd)
#    os.system(cmd)
##poe transpose bgt.txt and write an output in bimbam mgt format. Each SNP occupies three rows: genotype, paternal, maternal.
##output file name is *.snpab-bgt.txt

    kin = "~/data/rnaseq/geno/kinship.1477.rkm.gz"
    gt=f2+".snpab-bgt.txt"
    cov="~/data/rnaseq/pheno/age.sex.bmi.blood.cells.1477"
    phchr == chr 
    ph="~/data/rnaseq/pheno/rnaseq.1477.16824.chr"+str(phchr)+".qqnorm"
    out = f2+".16824.chr"+str(phchr)+".qqnorm"
    cmd = "~/work/idul/idul -b -g " + gt + " -c " + cov + " -p " + ph  + " -k " + kin + " -f 0.001 -t 100 -o " +  out
    print(cmd)
#    os.system(cmd)
##this command use idul to fit a linear mixed model for association test.
##IDUL can be found at github.com/haplotype/idul.
##the kinship matrix can be inferred using Kindred, which can be found at github.com/haplotype/kindred
##this does single marker test for genotype, paternal and maternal alleles for each SNP. 
##-b output Bayes factors; without -b output p-values. 

#   cmd = "~/work/idul/idul -g " + gt + " -c " + cov + " -p " + ph  + " -k " + kin + " -j 2 -b -a -f 0.001 -t 100 -o " +  out 
#   print(cmd) 
#   os.system(cmd)
##this does joint analysis, that is, ignoring genotypes, but do y=Wa+pb1+mb0+u+e; 
##-b output Bayes factors (recommended for 2 d.f. tests).
```
