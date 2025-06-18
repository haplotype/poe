/////////////////////////////////////////////////////////////////////////////////
//The MIT License (MIT)
//
//Copyright (c) 2023 Yongtao Guan
//  Bug report: ytguan@gmail.com 
//
//Permission is hereby granted, free of charge, to any person obtaining a copy
//of this software and associated documentation files (the "Software"), to deal
//in the Software without restriction, including without limitation the rights
//to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
//copies of the Software, and to permit persons to whom the Software is
//furnished to do so, subject to the following conditions:
//
//The above copyright notice and this permission notice shall be included in all
//copies or substantial portions of the Software.
//
//THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
//IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
//FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
//AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
//LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
//OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
//SOFTWARE.

#include <random> 
#include <iostream> 
#include <fstream> 
#include <sstream> 
#include <algorithm> 
#include <vector> 
#include <stdio.h> 
#include <stdlib.h> 
#include <cstring> 
#include <string.h> 
#include <unistd.h>
//#include <sys/stat.h>
//#include <sys/types.h>
#include "htslib/hts.h"
#include "htslib/sam.h"
#include "htslib/vcf.h"              
#include <map> 
#include <cmath>
#include <math.h>
#include <pthread.h>
#include <stdint.h>
#include <zlib.h>
//#include "thpool.h"

#define EIGEN_USE_OPENMP
//#define EIGEN_USE_LAPACKE

using namespace std; 

#define BUFFER_SIZE  4096 
#define VERSION "0.31"
//ver0.30 initiated on 17 Aug 2023
//ver0.31 on 19 Aug 2023


void print_progress_num(int last, int p)
{
    if(!last) {
	printf("##processed variants = %d \r", p); 
	fflush(stdout); 	
    } else {
	printf("##processed variants = %d \n", p); 
    }
}

void print_progress_bar(int last, long p, long total)
{
    char str[] = "##processed pairs:"; 
    int progress = (int) (100.0 * p / total); 
//    int barsize = (int) (progress / 2.0); 
//    char bar[100];
//    memset(bar, '\0', 100); 
    if(!last) {
//	    for (int i = 0; i < barsize; i++)
//		    bar[i] = '>'; 
	    printf("%s %ld or %d%%\r", str, p, progress); 
	    fflush(stdout); 	
    } else {
//	    for (int i = 0; i < barsize; i++)
//		    bar[i] = '>'; 
//	    printf("%s [%-50s] 100%%\n", str, bar); 
	    printf("%s %ld or %d%%\n", str, p, progress); 
    }
}


void read_from_gzipped_file(string fn, vector<string> &val, int &rows) {
    gzFile file = gzopen(fn.c_str(), "rb");
    if (file == NULL) {
        std::cerr << "Failed to open file: " << fn << std::endl;
        exit(1);
    }

    rows = 0;
    char buffer[BUFFER_SIZE];
    std::string decompressed_data;

    // Read and decompress the entire file into a string
    int bytes_read;
    while ((bytes_read = gzread(file, buffer, BUFFER_SIZE - 1)) > 0) {
        buffer[bytes_read] = '\0';
        decompressed_data.append(buffer);
    }
    gzclose(file);

    // Process the decompressed data line by line
    std::stringstream ss(decompressed_data);
    std::string line;
    while (std::getline(ss, line)) {
        // Ignore lines starting with #
        if (!line.empty() && line[0] == '#') {
            continue;
        }

	std::stringstream line_stream(line);
	string value;
	while (line_stream >> value) {
	    val.push_back(value); 
	}
        rows++;
    }
}

void read_from_gzipped_num_file(string fn, vector<float> &val, int &rows) {
    gzFile file = gzopen(fn.c_str(), "rb");
    if (file == NULL) {
        std::cerr << "Failed to open file: " << fn << std::endl;
        exit(1);
    }

    rows = 0;
    char buffer[BUFFER_SIZE];
    std::string decompressed_data;

    // Read and decompress the entire file into a string
    int bytes_read;
    while ((bytes_read = gzread(file, buffer, BUFFER_SIZE - 1)) > 0) {
        buffer[bytes_read] = '\0';
        decompressed_data.append(buffer);
    }
    gzclose(file);

    // Process the decompressed data line by line
    std::stringstream ss(decompressed_data);
    std::string line;
    while (std::getline(ss, line)) {
        // Ignore lines starting with #
        if (!line.empty() && line[0] == '#') {
            continue;
        }

	std::stringstream line_stream(line);
	float value;
	while (line_stream >> value) {
	    val.push_back(value); 
	}
        rows++;
    }
}

//bool startsWithHash(const std::string * str) {
//    for (char ch : str) {
//	if(!isspace(ch)) {
//	    return ch == '#';
//	}
//    }
//    return false; 
//}

bool getGzLine(gzFile &gzf, std::string &line) {
    line.clear(); 
//    char buffer[4096]; 

    while (true) {
	int ch = gzgetc(gzf); 
	if(ch == EOF) {
	    if(line.empty()) 
		return false; 
	    return true; 
	}
	if (ch == '\n')
	    return true; 
	line += static_cast<char>(ch); 
    }
}

void transpose(string fnab, string fn, string fout)
{
    int nrows = 0; 
    vector<string> vsnpab; 
    read_from_gzipped_file(fnab, vsnpab, nrows); 
    if(nrows == 0) {
	cout << "zero lines from " << fnab << endl;
	exit(0); 
    }
    cout << "number of lines = " << nrows << endl;

    int ncols = 0; 
    vector<int> val; 
    gzFile gzf = gzopen(fn.c_str(), "rb"); 
    if(!gzf) {
	std::cerr << "Failed to open file." << std::endl;
	return;
    }
    std::string line; 
    getGzLine(gzf, line); 
    int temp;
    while (!line.empty()) 
    {
	ncols++; 
	std::istringstream iss(line); 
	while( iss >> temp) { 
	   val.push_back(temp); 
	}
        getGzLine(gzf, line); 
    }
    gzclose(gzf); 
    cout << "number of cols = " << ncols << endl; 

    cout << val.size() << endl; 
    cout << (long) (nrows) * (long) ncols << endl; 

    FILE * fp = fopen(fout.c_str(), "w");
    if (fp == NULL) {
	    fprintf(stderr, "can't open file %s\n", fout.c_str()); 
	    exit(EXIT_FAILURE);
    }   
    //each SNP output three rows, genotypes, paternal haplotype, maternal haplotype, where genotype for heterozygous are imputed.  
    for (int m = 0; m < nrows; m++)
    {
	fprintf(fp, "%s:2 %s %s ", vsnpab.at(3*m).c_str(), vsnpab.at(3*m+1).c_str(), vsnpab.at(3*m+2).c_str()); 
//	long pos = m; 
	for (int i = 0; i < ncols; i+=2)   //ncols is number of haplotypes;
	{
	    long pos = i*nrows+m; 
	    fprintf(fp, "%d ", val.at(pos)+val.at(pos+nrows));  
//	    pos += nrows*2; 
	}                                       
	fprintf(fp, "\n"); 

	fprintf(fp, "%s:1 %s %s ", vsnpab.at(3*m).c_str(), vsnpab.at(3*m+1).c_str(), vsnpab.at(3*m+2).c_str()); 
//	pos = m; 
	for (int i = 0; i < ncols; i+=2)
	{
	    long pos = i*nrows+m; 
	    fprintf(fp, "%d ", val.at(pos));  
//	    pos += nrows*2; 
	}                                       
	fprintf(fp, "\n"); 

	fprintf(fp, "%s:0 %s %s ", vsnpab.at(3*m).c_str(), vsnpab.at(3*m+1).c_str(), vsnpab.at(3*m+2).c_str()); 
//	pos = m+nrows; 
	for (int i = 0; i < ncols; i+=2)
	{
	    long pos = i*nrows+m; 
	    fprintf(fp, "%d ", val.at(pos+nrows));  
//	    pos += nrows*2; 
	}                                       
	fprintf(fp, "\n"); 
    }
}

double get_nubeam(int *gt, int n)
{
    double ww[4] = {1, sqrt(2), sqrt(3), sqrt(5)}; 
    double prod[4] = {1,0,0,1}; 
    for(int i = 0; i < n; i++)
    {
	if(gt[i] == 0) {
	    prod[2] += prod[0]; 
	    prod[3] += prod[1]; 
	} else {
	    prod[1] += prod[0];  
	    prod[3] += prod[2]; 
	}
    }
//    for (int j = 0; j < 4; j++)
//	cout << prod[j] << " "; 
//    cout << endl;
    double r = prod[0]+prod[3]; 
//    for (int j = 0; j < 4; j++)
//	r += ww[j] * prod[j]; 
    return(log(r)); 
}

void nubeam(string fnab, string fn, string fout, int len)
{
    int nsnp = 0; 
    vector<string> vsnpab; 
    read_from_gzipped_file(fnab, vsnpab, nsnp); 
    if(nsnp == 0) {
	cout << "zero lines from " << fnab << endl;
	exit(0); 
    }
    cout << "number of lines = " << nsnp << endl;

    int nhap = 0; 
    vector<int> val; 
    gzFile gzf = gzopen(fn.c_str(), "rb"); 
    if(!gzf) {
	std::cerr << "Failed to open file." << std::endl;
	return;
    }
    std::string line; 
    getGzLine(gzf, line); 
    int temp;
    while (!line.empty()) 
    {
	nhap++; 
	std::istringstream iss(line); 
	while( iss >> temp) { 
	   val.push_back(temp); 
	}
        getGzLine(gzf, line); 
    }
    gzclose(gzf); 
    cout << "number of cols = " << nhap << endl; 

    cout << val.size() << endl; 
    cout << (long) (nsnp) * (long) nhap << endl; 

    FILE * fp = fopen(fout.c_str(), "w");
    if (fp == NULL) {
	    fprintf(stderr, "can't open file %s\n", fout.c_str()); 
	    exit(EXIT_FAILURE);
    }   
    //each SNP output three rows, genotypes, paternal haplotype, maternal haplotype, where genotype for heterozygous are imputed.  
    int * hap1 = new int[len*2+1]; 
    int * hap0 = new int[len*2+1]; 
    double * nbm1 = new double[nhap]; 
    double * nbm0 = new double[nhap]; 
    for (int m = 0; m < nsnp; m++)
    {
	long pos = 0;    
	for (int i = 0; i < nhap; i+=2) {
	    for(int j = -len; j <= len; j++)
	    {
	        long offset=0; 
		if(m+j <= 0) 
		    offset = 0; 
		else if(m+j >= nsnp-1) 
		    offset = nsnp-1; 
		else 
		    offset = m+j; 
		hap1[len+j] = val.at(pos+offset); 
		//pos+j = (2*i-1)*nsnp + m + j; 
		hap0[len+j] = val.at(pos+nsnp+offset); 
	    }
	    pos += 2*nsnp; 
//	    for(int j = -len; j < len; j++)
//	    {
//		if(m+j <= 0) 
//		    offset = 0; 
//		else if(m+j >= nsnp-1) 
//		    offset = nsnp-1; 
//		else 
//		    offset = pos+m+j; 
//		hap0[j] = val.at(pos1); 
//		//pos =j = 2*i * snp + m +j; 
//	    }                 
	    nbm1[i] = get_nubeam(hap1, 2*len+1); 
	    nbm0[i] = get_nubeam(hap0, 2*len+1); 
	}
//	    if(m == 3) {
//		for(int t = 0; t < len; t++)
//		    cout << hap1[t]; 
//		cout << "\t " << nbm1[i] << endl; 
//		for(int t = 0; t < len; t++)
//		    cout << hap0[t]; 
//		cout << "\t " << nbm0[i] << endl; 
//	    }

//	fprintf(fp, "%s:N %s %s ", vsnpab.at(3*m).c_str(), vsnpab.at(3*m+1).c_str(), vsnpab.at(3*m+2).c_str()); 
	for (int i = 0; i < nhap; i+=2)
	{
	    fprintf(fp, "%.5f ", nbm1[i]+nbm0[i]);  
	}                                       
	fprintf(fp, "\n"); 

//	fprintf(fp, "%s:1 %s %s ", vsnpab.at(3*m).c_str(), vsnpab.at(3*m+1).c_str(), vsnpab.at(3*m+2).c_str()); 
//	for (int i = 0; i < nhap; i+=2)
//	{
//	    fprintf(fp, "%.5f ", nbm1[i]);  
//	}                                       
//	fprintf(fp, "\n"); 
//
//	fprintf(fp, "%s:0 %s %s ", vsnpab.at(3*m).c_str(), vsnpab.at(3*m+1).c_str(), vsnpab.at(3*m+2).c_str()); 
//	for (int i = 0; i < nhap; i+=2)
//	{
//	    fprintf(fp, "%.5f ", nbm0[i]);  
//	}                                       
//	fprintf(fp, "\n"); 
    }
    delete[] nbm1; 
    delete[] nbm0; 
    delete[] hap1; 
    delete[] hap0; 
}

void transmission(string fn, string fout)
{
    FILE * fp = fopen(fout.c_str(), "w");
    if (fp == NULL) {
	    fprintf(stderr, "can't open file %s\n", fout.c_str()); 
	    exit(EXIT_FAILURE);
    }   

    int nlines = 0; 
    vector<int> val; 
    gzFile gzf = gzopen(fn.c_str(), "rb"); 
    if(!gzf) {
	std::cerr << "Failed to open file." << std::endl;
	return;
    }
    std::string line; 
    getGzLine(gzf, line); 
    int temp;
    while (!line.empty()) 
    {
	nlines++; 
	string rs, a, b; 
	std::istringstream iss(line); 
	iss >> rs; 
	iss >> a; 
	iss >> b; 
	while( iss >> temp) { 
	   val.push_back(temp); 
	}

	if(nlines % 3 == 0) {
	    int ni = val.size() / 3; 
	    int n2 = 0; 
	    int n1 = 0; 
	    int n0 = 0; 
	    int pat = 0; 
	    int mat = 0; 
	    for(int i = 0; i < ni; i++)
	    {
               int gt = val.at(i); 
	       int h1 = val.at(ni+i); 
	       int h0 = val.at(ni+ni+i); 
	       if(gt == 2 ) n2 ++; 
	       else if(gt == 0) n0++; 
	       else if(gt == 1) {
                   n1++; 
		   if(h1 == 1) pat++; 
		   if(h0 == 1) mat++; 
	       }
	    }
	    fprintf(fp, "%s %d %d %d %d %d\n", rs.c_str(), n2, n0, n1, pat, mat); 
	    vector<int>().swap(val); 
	}
        getGzLine(gzf, line); 
    }
    gzclose(gzf); 
    fclose(fp); 
}


void simulate_trios(string fn_vcf, string pref) 
{

    htsFile *fpv = hts_open(fn_vcf.c_str(), "r");   
    bcf_hdr_t *hdr = bcf_hdr_read(fpv);  
    bcf1_t* line = bcf_init();   
    int ni = bcf_hdr_nsamples(hdr); 
    int ntrios = ni / 2; 
    fprintf(stdout, "##numbers of samples and trios: %d  %d\n", ni, ntrios); 

    FILE * fp2 = NULL; 
    {
	string buf2; 
	buf2.assign(pref); 
	buf2.append(".bimbam.pos");  
	fp2 = fopen(buf2.c_str(), "w");
	if (fp2 == NULL) {
		fprintf(stderr, "can't open file %s\n", buf2.c_str()); 
		exit(EXIT_FAILURE);
	}   // for SNPs pvals etc. 
    }

    FILE * fp3 = NULL; 
    {
	string buf3; 
	buf3.assign(pref); 
	buf3.append(".bimbam.gt");  
	fp3 = fopen(buf3.c_str(), "w");
	if (fp3 == NULL) {
		fprintf(stderr, "can't open file %s\n", buf3.c_str()); 
		exit(EXIT_FAILURE);
	}   // for SNPs pvals etc. 
	fprintf(fp3, "%d = \n 1000 \n", ntrios); 
    }

    FILE * fp6 = NULL; 
    {
	string buf6; 
	buf6.assign(pref); 
	buf6.append(".truth.txt");  
	fp6 = fopen(buf6.c_str(), "w");
	if (fp6 == NULL) {
		fprintf(stderr, "can't open file %s\n", buf6.c_str()); 
		exit(EXIT_FAILURE);
	}    
	fprintf(fp6, "%d = \n 1000 \n", ntrios); 
    }

    int ns = 0; 
    int ns1 = 1; 

    float * snpgt = new float[ntrios]; 
    float * hap1 = new float[ntrios]; 
    float * hap0 = new float[ntrios]; 

    int nlines = 0; 
    while(bcf_read1(fpv, hdr, line) == 0) {   

	if(ns1 % 1000 == 0) 
	    print_progress_num(0, ns1);  
	ns1++; 

	bcf_unpack(line,BCF_UN_STR); 
	bcf_get_variant_types(line); 
	if(line->d.var_type != VCF_SNP) continue; 
	if(line->n_allele > 2) continue;       

	int32_t *gt_arr = NULL, ngt_arr = 0;
	int ngt;
	ngt = bcf_get_genotypes(hdr, line, &gt_arr, &ngt_arr);
//	    cout << "ngt = " << ngt << endl; 
	assert(ngt == ni*2); 
	if ( ngt<=0 ) 
	{
	    cout << "no genotypes at " << bcf_seqname(hdr, line) << " " << line->pos << endl; 
	    continue;
	    //fiter by gt presence. 
	}

	string allele_a(line->d.allele[0]); 
	string allele_b(line->d.allele[1]); 
	if(allele_a.size() > 1 || allele_b.size() > 1) {
	    allele_a.assign("+"); 
	    allele_b.assign("-"); 
	}

	string rs; 
	rs = std::to_string(line->rid+1)+ "_" + std::to_string(line->pos+1); 
	
	fprintf(fp2, "%s:%s:%s %ld %d\n", rs.c_str(), allele_a.c_str(), allele_b.c_str(), line->pos+1, line->rid+1); 
	fprintf(fp3, "%s:%s:%s ", rs.c_str(),allele_a.c_str(), allele_b.c_str()); 
	fprintf(fp6, "%s:%s:%s ", rs.c_str(),allele_a.c_str(), allele_b.c_str()); 

	for (int j=0; j<ntrios; j++)
	{
	    int i = j*2;  
	    int gt_fa = bcf_gt_allele(gt_arr[2*i+0]) + bcf_gt_allele(gt_arr[2*i+1]);
	    int gt_mo = bcf_gt_allele(gt_arr[2*i+2]) + bcf_gt_allele(gt_arr[2*i+3]);
	    int gt_ch = bcf_gt_allele(gt_arr[2*i+0]) + bcf_gt_allele(gt_arr[2*i+2]);
	    hap1[j] = bcf_gt_allele(gt_arr[2*i+0]); 
	    hap0[j] = bcf_gt_allele(gt_arr[2*i+2]); 
	    if(gt_ch == 0 || gt_ch == 2) 
		snpgt[j] = gt_ch; 
	    else if(gt_ch == 1) 
	    {
		//0 0; error; 
		//0 1; good;      0 1
		//0 2; good;      0 1
		//1 0; good;      1 0
		//1 1; mask; 
		//1 2; good;      0 1
		//2 0; good;      1 0
		//2 1; good;      1 0 
		//2 2; error; 
		if(gt_fa == 0  && gt_mo == 1)
		    snpgt[j] = 7; 
		else if(gt_fa == 0  && gt_mo == 2)
		    snpgt[j] = 7; 
		else if(gt_fa == 1  && gt_mo == 0)
		    snpgt[j] = 8; 
		else if(gt_fa == 1  && gt_mo == 2)
		    snpgt[j] = 7; 
		else if(gt_fa == 2  && gt_mo == 0)
		    snpgt[j] = 8; 
		else if(gt_fa == 2  && gt_mo == 1)
		    snpgt[j] = 8; 
		else if(gt_fa == 1  && gt_mo == 1)
		{   
		    snpgt[j] = 9; 
		}
		else {
		    snpgt[j] = 10; 
		}
	    }
	}
	free(gt_arr);

	ns++; 
	char A = allele_a[0]; 
	char B = allele_b[0]; 
	for (int i = 0; i < ntrios; i++) 
	{
	    if(int(hap1[i]) == 0 && int(hap0[i]) == 0) 
		fprintf(fp6, "%c%c ", A, A);  
	    if(int(hap1[i]) == 0 && int(hap0[i]) == 1) 
		fprintf(fp6, "%c%c ", A, B);  
	    if(int(hap1[i]) == 1 && int(hap0[i]) == 0) 
		fprintf(fp6, "%c%c ", B, A);  
	    if(int(hap1[i]) == 1 && int(hap0[i]) == 1) 
		fprintf(fp6, "%c%c ", B, B);  
	}
	fprintf(fp6, "\n");  
	

	for (int i = 0; i < ntrios; i++) 
	{
	    if(snpgt[i] == 0) 
		fprintf(fp3, " %c%c", A, A);  
	    else if(snpgt[i] == 2) 
		fprintf(fp3, " %c%c", B, B);  
	    else if(snpgt[i] == 7) 
		fprintf(fp3, " %c%c", A, B);  
	    else if(snpgt[i] == 8) 
		fprintf(fp3, " %c%c", B, A);  
	    else if(snpgt[i] == 9) 
		fprintf(fp3, " NN");  
	    else if(snpgt[i] == 10) 
		fprintf(fp3, " ??");  
	}
	fprintf(fp3, "\n");  

	nlines++; 
	//fill the missing gentoypes by the mean of the rest. 
    }
    delete[] snpgt; 
    delete[] hap1; 
    delete[] hap0; 


    fprintf(stdout, "##number of biallelic SNPs: %d \n", ns); 
    print_progress_num(1, ns1);  
    hts_close(fpv); 
    if(ns == 0) 
    {
       fprintf(stdout, "##no valid SNPs. abort. \n"); 
       exit(0); 
    }
}


int usage()
{
	fprintf(stderr, "\n");
	fprintf(stderr, "Version: %s\n", VERSION);
	fprintf(stderr, "Usage:   poe -i in.vcf.gz -f trio.txt -o pref\n");
	fprintf(stderr, "Usage:       this will produce partially phased child genotypes.\n");
	fprintf(stderr, "Usage:       output pref.bimbam.pos and pref.bimbam.gt\n");
	fprintf(stderr, "Usage:   poe -s snpab.txt -g mgt.txt -o pref\n");
	fprintf(stderr, "Usage:       this will transpose genotypes and add SNP ID for association.\n");
	fprintf(stderr, "Usage:       output pref.snpab-mgt.txt\n");
	fprintf(stderr, "Usage:   poe -s snpab.txt -g mgt.txt -o pref -n len\n");
	fprintf(stderr, "Usage:       this will produce nubeam input file.\n");
	fprintf(stderr, "Usage:       output pref.nubeam.len.txt\n");
	fprintf(stderr, "Usage:   poe -g bgt.txt -o pref -t\n");
	fprintf(stderr, "Usage:       this will output SNPID n2 n0 n1 n_pat n_mat \n");
	fprintf(stderr, "Usage:       output pref.trans.txt\n");
        fprintf(stderr, "Options: \n");
	fprintf(stderr, "         -f str        trio file name; three id in a row: fa mo child\n");    
	fprintf(stderr, "         -h            print usage\n");
	fprintf(stderr, "         -i str        (indexed) bgzipped vcf file\n");
	fprintf(stderr, "         -n len        nubeam half haplotype length\n");
	fprintf(stderr, "         -o str        output prefix [out]\n");
	fprintf(stderr, "         -g str        genotype file (each row is an individual)\n");
	fprintf(stderr, "         -s str        snp file (each row is a snp)\n");
	fprintf(stderr, "         -t            tally transmission counts per SNP\n");
	fprintf(stderr, "\n");
	fprintf(stderr, "Bug report: Yongtao Guan <ytguan@gmail.com>\n\n");
	return 1;
}

int main(int argc, char *argv[])
{
//    time_t time0, time1, time2; 
    string fn_vcf; 
    string fn_fam; 
    string pref("out");
    char c;
    int flag_vcf = 0; 
    int flag_fam = 0; 
    int flag_help = 0; 

    string fn_snpab; 
    string fn_gt; 
    int flag_gt = 0; 
    int flag_snpab = 0; 

    int flag_nubeam = 0; 
    int nbm_len = 0; 

    int flag_trans_pval = 0; 
    int flag_x = 0; //flag of simulating cross-over; 

    while ((c = getopt(argc, argv, "f:g:hi:n:o:s:tx")) >= 0) 
    {
	switch (c) {
	    case 'f': 
		flag_fam = 1; 
		fn_fam.assign(optarg); 
		break;
	    case 'g': 
		flag_gt = 1; 
		fn_gt.assign(optarg); 
		break;
	    case 'h': 
		flag_help = 1; 
		break; 
	    case 'i': 
		flag_vcf = 1; 
		fn_vcf.assign(optarg); 
		break;
	    case 'n': 
		flag_nubeam = 1; 
		nbm_len = atoi(optarg); 
		break;
	    case 'o': 
		pref.assign(optarg); 
		break;
	    case 't': 
		flag_trans_pval = 1; 
		break; 
	    case 's': 
		flag_snpab = 1; 
		fn_snpab.assign(optarg); 
		break;
	    case 'x': 
		flag_x = 1; 
		break; 
	    default: 
		break; 
	}
    }


    fprintf(stdout, "\n"); 
    fprintf(stdout, "POE %s by Yongtao Guan at Framingham Heart Study \n", VERSION); 
    fprintf(stdout, "  National Heart, Lung, and Blood Institute (C) 2023 \n"); 

    FILE * fplog = NULL;
    string buf; 
    buf.assign(pref); 
    buf.append(".log");  
    fplog = fopen(buf.c_str(), "w");
    if (fplog == NULL) {
	    fprintf(stderr, "can't open file %s\n", buf.c_str()); 
	    exit(EXIT_FAILURE);
    }   // for SNPs. 
    for (int i = 0; i < argc; i++)
	fprintf(fplog, "%s ", argv[i]); 
    fprintf(fplog, "\n"); 

    if(flag_trans_pval == 1) 
    {
	string buf; 
	buf.assign(pref); 
	buf.append(".trans.txt");  
	cout << buf << endl; 
	transmission(fn_gt, buf); 
	exit(EXIT_SUCCESS); 

    }
    if(flag_x == 1) 
    {
	simulate_trios(fn_vcf, pref); 
	exit(EXIT_SUCCESS); 
    }
    
    if(flag_gt == 1 && flag_snpab == 1) 
    {
	if(flag_nubeam == 0)  {
	    string buf; 
	    buf.assign(pref); 
	    buf.append(".snpab-bgt.txt");  
	    cout << buf << endl; 
	    transpose(fn_snpab, fn_gt, buf); 
	}
	else {
	    string buf(pref);        
	    buf = buf + ".nubeam."; 
            buf = buf + std::to_string(nbm_len); 
	    buf = buf + ".txt"; 
	    cout << buf << endl; 
	    {
		int hap[7]={1,1,1,1,1,1,1};
		cout << get_nubeam(hap,7) << endl; 
	    }
	    {
		int hap[7]={1,0,1,0,1,0,1};
		cout << get_nubeam(hap,7) << endl; 
	    }
	    nubeam(fn_snpab, fn_gt, buf, nbm_len); 
	}
	exit(EXIT_SUCCESS); 
    }

    if (flag_help == 1 || flag_vcf == 0 || flag_fam == 0 ) return usage();

    int ntrios = 0; 
    std::vector<string> vec_fam; 
    /////////////////////////////////////////////////////
    //read fam;  
    {
	read_from_gzipped_file(fn_fam, vec_fam, ntrios); 
	//a row begin with # will be ignored. 
	int ncols = vec_fam.size() / ntrios; 
	cout << "number of trios = " << ntrios << endl;  
	assert(ncols == 3); 
    }

    FILE * fp2 = NULL; 
    {
	string buf2; 
	buf2.assign(pref); 
	buf2.append(".bimbam.pos");  
	fp2 = fopen(buf2.c_str(), "w");
	if (fp2 == NULL) {
		fprintf(stderr, "can't open file %s\n", buf2.c_str()); 
		exit(EXIT_FAILURE);
	}   // for SNPs pvals etc. 
    }

    FILE * fp3 = NULL; 
    {
	string buf3; 
	buf3.assign(pref); 
	buf3.append(".bimbam.gt");  
	fp3 = fopen(buf3.c_str(), "w");
	if (fp3 == NULL) {
		fprintf(stderr, "can't open file %s\n", buf3.c_str()); 
		exit(EXIT_FAILURE);
	}   // for SNPs pvals etc. 
	fprintf(fp3, "%d = \n 1000 \n", ntrios); 
	fprintf(fp3, "IND "); 
	for (int i = 0; i < ntrios; i++)
	    fprintf(fp3, "%s ", vec_fam.at(i*3+2).c_str()); 
	fprintf(fp3, "\n");

    }

    FILE * fp4 = NULL; 
    {
	string buf4; 
	buf4.assign(pref); 
//	buf4.append(".untransmitted.gt.012");  
	buf4.append(".snpab-untransmitted.gt");  
	fp4 = fopen(buf4.c_str(), "w");
	if (fp4 == NULL) {
		fprintf(stderr, "can't open file %s\n", buf4.c_str()); 
		exit(EXIT_FAILURE);
	}    
    }

    FILE * fp5 = NULL; 
    {
	string buf5; 
	buf5.assign(pref); 
	buf5.append(".af.txt");  
	fp5 = fopen(buf5.c_str(), "w");
	if (fp5 == NULL) {
		fprintf(stderr, "can't open file %s\n", buf5.c_str()); 
		exit(EXIT_FAILURE);
	}    
    }

    FILE * fp6 = NULL; 
    {
	string buf6; 
	buf6.assign(pref); 
	buf6.append(".bimbam.gt.012");  
	fp6 = fopen(buf6.c_str(), "w");
	if (fp6 == NULL) {
		fprintf(stderr, "can't open file %s\n", buf6.c_str()); 
		exit(EXIT_FAILURE);
	}    
    }


    std::map<string, int> rs2yes; //if this SNP ID exits
    //////////////////////////////////////////////////////
    std::map<string, int> id2idx; //sample id to genotype index; 
    if(flag_vcf == 1) 
    {
	htsFile *fpv = hts_open(fn_vcf.c_str(), "r");   
	bcf_hdr_t *hdr = bcf_hdr_read(fpv);  
	bcf1_t* line = bcf_init();   
	int ni = bcf_hdr_nsamples(hdr); 
	fprintf(fplog, "##number of samples: %d \n", ni); 
	fprintf(stdout, "##number of samples: %d \n", ni); 
    	if(hdr->samples != NULL) {
    	    for (int i = 0; i < ni; i++)       
	    {
		string sa(hdr->samples[i]); 
		id2idx.insert(std::make_pair(sa, i)); 
	    }
    	} else {
	    cout << "no sample ID in vcf, quit. " << endl; 
	    exit(0); 
	}

	int ns = 0; 
	int ns1 = 1; 

	float * hap1 = new float[ntrios]; 
	float * hap0 = new float[ntrios]; 
	float * snpgt = new float[ntrios]; 
	float * snpgt2 = new float[ntrios]; 
	float * snpgt3 = new float[ntrios]; //child gt.012
	float * menderr = new float[ntrios]; 
	for(int i = 0; i < ntrios; i++)
	    menderr[i] = 0; 
	int nlines = 0; 
	while(bcf_read1(fpv, hdr, line) == 0) {   

	    if(ns1 % 1000 == 0) 
	        print_progress_num(0, ns1);  
	    ns1++; 

	    bcf_unpack(line,BCF_UN_STR); 
	    bcf_get_variant_types(line); 
	    if(line->d.var_type != VCF_SNP) continue; 
	    if(line->n_allele > 2) continue;       

	    int32_t *gt_arr = NULL, ngt_arr = 0;
            int ngt;
	    ngt = bcf_get_genotypes(hdr, line, &gt_arr, &ngt_arr);
//	    cout << "ngt = " << ngt << endl; 
	    assert(ngt == ni*2); 
	    if ( ngt<=0 ) 
	    {
		cout << "no genotypes at " << bcf_seqname(hdr, line) << " " << line->pos << endl; 
		continue;
		//fiter by gt presence. 
	    }

	    string allele_a(line->d.allele[0]); 
	    string allele_b(line->d.allele[1]); 
	    if(allele_a.size() > 1 || allele_b.size() > 1) {
		allele_a.assign("+"); 
		allele_b.assign("-"); 
	    }

	    string rs; 
	    if(line->d.id[0] == '.') {
		rs = std::to_string(line->rid+1)+ "_" + std::to_string(line->pos+1); 
	    } else {
		rs.assign(line->d.id); 
	    }
	    map<string, int>::iterator it; 
	    it = rs2yes.find(rs); 
	    while (it != rs2yes.end())   
	    {
		rs.append("_dup");
	        it = rs2yes.find(rs); 
	    } 
	    rs2yes[rs] = 1; 
	    fprintf(fp2, "%s:%s:%s %ld %d\n", rs.c_str(), allele_a.c_str(), allele_b.c_str(), line->pos+1, line->rid+1); 
	    fprintf(fp3, "%s:%s:%s ", rs.c_str(),allele_a.c_str(), allele_b.c_str()); 
	    fprintf(fp4, "%s:%s:%s:%d:%ld %s %s ", rs.c_str(),allele_a.c_str(), allele_b.c_str(), line->rid+1, line->pos+1, allele_a.c_str(), allele_b.c_str()); 

	    double sumf = 0; 
	    double summ = 0; 
	    double sumc = 0; 
	    int merr = 0; 
	    int trih = 0; 
	    for (int i=0; i<ntrios; i++)
	    {
		int fa_idx=-1, mo_idx=-1, ch_idx=-1; 

		map<string, int>::iterator it; 
		it = id2idx.find(vec_fam.at(i*3+0)); 
		if(it != id2idx.end()) 
		    fa_idx = it->second; 
		it = id2idx.find(vec_fam.at(i*3+1)); 
		if(it != id2idx.end()) 
		    mo_idx = it->second; 
		it = id2idx.find(vec_fam.at(i*3+2)); 
		if(it != id2idx.end()) 
		    ch_idx = it->second; 
		assert(fa_idx>=0); 
		assert(mo_idx>=0); 
		assert(ch_idx>=0); 

		int gt_fa = bcf_gt_allele(gt_arr[2*fa_idx+0]) + bcf_gt_allele(gt_arr[2*fa_idx+1]);
		int gt_mo = bcf_gt_allele(gt_arr[2*mo_idx+0]) + bcf_gt_allele(gt_arr[2*mo_idx+1]);
		int gt_ch = bcf_gt_allele(gt_arr[2*ch_idx+0]) + bcf_gt_allele(gt_arr[2*ch_idx+1]);
		hap1[i] = bcf_gt_allele(gt_arr[2*ch_idx+0]); 
		hap0[i] = bcf_gt_allele(gt_arr[2*ch_idx+1]); 
		sumf += gt_fa; 
		summ += gt_mo; 
		sumc += gt_ch; 
		snpgt2[i] = gt_fa + gt_mo - gt_ch; 
		snpgt3[i] = gt_ch; 
		if(gt_ch == 0 || gt_ch == 2) 
		    snpgt[i] = gt_ch; 
		else if(gt_ch == 1) 
		{
		    //0 0; error; 
		    //0 1; good;      0 1
		    //0 2; good;      0 1
		    //1 0; good;      1 0
		    //1 1; mask; 
		    //1 2; good;      0 1
		    //2 0; good;      1 0
		    //2 1; good;      1 0 
		    //2 2; error; 
		    if(gt_fa == 0  && gt_mo == 1)
			snpgt[i] = 7; 
		    else if(gt_fa == 0  && gt_mo == 2)
			snpgt[i] = 7; 
		    else if(gt_fa == 1  && gt_mo == 0)
			snpgt[i] = 8; 
		    else if(gt_fa == 1  && gt_mo == 2)
			snpgt[i] = 7; 
		    else if(gt_fa == 2  && gt_mo == 0)
			snpgt[i] = 8; 
		    else if(gt_fa == 2  && gt_mo == 1)
			snpgt[i] = 8; 
		    else if(gt_fa == 1  && gt_mo == 1)
		    {   
			trih ++; 
			snpgt[i] = 9; 
		    }
		    else {
			snpgt[i] = 10; 
			menderr[i]++;
			merr ++; 
		    }
		}
	    }
	    double n2=ntrios*2.0; 
	    sumf /= n2; 
	    summ /= n2; 
	    sumc /= n2; 
	    fprintf(fp5, "%s:%s:%s %d %ld %5.4f %5.4f %5.4f %d %d\n", rs.c_str(),allele_a.c_str(), allele_b.c_str(), line->rid+1, line->pos+1, sumf, summ, sumc, merr, trih); 

	    free(gt_arr);

	    ns++; 
	    char A = allele_a[0]; 
	    char B = allele_b[0]; 

	    for (int i = 0; i < ntrios; i++) 
	    {
		if(snpgt[i] == 0) 
		    fprintf(fp3, " %c%c", A, A);  
		else if(snpgt[i] == 2) 
		    fprintf(fp3, " %c%c", B, B);  
		else if(snpgt[i] == 7) 
		    fprintf(fp3, " %c%c", A, B);  
		else if(snpgt[i] == 8) 
		    fprintf(fp3, " %c%c", B, A);  
		else if(snpgt[i] == 9) 
		    fprintf(fp3, " NN");  
		else if(snpgt[i] == 10) 
		    fprintf(fp3, " ??");  
	    }
	    fprintf(fp3, "\n");  

	    for (int i = 0; i < ntrios; i++) 
	    {
		int temp = (int) snpgt2[i]; 
		temp = temp > 2 ? 2 : temp; 
		temp = temp < 0 ? 0 : temp; 
		fprintf(fp4, " %d", temp);  
	    }
	    fprintf(fp4, "\n");  

	    for (int i = 0; i < ntrios; i++) 
	    {
		int temp = (int) snpgt3[i]; 
		temp = temp > 2 ? 2 : temp; 
		temp = temp < 0 ? 0 : temp; 
		fprintf(fp6, " %d", temp);  
	    }
	    fprintf(fp6, "\n");  
	              
	    nlines++; 
	    //fill the missing gentoypes by the mean of the rest. 
        }
	delete[] snpgt; 
	delete[] snpgt2; 
	delete[] snpgt3; 
//	delete[] hap1; 
//	delete[] hap0; 


	fprintf(fplog, "##number of biallelic SNPs: %d \n", ns); 
	fprintf(stdout, "##number of biallelic SNPs: %d \n", ns); 
//	fprintf(fplog, "##number of biallelic SNPs without tag: %s = %d \n", af_tag.c_str(), ns0); 
//	fprintf(stdout, "##number of biallelic SNPs without tag: %s = %d \n", af_tag.c_str(), ns0); 
//	for (int i = 0; i < ntrios; i++)
//	   fprintf(fplog, "%f ", menderr[i]); 
	print_progress_num(1, ns1);  
	fflush(fplog); 
	hts_close(fpv); 
	if(ns == 0) 
	{
	   fprintf(fplog, "##no valid SNPs. abort. \n"); 
	   fprintf(stdout, "##no valid SNPs. abort. \n"); 
	   fflush(fplog); 
	   fclose(fplog); 
	   exit(0); 
	}
    }
    fclose(fp2); 
    fclose(fp3); 
    fclose(fp4); 
    fclose(fp5); 
    fclose(fp6); 
    fclose(fplog); 
///////////////////////////////////////////////////////

    return 0;
}




