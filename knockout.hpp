#include "akt.hpp"

#include "kseq.h"
#include "kstring.h"
#include "pedigree.h"
#include <iomanip>  
#include <vector>
#include <zlib.h>

typedef struct _args
{
    char *pedigree,*inputfile,*include,*gene_bedfile;
} args;


class Gene{ 
public:
    Gene(const string & gene_name,int rid,int start,int stop);
    int isOverlapping(int rid,int a,int b);
    int getStart() {return _start;};
    int getStop() {return _stop;};
    int getChrom() {return _rid;};    
    const string  & getGeneName() {return _gene_name;};
    void print() {cerr<<_rid<<"\t"<<_start<<"\t"<<_stop<<"\t"<<_gene_name<<endl;}
private:
    string _gene_name;
    int _rid,_start,_stop;
};

//list of genes with some basic query tools
class GeneList{ 
public:
    GeneList(char *fname,bcf_hdr_t *h);
    const string & getGeneName(int rid,int a,int b);
    vector< vector<Gene> > _genes; //vectors of Genes (one per rid)
};
