#include "knockout.hpp"


KSTREAM_INIT(gzFile, gzread, 16384)
using namespace std;


Gene::Gene(const string & gene_name,int rid,int start,int stop)
{
    _gene_name=gene_name;
    _rid=rid;
    _start=start;
    _stop=stop;
};

int Gene::isOverlapping(int rid,int a,int b)
{
    if(rid!=_rid)
    {
	return(0);
    }
    if(a >= _start && a<= _stop)
    {
	return(1);
    }
    if(b >= _start && b<= _stop)
    {
	return(2);
    }
    
    return(0);
};


GeneList::GeneList(char *fname,bcf_hdr_t *hdr)
{
    gzFile fp = gzopen(fname, "r");
    if(fp==NULL) {
	fprintf(stderr,"problem opening %s\n",fname);
	exit(1);
    }
    kstream_t *ks = ks_init(fp);
    kstring_t str = {0,0,0};
    ks_tokaux_t aux;
    int    count=0;
    cerr << "Reading genes...";
    while(    ks_getuntil(ks, '\n', &str, 0) >=0) {
	char *ptr = kstrtok(str.s,"\t",&aux);//chrom
	stringstream ss;
	ss << ptr;
	string chrom,start,stop,gene;
	ss >> chrom;
	ss >> start;
	ss >> stop;
	ss >> gene;
	int rid = bcf_hdr_name2id(hdr, chrom.c_str());
	if(rid == -1)
	{
	    cerr << "WARNING: gene contig was not in the VCF header (ignored): "<< ptr << endl;

	}
	else
	{
	    int a=atoi(start.c_str())-1;
	    int b=atoi(stop.c_str())-1;
	    count++;
	    assert(rid>=0);
	    while(_genes.size()<=rid)
	    {
		_genes.push_back(vector<Gene>());
	    }
	    assert(_genes.size()>rid);
	    Gene g(gene,rid,a,b);
//	    g.print();
	    _genes[rid].push_back(g);
	}
    }
    cerr<<"found "<<count<<endl;
}
	

/**
 * @name    usage
 * @brief   print out options
 *
 * List of input options
 *
 */
static void usage()
{ 
    fprintf(stderr, "\n");
    fprintf(stderr, "About:   akt knockout - profiles duo/trios\n");
    fprintf(stderr, "Usage:   akt knockout input.bcf -p pedigree.fam -G genes.bed.gz\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "Options:\n");
    fprintf(stderr, "    -p, --pedigree              pedigree information in plink .fam format\n");
    fprintf(stderr, "    -i, --include               your definition of LoF goes here eg. 'INFO/CSQ[*]~\"stop_gained\" | INFO/CSQ[*]~\"frameshift_variant\" | INFO/CSQ[*]~\"splice_acceptor_variant\" | INFO/CSQ[*]~\"splice_donor_variant\"\n");
	fprintf(stderr, "    -G, --genes        <file>   four column bed file containing your gene list\n");
    exit(1);
}

int knockout(args & a)
{
    //open a file.
    bcf_srs_t *sr =  bcf_sr_init() ; 

    if(a.gene_bedfile!=NULL)
    {
	if ( bcf_sr_set_targets(sr, a.gene_bedfile, 1 ,0)<0 )
	{
	    cerr << "ERROR: Failed to set targets " <<a.gene_bedfile<<endl;;
	    exit(1);
	}
    }

    if(bcf_sr_add_reader(sr,a.inputfile)!=1) exit(1);
    GeneList genes = GeneList(a.gene_bedfile,sr->readers[0].header);
    bcf_hdr_t *hdr=sr->readers[0].header;
    filter_t *filter=NULL;
    if(a.include!=NULL)   
	filter = filter_init(hdr, a.include);

    //note: this calls set_samples on the hdr
    sampleInfo ped(a.pedigree,hdr);

    int nsnp=0;

    bcf1_t *line=bcf_init1();;
    int  nsample = ped.N;
    int ngt = 2*nsample;
    int *gt_arr=(int *)malloc(ngt * sizeof(int));

    cerr << "Reading input from "<<a.inputfile<<endl;
    cout << "CHROM\tPOS\tREF\tALT\tGENE\tSAMPLE\tGT" <<endl;
    vector<int> lof_counts(nsample,0);
    vector< vector<string> > lof_gts(nsample, vector<string>() );
    int prev_rid = -1;
    vector<Gene>::iterator gene_it;
    while(bcf_sr_next_line (sr))
    {
	line =  bcf_sr_get_line(sr, 0);
	bcf_unpack(line,BCF_UN_ALL);
	assert(line->n_allele==2);	
	if(a.include==NULL||filter_test(filter,line,NULL))
	{
	    
	    if(prev_rid!=-1 && (line->pos > gene_it->getStop() ||line->rid != prev_rid) )//flush knockouts for last gene.
	    {
		for(int i=0;i<nsample;i++)
		{
		    if(lof_counts[i]>1)
		    {
			assert(lof_counts[i]>=lof_gts[i].size());
			for(size_t j=0;j<lof_gts[i].size();j++)
			{
			    cout << lof_gts[i][j];			    
			}			
//			cout << hdr->samples[i] <<"\t"<<bcf_hdr_int2id(hdr,BCF_DT_CTG,gene_it->getChrom()) <<"\t"<<gene_it->getStart()<<			    "\t"<<gene_it->getStop()<<"\t"<<gene_it->getGeneName()<<"\t"<<lof_counts[i] <<endl;
		    }
		}
		lof_counts.assign(nsample,0);
		lof_gts.assign(nsample, vector<string>() );
	    }
	    if(line->rid != prev_rid)
	    {
		prev_rid = line->rid;	    
		gene_it=genes._genes[prev_rid].begin();
	    }
	    
	    //move the gene iterator forward until we get past the position (then move it back one)
	    while(gene_it!=genes._genes[prev_rid].end() &&  gene_it->getStart()<line->pos)
	    {
		gene_it++;
	    }
	    gene_it--;

	    if(  gene_it->isOverlapping(line->rid,line->pos,line->pos+strlen(line->d.allele[1])) )
	    {
		int ret = bcf_get_genotypes(hdr, line, &gt_arr, &ngt);
		bool diploid = ret==2*nsample;
		for(int i=0;i<nsample;i++)
		{
		    if(!bcf_gt_is_missing(gt_arr[i*2]))
		    {
			lof_counts[i]+=bcf_gt_allele(gt_arr[i*2]);
		    }
		    if(!bcf_gt_is_missing(gt_arr[i*2+1]))
		    {
			lof_counts[i]+=bcf_gt_allele(gt_arr[i*2+1]);
		    }
		    if( !bcf_gt_is_missing(gt_arr[i*2]) && !bcf_gt_is_missing(gt_arr[i*2+1]) )
		    {
			if(bcf_gt_allele(gt_arr[i*2])>0||bcf_gt_allele(gt_arr[i*2+1])>0)
			{
			    stringstream ss;
			    ss  << hdr->samples[i] <<"\t"<<bcf_hdr_int2id(hdr,BCF_DT_CTG,gene_it->getChrom()) <<"\t"<<line->pos<<"\t"<<line->d.allele[0]<<"\t"<<line->d.allele[1]<<bcf_gt_allele(gt_arr[i*2])<<bcf_gt_allele(gt_arr[i*2+1])<<"\t"<<gene_it->getGeneName()<<endl;
			    lof_gts[i].push_back(ss.str());
//			    cerr << hdr->samples[i] <<"\t"<<bcf_hdr_int2id(hdr,BCF_DT_CTG,gene_it->getChrom()) <<"\t"<<line->pos<<"\t"<<line->d.allele[0]<<"\t"<<line->d.allele[1]<<bcf_gt_allele(gt_arr[i*2])<<bcf_gt_allele(gt_arr[i*2+1])<<"\t"<<gene_it->getGeneName()<<endl;
			}
		    }
		}
	    }
	}
    }

    free(gt_arr);
    if(a.include!=NULL)  filter_destroy(filter);
    bcf_sr_destroy(sr);	
    return(0);
}


int knockout_main(int argc,char **argv)
{
    int c;
    args a;
    if(argc<3) usage();
    static struct option loptions[] = {
	{"pedigree",1,0,'p'},
	{"include",1,0,'i'},
	{"genes-file",required_argument,NULL,'G'},
	{0,0,0,0}
    };
    a.pedigree=a.inputfile=a.include=a.gene_bedfile=NULL;

    while ((c = getopt_long(argc, argv, "p:i:G:",loptions,NULL)) >= 0)
    {
	switch (c)
	{
	case 'p': a.pedigree = optarg; break;
	case 'i': a.include = optarg; break;
	case 'G': a.gene_bedfile = optarg; break;
	}
    }
    optind++;
    a.inputfile=argv[optind];
    if(a.inputfile==NULL) die("no input provided");
    if(a.pedigree==NULL) die("the -p argument is required");
    if(a.gene_bedfile==NULL) die("the -G argument is required");
    knockout(a);
    return(0);
}

