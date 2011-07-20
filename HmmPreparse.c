/* Matt.c -- Multiple alignment program (Multiple Alignment with Translations and Twists).
 *
 * Author: Matt Menke (2007)
 *
 * Matt is licensed under the GNU public license version 2.0. 
 *
 * If you would like to license Matt in an enviroment where the GNU public license is 
 * unacceptable (such as inclusion in a non-GPL software package) comercial Matt
 * licensing is available through the MIT and Tufts offices of Technology Transfer. 
 * Contact betawrap@csail.mit.edu or cowen@cs.tufts.edu for more information.
 */

#ifdef _DEBUG
/* Memory leak detection library. */
/* #include <vld.h> */
#endif
#include <stdarg.h>
#include <stdio.h>
// The following library is included to use options
#include <unistd.h>

#ifdef _OPENMP
#include <omp.h>
#define THREADING "OpenMP"
#else
#define THREADING "Single Threaded"
#endif

#include "Score.h"
#include "Protein.h"
#include "pdb.h"
//#include "MultipleAlignment.h"
//#include "MultipleAlignmentOutput.h"
//#include "Extend.h"
/*
#define OUTPUT_PDB		0x01
#define OUTPUT_TXT		0x02
#define OUTPUT_FASTA	0x04
#define OUTPUT_SPT		0x08
#define OUTPUT_MSF		0x10
#define OUTPUT_PIR		0x20
#define OUTPUT_SPDBV	0x40

const char formats[][6] = {"pdb", "txt", "fasta", "spt", "msf", "pir", "spdbv"};
//*/

struct HELIX_COLLECTION {
	AlphaHelix *helices;
	int length;
};
typedef struct HELIX_COLLECTION Helices;


#define VERSION "1.00"

void Usage() {
	printf("\nMatt version " VERSION " " THREADING " build.\n"
		"See http://matt.csail.mit.edu for more information.\n\n"
		"Usage: Matt -o outprefix [-c cutoff] [-t threads] [-[rslbdvp][01]]*\n"
		"            [-f <extension list>] [-L listfile]* [file[:<chain list>]]*\n"
	);
}

void Feedback(FILE *out, char *prefix, char *string, va_list arg, char *suffix) {
	fprintf(out, "%s", prefix);

	vfprintf(out, string, arg);

	fprintf(out, "%s", suffix);
}

/* Not necessarily accurate counts in OpenMP mode, but doesn't really matter.
 */
int errors = 0;
int warnings = 0;

void WarningLog(int console, FILE *f2, char *string, ...) {
	va_list arg;
	warnings ++;

	if (console) {
		va_start(arg, string);
		Feedback(stdout, " * ", string, arg, "\n");
		va_end(arg);
	}

	if (f2) {
		va_start(arg, string);
		Feedback(f2, " * ", string, arg, "\n");
		va_end(arg);
	}
}

void StatusLog(int console, FILE *f2, char *string, ...) {
	va_list arg;

	if (console) {
		va_start(arg, string);
		Feedback(stdout, "   ", string, arg, "\n");
		va_end(arg);
	}

	if (f2) {
		va_start(arg, string);
		Feedback(f2, "   ", string, arg, "\n");
		va_end(arg);
	}
}

void ErrorLog(int console, FILE *f2, char *string, ...) {
	va_list arg;
	errors++;

	if (console) {
		va_start(arg, string);
		Feedback(stdout, "** ", string, arg, "\n");
		va_end(arg);
	}

	if (f2) {
		va_start(arg, string);
		Feedback(f2, "** ", string, arg, "\n");
		va_end(arg);
	}
}

void RootStatusLog(int console, FILE *f2, char *string, ...) {
	va_list arg;

	if (console) {
		va_start(arg, string);
		Feedback(stdout, "", string, arg, "\n");
		va_end(arg);
	}

	if (f2) {
		va_start(arg, string);
		Feedback(f2, "", string, arg, "\n");
		va_end(arg);
	}
}

struct ALIGNED_RESIDUE_PAIR {
	int res1, res2;
	int dir;
	int count;
};
typedef struct ALIGNED_RESIDUE_PAIR AlignedResiduePair;

struct ALIGNED_STRAND_PAIR {
	int res1, res2;
	int dir;
	int length;
};
typedef struct ALIGNED_STRAND_PAIR AlignedStrandPair;

// Define Table structure
typedef struct {
	char id;
	double scores[21][21];
} Table;
Table *tables = 0;
int numTables = 0;


// Load Beta-mutation table
int LoadTables() {
	int res = 1;
	char *line = 0;
	int len = 0;
	double bg[21] = {
		0.0787945,		/* A */
		0.0151600,		/* C */
		0.0535222,		/* D */
		0.0668298,		/* E */
		0.0397062,		/* F */
		0.0695071,		/* G */
		0.0229198,		/* H */
		0.0590092,		/* I */
		0.0594422,		/* K */
		0.0963728,		/* L */
		0.0237718,		/* M */
		0.0414386,		/* N */
		0.0482904,		/* P */
		0.0395639,		/* Q */
		0.0540978,		/* R */
		0.0683364,		/* S */
		0.0540687,		/* T */
		0.0673417,		/* V */
		0.0114135,		/* W */
		0.0304133,		/* Y */
		1
	};
	FileReader *in;
	int i,j=0;
	numTables = 0;
	if (tables) free(tables);
	tables = (Table*) calloc(3, sizeof(Table));
	in = readerOpen("betaDataBuriedRevised.dat");
	if (!in) res = 0;
	else {
		Table *table = &tables[numTables++];
		table->id = 'i';
		while(readerGetLine(in, &line, &len)) {
			char r1, r2;
			if (line[0] < 'A' || line[0] >'Z' ||
				line[1] < 'A' || line[1] >'Z') {
					res = 0;
					continue;
			}
			r1 = lookupTable[line[0]-'A'];
			r2 = lookupTable[line[1]-'A'];
			if (r1 > 20) r1 = 20;
			if (r2 > 20) r2 = 20;
			if (!sscanf(line+2, "%lf", &table->scores[r1][r2])) break;
		}
		readerClose(in);
	}

	in = readerOpen("betaDataExposedRevised.dat");
	if (!in) res = 0;
	else {
		Table *table = &tables[numTables++];
		table->id = 'o';
		while(readerGetLine(in, &line, &len)) {
			char r1, r2;
			if (line[0] < 'A' || line[0] >'Z' ||
				line[1] < 'A' || line[1] >'Z') {
					res = 0;
					continue;
			}
			r1 = lookupTable[line[0]-'A'];
			r2 = lookupTable[line[1]-'A'];
			if (r1 > 20) r1 = 20;
			if (r2 > 20) r2 = 20;
			if (!sscanf(line+2, "%lf", &table->scores[r1][r2])) break;
		}
		readerClose(in);
	}

	free(line);
	if (res) {
		int q;
		for (q=0; q<numTables; q++)
		{
			Table *table = tables+q;
			double r[21];
			double t = 0;
			for (i=0; i<20; i++) {
				r[i] = exp(table->scores[i][0])/exp(table->scores[0][i]);
				t += r[i];
			}
			for (i=0; i<20; i++) {
				r[i] /= t;
			}
			for (i=0; i<20; i++) {
				for (j=0; j<20; j++) {
					table->scores[i][j] = exp(table->scores[i][j]) * r[j];
				}
			}
			for (i=0; i<20; i++) {
				for (j=i; j<20; j++) {
					table->scores[j][i] = table->scores[i][j] =
						log((table->scores[i][j] + table->scores[j][i]) / 2 / bg[i]/bg[j]);
				}
			}
		}
	}
	return res;
}


// Cleaning up pool of alpha helices
void CleanupHelices(PDBModel *model) {
	// If more than two alpha helices share same region, manipulate indices and extend lengths
	// Total size of pool will be reduced if duplicated alpha helices are found
	int nSeqs;
	for(nSeqs = 0; nSeqs < model->numChains ; nSeqs++) {
		PDBChain *chain = model->chains[nSeqs];
		int i=0;
		while(i < chain->numAlphaHelices) {
			AlphaHelix helix1 = chain->alphaHelices[i];
			AlphaHelix helix2 = chain->alphaHelices[i+1];
			
			if( ((helix1.start <= helix2.start) && (helix2.start <= helix1.start + helix1.length) && (helix1.type == helix2.type)) 
			 || ((helix2.start <= helix1.start) && (helix1.start <= helix2.start + helix2.length) && (helix1.type == helix2.type)) ) { 
		
				int end = ((helix1.start + helix1.length) > (helix2.start + helix2.length)) ? (helix1.start + helix1.length) : (helix2.start + helix2.length);
				AlphaHelix helix;
				helix.start = (helix1.start < helix2.start) ? helix1.start : helix2.start;
				helix.length = end - helix.start;
	
				helix.type = helix1.type;	
				helix.id[0] = helix.comment[0] = 0;
				chain->alphaHelices[i] = helix;
				
				int j;
				for(j=i+1;j<chain->numAlphaHelices;j++) {
					chain->alphaHelices[j] = chain->alphaHelices[j+1];
				}
				
				chain->numAlphaHelices--;
			}
			else i++;
		}
	}
}


// Have same functionality as CleanupHelices(PDBChain *chain) 
// Defined for different intput data type
// Cleaning up pool of alpha helices
void TempCleanupHelices(Helices* h) {
	// Calculate total number of elements in helices

	// If more than two alpha helices share same region, manipulate indices and extend lengths
	// Total size of pool will be reduced if duplicated alpha helices are found
	int i=0;
	while(i < h->length) {
		AlphaHelix helix1 = h->helices[i];
		int j;
		for(j=i+1;j<h->length;j++) {
			AlphaHelix helix2 = h->helices[j];
			
			if( ((helix1.start <= helix2.start) && (helix2.start <= helix1.start + helix1.length) && (helix1.type == helix2.type))
			 || ((helix2.start <= helix1.start) && (helix1.start <= helix2.start + helix2.length) && (helix1.type == helix2.type)) ) { 
				int end = ((helix1.start + helix1.length) > (helix2.start + helix2.length)) ? (helix1.start + helix1.length) : (helix2.start + helix2.length);
				AlphaHelix helix;
				helix.start = (helix1.start < helix2.start) ? helix1.start : helix2.start;
				helix.length = end - helix.start;
				helix.type = helix1.type;
				helix.id[0] = helix.comment[0] = 0;
				h->helices[i] = helix;
				
				int k;
				for(k=j;k<h->length;k++) {
					h->helices[k] = h->helices[k+1];
					
					// To put helices[j+1] out of consideration later
					h->helices[k+1].start = 0;
					h->helices[k+1].length = 0;
					h->helices[k+1].type = 0;
				}
				h->length--;
			}
		}
		i++;
	}
	
	// Final cleanup: Remove all duplicated helices
	i=0;
	while( i < h->length ) {
		if( (h->helices[i].start == h->helices[i+1].start) && 
		    (h->helices[i].length == h->helices[i+1].length) && 
		    (h->helices[i].type == h->helices[i+1].type) ) {
			h->helices[i+1].start = 0;
			h->helices[i+1].length = 0;
			h->helices[i+1].type = 0;
			h->length--;	
		}
		i++;
	}
	h->length--;

}


// Reindex alpha annotations considering gaps from results of multiple alignment
void ReindexingHelices(PDBModel *model, SequenceAlignment *seqs) {
     int i,j;
     int action = 0;
     int s = 0;
     int e = 0;
     int len = 0;
     int acc_len = 0;
     char previous = 1;

	 int nSeqs;
	 for(nSeqs=0 ; nSeqs<model->numChains ; nSeqs++) {
		 PDBChain *chain = model->chains[nSeqs];
	
		 for(i=0;i<seqs->seqs[nSeqs]->length;i++) {
			char c = seqs->seqs[nSeqs]->seq[i];
			if(c<0) { // Found gap
				if(previous < 0) { // Middle of gap
					len++;
					e = s + len;
					previous = c;	  
				}
				else { // first gap
				   s = i;
				   previous = c;
				   len++;	   	 
				}
			}
			else { // Found residue	
				// After identifying gaps, push indices and lengths of alpha helices accordingly
				if(previous < 0) { // This residue is right after a gap
					for(j=0;j<chain->numAlphaHelices;j++) {
						if(s <= chain->alphaHelices[j].start) {
							chain->alphaHelices[j].start += len;            		
						}
						else if(s < (chain->alphaHelices[j].start + chain->alphaHelices[j].length) ) {
							chain->alphaHelices[j].length += len;
						}
						
					}
				}
				len = 0;
				previous = c;
			}
		}

	}
	CleanupHelices(model);
}


// Filtering alpha helices by deciding residues in consensus using ratio of amino acids
void HelixResidueRatio(PDBModel *model, SequenceAlignment *seqs, Helices *hout, float threshold) {
	int i,j;
	int max_len = 0;
	
	// Get maximum length of chain
	for(i=0 ; i<seqs->numSeqs ; i++) {
		if(seqs->seqs[i]->length > max_len) max_len = seqs->seqs[i]->length;
	}
	
	// Construct matrix
	int* matrix = (int*) calloc(1, max_len*sizeof(int));

	// If residue is found at index i, increment value in the matrix[i]
	// Value in matrix[i] cannot exceed total number of sequences
	// Value in matrix[i] indicates how many chains have residue at index of i
	for(i=0;i<seqs->numSeqs;i++) {
		for(j=0;j<seqs->seqs[i]->length;j++) {
			char c = seqs->seqs[i]->seq[j];
			if( 'A' <= ShortNames[c] && ShortNames[c] <= 'Z' ) {
				matrix[j]++;
			}		
		}
	}
	
	int hnum = 0;
	hout->helices = (AlphaHelix*) malloc(sizeof(AlphaHelix*));
	
	for(i=0;i<seqs->numSeqs;i++) {
		PDBChain *chain = model->chains[i];
		for(j=0;j<chain->numAlphaHelices;j++) {
			int ind1 = chain->alphaHelices[j].start;
			int ind2 = chain->alphaHelices[j].start + chain->alphaHelices[j].length;

			int k;
			int s = ind1; 
			int len = 0;

			// For each region that alpha helix lies, check if enough number of amino acids (determined by threshold) exist
			// If so, add the alpha helix to result
			for(k=ind1;k<ind2;k++) {
				float t = (float)(matrix[k]) / (float)(seqs->numSeqs);
				if(t >= threshold) {
					len++;	
					// Whole region of alpha helix is in consensus
					if(len == chain->alphaHelices[j].length) {
						AlphaHelix helix;
						helix.start = s;
						helix.length = len;
						helix.type = chain->alphaHelices[j].type;
						helix.id[0] = helix.comment[0] = 0;
																		
						hnum++;
						//out = (AlphaHelix*) realloc(out, hnum*sizeof(AlphaHelix));
						hout->helices = (AlphaHelix*) realloc(hout->helices, hnum*sizeof(AlphaHelix));
						hout->helices[hnum-1] = helix;
						len = 0;
					}
				}
				else {
					// Part of alpha helix is in consensus
					// Only store ones with length > 3
					if(len>3) {
						AlphaHelix helix;
						helix.start = s;
						helix.length = len;
						helix.type = chain->alphaHelices[j].type;
						helix.id[0] = helix.comment[0] = 0;
						
						hnum++;
						hout->helices = (AlphaHelix*) realloc(hout->helices, hnum*sizeof(AlphaHelix));
						hout->helices[hnum-1] = helix;
					}
					s = k + 1;
					len = 0;
				}
			}
		}
	}
	
	hout->length = hnum;
	free(matrix);
}

// Filtering alpha helices by deciding residues in consensus using ratio of alpha helices
void HelixAlphaRatio(PDBModel *model, SequenceAlignment *seqs, Helices *hout, float threshold) {
	int i,j,k;
	int max_len = 0;
	
	// Get maximum length of chain
	for(i=0 ; i<seqs->numSeqs ; i++) {
		if(seqs->seqs[i]->length > max_len) max_len = seqs->seqs[i]->length;
	}

	// Construct matrix
	int* matrix = (int*) calloc(1, max_len*sizeof(int));

	// For each chain,
	// increment value in the matrix[i], if one residue of alpha helix exists at index i  
	// Value in matrix[i] cannot exceed total number of sequences
	// Value in matrix[i] indicates how many chains have part of alpha helix at index of i
	for(i=0;i<seqs->numSeqs;i++) {
		PDBChain *chain = model->chains[i];
		for(j=0;j<chain->numAlphaHelices;j++) {
			int ind1 = chain->alphaHelices[j].start;
			int ind2 = chain->alphaHelices[j].start + chain->alphaHelices[j].length;

			for(k=ind1;k<ind2;k++) {
				char c = seqs->seqs[i]->seq[k];
				if( 'A' <= ShortNames[c] && ShortNames[c] <= 'Z' ) {
					matrix[k]++;
				}		
			}			
		}
	}
	
	int hnum = 0;
	hout->helices = (AlphaHelix*) malloc(sizeof(AlphaHelix*));
	
	int beginning = 0;
	int s = 0;
	int len = 0;
	for(i=0;i<max_len;i++) {
		float t = (float)(matrix[i]) / (float)(seqs->numSeqs); 
		if(t >= threshold) { // Found first residue of potential alpha helix in consensus
			if(beginning == 0) {
				s = i;
				beginning = 1;
			}
			len++;	
		}			
		else { // Residue at index i of any chain cannot be part of potential alpha helix in consensus
			beginning = 0;
			// Store alpha helix with length > 3
			if(len > 3) {
				AlphaHelix helix;
				helix.start = s;
				helix.length = len;
				
				// Find correct type of alpha helix
				int type = 0;
				j = 0;
				
				
				// Inifinite loop here
				while(type == 0 && j<seqs->numSeqs) {
					PDBChain *chain = model->chains[j];
					for(k=0;k<chain->numAlphaHelices-1;k++) {
						if( (s >= chain->alphaHelices[k].start) && (s+len <= (chain->alphaHelices[k+1].start + chain->alphaHelices[k+1].length)) ) 
						{
							type = chain->alphaHelices[k].type;
							break;
						}
					}	
				}
				
				helix.type = type;
				helix.id[0] = helix.comment[0] = 0;
																
				hnum++;
				hout->helices = (AlphaHelix*) realloc(hout->helices, hnum*sizeof(AlphaHelix));
				hout->helices[hnum-1] = helix;
			}
			len = 0;
		}
	}
	hout->length = hnum;
	free(matrix);
}


// Print error messages when wrong options are used
void printOption() 
{
	printf("Usage: [executable] <protein name>\n");
	printf("Usage: [executable] -m [alpha|residue] <protein name>\n");
	printf("Usage: [executable] -t [consensus threshold =< 1 ] <protein name>\n");
	printf("Usage: [executable] -o [alpha|beta] <protein name>\n");
	printf("Multiple options can be used. \n");
}

// CDECL is needed to avoid a warning under MSVC when using fastcall.
// Defined to be nothing in other environments.
int CDECL main(int realArgc, char **realArgv) {

	if (realArgc <= 1) {
	// When the program has no input
		printf("Nothing to do.\n");
		return 0;
	}


////////////////////////////////////////////////////////////////////////////////////////////////	

	// Flags for printing annotations 
	// If alpha=1, alpha annotation will be printed. If beta=1, beta annotation will be printed
	// Default is to print both annotations
	int alpha = 1;
	int beta = 1;
	
	// Threshold for consensus between sequences
	float threshold = 0.5;
	
	// Flags for method that calculates helices in consensus
	int method = 0; // when method is 0, calculate alpha helices in consensus using ratio of alpha helices
					// when method is 1, calculate alpha helices in consensus using ratio of amino acids
	
	// Get option for alpha/beta/none
	int c = 0;
	char *option = "";
	
	while ( (c = getopt(realArgc, realArgv, "m:o:t:") ) != -1 ) {
		switch (c) {
		case 'm': 
			option = optarg;
			if(	strcmp(option, "alpha") == 0 ) method = 0;
			else if( strcmp(option, "residue") == 0 ) method = 1;
			else {
				printOption();
				return 1;
			}
			break;
		case 'o':
			option = optarg;
			if( strcmp(option, "alpha") == 0 ) beta = 0;			
			else if( strcmp(option, "beta") == 0 )  alpha = 0;
			else {
				printOption();
				return 1;
			}
			break;
		case 't':
			threshold = atof(optarg);
			if(threshold > 1) {
				printOption();
				return 1;
			}
			break;
		case '?':
			printf("Unrecognizable option \n");
			return 1;
		default:
			return 1;
		}
	}
	

///////////////////////////////////////////////////////////////////////////////////////////////

	//test\\MattAlignment 
	int a;
	for (a=optind; a<realArgc; a++) {
		SequenceAlignment *seqs;
		PDBData *pdb;
		PDBModel *model;
		int **seqIndices;
		int **chainIndices;
		int residues = 0;
		int i, j, k;
		int *indices=0;
		int p1, p2;

		AlignedResiduePair * pairs = 0;
		int numPairs = 0;

		AlignedStrandPair * strands = 0;
		int numStrands = 0;

		int MIN_CONTINUE;
		int MIN_START;

		char *fileName = (char*) malloc(strlen(realArgv[a]) + 50);
		sprintf(fileName, "%s.fasta", realArgv[a]);
	    if(!(seqs = LoadMultipleAlignment(fileName))) {
			printf("%s.fasta does not exist.\n", realArgv[a]);
			return 0;
		}

		sprintf(fileName, "%s.pdb", realArgv[a]);
        if(!(pdb = LoadPDBFile(fileName, 1, 2, 0, 1))) {
			printf("%s.pdb does not exist.\n", realArgv[a]);
			return 0;
		}
		model = pdb->models[0];

		seqIndices = (int**) malloc(2*seqs->numSeqs * sizeof(int*));
		chainIndices = seqIndices+seqs->numSeqs;

		MIN_CONTINUE = seqs->numSeqs/2+1;
		MIN_START = MIN_CONTINUE;

		LoadTables();

		for (i=0; i<seqs->numSeqs; i++) {
			CalcChainSecondary(model->chains[i]);
			residues += seqs->seqs[i]->length;
			residues += model->chains[i]->length;
		}
		indices = (int *)malloc(residues * sizeof(int));
		for (i=0; i<seqs->numSeqs; i++) {
			seqIndices[i] = indices;
			chainIndices[i] = indices += seqs->seqs[i]->length;
			indices += model->chains[i]->length;
		}


////////////////////////////////////////////////////////////////////////////////////////////////
		
		// Reindexing alpha helices considering gaps in multiple alignment
		ReindexingHelices(model, seqs);
		
		// Deciding alpha helices in consensus
		Helices *halpha = (Helices*) malloc(sizeof(Helices*));
		if(method == 0) {
			HelixAlphaRatio(model, seqs, halpha, threshold);
		}
		else if(method == 1) {
			HelixResidueRatio(model, seqs, halpha, threshold);
			TempCleanupHelices(halpha);
		}
		
		
///////////////////////////////////////////////////////////////////////////////////		

		// Output Alpha helix information
		// Right-handed Alpha, type = 1
		// Left-handed Alpha, type = 6
		// 310, type = 5

		FILE *out;
		sprintf(fileName, "%s.ssi", realArgv[a]);
		out = fopen(fileName, "wb");

		if (out == NULL) {
		  fprintf(stderr, "Can't open output file %s!\n",
				  fileName);
		  return 1;
		}

		fprintf(out, "# STOCKHOLM 1.0\n\n");
		fprintf(out, "#=ALPHA type description: RIGHT_ALPHA = 1;\n");
		fprintf(out, "#=ALPHA LEFT_ALPHA.ll = 6; RIGHT_310 = 5;\n");
		fprintf(out, "#=ALPHA format: start length type\n\n");
		
		// Print alpha annotations only if alpha = 1
		for(i=0;i<halpha->length;i++) {
			fprintf(out, "#=ALPHA %d %d %d\n", halpha->helices[i].start, halpha->helices[i].length, halpha->helices[i].type);			
		}
		fprintf(out,"\n\n");
		
		
		// Free memories
		if(halpha) {
			free(halpha->helices);
			free(halpha);
		}


////////////////////////////////////////////////////////////////////////////////////		
				
		// Beta strands annotation starts from here
		// Adding modularity can be done later
		for (i=0; i<seqs->numSeqs; i++) {
			Sequence *seq = seqs->seqs[i];
			PDBChain *chain = model->chains[i];
			p1 = p2 = 0;
			while (1) {
				if (p1 == seq->length) break;
				if (seq->seq[p1] == -1) {
					seqIndices[i][p1] = -1;
					p1++;
					continue;
				}
				if (p2 == chain->length || seq->seq[p1] != chain->seq->seq[p2]) break;
				seqIndices[i][p1] = p2;
				chainIndices[i][p2] = p1;
				p1++;
				p2++;
			}
			if (p1 != seq->length || p2 != chain->length) {
				p1=p1;
			}
		}
		while (1) {
			char *markForDeath = 0;
			BetaResiduePair *betaPairs = (BetaResiduePair*) malloc(2 * sizeof(BetaResiduePair) * seqs->numSeqs);
			for (p1=0; p1<seqs->seqs[0]->length; p1++) {
				int count;
				int numBetaPairs = 0;
				for (i=0; i<seqs->numSeqs; i++) {
					Sequence *seq = seqs->seqs[i];
					PDBChain *chain = model->chains[i];
					count = 0;
					if (seqIndices[i][p1] == -1) continue;				
					for (j=0; j<chain->numBetaPairs; j++) {
						BetaResiduePair pair;
						if (FindMatch(chain->betaPairs+j, seqIndices[i][p1], &pair)) {
							// res1 is always seqIndices[i][p1]
							if (pair.res1 < pair.res2 && chainIndices[i][pair.res2] != -1) {
								k = numBetaPairs++;
								pair.res1 = chainIndices[i][pair.res1];
								pair.res2 = chainIndices[i][pair.res2];
								while (k &&
									(pair.res2 < betaPairs[k-1].res2 ||
									(pair.res2 == betaPairs[k-1].res2 && pair.dir < betaPairs[k-1].dir))) {
									betaPairs[k] = betaPairs[k-1];
									k--;
								}
								betaPairs[k] = pair;
							}
							// Should never happen.  Check just in case, because of memory allocation.
							if (++count == 2) break;
						}
					}
				}
				for (i=0; i<numBetaPairs; i++) {
					count = 1;
					while ( i+1 < numBetaPairs &&
						betaPairs[i].res2 == betaPairs[i+1].res2 &&
						betaPairs[i].dir == betaPairs[i+1].dir ) {
							count ++;
							i++;
					}
					if (count < MIN_START) continue;
					//if (betaPairs[i].res1 == 974 || betaPairs[i].res1 == 976 || betaPairs[i].res2 == 974 || betaPairs[i].res2 == 976) {
					//	i=i;
					//}
					if (numPairs % 256 == 0) {
						pairs = (AlignedResiduePair*) realloc(pairs, (numPairs+256)*sizeof(AlignedResiduePair));
					}
					k = numPairs++;
					while (k > 0 && pairs[k-1].count < count) {
						pairs[k] = pairs[k-1];
						k--;
					}
					pairs[k].res1 = betaPairs[i].res1;
					pairs[k].res2 = betaPairs[i].res2;
					pairs[k].dir = betaPairs[i].dir;
					pairs[k].count = count;
				}
			}
			free(betaPairs);
			for (i=0; i<numPairs; i++) {
				int dir;
				int *newPos;
				AlignedStrandPair * strand;
				for (j=0; j<numStrands; j++) {
					if (strands[j].res1 <= pairs[i].res1 && pairs[i].res1 < strands[j].res1 + strands[j].length) {
						if (strands[j].dir == 1 && strands[j].res2 <= pairs[i].res2 && pairs[i].res2 < strands[j].res2 + strands[j].length) {
							break;
						}
						if (strands[j].dir == -1 && strands[j].res2-strands[j].length < pairs[i].res2 && pairs[i].res2 <= strands[j].res2) {
							break;
						}
					}
				}
				if (j<numStrands) continue;
				newPos = (int*)malloc(2 * model->numChains * sizeof(int));
				strands = (AlignedStrandPair*)realloc(strands, (numStrands+1) * sizeof(AlignedStrandPair));
				strand = strands + numStrands;
				strand->length = 1;
				strand->res1 = pairs[i].res1;
				strand->res2 = pairs[i].res2;
				strand->dir = pairs[i].dir;
				for (dir = 1; dir > -2; dir -= 2) {
					int p1old = strand->res1;
					int p2old = strand->res2;
					while (1) {
						int count;
						int best;
						int p1, p2;

						/*
						p1 = p1old + dir;
						p2 = p1old + dir*strand->dir;

						if (p1 < 0 || p2 < 0 || p1 >= seqs->seqs[0]->length || p2 >= seqs->seqs[0]->length) break;
						//*/

						for (j=0; j<model->numChains; j++) {
							PDBChain *chain = model->chains[j];
							newPos[2*j] = -1;
							if (seqIndices[j][p1old] == -1 || seqIndices[j][p2old] == -1) continue;
							p1 = p1old;
							p2 = p2old;
							while (1) {
								p1 += dir;
								if (p1 < 0 || p1 >= seqs->seqs[0]->length) break;
								if (seqIndices[j][p1] != -1) break;
							}
							if (seqIndices[j][p1] == -2) continue;
							while (1) {
								p2 += dir*pairs[i].dir;
								if (p2 < 0 || p2 >= seqs->seqs[0]->length) break;
								if (seqIndices[j][p2] != -1) break;
							}
							if (seqIndices[j][p2] == -2) continue;
							if (p1 < 0 || p2 < 0 || p1 >= seqs->seqs[0]->length || p2 >= seqs->seqs[0]->length) continue;
							if (seqIndices[j][p1old]+dir != seqIndices[j][p1] ||
								seqIndices[j][p2old]+dir * pairs[i].dir != seqIndices[j][p2]) continue;

							for (k=0; k<chain->numBetaPairs; k++) {
								BetaResiduePair pair;
								if (FindMatch(chain->betaPairs+k, seqIndices[j][p1old], &pair) &&
									p2old == chainIndices[j][pair.res2] &&
									pair.dir == pairs[i].dir &&
									FindMatch(chain->betaPairs+k, seqIndices[j][p1], &pair) &&
									p2 == chainIndices[j][pair.res2]) {
										newPos[j*2] = p1;
										newPos[j*2+1] = p2;
										//count ++;
										break;
								}
							}
						}
						//if (count < MIN_CONTINUE) break;
						best = 0;
						count = 0;
						for (j = 0; j<model->numChains; j++) {
							int current = 1;
							if (newPos[2*j] == -1) continue;
							for (k = j+1; k<model->numChains; k++) {
								if (newPos[2*j] == newPos[2*k] &&
									newPos[2*j+1] == newPos[2*k+1]) current ++;
							}
							if (current > count) {
								count = current;
								best = j;
							}
						}
						if (count < MIN_CONTINUE) break;
						p1 = newPos[best*2];
						p2 = newPos[best*2+1];
						if (p1old + dir != p1 || p2old + dir * pairs[i].dir != p2) {
							if (!markForDeath) {
								markForDeath = (char*) calloc(seqs->seqs[0]->length, sizeof(char));
							}
							for (j=p1old+dir; j != p1; j += dir) {
								markForDeath[j] = 1;
							}
							for (j=p2old+dir*pairs[i].dir; j != p2; j += dir*pairs[i].dir) {
								markForDeath[j] = 1;
							}
						}//*/
						strand->length++;
						if (dir < 0) {
							strand->res1 --;
							strand->res2 -= strand->dir;
						}
						p1old = p1;
						p2old = p2;
					}
				}
				if (strand->length > 1) {
					numStrands++;
				}
				free(newPos);
				if (markForDeath) break;
			}
			free(pairs);
			pairs = 0;
			numPairs = 0;
			if (markForDeath) {
			// Print nothing
				for (i=seqs->seqs[0]->length-1; i>=0; i--) {
					if (!markForDeath[i]) continue;
					for (j=0; j<model->numChains; j++) {
						int start = seqIndices[j][i];
						if (start >= 0) {
							chainIndices[j][start] = -1;
						}
						else {
							int w = i+1;
							start = model->chains[j]->length;
							while (w < seqs->seqs[j]->length) {
								if (seqIndices[j][w] >= 0) {
									start = seqIndices[j][w];
									break;
								}
								w++;
							}
						}
						for (k = start; k < model->chains[j]->length; k++) {
							if (chainIndices[j][k] >= 0) chainIndices[j][k]--;
						}
						for (k=i+1; k<seqs->seqs[j]->length; k++) {
							seqs->seqs[j]->seq[k-1] = seqs->seqs[j]->seq[k];
							seqIndices[j][k-1] = seqIndices[j][k];
						}
						seqs->seqs[j]->length--;
					}
				}
				//free(markForDeath);
				numStrands = 0;
				free(strands);
				strands = 0;
			}
			else {
				//FILE *out;
				//sprintf(fileName, "%s.ssi", realArgv[a]);
				//out = fopen(fileName, "wb");
				//fprintf(out, "# STOCKHOLM 1.0\n\n");

				
				for (i=0; i<numStrands; i++) {
					int j = i;
					AlignedStrandPair temp = strands[i];
					while (j>0 && strands[j-1].res1 > temp.res1) {
						strands[j] = strands[j-1];
						j--;
					}
					strands[j] = temp;
				}	
				
				// Print beta annotations only if beta = 1
				if(beta == 1) {
				
					// Print out beta strand annotation
					for (i=0; i<numStrands; i++) {
						int p1 = strands[i].res1;
						int p2 = strands[i].res2;
						int start2 = p2;
						int longest = 0;
						int chain;
						if (strands[i].dir == -1) start2 -= strands[i].length-1;
						for (chain=0; chain<model->numChains; chain++) {
							int len = 0;
							for (j=p1+strands[i].length; j<start2; j++) {
								if (seqs->seqs[chain]->seq[j] >= 0) len++;
							}
							if (len > longest) longest = len;
						}
						fprintf(out, "#=BETA %i\t%i\t%i\t%i\t%i\t", 1+strands[i].res1, 1+start2, strands[i].length, longest+20, strands[i].dir);
						for (j=0; j<strands[i].length; j++) {
							int chain;
							double bestScore = -90000;
							int bestTable = 0;
							int table;
							for (table=0; table<numTables; table++) {
								double score = 0;
								for (chain=0; chain<model->numChains; chain++) {
									int r1 = seqs->seqs[chain]->seq[p1];
									int r2 = seqs->seqs[chain]->seq[p2];
									// TODO: FIX
									if (r1 >= 20) r1 = 19;
									if (r2 >= 20) r2 = 19;
	
									if (r1 >= 0 && r2 >= 0) {
										score += tables[table].scores[r1][r2];
									}
								}
								if (score > bestScore) {
									bestTable = table;
									bestScore = score;
								}
							}
							fputc(tables[bestTable].id, out);
							p1++;
							p2+=strands[i].dir;
						}
						fputc('\n', out);
					}

				}
						
				for (i=0; i<model->numChains; i++) {
					fprintf(out, "\n%-20s ", seqs->seqs[i]->name);
					for (j=0; j<seqs->seqs[i]->length; j++) {
						char c = seqs->seqs[i]->seq[j];
						if (c<0) fputc('-', out);
						else fputc(ShortNames[c], out);
					}
				}
				fprintf(out, "\n//\n");
				fclose(out);
				numStrands = 0;
				free(strands);
				break;
				/*
				FILE *out;
				sprintf(fileName, "%s.fasta+", realArgv[a]);
				out = fopen(fileName, "wb");
				fprintf(out, "Fasta algnment augmented with beta strand info\n");
				for (i=0; i<numStrands; i++) {
					int j = i;
					AlignedStrandPair temp = strands[i];
					while (j>0 && strands[j-1].res1 > temp.res1) {
						strands[j] = strands[j-1];
						j--;
					}
					strands[j] = temp;
				}
				for (i=0; i<numStrands; i++) {
					int p1 = strands[i].res1;
					int p2 = strands[i].res2;
					int start2 = p2;
					int longest = 0;
					int chain;
					if (strands[i].dir == -1) start2 -= strands[i].length-1;
					for (chain=0; chain<model->numChains; chain++) {
						int len = 0;
						for (j=p1+strands[i].length; j<start2; j++) {
							if (seqs->seqs[chain]->seq[j] >= 0) len++;
						}
						if (len > longest) longest = len;
					}
					fprintf(out, "%i\t%i\t%i\t%i\t%i\t", 1+strands[i].res1, 1+start2, strands[i].length, longest+20+strands[i].length, strands[i].dir);
					for (j=0; j<strands[i].length; j++) {
						int chain;
						double bestScore = -90000;
						int bestTable = 0;
						int table;
						for (table=0; table<numTables; table++) {
							double score = 0;
							for (chain=0; chain<model->numChains; chain++) {
								int r1 = seqs->seqs[chain]->seq[p1];
								int r2 = seqs->seqs[chain]->seq[p2];
								// TODO: FIX
								if (r1 >= 20) r1 = 19;
								if (r2 >= 20) r2 = 19;

								if (r1 >= 0 && r2 >= 0) {
									score += tables[table].scores[r1][r2];
								}
							}
							if (score > bestScore) {
								bestTable = table;
								bestScore = score;
							}
						}
						fputc(tables[bestTable].id, out);
						p1++;
						p2+=strands[i].dir;
					}
					fputc('\n', out);
				}
				for (i=0; i<model->numChains; i++) {
					fprintf(out, "\n>%s\n", seqs->seqs[i]->name);
					for (j=0; j<seqs->seqs[i]->length; j++) {
						char c = seqs->seqs[i]->seq[j];
						if (j && !(j%60)) {
							fputc('\n', out);
						}
						if (c<0) fputc('-', out);
						else fputc(ShortNames[c], out);
					}
				}
				fclose(out);
				numStrands = 0;
				free(strands);
				break;
				//*/
			}
		}
		if (seqIndices) {
			free(seqIndices[0]);
			free(seqIndices);
		}
		FreeSequenceAlignment(seqs);
		CleanupPDB(pdb);
		free(fileName);
	}
	return 0;
}
