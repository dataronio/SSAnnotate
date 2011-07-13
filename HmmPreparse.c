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

#define VERSION "1.00"

void Usage() {
	printf("\nMatt version " VERSION " " THREADING " build.\n"
		"See http://matt.csail.mit.edu for more information.\n\n"
		"Usage: Matt -o outprefix [-c cutoff] [-t threads] [-[rslbdvp][01]]*\n"
		"            [-f <extension list>] [-L listfile]* [file[:<chain list>]]*\n"
	);
}

void Feedback(FILE *out, char *prefix, char *string, va_list arg, char *suffix) {
	fprintf(out, prefix);

	vfprintf(out, string, arg);

	fprintf(out, suffix);
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
	
	// Get option for alpha/beta/none
	int c;
	char *option;
	
	while ( (c = getopt(realArgc, realArgv, "o:") ) != -1 ) {
		switch (c) {
		case 'o':
			option = optarg;
			break;
		case '?':
			printf("Unrecognizable option \n");
			return 1;
		default:
			return 1;
		}
	}

	// To add more options, this part should be modified, too
	if( strcmp(option, "alpha") == 0 ) beta = 0;			
	else if( strcmp(option, "beta") == 0 )  alpha = 0;
	else {
		if(optind != 1) {
			printf("Usage: [executable] <protein name>\n");
			printf("Usage: [executable] -o [alpha|beta] <protein name>\n");
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
		seqs = LoadMultipleAlignment(fileName);

		sprintf(fileName, "%s.pdb", realArgv[a]);
		pdb = LoadPDBFile(fileName, 1, 2, 0, 1);
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

		
		
///////////////////////////////////////////////////////////////////////////////////		

		// Output Alpha helix information
		// Right-handed Alpha, type = 1
		// Left-handed Alpha, type = 6
		// 310, type = 5

		FILE *out;
		sprintf(fileName, "%s.ssi", realArgv[a]);
		out = fopen(fileName, "wb");

		int nSeqs;
		int nHelix;

		fprintf(out, "# STOCKHOLM 1.0\n\n");
		
		// Print alpha annotations only if alpha = 1
	    if(alpha == 1) {
	    
		// To be removed later
		// fprintf(out, "# ALPHA type description: RIGHT_ALPHA = 1; LEFT_ALPHA = 6;\n"),
		// fprintf(out, "# ALPHA format: start length type\n");
		// 
		
		// Remove redundant alpha helices
		PDBChain *chainA = model->chains[0];
		for(nSeqs = 1; nSeqs < seqs->numSeqs ; nSeqs++) 
		{
			PDBChain *chainB = model->chains[nSeqs];	
			int loop_start = chainA->numAlphaHelices - 1;
			int ia, ib;
			for(ib=chainB->numAlphaHelices-1 ; ib >= 0 ; ib--) 
			{
				for(ia=loop_start ; ia >=0 ; ia--) 
				{
					if(chainA->alphaHelices[ia].start == chainB->alphaHelices[ib].start) 
					{
						int j;
						for(j=ib+1 ; j < chainB->numAlphaHelices ; j++)
						{
							chainB->alphaHelices[j-1] = chainB->alphaHelices[j];
						}
						chainB->numAlphaHelices--;
						break;
					}
				}
			}
		}

		// Print out each alpha helices detected
		for (nSeqs = 0 ; nSeqs < seqs->numSeqs; nSeqs++ )
		{
			PDBChain *chain = model->chains[nSeqs];
			for (nHelix = 0; nHelix < chain->numAlphaHelices; nHelix++ )
			{
				fprintf(out, "#=ALPHA %d %d %d\n",
					chain->alphaHelices[nHelix].start, 
					chain->alphaHelices[nHelix].length, 
					chain->alphaHelices[nHelix].type);
			}
		}

		fprintf(out, "\n");		
		}

////////////////////////////////////////////////////////////////////////////////////		
		
		

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
					while (i+1 < numBetaPairs &&
						betaPairs[i].res2 == betaPairs[i+1].res2 &&
						betaPairs[i].dir == betaPairs[i+1].dir) {
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
				free(markForDeath);
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
