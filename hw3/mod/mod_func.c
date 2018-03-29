#include <stdio.h>
#include "hocdec.h"
#define IMPORT extern __declspec(dllimport)
IMPORT int nrnmpi_myid, nrn_nobanner_;

extern void _burste_reg();
extern void _bursti_reg();
extern void _fast_reg();
extern void _fasta_reg();
extern void _regular_reg();
extern void _regulara_reg();

modl_reg(){
	//nrn_mswindll_stdio(stdin, stdout, stderr);
    if (!nrn_nobanner_) if (nrnmpi_myid < 1) {
	fprintf(stderr, "Additional mechanisms from files\n");

fprintf(stderr," burste.mod");
fprintf(stderr," bursti.mod");
fprintf(stderr," fast.mod");
fprintf(stderr," fasta.mod");
fprintf(stderr," regular.mod");
fprintf(stderr," regulara.mod");
fprintf(stderr, "\n");
    }
_burste_reg();
_bursti_reg();
_fast_reg();
_fasta_reg();
_regular_reg();
_regulara_reg();
}
