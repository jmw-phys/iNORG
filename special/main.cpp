
#include "specs.h"
#include "prmtr.h"
#include "model.h"
#include "apizen.h"
#include "dmft_bethe.h"

int main(int argc, char* argv[])
{
	using namespace std;
	MPI_Init(&argc, &argv);
	MyMpi mm;
	if (mm) cout << "\n\nVersion: v1.2.01 @ 2023.07.04\
(running "<< present() <<")\n\n" << endl;
	if (mm) cout << NAV(pwd()) << endl; 
	use_mkl(mm);
	// if (mm) io_init();
	clock_t t_program_bgn;
	if (mm) TIME_BGN("program", t_program_bgn);
	Prmtr prmtr(mm);

	APIzen norg(mm, prmtr, "solver");
	// DMFT dmft(mm, prmtr, 1);
	
    if (mm)	TIME_END("program", t_program_bgn);
    MPI_Finalize();
	return 0;
}