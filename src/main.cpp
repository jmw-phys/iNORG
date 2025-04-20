/*
code developed and maintained by (jmw@ruc.edu.cn, RUC, China) date 2022 - 2024
*/

#include "specs.h"
#include "prmtr.h"
// #include "model.h"
// #include "api_zen.h"
#include "api_edmft.h"
// #include "dmft_bethe.h"

// Function to read the run mode from PARAMS.norg
// Returns the mode (0 or 1) or -1 on error.
int get_run_mode(const MyMpi& mm) {
    int mode = -1; // Default value, indicates error or not set

    if (mm) { // Only master process reads the file
        std::ifstream params_file("PARAMS.norg");
        if (!params_file.is_open()) {
            std::cerr << "Error: Could not open PARAMS.norg file." << std::endl;
            // mode remains -1
        } else {
            std::string line;
            if (std::getline(params_file, line)) {
                try {
                    // Assuming the line is "mode = value" or just "value"
                    size_t eq_pos = line.find('=');
                    std::string value_str;
                    if (eq_pos != std::string::npos) {
                        value_str = line.substr(eq_pos + 1);
                    } else {
                        // If '=' not found, assume the whole line is the value
                        value_str = line;
                    }
                    // Trim whitespace (optional, but good practice)
                    value_str.erase(0, value_str.find_first_not_of(" \t\n\r\f\v"));
                    value_str.erase(value_str.find_last_not_of(" \t\n\r\f\v") + 1);

                    mode = std::stoi(value_str);
                    if (mode != 0 && mode != 1) {
                         std::cerr << "Error: Invalid mode value (" << mode << ") in PARAMS.norg. Must be 0 or 1." << std::endl;
                         mode = -1; // Set to error value if not 0 or 1
                    }
                } catch (const std::invalid_argument& ia) {
                    std::cerr << "Error: Invalid format in PARAMS.norg. Could not parse mode value from line: " << line << std::endl;
                    mode = -1; // Set mode to error value
                } catch (const std::out_of_range& oor) {
                    std::cerr << "Error: Mode value out of range in PARAMS.norg from line: " << line << std::endl;
                    mode = -1; // Set mode to error value
                }
            } else {
                std::cerr << "Error: PARAMS.norg is empty or could not read the first line." << std::endl;
                mode = -1; // Set mode to error value
            }
            params_file.close();
        }
    }

    // Broadcast the mode value from the master process to all other processes
    MPI_Bcast(&mode, 1, MPI_INT, 0, MPI_COMM_WORLD);

    return mode;
}

int main(int argc, char* argv[])
{
	using namespace std;
	MPI_Init(&argc, &argv);
	MyMpi mm;
	if (mm) cout << "\n\nVersion: v1.9.16@2025.04.20\
(running "<< present() <<")\n\n" << endl;
	if (mm) cout << NAV(pwd()) << endl; 
	use_mkl(mm);
	// if (mm) io_init();
	clock_t t_program_bgn;
	if (mm) TIME_BGN("program", t_program_bgn);
	Prmtr prmtr(mm);

    // Get the run mode from the dedicated function
    int mode = get_run_mode(mm);

	if (mode == 0) {
		if (mm) cout << "Mode 0 selected: Running DMFT..." << endl;
		// DMFT dmft(mm, prmtr, 1);
        WRN("DMFT is not implemented yet.");
	} else if (mode == 1) {
		if (mm) cout << "Mode 1 selected: Running APIedmft..." << endl;
		APIedmft norg(mm, prmtr, "solver");
		// if(mm) WRN(NAV(prmtr.control_divs));
	} else {
		if (mm) { // Only master prints the final error message (it already printed details)
			cerr << "Error: Could not determine valid run mode from PARAMS.norg. Aborting." << endl;
		}
		MPI_Abort(MPI_COMM_WORLD, 1); // Abort MPI execution
	}

	// APIzen norg(mm, prmtr, "solver"); // Original line commented or removed
	// if(mm) WRN(NAV(prmtr.control_divs)); // Original line commented or removed
	// DMFT dmft(mm, prmtr, 1); // Original line commented or removed

    if (mm)	TIME_END("program", t_program_bgn);
    MPI_Finalize();
	return 0;
}