/*
code developed and maintained by (jmw@ruc.edu.cn, RUC, China) date 2022 - 2025
*/

#include "specs.h"
#include "prmtr.h"
// #include "model.h"
// #include "api_zen.h"
#include "api_edmft.h"
#include "dmft_bethe.h"

int get_run_mode(const MyMpi& mm); // Function to read the run mode from PARAMS.norg
void read_params_norg(const MyMpi& mm, Prmtr& p, const Str& filename);
void set_control_divs_from_params(Prmtr& p, const std::vector<int>& restrain_values, const std::vector<int>& distribute_values);

int main(int argc, char* argv[])
{
	using namespace std;
	MPI_Init(&argc, &argv);
	MyMpi mm;
	if (mm) cout << "\n\nVersion: v1.10.01@2025.04.30\
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
        read_params_norg(mm, prmtr, "PARAMS.norg");
		DMFT dmft(mm, prmtr, 1);
        // WRN("DMFT is not implemented yet.");
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
            bool found_mode = false;
            
            // Read through file lines looking for "mode" entry
            while (std::getline(params_file, line) && !found_mode) {
                try {
                    // Remove any comments (anything after #)
                    size_t comment_pos = line.find('#');
                    if (comment_pos != std::string::npos) {
                        line = line.substr(0, comment_pos);
                    }
                    
                    // Trim whitespace from both ends
                    line.erase(0, line.find_first_not_of(" \t\n\r\f\v"));
                    line.erase(line.find_last_not_of(" \t\n\r\f\v") + 1);
                    
                    // Skip empty lines
                    if (line.empty()) continue;
                    
                    // Split by whitespace to get key and value
                    std::istringstream iss(line);
                    std::string key, value_str;
                    iss >> key;
                    
                    // If this is the mode line
                    if (key == "mode") {
                        found_mode = true;
                        iss >> value_str;
                        
                        // Convert to integer
                        mode = std::stoi(value_str);
                        
                        if (mode != 0 && mode != 1) {
                            std::cerr << "Error: Invalid mode value (" << mode << ") in PARAMS.norg. Must be 0 or 1." << std::endl;
                            mode = -1; // Set to error value if not 0 or 1
                        }
                    }
                } catch (const std::invalid_argument& ia) {
                    std::cerr << "Error: Invalid format in PARAMS.norg. Could not parse mode value from line: " << line << std::endl;
                    mode = -1; // Set mode to error value
                } catch (const std::out_of_range& oor) {
                    std::cerr << "Error: Mode value out of range in PARAMS.norg from line: " << line << std::endl;
                    mode = -1; // Set mode to error value
                }
            }
            
            if (!found_mode) {
                std::cerr << "Error: Could not find 'mode' entry in PARAMS.norg file." << std::endl;
                mode = -1;
            }
            
            params_file.close();
        }
    }

    // Broadcast the mode value from the master process to all other processes
    MPI_Bcast(&mode, 1, MPI_INT, 0, MPI_COMM_WORLD);

    return mode;
}


// Function to read parameters from PARAMS.norg file
void read_params_norg(const MyMpi& mm, Prmtr& p, const Str& filename) {
    std::ifstream file(filename);
    if (!file) {
        if (mm) WRN("Cannot open file: " + filename);
        return;
    }

    std::string line;
    int mode = 0;
    double mu = 0.0;
    int band = 0;
    int pred_gs_deg = 0;
    std::vector<double> t_values;
    double u = 0.0;
    double uprim = 0.0;
    std::vector<int> fit_nbaths_values;
    std::vector<int> restrain_values;
    std::vector<int> distribute_values;

    while (std::getline(file, line)) {
        std::istringstream iss(line);
        std::string key;
        iss >> key;

        if (key == "mode") {
            iss >> mode;
        } else if (key == "bethe_mu") {
            iss >> mu;
        } else if (key == "bethe_band") {
            iss >> band;
        } else if (key == "pred_gs_deg") {
            iss >> pred_gs_deg;
        } else if (key == "bethe_t") {
            char ch;
            while (iss >> ch && ch != ']') {
                if (ch != ',' && ch != '[') {
                    double value;
                    iss.unget();
                    iss >> value;
                    t_values.push_back(value);
                }
            }
        } else if (key == "bethe_u") {
            iss >> u;
        } else if (key == "bethe_uprim") {
            iss >> uprim;
        } else if (key == "fit_nbaths") {
            char ch;
            while (iss >> ch && ch != ']') {
                if (ch != ',' && ch != '[') {
                    int value;
                    iss.unget();
                    iss >> value;
                    fit_nbaths_values.push_back(value);
                }
            }
        } else if (key == "restrain") {
            char ch;
            while (iss >> ch && ch != ']') {
                if (ch == '0') {
                    restrain_values.push_back(0);
                }
                else if (ch == '-' || isdigit(ch)) {
                    iss.unget();
                    int value;
                    iss >> value;
                    restrain_values.push_back(value);
                }
            }
        } else if (key == "distribute") {
            char ch;
            while (iss >> ch && ch != ']') {
                if (ch == '0') {
                    distribute_values.push_back(0);
                }
                else if (ch == '-' || isdigit(ch)) {
                    iss.unget();
                    int value;
                    iss >> value;
                    distribute_values.push_back(value);
                }
            }
        }
    }

    // Update parameters in p

    p.mu = mu;
    p.nband = band;
    p.norbs = 2 * band;
    p.degel = pred_gs_deg;
    p.if_norg_degenerate = pred_gs_deg;
    p.bethe_t.reset(VecReal(t_values));     // Convert std::vector to VecReal
    // bethe_t = VecReal(t_values);         // Convert std::vector to VecReal
    p.U = u;
    p.delta = p.U - uprim;
    Real bethe_u = u;
    Real bethe_u12 = uprim;
    VecReal bethe_t = VecReal(t_values);
    
    // Update dependent parameters
    p.norg_sets = p.norbs;
    p.bandw = 2 * (2 * SQRT(SQR(bethe_u) + SQR(bethe_u12) + SQR(bethe_u12) + SUM(bethe_t * bethe_t)));

    p.derive_ImGreen();
    
    // Convert std::vector to VecReal/VecInt for output
    VecReal t_vals = VecReal(t_values);
    VecInt fit_nbaths_vals = VecInt(fit_nbaths_values);
    
    // Set restrain and distribute directly
    p.templet_restrain = VecInt(restrain_values);
    // p.distribute = VecInt(distribute_values);
    
    if (mm) WRN(NAV9(mu, band, t_vals, u, uprim, pred_gs_deg, fit_nbaths_vals, VecInt(restrain_values), VecInt(distribute_values)));

    // Call the function to set control_divs if distribute_values is not empty before calling after_modify_prmtr
    if (!distribute_values.empty()) {
        set_control_divs_from_params(p, restrain_values, distribute_values);
    }
    
    p.after_modify_prmtr(fit_nbaths_vals);
    p.recalc_partical_number(); p.derive();
}


// Function to set control_divs based on restrain and distribute values
void set_control_divs_from_params(Prmtr& p, const std::vector<int>& restrain_values, const std::vector<int>& distribute_values) {
    if (distribute_values.empty()) return;

    // Convert std::vector to VecInt using constructor that takes StdVec
    VecInt restrain = VecInt(restrain_values);
    VecInt distribute = VecInt(distribute_values);
		
    // Save the values to templet_restrain
    p.templet_restrain = restrain;
    
    // Setup control_divs based on input parameters
    p.control_divs.reset(p.norg_sets + 1, distribute.size(), 0);
    
    // First row is the restrain values
    p.control_divs[0] = restrain;
    
    // Other rows are the distribute values
    for_Int(i, 0, p.norg_sets) {
        p.control_divs[i + 1] = distribute;
    }
    
    // Set ndiv to the size of distribute
    p.ndiv = distribute.size();
}