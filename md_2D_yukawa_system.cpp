//19.04.2021
//2D Yukawa system

#define _USE_MATH_DEFINES
#include <iostream>
#include <fstream>
#include <cmath>
#include <random>
#include <chrono>
#include <string>
#include <complex>
#include <vector>
// #include <windows.h>
using namespace std;

//set parameters
const int num_part = 10000;
int norm_thermo_time = 100000;
int norm_free_time = 10000;
int norm_meas_time = 100000;
int write_period = 100;
double time_step = 0.01;
double coupl_param = 200.0;
double screen_param = 2.0;
double friction_coef = 0.0;

double max_wave_num = 3.0;
double max_disp_freq = 1.5;

//file names
string title_of_coordinates = "XYZ ().xyz";
string title_of_energy = "ENERGY ().dat";
string title_of_rdf = "RDF ().dat";
string title_of_oop = "OOP ().dat";
string title_of_dispersion = "DISPERSION ().dat";
string title_of_msd = "MSD ().dat";
string title_of_vacf = "VACF ().dat";
string title_of_sacf = "SACF ().dat";

//functions
double Cutoff_Radius();
void Initialization();
void Force(int norm_time);
void Integration(int switch_on, int norm_time);
void Radial_Distribution_Function(int switch_on);
void Orientational_Order_Parameter(int norm_time);
void Dispersion(int switch_on);
void MSD_VACF_SACF(int switch_on);

//random number engine
random_device random_dev;
mt19937_64 Mersenne_Twister(random_dev());

//system parameters
double box_length = sqrt(M_PI * num_part);
int norm_total_time = norm_thermo_time + norm_free_time + norm_meas_time;
double r_c = Cutoff_Radius();
int num_cells_1D = floor(box_length / r_c);
double cell_length = box_length / num_cells_1D;
int num_cells_2D = num_cells_1D * num_cells_1D;
double random_force = sqrt(friction_coef / (coupl_param * time_step));
normal_distribution<double> normal_dist(0.0, 1.0);

double wave_num_step = 2.0 * M_PI / box_length;
double disp_freq_step = 2.0 * M_PI / (norm_meas_time * time_step);
int norm_max_wave_num = round(max_wave_num / wave_num_step);
int norm_max_disp_freq = round(max_disp_freq / disp_freq_step);

//global variables
//coordinates and velocity components
double x[num_part], y[num_part], v_x[num_part], v_y[num_part];

//previous, current, and next time force components
double F_x_prev[num_part], F_x_curr[num_part], F_x_next[num_part],
F_y_prev[num_part], F_y_curr[num_part], F_y_next[num_part];

//potential and kinetic energy functions
double kin_en, pot_en;

//linked list variables
int linked_list[num_part];
vector<int> head(num_cells_2D);

//RDF variables
int rdf_samples;
double ring_width;
vector<double> rdf;

//OOP variables
vector<int> num_neighbor_part;
vector<complex<double>> oop;

//Dispersion variables
vector<vector<complex<double>>> micro_long_current, micro_tran_current;

//MSD, VACF and SACF variables
int time_counter, num_blocks, num_block_elem, num_points_log10;
vector<int> block_length;
vector<vector<int>> msd_vacf_safc_samples;

vector<int> cross_x, cross_y;
vector<vector<double>> msd;
vector<vector<vector<double>>> x_time_0, y_time_0;

vector<vector<double>> vacf;
vector<vector<vector<double>>> v_x_time_0, v_y_time_0;

double kin_micro_stress_tensor, pot_micro_stress_tensor;
vector<vector<double>> sacf, micro_stress_tensor_0;

//MD CORE
int main()
{
	//keep track of time (start)
	auto start = chrono::high_resolution_clock::now();

	//initialize the system
	Initialization();

	//calc. previous and current time forces
	Force(-1);
	for (int i = 0; i < num_part; i++)
	{
		F_x_prev[i] = F_x_next[i];
		F_x_curr[i] = F_x_next[i];
		F_y_prev[i] = F_y_next[i];
		F_y_curr[i] = F_y_next[i];
	}

	//start the loop
	for (int norm_time = 0; norm_time < norm_total_time; norm_time++)
	{
		//calc. next time coordinates
		Integration(1, norm_time);

		//calc. next time forces
		Force(norm_time);

		//calc. next time velocities
		Integration(2, norm_time);

		//calc. properties of the system
		if (norm_time == norm_thermo_time + norm_free_time - 1)
		{
			Radial_Distribution_Function(0);
			Dispersion(0);
			MSD_VACF_SACF(0);
		}
		if (norm_time >= norm_thermo_time + norm_free_time)
		{
			Radial_Distribution_Function(1);
			Dispersion(1);
			MSD_VACF_SACF(1);
		}

		if (norm_time % write_period == 0)
		{
			//show energies on the screen
			cout << "time: " << (norm_time - (norm_thermo_time + norm_free_time)) * time_step << '\n'
				<< "pot_en: " << pot_en << '\n' << "coupl_param: " << 1.0 / kin_en << '\n';

			//write energies to the file
			ofstream energy_file(title_of_energy, ios::app);
			energy_file << (norm_time - (norm_thermo_time + norm_free_time)) * time_step << '\t'
				<< pot_en << '\t' << kin_en << '\t' << 1.0 / kin_en << '\n';
			energy_file.close();

			if (norm_time >= norm_thermo_time + norm_free_time)
			{
				//write coordinates to the file
				ofstream coordinate_file(title_of_coordinates, ios::app);
				coordinate_file << num_part << '\n' << "time: " << (norm_time - (norm_thermo_time + norm_free_time)) * time_step << '\n';
				for (int i = 0; i < num_part; i++)
					coordinate_file << i << '\t' << x[i] << '\t' << y[i] << '\n';
				coordinate_file.close();

				//calc. OOP
				Orientational_Order_Parameter(norm_time);
			}
		}
	}

	//write properties of the system to the file
	Radial_Distribution_Function(2);
	Dispersion(2);
	MSD_VACF_SACF(2);

	//keep track of time (finish)
	auto finish = chrono::high_resolution_clock::now();
	chrono::duration<double> duration = finish - start;
	ofstream simulation_time("simulation time.dat", ios::app);
	simulation_time << "simulation time is:\n" << duration.count() << " s\n"
		<< duration.count() / 60.0 << " min\n" << duration.count() / 3600.0 << " h\n";
	simulation_time.close();

	/*Beep(1568, 200);
	Beep(1568, 200);
	Beep(1568, 200);
	Beep(1245, 1000);
	Beep(1397, 200);
	Beep(1397, 200);
	Beep(1397, 200);
	Beep(1175, 1000);*/
}

//CUTOFF RADIUS
double Cutoff_Radius()
{
	double r = 1.0, accuracy = pow(10.0, -6.0);
	while (exp(-screen_param * r) / r >= accuracy)
		r += 0.01;
	cout << "r_c = " << r << '\n';
	return r;
}

//INITIALIZATION
void Initialization()
{
	//place particles in a 2D grid
	int n1 = 0;
	int n2 = 0;
	double stand_dev = sqrt(0.5 / coupl_param);
	normal_distribution<double> velocity_dist(0.0, stand_dev);
	for (int i = 0; i < num_part; i++)
	{
		if (n1 >= static_cast<int>(sqrt(num_part)))
		{
			n1 = 0;
			n2++;
		}
		x[i] = (0.5 + n1) * sqrt(M_PI);
		y[i] = (0.5 + n2) * sqrt(M_PI);
		n1++;
		//set velocity values
		v_x[i] = velocity_dist(Mersenne_Twister);
		v_y[i] = velocity_dist(Mersenne_Twister);
	}
}

//FORCE
void Force(int norm_time)
{
	//potential energy
	pot_en = 0.0;

	if (norm_time >= norm_thermo_time + norm_free_time)
	{
		if (norm_time % write_period == 0)
		{
			//the number of neighbor particles
			num_neighbor_part.assign(num_part, 0);
			//the orientational order parameter
			oop.assign(num_part, complex<double>(0.0, 0.0));
		}

		//potential microscopic stress tensor
		pot_micro_stress_tensor = 0.0;
	}

	//reset headers in linked-list method
	for (int cell = 0; cell < num_cells_2D; cell++)
		head[cell] = -1;
	//scan particles to construct headers and linked-lists
	for (int i = 0; i < num_part; i++)
	{
		//periodic boundary conditions
		if (x[i] < 0)
		{
			x[i] += box_length;
			if (norm_time >= norm_thermo_time + norm_free_time)
				cross_x[i]--;
		}
		if (x[i] >= box_length)
		{
			x[i] -= box_length;
			if (norm_time >= norm_thermo_time + norm_free_time)
				cross_x[i]++;
		}
		if (y[i] < 0)
		{
			y[i] += box_length;
			if (norm_time >= norm_thermo_time + norm_free_time)
				cross_y[i]--;
		}
		if (y[i] >= box_length)
		{
			y[i] -= box_length;
			if (norm_time >= norm_thermo_time + norm_free_time)
				cross_y[i]++;
		}

		//vector cell index to which the particle belongs
		int cell_x = floor(x[i] / cell_length);
		int cell_y = floor(y[i] / cell_length);
		//translate the vector cell index to the scalar cell index
		int cell = cell_x + cell_y * num_cells_1D;
		//link to the previous occupant (or EMPTY if it's the 1st)
		linked_list[i] = head[cell];
		//the last one goes to the header
		head[cell] = i;

		//Langevin dynamics
		F_x_next[i] = -friction_coef * v_x[i] + random_force * normal_dist(Mersenne_Twister);
		F_y_next[i] = -friction_coef * v_y[i] + random_force * normal_dist(Mersenne_Twister);
	}

	//scan inner cells
	for (int cell_y = 0; cell_y < num_cells_1D; cell_y++)
		for (int cell_x = 0; cell_x < num_cells_1D; cell_x++)
		{
			//calc. the scalar cell index
			int cell = cell_x + cell_y * num_cells_1D;
			//scan the neighbor cells (including itself) of the cell
			for (int neighbor_cell_y = cell_y - 1; neighbor_cell_y <= cell_y + 1; neighbor_cell_y++)
				for (int neighbor_cell_x = cell_x - 1; neighbor_cell_x <= cell_x + 1; neighbor_cell_x++)
				{
					//PBCs by shifting coordinates
					int shift_x;
					int shift_y;
					if (neighbor_cell_x < 0)
						shift_x = -1;
					else if (neighbor_cell_x >= num_cells_1D)
						shift_x = 1;
					else
						shift_x = 0;
					if (neighbor_cell_y < 0)
						shift_y = -1;
					else if (neighbor_cell_y >= num_cells_1D)
						shift_y = 1;
					else
						shift_y = 0;
					//calc. the scalar cell index of the neighbor cell
					int neighbor_cell = (neighbor_cell_x + num_cells_1D) % num_cells_1D
						+ ((neighbor_cell_y + num_cells_1D) % num_cells_1D) * num_cells_1D;
					//scan particle i in the cell
					int i = head[cell];
					while (i != -1)
					{
						//scan atom j in neighbor cell
						int j = head[neighbor_cell];
						while (j != -1)
						{
							//avoid double counting of pair (i, j)
							if (i < j)
							{
								//image-corrected relative pair distance
								double r_x = x[i] - (x[j] + shift_x * box_length);
								double r_y = y[i] - (y[j] + shift_y * box_length);
								double r = sqrt(r_x * r_x + r_y * r_y);
								//calc. within the cutoff radius
								if (r < r_c)
								{
									//Yukawa potential gradient
									double F = exp(-screen_param * r) * (1.0 + screen_param * r) / (2.0 * r * r * r);
									//calc. forces
									F_x_next[i] += F * r_x;
									F_x_next[j] -= F * r_x;
									F_y_next[i] += F * r_y;
									F_y_next[j] -= F * r_y;

									//calc. potential energy
									pot_en += exp(-screen_param * r) / r;

									//-------------------------------------------------------------------------------
									if (norm_time >= norm_thermo_time + norm_free_time)
									{
										//calc. RDF with contribution for particle i and j
										int bin = floor(r / ring_width);
										rdf[bin] += 2;

										//calc. OOP with contribution for particle i and j
										if (norm_time % write_period == 0 && r <= 0.5 * (sqrt(2.0) + 1.0) * sqrt(2.0 * M_PI / sqrt(3.0)))
										{
											double theta = atan(r_y / r_x);
											double Re = cos(6.0 * theta);
											double Im = sin(6.0 * theta);
											oop[i] += complex<double>(Re, Im);
											num_neighbor_part[i]++;
											oop[j] += complex<double>(Re, Im);
											num_neighbor_part[j]++;
										}

										//calc. potential microscopic stress tensor
										pot_micro_stress_tensor += 2.0 * F * r_x * r_y;
									}
									//-------------------------------------------------------------------------------
								}
							}
							j = linked_list[j];
						}
						i = linked_list[i];
					}
				}
		}
	//potential energy per particle
	pot_en /= num_part;
}

//INTEGRATING THE EQ. OF MOTION
void Integration(int switch_on, int norm_time)
{
	switch (switch_on)
	{
	case 1:
		//calc. next time coordinates
		for (int i = 0; i < num_part; i++)
		{
			x[i] += time_step * v_x[i] + time_step * time_step * (4.0 * F_x_curr[i] - F_x_prev[i]) / 6.0;
			y[i] += time_step * v_y[i] + time_step * time_step * (4.0 * F_y_curr[i] - F_y_prev[i]) / 6.0;
		}
		break;
	case 2:
	{
		double sum_v_x = 0.0;
		double sum_v_y = 0.0;
		double sum_v2_x = 0.0;
		double sum_v2_y = 0.0;

		//kinetic microscopic stress tensor
		if (norm_time >= norm_thermo_time + norm_free_time)
			kin_micro_stress_tensor = 0.0;

		for (int i = 0; i < num_part; i++)
		{
			//calc. next time velocities
			v_x[i] += time_step * (2.0 * F_x_next[i] + 5.0 * F_x_curr[i] - F_x_prev[i]) / 6.0;
			v_y[i] += time_step * (2.0 * F_y_next[i] + 5.0 * F_y_curr[i] - F_y_prev[i]) / 6.0;

			if (norm_time < norm_thermo_time)
			{
				sum_v_x += v_x[i];
				sum_v_y += v_y[i];
			}

			sum_v2_x += v_x[i] * v_x[i];
			sum_v2_y += v_y[i] * v_y[i];

			//calc. previous and current time forces
			F_x_prev[i] = F_x_curr[i];
			F_x_curr[i] = F_x_next[i];
			F_y_prev[i] = F_y_curr[i];
			F_y_curr[i] = F_y_next[i];

			//calc. kinetic microscopic stress tensor
			if (norm_time >= norm_thermo_time + norm_free_time)
				kin_micro_stress_tensor += 2.0 * v_x[i] * v_y[i];
		}

		//kinetic energy per particle
		kin_en = (sum_v2_x + sum_v2_y) / num_part;

		//thermostat
		if (norm_time < norm_thermo_time)
		{
			//velocity of the center of mass
			sum_v_x /= num_part;
			sum_v_y /= num_part;
			//scale factor of velocities
			double scale_factor_x = sqrt(0.5 * num_part / (coupl_param * sum_v2_x));
			double scale_factor_y = sqrt(0.5 * num_part / (coupl_param * sum_v2_y));
			//rescale velocities
			for (int i = 0; i < num_part; i++)
			{
				v_x[i] = (v_x[i] - sum_v_x) * scale_factor_x;
				v_y[i] = (v_y[i] - sum_v_y) * scale_factor_y;
			}
		}
		break;
	}
	default:
		break;
	}
}

//RADIAL DISTRIBUTION FUNCTION
void Radial_Distribution_Function(int switch_on)
{
	switch (switch_on)
	{
	case 0:
	{
		//number of samples
		rdf_samples = 0;
		//number of bins
		int num_bins = floor(10.0 * r_c);
		//bin size
		ring_width = r_c / num_bins;
		//Radial Distribution Function (RDF)
		rdf.assign(num_bins, 0.0);
		break;
	}
	case 1:
		rdf_samples++;
		break;
	case 2:
	{
		//constant multiplier
		double multiplier = num_part * rdf_samples * ring_width * ring_width;
		ofstream rdf_file(title_of_rdf, ios::app);
		for (int bin = 0; bin < rdf.size(); bin++)
		{
			//RDF averaged over number of particles and time samples
			rdf[bin] /= (multiplier * (2.0 * bin + 1));
			//write RDF to the file
			rdf_file << (bin + 0.5) * ring_width << '\t' << rdf[bin] << '\n';
		}
		rdf_file.close();
		break;
	}
	default:
		break;
	}
}

//ORIENTATIONAL ORDER PARAMETER
void Orientational_Order_Parameter(int norm_time)
{
	double local_oop = 0.0;
	complex<double> overall_oop(0.0, 0.0);
	for (int i = 0; i < num_part; i++)
	{
		local_oop += abs(oop[i] / static_cast<double>(num_neighbor_part[i]));
		overall_oop += oop[i] / static_cast<double>(num_neighbor_part[i]);
	}
	//local OOP of the system
	local_oop /= num_part;
	//overall OOP of the system
	overall_oop /= static_cast<double>(num_part);
	//write OOP to the file
	ofstream oop_file(title_of_oop, ios::app);
	oop_file << (norm_time - (norm_thermo_time + norm_free_time)) * time_step << '\t' << local_oop << '\t' << abs(overall_oop) << '\n';
	oop_file.close();
}

//DISPERSION
void Dispersion(int switch_on)
{
	switch (switch_on)
	{
	case 0:
		//microscopic longitudinal and transverse currents
		micro_long_current.reserve(norm_meas_time);
		micro_tran_current.reserve(norm_meas_time);
		break;
	case 1:
	{
		//temporary arrays
		vector<complex<double>> temp_micro_long_current(norm_max_wave_num, complex<double>(0.0, 0.0));
		vector<complex<double>> temp_micro_tran_current(norm_max_wave_num, complex<double>(0.0, 0.0));
		for (int norm_wave_num = 0; norm_wave_num < norm_max_wave_num; norm_wave_num++)
		{
			for (int i = 0; i < num_part; i++)
			{
				double Re = cos(norm_wave_num * wave_num_step * x[i]);
				double Im = sin(norm_wave_num * wave_num_step * x[i]);
				temp_micro_long_current[norm_wave_num] += v_x[i] * complex<double>(Re, Im);
				temp_micro_tran_current[norm_wave_num] += v_y[i] * complex<double>(Re, Im);
			}
			temp_micro_long_current[norm_wave_num] *= norm_wave_num * wave_num_step;
			temp_micro_tran_current[norm_wave_num] *= norm_wave_num * wave_num_step;
		}
		micro_long_current.emplace_back(temp_micro_long_current);
		micro_tran_current.emplace_back(temp_micro_tran_current);
		break;
	}
	case 2:
		for (int norm_wave_num = 0; norm_wave_num < norm_max_wave_num; norm_wave_num++)
			for (int norm_disp_freq = 0; norm_disp_freq < norm_max_disp_freq; norm_disp_freq++)
			{
				//Fourier Transform of the microscopic longitudinal and transverse currents
				complex<double> FT_micro_long_current(0.0, 0.0);
				complex<double> FT_micro_tran_current(0.0, 0.0);
				//integrate by the Simpson's rule
				for (int norm_time = 0; norm_time < norm_meas_time; norm_time++)
				{
					//Hann window function
					double window_func = 0.5 * (1.0 - cos(2.0 * M_PI * norm_time / (norm_meas_time - 1)));
					//real and imaginary parts of FT
					double Re = cos(norm_disp_freq * disp_freq_step * norm_time * time_step);
					double Im = -sin(norm_disp_freq * disp_freq_step * norm_time * time_step);
					//integrands
					complex<double> integrand_micro_long_current = micro_long_current[norm_time][norm_wave_num] * complex<double>(Re, Im) * window_func;
					complex<double> integrand_micro_tran_current = micro_tran_current[norm_time][norm_wave_num] * complex<double>(Re, Im) * window_func;
					if (norm_time == 0 || norm_time == (norm_meas_time - 1))
					{
						FT_micro_long_current += integrand_micro_long_current;
						FT_micro_tran_current += integrand_micro_tran_current;
					}
					else if (norm_time % 2 != 0)
					{
						FT_micro_long_current += 4.0 * integrand_micro_long_current;
						FT_micro_tran_current += 4.0 * integrand_micro_tran_current;
					}
					else
					{
						FT_micro_long_current += 2.0 * integrand_micro_long_current;
						FT_micro_tran_current += 2.0 * integrand_micro_tran_current;
					}
				}
				//fluctuation spectra
				double long_spectrum = abs(FT_micro_long_current) * abs(FT_micro_long_current) / (6.0 * M_PI * num_part * norm_meas_time);
				double tran_spectrum = abs(FT_micro_tran_current) * abs(FT_micro_tran_current) / (6.0 * M_PI * num_part * norm_meas_time);
				//write dispersion to the file
				ofstream dispersion_file(title_of_dispersion, ios::app);
				dispersion_file << norm_wave_num * wave_num_step << '\t' << norm_disp_freq * disp_freq_step << '\t' << log10(long_spectrum) << '\t' << log10(tran_spectrum) << '\n';
				dispersion_file.close();
			}
		break;
	default:
		break;
	}
}

//MEAN-SQUARED DISPLACEMENT, VELOCITY AND STRESS AUTOCORRELATION FUNCTIONS
void MSD_VACF_SACF(int switch_on)
{
	switch (switch_on)
	{
	case 0:
	{
		//time counter
		time_counter = 0;
		//number of blocks
		num_blocks = 4;
		//number of block elements
		num_block_elem = 100;
		//number of points in one unit of "log10" scale
		num_points_log10 = 10;
		//length of a block
		block_length.assign(num_blocks, 0);
		//temporary arrays
		vector<int> temp_int_v(num_block_elem, 0);
		vector<double> temp_double_v(num_block_elem, 0.0);
		vector<vector<double>> temp_double_v2(num_block_elem, vector<double>(num_part, 0.0));
		//number of samples
		msd_vacf_safc_samples.assign(num_blocks, temp_int_v);

		//Mean-Squared Displament (MSD)
		msd.assign(num_blocks, temp_double_v);
		//coordinates at initial time
		x_time_0.assign(num_blocks, temp_double_v2);
		y_time_0.assign(num_blocks, temp_double_v2);
		//number of boundary crossings in each direction
		cross_x.assign(num_part, 0);
		cross_y.assign(num_part, 0);

		//Velocity Autocorrealation Function (VACF)
		vacf.assign(num_blocks, temp_double_v);
		//velocities at initial time
		v_x_time_0.assign(num_blocks, temp_double_v2);
		v_y_time_0.assign(num_blocks, temp_double_v2);

		//Stress Autocorrealation Function (SACF)
		sacf.assign(num_blocks, temp_double_v);
		//microscopic stress tensor at initial time
		micro_stress_tensor_0.assign(num_blocks, temp_double_v);
		break;
	}
	case 1:
		//loop over all the blocks to test which blocks need sampling
		for (int block = 0; block < num_blocks; block++)
			//test for blocking operation, i.e. when "time_counter" is a multiple of "num_block_elem^block"
			if (time_counter % static_cast<int>(pow(num_block_elem / num_points_log10, block)) == 0)
			{
				//increase the current block-length
				block_length[block]++;
				//compute the current length of the block, limited to size "num_block_elem"
				int curr_block_length = min(block_length[block], num_block_elem);
				//loop over the block elements
				for (int block_elem = 0; block_elem < curr_block_length; block_elem++)
				{
					//loop over the particles
					for (int i = 0; i < num_part; i++)
					{
						//set last index to the correlation value
						if (block_elem == curr_block_length - 1)
						{
							//store positions of initial time, take into account the boundary crossings
							x_time_0[block][curr_block_length - 1][i] = x[i] + cross_x[i] * box_length;
							y_time_0[block][curr_block_length - 1][i] = y[i] + cross_y[i] * box_length;

							//store velocities of initial time
							v_x_time_0[block][curr_block_length - 1][i] = v_x[i];
							v_y_time_0[block][curr_block_length - 1][i] = v_y[i];
						}
						//shift to the right if the block is full
						else if (block_length[block] > num_block_elem)
						{
							x_time_0[block][block_elem][i] = x_time_0[block][block_elem + 1][i];
							y_time_0[block][block_elem][i] = y_time_0[block][block_elem + 1][i];

							v_x_time_0[block][block_elem][i] = v_x_time_0[block][block_elem + 1][i];
							v_y_time_0[block][block_elem][i] = v_y_time_0[block][block_elem + 1][i];
						}

						//calc. MSD
						double dx = (x[i] + cross_x[i] * box_length) - x_time_0[block][block_elem][i];
						double dy = (y[i] + cross_y[i] * box_length) - y_time_0[block][block_elem][i];
						msd[block][curr_block_length - 1 - block_elem] += dx * dx + dy * dy;

						//calc. VACF
						vacf[block][curr_block_length - 1 - block_elem] += v_x[i] * v_x_time_0[block][block_elem][i]
							+ v_y[i] * v_y_time_0[block][block_elem][i];
					}

					//set last index to the correlation value
					if (block_elem == curr_block_length - 1)
						//store values of microscopic stress tensor for initial time
						micro_stress_tensor_0[block][curr_block_length - 1] = kin_micro_stress_tensor + pot_micro_stress_tensor;
					//shift to the right if the block is full
					else if (block_length[block] > num_block_elem)
						micro_stress_tensor_0[block][block_elem] = micro_stress_tensor_0[block][block_elem + 1];
					//calc. SACF
					sacf[block][curr_block_length - 1 - block_elem] += (kin_micro_stress_tensor + pot_micro_stress_tensor)
						* micro_stress_tensor_0[block][block_elem];

					//count the number of sampling
					msd_vacf_safc_samples[block][curr_block_length - 1 - block_elem]++;
				}
			}
		//count the current sampling
		time_counter++;
		break;
	case 2:
	{
		ofstream msd_file(title_of_msd, ios::app);
		ofstream vacf_file(title_of_vacf, ios::app);
		ofstream sacf_file(title_of_sacf, ios::app);
		for (int block = 0; block < num_blocks; block++)
		{
			int start_block_elem;
			if (block > 0)
				start_block_elem = num_points_log10;
			else
				start_block_elem = 0;
			int curr_block_length = min(block_length[block], num_block_elem);
			for (int block_elem = start_block_elem; block_elem < curr_block_length; block_elem++)
			{
				//write MSD to the file
				msd_file << block_elem * time_step * pow(num_block_elem / num_points_log10, block) << '\t'
					<< msd[block][block_elem] / (num_part * msd_vacf_safc_samples[block][block_elem]) << '\n';

				//write VACF to the file
				vacf_file << block_elem * time_step * pow(num_block_elem / num_points_log10, block) << '\t'
					<< vacf[block][block_elem] / (num_part * msd_vacf_safc_samples[block][block_elem]) << '\n';

				//write SACF to the file
				sacf_file << block_elem * time_step * pow(num_block_elem / num_points_log10, block) << '\t'
					<< sacf[block][block_elem] * coupl_param / (2.0 * num_part * msd_vacf_safc_samples[block][block_elem]) << '\n';
			}
		}
		msd_file.close();
		vacf_file.close();
		sacf_file.close();
		break;
	}
	default:
		break;
	}
}

//DISPERSION
//void Dispersion(int switch_on)
//{
//	switch (switch_on)
//	{
//	case 1:
//		//clear the temporary arrays
//		temp_lambda_x.clear();
//		temp_tau_x.clear();
//
//		temp_lambda_y.clear();
//		temp_tau_y.clear();
//
//		for (double k = dk; k <= k_max; k += dk)
//		{
//			int Nk = round(k / dk) - 1;
//			temp_lambda_x.push_back(complex<double>(0.0, 0.0));
//			temp_tau_x.push_back(complex<double>(0.0, 0.0));
//
//			temp_lambda_y.push_back(complex<double>(0.0, 0.0));
//			temp_tau_y.push_back(complex<double>(0.0, 0.0));
//
//			for (int i = 0; i < N; i++)
//			{
//				double Re_x = cos(k * x[i]);
//				double Im_x = sin(k * x[i]);
//				temp_lambda_x[Nk] += v_x[i] * complex<double>(Re_x, Im_x);
//				temp_tau_x[Nk] += v_y[i] * complex<double>(Re_x, Im_x);
//
//				double Re_y = cos(k * y[i]);
//				double Im_y = sin(k * y[i]);
//				temp_lambda_y[Nk] += v_x[i] * complex<double>(Re_y, Im_y);
//				temp_tau_y[Nk] += v_y[i] * complex<double>(Re_y, Im_y);
//			}
//			temp_lambda_x[Nk] *= k;
//			temp_tau_x[Nk] *= k;
//
//			temp_lambda_y[Nk] *= k;
//			temp_tau_y[Nk] *= k;
//		}
//		//longitudinal & transverse currents
//		lambda_x.push_back(temp_lambda_x);
//		tau_x.push_back(temp_tau_x);
//
//		lambda_y.push_back(temp_lambda_y);
//		tau_y.push_back(temp_tau_y);
//
//		break;
//	case 2:
//		for (double k = dk; k <= k_max; k += dk)
//		{
//			int Nk = round(k / dk) - 1;
//			for (double omega = domega; omega <= omega_max; omega += domega)
//			{
//				//Fourier Transform of longitudinal & transverse currents
//				complex<double> FT_lambda_x(0.0, 0.0);
//				complex<double> FT_tau_x(0.0, 0.0);
//
//				complex<double> FT_lambda_y(0.0, 0.0);
//				complex<double> FT_tau_y(0.0, 0.0);
//
//				for (int Nt = 0; Nt < Nt_meas; Nt++)
//				{
//					//Hann window function
//					double window_func = 0.5 * (1.0 - cos(2.0 * M_PI * Nt / (Nt_meas - 1)));
//					//real and imaginary parts of FT
//					double Re = cos(omega * Nt * dt);
//					double Im = -sin(omega * Nt * dt);
//					//integrands
//					complex<double> integrand_lambda_x = lambda_x[Nt][Nk] * window_func * complex<double>(Re, Im);
//					complex<double> integrand_tau_x = tau_x[Nt][Nk] * window_func * complex<double>(Re, Im);
//
//					complex<double> integrand_lambda_y = lambda_y[Nt][Nk] * window_func * complex<double>(Re, Im);
//					complex<double> integrand_tau_y = tau_y[Nt][Nk] * window_func * complex<double>(Re, Im);
//
//					//integrate by the Simpson's rule
//					if (Nt == 0 || Nt == (Nt_meas - 1))
//					{
//						FT_lambda_x += integrand_lambda_x;
//						FT_tau_x += integrand_tau_x;
//
//						FT_lambda_y += integrand_lambda_y;
//						FT_tau_y += integrand_tau_y;
//					}
//					else if (Nt % 2 != 0)
//					{
//						FT_lambda_x += 4.0 * integrand_lambda_x;
//						FT_tau_x += 4.0 * integrand_tau_x;
//
//						FT_lambda_y += 4.0 * integrand_lambda_y;
//						FT_tau_y += 4.0 * integrand_tau_y;
//					}
//					else
//					{
//						FT_lambda_x += 2.0 * integrand_lambda_x;
//						FT_tau_x += 2.0 * integrand_tau_x;
//
//						FT_lambda_y += 2.0 * integrand_lambda_y;
//						FT_tau_y += 2.0 * integrand_tau_y;
//					}
//				}
//				//fluctuation spectra
//				double Longitudinal_x = abs(FT_lambda_x) * abs(FT_lambda_x) / (6.0 * M_PI * N * Nt_meas);
//				double Transverse_x = abs(FT_tau_x) * abs(FT_tau_x) / (6.0 * M_PI * N * Nt_meas);
//
//				double Longitudinal_y = abs(FT_lambda_y) * abs(FT_lambda_y) / (6.0 * M_PI * N * Nt_meas);
//				double Transverse_y = abs(FT_tau_y) * abs(FT_tau_y) / (6.0 * M_PI * N * Nt_meas);
//
//				//write dispersion to the file
//				ofstream dispersion(title_of_dispersion, ios::app);
//				dispersion << k << '\t' << omega << '\t' << log10(Longitudinal_x) << '\t' << log10(Transverse_x)
//					<< '\t' << log10(Longitudinal_y) << '\t' << log10(Transverse_y) << '\n';
//				dispersion.close();
//			}
//		}
//		break;
//	default:
//		break;
//	}
//}