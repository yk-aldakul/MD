//10.05.2021
//3D Yukawa system

#define _USE_MATH_DEFINES
#include <iostream>
#include <fstream>
#include <cmath>
#include <complex>
#include <random>
#include <string>
#include <vector>
#include <chrono>
//#include <windows.h>

//set parameters
int num_part = 103823;
int norm_thermo_time = 100000;
int norm_free_time = 10000;
int norm_meas_time = 100000;
int write_period = 100;
double time_step = 0.01;
double fric_coef = 0.0;
double coupl_param = 200.0;
double screen_param = 2.0;

double max_wave_num = 3.0;
double max_disp_freq = 1.5;

//file names
std::string title_of_coordinates = "XYZ ().xyz";
std::string title_of_energy = "ENERGY ().dat";
std::string title_of_rdf = "RDF ().dat";
std::string title_of_dispersion = "DISPERSION ().dat";
std::string title_of_msd = "MSD ().dat";
std::string title_of_vacf = "VACF ().dat";
std::string title_of_sacf = "SACF ().dat";

//random number engine
std::random_device rand_dev;
std::mt19937_64 Mersenne_Twister(rand_dev());

//system parameters
double box_length = std::cbrt(4.0 * M_PI * static_cast<double>(num_part) / 3.0);
double Cutoff_Radius()
{
	double r = 1.0, accuracy = std::pow(10.0, -6.0);
	while (std::exp(-screen_param * r) / r >= accuracy)
		r += 0.01;
	std::cout << "r_c = " << r << '\n';
	return r;
}
double r_c = Cutoff_Radius();
int num_cells_1D = std::floor(box_length / r_c);
double cell_length = box_length / static_cast<double>(num_cells_1D);
int num_cells_2D = num_cells_1D * num_cells_1D;
int num_cells_3D = num_cells_2D * num_cells_1D;
double rand_force = std::sqrt(2.0 * fric_coef / (3.0 * coupl_param * time_step));
std::normal_distribution<double> norm_dist(0.0, 1.0);

double wave_num_step = 2.0 * M_PI / box_length;
double disp_freq_step = 2.0 * M_PI / (static_cast<double>(norm_meas_time) * time_step);
int norm_max_wave_num = std::round(max_wave_num / wave_num_step);
int norm_max_disp_freq = std::round(max_disp_freq / disp_freq_step);

//global variables
//coordinates & velocity components
std::vector<double> x(num_part), y(num_part), z(num_part);
std::vector<double> v_x(num_part), v_y(num_part), v_z(num_part);

//previous, current, & next time force components
std::vector<double> F_prev_x(num_part), F_prev_y(num_part), F_prev_z(num_part);
std::vector<double> F_curr_x(num_part), F_curr_y(num_part), F_curr_z(num_part);
std::vector<double> F_next_x(num_part), F_next_y(num_part), F_next_z(num_part);

//potential & kinetic energy functions
double pot_en;
double kin_en;

//linked list variables
std::vector<int> linked_list(num_part);
std::vector<int> head(num_cells_3D);

//RDF variables
int rdf_samples;
double ring_width;
std::vector<double> rdf;

//Dispersion variables
std::vector<std::vector<std::complex<double>>> long_mc;
std::vector<std::vector<std::complex<double>>> tran_mc;

//MSD, VACF & SACF variables
int time_counter;
int num_blocks;
int num_block_elem;
int num_points_log10;
std::vector<int> block_length;
std::vector<std::vector<int>> msd_vacf_safc_samples;

std::vector<double> cross_x, cross_y, cross_z;
std::vector<std::vector<std::vector<double>>> x_time_0, y_time_0, z_time_0;
std::vector<std::vector<double>> msd;

std::vector<std::vector<std::vector<double>>> v_x_time_0, v_y_time_0, v_z_time_0;
std::vector<std::vector<double>> vacf;

double kin_mst_xy, kin_mst_yz, kin_mst_zx;
double pot_mst_xy, pot_mst_yz, pot_mst_zx;
std::vector<std::vector<double>> mst_xy_time_0, mst_yz_time_0, mst_zx_time_0;
std::vector<std::vector<double>> sacf_xy, sacf_yz, sacf_zx;

//INITIALIZATION
void Initialization()
{
	//place particles in a 3D grid
	int n1 = 0;
	int n2 = 0;
	int n3 = 0;
	//edge limit
	int edge = std::cbrt(num_part);
	double separation = std::cbrt(4.0 * M_PI / 3.0);
	double stand_dev = 1.0 / std::sqrt(3.0 * coupl_param);
	std::normal_distribution<double> velocity_dist(0.0, stand_dev);
	for (int i = 0; i < num_part; i++)
	{
		if (n1 >= edge)
		{
			n1 = 0;
			n2++;
		}
		if (n2 >= edge)
		{
			n2 = 0;
			n3++;
		}
		if (n3 >= edge)
			std::cout << "Looks like you are out of the main box!" << '\n';
		x[i] = (0.5 + static_cast<double>(n1)) * separation;
		y[i] = (0.5 + static_cast<double>(n2)) * separation;
		z[i] = (0.5 + static_cast<double>(n3)) * separation;
		n1++;
		//set velocity values
		v_x[i] = velocity_dist(Mersenne_Twister);
		v_y[i] = velocity_dist(Mersenne_Twister);
		v_z[i] = velocity_dist(Mersenne_Twister);
	}

	//read MD data from a file
	/*int i;
	double x_i;
	double y_i;
	double z_i;
	double v_x_i;
	double v_y_i;
	double v_z_i;
	std::ifstream md_data("");
	while (md_data >> i >> x_i >> y_i >> z_i >> v_x_i >> v_y_i >> v_z_i)
	{
		x[i] = x_i;
		y[i] = y_i;
		z[i] = z_i;
		v_x[i] = v_x_i;
		v_y[i] = v_y_i;
		v_z[i] = v_z_i;
	}
	md_data.close();*/
}

//FORCE
void Force(int& norm_time)
{
	//potential energy
	pot_en = 0.0;

	//potential microscopic stress tensor
	if (norm_time >= norm_thermo_time + norm_free_time)
	{
		pot_mst_xy = 0.0;
		pot_mst_yz = 0.0;
		pot_mst_zx = 0.0;
	}

	//reset headers in linked-list method
	for (int cell = 0; cell < num_cells_3D; cell++)
		head[cell] = -1;
	//scan particles to construct headers & linked-lists
	for (int i = 0; i < num_part; i++)
	{
		//periodic boundary conditions
		if (x[i] < 0)
		{
			x[i] += box_length;
			if (norm_time >= norm_thermo_time + norm_free_time)
				cross_x[i] -= box_length;
		}
		if (x[i] >= box_length)
		{
			x[i] -= box_length;
			if (norm_time >= norm_thermo_time + norm_free_time)
				cross_x[i] += box_length;
		}
		if (y[i] < 0)
		{
			y[i] += box_length;
			if (norm_time >= norm_thermo_time + norm_free_time)
				cross_y[i] -= box_length;
		}
		if (y[i] >= box_length)
		{
			y[i] -= box_length;
			if (norm_time >= norm_thermo_time + norm_free_time)
				cross_y[i] += box_length;
		}
		if (z[i] < 0)
		{
			z[i] += box_length;
			if (norm_time >= norm_thermo_time + norm_free_time)
				cross_z[i] -= box_length;
		}
		if (z[i] >= box_length)
		{
			z[i] -= box_length;
			if (norm_time >= norm_thermo_time + norm_free_time)
				cross_z[i] += box_length;
		}

		//vector cell index to which the particle belongs
		int cell_x = std::floor(x[i] / cell_length);
		int cell_y = std::floor(y[i] / cell_length);
		int cell_z = std::floor(z[i] / cell_length);
		//translate the vector cell index to the scalar cell index
		int cell = cell_x + cell_y * num_cells_1D + cell_z * num_cells_2D;
		//link to the previous occupant (or EMPTY if it's the 1st)
		linked_list[i] = head[cell];
		//the last one goes to the header
		head[cell] = i;

		//Langevin dynamics
		F_next_x[i] = -fric_coef * v_x[i] + rand_force * norm_dist(Mersenne_Twister);
		F_next_y[i] = -fric_coef * v_y[i] + rand_force * norm_dist(Mersenne_Twister);
		F_next_z[i] = -fric_coef * v_z[i] + rand_force * norm_dist(Mersenne_Twister);
	}

	//scan inner cells
	for (int cell_z = 0; cell_z < num_cells_1D; cell_z++)
		for (int cell_y = 0; cell_y < num_cells_1D; cell_y++)
			for (int cell_x = 0; cell_x < num_cells_1D; cell_x++)
			{
				//calc. the scalar cell index
				int cell = cell_x + cell_y * num_cells_1D + cell_z * num_cells_2D;
				//scan the neighbor cells (including itself) of the cell
				for (int neighbor_cell_z = cell_z - 1; neighbor_cell_z <= cell_z + 1; neighbor_cell_z++)
					for (int neighbor_cell_y = cell_y - 1; neighbor_cell_y <= cell_y + 1; neighbor_cell_y++)
						for (int neighbor_cell_x = cell_x - 1; neighbor_cell_x <= cell_x + 1; neighbor_cell_x++)
						{
							//PBCs by shifting coordinates
							double shift_x, shift_y, shift_z;
							if (neighbor_cell_x < 0)
								shift_x = -box_length;
							else if (neighbor_cell_x >= num_cells_1D)
								shift_x = box_length;
							else
								shift_x = 0;
							if (neighbor_cell_y < 0)
								shift_y = -box_length;
							else if (neighbor_cell_y >= num_cells_1D)
								shift_y = box_length;
							else
								shift_y = 0;
							if (neighbor_cell_z < 0)
								shift_z = -box_length;
							else if (neighbor_cell_z >= num_cells_1D)
								shift_z = box_length;
							else
								shift_z = 0;
							//calc. the scalar cell index of the neighbor cell
							int neighbor_cell = (neighbor_cell_x + num_cells_1D) % num_cells_1D
								+ ((neighbor_cell_y + num_cells_1D) % num_cells_1D) * num_cells_1D
								+ ((neighbor_cell_z + num_cells_1D) % num_cells_1D) * num_cells_2D;
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
										double r_x = x[i] - (x[j] + shift_x);
										double r_y = y[i] - (y[j] + shift_y);
										double r_z = z[i] - (z[j] + shift_z);
										double r = std::hypot(r_x, r_y, r_z);
										//calc. within the cutoff radius
										if (r < r_c)
										{
											//Yukawa potential gradient
											double F = std::exp(-screen_param * r) * (1.0 + screen_param * r) / (3.0 * r * r * r);
											//calc. forces
											F_next_x[i] += F * r_x;
											F_next_x[j] -= F * r_x;
											F_next_y[i] += F * r_y;
											F_next_y[j] -= F * r_y;
											F_next_z[i] += F * r_z;
											F_next_z[j] -= F * r_z;

											//calc. potential energy
											pot_en += std::exp(-screen_param * r) / r;

											//-----------------------------------------------------------------------------------
											if (norm_time >= norm_thermo_time + norm_free_time)
											{
												//calc. RDF with contribution for particle i & j
												int bin = std::floor(r / ring_width);
												rdf[bin] += 2;

												//calc. potential microscopic stress tensor
												pot_mst_xy += 3.0 * F * r_x * r_y;
												pot_mst_yz += 3.0 * F * r_y * r_z;
												pot_mst_zx += 3.0 * F * r_z * r_x;
											}
											//-----------------------------------------------------------------------------------
										}
									}
									j = linked_list[j];
								}
								i = linked_list[i];
							}
						}
			}
	//potential energy per particle
	pot_en /= static_cast<double>(num_part);
}

//INTEGRATING THE EQ. OF MOTION
void Integration(int switch_on, int& norm_time)
{
	switch (switch_on)
	{
	case 1:
		//calc. next time coordinates
		for (int i = 0; i < num_part; i++)
		{
			x[i] += time_step * v_x[i] + time_step * time_step * (4.0 * F_curr_x[i] - F_prev_x[i]) / 6.0;
			y[i] += time_step * v_y[i] + time_step * time_step * (4.0 * F_curr_y[i] - F_prev_y[i]) / 6.0;
			z[i] += time_step * v_z[i] + time_step * time_step * (4.0 * F_curr_z[i] - F_prev_z[i]) / 6.0;
		}
		break;
	case 2:
	{
		double sum_v_x = 0.0;
		double sum_v_y = 0.0;
		double sum_v_z = 0.0;
		double sum_v2_x = 0.0;
		double sum_v2_y = 0.0;
		double sum_v2_z = 0.0;

		//kinetic microscopic stress tensor
		if (norm_time >= norm_thermo_time + norm_free_time)
		{
			kin_mst_xy = 0.0;
			kin_mst_yz = 0.0;
			kin_mst_zx = 0.0;
		}

		for (int i = 0; i < num_part; i++)
		{
			//calc. next time velocities
			v_x[i] += time_step * (2.0 * F_next_x[i] + 5.0 * F_curr_x[i] - F_prev_x[i]) / 6.0;
			v_y[i] += time_step * (2.0 * F_next_y[i] + 5.0 * F_curr_y[i] - F_prev_y[i]) / 6.0;
			v_z[i] += time_step * (2.0 * F_next_z[i] + 5.0 * F_curr_z[i] - F_prev_z[i]) / 6.0;

			if (norm_time < norm_thermo_time)
			{
				sum_v_x += v_x[i];
				sum_v_y += v_y[i];
				sum_v_z += v_z[i];
			}

			sum_v2_x += v_x[i] * v_x[i];
			sum_v2_y += v_y[i] * v_y[i];
			sum_v2_z += v_z[i] * v_z[i];

			//calc. previous & current time forces
			F_prev_x[i] = F_curr_x[i];
			F_curr_x[i] = F_next_x[i];
			F_prev_y[i] = F_curr_y[i];
			F_curr_y[i] = F_next_y[i];
			F_prev_z[i] = F_curr_z[i];
			F_curr_z[i] = F_next_z[i];

			//calc. kinetic microscopic stress tensor
			if (norm_time >= norm_thermo_time + norm_free_time)
			{
				kin_mst_xy += 3.0 * v_x[i] * v_y[i];
				kin_mst_yz += 3.0 * v_y[i] * v_z[i];
				kin_mst_zx += 3.0 * v_z[i] * v_x[i];
			}
		}

		//kinetic energy per particle
		kin_en = 1.5 * (sum_v2_x + sum_v2_y + sum_v2_z) / static_cast<double>(num_part);

		//thermostat
		if (norm_time < norm_thermo_time)
		{
			//velocity of the center of mass
			sum_v_x /= static_cast<double>(num_part);
			sum_v_y /= static_cast<double>(num_part);
			sum_v_z /= static_cast<double>(num_part);
			//scale factor of velocities
			double scale_factor_x = 1.0 / std::sqrt(3.0 * coupl_param * sum_v2_x / static_cast<double>(num_part));
			double scale_factor_y = 1.0 / std::sqrt(3.0 * coupl_param * sum_v2_y / static_cast<double>(num_part));
			double scale_factor_z = 1.0 / std::sqrt(3.0 * coupl_param * sum_v2_z / static_cast<double>(num_part));
			//rescale velocities
			for (int i = 0; i < num_part; i++)
			{
				v_x[i] = (v_x[i] - sum_v_x) * scale_factor_x;
				v_y[i] = (v_y[i] - sum_v_y) * scale_factor_y;
				v_z[i] = (v_z[i] - sum_v_z) * scale_factor_z;
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
		int num_bins = std::floor(10.0 * r_c);
		//bin size
		ring_width = r_c / static_cast<double>(num_bins);
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
		double multiplier = static_cast<double>(num_part) * static_cast<double>(rdf_samples) * std::pow(ring_width, 3.0);
		std::ofstream rdf_file(title_of_rdf, std::ios::app);
		for (int bin = 0; bin < rdf.size(); bin++)
		{
			//RDF averaged over number of particles & time samples
			rdf[bin] /= (multiplier * (3.0 * static_cast<double>(bin) * (static_cast<double>(bin) + 1) + 1));
			//write RDF to the file
			rdf_file << (static_cast<double>(bin) + 0.5) * ring_width << '\t' << rdf[bin] << '\n';
		}
		rdf_file.close();
		break;
	}
	default:
		break;
	}
}

//DISPERSION
void Dispersion(int switch_on)
{
	switch (switch_on)
	{
	case 0:
		//microscopic longitudinal & transverse currents
		long_mc.reserve(norm_meas_time);
		tran_mc.reserve(norm_meas_time);
		break;
	case 1:
	{
		//temporary arrays
		std::vector<std::complex<double>> temp_long_mc(norm_max_wave_num, std::complex<double>(0.0, 0.0));
		std::vector<std::complex<double>> temp_tran_mc(norm_max_wave_num, std::complex<double>(0.0, 0.0));
		for (int norm_wave_num = 0; norm_wave_num < norm_max_wave_num; norm_wave_num++)
		{
			for (int i = 0; i < num_part; i++)
			{
				double Re = cos(static_cast<double>(norm_wave_num) * wave_num_step * x[i]);
				double Im = sin(static_cast<double>(norm_wave_num) * wave_num_step * x[i]);
				temp_long_mc[norm_wave_num] += v_x[i] * std::complex<double>(Re, Im);
				temp_tran_mc[norm_wave_num] += v_y[i] * std::complex<double>(Re, Im);
			}
			temp_long_mc[norm_wave_num] *= static_cast<double>(norm_wave_num) * wave_num_step;
			temp_tran_mc[norm_wave_num] *= static_cast<double>(norm_wave_num) * wave_num_step;
		}
		//microscopic longitudinal & transverse currents
		long_mc.emplace_back(temp_long_mc);
		tran_mc.emplace_back(temp_tran_mc);
		break;
	}
	case 2:
		for (int norm_wave_num = 0; norm_wave_num < norm_max_wave_num; norm_wave_num++)
			for (int norm_disp_freq = 0; norm_disp_freq < norm_max_disp_freq; norm_disp_freq++)
			{
				//Fourier Transform of the microscopic longitudinal & transverse currents
				std::complex<double> FT_long_mc(0.0, 0.0);
				std::complex<double> FT_tran_mc(0.0, 0.0);
				//integrate by the Simpson's rule
				for (int norm_time = 0; norm_time < norm_meas_time; norm_time++)
				{
					//Hann window function
					double window_func = 0.5 * (1.0 - cos(2.0 * M_PI * norm_time / (static_cast<double>(norm_meas_time) - 1.0)));
					//real & imaginary parts of FT
					double Re = cos(static_cast<double>(norm_disp_freq) * disp_freq_step * static_cast<double>(norm_time) * time_step);
					double Im = -sin(static_cast<double>(norm_disp_freq) * disp_freq_step * static_cast<double>(norm_time) * time_step);
					//integrands
					std::complex<double> integ_long_mc = long_mc[norm_time][norm_wave_num] * std::complex<double>(Re, Im) * window_func;
					std::complex<double> integ_tran_mc = tran_mc[norm_time][norm_wave_num] * std::complex<double>(Re, Im) * window_func;
					if (norm_time == 0 || norm_time == (norm_meas_time - 1))
					{
						FT_long_mc += integ_long_mc;
						FT_tran_mc += integ_tran_mc;
					}
					else if (norm_time % 2 != 0)
					{
						FT_long_mc += 4.0 * integ_long_mc;
						FT_tran_mc += 4.0 * integ_tran_mc;
					}
					else
					{
						FT_long_mc += 2.0 * integ_long_mc;
						FT_tran_mc += 2.0 * integ_tran_mc;
					}
				}
				//fluctuation spectra
				double long_spectrum = std::abs(FT_long_mc) * std::abs(FT_long_mc) / (6.0 * M_PI * num_part * norm_meas_time);
				double tran_spectrum = std::abs(FT_tran_mc) * std::abs(FT_tran_mc) / (6.0 * M_PI * num_part * norm_meas_time);
				//write dispersion to the file
				std::ofstream dispersion_file(title_of_dispersion, std::ios::app);
				dispersion_file << static_cast<double>(norm_wave_num) * wave_num_step << '\t'
					<< static_cast<double>(norm_disp_freq) * disp_freq_step << '\t'
					<< std::log10(long_spectrum) << '\t' << std::log10(tran_spectrum) << '\n';
				dispersion_file.close();
			}
		break;
	default:
		break;
	}
}

//MEAN-SQUARED DISPLACEMENT, VELOCITY & STRESS AUTOCORRELATION FUNCTIONS
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
		std::vector<int> temp_int_v(num_block_elem, 0);
		std::vector<double> temp_double_v(num_block_elem, 0.0);
		std::vector<std::vector<double>> temp_double_v2(num_block_elem, std::vector<double>(num_part, 0.0));
		//number of samples
		msd_vacf_safc_samples.assign(num_blocks, temp_int_v);

		//Mean-Squared Displament (MSD)
		msd.assign(num_blocks, temp_double_v);
		//coordinates at initial time
		x_time_0.assign(num_blocks, temp_double_v2);
		y_time_0.assign(num_blocks, temp_double_v2);
		z_time_0.assign(num_blocks, temp_double_v2);
		//number of boundary crossings in each direction
		cross_x.assign(num_part, 0);
		cross_y.assign(num_part, 0);
		cross_z.assign(num_part, 0);

		//Velocity Autocorrealation Function (VACF)
		vacf.assign(num_blocks, temp_double_v);
		//velocities at initial time
		v_x_time_0.assign(num_blocks, temp_double_v2);
		v_y_time_0.assign(num_blocks, temp_double_v2);
		v_z_time_0.assign(num_blocks, temp_double_v2);

		//Stress Autocorrealation Function (SACF)
		sacf_xy.assign(num_blocks, temp_double_v);
		sacf_yz.assign(num_blocks, temp_double_v);
		sacf_zx.assign(num_blocks, temp_double_v);
		//microscopic stress tensors at initial time
		mst_xy_time_0.assign(num_blocks, temp_double_v);
		mst_yz_time_0.assign(num_blocks, temp_double_v);
		mst_zx_time_0.assign(num_blocks, temp_double_v);
		break;
	}
	case 1:
		//loop over all the blocks to test which blocks need sampling
		for (int block = 0; block < num_blocks; block++)
			//test for blocking operation, i.e. when "time_counter" is a multiple of "num_block_elem^block"
			if (time_counter % static_cast<int>(std::pow(num_block_elem / num_points_log10, block)) == 0)
			{
				//increase the current block-length
				block_length[block]++;
				//compute the current length of the block, limited to size "num_block_elem"
				int curr_block_length = std::min(block_length[block], num_block_elem);
				//loop over the block elements
				for (int block_elem = 0; block_elem < curr_block_length; block_elem++)
				{
					//loop over the particles
					for (int i = 0; i < num_part; i++)
					{
						//set last index to the correlation value
						if (block_elem == curr_block_length - 1)
						{
							//store coordinates of initial time, take into account the boundary crossings
							x_time_0[block][curr_block_length - 1][i] = x[i] + cross_x[i];
							y_time_0[block][curr_block_length - 1][i] = y[i] + cross_y[i];
							z_time_0[block][curr_block_length - 1][i] = z[i] + cross_z[i];

							//store velocities of initial time
							v_x_time_0[block][curr_block_length - 1][i] = v_x[i];
							v_y_time_0[block][curr_block_length - 1][i] = v_y[i];
							v_z_time_0[block][curr_block_length - 1][i] = v_z[i];
						}
						//shift to the right if the block is full
						else if (block_length[block] > num_block_elem)
						{
							x_time_0[block][block_elem][i] = x_time_0[block][block_elem + 1][i];
							y_time_0[block][block_elem][i] = y_time_0[block][block_elem + 1][i];
							z_time_0[block][block_elem][i] = z_time_0[block][block_elem + 1][i];

							v_x_time_0[block][block_elem][i] = v_x_time_0[block][block_elem + 1][i];
							v_y_time_0[block][block_elem][i] = v_y_time_0[block][block_elem + 1][i];
							v_z_time_0[block][block_elem][i] = v_z_time_0[block][block_elem + 1][i];
						}

						//calc. MSD
						double dx = x[i] + cross_x[i] - x_time_0[block][block_elem][i];
						double dy = y[i] + cross_y[i] - y_time_0[block][block_elem][i];
						double dz = z[i] + cross_z[i] - z_time_0[block][block_elem][i];
						msd[block][curr_block_length - 1 - block_elem] += dx * dx + dy * dy + dz * dz;

						//calc. VACF
						vacf[block][curr_block_length - 1 - block_elem] += v_x[i] * v_x_time_0[block][block_elem][i]
							+ v_y[i] * v_y_time_0[block][block_elem][i] + v_z[i] * v_z_time_0[block][block_elem][i];
					}

					//set last index to the correlation value
					if (block_elem == curr_block_length - 1)
					{
						//store values of microscopic stress tensor for initial time
						mst_xy_time_0[block][curr_block_length - 1] = kin_mst_xy + pot_mst_xy;
						mst_yz_time_0[block][curr_block_length - 1] = kin_mst_yz + pot_mst_yz;
						mst_zx_time_0[block][curr_block_length - 1] = kin_mst_zx + pot_mst_zx;
					}
					//shift to the right if the block is full
					else if (block_length[block] > num_block_elem)
					{
						mst_xy_time_0[block][block_elem] = mst_xy_time_0[block][block_elem + 1];
						mst_yz_time_0[block][block_elem] = mst_yz_time_0[block][block_elem + 1];
						mst_zx_time_0[block][block_elem] = mst_zx_time_0[block][block_elem + 1];
					}
					//calc. SACF
					sacf_xy[block][curr_block_length - 1 - block_elem] += (kin_mst_xy + pot_mst_xy)
						* mst_xy_time_0[block][block_elem];
					sacf_yz[block][curr_block_length - 1 - block_elem] += (kin_mst_yz + pot_mst_yz)
						* mst_yz_time_0[block][block_elem];
					sacf_zx[block][curr_block_length - 1 - block_elem] += (kin_mst_zx + pot_mst_zx)
						* mst_zx_time_0[block][block_elem];

					//count the number of sampling
					msd_vacf_safc_samples[block][curr_block_length - 1 - block_elem]++;
				}
			}
		//count the current sampling
		time_counter++;
		break;
	case 2:
	{
		std::ofstream msd_file(title_of_msd, std::ios::app);
		std::ofstream vacf_file(title_of_vacf, std::ios::app);
		std::ofstream sacf_file(title_of_sacf, std::ios::app);
		for (int block = 0; block < num_blocks; block++)
		{
			int start_block_elem;
			if (block > 0)
				start_block_elem = num_points_log10;
			else
				start_block_elem = 0;
			int curr_block_length = std::min(block_length[block], num_block_elem);
			for (int block_elem = start_block_elem; block_elem < curr_block_length; block_elem++)
			{
				double time = static_cast<double>(block_elem) * std::pow(num_block_elem / num_points_log10, block) * time_step;
				double divisor = static_cast<double>(num_part) * static_cast<double>(msd_vacf_safc_samples[block][block_elem]);

				//write MSD to the file
				msd_file << time << '\t' << msd[block][block_elem] / divisor << '\n';

				//write VACF to the file
				vacf_file << time << '\t' << vacf[block][block_elem] / divisor << '\n';

				//write SACF to the file
				sacf_file << time << '\t' << (sacf_xy[block][block_elem] + sacf_yz[block][block_elem] + sacf_zx[block][block_elem])
					* coupl_param / (9.0 * divisor) << '\n';
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

//MD CORE
int main()
{
	//keep track of time (start)
	auto start = std::chrono::high_resolution_clock::now();

	//initialization step
	int norm_time = -1;

	//initialize the system
	Initialization();

	//calc. previous & current time forces
	Force(norm_time);
	for (int i = 0; i < num_part; i++)
	{
		F_prev_x[i] = F_next_x[i];
		F_curr_x[i] = F_next_x[i];
		F_prev_y[i] = F_next_y[i];
		F_curr_y[i] = F_next_y[i];
		F_prev_z[i] = F_next_z[i];
		F_curr_z[i] = F_next_z[i];
	}

	//start the loop
	for (norm_time = 0; norm_time < norm_thermo_time + norm_free_time + norm_meas_time; norm_time++)
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
			std::cout << "time: " << static_cast<double>((norm_time - (norm_thermo_time + norm_free_time))) * time_step << '\n'
				<< "pot_en: " << pot_en << '\n' << "coupl_param: " << 1.5 / kin_en << '\n';

			//write energies to the file
			std::ofstream energy_file(title_of_energy, std::ios::app);
			energy_file << static_cast<double>((norm_time - (norm_thermo_time + norm_free_time))) * time_step << '\t'
				<< pot_en << '\t' << kin_en << '\t' << 1.5 / kin_en << '\n';
			energy_file.close();

			//write coordinates to the file
			if (norm_time >= norm_thermo_time + norm_free_time)
			{
				std::ofstream coordinate_file(title_of_coordinates, std::ios::app);
				coordinate_file << num_part << '\n'
					<< static_cast<double>((norm_time - (norm_thermo_time + norm_free_time))) * time_step << '\n';
				for (int i = 0; i < num_part; i++)
					coordinate_file << i << '\t' << x[i] << '\t' << y[i] << '\t' << z[i] << '\n';
				coordinate_file.close();
			}
		}
	}

	//write properties of the system to the file
	Radial_Distribution_Function(2);
	Dispersion(2);
	MSD_VACF_SACF(2);

	//keep track of time (finish)
	auto finish = std::chrono::high_resolution_clock::now();
	std::chrono::duration<double> duration = finish - start;
	std::ofstream simulation_time("simulation time.dat", std::ios::app);
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