#ifndef POST_PROCESSING_HPP
#define POST_PROCESSING_HPP

#include "../rigid_molecules/rigid_molecules.hpp"

#include <array>

#include <fstream>
#include <iostream>

std::string white_spaces(int spaces);
std::string to_string(double angstroms);
void water_pdb(water_site_positions* molecule, int& site_counter, std::string path);
void water_buffer_pdb(h2o_buffer* water_molecules, std::string path);

void terminate_pbd(std::string path);

#endif
