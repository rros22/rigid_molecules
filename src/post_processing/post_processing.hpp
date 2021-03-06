#ifndef POST_PROCESSING_HPP
#define POST_PROCESSING_HPP

#include "../box/box.hpp"
#include "../rigid_molecules/rigid_molecules.hpp"

#include <array>

#include <fstream>
#include <iostream>

std::string white_spaces(int spaces);
std::string to_string(double angstroms);
void molecule_pdb(const std::shared_ptr<rigid_molecule>& molecule, int &site_counter, std::string path);
void box_pdb(Box& box, std::string path);

void terminate_pbd(std::string path);

#endif
