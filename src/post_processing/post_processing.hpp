#ifndef POST_PROCESSING_HPP
#define POST_PROCESSING_HPP

#include "../box/box.hpp"
#include "../rigid_molecules/rigid_molecules.hpp"

#include <array>

#include <fstream>
#include <iostream>

std::string white_spaces(int spaces);
std::string to_string(double angstroms);
void box_pdb(Box *box, std::string path);

#endif
