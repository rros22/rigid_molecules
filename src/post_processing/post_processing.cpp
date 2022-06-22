#include "post_processing.hpp"
#include <iomanip>
#include <memory>
#include <cmath>

std::string to_string(double angstroms){

    double temp = angstroms;

    std::string str(std::to_string(temp));

    int i = str.find('.');

    std::string str_1(str.substr(0, i));
    std::string str_2(str.substr(i, 4));

    return str_1 + str_2;
}

std::string white_spaces(int spaces){

    std::string str(spaces, ' ');

    return str;

}



void water_pdb(water_site_positions* molecule, int& site_counter, std::string path)
{
  std::ofstream file;
  file.open(path, std::fstream::app);
  file << std::fixed;
  file << std::setprecision(3);
  //number of interaction sites of molecule
  int site_no = 4;
  //molecule sites pointer, to current site
  site_positions* site;
  //lenght of field parameters
  int no_len;
  int name_len;
  int x_len;
  int y_len;
  int z_len;
  int sym_len;
  //for O,H1,H2
  //insert first column, element designator
  file << "ATOM";
  //calculate no_length (site number length)
  no_len = floor(log10(site_counter)) + 1;
  //insert white spaces for justification, then site number
  file << white_spaces(11 - 4 - no_len);
  file << std::to_string(site_counter);
  //calculate name_length (S for site + site number)
  name_len = no_len + 1;
  //insert with spaces for justification, then site_name
  file << white_spaces(2);
  file << "S" + std::to_string(site_counter);
  //get site from site list
  site = &(molecule->O);
  //get length of each coordinate
  x_len = to_string(site->x).size();
  y_len = to_string(site->y).size();
  z_len = to_string(site->z).size();
  //insert white space, right justification
  file << white_spaces(38 - 13 - x_len - name_len);
  file << to_string(site->x);
  file << white_spaces(46 - 38 - y_len);
  file << to_string(site->y);
  file << white_spaces(54 - 46 - z_len);
  file << to_string(site->z);
  //element symbol
  sym_len = 1;
  file << white_spaces(78 - 54 - sym_len);
  file << 'O';
  //increment element no
  site_counter += 1;
  //insert new line character
  file << std::endl;
  //insert first column, element designator
  file << "ATOM";
  //calculate no_length (site number length)
  no_len = floor(log10(site_counter)) + 1;
  //insert white spaces for justification, then site number
  file << white_spaces(11 - 4 - no_len);
  file << std::to_string(site_counter);
  //calculate name_length (S for site + site number)
  name_len = no_len + 1;
  //insert with spaces for justification, then site_name
  file << white_spaces(2);
  file << "S" + std::to_string(site_counter);
  //get site from site list
  site = &(molecule->H1);
  //get length of each coordinate
  x_len = to_string(site->x).size();
  y_len = to_string(site->y).size();
  z_len = to_string(site->z).size();
  //insert white space, right justification
  file << white_spaces(38 - 13 - x_len - name_len);
  file << to_string(site->x);
  file << white_spaces(46 - 38 - y_len);
  file << to_string(site->y);
  file << white_spaces(54 - 46 - z_len);
  file << to_string(site->z);
  //element symbol
  sym_len = 1;
  file << white_spaces(78 - 54 - sym_len);
  file << 'H';
  //increment element no
  site_counter += 1;
  //insert new line character
  file << std::endl;
  //insert first column, element designator
  file << "ATOM";
  //calculate no_length (site number length)
  no_len = floor(log10(site_counter)) + 1;
  //insert white spaces for justification, then site number
  file << white_spaces(11 - 4 - no_len);
  file << std::to_string(site_counter);
  //calculate name_length (S for site + site number)
  name_len = no_len + 1;
  //insert with spaces for justification, then site_name
  file << white_spaces(2);
  file << "S" + std::to_string(site_counter);
  //get site from site list
  site = &(molecule->H2);
  //get length of each coordinate
  x_len = to_string(site->x).size();
  y_len = to_string(site->y).size();
  z_len = to_string(site->z).size();
  //insert white space, right justification
  file << white_spaces(38 - 13 - x_len - name_len);
  file << to_string(site->x);
  file << white_spaces(46 - 38 - y_len);
  file << to_string(site->y);
  file << white_spaces(54 - 46 - z_len);
  file << to_string(site->z);
  //element symbol
  sym_len = 1;
  file << white_spaces(78 - 54 - sym_len);
  file << 'H';
  //increment element no
  site_counter += 1;
  //insert new line character
  file << std::endl;
}

void terminate_pbd(std::string path)
{
  std::ofstream file;
  file.open(path, std::fstream::app);

  file << "ENDMDL";

  file << std::endl;
}

void water_buffer_pdb(h2o_buffer* water_molecules, std::string path)
{
  //get vector of pointers to molecules, and amount of molecules in vector
  int molecule_no = water_molecules->n;
  water_site_positions* molecules = water_molecules->water_site_pos;
  //define site counter
  int k = 1;
  //write to file all molecules
  for (int i = 0; i < molecule_no; i++)
  {
    water_pdb(&molecules[i], k, path);
  }
  //terminate frame
  terminate_pbd(path);
}
