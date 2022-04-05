#include "post_processing.hpp"
#include <iomanip>

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

void molecule_pdb(rigid_molecule *molecule, int &site_counter, std::string path){

  std::ofstream file;
  file.open(path, std::fstream::app);

  file << std::fixed;
  file << std::setprecision(3);

  //number of interaction sites of molecule
  int site_no = molecule->return_sites_list().size();

  //molecule sites pointer, to current site
  site *site;

  //lenght of field parameters
  int no_len;
  int name_len;
  int x_len;
  int y_len;
  int z_len;
  int sym_len;

  //coordinates of current site
  std::array<double, 3> coordinates;

  //cartesian coordinates
  //loop through site list

  for (int j = 0; j < site_no - 1; j++){

    //insert first column, element designator
    file << "ATOM";

    //calculate no_length (site number length)
    no_len = std::to_string(site_counter).size();

    //insert white spaces for justification, then site number
    file << white_spaces(11 - 4 - no_len);
    file << std::to_string(site_counter);

    //calculate name_length (S for site + site number)
    name_len = ("S" + std::to_string(site_counter)).size();

    //insert with spaces for justification, then site_name
    file << white_spaces(2);
    file << "S" + std::to_string(site_counter);

    //get site from site list
    site = (molecule->return_sites_list())[j];

    //get global coordinates from site
    coordinates = molecule->return_coordinates_site(j);

    //get length of each coordinate
    x_len = to_string(coordinates[0]).size();
    y_len = to_string(coordinates[1]).size();
    z_len = to_string(coordinates[2]).size();

    //insert white space, right justification
    file << white_spaces(38 - 13 - x_len - name_len);
    file << to_string(coordinates[0]);

    file << white_spaces(46 - 38 - y_len);
    file << to_string(coordinates[1]);

    file << white_spaces(54 - 46 - z_len);
    file << to_string(coordinates[2]);

    //element symbol
    sym_len = site->get_symbol().size();

    file << white_spaces(78 - 54 - sym_len);
    file << site->get_symbol();


    //increment element no
    site_counter += 1;

    //insert new line character
    file << std::endl;

  }

}

void terminate_pbd(std::string path){

  std::ofstream file;
  file.open(path, std::fstream::app);

  file << "ENDMDL";

  file << std::endl;

}


void box_pdb(Box *box, std::string path){

    //get vector of pointers to molecules, and amount of molecules in vector
    std::vector<rigid_molecule*> molecules = box->get_molecules();
    int molecule_no = molecules.size();

    //define site counter
    int k = 1;

    //write to file all molecules
    for (int i = 0; i < molecule_no; i++){

          molecule_pdb(molecules[i], k, path);


    }

    //terminate frame
    terminate_pbd(path);

    //close file
    //file.close();
}
