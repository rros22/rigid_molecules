#include "parse.hpp"
#include <iostream>

bool parse_molecules(int argc, char* argv[], unsigned &molecule_no)
{
  if (argv[1] == NULL)
  {
    std::cout << "Program executed without providing enough arguments"
    << std::endl << "Number of molecules unspecified" << std::endl;
    return 1;
  }

  molecule_no = std::stoi(argv[1]);
  return 0;
}
