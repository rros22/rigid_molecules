#include "sites.hpp"
#include <cmath>

void force_input_data::force_input_component_allocator(int no_a, int no_b)
{
    //a and b represent number of molecules of each type
    int lj_no_a = 1;
    int charge_no_a = 3;
    int lj_no_b = 2;
    int charge_no_b = 0;
    //size of memory chunks for a and be
    int size_a = no_a*lj_no_a*sizeof(lj_force_input) + no_a*charge_no_a*sizeof(charge_force_input);
    int size_b = no_b*lj_no_b*sizeof(lj_force_input) + no_b*charge_no_b*sizeof(charge_force_input);
    //allocate all memory upfront in a contiguous array
    void* buffer = malloc(size_a + size_b);
    //push pointers to mark memory areas
    lj_a = (lj_force_input*) buffer;
    lj_b = lj_a + no_a*lj_no_a; //shift by number of lj_sites a
    charge_a = (charge_force_input*) (lj_b + no_b*lj_no_b); //shift by number of lj_sites b
    charge_b = charge_a + no_a*charge_no_a;
}
