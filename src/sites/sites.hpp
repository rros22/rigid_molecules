#ifndef SITES_HPP
#define SITES_HPP

#include <iostream>
#include <array>
#include <memory>



struct lj_force_input
{
    //coordinates
    double x;
    double y;
    double z;
    //force parameters specified by calling function
};

struct charge_force_input
{
    //coordinates
    double x;
    double y;
    double z;
    //charge
    double q;
};

struct force_input_data
{
    //pointers to separate parts of memory
    lj_force_input* lj_a;
    lj_force_input* lj_b;
    charge_force_input* charge_a;
    charge_force_input* charge_b;

    void force_input_component_allocator(int no_a, int no_b);
};



#endif
