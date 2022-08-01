#include <iostream>
#include "smoke.h"

int main()
{
    // instantiate the class with number of grid point per side "grid"
    // and time step "dt" as well as diffusion constant
    int grid = 10;
    double dt = 0.1;
    double diff = 2;

    ClassSolverNS solver(grid, dt, diff);

    int x = 100 % 12;
    int y = 100 - x*(grid+2);
    std::cout << "x = " << x << std::endl; 
    std::cout << "y = " << y << std::endl; 
    
    //solver.reverseIndices(99, x, y);

    //solver.initialize_density();
    //solver.print();
    //solver.diffuse(20);
    //solver.print();
    //solver.comparison();

    return 0;
}