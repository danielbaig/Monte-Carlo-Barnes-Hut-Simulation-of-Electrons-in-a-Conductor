#include "Random.h"
#include "Timer.h"
#include "ConductorElectronSystem_Node.h"

#define _USE_MATH_DEFINES

#include <cmath>
#include <iostream>
#include <vector>
#include <array>
#include <fstream> // Saving to file
#include <algorithm> // std::remove
#include <numeric> // std::iota
#include <tuple> // std::tuple, std::make_tuple

#include <math.h> // exp
#include <stdlib.h> // abs


double Node::m_theta = 1;

int main()
{

    constexpr int numElectrons{ 1000 };
    constexpr double temperature{ 300 }; // [K]
    constexpr unsigned int Nsteps{ static_cast<unsigned int>(1e+5) };



    ConductorElectronSystem simulation{numElectrons, temperature, Nsteps};

    Timer t{};
    simulation.performMCSimulation();
    std::cout << "Time taken: " << t.elapsed() 
        << " seconds\n";

    simulation.savePositions();




    return 0;
}
