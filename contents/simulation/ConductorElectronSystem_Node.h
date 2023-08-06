#ifndef CONDUCTORELECTRONSYSTEM_NODE
#define	CONDUCTORELECTRONSYSTEM_NODE

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


constexpr double ELECTRONMASS{ 9.11e-31 };
constexpr double ELECTRONCHARGE{ 1.6e-19 };
constexpr double BOLTZMANN{ 1.380649e-23 };
constexpr double COULOMB{ 8.9875517923e+9 };


class Node
{
    /*
    Creates a node for the quadtree.

    Child nodes are at the next level.
    If child==None then it must be a particle.

    Member variables:
        - double m_x: x-coord of COM
        - double m_y: y-coord of COM
        - std::array<double,2> m_r: Vector of COM.
        - const double m_centrex: Geometric x coord of the centre of the node.
        - const double m_centrey: Geometric y coord of the centre of the node.
        - double m_mass: Mass of the node.
        - double m_charge: Charge of the node.
        - const double m_sideLength: Length of a side of the square the node represents.
        - bool m_COM_calculated: Whether the COM has been calculated.
        - int m_numMembers: The number of members of the node.
        - Node* m_member: Single direct member of the node.
        - Node* m_children[4]: The child nodes.
        - static double m_theta: Distance sidelength threshold ratio.
        - Node* m_parent: Parent node.
    */
private:

    double m_x{};
    double m_y{};
    std::array<double, 2> m_r{};

    const double m_centrex{};
    const double m_centrey{};

    double m_mass{};
    double m_charge{};
    const double m_sideLength{};
    bool m_COM_calculated{};
    int m_numMembers{ 0 };
    Node* m_member{nullptr};
    Node* m_children[4]{nullptr, nullptr, nullptr, nullptr};
    static double m_theta;
    
    Node* m_parent{};



public:
    Node(Node* parent, double x = 0, double y = 0, double mass = 0,
        double charge = 0,
        double sideLength = 0,
        bool COM_calculated = false);


    ~Node();

    void calcCOM();

    Node* addChild(int quadrant);

    int findQuadrant(Node* particle) const;


    void insertQuadtree(Node* particle);


    double BH_potential(Node* node) const;


    friend class ConductorElectronSystem;
};


class ConductorElectronSystem
{
    /*
    Creates the system to model the electrons in a
    conductor and the external charges that impact the
    internal electrons.

    Member variables:
        - const unsigned int m_Nsteps: The number of MC steps to perform
        - const unsigned int m_numElectrons: The number of electrons used.
        - static inline constexpr double m_sphereRadius: The radius of the
        conducting sphere (circle used instead).
        - const double m_temperature: The temperature of the system.
        - const std::string m_shape: The shape of the conductor (not used in code).
        - const Node* m_externalCharge: Pointer to the node representing the external
        charge.
        - std::vector<Node*> m_particles: Vector of node pointers of the particle instances.
        - Node* m_root: The root node for the quadtree.
    */


private:
    const unsigned int m_Nsteps{};
    const unsigned int m_numElectrons{};
    static inline constexpr double m_sphereRadius{ 1e-7 }; // [m]
    const double m_temperature{}; // [K]

    const std::string m_shape{ "sphere" };
 
    const Node* m_externalCharge{};


    std::vector<Node*> m_particles{};
    Node* m_root{};


public:


    ConductorElectronSystem(unsigned int numParticles,
        double temperature, unsigned int Nsteps);

    ~ConductorElectronSystem();


    void savePositions() const;

    void saveEnergyHistory(std::vector<double> energyHistory,
        unsigned int Nsteps) const;

    static double calculatePairPotential(const Node* obj1, const Node* obj2);

    void performMCSimulation();

    double electron_electronPotential(bool recreateQuadtree=true);

    double electron_externalPotential() const;

    double conductorBarrierPotential() const;



    static double getSphereRadius();



    friend class Node;
};



#endif

