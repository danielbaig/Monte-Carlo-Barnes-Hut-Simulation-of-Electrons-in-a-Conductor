#include "ConductorElectronSystem_Node.h"
#include "Random.h"

#define _USE_MATH_DEFINES

#include <cmath>
#include <iostream>
#include <vector>
#include <array>
#include <fstream> // Saving to file
#include <algorithm> // std::remove, std::fill_n
#include <numeric> // std::iota
#include <tuple> // std::tuple, std::make_tuple

#include <math.h> // exp
#include <stdlib.h> // abs


Node* buildQuadTree(std::vector<Node*> particles, int numParticles)
{
    /*
    Generate a quadtree for all particles.

    Inputs:
        - std::vector<Node*> particles: Particles in the system.
        - int numParticles: The number of particles in the system.

    Outputs:
        - Node* root: The initial node.
    */

    // Create initial node.
    Node* root{ new Node{nullptr, 0,0,0,0,
        2 * ConductorElectronSystem::getSphereRadius() } };

    for (int i{ 0 }; i < numParticles; ++i)
    {
        root->insertQuadtree(particles[i]);
    }

    root->calcCOM();


    return root;
}




Node::Node(Node* parent, double x, double y, double mass,
    double charge,
    double sideLength,
    bool COM_calculated) :
    m_centrex{ x },
    m_centrey{ y },
    m_mass{ mass },
    m_charge{ charge },
    m_sideLength{ sideLength },
    m_COM_calculated{ COM_calculated },
    m_parent{parent}
{
    /*
    Constructor for the node class.

    Inputs:
        - Node* parent: Node's parent.
        - double x: x-coord of position of node.
        - double y: y-coord of position of node.
        - double mass: mass of the node.
        - double charge: Charge of the node.
        - double sideLength: Length of the side of the square the node represents.
        - bool m_COM_calculated: Whether the position is also the COM.
        */

    m_r = { m_centrex,m_centrey };
    if (m_COM_calculated)
    {
        m_x = m_centrex;
        m_y = m_centrey;
    }


}




Node::~Node()
    {
    /*
    Destructor for the node class
    */
    
        for (auto& child : m_children)
        {
            if (child != nullptr)
            {
                delete child;
                child = nullptr;
            }
        }

        

    }

void Node::calcCOM()
{
    /*
    Calculates the centre of mass of a node using recursion.
    */

    if (!m_COM_calculated)
    {
        if (m_numMembers == 1)
        {
            // Sets proporties to that of its only
            // member
            m_x = m_member->m_x;
            m_y = m_member->m_y;
            m_mass = m_member->m_mass;
            m_charge = m_member->m_charge;
        }
        else
        {
            // Calculate COM of child nodes
            for (auto& child : m_children)
            {
                if (child!=nullptr)
                {
                    child->calcCOM();

                    m_x += child->m_x
                        * child->m_mass;
                    m_y += child->m_y
                        * child->m_mass;
                    m_mass += child->m_mass;
                    m_charge += child->m_charge;
                }
            }
            m_x /= m_mass;
            m_y /= m_mass;
        }

        m_COM_calculated = true;
    }
}

Node* Node::addChild(int quadrant)
{
    /*
    Adds a new child node.

    Inputs:
        - int quadrant: The quadrant the child
        resides in.

    Outputs:
        - Node* child: The child node.
    */



    double dx{ quadrant % 2 == 0 ? -1. : 1. };
    double dy{ quadrant < 2 ? -1. : 1. };
    dx *= m_sideLength / 4;
    dy *= m_sideLength / 4;

    Node* node_ptr{ new Node(this, m_centrex + dx, m_centrey + dy, 0, 0,
        m_sideLength / 2) };

    return node_ptr;
}

int Node::findQuadrant(Node* particle) const
    {
        /*
        Find the quadrant the particle is in.
        Quadrant indexing as:
            2 3
            0 1

        Inputs:
            - Node* particle: The particle to locate the
            quadrant of.
            - int quadrant: The quadrant the particle is in.
        */

        int alpha{ particle->m_centrex < m_centrex ? 0 : 1 };
        int beta{ particle->m_centrey < m_centrey ? 0 : 1 };
        return alpha + 2 * beta;


    }


void Node::insertQuadtree(Node* particle)
{
    /*
    Insert particle at a node.

    Inputs:
        - Node* particle: The particle to insert.
    */

    m_numMembers += 1;

    if (m_numMembers > 2)
    {
        int quadrant{ findQuadrant(particle) };

        // Go into child quadrant
        if (m_children[quadrant] == nullptr)
        {
            m_children[quadrant] = addChild(quadrant);
        }
        m_children[quadrant]->insertQuadtree(particle);

    }

    else if (m_numMembers == 2)
    {

        // Original particle reallocation
        int quadrant{ findQuadrant(m_member) };
        m_children[quadrant] = addChild(quadrant);
        m_children[quadrant]->m_member = m_member;
        m_children[quadrant]->m_numMembers += 1;



        // Reset member
        m_member = nullptr;
        // New particle allocation
        quadrant = findQuadrant(particle);
        // Go into child quadrant
        if (m_children[quadrant] == nullptr)
        {
            m_children[quadrant] = addChild(quadrant);
        }

        m_children[quadrant]->insertQuadtree(particle);
    }

    else if (m_numMembers == 1)
    {
        m_member = particle;
        particle->m_parent = this;

    }
}



double Node::BH_potential(Node* node) const
{
    /*
    Determine the potential on a particle using the Barnes-Hut method.
    
    Inputs:
        - Node* node: Node to calculate the potential from.

    Outputs:
        - double energy: The energy due to the particle interacting
        with the potential.
    */

    double energy{ 0 };
    double dx{ m_x - node->m_x };
    double dy{ m_y - node->m_y };
    double dr{ sqrt(dx * dx + dy * dy) };

    if (dr != 0)
    {
        // Only one member or can be considered to be so.
        if (node->m_numMembers == 1 || node->m_sideLength / dr
            < m_theta)
        {
            energy = ConductorElectronSystem::calculatePairPotential(this, node);
        }
        else
        {
            // Apply algorithm to the subnodes.
            for (int i{ 0 }; i < 4; ++i)
            {
                if (node->m_children[i] != nullptr)
                {
                    energy += BH_potential(node->m_children[i]);
                }
            }
        }
    }

    return energy;
}




ConductorElectronSystem::ConductorElectronSystem(
    unsigned int numParticles,
    double temperature, unsigned int Nsteps)
    : m_temperature{ temperature },
    m_numElectrons{numParticles},
    m_Nsteps{Nsteps}
{
    /*
    Constructor.
    Sets the temperature and assigns the electron positions as well as creates 
    the node root.

    Inputs:
        - std::vector<Node*> particles: Instance of the particles.
        - unsigned int numParticles: The number of particles.
        - double temperature: temperature of the system.
        - unsigned int Nsteps: The number of Monte Carlo steps to perform.
    */

    m_particles.resize(m_numElectrons);

    double x{};
    double y{};
    double angle{};
    double distance{};

    for (unsigned int i{ 0 }; i < m_numElectrons; ++i)
    {
        angle = 2 * M_PI * Random::uniform();
        // http://www.anderswallin.net/2009/05/uniform-random-points-in-a-circle-using-polar-coordinates/
        distance = ConductorElectronSystem::getSphereRadius() * sqrt(Random::uniform());

        x = distance * cos(angle);
        y = distance * sin(angle);

        Node* node_ptr{ new Node(nullptr, x, y, ELECTRONMASS, -ELECTRONCHARGE, 0, true) };
        m_particles[i] = node_ptr;

    }

    m_root = buildQuadTree(m_particles, m_numElectrons);

    m_externalCharge =  new Node{ nullptr, -1.5 * m_sphereRadius,0,
        ELECTRONMASS * m_numElectrons, 1e+2*ELECTRONCHARGE *
        m_numElectrons, 0, true};


}

ConductorElectronSystem::~ConductorElectronSystem()
{
    /*
    Destructor.

    Deletes dynamic arrays.
    */

    for (unsigned int i{ 0 }; i < m_numElectrons; ++i)
    {
        delete m_particles[i];
        m_particles[i] = nullptr;
    }

    delete m_externalCharge;
    m_externalCharge = nullptr;
    delete m_root;
    m_root = nullptr;
}


void ConductorElectronSystem::savePositions() const
    {
        /*
        Saves the positions of the electrons.
        */
        std::ofstream myfile;
        myfile.open("../data/positions.txt");
        for (unsigned int i{ 0 }; i < m_numElectrons; ++i)
        {
            myfile << m_particles[i]->m_x << ','
                << m_particles[i]->m_y << '\n';
        }
        myfile.close();
    }

void ConductorElectronSystem::saveEnergyHistory(std::vector<double> energyHistory,
    unsigned int Nsteps) const
    {
        /*
        Saves the positions of the electrons.

        Inputs:
            - std::vector<double> energyHistory: The history of
            all the energies of the system
            - unsigned int Nsteps: The number of MC steps that
            were performed.
        */

        std::ofstream myfile;
        myfile.open("../data/energyHistory.txt");
        for (const auto& i : energyHistory)
        {
            myfile << i << ',';
        }
        myfile.close();
    }


double ConductorElectronSystem::calculatePairPotential(const Node* obj1, const Node* obj2)
    {
        /*
        Calculates the potential between two objects.

        Inputs:
            - const Node* obj1: Instance of first object.
            - const Node* obj2: Instance of second object.
        */

        double dx{ obj1->m_x - obj2->m_x };
        double dy{ obj1->m_y - obj2->m_y };
        double dr{ sqrt(dx * dx + dy * dy) };

        return COULOMB * obj1->m_charge * obj2->m_charge / dr;
    }

void ConductorElectronSystem::performMCSimulation()
{
    /*
    Performs a Monte Carlo simulation to find the equilibrium
    positions of the electrons in the conductor.
    */

    static const double moveMax{ 0.45*m_sphereRadius };
    unsigned int numAccepted{ 0 };


    double energyTotal{ electron_electronPotential(false)
        + electron_externalPotential()
        + conductorBarrierPotential() };

    // Loop variables
    std::vector<double> energyHistory(m_Nsteps);
    double movementAmount{};
    double direction{};
    int electronToMove{};
    double xprev{};
    double yprev{};
    double newEnergy{};
    double energyChange{};
    bool recreateQuadtree{true};

    for (unsigned int n{ 0 }; n < m_Nsteps; ++n)
    {
        energyHistory[n] = energyTotal;

        if ((n + 1) % (m_Nsteps / 10) == 0)
        {
            std::cout << 100 * (n + 1) / m_Nsteps << "% \n";
        }

        // Determine where the electron should move
        // http://www.anderswallin.net/2009/05/uniform-random-points-in-a-circle-using-polar-coordinates/
        movementAmount = moveMax * sqrt(Random::uniform());
        direction = 2 * M_PI * Random::uniform();

        electronToMove = Random::intUniform(0, m_numElectrons);

        xprev = m_particles[electronToMove]->m_x;
        yprev = m_particles[electronToMove]->m_y;

        m_particles[electronToMove]->m_x += movementAmount * cos(direction);
        m_particles[electronToMove]->m_y += movementAmount * sin(direction);

        // Calculate the new energy
        if (m_particles[electronToMove]->m_x - m_particles[electronToMove]->m_parent->m_centrex
            < m_particles[electronToMove]->m_parent->m_sideLength
            && m_particles[electronToMove]->m_y - m_particles[electronToMove]->m_parent->m_centrey
            < m_particles[electronToMove]->m_parent->m_sideLength)
        {
            recreateQuadtree = false;
            m_particles[electronToMove]->m_parent->m_x = m_particles[electronToMove]->m_x;
            m_particles[electronToMove]->m_parent->m_y = m_particles[electronToMove]->m_y;
        }
        else recreateQuadtree = true;


        newEnergy = electron_electronPotential(recreateQuadtree)
            + electron_externalPotential()
            + conductorBarrierPotential();

        energyChange = newEnergy - energyTotal;

        // Determine if the change should be accepted
        if (energyChange < 0 || Random::uniform()
            < exp(-energyChange / (BOLTZMANN * m_temperature)))
        {
            energyTotal = newEnergy;
            ++numAccepted;
        }

        else
        {
            recreateQuadtree = false;
            m_particles[electronToMove]->m_parent->m_x = xprev;
            m_particles[electronToMove]->m_parent->m_y = yprev;
            m_particles[electronToMove]->m_x = xprev;
            m_particles[electronToMove]->m_y = yprev;
        }
    }

    saveEnergyHistory(energyHistory, m_Nsteps);
    std::cout << "Acceptance probability: " 
        << 100*static_cast<double>(numAccepted) / static_cast<double>(m_Nsteps) << "%\n";
}





double ConductorElectronSystem::electron_electronPotential(bool recreateQuadtree)
{
    /*
    Calculates the potential between all the electrons.

    Inputs:
        - bool recreateQuadtree: If the quadtree needs to be recreated or not.

    Outputs:
        - double energy: Energy between the electrons.
    */


    double energy{ 0. };

    if (recreateQuadtree)
    {
        // Create the quadtree
        delete m_root;
        m_root = buildQuadTree(m_particles, m_numElectrons);
    }


    // Calculate energy between electrons
        
    for (unsigned int i{ 0 }; i < m_numElectrons; ++i)
    {
        energy += m_particles[i]->BH_potential(m_root);
    }

    return energy;
}

double ConductorElectronSystem::electron_externalPotential() const
    {
        /*
        Calculates the potential between all the electrons and
        the external charges.

        Outputs:
            - double energy: The energy due to the external charge.
        */
        double energy{ 0 };

        // Calculate energy between conductor and external charge
        for (unsigned int i{ 0 }; i < m_numElectrons; ++i)
        {
            energy += calculatePairPotential(
                m_particles[i], m_externalCharge);
        }

        return energy;
    }

double ConductorElectronSystem::conductorBarrierPotential() const
    {
        /*
        Calculates the potential of electrons "outside" of the
        conductor.

        Outputs:
            - double energy: Energy due to electrons being outside of the conductor.
        */
        static const double externalPotentialGrad{ m_sphereRadius*1e+4 };

        double energy{ 0. };
        double distance{};

        for (unsigned int i{ 0 }; i < m_numElectrons; ++i)
        {
            distance = sqrt(m_particles[i]->m_x * m_particles[i]->m_x
                + m_particles[i]->m_y * m_particles[i]->m_y);

            if (distance > m_sphereRadius)
            {
                energy += COULOMB*ELECTRONCHARGE*(exp((distance
                    - m_sphereRadius) / externalPotentialGrad) - 1);
            }
        }

        return energy;
    }



double ConductorElectronSystem::getSphereRadius()
{
    /*
    Obtains the sphere radius.
    */

    return m_sphereRadius;

}



