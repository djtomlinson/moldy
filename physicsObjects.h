#ifndef PHYSOBJ_H
#define PHYSOBJ_H

//Requires typeAliases.h,<cmath>,<iostream>
#include <string>

class Bead
{
private:
    int m_id {};
    vec3 m_position {};
    vec3 m_velocity {};
    double m_mass {};
    double m_radius {};
    double m_charge {};
    int m_species {};
    static inline int s_idGenerator {};
public:
    //Constructor (currently constructs too large bead)
    Bead(vec3 position={0,0,0}, vec3 velocity={0,0,0}, double mass=1.0, double radius=1.0, double charge=0.0,
        int species = -1)
        : m_position {position}
        , m_velocity {velocity}
        , m_mass {mass}
        , m_radius {radius}
        , m_charge {charge}
        , m_species {species}
        , m_id {s_idGenerator++}
        {
            //std::cout << "ID: " << m_id << '\n';
        }
    
    void print()
    {
        std::cout << "Position\t" << '(' << m_position[0] << ',' << m_position[1] << ',' << m_position[2] << ")\n"
                  << "Velocity\t" << '(' << m_velocity[0] << ',' << m_velocity[1] << ',' << m_velocity[2] << ")\n"
                  << "Mass\t\t" << m_mass << '\n'
                  << "Radius\t\t" << m_radius << '\n'
                  << "Charge\t\t" << m_charge << '\n';
    }
    
    //Getters
    int getId() const { return m_id; }
    const vec3& getPosition() const { return m_position; }
    const vec3& getVelocity() const { return m_velocity; }
    double getMass() const { return m_mass; }
    double getRadius() const { return m_radius; }
    double getCharge() const { return m_charge; }
    int getSpecies() const { return m_species; }

    //Setters
    void setPosition(vec3& position){ m_position = position; }
    void setVelocity(vec3& velocity){ m_velocity = velocity; }
    void setMass(double mass){ m_mass = mass; }
    void setRadius(double radius){ m_radius = radius; }
    void setCharge(double charge){ m_charge = charge; }
    void setSpecies(int species){ m_species = species;}

    //Physics steps
    void addVelocity(vec3 velocityChange)
    {
        m_velocity[0] += velocityChange[0];
        m_velocity[1] += velocityChange[1];
        m_velocity[2] += velocityChange[2];

    }
    void updatePosition(double timestep, char boundaryConditions='p')
    {
        m_position[0] += (m_velocity[0] * timestep);
        m_position[1] += (m_velocity[1] * timestep);
        m_position[2] += (m_velocity[2] * timestep);
    }
    vec3 binaryOverlapForce(Bead& otherBead, double repulsionStrength)
    {
        double beadSeparation {sqrt((m_position[0]-otherBead.m_position[0])
                                   *(m_position[0]-otherBead.m_position[0])
                                   +(m_position[1]-otherBead.m_position[1])
                                   *(m_position[1]-otherBead.m_position[1])
                                   +(m_position[2]-otherBead.m_position[2])
                                   *(m_position[2]-otherBead.m_position[2]))};
                                   
        if (beadSeparation < (m_radius+otherBead.m_radius))
        {
            //std::cout << beadSeparation << '\n';
            double invSeparation {1/beadSeparation};
            return {repulsionStrength*(m_position[0]-otherBead.m_position[0])*invSeparation,
                    repulsionStrength*(m_position[1]-otherBead.m_position[1])*invSeparation,
                    repulsionStrength*(m_position[2]-otherBead.m_position[2])*invSeparation};
        }
        else
        {
            return {0,0,0};
        }
    }
    vec3 lennardJonesForce(Bead& otherBead, double repulsionStrength)
    {
        double beadSeparation {sqrt((m_position[0]-otherBead.m_position[0])
                                   *(m_position[0]-otherBead.m_position[0])
                                   +(m_position[1]-otherBead.m_position[1])
                                   *(m_position[1]-otherBead.m_position[1])
                                   +(m_position[2]-otherBead.m_position[2])
                                   *(m_position[2]-otherBead.m_position[2]))};
        double invSeparation {(m_radius+otherBead.m_radius)/beadSeparation}; //includes numer.
        double invSeparation5 {invSeparation*invSeparation*invSeparation*invSeparation
                              *invSeparation};
        double invSeparation11 {invSeparation5*invSeparation5*invSeparation};
        
        double forceStrength {4*repulsionStrength*(-12*invSeparation11+6*invSeparation5)}; //LJ strength
        
        invSeparation = 1/beadSeparation; //does not include numer.

        return {-forceStrength*(m_position[0]-otherBead.m_position[0])*invSeparation,
                -forceStrength*(m_position[1]-otherBead.m_position[1])*invSeparation,
                -forceStrength*(m_position[2]-otherBead.m_position[2])*invSeparation};        
    }

    //Friends
    friend class Spring;
};

class Spring
{
private:
    int m_id {};
    int m_start {};
    int m_end {};
    double m_naturalLength {};
    double m_springConstant {};
    static inline int s_idGenerator {};
public:
    //Constructor
    Spring(int start=0, int end=1, double naturalLength=1.0, double springConstant=1.0)
        : m_start {start}
        , m_end {end}
        , m_naturalLength {naturalLength}
        , m_springConstant {springConstant}
        , m_id {s_idGenerator++}
        {
            //std::cout << "ID: " << m_id << '\n';
        }

    //Getters
    int getId() const { return m_id; }
    int getStart() const { return m_start; }
    int getEnd() const { return m_end; }
    double getNaturalLength() const { return m_naturalLength; }
    double getSpringConstant() const { return m_springConstant; }

    //Setters
    void setStart(int start) { m_start = start; }
    void setEnd(int end) { m_end = end; }
    void setNaturalLength(double naturalLength) { m_naturalLength = naturalLength; }
    void setSpringConstant(double springConstant) { m_springConstant = springConstant; }

    //Physics steps
    vec3 springForce(Bead& beadStart, Bead& beadEnd)
    {
        vec3 startPosition {beadStart.m_position};
        vec3 endPosition {beadEnd.m_position};
        double beadSeparation {sqrt((startPosition[0]-endPosition[0])*(startPosition[0]-endPosition[0])
                                   +(startPosition[1]-endPosition[1])*(startPosition[1]-endPosition[1])
                                   +(startPosition[2]-endPosition[2])*(startPosition[2]-endPosition[2]))};
        double forceMagnitude {m_springConstant*(m_naturalLength-beadSeparation)};
        
        double invBeadSeparation {1/beadSeparation};

        return {forceMagnitude * (startPosition[0]-endPosition[0]) * invBeadSeparation,
                forceMagnitude * (startPosition[1]-endPosition[1]) * invBeadSeparation,
                forceMagnitude * (startPosition[2]-endPosition[2]) * invBeadSeparation};
    }

};

#endif

