#include <iostream>
#include <fstream>
#include <sstream>
#include <array>
#include <random>
#include <cmath>

#include "typeAliases.h"
#include "constants.h"
#include "containers.h"
#include "physicsObjects.h"

/*
using vec3 = std::array<double,3>;

class Bead
{
private:
    int m_id {};
    vec3 m_position {};
    vec3 m_velocity {};
    double m_mass {};
    double m_radius {};
    double m_charge {};
    static inline int s_idGenerator {};
public:
    //Constructor (currently constructs too large bead)
    Bead(vec3 position={0,0,0}, vec3 velocity={0,0,0}, double mass=1.0, double radius=1.0, double charge=0.0)
        : m_position {position}
        , m_velocity {velocity}
        , m_mass {mass}
        , m_radius {radius}
        , m_charge {charge}
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

    //Setters
    void setPosition(vec3& position){ m_position = position; }
    void setVelocity(vec3& velocity){ m_velocity = velocity; }
    void setMass(double mass){ m_mass = mass; }
    void setRadius(double radius){ m_radius = radius; }
    void setCharge(double charge){ m_charge = charge; }

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
    //vec3 lennardJones()

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
*/
//force functions START

vec3 stokesDrag(vec3 velocity, double radius, double viscosity) //pass velocity by referen
{
    return {-6*k::pi*radius*viscosity*velocity[0],
            -6*k::pi*radius*viscosity*velocity[1],
            -6*k::pi*radius*viscosity*velocity[2]};
}

vec3 coulombForce(Bead& bead1, Bead& bead2)
{
    double q1 {bead1.getCharge()};
    double q2 {bead2.getCharge()};
    vec3 p1 {bead1.getPosition()};
    vec3 p2 {bead2.getPosition()};

    double invBeadSeparationSq {(p1[0]-p2[0])*(p1[0]-p2[0])
                               +(p1[1]-p2[1])*(p1[1]-p2[1])
                               +(p1[2]-p2[2])*(p1[2]-p2[2])};
    double coeff {k::k_e*invBeadSeparationSq*q1*q2};
    return {coeff*(p1[0]-p2[0]),coeff*(p1[1]-p2[1]),coeff*(p1[2]-p2[2])};
}

vec3 forceToVelocityChange(vec3 force, double mass, double timestep)
{
    /* F = ma
     * a = F/m
     * dv = F/m * dt
     */
    double invMass = 1/mass;
    return {force[0]*invMass*timestep,
            force[1]*invMass*timestep,
            force[2]*invMass*timestep};
}

vec3 generateThermalNoise(double mass, double radius, double temperature, double viscosity)
{
    constexpr double k_B2 {2*k::k_B};
    double stokesCoeff {6*k::pi*radius*viscosity};
    
    //generate gaussian using Box-Muller
    static std::mt19937 rng(std::random_device{}()); 
    static std::uniform_real_distribution<> runif(0.0, 1.0);
    double uniform1;
    double uniform3;
    do
    {
        uniform1 = runif(rng);
        uniform3 = runif(rng);
    }
    while(uniform1 == 0 || uniform3 == 0);
    
    double uniform2 {runif(rng)};
    double uniform4 {runif(rng)};

    constexpr double pi2 {2*k::pi};
    double R12 {sqrt(-2.0*log(uniform1))};
    double R34 {sqrt(-2.0*log(uniform3))};
    
    double normal1 {R12*cos(pi2*uniform2)};
    double normal2 {R12*sin(pi2*uniform2)};
    double normal3 {R34*cos(pi2*uniform4)};
    //4th result not calculated as we only need 3
    
    double magnitude {sqrt(stokesCoeff*k_B2*temperature)};
    return {magnitude*normal1,
            magnitude*normal2,
            magnitude*normal3};
}

void outputBeadPosition(Bead& bead, double timestep)
{
    vec3 position {bead.getPosition()};
    std::cout << '(' << position[0] << ',' << position[1] << ',' << position[2] << ")\n";
}


int main(){

    double maxTime {}; //try to keep in SI
    double timestep {};
    double viscosity {};
    double temperature {};
    double overlapStrength {};
    
    double currentTime {0.0};
    
    std::ifstream simParams {"simulationParameters.txt"};
    std::string params {};
    if (!simParams)
    {
        return 1;
    }
    while(std::getline(simParams, params))
    {
        std::istringstream in {params};
        in >> maxTime >> timestep >> viscosity >> temperature >> overlapStrength;
    }

    std::ofstream beadOut {"beadOutput.txt"};
    if (!beadOut)
    {
        return 1;
    }
    std::ofstream springOut {"springOutput.txt"};
    if (!springOut)
    {
        return 1;
    }
    //copy output to outputBackup folder, numerate


    //some sort of construction loop that loads from file written by python program
    std::ifstream beadInit {"beadInit.txt"};
    std::ifstream beadInitCopy {"beadInit.txt"};
    if (!beadInit || !beadInitCopy)
    {
        std::cout << !beadInit << !beadInitCopy
                  <<" Bead initialisation file not found, exiting...\n";
        return 1;
    }
    int numberBeads {};
    std::string beadLine;
    
    //Dislike this workaround, ineffiencient
    while(std::getline(beadInit,beadLine)) { ++numberBeads; }
    Array<Bead> beadArray(numberBeads);

    int lineNo {};
    while(std::getline(beadInitCopy,beadLine))    
    {
        std::istringstream in {beadLine};
        vec3 position {0,0,0};
        vec3 velocity {0,0,0};
        double mass {};
        double radius {};
        double charge {};
        in >> position[0] >> position[1] >> position[2]
           >> velocity[0] >> velocity[1] >> velocity[2]
           >> mass >> radius >> charge;
        beadArray[lineNo].setPosition(position);
        beadArray[lineNo].setVelocity(velocity);
        beadArray[lineNo].setMass(mass);
        beadArray[lineNo].setRadius(radius);
        beadArray[lineNo].setCharge(charge);
        ++lineNo;
    }

    std::ifstream springInit {"springInit.txt"};
    std::ifstream springInitCopy {"springInit.txt"};
    if (!springInit || !springInitCopy)
    {
        std::cout << !springInit << !springInitCopy
                  <<" Spring initialisation file not found, exiting...\n";
        return 1;
    }
    int numberSprings {};
    std::string springLine;
    while(std::getline(springInit,springLine)) { ++numberSprings; }
    Array<Spring> springArray(numberSprings);
    
    int springNo {};
    while(std::getline(springInitCopy,springLine))
    {
        std::istringstream in {springLine};
        int start {};
        int end {};
        double naturalLength {};
        double springConstant {};
        in >> start >> end >> naturalLength >> springConstant;
        springArray[springNo].setStart(start);
        springArray[springNo].setEnd(end);
        springArray[springNo].setNaturalLength(naturalLength);
        springArray[springNo].setSpringConstant(springConstant);
        ++springNo;
    }
    
    beadOut << "Time"  << ' ' << "X" 
                       << ' ' << "Y"
                       << ' ' << "Z" 
                       << ' ' << "ID" 
                       << ' ' << "Radius" << '\n';
    
    springOut << "Time"  << ' ' << "Start" 
                         << ' ' << "End"
                         << ' ' << "Force" 
                         << ' ' << "ID" << '\n'; 
        
    while(currentTime <= maxTime)
    {
               
        vec3 overlapForce {};
        vec3 overlapVelocity {};
        vec3 chargeForce {};
        vec3 chargeVelocity {};
        
        vec3 springForce {};
        vec3 springVelocity {};
        int startInd {};
        int endInd {};

        vec3 noiseForce {};
        vec3 noiseVelocity {};

        vec3 currentVelocity {};
        vec3 dragForce {};
        vec3 dragVelocity {};

        vec3 outputPosition {};

        for (int i=0; i<springArray.getLength(); ++i)
        {   
            startInd = springArray[i].getStart();
            endInd = springArray[i].getEnd();
            double beadMassStart = beadArray[startInd].getMass();
            double beadMassEndRatio = beadMassStart/beadArray[endInd].getMass();
            springForce = springArray[i].springForce(beadArray[startInd],
                                                     beadArray[endInd]);
            springVelocity = forceToVelocityChange(springForce,beadMassStart,timestep);
            beadArray[startInd].addVelocity(springVelocity);
            beadArray[endInd].addVelocity({-beadMassEndRatio * springVelocity[0],
                                           -beadMassEndRatio * springVelocity[1],
                                           -beadMassEndRatio * springVelocity[2]});
            springOut << currentTime  << ' ' << startInd 
                                      << ' ' << endInd
                                      << ' ' << sqrt(springForce[0]*springForce[0]
                                                    +springForce[1]*springForce[1]
                                                    +springForce[2]*springForce[2])
                                      << ' ' << springArray[i].getId() << '\n'; 
     
        }

        for (int j=0; j<beadArray.getLength(); ++j)
        {
            double beadMass = beadArray[j].getMass();
            double beadRadius = beadArray[j].getRadius();
            currentVelocity = beadArray[j].getVelocity();
            outputPosition = beadArray[j].getPosition();
            beadOut << currentTime << ' ' << outputPosition[0] 
                                   << ' ' << outputPosition[1]
                                   << ' ' << outputPosition[2] 
                                   << ' ' << beadArray[j].getId() 
                                   << ' ' << beadRadius << '\n';

            for (int k=0; k<beadArray.getLength(); ++k)
            {
                if (j != k)
                {
                    overlapForce = beadArray[j].binaryOverlapForce(beadArray[k],
                                                                   overlapStrength);
                    overlapVelocity = forceToVelocityChange(overlapForce,beadMass,timestep);
                    beadArray[j].addVelocity(overlapVelocity);
                    chargeForce = coulombForce(beadArray[j],beadArray[k]);
                    chargeVelocity = forceToVelocityChange(chargeForce,beadMass,timestep);
                    beadArray[j].addVelocity(chargeVelocity);
                }
            }
            
            noiseForce = generateThermalNoise(beadMass,beadRadius,temperature,viscosity);
            noiseVelocity = forceToVelocityChange(noiseForce,beadMass,timestep);
            beadArray[j].addVelocity(noiseVelocity);

            dragForce = stokesDrag(currentVelocity,beadRadius,viscosity);
            dragVelocity = forceToVelocityChange(dragForce,beadMass,timestep);
            beadArray[j].addVelocity(dragVelocity);

            beadArray[j].updatePosition(timestep);
        }

       currentTime += timestep;
    }
    return 0;
}
