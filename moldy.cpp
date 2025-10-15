#include <iostream>
#include <fstream>
#include <sstream>
#include <array>
#include <vector>
#include <random>
#include <cmath>
#include <string>

#include "typeAliases.h"
#include "constants.h"
#include "containers.h"
#include "physicsObjects.h"
#include "forces.h"
#include "actionManager.h"

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
    double overlapStrength {};//can be calculated 3*pi*eta*radius/2*timestep
    
    double currentTime {0.0};
    
    //simulation parameters input here
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
    
    //simulation script loaded here
    std::ifstream simScript {"simulationScript.txt"};
    if (!simScript)
    {
        return 1;
    }
    
    //interaction matrix input happens here
    std::ifstream interactionMatrixFile {"interactionMatrix.txt"};
    std::ifstream interactionMatrixFileCopy {"interactionMatrix.txt"};

    std::string matrixRows {};
    if (!interactionMatrixFile)
    {
        return 1;
        std::cout << "Error reading interaction matrix\n";
    }
    int numberRows {};
    while(std::getline(interactionMatrixFile,matrixRows)) { ++numberRows; }
    
    mat<double> interactionMatrix (numberRows, std::vector<double>(numberRows));

    int rowNumber {};
    while(std::getline(interactionMatrixFileCopy,matrixRows))
    {
        std::istringstream in {matrixRows};
        if (rowNumber!=0) //header with species names is row 0
        {
            
            for (int i=0; i<numberRows-1; i++)
            {
                in >> interactionMatrix[rowNumber-1][i]; //these indices may need swapping
                std::cout << interactionMatrix[rowNumber-1][i] << " currently writing: "
                          << "row " << rowNumber-1 << " column " << i << '\n';
            }
        }
        else
        {
            std::string in1;
            std::string in2;
            in >> in1 >> in2;
            std::cout << in1 << ' ' << in2 << '\n';
        }
        ++rowNumber;
    }

    //check/make output files
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
        int species {};
        in >> position[0] >> position[1] >> position[2]
           >> velocity[0] >> velocity[1] >> velocity[2]
           >> mass >> radius >> charge;
        beadArray[lineNo].setPosition(position);
        beadArray[lineNo].setVelocity(velocity);
        beadArray[lineNo].setMass(mass);
        beadArray[lineNo].setRadius(radius);
        beadArray[lineNo].setCharge(charge);
        beadArray[lineNo].setSpecies(species);
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
                       << ' ' << "Radius" 
                       << ' ' << "Species" << '\n';
    
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
                                   << ' ' << beadRadius 
                                   << ' ' << beadArray[j].getSpecies() << '\n';

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
       
       checkActions(simScript, beadArray, springArray);
       currentTime += timestep;
    }
    return 0;
}
