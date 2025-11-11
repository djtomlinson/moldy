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
#include "sharedFunctions.h" //includes mersenne twister seeding, CoM, StDev
#include "forces.h"
#include "actionManager.h"

int main(){

    double maxTime {}; //try to keep in SI
    double timestep {};
    double viscosity {};
    double temperature {};
    double overlapStrength {};//can be calculated 3*pi*eta*radius/2*timestep for harmonic
    
    double currentTime {0.0};
    
    vec3 zeroForce {0,0,0}; //used to reset forces
    
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
    //to do: copy output to outputBackup folder, numerate


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
           >> mass >> radius >> charge >> species;
        beadArray[lineNo].setPosition(position);
        beadArray[lineNo].setVelocity(velocity);
        beadArray[lineNo].setMass(mass);
        beadArray[lineNo].setRadius(radius);
        beadArray[lineNo].setCharge(charge);
        beadArray[lineNo].setSpecies(species);
        //beadArray[lineNo].print();
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
        /* CHECK SCRIPT, DO ACTIONS */
        checkActions(simScript, beadArray, springArray);
               
        /* DECLARE VARIABLES */
        vec3 overlapForce {};
        vec3 chargeForce {};
        vec3 speciesForce {};
        vec3 springForce {};
        int startInd {};
        int endInd {};

        vec3 noiseForce {};

        vec3 currentVelocity {};
        vec3 dragForce {};

        vec3 outputPosition {};
        
        /* SPRING LOOP */
        for (int i=0; i<springArray.getLength(); ++i)
        {   
            /* SPRING FORCE */
            startInd = springArray[i].getStart();
            endInd = springArray[i].getEnd();
            springForce = springArray[i].springForce(beadArray[startInd],
                                                     beadArray[endInd]);
            beadArray[startInd].addForce(springForce);
            beadArray[endInd].addForce({-springForce[0],-springForce[1],-springForce[2]});
            
            /* WRITING TO SPRING OUTPUT */
            springOut << currentTime  << ' ' << startInd 
                                      << ' ' << endInd
                                      << ' ' << sqrt(springForce[0]*springForce[0]
                                                    +springForce[1]*springForce[1]
                                                    +springForce[2]*springForce[2])
                                      << ' ' << springArray[i].getId() << '\n'; 
        }
        
        /* BEAD LOOP */
        for (int j=0; j<beadArray.getLength(); ++j)
        {
            /* GET BEAD PARAMETERS */
            double beadMass = beadArray[j].getMass();
            double beadRadius = beadArray[j].getRadius();
            
            /* WRITING TO BEAD OUTPUT */
            outputPosition = beadArray[j].getPosition();
            beadOut << currentTime << ' ' << outputPosition[0] 
                                   << ' ' << outputPosition[1]
                                   << ' ' << outputPosition[2] 
                                   << ' ' << beadArray[j].getId() 
                                   << ' ' << beadRadius 
                                   << ' ' << beadArray[j].getSpecies() << '\n';

            /* EVALUATE FORCES FROM OTHER BEADS */
            for (int k=0; k<beadArray.getLength(); ++k)
            {
                if (j != k)
                {
                    /* OVERLAP FORCE: constant force separates beads if overlapping. */
                    overlapForce = beadArray[j].binaryOverlapForce(beadArray[k],
                                                                   overlapStrength);
                    beadArray[j].addForce(overlapForce);

                    /* CHARGE FORCE: explicitly Coulomb (1/r^2) forces evaluated. */
//                    chargeForce = coulombForce(beadArray[j],beadArray[k]);
//                    chargeVelocity = forceToVelocityChange(chargeForce,beadMass,timestep);
//                    beadArray[j].addVelocity(chargeVelocity);
                    
                    /* SPECIES FORCE: Lennard-Jones interaction, mediated by interaction matrix. */
                    speciesForce = interactionForceLJ(beadArray[j],beadArray[k],interactionMatrix);
                    beadArray[j].addForce(speciesForce);
                }
            }
            
            /* NOISE FORCE: thermal force, Brownian motion. */
            noiseForce = generateThermalNoise(beadMass,beadRadius,temperature,viscosity,timestep);
            beadArray[j].addForce(noiseForce);
           
            /* DRAG FORCE: Stokes drag. */
            currentVelocity = beadArray[j].estimateVelocity(timestep); 
            dragForce = stokesDrag(currentVelocity,beadRadius,viscosity);
            beadArray[j].addForce(dragForce);
            
            /* UPDATE BEAD POSITION */
            beadArray[j].updatePositionVerlet(timestep);
            beadArray[j].setForce(zeroForce);
       }

       vec3 stdev {0,0,0};
       vec3 mean {0,0,0};
//       std::cout << mean[0] << ' ' << stdev[0] << '\n';
       getStDevAndCoM(beadArray,stdev,mean);
//       std::cout << mean[0] << ' ' << stdev[0] << '\n';


       currentTime += timestep;
    }
    return 0;
}
