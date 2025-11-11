void addMoleculeFromFile(Array<Bead>& beadArr, Array<Spring>& springArr, 
                         std::ifstream& beadScript, std::ifstream& springScript)
{
    //to do: check whether counting lines in beads+springs then resizing arrays once before reading lines again is quicker
    
    std::string currBead {};
    std::string currSpring {};
    int springIndOffset {beadArr[beadArr.getLength()-1].getId()+1}; //spring ind from files offset

    while (std::getline(beadScript,currBead))
    {
        //std::cout << "Adding bead... \n";
        vec3 position {};
        vec3 velocity {};
        double mass {};
        double radius {};
        double charge {};
        int species {};

        std::istringstream in {currBead};
        in >> position[0] >> position[1] >> position[2]
           >> velocity[0] >> velocity[1] >> velocity[2]
           >> mass >> radius >> charge >> species;
        beadArr.extend1();
        int arrLength {beadArr.getLength()};
        beadArr[arrLength-1].setPosition(position);
        beadArr[arrLength-1].setVelocity(velocity);
        beadArr[arrLength-1].setMass(mass);
        beadArr[arrLength-1].setRadius(radius);
        beadArr[arrLength-1].setCharge(charge);
        beadArr[arrLength-1].setSpecies(species);
        //std::cout << "New array length " << arrLength << '\n';
    }
    while (std::getline(springScript, currSpring))
    {
        int start {};
        int end {};
        double naturalLength {};
        double springConstant {};

        std::istringstream in {currSpring};
        in >> start >> end >> naturalLength >> springConstant;
        springArr.extend1();
        int arrLength {springArr.getLength()};
        springArr[arrLength-1].setStart(start+springIndOffset);
        springArr[arrLength-1].setEnd(end+springIndOffset);
        springArr[arrLength-1].setNaturalLength(naturalLength);
        springArr[arrLength-1].setSpringConstant(springConstant);
    
    }

}

void addMoleculeNearCoM(Array<Bead>& beadArr, Array<Spring>& springArr, 
                         std::ifstream& beadScript, std::ifstream& springScript)
{
    //to do: check whether counting lines in beads+springs then resizing arrays once before reading lines again is quicker
    
    std::string currBead {};
    std::string currSpring {};
    int springIndOffset {beadArr[beadArr.getLength()-1].getId()+1}; //spring ind from files offset

    static std::uniform_real_distribution<> thetaGen(0, k::pi);
    static std::uniform_real_distribution<> phiGen(0, 2*k::pi);

    vec3 beadCoM {};
    vec3 beadStDev {};
    getStDevAndCoM(beadArr, beadStDev, beadCoM);
    vec3 safeDistance {3*beadStDev[0], 3*beadStDev[1], 3*beadStDev[2]};
    safeDistance = {2e-8,2e-8,2e-8};

    double theta = thetaGen(rng);
    double phi = phiGen(rng);
//    std::cout << "theta=" << theta << " phi=" << phi << '\n';
    double xAdd {safeDistance[0]*sin(theta)*cos(phi)+beadCoM[0]};
    double yAdd {safeDistance[1]*sin(theta)*sin(phi)+beadCoM[1]};
    double zAdd {safeDistance[2]*cos(theta)+beadCoM[2]};

    while (std::getline(beadScript,currBead))
    {
        //std::cout << "Adding bead... \n";
        vec3 position {};
        vec3 velocity {};
        double mass {};
        double radius {};
        double charge {};
        int species {};

        std::istringstream in {currBead};
        in >> position[0] >> position[1] >> position[2]
           >> velocity[0] >> velocity[1] >> velocity[2]
           >> mass >> radius >> charge >> species;
        beadArr.extend1();

        position = {position[0]+xAdd,position[1]+yAdd,position[2]+zAdd};
        //std::cout << "Added bead of species: " << species << '\n';
        int arrLength {beadArr.getLength()};
        beadArr[arrLength-1].setPosition(position);
        beadArr[arrLength-1].setVelocity(velocity);
        beadArr[arrLength-1].setMass(mass);
        beadArr[arrLength-1].setRadius(radius);
        beadArr[arrLength-1].setCharge(charge);
        beadArr[arrLength-1].setSpecies(species);
        //std::cout << "New array length " << arrLength << '\n';
    }
    while (std::getline(springScript, currSpring))
    {
        int start {};
        int end {};
        double naturalLength {};
        double springConstant {};

        std::istringstream in {currSpring};
        in >> start >> end >> naturalLength >> springConstant;
        springArr.extend1();
        int arrLength {springArr.getLength()};
        springArr[arrLength-1].setStart(start+springIndOffset);
        springArr[arrLength-1].setEnd(end+springIndOffset);
        springArr[arrLength-1].setNaturalLength(naturalLength);
        springArr[arrLength-1].setSpringConstant(springConstant);
    
    }

}

void checkActions(std::ifstream& simScript, Array<Bead>& beadArr, Array<Spring>& springArr)
{
    std::string currentAction;

    double actionTime {};
    std::string actionType {};
    int actionNumber {};
    std::string actionInfo {};
    std::string actionInfo2 {};

    if (std::getline(simScript, currentAction))
    {
        std::istringstream in {currentAction};
        in >> actionTime >> actionType >> actionNumber >> actionInfo >> actionInfo2;
        /*std::cout << "Completing action " << actionType << " at time "
                  << actionTime << '\n';*/
    }

    else
    {
        
        //std::cout << "End of simulation script.\n";
    }

    /*if ((actionTime-currTime)!=0)
    {
        std::cout << "Simulation script encountered time mismatch. "
                  << actionTime << " vs " << currTime << '\n';
        std::cout << (actionTime-currTime) << '\n' ;
    }*/
    
    if (actionType=="add")
    {
        for (int i=0; i<actionNumber; i++)
        {
            std::ifstream beadFile {actionInfo};
            std::ifstream springFile {actionInfo2};
            
            // plus checks to see if opened correctly
            //addMoleculeFromFile(beadArr, springArr, beadFile, springFile);
            addMoleculeNearCoM(beadArr, springArr, beadFile, springFile); 
        }
    }
}


