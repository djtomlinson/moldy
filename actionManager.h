void addMoleculeFromFile(Array<Bead>& beadArr, Array<Spring>& springArr, 
                         std::ifstream& beadScript, std::ifstream& springScript)
{
    //check whether counting lines in beads+springs then resizing arrays once before reading lines again is quicker
    
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
        std::cout << "End of simulation script.\n";
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
            addMoleculeFromFile(beadArr, springArr, beadFile, springFile);
        }
    }
}


