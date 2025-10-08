#ifndef FORCES_H
#define FORCES_H

//requires <cmath>, <random>

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

vec3 interactionForce(Bead& bead1, Bead& bead2, auto interactionMatrix)
{
    int s1 {bead1.getSpecies()};
    int s2 {bead2.getSpecies()};
    double interactionStrength {interactionMatrix[s1][s2]};
    double coeff {interactionStrength};
    vec3 p1 {bead1.getPosition()};
    vec3 p2 {bead2.getPosition()};

    double beadSeparation {sqrt((p1[0]-p2[0])*(p1[0]-p2[0])
                               +(p1[1]-p2[1])*(p1[1]-p2[1])
                               +(p1[2]-p2[2])*(p1[2]-p2[2]))};
    if (beadSeparation<2.0*bead1.getRadius())
    {
        return {coeff*(p1[0]-p2[0]),coeff*(p1[1]-p2[1]),coeff*(p1[2]-p2[2])};
    }
    else
    {
        return {0.0, 0.0, 0.0};
    }

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

#endif
