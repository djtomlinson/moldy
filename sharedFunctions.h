static std::mt19937 rng(std::random_device{}());

void outputBeadPosition(Bead& bead, double timestep)
 {
     vec3 position {bead.getPosition()};
     std::cout << '(' << position[0] << ',' << position[1] << ',' << position[2] << ")\n";
 }
 
 vec3 getCoM(Array<Bead>& beadArr)
 {
     int noBeads {beadArr.getLength()};
     vec3 CoM {0,0,0};
     for (int b=0; b<noBeads; b++)
     {
         vec3 beadPosition {beadArr[b].getPosition()};
         CoM[0] += beadPosition[0];
         CoM[1] += beadPosition[1];
         CoM[2] += beadPosition[2];
     }
     CoM[0] /= noBeads;
     CoM[1] /= noBeads;
     CoM[2] /= noBeads;
 
     return CoM;
 }

 void getStDevAndCoM(Array<Bead>& beadArr, vec3& StDev, vec3& CoM)
 {
     //subject to catastrophic cancellation
     int noBeads {beadArr.getLength()};
     vec3 Ex {0,0,0}; // E(x)
     vec3 Ex2 {0,0,0}; // E(x^2)
     for (int b=0; b<noBeads; ++b)
     {
         vec3 beadPosition {beadArr[b].getPosition()};
         
         Ex[0] += beadPosition[0];
         Ex[1] += beadPosition[1];
         Ex[2] += beadPosition[2];
     
         Ex2[0] += beadPosition[0]*beadPosition[0];
         Ex2[1] += beadPosition[1]*beadPosition[1];
         Ex2[2] += beadPosition[2]*beadPosition[2];
     }

     double invNoBeads {1.0/noBeads};
     
     Ex[0] *= invNoBeads;
     Ex[1] *= invNoBeads;
     Ex[2] *= invNoBeads;

     Ex2[0] *= invNoBeads;
     Ex2[1] *= invNoBeads;
     Ex2[2] *= invNoBeads;

     CoM[0] = Ex[0];
     CoM[1] = Ex[1];
     CoM[2] = Ex[2];

     StDev[0] = sqrt(Ex2[0] - Ex[0]*Ex[0]);
     StDev[1] = sqrt(Ex2[1] - Ex[1]*Ex[1]);
     StDev[2] = sqrt(Ex2[2] - Ex[2]*Ex[2]);
}

