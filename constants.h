#ifndef CONSTANTS_H
#define CONSTANTS_H

/* Constants, k used for brevity.
 * Symbols follow standard physics conventions.
 * Greek letters are written explicity, where capitalisation
 * indicates use of the capital greek letter.
 * an underscore (_) indicates subscripts.
 * Units are given in square brackets.
 * Source is NIST: https://pml.nist.gov/cuu/pdf/wall_2022.pdf
*/

namespace k
{
    inline constexpr double pi {3.14159}; //[] circle constant
    inline constexpr double k_B {1.380649e-23}; //[JK^-1] Boltzmann constant
    inline constexpr double epsilon_0 {8.8541878188e-12}; //[Fm^-1]  vacuum electric permittivity
    inline constexpr double k_e {1/(4*pi*epsilon_0)}; //[mF^-1] Coulomb constant
}

#endif
