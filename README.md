moldy is a <ins>**mol**</ins>ecular <ins>**dy**</ins>namics simulation, using bead and spring dynamics and physically derived forces to simulate polymers etc.

- The simulation is run from compiled C++ code.

- The simulation can be initialised via python functions provided in [moldy.py](moldy.py)

- The aim is to provide a user friendly front-end utilising python scripting whilst retaining the speed of the C++ program.

- Examples of python function usage included in [moldy.py](moldy.py)

- Example graphical outputs included in [Examples](examples)

- Compiled with c++20, moldy executable pre-compiled with `clang++ -std=c++20 -o moldy moldy.cpp`
