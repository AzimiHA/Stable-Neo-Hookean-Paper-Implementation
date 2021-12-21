TO build the executable, remove the build folder,
then use:
mkdir build
cd build
cmake ..

to make the build folder. 
Then use Visual Studio to build the solution file.

However, an executable is already provided in the build directory by the name of "a3-finite-elements-3d.exe"
so there is no need to do this.
To run the simulations, simply execute the following commands from inside the build directory:

"./a3-finite-elements-3d" - to run the bunny simulation


"./a3-finite-elements-3d arma" - to run the armadillo simulation

"./a3-finite-elements-3d cube" - to run the cube simulation. Note that this mode has no fixed points and is meant to test
the degeneracies.
	While the simulation is running, you can execute the following:
	Press "J" - Will flatten the cube to a plane
	Press "L" - Will flatten the cube to a line
	Press "P" - Will flatten the cube to a point

"./a3-finite-elements-3d cube2" - to run the cube simulation with 4 fixed points, similar to how the armadillo and bunny meshes work.
In this case there are no options to crush the cube.