#include <math.h>
#include <iostream>
#include <iomanip>
#include <vector>
#include <chrono>
#include "World.h"
#include "PotentialSolver.h"
#include "Species.h"
#include "Output.h"
#include "Source.h"
/* 
*  CC: NOTE THAT FIELD.H IS NOT INCLUDED HERE
*  FIELD.H IS ONLY FOR THE OTHER HEADER FILES 
*/

using namespace std;		//to avoid having to write std::cout
using namespace Const;		//to avoid having to write Const::ME

/*program execution starts here*/
int main(int argc, char *args[])
{
	/*initialize domain*/
	double AR = 0.25; // Aspect ratio (0, 1]
	double lz = 1.0;
	double lx = AR * lz;
	double ly = AR * lz;
	World world(lz * 10 + 1, ly * 10 + 1, lz * 10 + 1);
	world.setExtents({ -lx / 2,-ly / 2,0 }, { lx / 2,ly / 2,lz }); // why are the extents set at a negative number
	world.setTime(1e-7, 400);

	
	/*set objects
	double phi_sphere = 0;		//set default
	if (argc>1)
		phi_sphere = atof(args[1]);	//convert argument to float
	cout << "Sphere potential: " << phi_sphere << " V" << endl;
	world.addSphere({ 0,0,0.15 }, 0.05, phi_sphere);
	*/
    world.addInlet();

	/*set up particle species*/
	vector<Species> species;
	species.push_back(Species("O+", 16*AMU, QE, 1e2, world)); // CC: (name, mass, QE, mpw0, world information)

	/*setup injection sources*/
	vector<ColdBeamSource> sources; 
	sources.push_back(ColdBeamSource(species[0], world, 7000, 1e10));	// ion source, CC: (which species, world, velocity, density)

	/*initialize potential solver and solve initial potential*/
    PotentialSolver solver(world,10000,1e-4); // CC: (world, max_iteration, tolerance)
    solver.setReferenceValues(0,1.5,1e10);
    solver.solve();

    /*obtain initial electric field*/
    solver.computeEF();

    /*main loop*/
	while(world.advanceTime())
    {
	  	/*inject particles*/
    	for (ColdBeamSource &source:sources) // right now there is only 1 source
    		source.sample();

        /*move particles*/
		for (Species &sp:species)
		{
			sp.advance();
			sp.computeNumberDensity();
		}

		/*compute charge density*/
		world.computeChargeDensity(species);

        /*update potential*/
        solver.solve();

        /*obtain electric field*/
        solver.computeEF();

		/*screen and file output*/
        Output::screenOutput(world,species);
        Output::diagOutput(world,species);

		/*periodically write out results*/
        if (world.getTs()%20==0 || world.isLastTimeStep())
			Output::fields(world, species);
    }
	
	/* grab starting time*/
	cout<<"Simulation took "<<world.getWallTime()<<" seconds"<<endl;
	return 0;		//indicate normal exit
}
