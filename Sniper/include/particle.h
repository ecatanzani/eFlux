#ifndef PARTICLE_H
#define PARTICLE_H

#include <iostream>

class mc_particle
{
public:
	mc_particle(const long pdg_code);
	mc_particle(){};
	~mc_particle(){};
	void Load(const long pdg_code);
	const bool IsElectron();
	const bool IsProton();

private:
	bool is_electron = false;
	bool is_proton = false;
};

#endif