#include "particle.h"

mc_particle::mc_particle(const long pdg_code)
{
	switch (pdg_code)
	{
	case 11:
		is_electron = true;
		is_proton = false;
		break;
	case 2212:
		is_electron = false;
		is_proton = true;
		break;
	default:
		std::cerr << "\n\nWrong particle ID\n\n";
	}
}

void mc_particle::Load(const long pdg_code)
{
	switch (pdg_code)
	{
	case 11:
		is_electron = true;
		is_proton = false;
		break;
	case 2212:
		is_electron = false;
		is_proton = true;
		break;
	default:
		std::cerr << "\n\nWrong particle ID\n\n";
	}
}

const bool mc_particle::IsElectron()
{
	return is_electron;
}

const bool mc_particle::IsProton()
{
	return is_proton;
}
