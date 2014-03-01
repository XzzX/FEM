#ifndef CONFIGURATION_H
#define CONFIGURATION_H

#include	<iostream>
#include	<sstream>
#include	<iomanip>
#include	<string>
#include    <stdexcept>

class Configuration{
	public:
		std::string		mFilename;
		unsigned int	mNumberOfRefinements;

		///load basic configuration
		Configuration();
		///read out the command line parameters
		void	ReadCommandLineParameters( unsigned int argc, char **argv );
		void	SaveConfiguration();
};

extern	Configuration	gConfig;

///overloaded << operator to show easily the configuration on the screen
std::ostream &operator << (std::ostream &stream, const Configuration& config);

#endif
