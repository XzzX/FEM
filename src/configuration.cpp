#include	"configuration.h"

#include	<fstream>

#include	"util.h"

Configuration	gConfig;

std::ostream &operator << (std::ostream &stream, const Configuration& config){
	stream << "#CONFIGURATION" << std::endl;
	stream << "#filename:                " << config.mFilename << std::endl;
	stream << "#number of refinements:   " << config.mNumberOfRefinements << std::endl;
    return stream;
}

Configuration::Configuration() :
mFilename("data"), mNumberOfRefinements(2)
{}

void	Configuration::ReadCommandLineParameters( unsigned int argc, char **argv ){
    for (unsigned int i=1; i<argc; i++){

        std::string cmd = argv[i];

		if (cmd.compare("-data") == 0){
			i++;
            mFilename = argv[i];
		} else if (cmd.compare("-refine") == 0) {
			i++;
			mNumberOfRefinements = StringTo<unsigned int>(argv[i]);
		}
    }
}

void Configuration::SaveConfiguration(){
    std::ofstream fout(gConfig.mFilename + "_conf.txt");
	fout << gConfig;

	fout << mFilename << std::endl;
	fout << mNumberOfRefinements << std::endl;

    fout.close();
}
