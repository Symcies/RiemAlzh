#include "RandomVariableRealizations.h"

void RealizationsEvolution(std::map< std::string, std::vector<std::vector<double>> > Realizations, std::ofstream& File )
{
    for(auto it = Realizations.begin(); it != Realizations.end(); ++it)
    {
        File << std::endl << "Realization : " << it->first << ", Nb components : " << it->second[0].size() << std::endl;
        for(int NbComponents = 0; NbComponents < it->second[0].size() ; ++NbComponents)
        {
            for (int TimePoint = 0; TimePoint != it->second.size(); ++TimePoint)
            {
                File << it->second[TimePoint][NbComponents] << ", ";

            }
            File << std::endl;
        }
    }
}