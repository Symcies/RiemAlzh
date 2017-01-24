#include "ReadData.h"


ReadData::Data
ReadData
::OpenFilesMultivariate(std::string FilesPath, int NbMaxOfSubjects = -1)
{
    Data D;
    
    std::ifstream IndivID(FilesPath + "group.csv");
    std::ifstream DataX(FilesPath + "X.csv");
    std::ifstream DataY(FilesPath + "Y.csv");
    
    /// Open the Group file;
    if(IndivID.is_open())
    {
        int i = 0;
        std::string line;
        std::vector< std::pair< LinearAlgebra<ScalarType>::VectorType, double> > IndivData;
        std::pair< LinearAlgebra<ScalarType>::VectorType, double> Observations;
        
        while(getline(IndivID, line))
        {
            int j = std::stoi(line);
            if(j == i)
            {
                IndivData.push_back(Observations);
            }
            else
            {
                D.push_back(IndivData);
                IndivData.clear();
                IndivData.push_back(Observations);
                i = j;
            }
            if(i >= NbMaxOfSubjects && NbMaxOfSubjects != -1)
            {
                break;
            }
        }
        D.push_back(IndivData);
    }
    else { std::cout << "Unable to open indiv id's"; }
    
    /// Open the DATA : timepoints
    if(DataX.is_open())
    {
        std::string line;
        getline(DataX, line);
        for(auto it = D.begin(); it != D.end(); it++)
        {
            for(auto it2 = it->begin(); it2 != it->end(); ++it2)
            {
                it2->second = std::stod(line);
                getline(DataX, line);
            }
        }
    }
    else { std::cout << "Unable to open timepoints"; }
    
    /// Open the DATA : observations
    if(DataY.is_open())
    {
        std::string line;
        getline(DataY, line);
        for(auto it = D.begin(); it != D.end(); ++it)
        {
            for(auto it2 = it->begin(); it2 != it->end(); ++it2)
            {
                LinearAlgebra<ScalarType>::VectorType X(4);
                // TODO : check if everything was parsed
                int i = 0;
                std::stringstream LineStream(line);
                std::string cell;
                while(std::getline(LineStream, cell, ','))
                {
                    X(i) = std::stod(cell);
                    i++;
                }
                it2->first = X;
                getline(DataY, line);
            }
        }
    }
    else { std::cout << "Unable to open the observation values"; }
    
    if( D[0].size() == 0) 
    {
        D.erase(D.begin());
    }
 
    return D;
}


ReadData::MatrixType
ReadData
::OpenKernel(std::string FilePath)
{
    std::ifstream file (FilePath); 
    std::string line;
    int NbLines = 0;
    while(std::getline(file, line))
    {
        ++NbLines;
    }
    file.clear();
    file.seekg(00, std::ios::beg);
    
    
    MatrixType Kernel(NbLines, 1827, 0);

   
    if(file.is_open())
    {
        std::string line;
        int i = 0;
        while(std::getline(file, line))
        {
            int j = 0;
            std::stringstream LineStream(line);
            std:: string cell;
            while(std::getline(LineStream, cell, ','))
            {
                auto a = std::stold(cell);
                if(fabs(a) < 10e-13)
                { 
                    Kernel(i, j) = 0; 
                } 
                else
                {
                   Kernel(i, j) = std::stold(cell);  
                }
                
                ++j;
            }
            ++i;
        }
    }
    else { std::cout << "Unable to open the kernel located in " + FilePath; }
 
    
    return Kernel;
}