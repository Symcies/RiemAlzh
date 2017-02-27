#ifndef m_ModelSettings_h
#define m_ModelSettings_h

typedef double ScalarType;

#include <string>
#include <iostream>

#include "tinyxml2.h"

namespace io {

class ModelSettings {

public:

    ////////////////////////////////////////////////////////////////////////////////////////////////
    // typedef :
    ////////////////////////////////////////////////////////////////////////////////////////////////

    ////////////////////////////////////////////////////////////////////////////////////////////////
    // Constructor(s) / Destructor :
    ////////////////////////////////////////////////////////////////////////////////////////////////

    ModelSettings(const char *XMLFile);

    ~ModelSettings();

    ////////////////////////////////////////////////////////////////////////////////////////////////
    /// Encapsulation method(s) :
    ////////////////////////////////////////////////////////////////////////////////////////////////

    std::string GetType() { return m_Type; }

    unsigned int GetManifoldDimension() const { return m_ManifoldDimension; }

    unsigned int GetNumberOfIndependentSources() const { return m_NbIndependentSources; }

    std::string GetInvertKernelPath() const { return m_InvertKernelMatrixPath; }

    std::string GetInterpolationKernelPath() const { return m_InterpolationMatrixPath; }

    ////////////////////////////////////////////////////////////////////////////////////////////////
    /// Other method(s) :
    ////////////////////////////////////////////////////////////////////////////////////////////////

private:

    ////////////////////////////////////////////////////////////////////////////////////////////////
    /// Attribute(s) :
    ////////////////////////////////////////////////////////////////////////////////////////////////

    /// Type of the model
    std::string m_Type;

    /// Name of the output file
    std::string m_OutputFileName;

    /// Dimension of the manifold
    unsigned int m_ManifoldDimension;

    /// Number of sources
    unsigned int m_NbIndependentSources;

    /// Path to the kernel matrix
    std::string m_InvertKernelMatrixPath;

    /// Path to the interpolation matrix
    std::string m_InterpolationMatrixPath;


    ////////////////////////////////////////////////////////////////////////////////////////////////
    /// Model specific methods(s) :
    ////////////////////////////////////////////////////////////////////////////////////////////////

    /// Path to kernel Kxd

    /// Path to kernel invKd

    ////////////////////////////////////////////////////////////////////////////////////////////////
    /// Methods(s) :
    ////////////////////////////////////////////////////////////////////////////////////////////////

    /// Load the fast network model
    void LoadFastNetwork(const char *XMLFile);

    /// Load the meshwork model
    void LoadMeshworkModel(const char *XMLFile);

    /// Load the network model
    void LoadNetworkModel(const char *XMLFile);
};


}    
#endif //m_ModelSettings_h
