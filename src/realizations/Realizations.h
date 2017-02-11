#ifndef _Realizations_h
#define _Realizations_h

typedef double ScalarType;

#include <unordered_map>

#include "LinearAlgebra.h"


class Realizations {
    
    
public:
    
    ////////////////////////////////////////////////////////////////////////////////////////////////////
    /// typedef :
    ////////////////////////////////////////////////////////////////////////////////////////////////////
    
    typedef typename LinearAlgebra<ScalarType>::VectorType VectorType;
    typedef typename std::unordered_map<int, VectorType> IntVectorHash;
    typedef typename std::unordered_map<std::string, int> StringIntHash;
    typedef typename std::unordered_map<int, std::string> IntStringHash;
    
    
    ////////////////////////////////////////////////////////////////////////////////////////////////////
    /// Constructor(s) / Destructor :
    ////////////////////////////////////////////////////////////////////////////////////////////////////
    
    Realizations();
    ~Realizations();
    
    ////////////////////////////////////////////////////////////////////////////////////////////////////
    /// Encapsulation method(s) :
    ////////////////////////////////////////////////////////////////////////////////////////////////////

    /// Efficient way to get the vector of realizations
    inline       VectorType& at(const int Key)       { return m_Realizations.at(Key); }
    inline const VectorType& at(const int Key) const { return m_Realizations.at(Key); }
    
    /// Less efficient way to get the vector of realizations
    inline       VectorType& at(const std::string Name)       { return m_Realizations.at(m_StringToIntKey.at(Name)); };
    inline const VectorType& at(const std::string Name) const { return m_Realizations.at(m_StringToIntKey.at(Name)); };
    
    /// Efficient way to get a particular realization
    inline       ScalarType& at(const int Key, const unsigned int Number)       { return m_Realizations.at(Key)(Number); }
    inline const ScalarType& at(const int Key, const unsigned int Number) const { return m_Realizations.at(Key)(Number); }
    
    /// Less efficient way to get a particular realization
    inline       ScalarType& at(const std::string Name, const unsigned int Number)       { return m_Realizations.at(m_StringToIntKey.at(Name))(Number); }
    inline const ScalarType& at(const std::string Name, const unsigned int Number) const { return m_Realizations.at(m_StringToIntKey.at(Name))(Number); }
    
    /*
        /// Efficient way to get the vector of realizations
    inline VectorType at(int Key) const { return m_Realizations.at(Key); }
    
    /// Less efficient way to get the vector of realizations
    inline VectorType at(std::string Name) const { return m_Realizations.at(m_StringToIntKey.at(Name)); };
    
    /// Efficient way to get a particular realization
    inline ScalarType at(int Key, unsigned int Number) const { return m_Realizations.at(Key)(Number); }
    
    /// Less efficient way to get a particular realization
    inline ScalarType at(std::string Name, unsigned int Number) const { return m_Realizations.at(m_StringToIntKey.at(Name))(Number); }
    
    /// Efficient way to set a particular realization
    inline void set(int Key, unsigned int Number, ScalarType Value) { m_Realizations.at(Key)(Number) = Value; }
    
    /// Less efficient way to set a particular realization
    inline void set(std::string Name, unsigned int Number, ScalarType Value) { m_Realizations.at(m_StringToIntKey.at(Name))(Number) = Value; }
    */
    ////////////////////////////////////////////////////////////////////////////////////////////////////
    /// Other method(s) :
    ////////////////////////////////////////////////////////////////////////////////////////////////////
        
    void AddRealizations(std::string Name, int Key, VectorType& Values);
    
    




private:
    ////////////////////////////////////////////////////////////////////////////////////////////////////
    /// Method(s) :
    ////////////////////////////////////////////////////////////////////////////////////////////////////
    

    
    ////////////////////////////////////////////////////////////////////////////////////////////////////
    /// Attribute(s) :
    ////////////////////////////////////////////////////////////////////////////////////////////////////
    
    /// Hash table of all the attributes
    IntVectorHash m_Realizations;
    
    /// Convert the names into the int key
    StringIntHash m_StringToIntKey;
    
    /// Convert the int into their names
    IntStringHash m_IntToStringKey;
};


#endif //_Realizations_h
