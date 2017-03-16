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
    inline       VectorType& at(const int key)       { return realizations_.at(key); }
    inline const VectorType& at(const int key) const { return realizations_.at(key); }

    /// Less efficient way to get the vector of realizations
    inline       VectorType& at(const std::string name)       { return realizations_.at(string_to_int_key_.at(name)); };
    inline const VectorType& at(const std::string name) const { return realizations_.at(string_to_int_key_.at(name)); };

    /// Efficient way to get a particular realization
    inline       ScalarType& at(const int key, const unsigned int number)       { return realizations_.at(key)(number); }
    inline const ScalarType& at(const int key, const unsigned int number) const { return realizations_.at(key)(number); }

    /// Less efficient way to get a particular realization
    inline       ScalarType& at(const std::string name, const unsigned int number)       { return realizations_.at(string_to_int_key_.at(name))(number); }
    inline const ScalarType& at(const std::string name, const unsigned int number) const { return realizations_.at(string_to_int_key_.at(name))(number); }


    ////////////////////////////////////////////////////////////////////////////////////////////////////
    /// Other method(s) :
    ////////////////////////////////////////////////////////////////////////////////////////////////////

    void AddRealizations(std::string name, int key, VectorType& values);

    inline int ReverseNameToKey(std::string name) const { return string_to_int_key_.at(name); }

    inline std::string ReverseKeyToName(int key) const { return int_to_string_key_.at(key); }

    inline IntVectorHash::iterator begin() { return realizations_.begin(); }
    inline IntVectorHash::iterator   end() { return realizations_.end(); }

    inline IntVectorHash::const_iterator begin() const { return realizations_.begin(); }
    inline IntVectorHash::const_iterator   end() const { return realizations_.end(); }

    inline VectorType::iterator begin(int key) { return realizations_.at(key).begin(); }
    inline VectorType::iterator   end(int key) { return realizations_.at(key).end(); }

    inline VectorType::const_iterator begin(int key) const { return realizations_.at(key).begin(); }
    inline VectorType::const_iterator   end(int key) const { return realizations_.at(key).end(); }

    inline VectorType::iterator begin(std::string name)  { return realizations_.at(string_to_int_key_.at(name)).begin(); }
    inline VectorType::iterator   end(std::string name) { return realizations_.at(string_to_int_key_.at(name)).end(); }

    inline VectorType::const_iterator begin(std::string name) const { return realizations_.at(string_to_int_key_.at(name)).begin(); }
    inline VectorType::const_iterator   end(std::string name) const { return realizations_.at(string_to_int_key_.at(name)).end(); }




private:
    ////////////////////////////////////////////////////////////////////////////////////////////////////
    /// Method(s) :
    ////////////////////////////////////////////////////////////////////////////////////////////////////



    ////////////////////////////////////////////////////////////////////////////////////////////////////
    /// Attribute(s) :
    ////////////////////////////////////////////////////////////////////////////////////////////////////

    /// Hash table of all the attributes
    IntVectorHash realizations_;

    /// Convert the names into the int key
    StringIntHash string_to_int_key_;

    /// Convert the int into their names
    IntStringHash int_to_string_key_;
};


#endif //_Realizations_h
