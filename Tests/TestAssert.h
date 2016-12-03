#ifndef TestAssert_h
#define TestAssert_h

typedef double ScalarType;

#include <iostream>
#include <vector>
#include <functional>
#include <cmath>
#include "../LinearAlgebra/LinearAlgebra.h"


class TestAssert {

public:
    
    ////////////////////////////////////////////////////////////////////////////////////////////////////
    // typedef :
    ////////////////////////////////////////////////////////////////////////////////////////////////////

    
    ////////////////////////////////////////////////////////////////////////////////////////////////////
    // Constructor(s) / Destructor :
    ////////////////////////////////////////////////////////////////////////////////////////////////////
    
    /// Act as a constructor by initializing the m_Singleton attribut
    static TestAssert* Instance();

    ////////////////////////////////////////////////////////////////////////////////////////////////////
    // Getter(s) and Setter(s) :
    ////////////////////////////////////////////////////////////////////////////////////////////////////
    
    /// Initialize the activity of the singleton
    static void Init(bool Active);
    

    ////////////////////////////////////////////////////////////////////////////////////////////////////
    // Other method(s) :
    ////////////////////////////////////////////////////////////////////////////////////////////////////
    
    /// Assert that the variables are equal
    template <typename T>
    static void WarningEquality_Object(T a, T b, std::string msg);
    
    
    /// Assert that functions returns the same output
    template <typename F>
    static void WarningEquality_Function(std::function<F> a, std::function<F> b, std::string msg);
    
    /// Assert that the greater_than inequality is fulfilled
    template <typename T>
    static void WarningInequality_GreaterThan(T a, T b, std::string msg); 
    
    /// Assert that the double greater_than inequality is fulfiled
    template <typename T>
            static void WarningInequality_GreaterThan(T a, T b, T c, std::string msg);
    
    
protected:
    ////////////////////////////////////////////////////////////////////////////////////////////////////
    // Constructor(s) / Destructor :
    ////////////////////////////////////////////////////////////////////////////////////////////////////
    
    TestAssert();
    ~TestAssert();
    
    
    ////////////////////////////////////////////////////////////////////////////////////////////////////
    // Method(s) :
    ////////////////////////////////////////////////////////////////////////////////////////////////////

    /// Set the activity
    void SetActive(bool Active);
    
    /// Get the activity
    bool GetActive();

    ////////////////////////////////////////////////////////////////////////////////////////////////////
    // Attribute(s)
    ////////////////////////////////////////////////////////////////////////////////////////////////////
    
    /// Check if the test object is active
    bool m_IsActive = false;
    
    
};


template <typename T>
void 
TestAssert
::WarningEquality_Object(T a, T b, std::string msg) 
{
    if(!Instance()->GetActive()) return;
    
    if(a != b)
    {
        std::cout << msg << std::endl;
    }
}

template <typename F>
void
TestAssert
::WarningEquality_Function(std::function<F> a, std::function<F> b, std::string msg) 
{
    if(!Instance()->GetActive()) return;
    
    double Val = fabs(a() - b());

    if(fabs(Val) > 10e-5)
    {
        std::cout << std::endl << "Value : |" << a() << " - " << b() << "| = " << Val << std::endl;
        //throw std::invalid_argument( msg );
    }
    if(fabs(Val) > 10e-6 && a() < 10e3 && b() < 10e3)
    {
        std::cout << std::endl << "Value : |" << a() << " - " << b() << "| = " << Val << std::endl;
        //throw std::invalid_argument( msg );
    }
}

template <typename T>
void
TestAssert
::WarningInequality_GreaterThan(T a, T b, std::string msg) 
{
    if(!Instance()->GetActive()) return;
    
    if(a <= b)
    {
        std::cout << std::endl << "Should be a: " << a << " > " << b << " :b" << std::endl;
        throw std::invalid_argument( msg );
    }
}


template <typename T>
void
TestAssert
::WarningInequality_GreaterThan(T a, T b, T c, std::string msg) 
{
    if(!Instance()->GetActive()) return;
    
    if(a <= b || b <= c )
    {
        std::cout << std::endl << "Should be :" << a << " > " << b << " > " << c << std::endl;
        throw std::invalid_argument( msg );
    }
    
}





#endif //TestAssert_h
