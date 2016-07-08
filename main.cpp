#include <iostream>

#include <GaussianRandomVariable.h>
#include <LongitudinalModel.h>
#include <Algorithm.h>
#include "Samplers/HastingMetropolisWithinGibbs.h"

using namespace std;

int main() {
    ////////////////////////////////////////////////////////////////////////////////////////////////////
    // Initialization Parameters
    ////////////////////////////////////////////////////////////////////////////////////////////////////
    int DimensionNumber = 5;
    int NumberOfIndependentComponents = 3;


    ////////////////////////////////////////////////////////////////////////////////////////////////////
    // Initialisation Manifold:
    ////////////////////////////////////////////////////////////////////////////////////////////////////
    auto Manifold = std::make_shared<MultivariateLogisticManifold>(DimensionNumber, NumberOfIndependentComponents);


    ////////////////////////////////////////////////////////////////////////////////////////////////////
    // Initialisation Model:
    ////////////////////////////////////////////////////////////////////////////////////////////////////
    auto LM = std::make_shared<LongitudinalModel>();
    LM->SetManifold(Manifold);


    ////////////////////////////////////////////////////////////////////////////////////////////////////
    // Initialisation Data:
    ////////////////////////////////////////////////////////////////////////////////////////////////////

    std::vector<std::vector<double>> TimePoint;

    // First subject
    std::vector<double> T1;
    T1.push_back(70.2);
    T1.push_back(72.3);
    T1.push_back(77.8);
    T1.push_back(80.2);

    // Second subject
    std::vector<double> T2;
    T2.push_back(55.2);
    T2.push_back(59.8);

    // Second subject
    std::vector<double> T3;
    T3.push_back(60.9);
    T3.push_back(65.3);
    T3.push_back(67.8);
    T3.push_back(69.0);
    T3.push_back(71.0);
    T3.push_back(75.8);

    TimePoint.push_back(T1);
    TimePoint.push_back(T2);
    TimePoint.push_back(T3);


    Data *D =  new Data;
    *D = LM->SimulateSpecificData(TimePoint);

    ////////////////////////////////////////////////////////////////////////////////////////////////////
    // Initialisation Sampler:
    ////////////////////////////////////////////////////////////////////////////////////////////////////
    auto Sampler = std::make_shared<HastingMetropolisWithinGibbs>();


    ////////////////////////////////////////////////////////////////////////////////////////////////////
    // Initialisation Algorithm:
    ////////////////////////////////////////////////////////////////////////////////////////////////////
    Algorithm *A = new Algorithm();
    A->SetModel(LM);
    A->SetSampler(Sampler);
    A->SetData(D);
    A->Initialize();



    ////////////////////////////////////////////////////////////////////////////////////////////////////
    // Initialisation Algorithm:
    ////////////////////////////////////////////////////////////////////////////////////////////////////


    A->ComputeMCMCSAEM(500);

    std::vector<double> Parameters = A->GetParameters();
    cout << Parameters.size() << endl;
    cout << Parameters[0] << endl;
    cout << Parameters[1] << endl;
    cout << Parameters[2] << endl;
    cout << Parameters[3] << endl;


    ////////////////////////////////////////////////////////////////////////////////////////////////////
    // Deletes
    ////////////////////////////////////////////////////////////////////////////////////////////////////
    delete D;








    //
    //                       _oo0oo_
    //                      o8888888o
    //                      88" . "88
    //                      (| -_- |)
    //                      0\  =  /0
    //                    ___/`---'\___
    //                  .' \\|     |// '.
    //                 / \\|||  :  |||// \
    //                / _||||| -:- |||||- \
    //               |   | \\\  -  /// |   |
    //               | \_|  ''\---/''  |_/ |
    //               \  .-\__  '-'  ___/-. /
    //             ___'. .'  /--.--\  `. .'___
    //          ."" '<  `.___\_<|>_/___.' >' "".
    //         | | :  `- \`.;`\ _ /`;.`/ - ` : | |
    //         \  \ `_.   \_ __\ /__ _/   .-` /  /
    //     =====`-.____`.___ \_____/___.-`___.-'=====
    //                       `=---='
    //
    //
    //     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    //
    //               佛祖保佑         永无BUG
    //
    return 0;
}