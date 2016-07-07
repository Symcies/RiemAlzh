#include <iostream>

#include <GaussianRandomVariable.h>
#include <LongitudinalModel.h>
#include <MultivariateLogisticManifold.h>
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
    // Initialisation Data:
    ////////////////////////////////////////////////////////////////////////////////////////////////////

    pair <double, vector<double>> Obs1S1 (70.2, vector<double> (5,10));
    pair <double, vector<double>> Obs2S1 (84.7, vector<double> (5,9));
    vector<std::pair<double, std::vector<double>>> S1;
    S1.push_back(Obs1S1);
    S1.push_back(Obs2S1);

    pair <double, vector<double>> Obs1S2 (81.1, vector<double> (5,8));
    pair <double, vector<double>> Obs2S2 (83.6, vector<double> (5,6));
    pair <double, vector<double>> Obs3S2 (87.5, vector<double> (5,3));
    vector<pair<double, vector<double>>> S2;
    S2.push_back(Obs1S2);
    S2.push_back(Obs2S2);
    S2.push_back(Obs3S2);

    Data *D = new Data;
    D->push_back(S1);
    D->push_back(S2);


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


    A->ComputeMCMCSAEM(300);

    std::vector<double> Parameters = A->GetParameters();
    cout << Parameters.size();
    cout << Parameters[0];
    cout << Parameters[1];
    cout << Parameters[2];
    cout << Parameters[3];



    cout << "END";







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