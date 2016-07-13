#include <iostream> // cout
#include <limits>

#include <GaussianRandomVariable.h>
#include <LongitudinalModel.h>
#include <Algorithm.h>
#include <algorithm>
#include "Samplers/HastingMetropolisWithinGibbs.h"

using namespace std;

std::vector< std::vector< double >>
SimulateTimePoint(int NumberOfSubjects, int MinObs, int MaxObs)
{
    random_device RD;
    mt19937 RNG(RD());
    uniform_int_distribution<int> Uni(MinObs, MaxObs);
    normal_distribution<double> Normal(70, 6);

    vector< vector< double >> TimePoint;
    for(int i = 0; i<NumberOfSubjects; ++i)
    {
        vector<double> SubjectTimePoint;
        for(int j = 0; j < Uni(RNG) ; ++j)
        {
            SubjectTimePoint.push_back( Normal(RNG));
        }
        sort(SubjectTimePoint.begin(), SubjectTimePoint.end());
        TimePoint.push_back(SubjectTimePoint);
    }

    return TimePoint;
}


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

    std::vector<std::vector<double>> TimePoint = SimulateTimePoint(4, 2, 6);

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


    A->ComputeMCMCSAEM(5);

    std::vector<double> Parameters = A->GetParameters();
    std::cout << "P0_Mean : " << Parameters[0] << std::endl;
    std::cout << "T0_Mean : " << Parameters[1] << std::endl;
    std::cout << "V0_Mean : " << Parameters[2] << std::endl;
    std::cout << "Ksi_0_Mean : " << Parameters[3] << std::endl;
    std::cout << "Tau_0_Mean : " << Parameters[4] << std::endl;
    std::cout << "Uncertainty_Var : " << Parameters[5] << std::endl;


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
