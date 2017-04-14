#include "GaussianModel.h"

////////////////////////////////////////////////////////////////////////////////////////////////////
/// Constructor(s) / Destructor :
////////////////////////////////////////////////////////////////////////////////////////////////////

GaussianModel::GaussianModel(io::ModelSettings &model_settings) {
  
}

GaussianModel::~GaussianModel() {
  
}


////////////////////////////////////////////////////////////////////////////////////////////////////
/// Other method(s) :
////////////////////////////////////////////////////////////////////////////////////////////////////

void GaussianModel::Initialize(const Observations &obs) {
  
}

void GaussianModel::InitializeValidationDataParameters(const io::SimulatedDataSettings &data_settings,
                                                       const io::ModelSettings &model_settings) {
  
}

void GaussianModel::UpdateModel(const Realizations &reals,
                                const MiniBlock &block_info,
                                const std::vector<std::string,
                                                  std::allocator<std::string>> names) {
  
}

AbstractModel::SufficientStatisticsVector GaussianModel::GetSufficientStatistics(const Realizations &reals,
                                                                                 const Observations &obs) {
  
}

void GaussianModel::UpdateRandomVariables(const SufficientStatisticsVector &stoch_sufficient_stats) {
  
}

Observations GaussianModel::SimulateData(io::SimulatedDataSettings &data_settings) {
  
}

std::vector<AbstractModel::MiniBlock> GaussianModel::GetSamplerBlocks() const {
  
}

////////////////////////////////////////////////////////////////////////////////////////////////////
/// Log-likelihood related method(s) :
////////////////////////////////////////////////////////////////////////////////////////////////////


AbstractModel::VectorType GaussianModel::ComputeLogLikelihood(const Observations &obs,
                                                              const MiniBlock &block_info) {
  
}

ScalarType GaussianModel::ComputeIndividualLogLikelihood(const IndividualObservations &obs,
                                                         const int subjects_tot_num_) {
  
}

ScalarType  GaussianModel::GetPreviousLogLikelihood(const MiniBlock &block_info) {
  
}

void GaussianModel::SetPreviousLogLikelihood(VectorType &log_likelihood,
                                             const MiniBlock &block_info) {
  
}


////////////////////////////////////////////////////////////////////////////////////////////////////
/// Outputs
////////////////////////////////////////////////////////////////////////////////////////////////////


void GaussianModel::DisplayOutputs(const Realizations &reals) {
  
}

void GaussianModel::SaveData(unsigned int IterationNumber, const Realizations &reals) {
  
}