#ifndef COMPONENTS_RELATEDNESS_COMPONENT_PCA_COMPUTER_HPP
#define COMPONENTS_RELATEDNESS_COMPONENT_PCA_COMPUTER_HPP

#include <Eigen/Core>
#include "genfile/SingleSNPGenotypeProbabilities.hpp"
#include "genfile/CohortIndividualSource.hpp"
#include "appcontext/OptionProcessor.hpp"
#include "appcontext/UIContext.hpp"
#include "components/RelatednessComponent/KinshipCoefficientManager.hpp"

struct PCAComputer: public KinshipCoefficientManager
{
	typedef std::auto_ptr< PCAComputer > UniquePtr ;
	~PCAComputer() throw() {}
	PCAComputer(
		appcontext::OptionProcessor const& options,
		genfile::CohortIndividualSource const& samples,
		worker::Worker* worker,
		appcontext::UIContext& ui_context
	) ;

	void compute_PCA() ;

	void begin_processing_snps( std::size_t number_of_samples ) ;
	void processed_snp( genfile::SNPIdentifyingData const&, genfile::VariantDataReader& ) ;
	void end_processing_snps() ;
	
	std::string get_PCA_prefix() const ;
	
private:
	appcontext::OptionProcessor const& m_options ;
	appcontext::UIContext& m_ui_context ;
	genfile::CohortIndividualSource const& m_samples ;
	std::string m_filename ;
	std::size_t m_number_of_samples ;
	std::size_t m_number_of_snps ;
	std::size_t m_number_of_snps_processed ;
	Eigen::MatrixXd m_kinship_matrix ;
	Eigen::MatrixXd m_kinship_eigendecomposition ;
	Eigen::VectorXd m_PCA_eigenvalues ;
	Eigen::MatrixXd m_PCA_components ;

	std::size_t m_number_of_PCAs_to_compute ;
	double m_threshhold ;
	genfile::SingleSNPGenotypeProbabilities m_genotype_probabilities ;
	Eigen::VectorXd m_genotype_calls ;
	Eigen::VectorXd m_non_missingness ;
	Eigen::VectorXd m_loading_vectors ;
	
	void load_matrix( std::string const& filename, Eigen::MatrixXd* matrix ) const ;
} ;

#endif