#include <boost/math/distributions.hpp>
#include "genfile/SingleSNPGenotypeProbabilities.hpp"
#include "genfile/Error.hpp"
#include "integration/NewtonRaphson.hpp"
#include "integration/Derivative.hpp"
#include "snptest/case_control/FrequentistTest.hpp"
#include "snptest/case_control/NullModelLogLikelihood.hpp"
#include "snptest/case_control/AlternativeModelLogLikelihood.hpp"
#include "snptest/FinitelySupportedFunctionSet.hpp"

//#define DEBUG_TEST 1

namespace snptest {
	namespace case_control {
		FrequentistTest::FrequentistTest(
			double sample_inclusion_threshhold,
			bool mimic_snptest
		):
			m_sample_inclusion_threshold( sample_inclusion_threshhold ),
			m_mimic_snptest( mimic_snptest ),
			m_chi_squared( 1.0 )
		{}
		
		FrequentistTest::Results FrequentistTest::test(
			genfile::SNPIdentifyingData const& snp,
			Vector const& phenotypes,
			FinitelySupportedFunctionSet const& genotypes,
			Matrix const& covariates,
			std::vector< std::size_t > const& indices_of_samples_to_exclude
		) const {
			if( covariates.rows() != 0 || covariates.cols() != 0 ) {
				throw genfile::BadArgumentError( "CaseControlFrequentistTest::test()", "covariates (nonempty)" ) ;
			}
			
			assert( phenotypes.size() == genotypes.get_values().rows() ) ;
			
			using integration::derivative ;
			snptest::case_control::NullModelLogLikelihood null_loglikelihood(
				phenotypes,
				genotypes,
				covariates,
				indices_of_samples_to_exclude
			) ;

			integration::Derivative< snptest::case_control::NullModelLogLikelihood > null_loglikelihood_derivative = derivative( null_loglikelihood ) ;
			Vector null_parameters = integration::find_root_by_newton_raphson(
				null_loglikelihood_derivative,
				Vector::Zero( 1 ),
				0.00001
			) ;

#if DEBUG_TEST
			std::cerr << "AssociationTester:        ################ testing SNP " << snp.get_rsid() << " ################\n" ;
			std::cerr << "AssociationTester:        null: MLE is " << null_parameters << ".\n" ;
			std::cerr << "AssociationTester:        null: loglikelihood is " << null_loglikelihood.get_value_of_function() << ".\n" ;
			std::cerr << "AssociationTester: alternative: Finding MLE...\n" ;
#endif
			snptest::case_control::AlternativeModelLogLikelihood alternative_loglikelihood(
				phenotypes,
				genotypes,
				covariates,
				indices_of_samples_to_exclude
			) ;
#if DEBUG_TEST
		//alternative_loglikelihood.evaluate_at( Vector::Zero( 2 ) ) ;
		//std::cerr << "AssociationTester: alternative: Initial derivative of LL is " << alternative_loglikelihood.get_value_of_first_derivative() << ".\n" ;
#endif
			integration::Derivative< snptest::case_control::AlternativeModelLogLikelihood > alternative_loglikelihood_derivative = derivative( alternative_loglikelihood ) ;
			Vector alternative_parameters = integration::find_root_by_newton_raphson(
				alternative_loglikelihood_derivative,
				Vector::Zero( 2 ),
				0.00001
			) ;
#if DEBUG_TEST
			std::cerr << "AssociationTester: alternative: MLE is " << alternative_parameters << ".\n" ;
			std::cerr << "AssociationTester: alternative: loglikelihood is " << alternative_loglikelihood.get_value_of_function() << ".\n" ;
			std::cerr << "AssociationTester: -2LR = " << -2.0 * ( null_loglikelihood.get_value_of_function() - alternative_loglikelihood.get_value_of_function() ) << ".\n" ;
#endif
			Results result ;
			result.test_statistic = result.p_value = result.beta = result.standard_error = std::numeric_limits< double >::quiet_NaN() ;
			result.null_loglikelihood = null_loglikelihood.get_value_of_function() ;
			result.alternative_loglikelihood = alternative_loglikelihood.get_value_of_function() ;
			result.test_statistic = -2.0 * ( null_loglikelihood.get_value_of_function() - alternative_loglikelihood.get_value_of_function() ) ;
			if( result.test_statistic > 0.0 ) {
				result.p_value = boost::math::cdf(
					boost::math::complement(
						m_chi_squared,
						-2.0 * ( null_loglikelihood.get_value_of_function() - alternative_loglikelihood.get_value_of_function() )
					)
				) ;
			}
			result.beta = alternative_parameters( 1 ) ;
			result.variance_covariance = (-alternative_loglikelihood.get_value_of_second_derivative()).inverse() ;
			result.standard_error = std::sqrt( result.variance_covariance( 1, 1 ) ) ;
			return result ;
		}
	}
}

