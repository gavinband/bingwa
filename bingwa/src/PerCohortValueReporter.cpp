
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#include <memory>
#include <string>
#include <vector>
#include <boost/function.hpp>
#include <boost/foreach.hpp>
#include <Eigen/Core>
#include "genfile/VariantIdentifyingData.hpp"
#include "genfile/VariantEntry.hpp"
#include "genfile/string_utils.hpp"
#include "appcontext/OptionProcessor.hpp"
#include "bingwa/BingwaComputation.hpp"
#include "bingwa/PerCohortValueReporter.hpp"

namespace bingwa {
	PerCohortValueReporter::PerCohortValueReporter( std::vector< std::string > const& cohort_names ):
		m_cohort_names( cohort_names )
	{}
	
	void PerCohortValueReporter::add_variable( std::string const& variable ) {
		m_extra_variables.insert( variable ) ;
	}
	
	void PerCohortValueReporter::set_effect_parameter_names( EffectParameterNamePack const& names ) {
		m_effect_parameter_names = names ;
	}
	
	void PerCohortValueReporter::get_variables( boost::function< void ( std::string, std::string ) > callback ) const {
		using genfile::string_utils::to_string ;
		
		std::size_t const N = m_cohort_names.size() ;
		for( std::size_t i = 0; i < N; ++i ) {
			std::string prefix = m_cohort_names[ i ] + ":" ;

			for( std::size_t i = 0; i < m_effect_parameter_names.size(); ++i ) {
//				callback( prefix + "beta_" + to_string( i+1 ) ) ;
//				callback( prefix + "se_" + to_string( i+1 ) ) ;
				callback( prefix + m_effect_parameter_names.parameter_name(i), "FLOAT" ) ;
				callback( prefix + m_effect_parameter_names.se_name(i), "FLOAT" ) ;
			}
			for( std::size_t i = 0; i < m_effect_parameter_names.size(); ++i ) {
				for( std::size_t j = i+1; j < m_effect_parameter_names.size(); ++j ) {
					callback( prefix + m_effect_parameter_names.covariance_name(i,j), "FLOAT" ) ;
				}
			}
			callback( prefix + "pvalue", "FLOAT" ) ;
		}
	}
	
	void PerCohortValueReporter::operator()(
		VariantIdentifyingData const&,
		DataGetter const& data_getter,
		ResultCallback callback
	) {
		std::size_t const N = m_cohort_names.size() ;
		for( std::size_t i = 0; i < N; ++i ) {
			if( data_getter.is_non_missing( i ) ) {
				Eigen::VectorXd betas ;
				Eigen::VectorXd ses ;
				Eigen::VectorXd covariance ;
				double pvalue ;
				data_getter.get_betas( i, &betas ) ;
				data_getter.get_ses( i, &ses ) ;
				data_getter.get_covariance_upper_triangle( i, &covariance ) ;
				data_getter.get_pvalue( i, &pvalue ) ;

				using genfile::string_utils::to_string ;
				std::string prefix = m_cohort_names[ i ] + ":" ;
				assert( betas.size() == ses.size() ) ;
				assert( covariance.size() == ( betas.size() - 1 ) * betas.size() / 2 ) ;
				for( int j = 0; j < betas.size(); ++j ) {
//					callback( prefix + "beta_" + to_string( j+1 ), betas(j) ) ;
					callback( prefix + m_effect_parameter_names.parameter_name(j), betas(j) ) ;
				}
				for( int j = 0; j < betas.size(); ++j ) {
					callback( prefix + m_effect_parameter_names.se_name(j), ses(j) ) ;
				}
				{
					int index = 0 ;
					for( int j = 0; j < betas.size(); ++j ) {
						for( int k = j+1; k < betas.size(); ++k, ++index ) {
							callback( prefix + m_effect_parameter_names.covariance_name(j,k), covariance( index ) ) ;
						}
					}
					assert( index == covariance.size() ) ;
				}
				callback( prefix + "pvalue", pvalue ) ;
			}
		}
	}
	
	std::string PerCohortValueReporter::get_spec() const {
		return "PerCohortValueReporter" ;
	}
	
	std::string PerCohortValueReporter::get_summary( std::string const& prefix, std::size_t column_width ) const {
		return prefix + get_spec() ;
	}
}
