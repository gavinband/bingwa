
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
	PerCohortValueReporter::PerCohortValueReporter(
		std::vector< std::string > const& cohort_names,
		std::vector< std::vector< std::string > > const& cohort_variables,
		bool report_counts
	):
		m_cohort_names( cohort_names ),
		m_variables( cohort_variables ),
		m_report_counts( report_counts )
	{
		std::vector< std::string > defaults ;
		defaults.push_back( "A" ) ;
		defaults.push_back( "B" ) ;
		defaults.push_back( "AA" ) ;
		defaults.push_back( "AB" ) ;
		defaults.push_back( "BB" ) ;
		defaults.push_back( "NULL" ) ;
		defaults.push_back( "N" ) ;
		defaults.push_back( "B_allele_frequency" ) ;
		defaults.push_back( "trusted" ) ;
		defaults.push_back( "pvalue" ) ;
		for( std::size_t cohort = 0; cohort < m_variables.size(); ++cohort ) {
			for( std::size_t j = 0; j < defaults.size(); ++j ) {
				std::vector< std::string >::iterator where = std::find( m_variables[cohort].begin(), m_variables[cohort].end(), defaults[j] ) ;
				if( where != m_variables[cohort].end() ) {
					m_variables[cohort].erase( where ) ;
				}
			}
		}
	}
	
	void PerCohortValueReporter::set_effect_parameter_names( EffectParameterNamePack const& names ) {
		m_effect_parameter_names = names ;
	}
	
	void PerCohortValueReporter::get_variables( boost::function< void ( std::string, std::string ) > callback ) const {
		using genfile::string_utils::to_string ;
		
		std::size_t const N = m_cohort_names.size() ;
		for( std::size_t cohort = 0; cohort < N; ++cohort ) {
			std::string prefix = m_cohort_names[ cohort ] + ":" ;

			if( m_report_counts ) {
				callback( prefix + "A", "FLOAT" ) ;
				callback( prefix + "B", "FLOAT" ) ;
				callback( prefix + "AA", "FLOAT" ) ;
				callback( prefix + "AB", "FLOAT" ) ;
				callback( prefix + "BB", "FLOAT" ) ;
				callback( prefix + "NULL", "FLOAT" ) ;
			}

			callback( prefix + "N", "FLOAT" ) ;
			callback( prefix + "B_allele_frequency", "FLOAT" ) ;
			
			BOOST_FOREACH( std::string const& variable, m_variables[cohort] ) {
				callback( prefix + variable, "NULL" ) ;
			}

			callback( prefix + "trusted", "INTEGER" ) ;

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
			for( std::size_t k = 0; k < m_variables[cohort].size(); ++k ) {
				callback( prefix + m_variables[cohort][k], "NULL" ) ;
			}
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
				Eigen::VectorXd counts ;
				double pvalue ;
				double info ;
				data_getter.get_betas( i, &betas ) ;
				data_getter.get_ses( i, &ses ) ;
				data_getter.get_covariance_upper_triangle( i, &covariance ) ;
				data_getter.get_pvalue( i, &pvalue ) ;
				data_getter.get_info( i, &info ) ;
				data_getter.get_counts( i, &counts ) ;

				using genfile::string_utils::to_string ;
				std::string prefix = m_cohort_names[ i ] + ":" ;

				if( m_report_counts ) {
					callback( prefix + "A", counts(0) ) ;
					callback( prefix + "B", counts(1) ) ;
					callback( prefix + "AA", counts(2) ) ;
					callback( prefix + "AB", counts(3) ) ;
					callback( prefix + "BB", counts(4) ) ;
					callback( prefix + "NULL", counts(5) ) ;
				}

				callback( prefix + "N", counts.segment(0,5).sum() ) ;
				
				double B_allele_count = 0 ;
				double total_allele_count = 0 ;
				if( counts(0) == counts(0) ) {
					B_allele_count += counts(1) ;
					total_allele_count += counts(0) + counts(1) ;
				}
				if( counts(2) == counts(2) ) {
					B_allele_count += counts(3) + 2 * counts(4) ;
					total_allele_count += 2.0 * ( counts(2) + counts(3) + counts(4) ) ;
				}
				callback( prefix + "B_allele_frequency", B_allele_count / total_allele_count ) ;
				{
					std::string value ;
					BOOST_FOREACH( std::string const& variable, m_variables[i] ) {
						data_getter.get_variable( variable, i, &value ) ;
						try {
							callback( prefix + variable, genfile::string_utils::to_repr< double >( value ) ) ;
						} catch( ... ) {
							callback( prefix + variable, value ) ;
						}
					}
				}
				
				callback( prefix + "trusted", genfile::VariantEntry::Integer( data_getter.is_trusted( i ) )) ;
				
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
