
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
#include "bingwa/PerCohortCountsReporter.hpp"

namespace bingwa {
	PerCohortCountsReporter::PerCohortCountsReporter(
		std::vector< std::string > const& cohort_names,
		std::vector< std::vector< std::string > > const& variables
	):
		m_cohort_names( cohort_names ),
		m_variables( variables )
	{}
	
	void PerCohortCountsReporter::set_effect_parameter_names( EffectParameterNamePack const& names ) {
	}
	
	void PerCohortCountsReporter::get_variables( boost::function< void ( std::string, std::string ) > callback ) const {
		using genfile::string_utils::to_string ;
		
		std::size_t const N = m_cohort_names.size() ;
		for( std::size_t i = 0; i < N; ++i ) {
			std::string prefix = m_cohort_names[ i ] + ":" ;

			callback( prefix + "A", "FLOAT" ) ;
			callback( prefix + "B", "FLOAT" ) ;
			callback( prefix + "AA", "FLOAT" ) ;
			callback( prefix + "AB", "FLOAT" ) ;
			callback( prefix + "BB", "FLOAT" ) ;
			callback( prefix + "NULL", "FLOAT" ) ;
			
			BOOST_FOREACH( std::string const& variable, m_variables[i] ) {
				callback( prefix + variable, "NULL" ) ;
			}
		}
	}
	
	void PerCohortCountsReporter::operator()(
		VariantIdentifyingData const&,
		DataGetter const& data_getter,
		ResultCallback callback
	) {
		std::size_t const N = m_cohort_names.size() ;
		for( std::size_t i = 0; i < N; ++i ) {
			if( data_getter.is_non_missing( i ) ) {
				Eigen::VectorXd counts ;
				double info ;
				double maf ;
				
				data_getter.get_counts( i, &counts ) ;

				assert( counts.size() == 6 ) ;
				using genfile::string_utils::to_string ;
				std::string prefix = m_cohort_names[ i ] + ":" ;

				callback( prefix + "A", counts(0) ) ;
				callback( prefix + "B", counts(1) ) ;
				callback( prefix + "AA", counts(2) ) ;
				callback( prefix + "AB", counts(3) ) ;
				callback( prefix + "BB", counts(4) ) ;
				callback( prefix + "NULL", counts(5) ) ;

				{
					std::string value ;
					BOOST_FOREACH( std::string const& variable, m_variables[i] ) {
						data_getter.get_variable( variable, i, &value ) ;
						callback( prefix + variable, value ) ;
					}
				}
			}
		}
	}
	
	std::string PerCohortCountsReporter::get_spec() const {
		return "PerCohortCountsReporter" ;
	}
	
	std::string PerCohortCountsReporter::get_summary( std::string const& prefix, std::size_t column_width ) const {
		return prefix + get_spec() ;
	}
}
