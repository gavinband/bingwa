
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)


#ifndef BINGWA_PER_COHORT_COUNTS_REPORTER_HPP
#define BINGWA_PER_COHORT_COUNTS_REPORTER_HPP

#include <memory>
#include <string>
#include <vector>
#include <boost/function.hpp>
#include <Eigen/Core>
#include "genfile/VariantIdentifyingData.hpp"
#include "genfile/VariantEntry.hpp"
#include "appcontext/OptionProcessor.hpp"
#include "bingwa/BingwaComputation.hpp"

namespace bingwa {	
	struct PerCohortCountsReporter: public bingwa::BingwaComputation {
		public:
			typedef std::auto_ptr< PerCohortCountsReporter > UniquePtr ;
	
		public:
			PerCohortCountsReporter(
				std::vector< std::string > const& cohort_names,
				std::vector< std::vector< std::string > > const& variables
			) ;
			void add_variable( std::string const& variable ) ;
			void get_variables( boost::function< void ( std::string, std::string ) > callback ) const ;
			void set_effect_parameter_names( EffectParameterNamePack const& names ) ;

			void operator()(
				VariantIdentifyingData const&,
				DataGetter const& data_getter,
				ResultCallback callback
			) ;

			std::string get_spec() const ;
			std::string get_summary( std::string const& prefix = "", std::size_t column_width = 20 ) const ;

		private:
			std::vector< std::string > const m_cohort_names ;
			std::vector< std::vector< std::string > > m_variables ;
	} ;
}

#endif
