
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#ifndef BINGWA_BINGWA_COMPUTATION_HPP
#define BINGWA_BINGWA_COMPUTATION_HPP

#include <memory>
#include <string>
#include <vector>
#include <map>
#include <boost/function.hpp>
#include <boost/noncopyable.hpp>
#include <Eigen/Core>
#include "genfile/VariantIdentifyingData.hpp"
#include "genfile/VariantEntry.hpp"
#include "EffectParameterNamePack.hpp"
#include "appcontext/OptionProcessor.hpp"

namespace bingwa {
	struct BingwaComputation: public boost::noncopyable {
		typedef std::auto_ptr< BingwaComputation > UniquePtr ;
		static UniquePtr create(
			std::string const& name,
			std::vector< std::string > const& cohort_names,
			appcontext::OptionProcessor const&
		) ;
		virtual ~BingwaComputation() {}
		typedef genfile::VariantIdentifyingData VariantIdentifyingData ;
		typedef boost::function< void ( std::string const& value_name, genfile::VariantEntry const& value ) > ResultCallback ;
	
		struct DataGetter: public boost::noncopyable {
			virtual ~DataGetter() {} ;
			virtual std::size_t get_number_of_cohorts() const = 0 ;
			virtual bool is_non_missing( std::size_t i ) const = 0 ;
			virtual bool is_trusted( std::size_t i ) const = 0 ;
			virtual void get_counts( std::size_t, Eigen::VectorXd* result ) const = 0 ;
			virtual void get_betas( std::size_t i, Eigen::VectorXd* result ) const = 0 ;
			virtual void get_ses( std::size_t i, Eigen::VectorXd* result  ) const = 0 ;
			virtual void get_covariance_upper_triangle( std::size_t i, Eigen::VectorXd* result  ) const = 0 ;
			virtual void get_pvalue( std::size_t i, double* result ) const = 0 ;
			virtual void get_info( std::size_t i, double* result ) const = 0 ;
			virtual void get_variable( std::string const& variable, std::size_t i, std::string* result ) const = 0 ;
		} ;

		typedef boost::function< bool ( DataGetter const&, int i ) > Filter ;
	
		virtual void set_effect_parameter_names( EffectParameterNamePack const& names ) = 0 ;
		virtual void get_variables( boost::function< void ( std::string, std::string ) > ) const = 0 ;
		virtual void operator()(
			VariantIdentifyingData const&,
			DataGetter const& data_getter,
			ResultCallback callback
		) = 0 ;
		virtual std::string get_summary( std::string const& prefix = "", std::size_t column_width = 20 ) const = 0 ;
		virtual std::string get_spec() const = 0 ;
	} ;
}

#endif
