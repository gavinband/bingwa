#ifndef QCTOOL_SNP_SUMMARY_COMPUTATION_HPP
#define QCTOOL_SNP_SUMMARY_COMPUTATION_HPP

#include <string>
#include <memory>
#include <boost/noncopyable.hpp>
#include <boost/function.hpp>
#include <eigen/Core>
#include "genfile/SNPIdentifyingData.hpp"
#include "genfile/VariantEntry.hpp"
#include "appcontext/OptionProcessor.hpp"

struct SNPSummaryComputation: public boost::noncopyable {
	typedef std::auto_ptr< SNPSummaryComputation > UniquePtr ;
	virtual ~SNPSummaryComputation() {}
	static UniquePtr create( std::string const& name ) ;
	static void list_computations( boost::function< void ( std::string ) > callback ) ;

	typedef genfile::SNPIdentifyingData SNPIdentifyingData ;
	typedef Eigen::MatrixXd Genotypes ;
	typedef boost::function< void ( std::string const& value_name, genfile::VariantEntry const& value ) > ResultCallback ;
	virtual void operator()( SNPIdentifyingData const&, Genotypes const&, ResultCallback ) const = 0 ;
} ;

#endif