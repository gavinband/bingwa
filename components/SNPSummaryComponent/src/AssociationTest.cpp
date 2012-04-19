
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#include <vector>
#include <string>
#include "genfile/CohortIndividualSource.hpp"
#include "appcontext/OptionProcessor.hpp"
#include "components/SNPSummaryComponent/SNPSummaryComputation.hpp"
#include "components/SNPSummaryComponent/AssociationTest.hpp"
#include "components/SNPSummaryComponent/FrequentistCaseControlAssociationTest.hpp"
#include "components/SNPSummaryComponent/AutosomalFrequentistCaseControlAssociationTest.hpp"
#include "components/SNPSummaryComponent/XChromosomeFrequentistCaseControlAssociationTest.hpp"

namespace {
	AssociationTest::Vector get_sample_column(
		genfile::CohortIndividualSource const& samples,
		std::string const& column_name
	) {
		AssociationTest::Vector result = AssociationTest::Vector::Constant( samples.get_number_of_individuals(), std::numeric_limits< double >::quiet_NaN() ) ;
		for( std::size_t i = 0; i < samples.get_number_of_individuals(); ++i ) {
			genfile::CohortIndividualSource::Entry const& entry = samples.get_entry( i, column_name ) ;
			if( !entry.is_missing() ) {
				result(i) = entry.as< double >() ;
			}
		}
		return result ;
	}
}

SNPSummaryComputation::UniquePtr AssociationTest::create(
	std::string const& type,
	std::string const& phenotype_name,
	std::vector< std::string > const& covariate_names,
	genfile::CohortIndividualSource const& samples,
	appcontext::OptionProcessor const& options
) {
	Vector phenotypes = get_sample_column( samples, phenotype_name ) ;
	Matrix covariates( phenotypes.rows(), covariate_names.size() ) ;
	for( std::size_t i = 0; i < covariate_names.size(); ++i ) {
		covariates.col(i) = get_sample_column( samples, covariate_names[i] ) ;
	}
	
	SNPSummaryComputation::UniquePtr result ;
	if( type == "autosomal" ) {
		result.reset( new AutosomalFrequentistCaseControlAssociationTest( samples, phenotypes, covariates ) ) ;
	}
	else if( type == "X chromosome" ) {
		result.reset(
			new XChromosomeFrequentistCaseControlAssociationTest( samples, phenotypes, covariates, !options.check( "-no-X-inactivation" ) )
		) ;
	} else {
		throw genfile::BadArgumentError( "AssociationTest::create()", "type=\"" + type + "\"" ) ;
	}

	return result ;
}
