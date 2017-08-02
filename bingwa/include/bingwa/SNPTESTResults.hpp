
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#ifndef SNPTEST_RESULTS_HPP
#define SNPTEST_RESULTS_HPP

#include <boost/optional.hpp>
#include <boost/function.hpp>
#include "genfile/Chromosome.hpp"
#include "statfile/FilteringStatSource.hpp"
#include "FlatFileFrequentistGenomeWideAssociationResults.hpp"

namespace bingwa {
	struct SNPTESTResults: public FlatFileFrequentistGenomeWideAssociationResults {
		typedef boost::function< bool ( double info, double maf, Row const& betas, Row const& ses, Row const& counts ) > Filter ;

		SNPTESTResults(
			genfile::VariantIdentifyingDataTest::UniquePtr test,
			boost::optional< genfile::Chromosome > const chromosome_hint = boost::optional< genfile::Chromosome >()
		) ;

		void set_effect_size_column_regex( std::string const& beta_column_regex ) ;
		bingwa::EffectParameterNamePack get_effect_parameter_names() const ;
	
		std::string get_summary( std::string const& prefix, std::size_t target_column ) const ;

	private:
		genfile::VariantIdentifyingDataTest::UniquePtr m_exclusion_test ;
		std::string m_beta_column_regex ;
		boost::optional< genfile::Chromosome > m_chromosome_hint ;
		std::vector< std::string > m_beta_columns ;
		std::vector< std::string > m_se_columns ;
		std::vector< std::string > m_cov_columns ;
		std::string m_pvalue_column ;
		std::string m_info_column ;

	private:
		DesiredColumns setup_columns( std::vector< std::string > const& column_names ) ;
		bool read_snp( statfile::BuiltInTypeStatSource& source, genfile::SNPIdentifyingData& snp ) const ;
		bool check_if_snp_accepted( std::size_t snp_i ) const ;
		void store_value(
			int snp_index,
			std::string const& variable,
			std::string const& value
		) ;
	} ;
}

#endif

