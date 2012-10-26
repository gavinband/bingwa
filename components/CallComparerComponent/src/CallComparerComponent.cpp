
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#include <string>
#include <boost/signals2/signal.hpp>
#include <boost/function.hpp>
#include <boost/tuple/tuple.hpp>


#include "genfile/SNPDataSourceProcessor.hpp"
#include "genfile/SNPIdentifyingData.hpp"
#include "genfile/VariantEntry.hpp"
#include "genfile/Error.hpp"
#include "genfile/get_set_eigen.hpp"
#include "genfile/vcf/get_set_eigen.hpp"
#include "statfile/BuiltInTypeStatSink.hpp"
#include "db/Connection.hpp"
#include "db/SQLStatement.hpp"
#include "components/CallComparerComponent/PairwiseCallComparer.hpp"
#include "components/CallComparerComponent/PairwiseCallComparerManager.hpp"
#include "components/CallComparerComponent/CallComparerComponent.hpp"
#include "components/CallComparerComponent/FrequentistTestCallMerger.hpp"
#include "components/CallComparerComponent/CallComparerDBOutputter.hpp"
#include "components/CallComparerComponent/CallComparerFileOutputter.hpp"
#include "components/CallComparerComponent/LeastMissingConsensusCaller.hpp"
#include "components/SNPSummaryComponent/DBOutputter.hpp"
#include "components/SNPSummaryComponent/FlatFileOutputter.hpp"


CallComparerProcessor::UniquePtr CallComparerProcessor::create( PairwiseCallComparerManager::UniquePtr comparer, std::vector< std::string > const& call_fields ) {
	UniquePtr result(
		new CallComparerProcessor( comparer, call_fields )
	) ;
	return result ;
}

CallComparerProcessor::CallComparerProcessor( PairwiseCallComparerManager::UniquePtr call_comparer, std::vector< std::string > const& call_fields  ):
	m_call_comparer( call_comparer ),
	m_call_fields( call_fields )
{}

void CallComparerProcessor::begin_processing_snps( std::size_t number_of_samples ) {
	m_call_comparer->begin_processing_snps( number_of_samples ) ;
}

void CallComparerProcessor::processed_snp( genfile::SNPIdentifyingData const& snp, genfile::VariantDataReader& data_reader ) {
	// Add all the calls to the call comparer.
	m_call_comparer->begin_processing_snp( snp ) ;
	for( std::size_t i = 0; i < m_call_fields.size(); ++i ) {
		genfile::vcf::GenotypeSetter< Eigen::MatrixBase< Eigen::MatrixXd > > setter( m_calls ) ;
		data_reader.get( m_call_fields[i], setter ) ;
		m_call_comparer->add_calls( m_call_fields[i], m_calls ) ;
	}
	m_call_comparer->end_processing_snp() ;
}

void CallComparerProcessor::end_processing_snps() {}

void CallComparerComponent::declare_options( appcontext::OptionProcessor& options ) {
	options.declare_group( "Call comparison options" ) ;
	options[ "-compare-calls" ]
		.set_description( "Compare genotype calls from the given fields of a VCF or other file. "
		 	"The value should be a comma-separated list of genotype call fields in the data.")
		.set_takes_values_until_next_option() ;
	options[ "-compare-calls-file" ]
		.set_description( "Name of output file to put call comparisons in." )
		.set_takes_single_value()
		.set_default_value( "call-comparisons.qcdb" ) ;
	options[ "-consensus-call" ]
		.set_description( "Combine calls from several different genotypers into a single consensus callset. "
			"The calls to use are specified using the -compare-calls option. "
			"A consensus call at each SNP is made when no significant difference is found "
			"between N-1 of the N callsets." )
		.set_takes_single_value()
		.set_default_value( "call-comparisons.vcf" ) ;

	options[ "-comparison-model" ]
		.set_description( "Set the model used for call comparison.  Currently this must be \"AlleleFrequencyTest\". or \"AcceptAll\"." )
		.set_takes_single_value()
		.set_default_value( "AlleleFrequencyTest" ) ;

	options[ "-consensus-model" ]
		.set_description( "Model used to combine calls in the consensus set of calls. "
		 	"Currently this can be \"LeastMissing\" or \"QuangStyle\"." )
		.set_takes_single_value()
		.set_default_value( "LeastMissing" ) ;
	options[ "-consensus-call-pvalue-threshhold" ]
		.set_description( "Treat calls in call comparison treated as distinct if the p-value of an association test between them "
			"is less than or equal to this value." )
		.set_takes_single_value()
		.set_default_value( 0.001 ) ;

	options.option_implies_option( "-consensus-call", "-s" ) ;
	options.option_implies_option( "-compare-calls", "-s" ) ;
	options.option_implies_option( "-consensus-call", "-compare-calls" ) ;
	options.option_implies_option( "-consensus-call-pvalue-threshhold", "-compare-calls" ) ;
}

CallComparerComponent::UniquePtr CallComparerComponent::create(
	genfile::CohortIndividualSource const& samples,
	appcontext::OptionProcessor const& options,
	appcontext::UIContext& ui_context
) {
	CallComparerComponent::UniquePtr result(
		new CallComparerComponent( samples, options, ui_context )
	) ;
	return result ;
}

CallComparerComponent::CallComparerComponent(
	genfile::CohortIndividualSource const& samples,
	appcontext::OptionProcessor const& options,
	appcontext::UIContext& ui_context
):
	m_samples( samples ),
	m_options( options ),
	m_ui_context( ui_context )
{}

namespace {
	genfile::VariantEntry get_sample_entry( genfile::CohortIndividualSource const& samples, std::string const& name, std::size_t i ) {
		return samples.get_entry( i, name ) ;
	}
	
	void send_results_to_sink(
		genfile::SNPDataSink::SharedPtr sink,
		genfile::SNPIdentifyingData const& snp,
		Eigen::MatrixXd const& genotypes,
		std::map< std::string, std::vector< genfile::VariantEntry > > const& info
	) {
		sink->write_snp(
			genotypes.rows(),
			snp,
			genfile::GenotypeGetter< Eigen::MatrixXd >( genotypes, 0ul ),
			genfile::GenotypeGetter< Eigen::MatrixXd >( genotypes, 1ul ),
			genfile::GenotypeGetter< Eigen::MatrixXd >( genotypes, 2ul ),
			info
		) ;
	}
}

namespace {
	// Adapter to adapt SNPSummaryComponent outputters to be used here.
	struct SNPSummaryOutputter: public PairwiseCallComparerManager::ComparisonClient {
	public:
		static UniquePtr create( snp_summary_component::Storage::UniquePtr outputter ) {
			return UniquePtr( new SNPSummaryOutputter( outputter ) ) ;
		}

		static SharedPtr create_shared( snp_summary_component::Storage::UniquePtr outputter ) {
			return SharedPtr( new SNPSummaryOutputter( outputter ) ) ;
		}

	public:
		SNPSummaryOutputter( snp_summary_component::Storage::UniquePtr outputter ):
			m_outputter( outputter )
		{}

		void begin_comparisons( genfile::SNPIdentifyingData const& snp ) {
			m_snp = snp ;
		} ;

		void end_comparisons() {} 
		
		void set_result(
			std::string const& first_callset,
			std::string const& second_callset,
			std::string const& comparison,
			std::string const& variable,
			genfile::VariantEntry const& value
		) {
			m_outputter->store_per_variant_data(
				m_snp,
				comparison + ":" + first_callset + "/" + second_callset + ":" + variable,
				value
			) ;
		}

	private:
		snp_summary_component::Storage::UniquePtr m_outputter ;
		genfile::SNPIdentifyingData m_snp ;
	} ;
	
}

void CallComparerComponent::setup( genfile::SNPDataSourceProcessor& processor ) const {
	PairwiseCallComparerManager::UniquePtr manager( PairwiseCallComparerManager::create().release() ) ;

	std::string filename ;
	if( m_options.check( "-o" )) {
		filename = m_options.get< std::string >( "-o" ) ;
	} else {
		filename = genfile::strip_gen_file_extension_if_present( m_options.get< std::string >( "-g" ) ) + ( m_options.check( "-flat-file" ) ? ".tsv" : ".qcdb" ) ;
	}

	SNPSummaryOutputter::SharedPtr snp_outputter ;
	
	if( m_options.check( "-flat-file" ) ) {
		snp_outputter = SNPSummaryOutputter::create_shared(
			snp_summary_component::FlatFileOutputter::create(
				filename,
				m_options.get< std::string >( "-analysis-name" ),
				m_options.get_values_as_map()
			)
		) ;
		/*
		manager->send_merge_to(
			CallComparerFileOutputter::create_shared(
				filename,
				m_options.get< std::string >( "-analysis-name" )
			)
		) ;
		*/
	} else {
		snp_outputter = SNPSummaryOutputter::create_shared(
			snp_summary_component::DBOutputter::create(
				filename,
				m_options.get< std::string >( "-analysis-name" ),
				m_options.get_values_as_map()
			)
		) ;
	/*
		manager->send_merge_to(
			CallComparerDBOutputter::create_shared(
				filename,
				m_options.get< std::string >( "-analysis-name" ),
				m_options.get_values_as_map()
			)
		) ;
		*/
	}

	manager->send_comparisons_to( snp_outputter ) ;
	
	std::string const comparison_model = m_options.get< std::string >( "-comparison-model" ) ;
	if( comparison_model != "AcceptAll" ) {
		manager->add_comparer( comparison_model, PairwiseCallComparer::create( comparison_model ) ) ;
	}
	manager->add_comparer( "GenotypeTabulatingCallComparer", PairwiseCallComparer::create( "GenotypeTabulatingCallComparer" ) ) ;

	manager->set_merger( PairwiseCallComparerManager::Merger::create( comparison_model, m_options ) ) ;

	if( m_options.check( "-consensus-call" )) {
		ConsensusCaller::SharedPtr consensus_caller ;
		{
			genfile::SNPDataSink::SharedPtr sink(
				genfile::SNPDataSink::create(
					m_options.get< std::string >( "-consensus-call" )
				).release()
			) ;
			sink->set_sample_names(
				m_samples.get_number_of_individuals(),
				boost::bind( get_sample_entry, boost::ref( m_samples ), "ID_1", _1 )
			) ;
			consensus_caller = ConsensusCaller::create_shared( m_options.get< std::string >( "-consensus-model" ) ) ;
			consensus_caller->send_results_to(
				boost::bind(
					&send_results_to_sink,
					sink,
					_1,
					_2,
					_3
				)
			) ;
		}
		manager->send_merge_to( consensus_caller ) ;
	}
	

	CallComparerProcessor::UniquePtr ccc = CallComparerProcessor::create(
		manager,
		m_options.get_values< std::string >( "-compare-calls" )
	) ;
	
	processor.add_callback( genfile::SNPDataSourceProcessor::Callback::UniquePtr( ccc.release() ) ) ;
}

