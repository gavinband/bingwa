
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#include <vector>
#include <boost/foreach.hpp>
#include <boost/function.hpp>
#include <boost/tuple/tuple.hpp>
#include <boost/thread.hpp>
#include "genfile/SNPIdentifyingData.hpp"
#include "genfile/SNPDataSourceProcessor.hpp"
#include "genfile/Error.hpp"
#include "../../qctool_version_autogenerated.hpp"
#include "components/SNPSummaryComponent/SNPSummaryComponent.hpp"
#include "components/SNPSummaryComponent/SNPSummaryComputation.hpp"
#include "components/SNPSummaryComponent/SNPSummaryComputationManager.hpp"
#include "qcdb/Storage.hpp"
#include "qcdb/FlatFileOutputter.hpp"
#include "qcdb/FlatTableDBOutputter.hpp"
#include "components/SNPSummaryComponent/DBOutputter.hpp"
#include "components/SNPSummaryComponent/AssociationTest.hpp"
#include "components/SNPSummaryComponent/SequenceAnnotation.hpp"
#include "components/SNPSummaryComponent/DifferentialMissingnessComputation.hpp"
#include "components/SNPSummaryComponent/StratifyingSNPSummaryComputation.hpp"
#include "components/SNPSummaryComponent/CrossDataSetConcordanceComputation.hpp"
#include "components/SNPSummaryComponent/CrossDataSetHaplotypeComparisonComputation.hpp"
#include "components/SNPSummaryComponent/GeneticMapAnnotation.hpp"
#include "components/SNPSummaryComponent/IntensityReporter.hpp"
#include "components/SNPSummaryComponent/CallComparerComponent.hpp"

void SNPSummaryComponent::declare_options( appcontext::OptionProcessor& options ) {
	options.declare_group( "SNP computation options" ) ;
	options[ "-snp-stats" ]
		.set_description( "Calculate and output per-SNP statistics.  This implies that no SNP filtering options are used." ) ;
	options[ "-snp-stats-columns" ]
        .set_description( "Comma-seperated list of extra columns to output in the snp-wise statistics file." )
		.set_takes_single_value()
		.set_default_value( "alleles,HWE,missingness,information" ) ;

	options.declare_group( "Intensity computation options" ) ;
	options[ "-intensity-stats" ]
		.set_description( "Compute intensity means and (co)variances for each genotype class at each SNP." ) ;
	options[ "-intensities" ]
		.set_description( "Report per-sample intensities for each sample at each SNP." ) ;
	
	options.declare_group( "Association test options" ) ;
	options[ "-test" ]
		.set_description( "Perform an association test on the given phenotype." )
		.set_takes_single_value() ;
	options[ "-covariates" ]
		.set_description( "Specify a comma-separated list of covariates to use in the association test." )
		.set_takes_single_value()
		.set_default_value( "" ) ;
	options[ "-no-X-inactivation" ]
		.set_description( "Specify that X chromosome inactivation in females should not be modelled in the association test. "
			"If this option is specified, females have twice the maximum exposure that males do." )
		.set_takes_single_value()
		.set_default_value( "" ) ;
	
	options[ "-stratify" ]
		.set_description( "Compute all SNP summary statistics seperately for each level of the given variable in the sample file." )
		.set_takes_single_value() ;

	options[ "-differential" ]
		.set_description( "Test for differences in SNP summary statistics between the categories of the given variable."
			" Currently a test for differential missingness is performed." )
		.set_takes_single_value() ;
	
	options.declare_group( "Sequence annotation options" ) ;
	options[ "-annotate-ancestral" ]
		.set_description( "Specify a FASTA-formatted file containing ancestral alleles to annotate variants with." )
		.set_takes_single_value() ;
	options[ "-annotate-reference" ]
		.set_description( "Specify a FASTA-formatted file containing ancestral alleles to annotate variants with." )
		.set_takes_single_value() ;
	options[ "-annotate-genetic-map" ]
		.set_description( "Specify a genetic map file or files.  QCTOOL will interpolate the map "
			"to produce approximate positions in centiMorgans for each SNP in the data." )
		.set_takes_single_value() ;
	options[ "-flanking" ]
		.set_description( "Specify that flanking sequence annotations [ pos - a, pos + b ] should be output when using "
			"-annotate-reference and -annotate-ancestral" )
		.set_takes_values( 2 )
		.set_minimum_multiplicity( 0 )
		.set_maximum_multiplicity( 1 ) ;

	options.declare_group( "Callset comparison options" ) ;
	options[ "-compare-to" ]
		.set_description( "Compute a matrix of values indicating concordance of samples between the main dataset and the dataset given as argument to this option. "
		 	"Values must be the genotype and sample files (in that order).  Samples are matched using the first ID column; "
			"SNPs are matched based on all the identifying information fields." )
		.set_takes_values( 2 )
		.set_minimum_multiplicity( 0 )
		.set_maximum_multiplicity( 1 ) ;
	options[ "-haplotypic" ]
		.set_description( "Instruct QCTOOL to perform haplotypic computations.  Currently this affects the -compare-to option only "
			"and turns on computation of switch error for two sets of haplotypes." ) ;
	options.option_implies_option( "-snp-stats", "-g" ) ;
	options.option_implies_option( "-intensity-stats", "-g" ) ;
	options.option_implies_option( "-annotate-ancestral", "-g" ) ;
	options.option_implies_option( "-annotate-reference", "-g" ) ;
	options.option_implies_option( "-test", "-g" ) ;
	options.option_implies_option( "-test", "-s" ) ;
	options.option_implies_option( "-stratify", "-s" ) ;
	
	CallComparerComponent::declare_options( options ) ;
}

SNPSummaryComponent::SNPSummaryComponent(
	genfile::CohortIndividualSource const& samples,
	appcontext::OptionProcessor const& options,
	appcontext::UIContext& ui_context
):
	m_samples( samples ),
	m_options( options ),
	m_ui_context( ui_context )
{}

void SNPSummaryComponent::setup( genfile::SNPDataSourceProcessor& processor ) const {
	processor.add_callback( genfile::SNPDataSourceProcessor::Callback::UniquePtr( create_manager().release() ) ) ;
}

SNPSummaryComputationManager::UniquePtr SNPSummaryComponent::create_manager() const {
	SNPSummaryComputationManager::UniquePtr manager( new SNPSummaryComputationManager( m_samples ) ) ;
	using genfile::string_utils::to_string ;
	
	std::string filename ;
	if( m_options.check( "-o" )) {
		filename = m_options.get_value< std::string >( "-o" ) ;
	}
	else {
		std::vector< std::string > filenames = m_options.get_values< std::string >( "-g" ) ;
		if( filenames.size() == 1 ) {
			filename = genfile::strip_gen_file_extension_if_present( filenames[0] ) + ( m_options.check( "-flat-file" ) ? ".snp-stats.tsv" : ".qcdb" ) ;
		} else {
			filename = "qctool_cohort_1-" + to_string( filenames.size() ) + ( m_options.check( "-flat-file" ) ? ".snp-stats.tsv" : ".qcdb" ) ;
		}
	}

	qcdb::Storage::SharedPtr storage ;

	if( m_options.check( "-flat-file" )) {
		storage = qcdb::FlatFileOutputter::create_shared(
			filename,
			m_options.get< std::string >( "-analysis-name" ),
			m_options.get_values_as_map()
		) ;
	}
	else if( m_options.check( "-flat-table" )) {
		storage = qcdb::FlatTableDBOutputter::create_shared(
			filename,
			m_options.get< std::string >( "-analysis-name" ),
			m_options.get_values_as_map()
		) ;
	}
	else {
		storage = snp_summary_component::DBOutputter::create_shared(
			filename,
			m_options.get< std::string >( "-analysis-name" ),
			m_options.get_values_as_map()
		) ;
	}

	manager->add_result_callback(
		boost::bind(
			&qcdb::Storage::store_per_variant_data,
			storage,
			_1, _2, _3
		)
	) ;
	
	add_computations( *manager, storage ) ;
	m_ui_context.logger() << "SNPSummaryComponent: the following components are in place:\n" << manager->get_summary( "  " ) << "\n" ;
	
	return manager ;
}

namespace impl {
	
	typedef std::map< genfile::VariantEntry, std::vector< int > > StrataMembers ;

	StrataMembers compute_strata( genfile::CohortIndividualSource const& samples, std::string const& variable ) {
		StrataMembers result ;
		genfile::CohortIndividualSource::ColumnSpec const spec = samples.get_column_spec() ;
		if( !spec[ variable ].is_discrete() ) {
			throw genfile::BadArgumentError( "void impl::compute_strata()", "variable=\"" + variable + "\"" ) ;
		}

		for( std::size_t i = 0; i < samples.get_number_of_individuals(); ++i ) {
			genfile::VariantEntry const& entry = samples.get_entry( i, variable ) ;
			if( !entry.is_missing() ) {
				result[ samples.get_entry( i, variable ) ].push_back( int( i ) ) ;
			}
		}

		return result ;
	}
}

void SNPSummaryComponent::add_computations( SNPSummaryComputationManager& manager, qcdb::Storage::SharedPtr storage ) const {
	using genfile::string_utils::split_and_strip_discarding_empty_entries ;

	if( m_options.check( "-snp-stats" )) {
		std::vector< std::string > elts = split_and_strip_discarding_empty_entries( m_options.get_value< std::string >( "-snp-stats-columns" ), ",", " \t" ) ;
		BOOST_FOREACH( std::string const& elt, elts ) {
			manager.add_computation( elt, SNPSummaryComputation::create( elt )) ;
		}
	}

	if( m_options.check( "-intensity-stats" )) {
		manager.add_computation( "intensity-stats", SNPSummaryComputation::create( "intensity-stats" )) ;
	}

	if( m_options.check( "-intensities" )) {
		SNPSummaryComputation::UniquePtr computation(
			new snp_summary_component::IntensityReporter( m_samples )
		) ;
		manager.add_computation( "intensities", computation ) ;
	}

	if( m_options.check( "-test" )) {
		std::vector< std::string > phenotypes = split_and_strip_discarding_empty_entries( m_options.get_value< std::string >( "-test" ), ",", " \t" ) ;
		std::vector< std::string > covariates ;
		if( m_options.check( "-covariates" ))  {
			covariates = split_and_strip_discarding_empty_entries( m_options.get_value< std::string >( "-covariates" ), ",", " \t" ) ;
		}
		BOOST_FOREACH( std::string const& phenotype, phenotypes ) {
			manager.add_computation(
				"association_test",
				AssociationTest::create(
					"autosomal",
					phenotype,
					covariates,
					m_samples,
					m_options
				)
			) ;
			manager.add_computation(
				"X_chromosome_association_test",
				AssociationTest::create(
					"X chromosome",
					phenotype,
					covariates,
					m_samples,
					m_options
				)
			) ;
		}
	}

	if( m_options.check( "-annotate-reference" )) {
		appcontext::UIContext::ProgressContext progress = m_ui_context.get_progress_context( "Loading reference sequence" ) ;
		SequenceAnnotation::UniquePtr computation(
			new SequenceAnnotation( "reference", m_options.get< std::string >( "-annotate-reference" ), progress )
		) ;
		
		if( m_options.check( "-flanking" )) {
			std::vector< std::size_t > data = m_options.get_values< std::size_t >( "-flanking" ) ;
			assert( data.size() == 2 ) ;
			computation->set_flanking( data[0], data[1] ) ;
		}
		
		manager.add_computation(
			"reference_sequence",
			SNPSummaryComputation::UniquePtr(
				computation.release()
			)
		) ;
	}

	if( m_options.check( "-annotate-ancestral" )) {
		appcontext::UIContext::ProgressContext progress = m_ui_context.get_progress_context( "Loading ancestral sequence" ) ;
		SequenceAnnotation::UniquePtr computation(
			new SequenceAnnotation( "ancestral", m_options.get< std::string >( "-annotate-ancestral" ), progress )
		) ;
		
		if( m_options.check( "-flanking" )) {
			std::vector< std::size_t > data = m_options.get_values< std::size_t >( "-flanking" ) ;
			assert( data.size() == 2 ) ;
			computation->set_flanking( data[0], data[1] ) ;
		}

		manager.add_computation(
			"ancestral_sequence",
			SNPSummaryComputation::UniquePtr(
				computation.release()
			)
		) ;
	}

	if( m_options.check( "-annotate-genetic-map" )) {
		appcontext::UIContext::ProgressContext progress = m_ui_context.get_progress_context( "Loading genetic map" ) ;
		GeneticMapAnnotation::UniquePtr computation = GeneticMapAnnotation::create(
			genfile::wildcard::find_files_by_chromosome(
				m_options.get< std::string >( "-annotate-genetic-map" )
			),
			progress
		) ;
		
		manager.add_computation(
			"genetic_map",
			SNPSummaryComputation::UniquePtr(
				computation.release()
			)
		) ;
	}

	if( m_options.check( "-compare-to" )) {
		std::vector< std::string > filenames = m_options.get_values< std::string >( "-compare-to" ) ;

		std::vector< std::string > sample_id_columns = genfile::string_utils::split(
			m_options.get< std::string >( "-match-sample-ids" ),
			"~"
		) ;
		assert( sample_id_columns.size() == 2 ) ;
		
		snp_stats::CrossDataSetComparison::UniquePtr computation ;

		if( m_options.check( "-haplotypic" )) {
			computation.reset(
				new snp_stats::CrossDataSetHaplotypeComparisonComputation(
					m_samples,
					sample_id_columns[0]
				)
			) ;
		}
		else {
			computation.reset(
				new snp_stats::CrossDataSetConcordanceComputation(
					m_samples,
					sample_id_columns[0]
				)
			) ;
		}

		computation->set_comparer(
			genfile::SNPIdentifyingData::CompareFields(
				m_options.get_value< std::string >( "-snp-match-fields" ),
				m_options.check( "-match-alleles-to-cohort1" )
			)
		) ; 
	
		if( m_options.check( "-match-alleles-to-cohort1") ) {
			computation->set_match_alleles() ;
		}
	
		genfile::SNPDataSource::UniquePtr alternate_dataset = genfile::SNPDataSource::create_chain(
			genfile::wildcard::find_files_by_chromosome(
				filenames[0]
			),
			genfile::vcf::MetadataParser::Metadata(),
			m_options.get< std::string >( "-filetype" )
		) ;
		
		computation->set_alternate_dataset(
			genfile::CohortIndividualSource::create( filenames[1] ),
			sample_id_columns[1],
			alternate_dataset
		) ;

		manager.add_computation(
			"dataset_comparison",
			SNPSummaryComputation::UniquePtr(
				computation.release()
			)
		) ;
	}


	if( m_options.check( "-differential" )) {
		std::string const variable = m_options.get< std::string >( "-differential" ) ;
		impl::StrataMembers strata = impl::compute_strata( m_samples, variable ) ;
		manager.add_computation(
			"differential_missingness",
			DifferentialMissingnessComputation::create( variable, strata )
		) ;
	}

	if( m_options.check( "-stratify" )) {
		std::string const variable = m_options.get< std::string >( "-stratify" ) ;
		impl::StrataMembers strata = impl::compute_strata( m_samples, variable ) ;
		manager.stratify_by( strata, variable ) ;
	}
	
	if( m_options.check_if_option_was_supplied_in_group( "Call comparison options" ) ) {
		CallComparerComponent::UniquePtr cc = CallComparerComponent::create(
			m_samples,
			m_options,
			m_ui_context
		) ;

		cc->setup( manager, storage ) ;
	}
}

SNPSummaryComputation::UniquePtr SNPSummaryComponent::create_computation( std::string const& name ) const {
	if( name != "association_test" ) {
		return SNPSummaryComputation::UniquePtr( SNPSummaryComputation::create( name )) ;
	} else {
		assert(0) ;
	}
}
