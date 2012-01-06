#include <boost/foreach.hpp>
#define foreach BOOST_FOREACH
#include <boost/function.hpp>
#include <boost/tuple/tuple.hpp>
#include <boost/thread.hpp>
#include "genfile/SNPIdentifyingData.hpp"
#include "genfile/VariantEntry.hpp"
#include "genfile/vcf/get_set_eigen.hpp"
#include "genfile/SNPDataSourceProcessor.hpp"
#include "statfile/BuiltInTypeStatSink.hpp"
#include "db/Connection.hpp"
#include "db/SQLStatement.hpp"
#include "SNPSummaryComponent.hpp"
#include "../qctool_version_autogenerated.hpp"

void SNPSummaryComputationManager::add_computation( std::string const& name, SNPSummaryComputation::UniquePtr computation ) {
	m_computations.insert( name, computation ) ;
}

void SNPSummaryComputationManager::add_result_callback( ResultCallback callback ) {
	m_result_signal.connect( callback ) ;
}

void SNPSummaryComputationManager::begin_processing_snps( std::size_t number_of_samples, std::size_t number_of_snps ) {
	m_snp_index = 0 ;
}

void SNPSummaryComputationManager::processed_snp( genfile::SNPIdentifyingData const& snp, genfile::VariantDataReader& data_reader ) {
	{
		genfile::vcf::GenotypeSetter< Eigen::MatrixBase< SNPSummaryComputation::Genotypes > > setter( m_genotypes ) ;
		data_reader.get( "genotypes", setter ) ;
	}
	Computations::const_iterator i = m_computations.begin(), end_i = m_computations.end() ;
	for( ; i != end_i; ++i ) {
		i->second->operator()(
			snp,
			m_genotypes,
			boost::bind(
				boost::ref( m_result_signal ),
				m_snp_index,
				snp,
				i->first,
				_1,
				_2
			)
		) ;
	}
	++m_snp_index ;
}

void SNPSummaryComputationManager::end_processing_snps() {}

void SNPSummaryComponent::declare_options( appcontext::OptionProcessor& options ) {
	options.declare_group( "SNP computation options" ) ;
	options[ "-snp-stats" ]
		.set_description( "Calculate and output per-SNP statistics.  This implies that no SNP filtering options are used." ) ;
    options[ "-snp-stats-file" ]
        .set_description( 	"Override the auto-generated path(s) of the snp-stats file to use when outputting snp-wise statistics.  "
							"(By default, the paths are formed by adding \".snp-stats\" to the input gen filename(s).)  "
							"The '#' character can also be used here to specify one output file per chromosome." )
        .set_takes_values_per_use(1)
		.set_maximum_multiplicity(1)
		.set_default_value( "" ) ;

	options[ "-snp-stats-columns" ]
        .set_description( "Comma-seperated list of extra columns to output in the snp-wise statistics file.  "
 						"The standard columns are: "
						"SNPID, RSID, position, minor_allele, major_allele, MAF, HWE, missing, information."
						" Your choices here are old_information, jonathans_information, mach_r2, and entropy." )
		.set_takes_single_value()
		.set_default_value( "alleles,HWE,missingness,information" ) ;

	options.option_implies_option( "-snp-stats", "-cohort-name" ) ;
}

SNPSummaryComponent::SNPSummaryComponent( appcontext::OptionProcessor const& options ):
	m_options( options )
{}

namespace {
	struct FileOutputter: public boost::noncopyable {
		typedef std::auto_ptr< FileOutputter > UniquePtr ;
		typedef boost::shared_ptr< FileOutputter > SharedPtr ;
		
		static UniquePtr create( std::string const& filename ) { return UniquePtr( new FileOutputter( filename ) ) ; }
		static SharedPtr create_shared( std::string const& filename ) { return SharedPtr( new FileOutputter( filename ) ) ; }

		FileOutputter( std::string const& filename ):
			m_filename( filename ),
			m_sink( statfile::BuiltInTypeStatSink::open( filename ))
		{
			(*m_sink) | "SNPID" | "rsid" | "chromosome" | "position" | "alleleA" | "alleleB" | "computation_name" | "variable" | "value" ;
		}

		void operator()(
			std::size_t index,
			genfile::SNPIdentifyingData const& snp,
			std::string const& computation_name,
			std::string const& value_name,
			genfile::VariantEntry const& value
		) {
			(*m_sink)
				<< snp.get_SNPID()
				<< snp.get_rsid()
				<< std::string( snp.get_position().chromosome() )
				<< snp.get_position().position()
				<< snp.get_first_allele()
				<< snp.get_second_allele()
				<< computation_name
				<< value_name ;
			if( value == value ) {
				(*m_sink) << value ;
			}
			else {
				(*m_sink) << "NA" ;
			}
			(*m_sink) << statfile::end_row() ;
			;
		}

	private:
		std::string const m_filename ;
		statfile::BuiltInTypeStatSink::UniquePtr m_sink ;
	} ;

	struct DBOutputter {
		typedef std::auto_ptr< DBOutputter > UniquePtr ;
		typedef boost::shared_ptr< DBOutputter > SharedPtr ;

		static UniquePtr create( std::string const& filename, std::string const& cohort_name, std::string const& exclusions_name ) {
			return UniquePtr( new DBOutputter( filename, cohort_name, exclusions_name ) ) ;
		}
		static SharedPtr create_shared( std::string const& filename, std::string const& cohort_name, std::string const& exclusions_name ) {
			return SharedPtr( new DBOutputter( filename, cohort_name, exclusions_name ) ) ;
		}

		DBOutputter( std::string const& filename, std::string const& cohort_name, std::string const& exclusions_name ):
			m_connection( db::Connection::create( filename )),
			m_max_transaction_count( 100 ),
			m_cohort_name( cohort_name ),
			m_exclusions_name( exclusions_name )
		{
			db::Connection::ScopedTransactionPtr transaction = m_connection->open_transaction() ;
			m_connection->run_statement(
				"CREATE TABLE IF NOT EXISTS Variant ( id INTEGER PRIMARY KEY, snpid TEXT, rsid TEXT, chromosome TEXT, position INTEGER, alleleA TEXT, alleleB TEXT )"
			) ;
			m_connection->run_statement(
				"CREATE INDEX IF NOT EXISTS Variant_index ON Variant( chromosome, position, rsid )"
			) ;
			m_connection->run_statement(
				"CREATE TABLE IF NOT EXISTS Entity ( id INTEGER PRIMARY KEY, name TEXT, description TEXT )"
			) ;
			m_connection->run_statement(
				"CREATE TABLE IF NOT EXISTS EntityData ( "
				"entity_id INTEGER NOT NULL, "
				"variable_id INTEGER NOT NULL, "
				"value TEXT, "
				"FOREIGN KEY (entity_id) REFERENCES Entity( id ), "
				"FOREIGN KEY (variable_id) REFERENCES Entity( id ) "
				")"
			) ;
			m_connection->run_statement(
				"CREATE TABLE IF NOT EXISTS SummaryData ( "
				"variant_id INT, variable_id INT, value FLOAT, "
				"FOREIGN KEY( variant_id ) REFERENCES Variant( id ), "
				"FOREIGN KEY( variable_id ) REFERENCES Entity( id ), "
				"UNIQUE( variant_id, variable_id ) "
				")"
			) ;
			m_connection->run_statement(
				"CREATE INDEX IF NOT EXISTS SummaryDataIndex ON SummaryData( variant_id, variable_id )"
			) ;
			m_connection->run_statement(
				"CREATE INDEX IF NOT EXISTS EntityDataIndex ON EntityData( entity_id, variable_id )"
			) ;
			m_connection->run_statement(
				"CREATE VIEW IF NOT EXISTS SNPStatsView AS "
				"SELECT V.chromosome, V.position, V.rsid, cohort.value AS cohort, tool.value AS tool, variable.name AS variable, SD.value "
				"FROM SummaryData SD "
				"INNER JOIN Variant V ON V.id == SD.variant_id "
				"INNER JOIN Entity variable ON variable.id = SD.variable_id "
				"LEFT OUTER JOIN EntityData tool ON tool.entity_id = variable.id AND tool.variable_id = ( SELECT id FROM Entity WHERE name = 'tool' )"
				"LEFT OUTER JOIN EntityData cohort ON cohort.entity_id = variable.id AND cohort.variable_id = ( SELECT id FROM entity WHERE name = 'cohort') "
			) ;
			m_connection->run_statement(
				"CREATE VIEW SNPFilterView AS "
				"SELECT          V.chromosome, V.position, V.rsid, "
				"MAF.value AS 'MAF', "
				"HWE.value AS 'minus_log10_exact_HW_p_value', "
				"Missing.value AS 'missing_proportion', "
				"Info.value AS 'impute_info' "
				"FROM            Variant V "
				"LEFT OUTER JOIN      SummaryData Missing "
				"ON          Missing.variable_id = ( SELECT id FROM Entity WHERE name = 'missing proportion' ) "
				"AND         Missing.variant_id == V.id "
				"LEFT OUTER JOIN      SummaryData MAF "
				"ON          MAF.variable_id = ( SELECT id FROM Entity WHERE name = 'minor_allele_frequency' ) "
				"AND         MAF.variant_id == V.id "
				"LEFT OUTER JOIN      SummaryData HWE "
				"ON          HWE.variant_id == V.id "
				"AND         HWE.variable_id == ( SELECT id FROM Entity WHERE name == 'minus_log10_exact_HW_p_value' ) "
				"LEFT OUTER JOIN      SummaryData Info "
				"ON          Info.variant_id == V.id "
				"AND         Info.variable_id == ( SELECT id FROM Entity WHERE name == 'impute_info' ) "
			) ;
			construct_statements() ;
		}

		~DBOutputter() {
			write_data( m_data ) ;
		}

		void operator()(
			std::size_t index,
			genfile::SNPIdentifyingData const& snp,
			std::string const& computation_name,
			std::string const& variable,
			genfile::VariantEntry const& value
		) {
			m_data.resize( m_data.size() + 1 ) ;
			m_data.back().get<0>() = snp ;
			m_data.back().get<1>() = variable ;
			m_data.back().get<2>() = value ;

			if( m_data.size() == m_max_transaction_count ) {
				write_data( m_data ) ;
				m_data.clear() ;
			}
		}

	private:
		db::Connection::UniquePtr m_connection ;
		std::size_t const m_max_transaction_count ;
		std::string const m_cohort_name ;
		std::string const m_exclusions_name ;

		db::Connection::StatementPtr m_find_variant_statement ;
		db::Connection::StatementPtr m_insert_variant_statement ;
		db::Connection::StatementPtr m_find_entity_statement ;
		db::Connection::StatementPtr m_insert_entity_statement ;
		db::Connection::StatementPtr m_insert_summarydata_statement ;
		db::Connection::StatementPtr m_insert_entity_data_statement ;
		typedef std::vector< boost::tuple< genfile::SNPIdentifyingData, std::string, genfile::VariantEntry > > Data ;
		Data m_data ;

	private:
		void construct_statements() {
			m_find_variant_statement = m_connection->get_statement(
				"SELECT id FROM Variant WHERE rsid == ?1 AND chromosome == ?2 AND position == ?3"
			) ;
			m_insert_variant_statement = m_connection->get_statement(
				"INSERT INTO Variant ( snpid, rsid, chromosome, position, alleleA, alleleB) "
				"VALUES( ?1, ?2, ?3, ?4, ?5, ?6 )"
			) ;
			m_find_entity_statement = m_connection->get_statement( "SELECT * FROM Entity WHERE name == ?1" ) ;
			m_insert_entity_statement = m_connection->get_statement( "INSERT INTO Entity ( name ) VALUES ( ?1 )" ) ;
			m_insert_entity_data_statement = m_connection->get_statement( "INSERT INTO EntityData ( entity_id, variable_id, value ) VALUES ( ?1, ?2, ?3 )" ) ;
			m_insert_summarydata_statement = m_connection->get_statement(
				"INSERT OR REPLACE INTO SummaryData ( variant_id, variable_id, value ) "
				"VALUES( ?1, ?2, ?3 )"
			) ;
		}

		void write_data( Data const& data ) {
			db::Connection::ScopedTransactionPtr transaction ;

			for( std::size_t i = 0; i < 100; ++i ) {
				try {
					transaction = m_connection->open_transaction() ;
					break ;
				}
				catch( db::StatementStepError const& e ) {
					// wait a tenth of a second
					std::cerr << "SNPSummaryComponent::DBOutputter::write_data(): failed to open transaction, trying again in 0.1s...\n" ;
					boost::this_thread::sleep( boost::posix_time::milliseconds( 100 ) ) ;
				}
				catch( ... ) {
					std::cerr << "SNPSummaryComponent::write_data(): OMG, a strange exception was caught.\n" ;
					boost::this_thread::sleep( boost::posix_time::milliseconds( 100 ) ) ;
				}
			}
			if( !transaction.get() ) {
				throw genfile::OperationFailedError( "SNPSummaryComponent::write_data()", m_connection->get_spec(), "Opening transaction." ) ;
			}
			for( std::size_t i = 0; i < data.size(); ++i ) {
				store_data(
					data[i].get<0>(),
					data[i].get<1>(),
					data[i].get<2>()
				) ;
			}
		}

		db::Connection::RowId get_or_create_snp( genfile::SNPIdentifyingData const& snp ) const {
			m_find_variant_statement->reset()
				.bind( 1, snp.get_rsid() )
				.bind( 2, std::string( snp.get_position().chromosome() ))
				.bind( 3, snp.get_position().position() )
				.step()
			;
			if( m_find_variant_statement->empty() ) {
				m_insert_variant_statement
					->reset()
					.bind( 1, snp.get_SNPID() )
					.bind( 2, snp.get_rsid() )
					.bind( 3, std::string( snp.get_position().chromosome() ) )
					.bind( 4, snp.get_position().position() )
					.bind( 5, snp.get_first_allele())
					.bind( 6, snp.get_second_allele())
					.step()
				;

				return m_connection->get_last_insert_row_id() ;
			} else {
				return m_find_variant_statement->get< db::Connection::RowId >( 0 ) ;
			}
		}

		db::Connection::RowId get_or_create_variable( std::string const& name ) const {
			db::Connection::RowId result ;

			m_find_entity_statement
				->reset()
				.bind( 1, name ).step() ;

			if( m_find_entity_statement->empty() ) {
				m_insert_entity_statement
					->reset()
					.bind( 1, name )
					.step() ;
					
				result = m_connection->get_last_insert_row_id() ;
				
				m_insert_entity_data_statement
					->reset()
					.bind( 1, result )
					.bind( 2, get_or_create_entity( "tool" ))
					.bind( 3, "qctool revision " + std::string( globals::qctool_revision ) )
					.step()
				;

				m_insert_entity_data_statement
					->reset()
					.bind( 1, result )
					.bind( 2, get_or_create_entity( "cohort" ))
					.bind( 3, m_cohort_name )
					.step()
				;

				m_insert_entity_data_statement
					->reset()
					.bind( 1, result )
					.bind( 2, get_or_create_entity( "sample_exclusions" ))
					.bind( 3, m_exclusions_name )
					.step()
				;
			} else {
				result = m_find_entity_statement->get< db::Connection::RowId >( 0 ) ;
			}
			return result ;
		}

		db::Connection::RowId get_or_create_entity( std::string const& name ) const {
			db::Connection::RowId result ;

			m_find_entity_statement
				->reset()
				.bind( 1, name ).step() ;

			if( m_find_entity_statement->empty() ) {
				m_insert_entity_statement
					->reset()
					.bind( 1, name )
					.step() ;
					
				result = m_connection->get_last_insert_row_id() ;
			} else {
				result = m_find_entity_statement->get< db::Connection::RowId >( 0 ) ;
			}
			return result ;
		}

		void store_data(
			genfile::SNPIdentifyingData const& snp,
			std::string const& variable,
			genfile::VariantEntry const& value
		) {

			db::Connection::RowId snp_id = get_or_create_snp( snp ) ;
			db::Connection::RowId variable_id = get_or_create_variable( variable );

			assert( m_insert_summarydata_statement.get() ) ;
			m_insert_summarydata_statement
				->reset()
				.bind( 1, snp_id )
				.bind( 2, variable_id )
				.bind( 3, value )
				.step()
			;
		}
	} ;
}

genfile::SNPDataSourceProcessor::Callback::UniquePtr SNPSummaryComponent::create() const {
	genfile::SNPDataSourceProcessor::Callback::UniquePtr result( create_manager().release() ) ;
	return result ;
}
	

SNPSummaryComputationManager::UniquePtr SNPSummaryComponent::create_manager() const {
	SNPSummaryComputationManager::UniquePtr manager( new SNPSummaryComputationManager() ) ;
	std::vector< std::string > elts = genfile::string_utils::split_and_strip_discarding_empty_entries( m_options.get_value< std::string >( "-snp-stats-columns" ), ",", " \t" ) ;
	foreach( std::string const& elt, elts ) {
		manager->add_computation( elt, SNPSummaryComputation::create( elt )) ;
	}

	std::string filename = m_options.get_value< std::string >( "-snp-stats-file" ) ;
	if( filename.empty() ) {
		filename = m_options.get< std::string >( "-g" ) + ".snp-stats" ;
	}
	
	if( m_options.check( "-nodb" )) {
		manager->add_result_callback(
			boost::bind(
				&FileOutputter::operator(),
				FileOutputter::create_shared( filename ),
				_1, _2, _3, _4, _5
			)
		) ;
	}
	else {
			std::string sample_set_spec = "" ;
		if( m_options.check( "-excl-samples" )) {
			sample_set_spec += "excluded:" + genfile::string_utils::join( m_options.get_values< std::string >( "-excl-samples"  ), "," ) ;
		}
		if( m_options.check( "-incl-samples" )) {
			sample_set_spec += "included:" + genfile::string_utils::join( m_options.get_values< std::string >( "-incl-samples"  ), "," ) ;
		}
	
		DBOutputter::SharedPtr outputter = DBOutputter::create_shared(
			filename,
			m_options.get< std::string >( "-cohort-name" ), 
			sample_set_spec
		) ;
	
		manager->add_result_callback(
			boost::bind(
				&DBOutputter::operator(),
				outputter,
				_1, _2, _3, _4, _5
			)
		) ;
	}
	return manager ;
}

SNPSummaryComputation::UniquePtr SNPSummaryComponent::create_computation( std::string const& name ) const {
	return SNPSummaryComputation::UniquePtr( SNPSummaryComputation::create( name )) ;
}
