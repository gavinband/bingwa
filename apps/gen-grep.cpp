#include <limits>

#include "genfile/SNPDataSource.hpp"
#include "genfile/SNPDataSourceChain.hpp"
#include "genfile/get_set.hpp"
#include "genfile/wildcard.hpp"
#include "appcontext/CmdLineOptionProcessor.hpp"
#include "appcontext/ProgramFlow.hpp"
#include "PositionRange.hpp"
#include "wildcard.hpp"
#include "FileSet.hpp"

using namespace appcontext ;

namespace globals {
	std::string program_name = "gen-grep" ;
}

struct IDDataPrinterOptionProcessor: public CmdLineOptionProcessor
{
	std::string get_program_name() const { return globals::program_name ; }
	void declare_options( OptionProcessor& options ) {
		// Meta-options
		options.set_help_option( "-help" ) ;
		options.declare_group( "File-related options" ) ;
		options[ "-g" ]
			.set_description( "Specify the name (or names) of gen files to process." )
			.set_is_required()
			.set_takes_single_value()
			.set_takes_value_by_position( 1 ) ;
			
		options.declare_group( "SNP Selection options" ) ;
		options[ "-snp-interval" ]
			.set_description( "Only output SNPs in the given interval" )
			.set_takes_values_per_use( 2 )
			.set_maximum_multiplicity( 1 ) ;
			
		options[ "-snp-id-file" ]
			.set_description( "Specify a file containing SNP ids to match" )
			.set_takes_single_value() ;

		options[ "-rsid-file" ]
			.set_description( "Specify a file containing RSids to match" )
			.set_takes_single_value() ;
	}
	
	
	std::string gen_filename() const {
		return get_value< std::string > ( "-g" ) ;
	}

	std::set< std::string > snp_ids() const {
		std::set< std::string > result ;
		if( check_if_option_was_supplied( "-snp-id-file" )) {
			result = read_set_from_file< std::set< std::string > >( get_value< std::string >( "-snp-id-file" ) ) ;
		}
		return result ;
	}

	std::set< std::string > rsids() const {
		std::set< std::string > result ;
		if( check_if_option_was_supplied( "-rsid-file" )) {
			result = read_set_from_file< std::set< std::string > >( get_value< std::string >( "-rsid-file" ) ) ;
		}
		return result ;
	}
	
	bool have_snp_ids() const {
		return check_if_option_was_supplied( "-snp-id-file" ) ;
	}

	bool have_rsids() const {
		return check_if_option_was_supplied( "-rsid-file" ) ;
	}
	
	PositionRange position_range() const {
		if( check_if_option_was_supplied( "-snp-interval" )) {
			std::vector< std::size_t > values = get_values< std::size_t > ( "-snp-interval" ) ;
			assert( values.size() == 2 ) ;
			return PositionRange( values[0], values[1] ) ;
		}
		else {
			return PositionRange( 0, std::numeric_limits< uint32_t >::max() ) ;
		}
	}
} ;


struct IDDataPrinterContext
{
	IDDataPrinterContext( int argc, char** argv ) {
		//write_start_banner() ;
		m_options.process( argc, argv ) ;
		construct_snp_data_source() ;
		if( m_options.have_snp_ids() ) {
			m_snp_ids = m_options.snp_ids() ;
		}
		if( m_options.have_rsids() ) {
			m_rsids = m_options.rsids() ;
		}
		write_preamble() ;
	}

	~IDDataPrinterContext() {
		//write_end_banner() ;
	}

	genfile::SNPDataSource& snp_data_source() { return *m_snp_data_source ; }

	PositionRange const position_range() const { return m_options.position_range() ; }

 	bool have_snp_ids() const { return m_options.have_snp_ids() ; }
	std::set< std::string > const& snp_ids() const { return m_snp_ids ; }

 	bool have_rsids() const { return m_options.have_rsids() ; }
	std::set< std::string > const& rsids() const { return m_rsids ; }
	
private:

	void write_start_banner() {
		std::cout << "Welcome to " << globals::program_name << ".\n"
		 	<< "(C) 2009 University of Oxford\n\n";
	}

	void write_end_banner() {
		std::cout << "\nThank you for using " << globals::program_name << ".\n" ;
	}

	void write_preamble() {
		
	}

	void construct_snp_data_source() {
		genfile::SNPDataSourceChain::UniquePtr chain = genfile::SNPDataSourceChain::create(
			genfile::wildcard::find_files_by_chromosome(
				m_options.gen_filename(),
				genfile::wildcard::eALL_CHROMOSOMES
			)
		) ;
		
		m_snp_data_source.reset( chain.release() ) ;
	}
	
private:
	IDDataPrinterOptionProcessor m_options ;
	std::auto_ptr< genfile::SNPDataSource > m_snp_data_source ;
	std::set< std::string > m_snp_ids ;
	std::set< std::string > m_rsids ;
} ;

struct IDDataPrinter
{
	IDDataPrinter( IDDataPrinterContext& context )
		: m_context( context )
	{}

	void process() {
		std::cout << "SNPID RSID chromosome position allele1 allele2\n" ;
		
		std::string SNPID, RSID ;
		genfile::Chromosome chromosome ;
		uint32_t number_of_samples, SNP_position ;
		char allele1, allele2 ;
		while( m_context.snp_data_source().get_snp_identifying_data(
		 	genfile::set_value( number_of_samples ),
			genfile::set_value( SNPID ),
			genfile::set_value( RSID ),
			genfile::set_value( chromosome ),
			genfile::set_value( SNP_position ),
			genfile::set_value( allele1 ),
			genfile::set_value( allele2 )
		)) {
			if(
				m_context.position_range().contains( SNP_position )
				&& (
					!m_context.have_snp_ids()
					|| ( m_context.snp_ids().find( SNPID ) != m_context.snp_ids().end() )
					)
				&& (
					!m_context.have_rsids()
					|| ( m_context.rsids().find( RSID ) != m_context.rsids().end() )
				)
			) {
				std::cout << SNPID << " " << RSID << " " << chromosome << " " << SNP_position << " " << allele1 << " " << allele2 << "\n" ;
			}
			m_context.snp_data_source().ignore_snp_probability_data() ;
		}
	}
	
private:
	
	IDDataPrinterContext& m_context ;
} ;


int main( int argc, char** argv ) {
	try {
		IDDataPrinterContext context( argc, argv ) ;
		IDDataPrinter printer( context ) ;
		printer.process() ;
	}
	catch( HaltProgramWithReturnCode const& e ) {
		return e.return_code() ;
	}
	return 0 ;
}
