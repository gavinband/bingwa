#include <string>
#include <map>
#include <memory>
#include <iostream>
#include <boost/tuple/tuple.hpp>
#include "genfile/SNPDataSource.hpp"
#include "genfile/snp_data_utils.hpp"
#include "genfile/string_utils.hpp"
#include "genfile/Error.hpp"
#include "genfile/VCFFormatSNPDataSource.hpp"

using boost::tuples::tie ;

namespace genfile {
	VCFFormatSNPDataSource::VCFFormatSNPDataSource( std::auto_ptr< std::istream > stream_ptr ):
		m_spec( "(unnamed stream)" ),
		m_stream_ptr( stream_ptr ),
		m_metadata( vcf::MetadataParser( m_spec, *m_stream_ptr ).get_metadata() ),
		m_column_names( read_column_names( *m_stream_ptr )),
		m_number_of_lines( count_lines( *m_stream_ptr ))
	{
		reset_stream() ;
	}

	VCFFormatSNPDataSource::VCFFormatSNPDataSource( std::string const& filename ):
		m_spec( "file://" + filename ),
		m_stream_ptr( open_text_file_for_input( filename, get_compression_type_indicated_by_filename( filename ))),
		m_metadata( vcf::MetadataParser( m_spec, *m_stream_ptr ).get_metadata() ),
		m_column_names( read_column_names( *m_stream_ptr )),
		m_number_of_lines( count_lines( *m_stream_ptr ))
	{
	}

	void VCFFormatSNPDataSource::reset_stream() const {
		assert( 0 ) ;
	}

	std::vector< std::string > VCFFormatSNPDataSource::read_column_names( std::istream& stream ) const {
		std::string line ;
		if( !std::getline( stream, line ) ) {
			throw MalformedInputError( m_spec, m_metadata.size() ) ;
		}
		std::vector< std::string > elts = string_utils::split( line, "\t" ) ;
		if( elts.size() < 8 ) {
			throw MalformedInputError( m_spec, m_metadata.size() + 1 ) ;
		}
		else if(
			elts[0] != "#CHROM"
			|| elts[1] != "POS"
			|| elts[2] != "ID"
			|| elts[3] != "REF"
			|| elts[4] != "ALT"
			|| elts[5] != "QUAL"
			|| elts[6] != "FILTER"
			|| elts[7] != "INFO"
			|| elts[8] != "FORMAT"
		) {
			throw MalformedInputError( m_spec, m_metadata.size() ) ;
		}
		return elts ;
	}

	std::size_t VCFFormatSNPDataSource::count_lines( std::istream& str ) const {
		std::size_t count = 0 ;
		std::vector< char > buffer( 10000000 ) ;
		do {
			str.read( &(buffer[0]), 10000000 ) ;
			count += std::count( buffer.begin(), buffer.begin() + str.gcount(), '\n' ) ;
			// A vcf file can't contain a blank line.
			// Because popular editors (vim, nano, ..., but not emacs) typically add a trailing newline,
			// we might get in the situation where the file has two trailing newlines thus messing
			// with our count.
			// Therefore we check here for the special case where what we've read ends in two newlines.
			if( (str.gcount() > 1) && (buffer[ str.gcount() - 1] == '\n') && (buffer[ str.gcount() - 2] == '\n') ) {
				throw FileHasTwoTrailingNewlinesError( "(unknown)", count ) ;
			}
		}
		while( str ) ;
		return count ;
	}
	
	VCFFormatSNPDataSource::operator bool() const {
		return *m_stream_ptr ;
	}

	unsigned int VCFFormatSNPDataSource::number_of_samples() const {
		return m_column_names.size() - 8 ;
	}

	unsigned int VCFFormatSNPDataSource::total_number_of_snps() const {
		return m_number_of_lines ;
	}

	std::string VCFFormatSNPDataSource::get_source_spec() const {
		return m_spec ;
	}

	std::string VCFFormatSNPDataSource::get_summary( std::string const& prefix, std::size_t column_width ) const {
		return prefix + m_spec ;
	}

	void VCFFormatSNPDataSource::get_snp_identifying_data_impl( 
		IntegerSetter const& set_number_of_samples,
		StringSetter const& set_SNPID,
		StringSetter const& set_RSID,
		ChromosomeSetter const& set_chromosome,
		SNPPositionSetter const& set_SNP_position,
		AlleleSetter const& set_allele1,
		AlleleSetter const& set_allele2
	) {
		std::vector< std::string > ids ;
		GenomePosition pos ;
		std::vector< std::string > alleles ;
		//m_data_reader.read_id_data( set_value( pos ), set_value( ids ), set_value( alleles ), ignore(), ignore() ) ;
		//m_data_reader.read_info( ignore() ) ;
		
		set_number_of_samples( number_of_samples() ) ;
		set_chromosome( pos.chromosome() ) ;
		set_SNP_position( pos.position() ) ;

		// We currently only support SNPs
		if( alleles.size() == 2 && alleles[0].size() == 1 && alleles[1].size() == 1 ) {
			set_allele1( alleles[0][0] ) ;
			set_allele2( alleles[1][0] ) ;
		}
		else {
			set_allele1( '?' ) ;
			set_allele2( '?' ) ;
		}
		
		if( ids.size() == 0 ) {
			set_RSID( "?" ) ;
			set_SNPID( "?" ) ;
		}
		else if( ids.size() > 0 ) {
			set_RSID( ids[0] ) ;
		}
		if( ids.size() > 1 ) {
			set_SNPID( ids[1] ) ;
		}
	}

	void VCFFormatSNPDataSource::read_snp_probability_data_impl(
		GenotypeProbabilitySetter const& set_genotype_probabilities
	) {
		//m_call_data_reader( *m_stream_ptr, "GL", set_genotype_probabilities() ) ;
		assert(0) ;
	}

	void VCFFormatSNPDataSource::ignore_snp_probability_data_impl() {
		// m_call_data_reader( *m_stream_ptr, "", ignore() ) ;
		assert(0) ;
	}

	void VCFFormatSNPDataSource::reset_to_start_impl() {
		// seek back to the start and skip over metadata and column names.
		reset_stream() ;
	}
}
