#include <iostream>
#include <string>
#include "genfile/snp_data_utils.hpp"
#include "genfile/SNPDataSource.hpp"
#include "genfile/bgen.hpp"
#include "genfile/BGenFileSNPDataSource.hpp"

namespace genfile {
	BGenFileSNPDataSource::BGenFileSNPDataSource( std::string const& filename )
		: m_filename( filename )
	{
		setup( filename, get_compression_type_indicated_by_filename( filename )) ;
	}

	BGenFileSNPDataSource::BGenFileSNPDataSource( std::string const& filename, CompressionType compression_type )
		: m_filename( filename )
	{
		setup( filename, compression_type ) ;
	}

	void BGenFileSNPDataSource::reset_to_start_impl() {
		stream().clear() ;
		stream().seekg(0) ;

		// read offset again and skip to first SNP data block
		bgen::uint32_t offset ;
		bgen::read_offset( (*m_stream_ptr), &offset ) ;
		m_stream_ptr->ignore( offset ) ;
	}

	void BGenFileSNPDataSource::read_snp_identifying_data_impl( 
		uint32_t* number_of_samples,
		std::string* SNPID,
		std::string* RSID,
		Chromosome* chromosome,
		uint32_t* SNP_position,
		char* allele1,
		char* allele2
	) {
		unsigned char chr ;
		bgen::impl::read_snp_identifying_data( stream(), number_of_samples, SNPID, RSID, &chr, SNP_position, allele1, allele2 ) ;
		*chromosome = Chromosome( chr ) ;
	}

	void BGenFileSNPDataSource::read_snp_probability_data_impl(
		GenotypeProbabilitySetter const& set_genotype_probabilities
	) {
		if( m_flags & bgen::e_CompressedSNPBlocks ) {
			bgen::impl::read_compressed_snp_probability_data( stream(), number_of_samples(), set_genotype_probabilities ) ;
		}
		else {
			bgen::impl::read_snp_probability_data( stream(), number_of_samples(), set_genotype_probabilities ) ;
		}
	}

	void BGenFileSNPDataSource::ignore_snp_probability_data_impl() {
		if( m_flags & bgen::e_CompressedSNPBlocks ) {
			bgen::impl::read_compressed_snp_probability_data( stream(), number_of_samples(), ignore() ) ;
		}
		else {
			bgen::impl::read_snp_probability_data( stream(), number_of_samples(), ignore() ) ;
		}
	}

	void BGenFileSNPDataSource::setup( std::string const& filename, CompressionType compression_type ) {
		m_stream_ptr = open_binary_file_for_input( filename, compression_type ) ;
		bgen::uint32_t offset ;
		bgen::read_offset( (*m_stream_ptr), &offset ) ;
		bgen::uint32_t header_size = read_header_data() ;

		if( offset < header_size ) {
			throw FileStructureInvalidError() ;
		}
		// skip any remaining bytes before the first snp block
		m_stream_ptr->ignore( offset - header_size ) ;
		if( !*m_stream_ptr ) {
			throw FileStructureInvalidError() ;
		}
	}

	bgen::uint32_t BGenFileSNPDataSource::read_header_data() {
		bgen::uint32_t header_size ;

		bgen::read_header_block(
			(*m_stream_ptr),
			set_value( header_size ),
			set_value( m_total_number_of_snps ),
			set_value( m_number_of_samples ),
			ignore(),
			set_value( m_flags )
		) ;

		return header_size ;
	}
}	
