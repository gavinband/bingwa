
//			Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//	  (See accompanying file LICENSE_1_0.txt or copy at
//			http://www.boost.org/LICENSE_1_0.txt)

#ifndef GENBIN_REFERENCE_IMPLEMENTATION_HPP
#define GENBIN_REFERENCE_IMPLEMENTATION_HPP

#include <iostream>
#include <vector>
#include <cassert>
#include <cmath>
#include <stdint.h>
#include "../../config.hpp"
#ifdef HAVE_ZLIB
#include <sstream>
#include <zlib.h>
#endif
#include "genfile/endianness_utils.hpp"
#include "genfile/snp_data_utils.hpp"
#include "genfile/get_set.hpp"
#include "genfile/Error.hpp"
#include "genfile/zlib.hpp"
#include "genfile/bgen/impl.hpp"

/*
* This file contains a reference implementation of the BGEN file format
* specification described at:
* http://www.well.ox.ac.uk/~gav/bgen_format/bgen_format.html
*
* To use this file you will also need "endianness_utils.hpp".
*
*/

namespace genfile {
	namespace bgen {
		struct BGenError: public virtual std::exception {
			~BGenError() throw() {}
			char const* what() const throw() { return "BGenError" ; }
		} ;
		typedef ::uint32_t uint32_t ;
		typedef ::uint16_t uint16_t ;

		// v1.1 definitions were:
		// enum FlagsType { e_NoFlags = 0, e_CompressedSNPBlocks = 0x1, e_LongIds = 0x4 } ;
		// v1.2 definitions:
		enum FlagsType { e_NoFlags = 0, e_CompressedSNPBlocks = 0x1, e_Layout = 0x3C } ;
		enum Layout { e_v10Layout = 0x0, e_v11Layout = 0x4, e_v12Layout = 0x8 } ;

		// Read the offset from the start of the stream.
		void read_offset( std::istream& iStream, uint32_t* offset ) ;
		// Write an offset value to the stream.
		void write_offset( std::ostream& oStream, uint32_t const offset ) ;

		// Get the size of the header block with the given free data
		std::size_t get_header_block_size(
			std::string const& free_data
		) ;

		// Read a header block from the supplied stream, using
		// the supplied setter objects (which can be function pointers) to return
		// the data to the user.
		template<
			typename HeaderSizeSetter,
			typename NumberOfSNPBlocksSetter,
			typename NumberOfSamplesSetter,
			typename FlagsSetter,
			typename FreeDataSetter,
			typename VersionSetter
		>
		void read_header_block(
			std::istream& aStream,
			HeaderSizeSetter set_header_size,
			NumberOfSNPBlocksSetter set_number_of_snp_blocks,
			NumberOfSamplesSetter set_number_of_samples,
			FreeDataSetter set_free_data,
			FlagsSetter set_flags,
			VersionSetter set_version
		) {
			uint32_t
				header_size = 0,
				number_of_snp_blocks = 0,
				number_of_samples = 0,
				flags = 0 ;

			char version[4] ;

			std::size_t fixed_data_size = 20 ;
			std::vector<char> free_data ;

			genfile::read_little_endian_integer( aStream, &header_size ) ;
			assert( header_size >= fixed_data_size ) ;
			genfile::read_little_endian_integer( aStream, &number_of_snp_blocks ) ;
			genfile::read_little_endian_integer( aStream, &number_of_samples ) ;
			aStream.read( &(version[0]), 4 ) ;
			free_data.resize( header_size - fixed_data_size ) ;
			aStream.read( &(free_data[0]), free_data.size() ) ;
			genfile::read_little_endian_integer( aStream, &flags ) ;

			if( aStream ) {
				set_header_size( header_size ) ;
				set_number_of_snp_blocks( number_of_snp_blocks ) ;
				set_number_of_samples( number_of_samples ) ;
				set_free_data( std::string( free_data.begin(), free_data.end() )) ;
				set_flags( flags ) ;
				set_version( std::string( version, version + 4 )) ;
			} else {
				throw BGenError() ;
			}
		}

		// Write a bgen header block to the supplied stream.
		void write_header_block(
			std::ostream& aStream,
			uint32_t number_of_snp_blocks,
			uint32_t number_of_samples,
			std::string const& free_data,
			uint32_t flags
		) ;

		// Read a sample identifier block from the given stream
		template<
			typename NumberOfSamplesSetter,
			typename SampleIdentifierSetter
		>
		void read_sample_block(
			std::istream& aStream,
			NumberOfSamplesSetter set_number_of_samples,
			SampleIdentifierSetter set_identifier
		) {
			uint32_t block_size = 0 ;
			uint32_t number_of_samples = 0 ;
			genfile::read_little_endian_integer( aStream, &block_size ) ;
			genfile::read_little_endian_integer( aStream, &number_of_samples ) ;
			set_number_of_samples( number_of_samples ) ;

			uint16_t identifier_size ;
			std::string identifier ;

			for( std::size_t i = 0; i < number_of_samples; ++i ) {
				read_length_followed_by_data( aStream, &identifier_size, &identifier ) ;
				if( aStream ) {
					set_identifier( i, identifier ) ;
				} else {
					throw BGenError() ;
				}
			}
		}
		
		// Write a sample identifier block to the given stream
		template<
			typename SampleIdentifierGetter
		>
		void write_sample_block(
			std::ostream& aStream,
			uint32_t number_of_samples,
			SampleIdentifierGetter get_sample_identifier
		) {
			uint32_t block_size = 8 ;
			for( uint32_t i = 0; i < number_of_samples; ++i ) {
				block_size += 2 + get_sample_identifier(i).size() ;
			}
			write_little_endian_integer( aStream, block_size ) ;
			write_little_endian_integer( aStream, number_of_samples ) ;
			for( uint32_t i = 0; i < number_of_samples; ++i ) {
				std::string const& identifier = get_sample_identifier(i) ;
				assert( identifier.size() <= std::size_t( std::numeric_limits< uint16_t >::max() ) ) ;
				uint16_t const id_size = identifier.size() ;
				write_little_endian_integer( aStream, id_size ) ;
				write_length_followed_by_data( aStream, id_size, identifier ) ;
			}
		}

		// Read identifying data fields for the next variant in the file
		void read_snp_identifying_data(
			std::istream& aStream,
			uint32_t const flags,
			uint32_t* number_of_samples,
			std::string* SNPID,
			std::string* RSID,
			unsigned char* chromosome,
			uint32_t* SNP_position,
			std::string* first_allele,
			std::string* second_allele
		) ;

		// Write identifying data fields for the given variant.
		void write_snp_identifying_data(
			std::ostream& aStream,
			uint32_t const flags,
			uint32_t number_of_samples,
			unsigned char max_id_size,
			std::string SNPID,
			std::string RSID,
			unsigned char chromosome,
			uint32_t SNP_position,
			std::string first_allele,
			std::string second_allele
		) ;

		// Skip over probability data for the current variant
		void ignore_snp_probability_data(
			std::istream& aStream,
			uint32_t const flags,
			uint32_t number_of_samples
		) ;
		
		// Read the probability data for the current variant
		template<
			typename GenotypeProbabilitySetter
		>
		void read_snp_probability_data(
			std::istream& aStream,
			uint32_t const flags,
			uint32_t number_of_samples,
			GenotypeProbabilitySetter set_genotype_probabilities,
			std::vector< char >* buffer1,
			std::vector< char >* buffer2
		) {
			uint32_t const data_size = 6 * number_of_samples ;
			buffer1->resize( data_size ) ;
			if( flags & bgen::e_CompressedSNPBlocks ) {
				uint32_t compressed_data_size ;
				genfile::read_little_endian_integer( aStream, &compressed_data_size ) ;
				buffer2->resize( compressed_data_size ) ;
				aStream.read( &(*buffer2)[0], compressed_data_size ) ;
				zlib_uncompress( *buffer2, buffer1 ) ;
				assert( buffer1->size() == data_size ) ;
			}
			else {
				aStream.read( &(*buffer1)[0], data_size ) ;
			}
			impl::read_uncompressed_snp_probability_data(
				&(*buffer1)[0],
				&(*buffer1)[0] + data_size,
				impl::get_probability_conversion_factor( flags ),
				number_of_samples,
				set_genotype_probabilities
			) ;
		}

		// Write per-variant probability data to the given stream.
		// Genotype probabilities must be supplied by the given GenotypeProbabilityGetter
		// objects, which must be callable as
		// - get_AA_probability( index )
		// - get_AB_probability( index )
		// - get_BB_probability( index )
		// where index is the index of the individual in the SNP block.
		template< typename GenotypeProbabilityGetter >
		void write_snp_probability_data(
			std::ostream& aStream,
			uint32_t const flags,
			uint32_t const number_of_samples,
			GenotypeProbabilityGetter get_AA_probability,
			GenotypeProbabilityGetter get_AB_probability,
			GenotypeProbabilityGetter get_BB_probability
		) {
			if( flags & e_CompressedSNPBlocks ) {
				impl::write_compressed_snp_probability_data(
					aStream,
					flags,
					number_of_samples,
					get_AA_probability,
					get_AB_probability,
					get_BB_probability
				) ;
			}
			else {
				impl::write_uncompressed_snp_probability_data(
					aStream,
					flags,
					number_of_samples,
					get_AA_probability,
					get_AB_probability,
					get_BB_probability
				) ;
			}
		}

	}
}

#endif
