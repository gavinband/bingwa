#ifndef GENFILE_SHAPEIT_HAPLOTYPES_SNP_DATA_SOURCE_HPP
#define GENFILE_SHAPEIT_HAPLOTYPES_SNP_DATA_SOURCE_HPP

#include <string>
#include <vector>
#include "genfile/SNPDataSource.hpp"
#include "genfile/string_utils/slice.hpp"

namespace genfile {
	// Read genotypes from a ShapeIT-style output file.
	// It is assumed that haplotypes come in pairs, corresponding to haplotypes from the same individual.
	class ShapeITHaplotypesSNPDataSource: public IdentifyingDataCachingSNPDataSource
	{
	public:
		ShapeITHaplotypesSNPDataSource( std::string const& haplotypes_filename, Chromosome chromosome ) ;
		ShapeITHaplotypesSNPDataSource( std::string const& haplotypes_filename, Chromosome chromosome, CompressionType compression_type ) ;

		unsigned int number_of_samples() const { return m_number_of_samples ; }
		OptionalSnpCount total_number_of_snps() const { return OptionalSnpCount() ; }
		
		operator bool() const { return m_good ; }
		std::istream& stream() { return *m_stream_ptr ; }
		std::istream const& stream() const { return *m_stream_ptr ; }

		Chromosome chromosome() const { return m_chromosome ; }
		std::string get_source_spec() const { return m_haplotypes_filename ; }

	private:
		void reset_to_start_impl() ;
		
		void read_snp_identifying_data_impl( 
			uint32_t* number_of_samples,
			std::string* SNPID,
			std::string* RSID,
			Chromosome* chromosome,
			uint32_t* SNP_position,
			std::string* allele1,
			std::string* allele2
		) ;

		VariantDataReader::UniquePtr read_variant_data_impl() ;

		void ignore_snp_probability_data_impl() ;

	private:
		std::string const m_haplotypes_filename ;
		CompressionType m_compression_type ;
		unsigned int m_number_of_samples ;
		std::auto_ptr< std::istream > m_stream_ptr ;
		bool m_have_chromosome_column ;
		Chromosome m_chromosome ;
		bool m_good ;
		std::string m_current_line ;

		void setup( std::string const& filename, CompressionType compression_type ) ;
		void setup( std::auto_ptr< std::istream > stream_ptr ) ;
		void count_samples() ;
	} ;
}

#endif