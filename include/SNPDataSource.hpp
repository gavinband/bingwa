#ifndef SNPDATAPROVIDER_HPP
#define SNPDATAPROVIDER_HPP

#include <iostream>
#include <string>
#include "snp_data_utils.hpp"
#include "gen.hpp"
#include "bgen.hpp"

namespace genfile {
	struct SNPDataSourceError: public SNPDataError { char const* what() const throw() { return "SNPDataSourceError" ; } } ;
	struct FileIsInvalidError: public SNPDataSourceError { char const* what() const throw() { return "FileIsInvalidError" ; } } ;
	struct FileHasTrailingDataAfterLastSNP: public FileIsInvalidError { char const* what() const throw() { return "FileHasTrailingDataAfterLastSNP" ; } } ;
	struct FileContainsSNPsOfDifferentSizes: public FileIsInvalidError { char const* what() const throw() { return "FileContainsSNPsOfDifferentSizes" ; } } ;
	struct FileStructureInvalidError: public FileIsInvalidError { char const* what() const throw() { return "FileStructureInvalidError" ; } } ;

	// Base class for classes which provide SNP assay data, one snp at a time.
	// After the class is constructed, the intention is that
	// 1. the number_of_samples() and total_number_of_snps() functions return information
	// reflecting the data in the file or source, and
	// 2. The stream accessible through stream() is readable and pointing to the first snp in the file.
	class SNPDataSource
	{
	public:

		SNPDataSource() {} ;
		virtual ~SNPDataSource() {} ;

		// Factory functions
		static std::auto_ptr< SNPDataSource > create( std::string const& filename ) ;
		static std::auto_ptr< SNPDataSource > create( std::string const& filename, bool file_is_gzipped ) ;
		static std::auto_ptr< SNPDataSource > create( std::vector< std::string > const& filenames ) ;

		virtual unsigned int number_of_samples() const = 0;
		virtual unsigned int total_number_of_snps() const = 0 ;
		virtual operator bool() const { return stream() ; }
		virtual std::istream& stream() = 0 ;
		virtual std::istream const & stream() const = 0 ;
		enum FormatType { e_GenFormat = 0, e_BGenFormat = 1 } ;
		virtual FormatType format() const = 0;

		// Function read_snp().
		// Ideally this would also be a virtual member function.
		// However, a template member function can't be virtual.
		// Therefore, we dispatch to the correct implementation using the format()
		// and stream() members which implementations must provide.
		template<
			typename IntegerSetter,
			typename StringSetter,
			typename AlleleSetter,
			typename SNPPositionSetter,
			typename GenotypeProbabilitySetter
		>
		SNPDataSource& read_snp(
			IntegerSetter set_number_of_samples,
			StringSetter set_SNPID,
			StringSetter set_RSID,
			SNPPositionSetter set_SNP_position,
			AlleleSetter set_allele1,
			AlleleSetter set_allele2,
			GenotypeProbabilitySetter set_genotype_probabilities
		) {
			pre_read_snp() ;

			if( format() == e_GenFormat ) {
				gen::read_snp_block( stream(), set_number_of_samples, set_SNPID, set_RSID, set_SNP_position, set_allele1, set_allele2, set_genotype_probabilities ) ;
			}
			else if( format() == e_BGenFormat ){
				bgen::read_snp_block( stream(), set_number_of_samples, set_SNPID, set_RSID, set_SNP_position, set_allele1, set_allele2, set_genotype_probabilities ) ;
			}
			else {
				assert(0) ; // invalid format type.
			}
		
			post_read_snp() ;

			return *this ;
		};
	
	protected:
		virtual void pre_read_snp() {} ;
		virtual void post_read_snp() {} ;
	private:
	
		SNPDataSource( SNPDataSource const& other ) ;
		SNPDataSource& operator=( SNPDataSource const& other ) ;
	} ;
}

#endif