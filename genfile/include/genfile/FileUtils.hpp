#ifndef GENFILE_FILE_UTILS_HPP
#define GENFILE_FILE_UTILS_HPP

#include <string>
#include <memory>

namespace genfile {
	enum CompressionType { e_NoCompression = 0, e_GzipCompression = 1 } ;
	bool filename_indicates_gen_format( std::string const& filename ) ;
	bool filename_indicates_bgen_format( std::string const& filename ) ;
	bool filename_indicates_gen_or_bgen_format( std::string const& filename ) ;
	CompressionType get_compression_type_indicated_by_filename( std::string const& filename ) ;
	Chromosome get_chromosome_indicated_by_filename( std::string const& filename ) ;

	std::pair< std::string, std::string > uniformise( std::string filename ) ;

	std::string get_gen_file_extension_if_present( std::string const& filename ) ;
	std::string strip_gen_file_extension_if_present( std::string const& filename ) ;

	std::string create_temporary_filename() ;

	std::auto_ptr< std::istream > open_text_file_for_input( std::string filename, CompressionType compression_type ) ;
	std::auto_ptr< std::ostream > open_text_file_for_output( std::string filename, CompressionType compression_type ) ;
	std::auto_ptr< std::istream > open_binary_file_for_input( std::string filename, CompressionType compression_type ) ;
	std::auto_ptr< std::ostream > open_binary_file_for_output( std::string filename, CompressionType compression_type ) ;
}

#endif