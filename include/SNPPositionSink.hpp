#ifndef SNP_POSITION_SINK_HPP
#define SNP_POSITION_SINK_HPP

#include <string>
#include <iostream>
#include <memory>
#include <stddef.h>

#include "SNPDataSink.hpp"

struct SNPPositionSink: public genfile::SNPDataSink
{
	SNPPositionSink( std::string const& filename ) {
		setup( filename ) ;
	}
	
	operator bool() const { return *m_stream_ptr ; }

	void write_snp_impl(
		uint32_t,
		std::string,
		std::string,
		uint32_t SNP_position,
		char,
		char,
		GenotypeProbabilityGetter const&,
		GenotypeProbabilityGetter const&,
		GenotypeProbabilityGetter const&
	) {
		stream() << SNP_position << "\n" ;
	} ;
	
private:
	
	std::ostream& stream() { return *m_stream_ptr ; }
	
	void setup( std::string const& filename ) {
		m_stream_ptr = open_file_for_output( filename ) ;
	}
	
	std::auto_ptr< std::ostream > m_stream_ptr ;
} ;

#endif