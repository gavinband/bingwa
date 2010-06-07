#include <vector>
#include <string>
#include <iostream>
#include <sstream>
#include <cassert>
#include "SampleRow.hpp"
#include "string_utils/string_utils.hpp"
#include "Whitespace.hpp"
#include "genfile/CohortIndividualSource.hpp"

bool check_char_is_space( char c ) { return c == ' ' ; }
bool check_char_is_newline( char c ) { return c == '\n' ; }

SampleRow::SampleRow() {}

SampleRow::SampleRow( SampleRow const& other ) :
	m_id1( other.m_id1 ),
	m_id2( other.m_id2 ),
	m_further_data( other.m_further_data ),
	m_column_headings( other.m_column_headings ),
	m_column_types( other.m_column_types )
{
}

void SampleRow::reset( std::vector< genfile::CohortIndividualSource::SingleColumnSpec > const& column_spec ) {
	m_column_headings.resize( column_spec.size() ) ;
	m_column_types.resize( column_spec.size() ) ;
	for( std::size_t i = 0; i < column_spec.size(); ++i ) {
		m_column_headings[i] = column_spec[i].name() ;
		std::ostringstream ostr ;
		ostr << column_spec[i].type() ;
		assert( ostr.str().size() == 1 ) ;
		m_column_types[i] = ostr.str()[0] ;
	}
	m_id1 = m_id2 = "" ;
	m_further_data.clear() ;

	assert( m_column_headings.size() >= 3 ) ;
	assert( m_column_types.size() == m_column_headings.size()) ;
	assert( m_column_headings[2] == "missing" ) ;
	assert( m_column_types[0] == '0' ) ;
	assert( m_column_types[1] == '0' ) ;
	assert( m_column_types[2] == '0' ) ;
}

SampleRow::Entry SampleRow::further_data( std::string const& column_heading ) const {
	FurtherData::const_iterator further_data_i = m_further_data.find( column_heading ) ;
	assert( further_data_i != m_further_data.end() ) ;
	return further_data_i->second ;
}

void SampleRow::add_column( std::string const& heading, char type, Entry value ) {
	if( !has_column( heading )) {
		m_column_headings.push_back( heading ) ;
		m_column_types.push_back( type ) ;
		m_further_data[heading] = value ;
	}
}

bool SampleRow::has_column( std::string const& heading ) const {
	FurtherData::const_iterator where = m_further_data.find( heading ) ;
	return ( where != m_further_data.end()) ;
}

bool SampleRow::has_value( std::string const& name ) const {
	return SampleRow::has_column( name ) ;
}

void SampleRow::set_value( std::string const& heading, Entry value ) {
	assert( has_column( heading )) ;
	m_further_data[ heading ] = value ;
}

double SampleRow::get_double_value( std::string const& heading ) const {
	FurtherData::const_iterator where = m_further_data.find( heading ) ;
	assert( where != m_further_data.end() ) ;
	return where->second.as< double >() ;
}

std::string SampleRow::get_string_value( std::string const& name ) const {
	if( name == "ID1" ) {
		return m_id1 ;
	}
	else if( name == "ID2" ) {
		return m_id2 ;
	}
	else {
		double value = get_double_value( name ) ;
		std::ostringstream ostr ;
		ostr << value ;
		return ostr.str() ;
	}
}

void SampleRow::read_ith_sample_from_source( std::size_t sample_i, genfile::CohortIndividualSource const& source ) {
	assert( sample_i < source.get_number_of_individuals() ) ;
	std::vector< genfile::CohortIndividualSource::SingleColumnSpec > column_spec = source.get_column_spec() ;
	reset( column_spec ) ;
	for( std::size_t i = 0; i < source.get_number_of_columns(); ++i ) {
		set_value( column_spec[i].name(), source.get_entry( sample_i, column_spec[i].name() )) ;
	}
} 

std::ostream& operator<<( std::ostream& aStream, SampleRow const& row ) {
	aStream << row.ID1() << " " << row.ID2() ;
	for( std::size_t i = 2 ; i < row.column_headings().size(); ++i ) {
		if( i > 0 )
			aStream << " " ;
		aStream << row.further_data( row.column_headings()[i] );
	}
	return aStream << "\n" ;
}
