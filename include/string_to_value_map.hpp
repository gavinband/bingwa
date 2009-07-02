#ifndef __GTOOL_STRING_TO_VALUE_MAP_HPP
#define __GTOOL_STRING_TO_VALUE_MAP_HPP


#include <string>


// Base class for objects which map strings to values.
// The values can be returned either as string or double values.
class string_to_value_map
{
public:
	
	template< typename T >
	T get_value( std::string const& name ) const ;

protected:

	virtual double get_double_value( std::string const& name ) const = 0 ;
	virtual std::string get_string_value( std::string const& name ) const = 0 ;
} ;


#endif