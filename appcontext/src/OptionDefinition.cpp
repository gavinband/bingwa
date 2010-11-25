#include <string>
#include <sstream>
#include <cassert>
#include "appcontext/OptionDefinition.hpp"
#include "appcontext/OptionProcessor.hpp"

namespace appcontext {
	OptionDefinition::OptionDefinition()
		: m_description(""),
			m_is_required( false ),
			m_number_of_values_per_use( 0 ),
			m_lower_number_of_permitted_values( 0 ),
			m_upper_number_of_permitted_values( 1 ),
			m_default_values(),
			m_position( -1 )
	{}

	OptionDefinition::OptionDefinition( OptionDefinition const& other )
		: m_description( other.m_description ),
			m_is_required( other.m_is_required ),
			m_number_of_values_per_use( other.m_number_of_values_per_use ),
			m_lower_number_of_permitted_values( other.m_lower_number_of_permitted_values ),
			m_upper_number_of_permitted_values( other.m_upper_number_of_permitted_values ),
			m_default_values( other.m_default_values )
	{}

	void OptionDefinition::check_option_values( std::string const& option_name, std::vector< std::string > const& option_values ) const {

		// Check if number or values is too few...
		if( option_values.size() < m_lower_number_of_permitted_values ) {
			std::ostringstream ostr ;
			ostr << "Option \"" << option_name << "\" requires at least " << m_lower_number_of_permitted_values << " value" ;
			if( m_lower_number_of_permitted_values > 1 ) ostr << "s" ;
		
			throw OptionValueInvalidException( option_name, option_values, ostr.str() ) ;
		}

		// ...or too many.
		if( option_values.size() > m_upper_number_of_permitted_values ) {
			std::ostringstream ostr ;
			ostr << "Option \"" << option_name << "\" takes at most " << m_upper_number_of_permitted_values << " value" ;
			if( m_upper_number_of_permitted_values > 1 ) ostr << "s" ;

			throw OptionValueInvalidException( option_name, option_values, ostr.str() ) ;
		}

		std::vector< value_checker_t >::const_iterator
			i = m_value_checkers.begin(),
			end_i = m_value_checkers.end() ;
		for( ; i != end_i; ++i ) {
			(*i)( option_name, option_values ) ;
		}
	}
}