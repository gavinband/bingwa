#include <set>
#include <string>
#include <cassert>
#include "genfile/GenomePosition.hpp"
#include "genfile/SNPIdentifyingDataTest.hpp"
#include "genfile/SNPPositionInRangeTest.hpp"
#include "genfile/string_utils.hpp"

namespace genfile {
	SNPPositionInRangeTest::SNPPositionInRangeTest( GenomePositionRange const range ):
		m_range( range )
	{
	}
		
	bool SNPPositionInRangeTest::operator()(
		std::string,
		std::string,
		GenomePosition position,
		char,
		char
	) const {
		return m_range.check_if_contains( position ) ;
	}
	
	std::string SNPPositionInRangeTest::display() const {
		return "position in " + string_utils::to_string( m_range ) ;
	}
}