#ifndef GENFILE_CHROMOSOME_IN_SET_TEST_HPP
#define GENFILE_CHROMOSOME_IN_SET_TEST_HPP

#include <set>
#include <string>
#include "genfile/Chromosome.hpp"
#include "genfile/SNPIdentifyingDataTest.hpp"

namespace genfile {
	struct ChromosomeInSetTest: public SNPIdentifyingDataTest
	{
		ChromosomeInSetTest( std::set< Chromosome > const& chromosomes ) ;
		bool operator()(
			std::string SNPID,
			std::string,
			GenomePosition,
			char,
			char
		) const ;
		
		std::string display() const ;
	private:
		std::set< Chromosome > m_chromosomes ;
	} ;
}

#endif
