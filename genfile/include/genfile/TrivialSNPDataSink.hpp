#ifndef TrivialSNPDataSink_HPP
#define TrivialSNPDataSink_HPP

#include <iostream>
#include <string>
#include "SNPDataSink.hpp"

namespace genfile {
	// class TrivialSNPDataSink represents a SNPDataSink
	// which never fails and always discards its output.
	class TrivialSNPDataSink: public SNPDataSink
	{
	public:
		operator bool() const {
			return true ;
		}

	private:
		void write_snp_impl(
			uint32_t number_of_samples,
			std::string SNPID,
			std::string RSID,
			Chromosome chromosome,
			uint32_t SNP_position,
			char first_allele,
			char second_allele,
			GenotypeProbabilityGetter const& get_AA_probability,
			GenotypeProbabilityGetter const& get_AB_probability,
			GenotypeProbabilityGetter const& get_BB_probability
		) {
			return ;
		}
	} ;
}

#endif