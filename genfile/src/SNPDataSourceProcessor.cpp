//#include <boost/signal.hpp>
#include "genfile/SNPIdentifyingData.hpp"
#include "genfile/SNPDataSourceProcessor.hpp"
#include "genfile/get_set.hpp"

namespace genfile {
	SNPDataSourceProcessor::Callback::~Callback() {
	}
	
	SNPDataSourceProcessor::~SNPDataSourceProcessor() {}

	void SNPDataSourceProcessor::add_callback( Callback& callback ) {
		m_callbacks.push_back( &callback ) ;
	}

	void SNPDataSourceProcessor::call_begin_processing_snps( std::size_t const& number_of_samples, std::size_t const& number_of_snps ) const {
		for( std::size_t i = 0; i < m_callbacks.size(); ++i ) {
			m_callbacks[i]->begin_processing_snps( number_of_samples, number_of_snps ) ;
		}
	}
	
	void SNPDataSourceProcessor::call_processed_snp(  SNPIdentifyingData const& id_data, SingleSNPGenotypeProbabilities const& genotypes ) const {
		for( std::size_t i = 0; i < m_callbacks.size(); ++i ) {
			m_callbacks[i]->processed_snp( id_data, genotypes ) ;
		}
	}

	void SNPDataSourceProcessor::call_end_processing_snps() const {
		for( std::size_t i = 0; i < m_callbacks.size(); ++i ) {
			m_callbacks[i]->end_processing_snps() ;
		}
	}

	void SimpleSNPDataSourceProcessor::process( genfile::SNPDataSource& source, ProgressCallback progress_callback ) {
		SNPIdentifyingData id_data ;
		SingleSNPGenotypeProbabilities genotypes ;
		genotypes.resize( source.number_of_samples() ) ;
		
		call_begin_processing_snps( source.number_of_samples(), source.total_number_of_snps() ) ;

		while( source.read_snp(
				ignore(),
				set_value( id_data.SNPID() ),
				set_value( id_data.rsid() ),
				set_value( id_data.position().chromosome() ),
				set_value( id_data.position().position() ),
				set_value( id_data.first_allele() ),
				set_value( id_data.second_allele() ),
				set_genotypes( genotypes )
			)
		) {
			call_processed_snp( id_data, genotypes ) ;
			if( progress_callback ) {
				progress_callback( source.number_of_snps_read(), source.total_number_of_snps() ) ;
			}
		}

		call_end_processing_snps() ;
	}
}

