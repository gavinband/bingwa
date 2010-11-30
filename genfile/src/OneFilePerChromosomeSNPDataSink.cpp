#include <memory>
#include <map>
#include <string>
#include <boost/ptr_container/ptr_vector.hpp>
#include <boost/tuple/tuple.hpp>
#include "genfile/SNPDataSink.hpp"
#include "genfile/OneFilePerChromosomeSNPDataSink.hpp"

namespace genfile {
	namespace impl {
		boost::tuple< std::string, std::string, std::string > split_filename_template(
			std::string const& filename_template,
			std::string const& wildcard
		) {
			std::size_t pos = filename_template.find( wildcard ) ;
			if( pos == std::string::npos ) {
				return boost::make_tuple( filename_template, std::string( "" ), std::string( "" )) ;
			}
			else {
				return boost::make_tuple(
					filename_template.substr( 0, pos ),
					wildcard,
					filename_template.substr( pos + wildcard.size(), filename_template.size() )
				) ;
			}
		}
	}
	
	OneFilePerChromosomeSNPDataSink::UniquePtr OneFilePerChromosomeSNPDataSink::create(
		std::string const& filename_template,
		std::string const& free_data,
		std::string const& wildcard
	) {
		return UniquePtr(
			new OneFilePerChromosomeSNPDataSink(
				filename_template,
				free_data,
				wildcard
			)
		) ;
	}
	
	OneFilePerChromosomeSNPDataSink::OneFilePerChromosomeSNPDataSink(
		std::string const& filename_template,
		std::string const& free_data,
		std::string const& wildcard
	):
		m_filename_template( impl::split_filename_template( filename_template, wildcard ) ),
		m_free_data( free_data )
	{
	}
	
	OneFilePerChromosomeSNPDataSink::operator bool() const {
		for(
			std::map< Chromosome, SNPDataSink* >::const_iterator i = m_sinks.begin() ;
			i != m_sinks.end() ;
			++i
		) {
			if( !*(i->second) ) {
				return false ;
			}
		} 
		return true ;
	}
	void OneFilePerChromosomeSNPDataSink::write_snp_impl(
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
		get_sink_for_chromosome( chromosome ).write_snp(
			number_of_samples,
			SNPID,
			RSID,
			chromosome,
			SNP_position,
			first_allele,
			second_allele,
			get_AA_probability,
			get_AB_probability,
			get_BB_probability
		) ;
	}

	SNPDataSink& OneFilePerChromosomeSNPDataSink::get_sink_for_chromosome( Chromosome const& chromosome ) {
		std::map< Chromosome, SNPDataSink* >::const_iterator where = m_sinks.find( chromosome ) ;
		if( where != m_sinks.end() ) {
			// Already constructed this sink.
			return *where->second ;
		} else {
			if( m_filename_template.get<1>().size() == 0 ) {
				// no wildcard, everything goes to the same sink.
				assert( m_sink_storage.size() < 2 ) ;
				if( m_sink_storage.empty() ) {
					m_sink_storage.push_back( SNPDataSink::create( m_filename_template.get<0>(), m_free_data ).release() ) ;
				}
				m_sinks[ chromosome ] = &m_sink_storage[0] ;
				return m_sink_storage[0] ;
			}
			else {
				std::string const filename = m_filename_template.get<0>() + string_utils::to_string( chromosome ) + m_filename_template.get<2>() ;
				m_sink_storage.push_back( SNPDataSink::create( filename, m_free_data ).release() ) ;
				m_sinks[ chromosome ] = &m_sink_storage.back() ;
			}
			return m_sink_storage.back() ;
		}
	}
}
