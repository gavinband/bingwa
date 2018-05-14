
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#include <string>
#include <vector>
#include <cassert>
#include "genfile/string_utils/string_utils.hpp"
#include "genfile/string_utils/slice.hpp"
#include "bingwa/EffectParameterNamePack.hpp"

namespace bingwa {
	EffectParameterNamePack::EffectParameterNamePack() {}
	EffectParameterNamePack::EffectParameterNamePack(
		std::vector< std::string > betas,
		std::vector< std::string > ses,
		std::vector< std::string > cov
	):
		m_betas( betas ),
		m_ses( ses ),
		m_cov( cov )
	{
		assert( m_ses.size() == m_betas.size() ) ;
		assert( m_cov.size() == ( m_betas.size() * ( m_betas.size() - 1 )) / 2 ) ;
	}

	EffectParameterNamePack::EffectParameterNamePack( EffectParameterNamePack const& other ):
		m_betas( other.m_betas ),
		m_ses( other.m_ses ),
		m_cov( other.m_cov )
	{
	}

	EffectParameterNamePack& EffectParameterNamePack::operator=( EffectParameterNamePack const& other ) {
		m_betas = other.m_betas ;
		m_ses = other.m_ses ;
		m_cov = other.m_cov ;
		return *this ;
	}

	std::size_t EffectParameterNamePack::size() const {
		return m_betas.size() ;
	}

	std::size_t EffectParameterNamePack::find( std::string const& name ) const {
		std::vector< std::string >::const_iterator where = std::find( m_betas.begin(), m_betas.end(), name ) ;
		if( where == m_betas.end() ) {
			return npos ;
		} else {
			return std::size_t( where - m_betas.begin() ) ;
		}
	}

	std::string EffectParameterNamePack::parameter_name( std::size_t i ) const {
		return m_betas[i] ;
	}

	std::string EffectParameterNamePack::se_name( std::size_t i ) const {
		return m_ses[i] ;
	}

	std::string EffectParameterNamePack::wald_pvalue_name( std::size_t i ) const {
		std::string result = m_ses[i] ;
		std::size_t pos = result.find( "se_" ) ;
		assert( pos != std::string::npos ) ;
		result.replace( pos, 3, "wald_pvalue_" ) ;
		return result ;
	}

	std::string EffectParameterNamePack::covariance_name( std::size_t i, std::size_t j ) const {
		return m_cov[ get_index_of_covariance(i,j) ] ;
	}

	std::size_t EffectParameterNamePack::get_index_of_covariance( std::size_t i, std::size_t j ) const {
		assert( j > i ) ;
		// index in cov skips the lower diagonal.
		// Lower diagonal has ((i+1)*(i+2)/2) entries since i is a 0-based index.
		// index in full array would be i*N+j
		std::size_t const N = m_betas.size() ;
		return i*N + j - ( (i+1)*(i+2)/2 ) ;
	}

	void EffectParameterNamePack::rename(
		std::size_t parameter_i,
		std::string const& newName,
		std::string const& newSubscript,
		std::string const& identifier
	) {
		using genfile::string_utils::to_string ;
		using genfile::string_utils::slice ;
		std::size_t const N = m_betas.size() ;
		assert( parameter_i < N ) ;
		m_betas[parameter_i] = newName + "_" + newSubscript + ":" + identifier ;
		m_ses[parameter_i] = "se_" + newSubscript ;
		// Rename covariance parameters
		// They are expected to be of the form <name>_i,j
		{
			std::size_t j = parameter_i ;
			for( std::size_t i = 0; i < parameter_i; ++i ) {
				std::size_t index = get_index_of_covariance(i,j) ;
				std::vector< slice > elts = slice( m_cov[index] ).split( "_" ) ;
				assert( elts.size() == 2 ) ;
				std::vector< slice > elts2 = elts[1].split( "," ) ;
				assert( elts2.size() == 2 ) ;
				m_cov[index] = elts[0] + "_" + elts2[0] + "," + newSubscript ;
			}
		}
		{
			std::size_t i = parameter_i ;
			for( std::size_t j = parameter_i+1; j < N; ++j ) {
				std::size_t index = get_index_of_covariance(i,j) ;
				std::vector< slice > elts = slice( m_cov[index] ).split( "_" ) ;
				assert( elts.size() == 2 ) ;
				std::vector< slice > elts2 = elts[1].split( "," ) ;
				assert( elts2.size() == 2 ) ;
				m_cov[index] = elts[0] + "_" + newSubscript + "," + elts2[1] ;
			}
		}
	}

	bool EffectParameterNamePack::operator==( EffectParameterNamePack const& other ) const {
		return m_betas == other.m_betas && m_ses == other.m_ses && m_cov == other.m_cov ;
	}

	bool EffectParameterNamePack::operator!=( EffectParameterNamePack const& other ) const {
		return !operator==( other ) ;
	}

	std::string EffectParameterNamePack::get_summary() const {
		using genfile::string_utils::join ;
		std::ostringstream ostr ;
		ostr << "(" << join( m_betas, "," ) << "; " << join( m_ses, "," ) << "; " << join( m_cov, ", " ) << ")" ;
		return ostr.str() ;
	}
}
