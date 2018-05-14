
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#ifndef EFFECTPARAMETERNAMEPACK_HPP
#define EFFECTPARAMETERNAMEPACK_HPP

#include <string>
#include <vector>

namespace bingwa {
	struct EffectParameterNamePack {
	public:
		static std::size_t const npos = std::string::npos ;
		EffectParameterNamePack() ;
		EffectParameterNamePack(
			std::vector< std::string > betas,
			std::vector< std::string > ses,
			std::vector< std::string > cov
		) ;

		EffectParameterNamePack( EffectParameterNamePack const& other ) ;
		EffectParameterNamePack& operator=( EffectParameterNamePack const& other ) ;

		std::size_t size() const ;
		std::size_t find( std::string const& ) const ;
		std::string parameter_name( std::size_t i ) const ;
		std::string se_name( std::size_t i ) const ;
		std::string wald_pvalue_name( std::size_t i ) const ;
		std::string covariance_name( std::size_t i, std::size_t j ) const ;
	
		/* rename a parameter, e.g. as beta_1:add, se_1 */
		void rename( std::size_t parameter_i, std::string const& newName, std::string const& newNumber, std::string const& identifier ) ;

		bool operator==( EffectParameterNamePack const& other ) const ;
		bool operator!=( EffectParameterNamePack const& other ) const ;
	
		std::string get_summary() const ;

	private:
		std::vector< std::string > m_betas ;
		std::vector< std::string > m_ses ;
		std::vector< std::string > m_cov ;
		
		std::size_t get_index_of_covariance( std::size_t i, std::size_t j ) const ;
	} ;
}

#endif
