#ifndef GENFILE_VCF_GET_SET_HPP
#define GENFILE_VCF_GET_SET_HPP

#include "genfile/VariantEntry.hpp"
#include "genfile/Error.hpp"
#include "genfile/SingleSNPGenotypeProbabilities.hpp"
#include "genfile/VariantDataReader.hpp"

namespace genfile {
	namespace vcf {
		struct VectorSetter: public VariantDataReader::PerSampleSetter {
		public:
			VectorSetter( std::vector< std::vector< Entry > >& data ):
				m_data( data )
			{}

			void set_number_of_samples( std::size_t n ) { m_data.resize( n ) ; }
			void set_sample( std::size_t n ) { assert( n < m_data.size() ) ; m_sample = n ; }
			void set_number_of_entries( std::size_t n ) { m_data[ m_sample ].resize( n ) ; m_entry_i = 0 ; }

		private:
			template< typename T >
			void set( T value ) {
				assert( m_entry_i < m_data[ m_sample ].size() ) ;
				m_data[ m_sample ][ m_entry_i++ ] = value ;
			}
		public:
			void operator()( MissingValue const value ) { set( value ) ; }
			void operator()( std::string& value ) { set( value ) ; }
			void operator()( Integer const value ) { set( value ) ; }
			void operator()( double const value ) { set( value ) ; }

		private:
			std::vector< std::vector< Entry > >& m_data ;
			std::size_t m_number_of_samples ;
			std::size_t m_sample ;
			std::size_t m_entry_i ;
		} ;
		
		struct GenotypeSetterBase: public VariantDataReader::PerSampleSetter
		{
			virtual ~GenotypeSetterBase() throw() ;

			virtual void set_number_of_samples( std::size_t n ) ;
			virtual void set_sample( std::size_t n ) ;
			virtual void set_number_of_entries( std::size_t n ) ;
			virtual void operator()( MissingValue const value ) ;
			virtual void operator()( Integer const value ) ;
			virtual void operator()( double const value ) ;

		protected:
			virtual void set( std::size_t, double, double, double ) = 0 ;

		private:
			std::size_t m_number_of_samples ;
			std::size_t m_sample ;
			std::size_t m_number_of_entries ;
			std::size_t m_entry_i ;
			double m_store[3] ;
			double m_A, m_B ;
			bool m_missing ;

			void set() ;

			template< typename T >
			void store( T const value ) {
				if( m_number_of_entries == 2 ) {
					// Treat as two calls.
					assert( value == 0 || value == 1 ) ;
					m_A += ( value == 0 ) ? 1 : 0 ;
					m_B += ( value == 0 ) ? 0 : 1 ;
				}
				else if( m_number_of_entries == 3 || m_number_of_entries == 4 ) {
					// treat as probabilities.  Ignore the fourth probability, which we interpret as NULL call.
					if( m_entry_i < 3 ) {
						m_store[ m_entry_i ] = value ;
					}
				}
				else {
					assert(0) ;
				}
				++m_entry_i ;
				if( m_entry_i == m_number_of_entries ) {
					set() ;
				}
			}
		} ;

		template< typename Setter >
		struct GenotypeSetter: public GenotypeSetterBase
		{
			GenotypeSetter( Setter const& setter ): m_setter( setter ) {}
			void set( std::size_t sample_i, double AA, double AB, double BB ) {
				m_setter( sample_i, AA, AB, BB ) ;
			}
		private:
			Setter const& m_setter ;
		} ;
		
		template<>
		struct GenotypeSetter< SingleSNPGenotypeProbabilities >: public GenotypeSetterBase
		{
			GenotypeSetter( SingleSNPGenotypeProbabilities& result ) ;
			void set_number_of_samples( std::size_t n ) ;
			void set( std::size_t sample_i, double AA, double AB, double BB ) ;
		private:
			SingleSNPGenotypeProbabilities& m_result ;
		} ;

		template<>
		struct GenotypeSetter< std::vector< double > >: public GenotypeSetterBase
		{
			GenotypeSetter( std::vector< double >& result ) ;
			void set_number_of_samples( std::size_t n ) ;
			void set( std::size_t sample_i, double AA, double AB, double BB ) ;
		private:
			std::vector< double >& m_result ;
		} ;

		template<>
		struct GenotypeSetter< std::vector< VariantEntry > >: public GenotypeSetterBase
		{
			GenotypeSetter( std::vector< VariantEntry >& result, double threshhold ) ;
			void set_number_of_samples( std::size_t n ) ;
			void set( std::size_t sample_i, double AA, double AB, double BB ) ;
		private:
			std::vector< VariantEntry >& m_result ;
			double const m_threshhold ;
		} ;
		
		template< typename Matrix >
		struct MatrixSetter: public VariantDataReader::PerSampleSetter
		{
			MatrixSetter( Matrix& result ): m_result( result ) {}
			void set_number_of_samples( std::size_t n ) { m_number_of_samples = n ; }
			void set_sample( std::size_t n ) {
				assert( n < m_number_of_samples ) ; 
				m_sample = n ;
				m_entry_i = 0 ;
			}
			void set_number_of_entries( std::size_t n ) {
				if( m_sample == 0 ) {
					m_result.resize( n, m_number_of_samples ) ;
				}
				else if( !n == m_result.rows() ) {
					throw BadArgumentError( "genfile::vcf::MatrixSetter::set_number_of_entries()", "n" ) ;
				}
			}
			void operator()( MissingValue const value ) {
				m_result( m_entry_i++, m_sample ) = std::numeric_limits< double >::quiet_NaN() ;
			}
			void operator()( double const value ) {
				m_result( m_entry_i++, m_sample ) = value ;
			}

		private:
			Matrix& m_result ;
			std::size_t m_number_of_samples ;
			std::size_t m_sample ;
			std::size_t m_number_of_entries ;
			std::size_t m_entry_i ;
		} ;
	}
}

#endif
