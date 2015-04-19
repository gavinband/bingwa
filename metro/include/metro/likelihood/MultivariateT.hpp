
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#ifndef METRO_LIKELIHOOD_MULTIVARIATE_T_HPP
#define METRO_LIKELIHOOD_MULTIVARIATE_T_HPP

#include <cmath>
#include <boost/math/special_functions/gamma.hpp>
#include <Eigen/Core>
#include <Eigen/LU>
#include <Eigen/Cholesky>
#include "metro/LogLikelihood.hpp"
#include "metro/DataRange.hpp"
#include "metro/DataSubset.hpp"

// #define DEBUG_MULTIVARIATE_T 1

namespace metro {
	namespace likelihood {
		template< typename Scalar, typename Vector, typename Matrix >
		struct MultivariateT: public metro::LogLikelihood< Scalar, Vector, Matrix > {
		public:
			typedef typename Vector::SegmentReturnType Segment ;
			typedef typename Vector::ConstSegmentReturnType ConstSegment ;
		public:
			MultivariateT( double const degrees_of_freedom ):
				m_pi( 3.141592653589793238462643383279502884 ),
				m_nu( degrees_of_freedom ),
				m_data( 0 ),
				m_kappa( 0 )
			{
			}
			
			MultivariateT( Matrix const& data, double const degrees_of_freedom ):
				m_pi( 3.141592653589793238462643383279502884 ),
				m_nu( degrees_of_freedom ),
				m_data( 0 ),
				m_kappa( 0 )
			{
				set_data( data ) ;
			}

			MultivariateT( Matrix const& data, DataSubset const& data_range, double const degrees_of_freedom ):
				m_pi( 3.141592653589793238462643383279502884 ),
				m_nu( degrees_of_freedom ),
				m_data( 0 ),
				m_kappa( 0 )
			{
				set_data( data, data_range ) ;
			}

			MultivariateT( Matrix const& data, Vector const& weights, double const degrees_of_freedom ):
				m_pi( 3.141592653589793238462643383279502884 ),
				m_nu( degrees_of_freedom ),
				m_data( 0 ),
				m_kappa( 0 )
			{
				set_data( data, DataRange( 0, data.rows() ), weights ) ;
			}

			MultivariateT( Matrix const& data, DataSubset const& data_range, Vector const& weights, double const degrees_of_freedom ):
				m_pi( 3.141592653589793238462643383279502884 ),
				m_nu( degrees_of_freedom ),
				m_data( 0 ),
				m_kappa( 0 )
			{
				set_data( data, weights ) ;
			}

			void set_data( Matrix const& data ) {
				set_data( data, DataRange( 0, data.rows() ) ) ;
			}

			void set_data(
				Matrix const& data,
				DataSubset const& data_range
			) {
				set_data( data, data_range, Vector::Constant( data.rows(), 1 )) ;
			}

			void set_data(
				Matrix const& data,
				Vector const& weights
			) {
				set_data( data, DataSubset( 1, DataRange( 0, data.rows() ), weights ) ) ; 
			}

			void set_data( Matrix const& data, DataSubset const& data_range, Vector const& weights ) {
				m_data = &data ;
				m_ranges = data_range ;
				assert( weights.size() == m_data->rows() ) ;
				m_weights = weights ;
				m_p = m_data->cols() ;
				m_kappa = compute_constant_terms( m_nu, m_p ) ;
				m_parameters = Vector::Zero( m_p + ( m_p * ( m_p + 1 ) / 2 ) ) ;
			}

			// Evaluate at a packed set of parameters
			// These are:
			// The p values of the mean, followed by
			// The p (p+1)/2 entries of the lower diagonal of the sigma matrix (in column-major order.)
			void evaluate_at( Vector const& parameters ) {
				assert( parameters.size() == m_parameters.size() ) ;
				assert( m_data ) ;
				unpack_parameters( parameters, &m_mean, &m_sigma ) ;
				evaluate_at( m_mean, m_sigma ) ;
			}

			void evaluate_at( Vector const& mean, Matrix const& sigma ) {
				assert( mean.size() == m_data->cols() ) ;
				assert( sigma.rows() == sigma.cols() ) ;
				assert( sigma.rows() == m_data->cols() ) ;

				m_mean = mean ;
				m_sigma = sigma ;

				// Pack parameters into the parameter vector
				pack_parameters( mean, m_sigma, &m_parameters ) ;

				// We compute up-front some quantities that are useful
				// when computing the log-likelihood.

				m_ldlt.compute( m_sigma ) ;
				m_log_determinant = m_ldlt.vectorD().array().log().sum() ;
				m_mean_centred_data = m_data->rowwise() - m_mean.transpose() ;
				// Z is a vector of the terms ( x_i - mu )^t Sigma^-1 ( x_i - mu ).
				m_Z = (
					m_mean_centred_data.array()
					* ( m_ldlt.solve( m_mean_centred_data.transpose() ).transpose().array() )
				).rowwise().sum() ;

#if DEBUG_MULTIVARIATE_T
					std::cerr << "metro::likelihood::MultivariateT::evaluate_at():\n"
						<< " m_mean " << m_mean.transpose() << "\n"
						<< " m_sigma = \n" << m_sigma << "\n" 
						<< " m_log_determinant = " << m_log_determinant << "\n" ;
#endif				
			}

			double get_value_of_function() const {
				double result = 0 ;
				if( m_nu == std::numeric_limits< double >::infinity() ) {
					// Multivariate normal.
					for( std::size_t i = 0; i < m_ranges.number_of_subranges(); ++i ) {
						DataRange const& range = m_ranges[i] ;
						ConstSegment const segmentZ = m_Z.segment( range.begin(), range.size() ) ;
						ConstSegment const segmentWeights = m_weights.segment( range.begin(), range.size() ) ;
						double const weightSum = segmentWeights.sum() ;
						
						result +=
							( weightSum * m_kappa )
							- ( 0.5 * weightSum * m_log_determinant )
							- ( 0.5 * ( segmentZ.array() * segmentWeights.array() ).sum() ) ;
					}
				} else {
					for( std::size_t i = 0; i < m_ranges.number_of_subranges(); ++i ) {
						DataRange const& range = m_ranges[i] ;
						ConstSegment const segmentZ = m_Z.segment( range.begin(), range.size() ) ;
						ConstSegment const segmentWeights = m_weights.segment( range.begin(), range.size() ) ;
						double const weightSum = segmentWeights.sum() ;
						result += -(( m_nu + m_p )
							* (
								(( segmentZ + Vector::Constant( range.size(), m_nu ) ).array().log() )
									* segmentWeights.array()
							).sum() ) / 2.0 ;
						result += weightSum * ( m_kappa - 0.5 * m_log_determinant ) ;
#if DEBUG_MULTIVARIATE_T
					std::cerr << "metro::likelihood::MultivariateT::get_value_of_function():\n"
						<< " Adding " << range << "\n"
						<< " segmentZ = " << segmentZ.transpose() << "\n"
						<< " segmentWeights = " << segmentWeights.transpose() << "\n"
						<< " weightSum = " << weightSum << "\n"
						<< " m_kappa = " << m_kappa << "\n"
						<< " m_log_determinant = " << m_log_determinant << "\n"
						<< " term = " << (( segmentZ + Vector::Constant( range.size(), m_nu ) ).array().log() )
						<< " result = " << result << ".\n" ;
#endif				

					}
				}
				return result ;
			}

			Vector get_value_of_first_derivative() const {
				assert(0) ;
			}

			Matrix get_value_of_second_derivative() const {
				assert(0) ;
			}


			// Fit multivariate T by EM until the loglikelihood increases by less than the given amount.
			template< typename StoppingCondition >
			bool estimate_by_em(
				StoppingCondition& stopping_condition
			) {
				Matrix regularising_sigma = Matrix::Constant( m_p, m_p, 0 ) ;
				double regularising_weight = 0 ;
				return estimate_by_em(
					stopping_condition,
					regularising_sigma,
					regularising_weight
				) ;
			}

			// Fit multivariate T by EM until the stopping condition is satisfied.
			// Algorithm details are from Nadarajah & Kotz, "Estimation methods for the Multivariate t Distribution.", p.103
			// This function includes a regularising variance-covariance matrix, and weight
			// to prevent the fit from becoming degenerate.
			//
			// StoppingCondition must support operator() with result convertible to bool.
			// a value of true means stop, otherwise continue.
			// It must also be callable as stopping_condition.converged(), which this function uses
			// to return a value to the caller.
			template< typename StoppingCondition >
			bool estimate_by_em(
				StoppingCondition& stopping_condition,
				Matrix const& regularising_sigma,
				double regularising_weight
			) {
				assert( regularising_sigma.rows() == m_p ) ;
				assert( regularising_sigma.cols() == m_p ) ;
				assert( regularising_weight >= 0.0 ) ;

				// Start with unit weights, giving MVN estimate
				Vector iterationWeights = Vector::Constant( m_data->rows(), 1 ) ;
				Vector mean = compute_weighted_mean( iterationWeights, m_weights, m_ranges ) ;
				Matrix sigma = compute_weighted_regularised_sigma( iterationWeights, m_weights, mean, regularising_sigma, regularising_weight, m_ranges ) ;

				evaluate_at( mean, sigma ) ;
				double loglikelihood = get_value_of_function() ;
				
#if DEBUG_MULTIVARIATE_T
				std::cerr << "metro::likelihood::MultivariateT::estimate_by_em(): start: params = "
					<< get_parameters().transpose() << ", ll = " << get_value_of_function() << ".\n" ;
#endif				

				// If nu = ∞ we are at the MLE already, so bail out...
				if( m_nu == std::numeric_limits< double >::infinity() ) {
					// Multivariate normal.  Quit right now.
					return true ;
				}

				// ..otherwise let's EM it.
				std::size_t iteration = 0 ;
				while( !stopping_condition( loglikelihood ) ) {
					// compute weights
					// Vector of weights is given as
					// (nu+p) / nu + (x_i-mean)^t R^-1 ( x_i - mean ).
					// Our x_i - mean_i is stored in a single row of m_mean_centred_data.
					iterationWeights = (
						m_mean_centred_data.array() * m_ldlt.solve( m_mean_centred_data.transpose() ).transpose().array()
					).rowwise().sum() ;
					iterationWeights += Vector::Constant( m_data->rows(), m_nu ) ;
					iterationWeights.array() = iterationWeights.array().inverse() * ( m_nu + m_p ) ;

					// compute new parameter estimates
					// these are
					// mean = sum( w_i x_i ) / sum( w_i )
					mean = compute_weighted_mean( iterationWeights, m_weights, m_ranges ) ;
					sigma = compute_weighted_regularised_sigma( iterationWeights, m_weights, mean, regularising_sigma, regularising_weight, m_ranges ) ;
					evaluate_at( mean, sigma ) ;
					loglikelihood = get_value_of_function() ;

#if DEBUG_MULTIVARIATE_T
					std::cerr << "metro::likelihood::MultivariateT::estimate_by_em(): after iteration "
						<< iteration << ": params = " << get_parameters().transpose()
						<< ", ll = " << loglikelihood
						<< ", iterationWeights = " << iterationWeights.head( std::min( iterationWeights.size(), 10l )).transpose() << ".\n" ;
#endif				
					++iteration ;
				}
				return stopping_condition.converged() ;
			}

			Vector const& get_parameters() const { return m_parameters ; }
			double get_degrees_of_freedom() const { return m_nu ; }
			Vector const& get_mean() const { return m_mean ; }
			Matrix const& get_sigma() const { return m_sigma ; }
			Matrix const& get_data() const { return m_data ; }
			
			std::string get_spec() const { return "MultivariateT" ; }

		private:
			double const m_pi ;
			double const m_nu ;
			Matrix const* m_data ;
			DataSubset m_ranges ;
			Vector m_weights ;
			double m_p ;
			double m_kappa ;
			Vector m_parameters ;
			Matrix m_sigma ;
			Vector m_mean ;

			Eigen::LDLT< Matrix > m_ldlt ;
			double m_log_determinant ;
			Matrix m_mean_centred_data ;
			Vector m_Z ;
			
		private:
			
			double compute_constant_terms( double const nu, double const p ) const {
				return ( nu == std::numeric_limits< double >::infinity() )
					? (
						-0.5 * p * std::log( 2 * m_pi )
					) : (
						lgamma( ( nu + p ) / 2.0 )
						- ( p / 2 ) * std::log( m_pi * nu )
						- lgamma( nu / 2.0 )
						+ (( nu + p ) * std::log( nu ) / 2.0 )
					)
				;
			}

			void pack_parameters( Vector const& mean, Matrix const& sigma, Vector* result ) const {
				result->resize( mean.size() + ( sigma.rows() * ( sigma.rows() + 1 ) / 2 )) ;
				result->segment( 0, mean.size() ) = mean ;
				int parameter_i = mean.size() ;
				for( int col = 0; col < m_p; ++col ) {
					result->segment( parameter_i, m_p - col ) = sigma.col( col ).segment( col, m_p - col ) ;
					parameter_i += m_p - col ;
				}
			}

			void unpack_parameters( Vector const& parameters, Vector* mean, Matrix* sigma ) const {
				*mean = parameters.head( m_p ) ;
				sigma->resize( m_p, m_p ) ;
				int parameter_i = m_p ;
				for( int col = 0; col < m_data->cols(); ++col ) {
					sigma->col( col ).segment( col, m_p - col ) = parameters.segment( parameter_i, m_p - col ) ;
					parameter_i += m_p - col ;
				}
			}

			Vector compute_weighted_mean(
				Vector const& iterationWeights,
				Vector const& dataWeights,
				DataSubset const& data_range
			) const {
				Vector result = Vector::Zero( m_p ) ;
				double total_weight = 0 ;
				for( std::size_t i = 0; i < data_range.number_of_subranges(); ++i ) {
					DataRange const& range = data_range[i] ;
					ConstSegment segmentIW = iterationWeights.segment( range.begin(), range.size() ) ;
					ConstSegment segmentDW = dataWeights.segment( range.begin(), range.size() ) ;
					
					result += (
						(
							( segmentIW.array() * segmentDW.array() ).matrix().asDiagonal()
							* m_data->block( range.begin(), 0, range.size(), m_p )
						).colwise().sum()
					).transpose() ;

					total_weight += ( segmentDW.array() * segmentIW.array() ).sum() ;
				}
				result /= total_weight ;
				return result ;
			}
			
			Matrix compute_weighted_regularised_sigma(
				Vector const& iterationWeights,
				Vector const& dataWeights,
				Vector const& mean,
				Matrix const& regularising_sigma,
				double const regularising_weight,
				DataSubset const& data_range
			) const {
				Matrix mean_centred_data = m_data->rowwise() - mean.transpose() ;
				// sigma = 1/N sum ( w_i (x_i-mean)(x_i-mean)^t)
				// we store x_i as a row not a column, so transposes go the opposite way.
				Matrix result = Matrix::Zero( m_p, m_p ) ;
				double n = 0 ;
				for( std::size_t i = 0; i < data_range.number_of_subranges(); ++i ) {
					DataRange const& range = data_range[i] ;
					ConstSegment segmentIW = iterationWeights.segment( range.begin(), range.size() ) ;
					ConstSegment segmentDW = dataWeights.segment( range.begin(), range.size() ) ;

					result += (
						(
							mean_centred_data.block( range.begin(), 0, range.size(), m_p ).transpose()
							* segmentIW.asDiagonal()
							* segmentDW.asDiagonal()
							* mean_centred_data.block( range.begin(), 0, range.size(), m_p )
						)
					) ;
					n += segmentDW.sum() ;
				}
				result += regularising_weight * regularising_sigma ;
				result /= ( n + regularising_weight ) ;
				return result ;
			}
		} ;
	}
}

#endif
