
#include <cassert>
#include "RowCondition.hpp"

StatisticInInclusiveRange::StatisticInInclusiveRange( std::string const& statistic_name, double lower_bound, double upper_bound, double epsilon )
: m_statistic_name( statistic_name ),
	m_lower_bound( lower_bound ),
	m_upper_bound( upper_bound ),
	m_epsilon( epsilon )
{} ;

bool StatisticInInclusiveRange::check_if_satisfied( string_to_value_map const& row_genotype_statistics ) const {
	double statistic_value = row_genotype_statistics.get_value< double >( m_statistic_name ) ;
	return ( statistic_value >= ( m_lower_bound - m_epsilon ))
		&& ( statistic_value <= ( m_upper_bound - m_epsilon )) ;
}

void StatisticInInclusiveRange::format_to_stream( std::ostream& oStream ) const {
	oStream
		<< m_statistic_name
		<< " in ["
	 	<< m_lower_bound
		<< ","
		<< m_upper_bound
		<< "]" ;
}

StatisticInExclusiveRange::StatisticInExclusiveRange( std::string const& statistic_name, double lower_bound, double upper_bound, double epsilon )
: m_statistic_name( statistic_name ),
	m_lower_bound( lower_bound ),
	m_upper_bound( upper_bound ),
	m_epsilon( epsilon )
{} ;

bool StatisticInExclusiveRange::check_if_satisfied( string_to_value_map const& row_genotype_statistics ) const {
	double statistic_value = row_genotype_statistics.get_value< double >( m_statistic_name ) ;
	return ( statistic_value > ( m_lower_bound - m_epsilon ))
		&& ( statistic_value < ( m_upper_bound - m_epsilon )) ;
}

void StatisticInExclusiveRange::format_to_stream( std::ostream& oStream ) const {
	oStream
		<< m_statistic_name
		<< " in ("
	 	<< m_lower_bound
		<< ","
		<< m_upper_bound
		<< ")" ;
}

StatisticGreaterThan::StatisticGreaterThan( std::string const& statistic_name, double lower_bound, double epsilon )
: m_statistic_name( statistic_name ),
	m_lower_bound( lower_bound ),
	m_epsilon( epsilon )
{} ;

bool StatisticGreaterThan::check_if_satisfied( string_to_value_map const& row_genotype_statistics ) const {
	double statistic_value = row_genotype_statistics.get_value< double >( m_statistic_name ) ;
	return ( statistic_value > ( m_lower_bound - m_epsilon )) ;
}

void StatisticGreaterThan::format_to_stream( std::ostream& oStream ) const {
	oStream
		<< m_statistic_name
		<< " > "
		<< m_lower_bound ;
}

StatisticLessThan::StatisticLessThan( std::string const& statistic_name, double upper_bound, double epsilon )
: m_statistic_name( statistic_name ),
	m_upper_bound( upper_bound ),
	m_epsilon( epsilon )
{} ;

bool StatisticLessThan::check_if_satisfied( string_to_value_map const& row_genotype_statistics ) const {
	double statistic_value = row_genotype_statistics.get_value< double >( m_statistic_name ) ;
	return ( statistic_value < ( m_upper_bound - m_epsilon )) ;
}

void StatisticLessThan::format_to_stream( std::ostream& oStream ) const {
	oStream
		<< m_statistic_name
		<< " < "
		<< m_upper_bound ;
}

