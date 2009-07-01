
#ifndef __GTOOL_ROWCONDITION_HPP__
#define __GTOOL_ROWCONDITION_HPP__

#include <set>
#include "GenRow.hpp"
#include "GenotypeAssayStatistics.hpp"
#include "Condition.hpp"

typedef Condition< GenotypeAssayStatistics > RowCondition ;
typedef CompoundCondition< GenotypeAssayStatistics > CompoundRowCondition ;
typedef AndCondition< GenotypeAssayStatistics > AndRowCondition ;
typedef OrCondition< GenotypeAssayStatistics > OrRowCondition ;
typedef InvertCondition< GenotypeAssayStatistics > NotRowCondition ;
typedef InvertCondition< GenotypeAssayStatistics > InvertRowCondition ;

struct TrivialRowCondition: public RowCondition
{
	bool check_if_satisfied( GenotypeAssayStatistics const& row_genotype_statistics ) const {
		return true ;
	}
	
	void format_to_stream( std::ostream& oStream ) const ;
} ;

struct GenotypeAssayStatisticInInclusiveRange: public RowCondition
{
	GenotypeAssayStatisticInInclusiveRange( std::string const& statistic_name, double lower_bound, double upper_bound, double epsilon = 0.0 ) ;

	bool check_if_satisfied( GenotypeAssayStatistics const& row_genotype_statistics ) const ;
	
	void format_to_stream( std::ostream& oStream ) const ;

	private:
	
		std::string const m_statistic_name ;
		double m_lower_bound, m_upper_bound, m_epsilon ;
} ;

struct GenotypeAssayStatisticInExclusiveRange: public RowCondition
{
	GenotypeAssayStatisticInExclusiveRange( std::string const& statistic_name, double lower_bound, double upper_bound, double epsilon = 0.0 ) ;

	bool check_if_satisfied( GenotypeAssayStatistics const& row_genotype_statistics ) const ;
	
	void format_to_stream( std::ostream& oStream ) const ;
	
	private:
	
		std::string const m_statistic_name ;
		double m_lower_bound, m_upper_bound, m_epsilon ;
} ;

struct GenotypeAssayStatisticGreaterThan: public RowCondition
{
	GenotypeAssayStatisticGreaterThan( std::string const& statistic_name, double lower_bound, double epsilon = 0.0 ) ;

	bool check_if_satisfied( GenotypeAssayStatistics const& row_genotype_statistics ) const ;
	
	void format_to_stream( std::ostream& oStream ) const ;
	
	private:
	
		std::string const m_statistic_name ;
		double m_lower_bound, m_epsilon ;
} ;

struct GenotypeAssayStatisticLessThan: public RowCondition
{
	GenotypeAssayStatisticLessThan( std::string const& statistic_name, double upper_bound, double epsilon = 0.0 ) ;

	bool check_if_satisfied( GenotypeAssayStatistics const& row_genotype_statistics ) const ;
	
	void format_to_stream( std::ostream& oStream ) const ;
	
	private:
	
		std::string const m_statistic_name ;
		double m_upper_bound, m_epsilon ;
} ;


// Factory function for conditions.
std::auto_ptr< RowCondition > condition_factory( std::string condition_spec ) ;


#endif

