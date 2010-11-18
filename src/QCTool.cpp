#include <string>
#include <map>
#include <set>
#include <fstream>
#include <memory>
#include <numeric>
#include <boost/bind.hpp>
#include <boost/function.hpp>
#include <iostream>
#include <iomanip>
#include "genfile/SNPDataSourceProcessor.hpp"

#include "Timer.hpp"
#include "GenRow.hpp"
#include "SampleRow.hpp"
#include "AlleleProportions.hpp"
#include "GToolException.hpp"
#include "RowCondition.hpp"
#include "GenRowStatistics.hpp"
#include "SampleRowStatistics.hpp"
#include "ObjectSource.hpp"
#include "SimpleFileObjectSource.hpp"
#include "SimpleFileObjectSink.hpp"
#include "OstreamTee.hpp"
#include "ProgramFlow.hpp"
#include "QCTool.hpp"
#include "QCToolContext.hpp"

QCTool::QCTool( QCToolContext & context ):
	m_context( context ),
	m_number_of_snps_processed( 0 )
{
}

void QCTool::begin_processing_snps(
	std::size_t number_of_samples,
	std::size_t number_of_snps
) {
	m_number_of_samples = number_of_samples ;
	m_number_of_snps = number_of_snps ;
	m_number_of_snps_processed = 0 ;
	m_timer.restart() ;
	m_context.logger() << "Processing SNPs...\n" ;
}

void QCTool::processed_snp(
	genfile::SNPIdentifyingData const& id_data,
	genfile::SingleSNPGenotypeProbabilities const& genotypes
) {
	try {
		unsafe_call_processed_snp( id_data, genotypes ) ;
	}
	catch( StatisticNotFoundException const& e ) {
		std::cerr << "!! ERROR: " << e << ".\n" ;
		std::cerr << "Note: required statistics must be added using -statistics.\n" ;
		throw HaltProgramWithReturnCode( -1 ) ;
	}
	catch( GToolException const& e) {
		std::cerr << "!! ERROR: " << e << ".\n" ;
		throw HaltProgramWithReturnCode( -1 ) ;
	}
	catch( GenAndSampleFileMismatchException const e ) {
		std::cerr << "!! ERROR: " << e.what() << ":\n"
			<< "There is a mismatch in sample counts between the GEN and sample files.\n"
			<< "  (" << e.actual_number_of_samples() << " in the GEN files versus " << e.expected_number_of_samples() << " in the sample files).\n" ; 
		throw HaltProgramWithReturnCode( -1 ) ;
	}
}

void QCTool::end_processing_snps() {
	assert( m_number_of_snps == m_number_of_snps_processed ) ;
	m_context.logger() << "\nProcessed " << m_number_of_snps << " SNPs in "
		<< std::fixed << std::setprecision(1) << m_timer.elapsed() << " seconds.\n" ;

	if( m_context.snp_filter().number_of_subconditions() > 0 ) {
		m_context.logger() << "(" << m_context.fltrd_in_snp_data_sink().number_of_snps_written() << " of " << m_number_of_snps << " SNPs passed the filter.)\n" ;
	}

	try {
		process_sample_rows() ;
		construct_plots() ;
	}
	catch( StatisticNotFoundException const& e ) {
		std::cerr << "!! ERROR: " << e << ".\n" ;
		std::cerr << "Note: required statistics must be added using -statistics.\n" ;
		throw HaltProgramWithReturnCode( -1 ) ;
	}
}

void QCTool::unsafe_call_processed_snp(
	genfile::SNPIdentifyingData const& id_data,
	genfile::SingleSNPGenotypeProbabilities const& genotypes
) {
	InternalStorageGenRow row( id_data, genotypes ) ;
	process_gen_row( row, ++m_number_of_snps_processed ) ;
	accumulate_per_column_amounts( row, m_per_column_amounts ) ;
	m_context.print_progress_if_necessary() ;
}		

void QCTool::process_gen_row( GenRow const& row, std::size_t row_number ) {
	m_context.snp_statistics().process( row ) ;
	if( m_context.snp_filter().check_if_satisfied( m_context.snp_statistics() )) {
		row.write_to_sink( m_context.fltrd_in_snp_data_sink() ) ;
		output_gen_row_stats( m_context.snp_statistics(), row_number ) ;
	}
	else {
		row.write_to_sink( m_context.fltrd_out_snp_data_sink() ) ;
		do_snp_filter_diagnostics( m_context.snp_statistics(), row_number ) ;
	}
}

void QCTool::output_gen_row_stats( GenotypeAssayStatistics const& row_statistics, std::size_t row_number ) {
	if( m_context.snp_statistics().size() > 0 ) {
		m_context.snp_stats_sink() << row_number ;
		for( std::size_t i = 0 ; i < row_statistics.size(); ++i ) {
			m_context.snp_stats_sink() << row_statistics.get_value< std::string >( row_statistics.get_statistic_name( i )) ;
		}
		m_context.snp_stats_sink() << statfile::end_row() ;
	}
}

void QCTool::do_snp_filter_diagnostics( GenRowStatistics const& row_statistics, std::size_t const row_index ) {
	std::ostream& log = m_context.logger()[ "log" ] ;
	log
		<< "Filtered out snp #" << row_index << " (" << row_statistics.row().SNPID() << " " << row_statistics.row().RSID() << " " << row_statistics.row().SNP_position() << ")"
		<< " because it does not satisfy " ;
	for( std::size_t i = 0, failed_condition_count = 0; i < m_context.snp_filter().number_of_subconditions() ; ++i ) {
		if( !m_context.snp_filter().subcondition(i).check_if_satisfied( row_statistics )) {
			++m_context.snp_filter_failure_counts()[ i ] ;
			if( failed_condition_count > 0 ) {
				log << " or " ;
			}
			log << "\"" << m_context.snp_filter().subcondition(i) << "\"" ;
			++failed_condition_count ;
		}
	}
	log << ".\n" ;
}

void QCTool::accumulate_per_column_amounts( GenRow& row, std::vector< GenotypeProportions >& per_column_amounts ) {
	// Keep totals for per-column stats.
	if( per_column_amounts.empty() ) {
		per_column_amounts.reserve( row.number_of_samples() ) ;
		std::copy( row.begin_genotype_proportions(), row.end_genotype_proportions(), std::back_inserter( per_column_amounts )) ;
	}
	else {
		assert( per_column_amounts.size() == row.number_of_samples() ) ;
		std::transform( per_column_amounts.begin(), per_column_amounts.end(),
		 				row.begin_genotype_proportions(),
		 				per_column_amounts.begin(),
						std::plus< GenotypeProportions >() ) ;
	}
}

void QCTool::process_sample_rows() {
	m_context.logger() << "Processing samples...\n" ;
	Timer timer ;

	apply_sample_filter() ;
	assert( m_context.sample_rows().size() == m_per_column_amounts.size() ) ;

	for( std::size_t i = 0 ; i < m_per_column_amounts.size(); ++i ) {
		SampleRow& sample_row = m_context.sample_rows()[i] ;
		m_context.sample_statistics().process( sample_row, m_per_column_amounts[i], m_number_of_snps ) ;
		output_sample_stats( i + 1, m_context.sample_statistics() ) ;
		
		m_context.sample_statistics().add_to_sample_row( sample_row, "missing" ) ;
		m_context.sample_statistics().add_to_sample_row( sample_row, "heterozygosity" ) ;
		m_context.fltrd_in_sample_sink() << sample_row ;
	}

	m_context.logger() << "Processed " << m_number_of_samples << " samples in " << std::fixed << std::setprecision(1) << timer.elapsed() << " seconds.\n" ;
	if( m_context.sample_filter().number_of_subconditions() > 0 ) {
		m_context.logger() << "(" << m_context.sample_rows().size() << " of " << m_number_of_samples << " samples passed the filter.)\n" ;
	}
	m_context.logger() << "\n" ;
}

void QCTool::output_sample_stats( std::size_t index, GenotypeAssayStatistics const& stats ) {
	m_context.sample_stats_sink() << index ;
	for( std::size_t i = 0; i < m_context.sample_statistics().size(); ++i ) {
		m_context.sample_stats_sink() << stats.get_value< std::string >( stats.get_statistic_name( i )) ;
	}
	m_context.sample_stats_sink() << statfile::end_row() ;
}

void QCTool::apply_sample_filter() {
	for( std::size_t pre_filter_i = 0, post_filter_i = 0; post_filter_i < m_context.sample_rows().size(); ++pre_filter_i ) {
		if( sample_row_is_filtered_out( pre_filter_i ) ) {
			m_context.fltrd_out_sample_sink().write( m_context.sample_rows()[ post_filter_i] ) ;
			do_sample_filter_diagnostics( m_context.sample_rows()[ post_filter_i ], pre_filter_i ) ;
			m_context.sample_rows().erase( m_context.sample_rows().begin() + post_filter_i ) ;
		}
		else {
			++post_filter_i ;
		}
	}
}

bool QCTool::sample_row_is_filtered_out( std::size_t const sample_row_index ) {
	return std::binary_search( m_context.indices_of_filtered_out_samples().begin(), m_context.indices_of_filtered_out_samples().end(), sample_row_index ) ;
}

void QCTool::do_sample_filter_diagnostics( SampleRow const& sample_row, std::size_t const sample_row_index ) {
	std::ostream& log = m_context.logger()["log"] ;
	log
		<< "Filtered out sample row " << sample_row_index << " (" << sample_row.ID1() << " " << sample_row.ID2() << ")"
		<< " because it does not satisfy " ;
	for( std::size_t i = 0, failed_condition_count = 0; i < m_context.sample_filter().number_of_subconditions() ; ++i ) {
		if( !m_context.sample_filter().subcondition(i).check_if_satisfied( sample_row )) {
			++m_context.sample_filter_failure_counts()[ i ] ;
			if( failed_condition_count > 0 ) {
				log << " or " ;
			}
			log << "\"" << m_context.sample_filter().subcondition(i) << "\"" ;
			++failed_condition_count ;
		}
	}
	log << ".\n" ;
}

void QCTool::construct_plots() {
	// not implemented.
}
