
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#include <string>
#include <vector>
#include <utility>
#include <cmath>
#include <stdlib.h>
#include <fstream>
#include <ctime>
#include <iostream>
#include <deque>
#include <algorithm>

#include <boost/bind.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_int.hpp>
#include <boost/random/uniform_real_distribution.hpp>
#include <boost/iterator/counting_iterator.hpp>
#include <boost/filesystem.hpp>
#include <boost/format.hpp>

#include "appcontext/ProgramFlow.hpp"
#include "appcontext/CmdLineOptionProcessor.hpp"
#include "appcontext/ApplicationContext.hpp"
#include "appcontext/get_current_time_as_string.hpp"

#include "genfile/SNPDataSource.hpp"
#include "genfile/SNPDataSourceChain.hpp"
#include "genfile/SNPFilteringSNPDataSource.hpp"
#include "genfile/get_list_of_snps_in_source.hpp"
#include "genfile/utility.hpp"
#include "genfile/CommonSNPFilter.hpp"
#include "genfile/VariantEntry.hpp"
#include "genfile/string_utils.hpp"
#include "genfile/vcf/StrictMetadataParser.hpp"

#include "statfile/BuiltInTypeStatSourceChain.hpp"
#include "statfile/DelimitedStatSource.hpp"

#include "db/SQLite3Connection.hpp"
#include "genfile/FromFilesGeneticMap.hpp"
#include "FileUtil.hpp"
#include "qcdb/Storage.hpp"
#include "qcdb/FlatTableDBOutputter.hpp"
#include "qcdb/FlatFileOutputter.hpp"

// #define DEBUG_INTHINNERATOR 1

namespace globals {
	std::string const program_name = "inthinnerator" ;
}

namespace {
	std::size_t get_random_seed() {
		std::size_t seed ;
		if( boost::filesystem::exists( "/dev/random" )) {
			std::ifstream ifs( "/dev/random" ) ;
			if( !ifs.is_open() ) {
				throw genfile::ResourceNotOpenedError( "/dev/random" ) ;
			}
			char buf[ sizeof( std::size_t ) ] ;
			ifs.read( buf, sizeof( std::size_t )) ;
			ifs.close() ;
			seed = *(reinterpret_cast< std::size_t* >( buf )) ;
		}
		else {
			seed = std::time(0) ;
		}
		return seed ;
	}
	
	struct CompareByOrder {
		typedef bool result_type ;
		
		CompareByOrder( std::deque< std::size_t > const& list ) {
			m_positions_in_sorted_list.resize( list.size(), std::numeric_limits< std::size_t >::max() ) ;
			for( std::size_t i = 0; i < list.size(); ++i ) {
				assert( list[i] < m_positions_in_sorted_list.size() ) ;
				assert( m_positions_in_sorted_list[ list[i] ] == std::numeric_limits< std::size_t >::max() ) ;
				m_positions_in_sorted_list[ list[i] ] = i ;
			}
		}

		bool operator()( std::size_t i, std::size_t j ) const {
			assert( i < m_positions_in_sorted_list.size() ) ;
			assert( j < m_positions_in_sorted_list.size() ) ;
			return m_positions_in_sorted_list[i] < m_positions_in_sorted_list[j] ;
		}
		
		private:
			std::vector< std::size_t > m_positions_in_sorted_list ;
	} ;

	struct TaggedSnp {
		TaggedSnp()
		{}

		TaggedSnp( genfile::SNPIdentifyingData2 const& snp, boost::optional< std::string > const tag ):
			m_snp( snp ),
			m_tag( tag )
		{}

		TaggedSnp( TaggedSnp const& other ):
			m_snp( other.m_snp ),
			m_tag( other.m_tag )
		{}

		TaggedSnp& operator=( TaggedSnp const& other ) {
			m_snp = other.m_snp ;
			m_tag = other.m_tag ;
			return *this ;
		}
		
		genfile::SNPIdentifyingData2 const& snp() const { return m_snp ; }
		boost::optional< std::string > const tag() const { return m_tag ; }

	private:
		genfile::SNPIdentifyingData2 m_snp ;
		boost::optional< std::string > m_tag ;
	} ;

	bool operator<( TaggedSnp const& left, TaggedSnp const& right ) {
		return left.snp() < right.snp() || ( left.snp() == right.snp() && left.tag() < right.tag() ) ;
	}
	
	namespace impl {
		// Compare SNP indices using their sorted order.
		genfile::SNPIdentifyingData2 const& get_snp( std::vector< TaggedSnp > const& list, std::size_t i ) {
			return list[i].snp() ;
		}
	
		void add_snp( std::vector< TaggedSnp >* result, genfile::SNPIdentifyingData2 const& snp, boost::optional< std::string > const& tag ) {
			assert( result ) ;
			result->push_back( TaggedSnp( snp, tag )) ;
		}
	}
}

struct InthinneratorOptionProcessor: public appcontext::CmdLineOptionProcessor
{
	// Methods needed for CmdLineOptionProcessor::process()
	std::string get_program_name() const { return globals::program_name ; }
	
	void declare_options( OptionProcessor& options ) {
		options.set_help_option( "-help" ) ;

		options.declare_group( "Input file options" ) ;
		options[ "-map" ]
			.set_description( "Set the path of the genetic map panel to use." )
			.set_takes_single_value()
			.set_is_required()
		;
		options[ "-g" ]
			.set_description( "Specify a file containing the SNPs to operate on." )
			.set_is_required()
			.set_takes_values( 1 )
			.set_maximum_multiplicity( 100 ) ;
		options[ "-metadata" ]
			.set_description(
				"Specify the name of a file containing VCF metadata to be used to parse "
				"a VCF file.  Keys in this file must either not be present in the VCF file, or must have "
				"identical values."
			)
			.set_takes_single_value() ;
		options[ "-genes" ]
			.set_description(
				"Specify the name of a file containing genes (in UCSC table format).  If this is supplied, inthinnerator "
				"will annotate each output row with the nearest gene and the nearest gene in the region."
			)
			.set_takes_single_value()
		;
		options[ "-take-longest-transcript" ]
			.set_description(
				"When using -genes, tell inthinnerator to ignore all but the longest transcript for each gene in the file."
			)
			.set_default_value( false ) ;
		;
		options.option_implies_option( "-take-longest-transcript", "-genes" ) ;
			
		options.declare_group( "Inclusion / exclusion options" ) ;
		options[ "-excl-rsids" ]
			.set_description( "Specify a file containing a whitespace-separated list of SNP rsids."
				" SNPs with ids in this file will be excluded from the analysis." )
			.set_takes_values( 1 )
			.set_maximum_multiplicity( 100 ) ;
		options[ "-incl-rsids" ]
			.set_description( "Specify a file containing a whitespace-separated list of SNP rsids."
				" SNPs with ids not in this file will be excluded from the analysis." )
			.set_takes_single_value() ;
		options[ "-excl-snpids" ]
			.set_description( "Specify a file containing a whitespace-separated list of SNP SNPIDs."
			" SNPs with ids in this file will be excluded from the analysis." )
			.set_takes_values( 1 )
			.set_maximum_multiplicity( 100 ) ;
		options[ "-incl-snpids" ]
			.set_description( "Specify a file containing a whitespace-separated list of SNP SNPIDs."
			" SNPs with ids not in this file will be excluded from the analysis." )
			.set_takes_single_value() ;
		options[ "-incl-range" ]
			.set_description( "Specify a range of SNPs (or comma-separated list of ranges of SNPs) to operate on. "
				"Each range should be in the format CC:xxxx-yyyy where CC is the chromosome and xxxx and yyyy are the "
				"start and end coordinates, or just xxxx-yyyy which matches that range from all chromosomes. "
				"You can also omit either of xxxx or yyyy to get all SNPs from the start or to the end of a chromosome." )
			.set_takes_values_until_next_option()
			.set_maximum_multiplicity( 100 ) ;
		options[ "-excl-range" ]
			.set_description( "Specify a range of SNPs (or comma-separated list of ranges of SNPs) to exclude from operation. "
				"Each range should be in the format CC:xxxx-yyyy where CC is the chromosome and xxxx and yyyy are the "
				"start and end coordinates, or just xxxx-yyyy which matches that range from all chromosomes. "
				"You can also omit either of xxxx or yyyy to get all SNPs from the start or to the end of a chromosome." )
			.set_takes_values_until_next_option()
			.set_maximum_multiplicity( 100 ) ;

		options.declare_group( "SNP thinning options" ) ;
		options[ "-min-distance" ]
			.set_description( "Specify a minimum acceptable distance between SNPs."
				" This must be in one of the forms \"<X>cM\", \"<X>M\", \"<X>bp\", \"<X>kb\", or \"<X>Mb\"."
				" The thinned list of SNPs will contain no SNPs within this distance of each other."
				" For recombination-distance offsets, a physical margin can also be specified as in \"<X>cM+<Y>kb\"." )
			.set_takes_single_value()
			.set_default_value( "0.01cM" ) ;
		options[ "-strategy" ]
			.set_description( "Specify the SNP thinning strategy if not using -rank."
				" This can be \"random\", \"first\", or \"random_by_position\"." )
			.set_takes_single_value()
			.set_default_value( "random" ) ;
		options[ "-bin-size" ]
			.set_description( "Specify the size of bins when computing occupied genomic intervals." )
			.set_takes_single_value()
			.set_default_value( "1kb" ) ;
		options[ "-rank" ]
			.set_description( "Specify name of a file containing numerical ranks of SNPs. "
				"SNPs will be picked in nonincreasing order of rank. "
				"This file must have first five columns SNPID, rsid, chromosome, position, allele1, allele2."
			)
			.set_takes_single_value() ;
		options[ "-rank-column" ]
			.set_description( "Specify the name of the column in the file for -rank containing the ranks." )
			.set_takes_single_value() ;
		options[ "-missing-code" ]
			.set_description( "Specify a comma-separated list of strings to be treated as missing ranks. "
				"Missing ranks are always picked last." )
			.set_takes_single_value()
			.set_default_value( "NA" ) ;

		options[ "-max-picks" ]
			.set_description( "Specify a number of SNPs to pick in each thinning."
				" By default we choose as many SNPs as it's possible to choose subject to the minimum distance constraint." )
			.set_takes_single_value()
		;

		options.option_excludes_option( "-rank", "-strategy" ) ;
		options.option_excludes_option( "-strategy", "-rank" ) ;
		options.option_implies_option( "-rank", "-rank-column" ) ;
		options.option_implies_option( "-missing-code", "-rank" ) ;

		options[ "-match-tag" ]
			.set_description( "Specify a file (in the same format as the files supplied to -g) containing SNPs to match by tag."
				"When picking SNPs, the ith SNP picked will be from the set having the tag of the ith SNP in this file."
				"This option also implies that no more than the given number of SNPs will be picked."
			)
			.set_takes_single_value()
			.set_maximum_multiplicity(1)
		;
		options[ "-tag-column" ]
			.set_description( "Specify the name of a column in the file supplied to -g containing a tag for each variant." )
			.set_takes_values( 1 )
			.set_maximum_multiplicity( 1 ) ;
		;
		options.option_implies_option( "-match-tag", "-tag-column" ) ;
		options.option_implies_option( "-tag-column", "-match-tag" ) ;
		options.option_excludes_option( "-match-tag", "-max-picks" ) ;

		options.declare_group( "Repetition options" ) ;
		options[ "-N" ]
			.set_description( "Specify a number of thinnings to perform.  This must be at least 1."
				" If this is larger than one, output files will be numbered in the form <filename>.####,"
				" where <filename> is the filename passed to -o and #### is a number indexing the repetition." )
			.set_takes_single_value()
			.set_default_value( 1 ) ;
		options[ "-start-N" ]
			.set_description( "Specify the first index to be used when running mutliple thinnings.  This is"
							" useful when running jobs in parallel but giving output as though run in a single command." )
			.set_takes_single_value()
			.set_default_value(0) ;

		options.declare_group( "Output file options" ) ;
		options[ "-o" ]
			.set_description( "Specify the output filename stub." )
			.set_takes_single_value() ;
		options[ "-odb" ]
			.set_description( "Specify the name of a database file to output." )
			.set_takes_single_value() ;
		options[ "-output-cols" ]
			.set_description( "Specify a comma-separated list of columns that should appear in the output files."
				" Possible columns are \"SNPID\", \"rsid\", \"chromosome\", \"position\", \"allele1\", \"allele2\", and \"cM_from_start_of_chromosome\"."
				" The special value \"all\" indicates that all available columns will be output." )
			.set_takes_single_value()
			.set_default_value( "all" ) ;
		options[ "-table-name" ]
			.set_description( "Specify a name for the table to use when using -odb." )
			.set_takes_single_value() ;
		options[ "-headers" ]
			.set_description( "Specify this to force output of column headersomit column headers in the output files." ) ;
		options[ "-no-headers" ]
			.set_description( "Specify this to suppress output of column headers in the output files." ) ;

		options.option_excludes_option( "-o", "-odb" ) ;
		options.option_excludes_option( "-odb", "-o" ) ;
		options.option_excludes_option( "-headers", "-no-headers" ) ;
		options.option_implies_option( "-headers", "-o" ) ;
		options.option_implies_option( "-no-headers", "-o" ) ;
		
		options[ "-suppress-excluded" ]
			.set_description( "Specify that inthinnerator should not produce an exclusion lists." ) ;
		options[ "-suppress-included" ]
			.set_description( "Specify that inthinnerator should not produce an inclusion lists." ) ;
		
		options.declare_group( "Miscellaneous options" ) ;
		options[ "-log" ]
			.set_description( "Set the path of the log file to output." )
			.set_takes_single_value()
			.set_default_value( globals::program_name + ".log" ) ;
		options[ "-analysis-name" ]
			.set_description( "Specify a name to label results from this analysis with.  (This applies to modules which store their results in a qcdb file.)" )
			.set_takes_single_value()
			.set_default_value( globals::program_name + " analysis" ) ;
		options[ "-analysis-chunk" ]
			.set_description( "Specify a name denoting the current genomic region or chunk on which this is run.  This is intended for use in parallel environments." )
			.set_takes_single_value()
			.set_default_value( genfile::MissingValue() ) ;
	}
} ;

// SNPPicker encapsulates the operation of picking a SNP from a given list.
// Different strategies are possible and are implemented by base classes.
// Picking takes place in the context of a particular list of SNPs, passed in by the set_snps() method.
// The pick() function picks a single SNP (by its index) from the list passed in.
// Typically picking involves computation of useful quantities, such as recombination position.
// We allow SNPPicker to report these via the get_attributes() method.
class SNPPicker
{
public:
	typedef std::auto_ptr< SNPPicker > UniquePtr ;
	virtual ~SNPPicker() {}
	
	typedef boost::function< bool ( std::size_t, std::size_t ) > SNPComparator ;
	typedef boost::function< genfile::SNPIdentifyingData2 const& ( std::size_t ) > SNPGetter ;
	virtual void set_snps( std::size_t, SNPGetter, SNPComparator ) = 0 ;

	virtual std::size_t pick(
		std::deque< std::size_t > const& among_these
	) const = 0 ;
	// Tell the picker the order in which SNPs are sorted.
	// This permits pickers to do binary lookup in the list of SNPs.
	// TODO: we should rewrite this to take a boost::function which is the sort comparator.
	virtual std::string display() const = 0 ;
	virtual std::set< std::string > get_attribute_names() const = 0 ;
	virtual std::map< std::string, genfile::VariantEntry > get_attributes( std::size_t chosen_snp ) const = 0 ;
} ;

// Pick first available SNP
class FirstAvailableSNPPicker: public SNPPicker
{
public:
	FirstAvailableSNPPicker() {} ;

	void set_snps( std::size_t number_of_snps, SNPGetter getter, SNPComparator ) {
		// This picker just picks the first index available each time.
		// It does not care about what SNP that correspondeth to
	}

	std::size_t pick(
		std::deque< std::size_t > const& among_these
	) const {
		assert( among_these.size() > 0 ) ;
		return *(among_these.begin()) ;
	}
	
	std::string display() const {
		return "FirstAvailableSNPPicker" ;
	}
	
	std::set< std::string > get_attribute_names() const {
		return std::set< std::string >() ;
	}
	std::map< std::string, genfile::VariantEntry > get_attributes( std::size_t chosen_snp ) const {
		return std::map< std::string, genfile::VariantEntry >() ;
	} ;
	
	
} ;

// Pick a random SNP
class RandomSNPPicker: public SNPPicker
{
private:
	typedef boost::mt19937 RNG ;
	typedef boost::uniform_int< std::size_t > Distribution ;
public:
	RandomSNPPicker():
		m_rng( new RNG( get_random_seed() ) )
	{
	}

	RandomSNPPicker( std::size_t seed ):
		m_rng( new RNG( seed ) )
	{}

	void set_snps( std::size_t, SNPGetter, SNPComparator ) {
		// This picker just picks a random index each time.
		// It does not care about what SNP that correspondeth to
	}

	std::size_t pick(
		std::deque< std::size_t > const& among_these
	) const {
		assert( among_these.size() > 0 ) ;
		Distribution distribution( 0, among_these.size() - 1 ) ;
		std::size_t choice = distribution( *m_rng ) ;
		assert( choice < among_these.size() ) ;
		std::deque< std::size_t >::const_iterator i = among_these.begin() ;
		std::advance( i, choice ) ;
		return *i ;
	}

	std::string display() const {
		return "RandomSNPPicker" ;
	}
	
	std::set< std::string > get_attribute_names() const {
		return std::set< std::string >() ;
	}
	std::map< std::string, genfile::VariantEntry > get_attributes( std::size_t chosen_snp ) const {
		return std::map< std::string, genfile::VariantEntry >() ;
	} ;

private:
	std::auto_ptr< RNG > m_rng ;
} ;

// Pick a random SNP
class RandomPositionSNPPicker: public SNPPicker
{
private:
	typedef boost::mt19937 RNG ;
	typedef boost::random::uniform_real_distribution< double > Distribution ;
	typedef std::multimap< genfile::GenomePosition, std::size_t > SnpByPositionMap ;
	typedef std::map< std::pair< double, double >, genfile::GenomePositionRange > RangeMap ;
	
public:
	RandomPositionSNPPicker(
		std::vector< genfile::GenomePositionRange > ranges,
		std::size_t maximum_distance
	):
		m_rng( new RNG( get_random_seed() ) ),
		m_uniform_01( 0.0, 1.0 ),
		m_range_map( compute_range_map( ranges )),
		m_maximum_distance( maximum_distance )
	{
	}
	
	void set_snps( std::size_t number_of_snps, SNPGetter getter, SNPComparator comparator ) {
		assert( number_of_snps > 0 ) ;
		m_snps_by_position.clear() ;
		for( std::size_t i = 0; i < number_of_snps; ++i ) {
			m_snps_by_position.insert( std::make_pair( getter(i).get_position(), i )) ;
		}
		m_snp_comparator = comparator ;
	}
	
	std::size_t pick(
		std::deque< std::size_t > const& among_these
	) const {
		bool picked = false ;
		std::size_t result = 0 ;
		while( !picked ) {
			double in_01 = m_uniform_01( *m_rng ) ;
			genfile::GenomePosition pos = map_uniform_01_to_position_in_ranges( in_01 ) ;

#if DEBUG_INTHINNERATOR
			std::cerr << "RandomPositionSNPPicker::pick(): testing " << pos << "..." ;
			std::cerr << "RandomPositionSNPPicker::pick(): among " << among_these.size() << " SNPs...\n" ;
#endif
			// find the two SNPs flanking this position.
			SnpByPositionMap::const_iterator previous = m_snps_by_position.lower_bound( pos ) ;
			if( previous != m_snps_by_position.end() && previous->first == pos ) {
				// Have landed right on a SNP!  pick it.
				result = previous->second ;
				picked = std::binary_search( among_these.begin(), among_these.end(), result ) ;
			} else {
				SnpByPositionMap::const_iterator next = previous ;
				--previous ;
				genfile::Position distanceToPrevious = std::numeric_limits< genfile::Position >::max() ;
				genfile::Position distanceToNext = std::numeric_limits< genfile::Position >::max() ;
				if( previous->first.chromosome() == pos.chromosome() ) {
					distanceToPrevious = pos.position() - previous->first.position() ;
				}
				if( next != m_snps_by_position.end() && next->first.chromosome() == pos.chromosome() ) {
					distanceToNext = next->first.position() - pos.position() ;
				}

#if DEBUG_INTHINNERATOR
				std::cerr << "RandomPositionSNPPicker::pick(): previous SNP is at " << previous->first << ", next is at " << next->first << ".\n" ;
#endif

				if( distanceToPrevious <= m_maximum_distance || distanceToNext <= m_maximum_distance ) {
					if( distanceToNext < distanceToPrevious ) {
						result = next->second ;
					} else {
						result = previous->second ;
					}

					picked = std::binary_search( among_these.begin(), among_these.end(), result ) ;
				} else {
					// no pickable SNP within distance of the chosen position.
					// Go back and try again
					picked = false ;
				}
#if DEBUG_INTHINNERATOR
				std::cerr << (picked ? "picked!\n" : "not picked.\n") ;
#endif
			}
		}
#if DEBUG_INTHINNERATOR
		std::cerr << "Result is " << result << ".\n" ;
#endif
		return result ;
	}

	std::string display() const {
		std::ostringstream ostr ;
		ostr << "RandomPositionSNPPicker with map:\n" ;
		RangeMap::const_iterator i = m_range_map.begin(), end_i = m_range_map.end() ;
		for( ; i != end_i; ++i ) {
			ostr << "   " << i->first.first << "-" << i->first.second << " --> " << i->second << "\n" ;
		}
		return ostr.str() ;
	}
	
	std::set< std::string > get_attribute_names() const {
		return std::set< std::string >() ;
	}
	std::map< std::string, genfile::VariantEntry > get_attributes( std::size_t chosen_snp ) const {
		return std::map< std::string, genfile::VariantEntry >() ;
	} ;

private:
	std::auto_ptr< RNG > m_rng ;
	Distribution m_uniform_01 ;
	RangeMap const m_range_map ;
	SnpByPositionMap m_snps_by_position ;
	SNPComparator m_snp_comparator ;
	std::size_t const m_maximum_distance ;
	
private:
	RangeMap compute_range_map( std::vector< genfile::GenomePositionRange > const& ranges ) const {
		double totalLength = 0 ;
		for( std::size_t i = 0; i < ranges.size(); ++i ) {
			totalLength += ranges[i].end().position() - ranges[i].start().position() ;
		}
#if DEBUG_INTHINNERATOR
		std::cerr << "RandomPositionSNPPicker::compute_range_map(): Total length is " << totalLength << ".\n" ;
#endif
		
		RangeMap result ;
		double start = 0 ;
		for( std::size_t i = 0; i < ranges.size(); ++i ) {
			double end = start + double( ranges[i].end().position() - ranges[i].start().position() ) / totalLength ;
			result.insert( std::make_pair( std::make_pair( start, end ), ranges[i] ) ) ;
#if DEBUG_INTHINNERATOR
		std::cerr << "RandomPositionSNPPicker::compute_range_map(): inserted element " << boost::format( "%.3f-%.3f" ) % start % end << " --> " << ranges[i] << ".\n" ;
#endif
		start = end ;
		}
		return result ;
	}

	genfile::GenomePosition map_uniform_01_to_position_in_ranges( double in_01 ) const {
		using boost::format ;
		assert( in_01 >= 0.0 && in_01 < 1 ) ;
		std::pair< double, double > u( in_01, 1.0 ) ;
#if DEBUG_INTHINNERATOR
		std::cerr << "RandomPositionSNPPicker::map_uniform_01_to_position_in_ranges(): looking for " << format( "%.3f" ) % u.first << "..." ;
#endif
		RangeMap::const_iterator where = m_range_map.lower_bound(u) ;
		if( u.first != where->first.first ) {
			assert( where != m_range_map.begin() ) ;
			--where ;
		}
		double x = ( in_01 - where->first.first ) / ( where->first.second - where->first.first ) ; // interpoland
		genfile::Position p = where->second.start().position() + ( where->second.end().position() - where->second.start().position() ) * x ;
#if DEBUG_INTHINNERATOR
		std::cerr << "RandomPositionSNPPicker::map_uniform_01_to_position_in_ranges(): chose range " << boost::format( "%.3f-%.3f" ) % where->first.first % where->first.second << ", " << where->second << " " << format( "x = %.4f" ) % x << ", position = " << p << "\n" ;
#endif
		return genfile::GenomePosition( where->second.start().chromosome(), p ) ;
	}
} ;

class HighestValueSNPPicker: public SNPPicker
{
private:
	struct DoubleComparator {
		bool operator()( double const a, double const b ) const {
			// Compare doubles putting NaN's last of all.
			// From http://stackoverflow.com/questions/4816156/are-ieee-floats-valid-key-types-for-stdmap-and-stdset
		    if ((a == a) && (b == b)) {
		        return a < b ;
		    }
		    if ((a != a) && (b != b)) return false ;
		    // We have one NaN and one non-NaN.
		    // Let's say NaN is less than everything
			return (a != a) ;
		}
	} ;

public:
	typedef std::auto_ptr< SNPPicker > UniquePtr ;
	
	// TODO: use a boost::bimap here.
	typedef std::map< genfile::SNPIdentifyingData2, double > SnpToValueMap ;
	typedef std::vector< double > IndexToValueMap ;
	typedef std::multimap< double, std::size_t, DoubleComparator > ValueToIndexMap ;

	HighestValueSNPPicker(
		SnpToValueMap const& values
	):
		m_values_per_snp( values ),
		m_value_to_index_map( DoubleComparator() )
	{
	}

	void set_snps( std::size_t number_of_snps, SNPGetter getter, SNPComparator comparator ) {
		m_value_to_index_map.clear() ;
		m_index_to_value_map.resize( number_of_snps ) ;
		for( std::size_t i = 0; i < number_of_snps; ++i ) {
			SnpToValueMap::const_iterator where = m_values_per_snp.find( getter(i) ) ;
			if( where == m_values_per_snp.end() ) {
				throw genfile::BadArgumentError( "HighestValueSNPPicker::set_snps()", "snps", "SNP " + genfile::string_utils::to_string( getter(i) ) + " is not in the value map." ) ;
			}
			m_index_to_value_map[i] = where->second ;
			m_value_to_index_map.insert( std::make_pair( where->second, i )) ;
		}
		m_snp_comparator = comparator ;
	}

	std::size_t pick(
		std::deque< std::size_t > const& among_these
	) const {
		assert( among_these.size() > 0 ) ;
		ValueToIndexMap::reverse_iterator pick = m_value_to_index_map.rbegin() ;
		while(
			!std::binary_search(
				among_these.begin(),
				among_these.end(),
				pick->second,
				m_snp_comparator
			)
		) {
			++pick ;
			assert( pick != m_value_to_index_map.rend() ) ;
		}

        std::size_t chosen_snp = pick->second ;
		
		// The following line is the main optimisation which ensures we don't keep repeating work
		// we've already done.  As the algorithm proceeds, m_value_to_index_map gets smaller
		// and we just look at the top of it (rbegin) each time.
		// Warning! Use of this line assumes that on each invocation to this function,
		// the set of indices passed to this function does not contain any indices
		// that we picked on previous calls.
		m_value_to_index_map.erase( pick.base(), m_value_to_index_map.end() ) ;

        return chosen_snp ;
	}

	std::string display() const {
		return "HighestValueSNPPicker" ;
	} ;
	
	std::set< std::string > get_attribute_names() const {
		std::set< std::string > result ;
		result.insert( "rank" ) ;
		return result ;
	}
	
	std::map< std::string, genfile::VariantEntry > get_attributes( std::size_t chosen_snp ) const {
		std::map< std::string, genfile::VariantEntry > result ;
		result[ "rank" ] = m_index_to_value_map[ chosen_snp ] ;
		return result ;
	}
	
private:
	SnpToValueMap const m_values_per_snp ;
	IndexToValueMap m_index_to_value_map ;
	mutable ValueToIndexMap m_value_to_index_map ;
	SNPComparator m_snp_comparator ;
} ;


// ProximityTest encapsulates a criterion for two SNPs being 'too close together'.
// Given the context of a fixed list of SNPs (passed in using set_snps),
// ProximityTest supports the primitive remove_snps_too_close_to()
// which removes from the given list indices of SNPs that are deemed 'too close in the genome'
// to the first argument.
//
// For efficiency reasons, ProximityTest is endowed with a prepare() method
// which arranges a given list of indices in an order (e.g. by genomic position) that
// is most suitable for the implementation of this test.
class ProximityTest
{
public:
	typedef std::auto_ptr< ProximityTest > UniquePtr ;
	typedef boost::function< genfile::SNPIdentifyingData2 const& ( std::size_t ) > SNPGetter ;
public:
	virtual ~ProximityTest() {}
	virtual void set_snps( std::size_t number_of_snps, SNPGetter snps ) = 0 ;
	virtual void prepare( std::deque< std::size_t >* among_these ) const = 0 ;
	virtual void remove_snps_too_close_to( std::size_t chosen_snp_i, std::deque< std::size_t >* among_these ) const = 0 ;	
	virtual std::string display() const = 0 ;
	virtual std::set< std::string > get_attribute_names() const = 0 ;
	virtual std::map< std::string, genfile::VariantEntry > get_attributes( std::size_t ) const = 0 ;
} ;

class RecombinationDistanceProximityTest: public ProximityTest
{
public:
	RecombinationDistanceProximityTest(
		genfile::GeneticMap const& genetic_map,
		double minimum_distance_in_cM,
		std::size_t margin_in_bp
	):
		m_genetic_map( genetic_map ),
		m_minimum_distance_in_cM( minimum_distance_in_cM ),
		m_margin_in_bp( margin_in_bp )
	{
	}
	
	void set_snps( std::size_t number_of_snps, SNPGetter snps ) {
		m_snps.reserve( number_of_snps ) ;
		for( std::size_t i = 0; i < number_of_snps; ++i ) {
			m_snps.push_back( snps(i) ) ;
		}
	}
	
	std::string display() const {
		std::string result = "RecombinationDistanceProximityTest( " + genfile::string_utils::to_string( m_minimum_distance_in_cM ) + "cM" ;
		if( m_margin_in_bp > 0 ) {
			result += "+" + genfile::string_utils::to_string( m_margin_in_bp ) + "bp" ;
		}
		result += " )" ;
		return result ;
	}

	void prepare( std::deque< std::size_t >* among_these ) const {
		// Sort by chromosome / recombination distance
		//std::sort( among_these->begin(), among_these->end(), boost::bind( &RecombinationDistanceProximityTest::compare_recombination_positions, this, _1, _2 )) ;
		std::sort( among_these->begin(), among_these->end(), boost::bind( &RecombinationDistanceProximityTest::compare_physical_positions, this, _1, _2 )) ;
	}

	void remove_snps_too_close_to( std::size_t chosen_snp, std::deque< std::size_t >* among_these ) const {
		std::pair< genfile::Chromosome, double >
			lower_bound = get_recombination_position( m_snps[ chosen_snp ].get_position() ),
			upper_bound = lower_bound ;
		lower_bound.second -= m_minimum_distance_in_cM ;
		lower_bound.second = std::max( lower_bound.second, 0.0 ) ;
		upper_bound.second += m_minimum_distance_in_cM ;

		// Find an iterator to the lowest SNP in among_these
		// such that recombination position of the (SNP + margin)
		// is greater than or equal to lower_bound
		std::deque< std::size_t >::iterator
			lower_bound_i = std::lower_bound(
				among_these->begin(),
				among_these->end(),
				lower_bound,
				boost::bind(
			 		&RecombinationDistanceProximityTest::compare_to_recombination_position,
					this,	
					_1,
					_2,
					m_margin_in_bp
				)
		) ;

		// Find an iterator to the lowest SNP in among_these
		// such that recombination position of the (SNP - margin)
		// is greater than upper_bound
		std::deque< std::size_t >::iterator
			upper_bound_i = std::upper_bound(
				among_these->begin(),
				among_these->end(),
				upper_bound,
				boost::bind(
			 		&RecombinationDistanceProximityTest::compare_recombination_position_to,
					this,
					_1,
					_2,
					m_margin_in_bp
				)
		) ;
		// We should always have the chosen SNP itself.
		// assert( std::distance( lower_bound_i, upper_bound_i ) >= 1 ) ;

		among_these->erase( lower_bound_i, upper_bound_i ) ;
	}

	std::set< std::string > get_attribute_names() const {
		std::set< std::string > result ;
		result.insert( "region_lower_bp" ) ;
		result.insert( "region_upper_bp" ) ;
		result.insert( "region_lower_cm" ) ;
		result.insert( "region_upper_cm" ) ;
		return result ;
	}

	std::map< std::string, genfile::VariantEntry > get_attributes( std::size_t chosen_snp ) const {
		std::pair< genfile::Chromosome, double >
			lower_bound = get_recombination_position( m_snps[ chosen_snp ].get_position() ),
			upper_bound = lower_bound ;

		genfile::GenomePosition lower_bound_bp, upper_bound_bp ;

		lower_bound.second -= m_minimum_distance_in_cM ;
		if( lower_bound.second >= m_genetic_map.get_start_of_map_in_cM( lower_bound.first )) {
		 	lower_bound_bp = m_genetic_map.find_least_physical_position( lower_bound.first, lower_bound.second ) ;
		}
		else {
			lower_bound_bp = m_snps[ chosen_snp ].get_position() ;
		}
		lower_bound_bp.position() -= std::min( m_margin_in_bp, lower_bound_bp.position() ) ;

		upper_bound.second += m_minimum_distance_in_cM ;
		if( upper_bound.second <=  m_genetic_map.get_end_of_map_in_cM( lower_bound.first )) {
			upper_bound_bp = m_genetic_map.find_greatest_physical_position( upper_bound.first, upper_bound.second ) ;
			upper_bound_bp.position() += m_margin_in_bp ;
		}
		else {
			upper_bound_bp.chromosome() = upper_bound.first ;
			upper_bound_bp.position() = std::numeric_limits< int >::max() ;
		}

		double lower_bound_cM = get_recombination_position( lower_bound_bp ).second ;
		double upper_bound_cM = get_recombination_position( upper_bound_bp ).second ;

		std::map< std::string, genfile::VariantEntry > result ;
		result[ "region_lower_bp" ] = int( lower_bound_bp.position()  ) ;
		result[ "region_upper_bp" ] = int( upper_bound_bp.position() ) ;
		result[ "region_lower_cm" ] = lower_bound_cM ;
		result[ "region_upper_cm" ] = upper_bound_cM ;
		return result ;
	}

private:
	std::vector< genfile::SNPIdentifyingData2 > m_snps ;
	genfile::GeneticMap const& m_genetic_map ;
	double const m_minimum_distance_in_cM ;
	genfile::Position const m_margin_in_bp ;

private:
	std::pair< genfile::Chromosome, double > get_recombination_position(
		genfile::GenomePosition const& position
	) const {
		return std::make_pair(
			position.chromosome(),
			m_genetic_map.find_cM_from_beginning_of_chromosome_at_position( position )
		) ;
	}

	bool compare_recombination_positions( std::size_t a, std::size_t b ) const {
		return get_recombination_position( m_snps[a].get_position() ) < get_recombination_position( m_snps[b].get_position() ) ;
	}

	bool compare_physical_positions( std::size_t a, std::size_t b ) const {
		return m_snps[a].get_position() < m_snps[b].get_position() ;
	}
	
	// Return true if the recombination of the position of the snp with index a
	// plus the margin, specified in kb, has recombination position less than the given one.
	bool compare_to_recombination_position(
		std::size_t a,
		std::pair< genfile::Chromosome, double > const& chromosome_and_offset,
		genfile::Position margin_in_bp
	) const {
		genfile::GenomePosition position = m_snps[a].get_position() ;
		position.position() += margin_in_bp ;
		return get_recombination_position( position ) < chromosome_and_offset ;
	}

	// Return true if the given position is less than that of the position of
	// the snp with index b minus the margin.
	bool compare_recombination_position_to(
		std::pair< genfile::Chromosome, double > const& chromosome_and_offset,
		std::size_t b,
		genfile::Position margin_in_bp
	) const {
		genfile::GenomePosition position = m_snps[b].get_position() ;
		position.position() -= std::min( position.position(), margin_in_bp ) ;
		return chromosome_and_offset < get_recombination_position( position ) ;
	}

} ;

class PhysicalDistanceProximityTest: public ProximityTest
{
public:
	PhysicalDistanceProximityTest(
		double minimum_distance_in_base_pairs
	):
		m_minimum_distance_in_base_pairs( minimum_distance_in_base_pairs )
	{}
	
	void set_snps( std::size_t number_of_snps, SNPGetter snps ) {
		m_snps.reserve( number_of_snps ) ;
		for( std::size_t i = 0; i < number_of_snps; ++i ) {
			m_snps.push_back( snps(i) ) ;
		}
	}

	void prepare( std::deque< std::size_t >* among_these ) const {
		// Sort by chromosome / physical distance
		std::sort( among_these->begin(), among_these->end(), boost::bind( &PhysicalDistanceProximityTest::compare, this, _1, _2 )) ;
	}

	void remove_snps_too_close_to( std::size_t chosen_snp, std::deque< std::size_t >* among_these ) const {
		genfile::GenomePosition
			lower_bound = m_snps[ chosen_snp ].get_position(),
			upper_bound = lower_bound ;

		if( lower_bound.position() > m_minimum_distance_in_base_pairs ) {
			lower_bound.position() -= m_minimum_distance_in_base_pairs ;
		}
		else {
			lower_bound.position() = 0 ;
		}
		
		upper_bound.position() += m_minimum_distance_in_base_pairs ;

		std::deque< std::size_t >::iterator
			lower_bound_i = std::upper_bound(
				among_these->begin(),
				among_these->end(),
				lower_bound,
				boost::bind(
			 		&PhysicalDistanceProximityTest::compare_position_to_b,
					this,
					_1,
					_2
				)
		) ;

		std::deque< std::size_t >::iterator
			upper_bound_i = std::lower_bound(
				among_these->begin(),
				among_these->end(),
				upper_bound,
				boost::bind(
			 		&PhysicalDistanceProximityTest::compare_a_to_position,
					this,
					_1,
					_2
				)
		) ;

		among_these->erase( lower_bound_i, upper_bound_i ) ;
	}
	
	std::set< std::string > get_attribute_names() const {
		std::set< std::string > result ;
		result.insert( "region_lower_bp" ) ;
		result.insert( "region_upper_bp" ) ;
		return result ;
	}

	std::map< std::string, genfile::VariantEntry > get_attributes( std::size_t chosen_snp ) const {
		genfile::GenomePosition
			lower = m_snps[ chosen_snp ].get_position(),
			upper = m_snps[ chosen_snp ].get_position() ;
			
		lower.position() -= std::min( lower.position(), m_minimum_distance_in_base_pairs ) ;
		upper.position() += m_minimum_distance_in_base_pairs ;

		std::map< std::string, genfile::VariantEntry > result ;
		result[ "region_lower_bp" ] = int( lower.position() ) ;
		result[ "region_upper_bp" ] = int( upper.position() ) ;
		return result ;
	}

	
	std::string display() const {
		return "PhysicalDistanceProximityTest( " + genfile::string_utils::to_string( m_minimum_distance_in_base_pairs ) + "bp )" ;
	}
	
private:
	std::vector< genfile::SNPIdentifyingData2 > m_snps ;
	genfile::Position const m_minimum_distance_in_base_pairs ;
	
	bool compare( std::size_t a, std::size_t b ) const {
		return m_snps[ a ].get_position() < m_snps[ b ].get_position() ;
	}
	
	bool compare_a_to_position( std::size_t a, genfile::GenomePosition position ) const {
		return m_snps[ a ].get_position() < position ;
	}

	bool compare_position_to_b( genfile::GenomePosition position, std::size_t b ) const {
		return position < m_snps[ b ].get_position() ;
	}
	
} ;

namespace genes {
	struct Feature {
		Feature(
			std::string const& name,
			genfile::GenomePosition start,
			genfile::GenomePosition end
		):
			m_name( name ),
			m_start( start ),
			m_end( end )
		{
			assert( m_start.chromosome() == m_end.chromosome() ) ;
		}
			
		Feature( Feature const& other ):
			m_name( other.m_name ),
			m_start( other.m_start ),
			m_end( other.m_end )
		{}
		
		Feature& operator=( Feature const& other ) {
			m_name = other.m_name ;
			m_start = other.m_start ;
			m_end = other.m_end ;
			return *this ;
		}
		
		std::string const& name() const { return m_name ; }
		genfile::Chromosome const chromosome() const { return m_start.chromosome() ; }
		genfile::GenomePosition const start() const { return m_start ; }
		genfile::GenomePosition const end() const { return m_end ; }
		
	private:
		std::string m_name ;
		genfile::GenomePosition m_start ;
		genfile::GenomePosition m_end ;
	} ;

	// order from left-to-right along chromosomes by start position then by end position.
	struct CompareFeaturesByStart {
		bool operator()( Feature const& left, Feature const& right ) const {
			return( left.start() < right.start() ) ;
		}
	} ;
	
	bool operator<( Feature const& left, Feature const& right ) {
		return( left.start() < right.start() ) ;
	}

	namespace {
		struct start {} ;
		struct end {} ;
	}
	
	struct Genes: public boost::noncopyable {
	private:
		typedef std::set< Feature > ChromosomeGeneSet ;
		typedef std::map< genfile::Chromosome, ChromosomeGeneSet > GeneMap ;
	public:
		typedef std::auto_ptr< Genes > UniquePtr ;

	public:
		Genes() {}

		void add_gene(
			Feature const& feature,
			bool longest_transcript_only
		) {
			m_genes[ feature.chromosome() ].insert( feature ) ;
		}

		std::size_t number_of_genes() const { return m_genes.size() ; }
		std::vector< Feature const* > find_genes_in_region( genfile::Chromosome const chromosome, genfile::Position const lower, genfile::Position const upper ) {
			return find_genes_in_region( genfile::GenomePosition( chromosome, lower ), genfile::GenomePosition( chromosome, upper )) ;
		}

		std::vector< Feature const* > find_genes_in_region( genfile::GenomePosition const& lower, genfile::GenomePosition const& upper ) {
			assert( lower.chromosome() == upper.chromosome() ) ;
			// Genes are ordered by start then by end.
			// We first find the one-past-the end possible gene intersecting the region.
			// Then we walk leftwards until either
			// 1. we run out of genes
			// 2. we run off the end of the chromosome
			// 3. the end 
			
			std::vector< Feature const* > result ;
			if( m_genes.find( lower.chromosome() ) == m_genes.end() ) {
				return result ;
			}
			ChromosomeGeneSet const& chromosomeGenes = m_genes.at( lower.chromosome() ) ;

			std::vector< int > distances( chromosomeGenes.size() ) ;
			// make dummy region end feature for comparison.
			Feature regionEnd( "end", upper, upper ) ;
			// Find the first gene past-the end (or the end iterator)
			ChromosomeGeneSet::const_iterator i = std::upper_bound( chromosomeGenes.begin(), chromosomeGenes.end(), regionEnd ) ;
			if( i != chromosomeGenes.end() && i->start().position() <= upper.position() ) {
				result.push_back( &(*i) ) ;
			}
			if( i != chromosomeGenes.begin() ) {
				for( --i; i != chromosomeGenes.begin(); --i ) {
					if( i->end().position() > lower.position() ) {
						result.push_back( &(*i) ) ;
					}
				}
			}
			return result ;
		}
	private:
		GeneMap m_genes ;
	} ;
	
	Genes::UniquePtr load_genes_from_refGene(
		std::string const& filename,
		appcontext::UIContext::ProgressContext& progress_context,
		bool longest_transcript_only = false
	) {
		statfile::BuiltInTypeStatSource::UniquePtr source( new statfile::DelimitedStatSource( filename, "\t" ) ) ;
		
		std::size_t chromColumn = source->index_of_column( "chrom" ) ;
		std::size_t txStartColumn = source->index_of_column( "txStart" ) ;
		std::size_t txEndColumn = source->index_of_column( "txEnd" ) ;
		std::size_t name2Column = source->index_of_column( "name2" ) ;

		//std::cerr << "chromColumn = " << chromColumn << ", txStartColumn = " << txStartColumn << ", txEndColumn = " << txEndColumn << ", name2Column = " << name2Column << ", number of columns is " << source->number_of_columns() << ".\n" ;
		int bin ;
		genfile::Position txStart ;
		genfile::Position txEnd ;
		std::string chrom ;
		std::string name2 ;
		
		Genes::UniquePtr result( new Genes() ) ;
		while( (*source) >> bin ) {
			//std::cerr << "Reading a gene...\n" ;
			(*source)
				>> statfile::ignore( chromColumn - 1 )
				>> chrom
				>> statfile::ignore( txStartColumn - chromColumn - 1 )
				>> txStart
				>> statfile::ignore( txEndColumn - txStartColumn - 1 )
				>> txEnd
				>> statfile::ignore( name2Column - txEndColumn - 1 )
				>> name2 ;
			if( !(*source)) {
				throw genfile::MalformedInputError( source->get_source_spec(), "Malformed refGene-format file", source->number_of_rows_read() ) ;
			}
			
			// UCSC coordinates are 0-based.  Map them here.
			++txStart ;
			++txEnd ;

			// reformat the chromosome by getting rid of the 'chr'.
			if( chrom.size() > 3 && chrom.substr(0, 3 ) == "chr" ) {
				chrom = chrom.substr( 3, chrom.size() ) ;
			}
			chrom = genfile::string_utils::replace_all( chrom, "chr", "" ) ;

			result->add_gene(
				Feature( name2, genfile::GenomePosition( chrom, txStart ), genfile::GenomePosition( chrom, txEnd ) ),
				longest_transcript_only
			) ;

			(*source) >> statfile::ignore_all() ;
			progress_context.notify_progress( source->number_of_rows_read(), source->number_of_rows() ) ;
		}
		
		return result ;
	}
}

class InthinneratorApplication: public appcontext::ApplicationContext
{
public:
	InthinneratorApplication( int argc, char** argv ):
		ApplicationContext( globals::program_name, std::auto_ptr< appcontext::OptionProcessor >( new InthinneratorOptionProcessor() ), argc, argv, "-log" )
	{
		process() ;
	}

private:
	enum DistanceType { ePhysical = 0, eRecombination = 1 } ;
	ProximityTest::UniquePtr m_proximity_test ;
	SNPPicker::UniquePtr m_snp_picker ;
	bool m_write_incl_list, m_write_excl_list ;

private:
	void process() {
		try {
			unsafe_process() ;
		}
		catch( genfile::FileNotFoundError const& e ) {
			get_ui_context().logger() << "!! Error (" << e.what() << "): the file \"" << e.filespec() << "\" could not be found.\n" ;
			throw appcontext::HaltProgramWithReturnCode( -1 ) ;
		}
		catch( genfile::ChromosomeNotInMapError const& e ) {
			get_ui_context().logger() << "!! Error (" << e.what() << "): chromosome " << e.chromosome() << " was not in the genetic map.\n" ;
			throw appcontext::HaltProgramWithReturnCode( -1 ) ;
		}
		catch( genfile::MalformedInputError const& e ) {
			get_ui_context().logger() << "!! Error (" << e.what() << "): " << e.format_message() << ".\n" ;
			throw appcontext::HaltProgramWithReturnCode( -1 ) ;
		}
		catch( genfile::BadArgumentError const& e ) {
			get_ui_context().logger() << "!! Error (" << e.what() << "): Bad argument ("
				<< e.arguments() << ") to function " << e.function() << ": "
				<< e.format_message() << "\n" ;
			throw appcontext::HaltProgramWithReturnCode( -1 ) ;
		}
		catch( statfile::FileNotOpenedError const& e ) {
			get_ui_context().logger() << "!! Error (" << e.what() << "): The file \"" << e.filename() << "\" could not be opened.\n" ;
			throw appcontext::HaltProgramWithReturnCode( -1 ) ;
		}
		catch( FileNotOpenedError const& e ) {
			get_ui_context().logger() << "!! Error (" << e.what() << "): The file \"" << e.filename() << "\" could not be opened.\n" ;
			throw appcontext::HaltProgramWithReturnCode( -1 ) ;
		}
		catch( db::StatementStepError const& e ) {
			get_ui_context().logger() << "!! Error (" << e.what() << "): error with SQL: " << e.sql() << "\n"
				<< e.description() << ".\n" ;
			throw appcontext::HaltProgramWithReturnCode( -1 ) ;
		}
	}
	
	void unsafe_process() {
		if( !options().check( "-o" ) && ! options().check( "-odb" ) ) {
			get_ui_context().logger() << "You must supply either -o or -odb.\n" ;
			throw appcontext::HaltProgramWithReturnCode( -1 ) ;
		}
		genfile::GeneticMap::UniquePtr map ;
		if( options().check( "-map" )) {
			load_genetic_map() ;
			get_ui_context().logger() << "Loaded: " << map->get_summary() << "\n";
		}

		std::vector< TaggedSnp > snps = get_list_of_snps( options().get_value< std::string >( "-g" ), map ) ;
		boost::optional< std::vector< double > > recombination_offsets ;
		if( map.get() ) {
			recombination_offsets = get_recombination_offsets(
				snps.size(),
				boost::bind( &impl::get_snp, boost::cref( snps ), _1 ),
				*map
			) ;
		}
		m_proximity_test = get_proximity_test( snps, map ) ;
		m_snp_picker = get_snp_picker( snps ) ;

		// Write an inclusion list either if user specified -write-incl-list, or didn't specify any output.
		m_write_incl_list = !options().check_if_option_was_supplied( "-suppress-included" ) ;
		// Write an exclusion list either if user specified -write-excluded, or didn't specify any output.
		m_write_excl_list = !options().check_if_option_was_supplied( "-suppress-excluded" ) ;

		genes::Genes::UniquePtr genes ;
		if( options().check( "-genes" )) {
			std::string const filename = options().get< std::string > ( "-genes" ) ;
			appcontext::UIContext::ProgressContext progress_context = get_ui_context().get_progress_context( "Loading genes from \"" + filename + "\"" ) ;
			genes = genes::load_genes_from_refGene( filename, progress_context, options().get< bool >( "-take-longest-transcript" ) ) ;
		}
		
		write_summary( map, snps, *m_proximity_test, *m_snp_picker ) ;
		
		perform_thinnings( snps, recombination_offsets, genes ) ;
	}
	
	genfile::GeneticMap::UniquePtr load_genetic_map() {
		UIContext::ProgressContext progress_context = get_ui_context().get_progress_context( "Loading genetic map" ) ;
		genfile::GeneticMap::UniquePtr map(
			new genfile::FromFilesGeneticMap(
				genfile::wildcard::find_files_by_chromosome(
					options().get_value< std::string >( "-map" ),
					genfile::wildcard::eALL_CHROMOSOMES
				),
				boost::ref( progress_context )
			)
		) ;
		return map ;
	}
	
	genfile::CommonSNPFilter::UniquePtr get_filter() const {
		genfile::CommonSNPFilter::UniquePtr filter(
			new genfile::CommonSNPFilter
		) ;
		if( options().has_value( "-excl-rsids" )) {
			std::vector< std::string > filenames = options().get_values< std::string >( "-excl-rsids" ) ;
			for( std::size_t i = 0; i < filenames.size(); ++i ) {
				filter->exclude_snps_in_file( filenames[i], genfile::CommonSNPFilter::RSIDs ) ;
			}
		}
		if( options().has_value( "-incl-rsids" )) {
			filter->exclude_snps_not_in_file( options().get< std::string >( "-incl-rsids" ), genfile::CommonSNPFilter::RSIDs ) ;
		}
		if( options().has_value( "-excl-snpids" )) {
			std::vector< std::string > filenames = options().get_values< std::string >( "-excl-rsids" ) ;
			for( std::size_t i = 0; i < filenames.size(); ++i ) {
				filter->exclude_snps_in_file( filenames[i], genfile::CommonSNPFilter::SNPIDs ) ;
			}
		}
		if( options().has_value( "-incl-snpids" )) {
			filter->exclude_snps_not_in_file( options().get< std::string >( "-incl-snpids" ), genfile::CommonSNPFilter::SNPIDs ) ;
		}

		if( options().check_if_option_was_supplied( "-incl-range" )) {
			std::vector< std::string > specs = options().get_values< std::string >( "-incl-range" ) ;
			for ( std::size_t i = 0; i < specs.size(); ++i ) {
				filter->include_snps_in_range(
					genfile::GenomePositionRange::parse( specs[i] )
				) ;
			}
		}
		
		if( options().check_if_option_was_supplied( "-excl-range" )) {
			std::vector< std::string > specs = options().get_values< std::string >( "-excl-range" ) ;
			for ( std::size_t i = 0; i < specs.size(); ++i ) {
				filter->exclude_snps_in_range(
					genfile::GenomePositionRange::parse( specs[i] )
				) ;
			}
		}
		
		return filter ;
	}

	std::vector< TaggedSnp > get_list_of_snps( std::string const& filename, genfile::GeneticMap::UniquePtr const& map ) const {
		std::vector< TaggedSnp > filtered_snps ;
		try {
			filtered_snps = get_list_of_snps_from_gen_file( filename ) ;
		}
		catch( std::exception const& e ) {
			get_ui_context().logger() << "File \"" << filename << "\" is not a GEN-style file.  Trying text format with columns SNPID rsid chromosome position...\n" ;
			// Not a GEN format file.  Try a different format.
			boost::optional< std::string > tag_column = ( options().check( "-tag-column" ) ? options().get< std::string >( "-tag-column" ) : boost::optional< std::string >() ) ;
			filtered_snps = get_list_of_snps_from_text_file( filename, tag_column ) ;
		}
		
		genfile::CommonSNPFilter::UniquePtr filter = get_filter() ;
		if( map.get() ) {
			filter->exclude_chromosomes_not_in_set( map->get_chromosomes() ) ;
		}

		std::vector< std::size_t > indices_of_included_snps = filter->get_indices_of_filtered_in_snps(
			filtered_snps.size(),
			boost::bind(
				&impl::get_snp,
				boost::cref( filtered_snps ),
				_1
			)
		) ;

		filtered_snps = genfile::utility::select_entries( filtered_snps, indices_of_included_snps ) ;

		return filtered_snps ;
	}

	std::vector< TaggedSnp > get_list_of_snps_from_gen_file( std::string const& filename ) const {
		std::vector< TaggedSnp > filtered_snps ;
		genfile::SNPDataSource::UniquePtr source ;
		{
			UIContext::ProgressContext progress_context = get_ui_context().get_progress_context( "Opening genotype files" ) ;
			genfile::vcf::MetadataParser::Metadata metadata ;
			if( options().check( "-metadata" ) ) {
				metadata = genfile::vcf::StrictMetadataParser(
					options().get_value< std::string >( "-metadata" )
				).get_metadata() ;
			}
			source.reset(
				genfile::SNPDataSourceChain::create(
					genfile::wildcard::find_files_by_chromosome( filename, genfile::wildcard::eALL_CHROMOSOMES ),
					metadata,
					"guess",
					boost::ref( progress_context )
				).release()
			) ;
		}
		{
			UIContext::ProgressContext progress_context = get_ui_context().get_progress_context( "Loading SNP list" ) ;
			source->list_snps( boost::bind( &impl::add_snp, &filtered_snps, _1, boost::optional< std::string >() ), boost::ref( progress_context ) ) ;
		}
		return filtered_snps ;
	}
	
	std::vector< TaggedSnp > get_list_of_snps_from_text_file(
		std::string const& filename,
		boost::optional< std::string > tag_column
	) const {
		std::vector< TaggedSnp > filtered_snps ;
		UIContext::ProgressContext progress_context = get_ui_context().get_progress_context( "Loading SNP list" ) ;

		statfile::BuiltInTypeStatSourceChain::UniquePtr chain(
			statfile::BuiltInTypeStatSourceChain::open(
				genfile::wildcard::find_files_by_chromosome( filename, genfile::wildcard::eALL_CHROMOSOMES )
			)
		) ;
		
		boost::optional< std::size_t > tag_column_index
			= tag_column ? chain->index_of_column( tag_column.get() ): boost::optional< std::size_t >() ;
		;

		std::string SNPID ;
		std::string rsid ;
		std::string chromosome ;
		genfile::Position position ;
		std::string alleleA = "?" ;
		std::string alleleB = "?" ;
		std::string tag_string ;
		while(
			(*chain)
				>> SNPID
				>> rsid
				>> chromosome
				>> position
		) {
			alleleA.clear() ;
			alleleB.clear() ;
			boost::optional< std::string > tag ;

			if( chain->has_column( "alleleA" ) ) {
				(*chain) >> alleleA ;
			}

			if( chain->has_column( "alleleB" ) ) {
				(*chain) >> alleleB ;
			}

			if( tag_column_index ) {
				(*chain)
					>> statfile::ignore( tag_column_index.get() - chain->current_column() )
					>> tag_string ;
				tag = tag_string ;
			}
			
			// std::cerr << "Read SNP: " << SNPID << " " << rsid << " " << chromosome << " " << position << ".\n" ;
			filtered_snps.push_back(
				TaggedSnp(
					genfile::SNPIdentifyingData2(
						SNPID,
						rsid,
						genfile::GenomePosition( chromosome, position ),
						alleleA,
						alleleB	
					),
					tag
				)
			) ;
			(*chain) >> statfile::ignore_all() ;
			progress_context.notify_progress( chain->number_of_rows_read(), chain->number_of_rows() ) ;
		}
		
		return filtered_snps ;
	}

	std::vector< double > get_recombination_offsets(
		std::size_t number_of_snps,
		boost::function< genfile::SNPIdentifyingData2 const& ( std::size_t ) > snps,
		genfile::GeneticMap const& map
	) const {
		UIContext::ProgressContext progress_context = get_ui_context().get_progress_context( "Computing recombination positions" ) ;
		std::vector< double > offsets( number_of_snps ) ;
		for( std::size_t i = 0; i < number_of_snps; ++i ) {
			offsets[i] = map.find_cM_from_beginning_of_chromosome_at_position( snps(i).get_position() ) ;
			progress_context( i+1, number_of_snps ) ;
		}
		return offsets ;
	}

	ProximityTest::UniquePtr get_proximity_test(
		std::vector< TaggedSnp > const& snps,
		genfile::GeneticMap::UniquePtr const& genetic_map
	) const {
		ProximityTest::UniquePtr test ;
		std::string distance_spec = options().get< std::string >( "-min-distance" ) ;
		std::vector< std::string > elts = genfile::string_utils::split_and_strip( distance_spec, "+" ) ;

		if( elts.size() < 1 || elts.size() > 2 ) {
			genfile::BadArgumentError( "InthinneratorApplication::get_proximity_test", "-min-distance = \"" + distance_spec + "\"" ) ;
		}
		
		std::vector< std::pair< double, std::string > > parsed_elts( elts.size() ) ;
		for( std::size_t i = 0; i < elts.size(); ++i ) {
			parsed_elts[i] = parse_distance_piece( elts[i] ) ;
		}

		if( parsed_elts.size() == 2 ) {
			// Assume <x>cM+20bp format.
			if( parsed_elts[0].second != "cM" || parsed_elts[1].second != "bp" ) {
				throw genfile::BadArgumentError( "InthinneratorApplication::get_proximity_test", "-min-distance = \"" + distance_spec + "\"" ) ;
			}
			if( !genetic_map.get() ) {
				throw genfile::BadArgumentError( "InthinneratorApplication::get_proximity_test", "genetic_map", "No genetic map was specified, use -map to specify one." ) ;
			}
			test.reset(
				new RecombinationDistanceProximityTest(
					*genetic_map,
					parsed_elts[0].first,
					parsed_elts[1].first
				)
			) ;
		}
		else if( parsed_elts[0].second == "cM" ) {
			if( !genetic_map.get() ) {
				throw genfile::BadArgumentError( "InthinneratorApplication::get_proximity_test", "genetic_map", "No genetic map was specified, use -map to specify one." ) ;
			}
			test.reset(
				new RecombinationDistanceProximityTest(
					*genetic_map,
					parsed_elts[0].first,
					0
				)
			) ;
		}
		else if( parsed_elts[0].second == "bp" ) {
			test.reset(
				new PhysicalDistanceProximityTest(
					parsed_elts[0].first
				)
			) ;
		}
		else {
			throw genfile::BadArgumentError( "InthinneratorApplication::get_proximity_test", "-min-distance = \"" + distance_spec + "\"" ) ;
		}
		
		test->set_snps(
			snps.size(),
			boost::bind( &impl::get_snp, boost::cref( snps ), _1 )
		) ;
		return test ;
	}

	std::pair< double, std::string > parse_distance_piece( std::string const piece ) const {
		std::pair< double, std::string > result ;
		if( piece.compare( piece.size() - 2, 2, "cM" ) == 0 ) {
			result.first = genfile::string_utils::to_repr< double >( piece.substr( 0, piece.size() - 2 )) ;
			result.second = "cM" ;
		}
		else if( piece.compare( piece.size() - 2, 2, "bp" ) == 0 ) {
			result.first = genfile::string_utils::to_repr< double >( piece.substr( 0, piece.size() - 2 )) ;
			result.second = "bp" ;
		}
		else if( piece.compare( piece.size() - 2, 2, "kb" ) == 0 ) {
			result.first = 1000.0 * genfile::string_utils::to_repr< double >( piece.substr( 0, piece.size() - 2 )) ;
			result.second = "bp" ;
		}
		else if( piece.compare( piece.size() - 2, 2, "Mb" ) == 0 ) {
			result.first = 1000000.0 * genfile::string_utils::to_repr< double >( piece.substr( 0, piece.size() - 2 )) ;
			result.second = "bp" ;
		}
		return result ;
	}

	SNPPicker::UniquePtr get_snp_picker(
		std::vector< TaggedSnp > const& snps
	) const {
		SNPPicker::UniquePtr picker ;
		if( options().check_if_option_was_supplied( "-rank" )) {
			std::map< genfile::SNPIdentifyingData2, double > value_map = load_ranks(
				options().get_value< std::string >( "-rank" ),
				options().get_value< std::string >( "-rank-column" ),
				genfile::string_utils::split_and_strip_discarding_empty_entries( options().get< std::string >( "-missing-code" ), ",", " \t" )
			) ;
			picker.reset( new HighestValueSNPPicker( value_map )) ;
		}
		else {
			std::string strategy = options().get< std::string >( "-strategy" ) ;
			if( strategy == "random" ) {
				picker.reset( new RandomSNPPicker ) ;
			}
			else if( strategy == "random_by_position" ) {
				std::pair< double, std::string > bin_size = parse_distance_piece( options().get< std::string >( "-bin-size" ) ) ;
				if( bin_size.second != "bp" ) {
					throw genfile::BadArgumentError(
						"InthinneratorApplication::get_snp_picker()",
						"option -bin-size",
						"Value must be expressed in bp, kb or Mb, not in units of recombination distance."
					) ;
				}
				std::vector< genfile::GenomePositionRange > ranges = compute_genomic_coverage( snps, bin_size.first ) ;
				picker.reset( new RandomPositionSNPPicker( ranges, bin_size.first )) ;
			}
			else if( strategy == "first" ) {
				picker.reset( new FirstAvailableSNPPicker ) ;
			}
			else {
				throw genfile::BadArgumentError( "InthinneratorApplication::get_snp_picker", "-strategy = \"" + strategy + "\"" ) ;
			}
		}
		return picker ;
	}

	std::vector< genfile::GenomePositionRange > compute_genomic_coverage(
		std::vector< TaggedSnp > snps,
		std::size_t max_distance
	) const {
		std::sort( snps.begin(), snps.end() ) ;

		std::vector< genfile::GenomePositionRange > result ;
		genfile::Position currentStart, currentEnd ;
		genfile::Chromosome currentChromosome ;
		bool inRange = false ;
		for( std::size_t i = 0; i < snps.size(); ++i ) {
			if( inRange ) {
				if( snps[i].snp().get_position().chromosome() != currentChromosome || ( snps[i].snp().get_position().position() - currentEnd ) > max_distance ) {
					if( currentEnd > currentStart ) { // ignore 1-SNP ranges
						result.push_back(
							genfile::GenomePositionRange(
								currentChromosome,
								currentStart,
								currentEnd
							)
						) ;
					}
					inRange = false ;
				} else {
					currentEnd = snps[i].snp().get_position().position() ;
				}
			}
			if( !inRange ) {
				currentStart = currentEnd = snps[i].snp().get_position().position() ;
				currentChromosome = snps[i].snp().get_position().chromosome() ;
				inRange = true ;
			}
		}
		assert( inRange ) ;
		result.push_back(
			genfile::GenomePositionRange(
				currentChromosome,
				currentStart,
				currentEnd
			)
		) ;
		return result ;
	}

	std::map< genfile::SNPIdentifyingData2, double > load_ranks(
		std::string const& filename,
		std::string const& column_name,
		std::vector< std::string > const& missing_values
	) const {
		statfile::BuiltInTypeStatSourceChain::UniquePtr chain(
			statfile::BuiltInTypeStatSourceChain::open(
				genfile::wildcard::find_files_by_chromosome(
					filename
				)
			)
		) ;

		std::size_t rank_column_index ;
		double sign ;
		
		if( chain->has_column( column_name )) {
			rank_column_index = chain->index_of_column( column_name ) ;
			sign = 1.0 ;
		}
		else if( column_name.size() > 0 && column_name[0] == '-' ) {
			rank_column_index = chain->index_of_column( column_name.substr( 1, column_name.size() ) ) ;
			sign = -1.0 ;
		}
		else {
			throw genfile::BadArgumentError( "load_ranks()", "column_name=\"" + column_name + "\"" ) ;
		}
		
		return load_ranks( *chain, rank_column_index, sign, std::set< std::string >( missing_values.begin(), missing_values.end() ) ) ;
	}

	std::map< genfile::SNPIdentifyingData2, double > load_ranks(
		statfile::BuiltInTypeStatSourceChain& chain,
		std::size_t rank_column_index,
		double sign,
		std::set< std::string > const& missing_values
	) const {
		if( rank_column_index < 4 ) {
			throw genfile::BadArgumentError( "InthinneratorApplication::load_ranks()", "rank_column_index" ) ;
		}

		UIContext::ProgressContext progress_context = get_ui_context().get_progress_context( "Loading ranks" ) ;
		
		std::map< genfile::SNPIdentifyingData2, double > result ;
		std::string SNPID ;
		std::string rsid ;
		std::string chromosome ;
		genfile::Position position ;
		std::string alleleA = "?", alleleB = "?" ;
		std::string rank_string ;
		double rank ;
		while(
			chain
				>> SNPID
				>> rsid
				>> chromosome
				>> position
		) {
			if( chain.has_column( "alleleA" ) ) {
				chain >> alleleA ;
			}
			if( chain.has_column( "alleleB" ) ) {
				chain >> alleleB ;
			}
			chain >> statfile::ignore( rank_column_index - chain.current_column() ) >> rank_string ;
			if( missing_values.find( rank_string ) != missing_values.end() ) {
				rank = std::numeric_limits< double >::quiet_NaN() ;
			}
			else {
				rank = genfile::string_utils::to_repr< double >( rank_string ) ;
			}
			result[ genfile::SNPIdentifyingData2( SNPID, rsid, genfile::GenomePosition( chromosome, position ), alleleA, alleleB ) ] = sign * rank ;
			chain >> statfile::ignore_all() ;
			progress_context.notify_progress( chain.number_of_rows_read(), chain.number_of_rows() ) ;
		}
		return result ;
	}

	void write_summary(
		genfile::GeneticMap::UniquePtr const& map,
		std::vector< TaggedSnp > const& snps_in_source,
		ProximityTest const& proximity_test,
		SNPPicker const& snp_picker
	) const {
		get_ui_context().logger() << "Genetic map file(s): "
			<< ( options().check( "-map" ) ? options().get_value< std::string >( "-map" ) : std::string( "(none specified)" ) )
			<< ".\n" ;
		get_ui_context().logger() << "Genotype file(s): " << options().get_value< std::string >( "-g" ) << ".\n" ;

		if( map.get() ) {
			get_ui_context().logger() << "Genetic map covers these chromosome(s): " ;
			std::set< genfile::Chromosome > chromosomes = map->get_chromosomes() ;
			for( std::set< genfile::Chromosome >::const_iterator i = chromosomes.begin(); i != chromosomes.end(); ++i ) {
				get_ui_context().logger() << *i << " " ;
			}
			get_ui_context().logger() << "\n" ;
		}
		
		get_ui_context().logger() << "- There are " << snps_in_source.size() << " genotyped SNPs in these chromosomes." ;
		if( snps_in_source.empty() ) {
			get_ui_context().logger() << ".\n" ;
		}
		else {
			get_ui_context().logger() << "  The first few are:\n" ;
			for( std::size_t i = 0; i < std::min( snps_in_source.size(), std::size_t(5) ); ++i ) {
				get_ui_context().logger() << "  " << snps_in_source[i].snp() ;
				if( snps_in_source[i].tag() ) {
					get_ui_context().logger() << " (tag=\"" << snps_in_source[i].tag().get() << "\")...\n" ;
				}
			}
			get_ui_context().logger() << "  .\n  .\n  .\n" ;
		}
		
		get_ui_context().logger()
			<< "Throwing out SNPs based on proximity test: "
			<< proximity_test.display() << ".\n" ;
		
		get_ui_context().logger()
			<< "Picking SNPs using SNP picker: "
			<< snp_picker.display() << ".\n" ;
		
	}

	void perform_thinnings(
		std::vector< TaggedSnp > const& snps,
		boost::optional< std::vector< double > > const& recombination_offsets,
		genes::Genes::UniquePtr const& genes
	) const {
		std::size_t const N = options().get_value< std::size_t >( "-N" ) ;
		assert( N > 0 ) ;
		std::size_t const start_N = options().get_value< std::size_t >( "-start-N" ) ;
		std::size_t const max_num_picks = options().check( "-max-picks" ) ? options().get_value< std::size_t >( "-max-picks" ) : snps.size() ;
		std::size_t const number_of_digits = std::max( std::size_t( std::log10( N + start_N ) ) + 1, std::size_t( 3u )) ;

		std::vector< boost::optional< std::string > > pick_tags ;
		if( options().check( "-match-tag" )) {
			std::vector< TaggedSnp > snps_to_match = get_list_of_snps_from_text_file(
				options().get< std::string >( "-match-tag" ),
				options().get< std::string >( "-tag-column" )
			) ;

			pick_tags.resize( snps_to_match.size() ) ;
			for( std::size_t i = 0; i < snps_to_match.size(); ++i ) {
				pick_tags[i] = snps_to_match[i].tag() ;
			}
		} else {
			pick_tags.resize( max_num_picks, boost::optional< std::string >() ) ;
		}

		std::string filename_stub = options().get< std::string >( "-o" ) ;
		
		std::string const formatstring = ( boost::format( filename_stub + ".%%0%dd" ) % number_of_digits ).str() ;
		std::cerr << "FORMAT STRING: " << formatstring << ".\n" ;
		for( std::size_t i = start_N; i < (start_N+N); ++i ) {
			get_ui_context().logger() << "Picking " << (i+1) << " of " << N << "..." ;
			std::set< std::size_t > picked_snps = pick_snps( snps, pick_tags ) ;
			get_ui_context().logger() << picked_snps.size() << " SNPs picked.\n" ;
			write_output( i, snps, recombination_offsets, picked_snps, genes, ( boost::format( formatstring ) % i ).str() ) ;
		}
	}

	std::set< std::size_t > pick_snps(
		std::vector< TaggedSnp > const& snps,
		std::vector< boost::optional< std::string > > const& pick_tags
	) const {
		// Create a list of all SNPs and sort it in the way preferred by the proximity test
		std::deque< std::size_t > all_snp_indices(
			boost::counting_iterator< std::size_t >( 0 ),
			boost::counting_iterator< std::size_t >( snps.size() )
		) ;
		m_proximity_test->prepare( &all_snp_indices ) ;
		CompareByOrder comparator( all_snp_indices ) ;

		// Tell the SNP picker our full list of SNPs
		m_snp_picker->set_snps(
			snps.size(),
			boost::bind( &impl::get_snp, snps, _1 ),
			boost::cref( comparator )
		) ;

		// Create per-tag subsets of the full list
		typedef std::map< boost::optional< std::string >, std::deque< std::size_t > > SnpsByTag ;
		SnpsByTag remaining_snps_by_tag ;
		for( std::size_t i = 0; i < all_snp_indices.size(); ++i ) {
			remaining_snps_by_tag[ snps[all_snp_indices[i]].tag() ].push_back( all_snp_indices[i] ) ;
		}

		return pick_snps( snps, pick_tags, &remaining_snps_by_tag ) ;
	}

	std::set< std::size_t > pick_snps(
		std::vector< TaggedSnp > const& snps,
		std::vector< boost::optional< std::string > > const& pick_tags,
		std::map< boost::optional< std::string >, std::deque< std::size_t > >* remaining_snps_by_tag
	) const {
		typedef std::map< boost::optional< std::string >, std::deque< std::size_t > > SnpsByTag ;
		std::set< std::size_t > result ;
		if( pick_tags.size() > 0 ) {
			UIContext::ProgressContext progress_context = get_ui_context().get_progress_context( "Picking SNPs" ) ;

			std::size_t remaining_snp_count = 0 ;
			SnpsByTag::const_iterator next_tag_i = remaining_snps_by_tag->find( pick_tags[0] ) ;
			while(
				result.size() < pick_tags.size()
				&& ( next_tag_i = remaining_snps_by_tag->find( pick_tags[ result.size() ] ) ) != remaining_snps_by_tag->end()
				&& next_tag_i->second.size() > 0
			) {
				std::size_t picked_snp = m_snp_picker->pick( next_tag_i->second ) ;
				remaining_snp_count = 0 ;
				for( SnpsByTag::iterator i = remaining_snps_by_tag->begin(); i != remaining_snps_by_tag->end(); ++i ) {
					m_proximity_test->remove_snps_too_close_to( picked_snp, &(i->second) ) ;
					remaining_snp_count += i->second.size() ;
				}
				result.insert( picked_snp ) ;
				progress_context.notify_progress( snps.size() - remaining_snp_count, snps.size() ) ;
			}

			if( next_tag_i == remaining_snps_by_tag->end() ) {
				throw genfile::BadArgumentError(
					"InthinneratorApplication::pick_snps()",
					"pick_tags",
					"A tag (\"" + ( pick_tags[ result.size() ] ? pick_tags[ result.size() ].get() : "NA" ) + "\") was supplied for which no SNP was available."
				) ;
			}

			progress_context.notify_progress( snps.size() - remaining_snp_count, snps.size() ) ;
		}
		return result ;
	}
	
	void write_output(
		std::size_t const iteration,
		std::vector< TaggedSnp > const& snps,
		boost::optional< std::vector< double > > const& recombination_offsets,
		std::set< std::size_t > const& picked_snps,
		genes::Genes::UniquePtr const& genes,
		std::string const& filename
	) const {
		std::vector< std::string > output_columns = get_output_columns() ;
		// remove crud we don't want in the output.
		output_columns.erase( std::remove( output_columns.begin(), output_columns.end(), "rsid" ), output_columns.end() ) ;
		output_columns.erase( std::remove( output_columns.begin(), output_columns.end(), "SNPID" ), output_columns.end() ) ;
		output_columns.erase( std::remove( output_columns.begin(), output_columns.end(), "chromosome" ), output_columns.end() ) ;
		output_columns.erase( std::remove( output_columns.begin(), output_columns.end(), "position" ), output_columns.end() ) ;
		output_columns.erase( std::remove( output_columns.begin(), output_columns.end(), "allele1" ), output_columns.end() ) ;
		output_columns.erase( std::remove( output_columns.begin(), output_columns.end(), "allele2" ), output_columns.end() ) ;

		qcdb::Storage::SharedPtr storage = qcdb::FlatFileOutputter::create_shared(
			filename,
			options().get< std::string >( "-analysis-name" ),
			options().get_values_as_map()
		) ;

		// Setup output columns
		if( options().check( "-match-tag" )) {
			storage->add_variable( "tag" ) ;
		}
		storage->add_variable( "iteration" ) ;
		storage->add_variable( "result" ) ;
		for( std::size_t i = 0; i < output_columns.size(); ++i ) {
			storage->add_variable( output_columns[i] ) ;
		}
		if( genes.get() ) {
			storage->add_variable( "nearest_gene_in_region" ) ;
			storage->add_variable( "distance_to_nearest_gene_in_region" ) ;
			storage->add_variable( "all_genes_in_region" ) ;
		}

		if( m_write_incl_list ) {
			UIContext::ProgressContext progress_context = get_ui_context().get_progress_context( "Writing included SNPs to \"" + filename + "\"" ) ;
			write_output(
				"included",
				iteration,
				snps,
				recombination_offsets,
				picked_snps,
				genes,
				*storage,
				output_columns,
				progress_context
			) ;
		}

		if( m_write_excl_list ) {
			std::set< std::size_t > unpicked_snps ;
			std::set_difference(
				boost::counting_iterator< std::size_t >( 0 ),
				boost::counting_iterator< std::size_t >( snps.size() ),
				picked_snps.begin(),
				picked_snps.end(),
				std::inserter( unpicked_snps, unpicked_snps.end() )
			) ;

			UIContext::ProgressContext progress_context = get_ui_context().get_progress_context( "Writing excluded SNPs to \"" + filename + "\"" ) ;
			write_output(
				"excluded",
				iteration,
				snps,
				recombination_offsets,
				unpicked_snps,
				genes,
				*storage,
				output_columns,
				progress_context
			) ;
		}
	}

	void write_output(
		std::string const& result,
		std::size_t const& iteration,
		std::vector< TaggedSnp > const& snps,
		boost::optional< std::vector< double > > const& recombination_offsets,
		std::set< std::size_t > const& indices_of_snps_to_output,
		genes::Genes::UniquePtr const& genes,
		qcdb::Storage& storage,
		std::vector< std::string > const& output_columns,
		UIContext::ProgressContext& progress_context
	) const {
		using genfile::string_utils::to_string ;
		progress_context.notify_progress( 0, snps.size() ) ;
		bool const by_tag = options().check( "-match-tag" ) ;
		
		std::size_t snp_index = 0 ;
		for(
			std::set< std::size_t >::const_iterator i = indices_of_snps_to_output.begin() ;
			i != indices_of_snps_to_output.end() ; 
			++i, ++snp_index
		) {
			if( by_tag ) {
				storage.store_per_variant_data(
					snps[*i].snp(),
					"tag",
					snps[*i].tag() ? genfile::VariantEntry( snps[*i].tag().get() ) : genfile::VariantEntry()
				) ;
			}

			storage.store_per_variant_data(
				snps[*i].snp(),
				"iteration",
				genfile::VariantEntry::Integer( iteration )
			) ;

			storage.store_per_variant_data(
				snps[*i].snp(),
				"result",
				result
			) ;

			std::map< std::string, genfile::VariantEntry > test_attributes = m_proximity_test->get_attributes( *i ) ;
			std::map< std::string, genfile::VariantEntry > const picker_attributes = m_snp_picker->get_attributes( *i ) ;
			for( std::size_t j = 0; j < output_columns.size(); ++j ) {
				std::string const column_name = genfile::string_utils::to_lower( output_columns[j] ) ;
				genfile::VariantEntry value ;
				if( column_name == "cm_from_start_of_chromosome" && recombination_offsets ) {
					value = recombination_offsets.get()[*i] ;
				}
				else {
					std::map< std::string, genfile::VariantEntry >::const_iterator where = test_attributes.find( column_name ) ;
					if( where == test_attributes.end() ) {
						where = picker_attributes.find( column_name ) ;
					}
					if( where == picker_attributes.end() ) {
						value = genfile::MissingValue() ;
					}
					else {
						value = where->second ;
					}
				}
				
				storage.store_per_variant_data(
					snps[ *i].snp(),
					output_columns[j],
					value
				) ;
			}
			
			if( genes.get() ) {
				genfile::Position const lower_bp = test_attributes[ "region_lower_bp" ].as< int >() ;
				genfile::Position const upper_bp = test_attributes[ "region_upper_bp" ].as< int >() ;
				std::vector< genes::Feature const* > const genes_in_region = genes->find_genes_in_region( snps[*i].snp().get_position().chromosome(), lower_bp, upper_bp ) ;
				if( genes_in_region.size() > 0 ) {
					std::set< std::string > geneNames ;
					std::size_t wNearest = genes_in_region.size() ;
					std::size_t nearestDistance = std::numeric_limits< std::size_t >::max() ;
					//std::cerr << "For SNP " << snps[*i] << ":\n" ;
					for( std::size_t gene_i = 0; gene_i < genes_in_region.size(); ++gene_i ) {
						// distance is 0 if SNP is in gene.
						// Otherwise it's the distance to the nearest end.
						// Note genes treated as closed.
						std::size_t distance = std::max(
							std::max( int( snps[*i].snp().get_position().position() ) - int( genes_in_region[gene_i]->end().position() ), 0 ),
							std::max( int( genes_in_region[gene_i]->start().position() ) - int( snps[*i].snp().get_position().position() ), 0 )
						) ;
						//std::cerr << " - " << genes_in_region[gene_i]->name() << " has distance " << distance << ".\n" ;
						if( distance < nearestDistance ) {
							wNearest = gene_i ;
							nearestDistance = distance ;
						}
						geneNames.insert( genes_in_region[gene_i]->name() ) ;
					}

					storage.store_per_variant_data(
						snps[*i].snp(),
						"nearest_gene_in_region",
						genes_in_region[ wNearest ]->name()
					) ;

					storage.store_per_variant_data(
						snps[*i].snp(),
						"distance_to_nearest_gene_in_region",
						genfile::VariantEntry::Integer( nearestDistance )
					) ;

					std::ostringstream theseGenes ;
					for( std::set< std::string >::const_iterator name_i = geneNames.begin(); name_i != geneNames.end(); ++name_i ) {
						theseGenes << (name_i == geneNames.begin() ? "" : "," ) << *name_i ;
					}
					storage.store_per_variant_data(
						snps[*i].snp(),
						"all_genes_in_region",
						theseGenes.str()
					) ;
				}
			}

			progress_context.notify_progress( snp_index+1, snps.size() ) ;
		}
	}


	std::string get_current_time_as_string() const {
		time_t rawtime ;
		struct tm * timeinfo ;
		char buffer[30] ;
		std::time( &rawtime ) ;
		timeinfo = std::localtime( &rawtime ) ;
		std::strftime( buffer, 80, "%Y-%m-%d %H:%M:%S", timeinfo ) ;
		return std::string( buffer ) ;
	}
		
	std::vector< std::string > get_output_columns() const {
		std::set< std::string > admissible_cols ;
		admissible_cols.insert( "snpid" ) ;
		admissible_cols.insert( "rsid" ) ;
		admissible_cols.insert( "chromosome" ) ;
		admissible_cols.insert( "position" ) ;
		admissible_cols.insert( "allele1" ) ;
		admissible_cols.insert( "allele2" ) ;
		admissible_cols.insert( "cm_from_start_of_chromosome" ) ;

		{
			std::set< std::string > const& dynamic_cols = m_proximity_test->get_attribute_names() ;
			for(
				std::set< std::string >::const_iterator i = dynamic_cols.begin();
				i != dynamic_cols.end();
				++i
			) {
				admissible_cols.insert( genfile::string_utils::to_lower( *i ) ) ;
			}
		}

		{
			std::set< std::string > const& dynamic_cols = m_snp_picker->get_attribute_names() ;
			for(
				std::set< std::string >::const_iterator i = dynamic_cols.begin();
				i != dynamic_cols.end();
				++i
			) {
				admissible_cols.insert( genfile::string_utils::to_lower( *i ) ) ;
			}
		}

		std::vector< std::string > output_cols = genfile::string_utils::split_and_strip( options().get< std::string >( "-output-cols" ), ",", " \t") ;
		if( output_cols.size() == 1 && output_cols[0] == "all" ) {
			// construct output cols.  Do it as follows to get them in a particular order.
			output_cols.clear() ;
			output_cols = genfile::string_utils::split_and_strip(
				"SNPID,rsid,chromosome,position,allele1,allele2,cM_from_start_of_chromosome,rank,region_lower_bp,region_upper_bp,region_lower_cM,region_upper_cM",
				","
			) ;
			for( std::vector< std::string >::iterator i = output_cols.begin(); i != output_cols.end(); ) {
				if( admissible_cols.find( genfile::string_utils::to_lower( *i )) == admissible_cols.end() ) {
					output_cols.erase( i ) ;
				} else {
					++i ;
				}
			}
		}

		for( std::size_t j = 0; j < output_cols.size(); ++j ) {
			if( admissible_cols.find( genfile::string_utils::to_lower( output_cols[j] )) == admissible_cols.end() ) {
				get_ui_context().logger() << "!! You have specified the output column \"" << output_cols[j] << "\", which I don't recognise.\n" ;
				get_ui_context().logger() << "   Available output columns are:\n" ;
				for(
					std::set< std::string >::const_iterator i = admissible_cols.begin();
					i != admissible_cols.end();
					++i
				) {
					get_ui_context().logger() << "   \"" << *i << "\".\n" ;
				}
				throw genfile::BadArgumentError(
					"InthinneratorApplication::get_output_columns()",
					"output column=" + genfile::string_utils::to_string( output_cols[j] ) + "\""
				) ;
			}
		}
		return output_cols ;
	}

} ;

int main( int argc, char** argv ) {
	try {
		InthinneratorApplication( argc, argv ) ;
	}
	catch( appcontext::HaltProgramWithReturnCode const& e ) {
		return e.return_code() ;
	}
	return 0 ;
}


	
