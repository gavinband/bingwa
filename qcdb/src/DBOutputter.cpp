
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#include <string>
#include <memory>
#include <boost/optional.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/tuple/tuple.hpp>
#include <boost/thread/thread.hpp>
#include "genfile/SNPIdentifyingData.hpp"
#include "genfile/VariantEntry.hpp"
#include "genfile/Error.hpp"
#include "db/Connection.hpp"
#include "db/SQLStatement.hpp"
#include "appcontext/get_current_time_as_string.hpp"
#include "qcdb/DBOutputter.hpp"
#include "../../qctool_version_autogenerated.hpp"

namespace qcdb {
	DBOutputter::UniquePtr DBOutputter::create( std::string const& filename, std::string const& analysis_name, std::string const& analysis_description, Metadata const& metadata ) {
		return UniquePtr( new DBOutputter( filename, analysis_name, analysis_description, metadata ) ) ;
	}
	DBOutputter::SharedPtr DBOutputter::create_shared( std::string const& filename, std::string const& analysis_name, std::string const& analysis_description, Metadata const& metadata ) {
		return SharedPtr( new DBOutputter( filename, analysis_name, analysis_description, metadata ) ) ;
	}

	DBOutputter::~DBOutputter() {}

	DBOutputter::DBOutputter( std::string const& filename, std::string const& analysis_name, std::string const& analysis_description, Metadata const& metadata ):
		m_connection( db::Connection::create( filename )),
		m_analysis_name( analysis_name ),
		m_analysis_description( analysis_description ),
		m_metadata( metadata ),
		m_create_indices( true )
	{
		try {
			m_connection->run_statement( "PRAGMA journal_mode = OFF" ) ;
			m_connection->run_statement( "PRAGMA synchronous = OFF" ) ;
		}
		catch( db::Error const& ) {
			std::cerr << "qcdb::DBOutputter::DBOutputter(): unable to set PRAGMA synchronous=OFF, is another connection using this database?" ;
		}

		db::Connection::ScopedTransactionPtr transaction = m_connection->open_transaction( 1200 ) ;

		m_connection->run_statement(
			"CREATE TABLE IF NOT EXISTS Variant ( id INTEGER PRIMARY KEY, rsid TEXT, chromosome TEXT, position INTEGER, alleleA TEXT, alleleB TEXT )"
		) ;
		m_connection->run_statement(
			"CREATE INDEX IF NOT EXISTS Variant_position_index ON Variant( chromosome, position )"
		) ;
		m_connection->run_statement(
			"CREATE TABLE IF NOT EXISTS VariantIdentifier ( variant_id INTEGER NOT NULL, identifier TEXT, FOREIGN KEY( variant_id ) REFERENCES Variant( id ) ) "
		) ;
		m_connection->run_statement(
			"CREATE INDEX IF NOT EXISTS VariantIdentifierIdentifierIndex ON VariantIdentifier( identifier )"
		) ;
		m_connection->run_statement(
			"CREATE TABLE IF NOT EXISTS Entity ( "
				"id INTEGER PRIMARY KEY, name TEXT, description TEXT, "
				"UNIQUE( name, description ) "
			")"
		) ;
		m_connection->run_statement(
			"CREATE TABLE IF NOT EXISTS EntityData ( "
			"entity_id INTEGER NOT NULL, "
			"variable_id INTEGER NOT NULL, "
			"value TEXT, "
			"FOREIGN KEY (entity_id) REFERENCES Entity( id ), "
			"FOREIGN KEY (variable_id) REFERENCES Entity( id ) "
			")"
		) ;
		m_connection->run_statement(
			"CREATE TABLE IF NOT EXISTS EntityRelationship ( "
				"entity1_id INTEGER NOT NULL, "
				"relationship_id INTEGER NOT NULL, "
				"entity2_id INTEGER NOT NULL, "
				"FOREIGN KEY( entity1_id ) REFERENCES Entity( id ), "
				"FOREIGN KEY( relationship_id ) REFERENCES Entity( id ), "
				"FOREIGN KEY( entity2_id ) REFERENCES Entity( id ), "
				"UNIQUE( entity1_id, relationship_id, entity2_id ) "
			")"
		) ;

		m_connection->run_statement(
			"CREATE TABLE IF NOT EXISTS SummaryData ( "
			"variant_id INT, analysis_id INT, variable_id INT, value NONE, "
			"FOREIGN KEY( variant_id ) REFERENCES Variant( id ), "
			"FOREIGN KEY( analysis_id ) REFERENCES Entity( id ), "
			"FOREIGN KEY( variable_id ) REFERENCES Entity( id ) "
			")"
		) ;
		m_connection->run_statement(
			"CREATE INDEX IF NOT EXISTS EntityDataIndex ON EntityData( entity_id, variable_id )"
		) ;

		m_connection->run_statement(
			"CREATE VIEW IF NOT EXISTS EntityDataView AS "
			"SELECT ED.entity_id, E.name, ED.variable_id, V.name AS variable, ED.value "
			"FROM EntityData ED "
			"INNER JOIN Entity E "
			"ON E.id = ED.entity_id "
			"INNER JOIN Entity V "
			"ON V.id = ED.variable_id"
		) ;

		m_connection->run_statement(
			"CREATE VIEW IF NOT EXISTS EntityRelationshipView AS "
			"SELECT ER.entity1_id AS entity1_id, E1.name AS entity1, ER.relationship_id, R.name AS relationship, ER.entity2_id AS entity2_id, E2.name AS entity2 "
			"FROM EntityRelationship ER "
			"INNER JOIN Entity E1 "
			"  ON E1.id = ER.entity1_id " 
			"INNER JOIN Entity R "
			"  ON R.id = ER.relationship_id " 
			"INNER JOIN Entity E2 "
			"  ON E2.id = ER.entity2_id" 
		) ;

		m_connection->run_statement(
			"CREATE VIEW IF NOT EXISTS AnalysisView AS "
			"SELECT entity1_id AS id, A.name AS name, Tbl.value AS \"table\", Tool.value AS tool, A.description AS description "
			"FROM EntityRelationshipView ER "
			"INNER JOIN Entity A ON A.id = ER.entity1_id "
			"LEFT OUTER JOIN EntityDataView Tool ON Tool.entity_id = A.id AND Tool.variable == 'tool' "
			"LEFT OUTER JOIN EntityDataView Tbl ON Tbl.entity_id = ER.entity1_id AND Tbl.variable == 'table' "
			"WHERE ER.relationship == 'is_a' AND ER.entity2 == 'analysis'"
		) ;
		
		m_connection->run_statement(
			"CREATE VIEW IF NOT EXISTS VariantView AS "
			"SELECT          V.id AS id, V.rsid AS rsid, V.chromosome AS chromosome, V.position AS position, V.alleleA AS alleleA, V.alleleB AS alleleB, "
			"GROUP_CONCAT( VI.identifier ) AS alternate_identifier "
			"FROM Variant V "
			"LEFT OUTER JOIN VariantIdentifier VI "
			"  ON VI.variant_id = V.id "
			"GROUP BY V.id"
		) ;
		
		m_connection->run_statement(
			"CREATE VIEW IF NOT EXISTS SummaryDataView AS "
			"SELECT          V.id AS variant_id, V.chromosome, V.position, V.rsid, "
			"CASE WHEN length( V.alleleA < 10 ) THEN V.alleleA ELSE substr( V.alleleA, 1, 10 ) || '...' END AS alleleA, "
			"CASE WHEN length( V.alleleB < 10 ) THEN V.alleleB ELSE substr( V.alleleB, 1, 10 ) || '...' END AS alleleB, "
			"SD.analysis_id, Analysis.name AS analysis, Variable.id AS variable_id, Variable.name AS variable, "
			"SD.value AS value "
			"FROM SummaryData SD "
			"INNER JOIN Variant V ON V.id == SD.variant_id "
			"INNER JOIN Entity Analysis ON Analysis.id = SD.analysis_id "
			"INNER JOIN Entity Variable ON Variable.id = SD.variable_id "
			"WHERE analysis_id IN ( SELECT id FROM Entity )"
		) ;
		m_connection->run_statement(
			"CREATE TABLE IF NOT EXISTS AnalysisStatus ( "
				"analysis_id INTEGER NOT NULL REFERENCES Entity( id ), "
				"started TEXT NOT NULL, "
				"completed TEXT, "
				"status TEXT NOT NULL "
			")"
		) ;
		m_connection->run_statement(
			"CREATE VIEW IF NOT EXISTS AnalysisStatusView AS "
			"SELECT analysis_id, E.name AS analysis, started, completed, status "
			"FROM AnalysisStatus A "
			"INNER JOIN Entity E "
			"ON E.id == A.analysis_id"
		) ;
		construct_statements() ;
		store_metadata() ;
	}

	void DBOutputter::finalise( long options ) {
		if( options & eCreateIndices ) {
			db::Connection::ScopedTransactionPtr transaction = m_connection->open_transaction( 600 ) ;
			m_connection->run_statement(
				"CREATE INDEX IF NOT EXISTS Variant_rsid_index ON Variant( rsid )"
			) ;
			m_connection->run_statement(
				"CREATE INDEX IF NOT EXISTS VariantIdentifierVariantIndex ON VariantIdentifier( variant_id )"
			) ;
			m_connection->run_statement(
				"CREATE INDEX IF NOT EXISTS SummaryDataIndex ON SummaryData( variant_id, variable_id )"
			) ;
		}
		end_analysis( m_analysis_id ) ;
	}

	void DBOutputter::construct_statements() {
		m_insert_entity_statement = m_connection->get_statement( "INSERT INTO Entity ( name, description ) VALUES ( ?1, ?2 )" ) ;
		m_find_entity_data_statement = m_connection->get_statement( "SELECT * FROM EntityData WHERE entity_id == ?1 AND variable_id == ?2" ) ;
		m_find_entity_statement = m_connection->get_statement( "SELECT id FROM Entity E WHERE name == ?1 AND description == ?2" ) ;
		m_insert_entity_data_statement = m_connection->get_statement( "INSERT OR REPLACE INTO EntityData ( entity_id, variable_id, value ) VALUES ( ?1, ?2, ?3 )" ) ;
		m_insert_entity_relationship_statement = m_connection->get_statement( "INSERT OR REPLACE INTO EntityRelationship( entity1_id, relationship_id, entity2_id ) VALUES( ?1, ?2, ?3 )") ;
		m_find_variant_statement = connection().get_statement(
			"SELECT id, rsid FROM Variant WHERE chromosome == ?1 AND position == ?2 AND alleleA = ?3 AND alleleB = ?4"
		) ;
		m_insert_variant_statement = connection().get_statement(
			"INSERT INTO Variant ( rsid, chromosome, position, alleleA, alleleB ) "
			"VALUES( ?1, ?2, ?3, ?4, ?5 )"
		) ;
		m_find_variant_identifier_statement = m_connection->get_statement( "SELECT * FROM VariantIdentifier WHERE variant_id == ?1 AND identifier == ?2" ) ;
		m_insert_variant_identifier_statement = m_connection->get_statement( "INSERT INTO VariantIdentifier( variant_id, identifier ) VALUES ( ?1, ?2 )" ) ;
		m_insert_summarydata_statement = m_connection->get_statement(
			"INSERT INTO SummaryData ( variant_id, analysis_id, variable_id, value ) "
			"VALUES( ?1, ?2, ?3, ?4 )"
		) ;
	}

	void DBOutputter::store_metadata() {
		load_entities() ;
		m_is_a = get_or_create_entity_internal( "is_a", "is_a relationship" ) ;
		m_used_by = get_or_create_entity_internal( "used_by", "used_by relationship" ) ;

		try {
			m_analysis_id = create_entity_internal(
				m_analysis_name,
				m_analysis_description,
				get_or_create_entity_internal( "analysis", "class of analyses" )
			) ;
		} catch( db::StatementStepError const& e ) {
			throw genfile::BadArgumentError( "qcdb::DBOutputter::store_metadata()", "analysis_name=\"" + m_analysis_name + "\"", "An analysis with name \"" + m_analysis_name + "\" and description \"" + m_analysis_description + "\" already exists" ) ;
		} 

		start_analysis( m_analysis_id ) ;

		get_or_create_entity_data(
			m_analysis_id,
			get_or_create_entity( "tool", "Executable, pipeline, or script used to generate these results." ),
			"qctool revision " + std::string( globals::qctool_revision )
		) ;
		
		db::Connection::RowId const cmd_line_arg_id = get_or_create_entity_internal( "command-line argument", "Value supplied to a script or executable" ) ;
		for( Metadata::const_iterator i = m_metadata.begin(); i != m_metadata.end(); ++i ) {
			db::Connection::RowId key_id = get_or_create_entity( i->first, "command-line argument", cmd_line_arg_id ) ;
			get_or_create_entity_data(
				m_analysis_id,
				key_id,
				genfile::string_utils::join( i->second.first, "," ) + " (" + i->second.second + ")"
			) ;
		}
	}
	
	db::Connection::RowId DBOutputter::start_analysis( db::Connection::RowId const analysis_id ) const {
		db::Connection::StatementPtr stmnt = m_connection->get_statement( "INSERT INTO AnalysisStatus( analysis_id, started, status ) VALUES( ?, ?, ? )" ) ;
		stmnt->bind( 1, analysis_id ) ;
		stmnt->bind( 2, appcontext::get_current_time_as_string() ) ;
		stmnt->bind( 3, "incomplete" ) ;
		stmnt->step() ;
	}

	db::Connection::RowId DBOutputter::end_analysis( db::Connection::RowId const analysis_id ) const {
		db::Connection::StatementPtr stmnt = m_connection->get_statement( "UPDATE AnalysisStatus SET completed = ?, status = ? WHERE analysis_id == ?" ) ;
		stmnt->bind( 1, appcontext::get_current_time_as_string() ) ;
		stmnt->bind( 2, "successfully completed" ) ;
		stmnt->bind( 3, analysis_id ) ;
		stmnt->step() ;
	}
	
	void DBOutputter::load_entities() {
		db::Connection::StatementPtr stmnt = m_connection->get_statement( "SELECT id, name, description FROM Entity" ) ;
		stmnt->step() ;
		while( !stmnt->empty() ) {
			m_entity_map.insert( std::make_pair( std::make_pair( stmnt->get< std::string >( 1 ), stmnt->get< std::string >( 2 )), stmnt->get< db::Connection::RowId >( 0 ))) ;
			stmnt->step() ;
		}
		stmnt->reset() ;
	}

	db::Connection::RowId DBOutputter::get_or_create_entity_internal( std::string const& name, std::string const& description, boost::optional< db::Connection::RowId > class_id ) const {
		db::Connection::RowId result ;
		EntityMap::const_iterator where = m_entity_map.find( std::make_pair( name, description )) ;
		if( where != m_entity_map.end() ) {
			result = where->second ;
		}
		else {
			m_find_entity_statement
				->bind( 1, name )
				.bind( 2, description )
				.step() ;
			if( m_find_entity_statement->empty() ) {
				result = create_entity_internal( name, description, class_id ) ;
			} else {
				result = m_find_entity_statement->get< db::Connection::RowId >( 0 ) ;
			}
			m_find_entity_statement->reset() ;
		}
		return result ;
	}
	
	db::Connection::RowId DBOutputter::create_entity_internal( std::string const& name, std::string const& description, boost::optional< db::Connection::RowId > class_id ) const {
		db::Connection::RowId result ;
		m_insert_entity_statement
			->bind( 1, name )
			.bind( 2, description )
			.step() ;
			
		result = m_connection->get_last_insert_row_id() ;
		m_entity_map.insert( std::make_pair( std::make_pair( name, description ), result ) ) ;
		m_insert_entity_statement->reset() ;

		if( class_id ) {
			create_entity_relationship( result, m_is_a, *class_id ) ;
		}
		return result ;
	}

	db::Connection::RowId DBOutputter::get_or_create_entity( std::string const& name, std::string const& description, boost::optional< db::Connection::RowId > class_id ) const {
		db::Connection::RowId result ;
		EntityMap::const_iterator where = m_entity_map.find( std::make_pair( name, description )) ;
		if( where != m_entity_map.end() ) {
			result = where->second ;
		}
		else {
			m_find_entity_statement
				->bind( 1, name )
				.bind( 2, description )
				.step() ;
			if( m_find_entity_statement->empty() ) {
				result = create_entity_internal( name, description, class_id ) ;
				create_entity_relationship( result, m_used_by, m_analysis_id ) ;
			} else {
				result = m_find_entity_statement->get< db::Connection::RowId >( 0 ) ;
			}
			m_find_entity_statement->reset() ;
		}
		return result ;
	}

	void DBOutputter::create_entity_relationship( db::Connection::RowId entity1_id, db::Connection::RowId relationship_id, db::Connection::RowId entity2_id ) const {
		m_insert_entity_relationship_statement
			->bind( 1, entity1_id )
			.bind( 2, relationship_id )
			.bind( 3, entity2_id )
			.step() ;
			
		m_insert_entity_relationship_statement->reset() ;
	}

	db::Connection::RowId DBOutputter::get_or_create_entity_data( db::Connection::RowId const entity_id, db::Connection::RowId const variable_id, genfile::VariantEntry const& value ) const {
		db::Connection::RowId result ;

		m_find_entity_data_statement
			->bind( 1, entity_id )
			.bind( 2, variable_id ).step() ;

		if( m_find_entity_data_statement->empty() ) {
			m_insert_entity_data_statement
				->bind( 1, entity_id )
				.bind( 2, variable_id )
				.bind( 3, value )
				.step() ;
			result = m_connection->get_last_insert_row_id() ;
			m_insert_entity_data_statement->reset() ;
		} else {
			result = m_find_entity_data_statement->get< db::Connection::RowId >( 0 ) ;
		}
		m_find_entity_data_statement->reset() ;
		return result ;
	}

	void DBOutputter::add_alternative_variant_identifier( db::Connection::RowId const variant_id, std::string const& identifier, std::string const& rsid ) const {
		if( identifier != rsid ) {
			add_variant_identifier( variant_id, identifier ) ;
		}
	}

	void DBOutputter::add_variant_identifier( db::Connection::RowId const variant_id, std::string const& identifier ) const {
		m_find_variant_identifier_statement
			->bind( 1, variant_id )
			.bind( 2, identifier )
			.step() ;
		if( m_find_variant_identifier_statement->empty() ) {
			m_insert_variant_identifier_statement
				->bind( 1, variant_id )
				.bind( 2, identifier )
				.step() ;
			m_insert_variant_identifier_statement->reset() ;
		}
		m_find_variant_identifier_statement->reset() ;
	}

	db::Connection::RowId DBOutputter::get_or_create_variant( genfile::SNPIdentifyingData2 const& snp ) const {
		db::Connection::RowId result ;
		m_find_variant_statement
			->bind( 1, std::string( snp.get_position().chromosome() ) )
			.bind( 2, snp.get_position().position() )
			.bind( 3, snp.get_first_allele() )
			.bind( 4, snp.get_second_allele() )
			.step()
		;
		if( m_find_variant_statement->empty() ) {
			m_insert_variant_statement
				->bind( 1, snp.get_rsid() )
				.bind( 2, std::string( snp.get_position().chromosome() ) )
				.bind( 3, snp.get_position().position() )
				.bind( 4, snp.get_first_allele())
				.bind( 5, snp.get_second_allele())
				.step()
			;

			result = connection().get_last_insert_row_id() ;
			m_insert_variant_statement->reset() ;
			snp.get_alternative_identifiers(
				boost::bind(
					&DBOutputter::add_alternative_variant_identifier,
					this,
					result,
					_1,
					snp.get_rsid()
				)
			) ;
		} else {
			result = m_find_variant_statement->get< db::Connection::RowId >( 0 ) ;
			std::string const rsid = m_find_variant_statement->get< std::string >( 1 ) ;
			add_alternative_variant_identifier( result, snp.get_rsid(), rsid ) ;
			snp.get_alternative_identifiers(
				boost::bind(
					&DBOutputter::add_alternative_variant_identifier,
					this,
					result,
					_1,
					rsid
				)
			) ;
		}
		m_find_variant_statement->reset() ;
		return result ;
	}
	
	void DBOutputter::insert_summary_data(
		db::Connection::RowId snp_id,
		db::Connection::RowId variable_id,
		genfile::VariantEntry const& value
	) const {
		m_insert_summarydata_statement
			->bind( 1, snp_id )
			.bind( 2, m_analysis_id )
			.bind( 3, variable_id )
			.bind( 4, value )
			.step() ;
		m_insert_summarydata_statement->reset() ;
	}
}
