#include <iostream>
#include <cassert>
#include <string>
#include <exception>

#include "sqlite3/sqlite3.h"
#include "db/SQLite3Connection.hpp"
#include "db/SQLStatement.hpp"
#include "db/SQLite3Statement.hpp"
#include "db/SQLite3Error.hpp"

namespace db {
	SQLite3Statement::SQLite3Statement( SQLite3Connection* connection, std::string const& SQL ):
		m_connection( connection ),
		m_have_results( false )
	{
		assert( m_connection ) ;
		m_statement = m_connection->prepare_sql( SQL ) ;
		//std::cerr << "SQLite3Statement::SQLite3Statement(): statement is \"" + SQL + "\".\n" ;
		assert( m_statement != 0 ) ;
	}
	
	SQLite3Statement::~SQLite3Statement() {
		m_connection->finalise_statement( m_statement ) ;
	}

	bool SQLite3Statement::step() {
		m_have_results = m_connection->step_statement( m_statement ) ;
		return m_have_results ;
	}

	SQLite3Statement::operator void*() const {
		if( m_have_results ) {
			return reinterpret_cast< void* >( m_statement ) ;
		}
		return 0 ;
	}

	std::size_t SQLite3Statement::get_number_of_columns() const {
		return std::size_t( get_column_count() ) ;
	}

	std::string SQLite3Statement::get_name_of_column( std::size_t i ) const {
		assert( m_statement != 0 ) ;
		assert( i < get_number_of_columns() ) ;
		return std::string( sqlite3_column_origin_name( m_statement, i )) ;
	}

	void SQLite3Statement::bind( std::size_t i, int value ) const {
		assert( m_statement != 0 ) ;
		int error = sqlite3_bind_int( m_statement, i, value ) ;
		if( error != SQLITE_OK ) {
			throw SQLite3Error( "SQLite3Statement::bind()", error ) ;
		}
	}

	void SQLite3Statement::bind( std::size_t i, std::string const& value ) const {
		assert( m_statement != 0 ) ;
		int error = sqlite3_bind_text( m_statement, i, value.c_str(), value.size(), SQLITE_TRANSIENT ) ;
		if( error != SQLITE_OK ) {
			throw SQLite3Error( "SQLite3Statement::bind()", error ) ;
		}
	}

	void SQLite3Statement::bind( std::size_t i, char const* buffer, std::size_t n ) const {
		assert( m_statement != 0 ) ;
		int error = sqlite3_bind_blob( m_statement, i, reinterpret_cast< void const* >( buffer ), int( n ), SQLITE_TRANSIENT ) ;
		if( error != SQLITE_OK ) {
			throw SQLite3Error( "SQLite3Statement::bind()", error ) ;
		}
	}

	void SQLite3Statement::reset() const {
		assert( m_statement != 0 ) ;
		int error = sqlite3_reset( m_statement ) ;
		if( error != SQLITE_OK ) {
			throw SQLite3Error( "SQLite3Statement::reset_to_start()", error ) ;
		}
	}

	std::string SQLite3Statement::get_sql() const {
		assert( m_statement != 0 ) ;
		return std::string( sqlite3_sql( m_statement )) ;
	}

	int SQLite3Statement::get_column_int( int column_id ) const {
		assert( m_statement != 0 ) ;
		return sqlite3_column_int( m_statement, column_id ) ;
	}

	double SQLite3Statement::get_column_double( int column_id ) const {
		assert( m_statement != 0 ) ;
		return sqlite3_column_double( m_statement, column_id ) ;
	}
	
	std::string SQLite3Statement::get_column_string( int column_id ) const {
		assert( m_statement != 0 ) ;
		return reinterpret_cast< char const * >( sqlite3_column_text( m_statement, column_id ) ) ;
	}

	char SQLite3Statement::get_column_char( int column_id ) const {
		assert( m_statement != 0 ) ;
		char const* p = reinterpret_cast< char const * >( sqlite3_column_text( m_statement, column_id )) ;
		int bytes = sqlite3_column_bytes( m_statement, column_id ) ;
		assert( bytes == 1 ) ;
		return *p ;
	}
	
	int SQLite3Statement::get_column_count() const {
		assert( m_statement != 0 ) ;
		return sqlite3_column_count( m_statement ) ;
	}
}