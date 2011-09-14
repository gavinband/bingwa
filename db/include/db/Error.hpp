#ifndef DB_SQL_ERROR_HPP
#define DB_SQL_ERROR_HPP

namespace db {
	struct Error: public std::exception
	{
		Error( std::string const& caller, std::string const& db_spec, int error ) ;
		~Error() throw() ;
		char const* what() const throw() { return "db::Error" ; }

		int const& error() const { return m_error ; }
		std::string const& spec() const { return m_spec ; }
		std::string description() const ;
		
	private:
		std::string const m_caller ;
		std::string const m_spec ;
		int m_error ;
	} ;
	
	struct ConnectionError: public Error
	{
		ConnectionError( std::string const& caller, std::string const& db_spec, int error_code ): Error( caller, db_spec, error_code ) {}
		char const* what() const throw() { return "db::ConnectionError" ; }
	} ;

	struct StatementPreparationError: public Error
	{
		StatementPreparationError( std::string const& caller, std::string const& db_spec, int error_code, std::string SQL ): Error( caller, db_spec, error_code ), m_SQL( SQL ) {}
		~StatementPreparationError() throw() {}
		char const* what() const throw() { return "db::StatementPreparationError" ; }
		std::string const& SQL() const { return m_SQL ; }
	private:
		std::string const m_SQL ;
	} ;

	struct StatementStepError: public Error
	{
		StatementStepError( std::string const& caller, std::string const& db_spec, int error_code ): Error( caller, db_spec, error_code ) {}
		char const* what() const throw() { return "db::StatementStepError" ; }
	} ;

	struct ValueBindError: public Error
	{
		~ValueBindError() throw() {}
		ValueBindError(
			std::string const& caller,
			std::string const& db_spec,
			int error_code,
			std::string const& slot
		): Error( caller, db_spec, error_code ), m_slot( slot ) {}
		char const* what() const throw() { return "db::StatementStepError" ; }
		std::string const& slot() const { return m_slot ; }
	private:
		std::string const m_slot ;
	} ;
	
}

#endif
