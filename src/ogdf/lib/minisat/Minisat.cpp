// disable VC++ warnings when using fopen
#define _CRT_SECURE_NO_WARNINGS

#include <ogdf/basic/ArrayBuffer.h>
#include <ogdf/lib/minisat/Minisat.h>

namespace Minisat
{

void Clause::addMultiple(int Amount, ...)
{
	va_list params;
	va_start(params, Amount);
	for (int i = 0; i < Amount; ++i) {
		Internal::Var paramVar = va_arg(params, Internal::Var);
		Internal::Lit l;
		if (paramVar >= 0) {
			l = Internal::mkLit(paramVar-1, true);
		}
		else {
			l = Internal::mkLit(-(paramVar + 1), false);
		}
		m_ps.push(l);
	}
	va_end(params);
}

Clause *Formula::newClause()
{
	m_Clauses.push_back(new Clause);
	return (m_Clauses.back());

}

void Formula::finalizeClause(const clause cl)
{
	for (int i = 0; i < cl->m_ps.size(); ++i) {
		// if an variable does not exist, it will be generated (and all between the gap)
		if (!(Internal::var(cl->m_ps[i]) < Solver::nVars())) {
			int max = Solver::nVars();
			for (int j = 0; j < Internal::var(cl->m_ps[i]) + 1 - max; ++j) {
				Solver::newVar();
			}
		}
	}
	Solver::addClause(cl->m_ps);
}

bool Formula::finalizeNotExtensibleClause(const clause cl)
{
	//proves if variables from clause are valid (still known to formula)
	for (int i = 0; i < cl->m_ps.size(); ++i) {
		if (!(Internal::var(cl->m_ps[i]) < Solver::nVars())) {
			m_messages << "Variable " << i << " is not present.";
			return false;
		}
	}
	Solver::addClause(cl->m_ps);
	return true;
}

Clause *Formula::getClause( const int pos )
{
	if ( pos < (int)m_Clauses.size() )
		return m_Clauses[pos];
	else
		return nullptr;
}

bool Formula::solve( Model &ReturnModel )
{
	bool solv = Solver::solve();

	if ( solv )
		ReturnModel.setModel ( *this );

	return solv;
}


bool Formula::solve( Model &ReturnModel, double& timeLimit )
{
    SolverStatus st;
    bool solv = Solver::solve(timeLimit, st);

    if ( solv )
       ReturnModel.setModel ( *this );

    ReturnModel.solverStatus = st;

    return solv;
}


void Formula::removeClause(int i)
{
	Internal::CRef cr = Solver::clauses[i];
	Solver::removeClause(cr);
	int j, k;
	for (j = k = 0; j < Solver::clauses.size(); ++j) {
		if (j == !i) {
			clauses[k++] = clauses[j];
		}
	}
	clauses.shrink( j - k );
	delete &m_Clauses[i];
	m_Clauses.erase( m_Clauses.begin() + i );
}

void Formula::reset()
{
	free();
	Solver::assigns.clear();
	Solver::vardata.clear();
	Solver::activity.clear();
	Solver::seen.clear();
	Solver::polarity.clear();
	Solver::decision.clear();
	Solver::trail.clear();
	Solver::dec_vars = 0;
	Solver::model.clear();
}

void Formula::free()
{
	for (auto i = 0; i < Solver::clauses.size(); ++i) {
		Solver::removeClause(Solver::clauses[i]);
	}
	for (auto &clause : m_Clauses) {
		delete clause;
	}
	Solver::clauses.shrink(Solver::clauses.size());
	m_Clauses.clear();
}

bool Formula::readDimacs(const char *filename)
{
	std::ifstream is(filename);
	if (!is.is_open()) return false;
	return readDimacs(is);
}

bool Formula::readDimacs(const string &filename)
{
	std::ifstream is(filename);
	if (!is.is_open()) return false;
	return readDimacs(is);
}

bool Formula::readDimacs(std::istream &in)
{
	std::string currentString;

	while (!in.eof()) {
		in >> currentString;
		if (currentString == "p") {
			in >> currentString;
			if (currentString == "cnf") {
				break;
			}
		}
	}
	if (in.eof()) {
		return false;
	}

	int numVars = -1;
	int numClauses = -1;
	in >> numVars >> numClauses;
	if (numVars < 0 || numClauses < 0) {
		return false;
	}
	newVars(numVars);

	int clauseCount = 0;
	ogdf::ArrayBuffer<int> literals;
	for (int lit; in >> lit;) {
		if (lit) {
			if (lit > numVars
			 || -lit > numVars) {
				ogdf::Logger::slout() << "Literal does not represent a valid variable (index too high)" << endl;
				return false;
			}
			literals.push(lit);
		} else {
			addClause(literals);
			literals.clear();
			clauseCount++;
		}
	}
	if (!literals.empty()) {
		ogdf::Logger::slout(ogdf::Logger::LL_MINOR) << "Last clause is not terminated by 0 marker, but we accept it nonetheless" << endl;
		addClause(literals);
	}
	if(clauseCount != numClauses) {
		ogdf::Logger::slout(ogdf::Logger::LL_MINOR) << "Number of clauses differs from file header" << endl;
	}

	return true;
}

bool Formula::writeDimacs(const char *filename)
{
	std::ofstream os(filename);
	if (!os.is_open()) return false;
	return writeDimacs(os);
}

bool Formula::writeDimacs(const string &filename)
{
	std::ofstream os(filename);
	if (!os.is_open()) return false;
	return writeDimacs(os);
}

bool Formula::writeDimacs(std::ostream &f)
{
	f << "p cnf " << getVariableCount() << " " << getClauseCount() << std::endl;
	for (auto &clause : m_Clauses) {
		for (int j = 0; j < clause->m_ps.size(); ++j) {
#if defined(OGDF_DEBUG)
			std::cout
			  << "Sign : "
			  << Internal::sign(clause->m_ps[j])
			  << "Var : "
			  << Internal::var(clause->m_ps[j]) + 1
			  << std::endl;
#endif
			f << " "
			  << Clause::convertLitSign(clause->m_ps[j])
			  << Internal::var(clause->m_ps[j]) + 1;
		}
		f << " 0" << std::endl;
	}
	return true;
}

} // MinisatAdvanced
