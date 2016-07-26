package com.nextprot.api.annotation.builder.statement.dao;

import java.util.List;

import org.nextprot.commons.statements.Statement;
import org.nextprot.commons.statements.constants.AnnotationType;

public interface StatementDao {

	List<Statement> findNormalStatements(AnnotationType annotationType, String entryName);

	List<Statement> findProteoformStatements(AnnotationType annotationType, String entryName);

	List<Statement> findStatementsByAnnotIsoIds(AnnotationType annotationType, List<String> ids);

	List<Statement> findStatementsByAnnotEntryId(AnnotationType annotationType, String annotEntryId);

}
