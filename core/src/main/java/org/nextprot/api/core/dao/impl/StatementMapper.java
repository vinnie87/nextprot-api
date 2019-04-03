package org.nextprot.api.core.dao.impl;

import org.nextprot.commons.statements.Statement;
import org.nextprot.commons.statements.StatementBuilder;
import org.nextprot.commons.statements.specs.NXFlatTableSchema;
import org.nextprot.commons.statements.specs.StatementField;
import org.springframework.jdbc.core.RowMapper;

import java.sql.ResultSet;
import java.sql.SQLException;

public class StatementMapper implements RowMapper<Statement> {
	public Statement mapRow(ResultSet rs, int rowNum) throws SQLException {

		NXFlatTableSchema schema = NXFlatTableSchema.fromResultSet(rs);

		StatementBuilder sfbuilder = new StatementBuilder()
				.withSpecifications(schema);

		for(StatementField key : schema.getFields()) {

			String value = rs.getString(key.getName());
			if (value != null) {
				sfbuilder.addField(key, value);
			}
		}
		return sfbuilder.build();
	}
	
}