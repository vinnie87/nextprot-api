package org.nextprot.api.core.dao.impl;

import java.sql.ResultSet;
import java.sql.SQLException;

import org.nextprot.api.commons.spring.jdbc.DataSourceServiceLocator;
import org.nextprot.api.commons.utils.SQLDictionary;
import org.nextprot.api.core.dao.EntryPropertiesDao;
import org.nextprot.api.core.domain.EntryProperties;
import org.springframework.beans.factory.annotation.Autowired;
import org.springframework.jdbc.core.namedparam.MapSqlParameterSource;
import org.springframework.jdbc.core.namedparam.NamedParameterJdbcTemplate;
import org.springframework.jdbc.core.namedparam.SqlParameterSource;
import org.springframework.jdbc.core.simple.ParameterizedRowMapper;
import org.springframework.stereotype.Repository;

@Repository
public class EntryPropertiesDaoImpl implements EntryPropertiesDao {

	@Autowired	private SQLDictionary sqlDictionary;
	@Autowired	private DataSourceServiceLocator dsLocator;
	
	@Override
	public EntryProperties findEntryProperties(String uniqueName) {
		SqlParameterSource namedParameters = new MapSqlParameterSource("uniqueName", uniqueName);
		return new NamedParameterJdbcTemplate(dsLocator.getDataSource()).queryForObject(sqlDictionary.getSQLQuery("entry-properties"), namedParameters, new EntryPropertyRowMapper());

	}
	
	static class EntryPropertyRowMapper implements ParameterizedRowMapper<EntryProperties> {

		@Override
		public EntryProperties mapRow(ResultSet resultSet, int row) throws SQLException {
			int res;
			EntryProperties properties = new EntryProperties();
			properties.setProteinExistence(resultSet.getString("pe"));
			res = resultSet.getInt("ptmcnt");
			properties.setPtmCount(res == -1? 0:res);
			properties.setIsoformCount(resultSet.getInt("isocnt"));
			res = resultSet.getInt("varcnt");
			properties.setVarCount(res == -1? 0:res);
			res = resultSet.getInt("mutcnt");
			properties.setMutagenesisCount(res == -1? 0:res);
			res = resultSet.getInt("intcnt");
			properties.setInteractionCount(res == -1? 0:res);
			properties.setMaxSeqLen(resultSet.getInt("maxlen"));
			properties.setFilterstructure(resultSet.getInt("structure") != -1);
			properties.setFilterdisease(resultSet.getInt("disease") != -1);
			properties.setFilterproteomics(resultSet.getInt("proteomic") != -1);
			properties.setFilterexpressionprofile(resultSet.getInt("expression") != -1);
			return properties;
		}
		
	}

}