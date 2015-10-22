package org.nextprot.api.core.dao.impl;

import org.nextprot.api.commons.constants.AnnotationCategory;
import org.nextprot.api.commons.spring.jdbc.DataSourceServiceLocator;
import org.nextprot.api.commons.utils.SQLDictionary;
import org.nextprot.api.core.dao.AntibodyMappingDao;
import org.nextprot.api.core.domain.annotation.Annotation;
import org.nextprot.api.core.domain.annotation.AnnotationEvidence;
import org.nextprot.api.core.domain.annotation.AnnotationIsoformSpecificity;
import org.springframework.beans.factory.annotation.Autowired;
import org.springframework.jdbc.core.RowMapper;
import org.springframework.jdbc.core.namedparam.MapSqlParameterSource;
import org.springframework.jdbc.core.namedparam.NamedParameterJdbcTemplate;
import org.springframework.jdbc.core.namedparam.SqlParameterSource;
import org.springframework.stereotype.Repository;

import java.sql.ResultSet;
import java.sql.SQLException;
import java.util.*;

@Repository
public class AntibodyMappingDaoImpl implements AntibodyMappingDao {

	@Autowired private SQLDictionary sqlDictionary;
	@Autowired private DataSourceServiceLocator dsLocator;

	@Override
	public List<Annotation> findAntibodyMappingAnnotationsById(long masterId) {
		
		SqlParameterSource namedParams = new MapSqlParameterSource("id", masterId);

		AntibodyRowMapper rowMapper =  new AntibodyRowMapper(); 
		new NamedParameterJdbcTemplate(dsLocator.getDataSource()).query(sqlDictionary.getSQLQuery("antibodies-by-id"), namedParams, rowMapper);
		return rowMapper.getAnnotations();
	}
	
	private static class AntibodyRowMapper implements RowMapper<Annotation>{
		
		Map<Long, Annotation> annotationsMap = new HashMap<>();
		
		@Override
		public Annotation mapRow(ResultSet resultSet, int row) throws SQLException {

			Long annotationId = resultSet.getLong("annotation_id");

			if (!annotationsMap.containsKey(annotationId)){

                Annotation annotation = new Annotation();
				annotation.setAnnotationId(annotationId);
				annotation.setCategory(AnnotationCategory.ANTIBODY_MAPPING);
                annotation.setQualityQualifier("GOLD"); // TODO: IS THIS KIND OF INFO ACCESSIBLE ?
				annotation.setTargetingIsoforms(new ArrayList<AnnotationIsoformSpecificity>());

                // From template antibody-list.xml.vm:
                //  <evidence is-negative='false' evidence-code='ECO:0000154' resource-assoc-type='evidence' resource-internal-ref='$xref.dbXrefId'
                //   resource-type='database' source-internal-ref='$antibody.assignedBy' #set($UNESCAPED_CDATA = $xref.resolvedUrl) />
                AnnotationEvidence evidence = new AnnotationEvidence();
                evidence.setAnnotationId(annotationId);
                evidence.setNegativeEvidence(false);
                evidence.setEvidenceCodeAC("ECO:0000154");
                evidence.setResourceAssociationType("evidence");
                evidence.setResourceType("database");
                evidence.setResourceId(resultSet.getLong("db_xref_id"));
                evidence.setAssignedBy(resultSet.getString("antibody_src"));

                annotation.setEvidences(Collections.singletonList(evidence));

                annotationsMap.put(annotationId, annotation);
			}

			String isoName = resultSet.getString("iso_unique_name");
			Annotation annotation = annotationsMap.get(annotationId);

			AnnotationIsoformSpecificity isoSpecificity = new AnnotationIsoformSpecificity();

            isoSpecificity.setIsoformName(isoName);
			isoSpecificity.setFirstPosition((Integer)resultSet.getObject("first_pos")); // better than using getInt which returns a primitive and can not be null
			isoSpecificity.setLastPosition((Integer)resultSet.getObject("last_pos"));
            isoSpecificity.setSpecificity("SPECIFIC");

			annotation.getTargetingIsoformsMap().put(isoName, isoSpecificity);

			return annotation;
		}
		
		public List<Annotation> getAnnotations() {
			return new ArrayList<>(annotationsMap.values());
		}
	}
}
