package org.nextprot.api.core.dao;

import static org.junit.Assert.assertTrue;

import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;

import org.junit.Test;
import org.nextprot.api.commons.service.MasterIdentifierService;
import org.nextprot.api.core.dao.impl.PeptideMappingDaoImpl;
import org.nextprot.api.core.domain.Family;
import org.nextprot.api.core.test.base.CoreUnitBaseTest;
import org.springframework.beans.factory.annotation.Autowired;
import org.springframework.test.context.ActiveProfiles;

@ActiveProfiles({ "dev" })
public class PeptideMappingDaoIntegrationTest extends CoreUnitBaseTest {

	@Autowired private PeptideMappingDao peptideMappingDao;
	@Autowired private MasterIdentifierService masterIdentifierService;

	
	@Test
	public void shouldReturn_Some_Expected_Mappings() {
		String entryName = "NX_Q9UGM3";
		Long masterId = this.masterIdentifierService.findIdByUniqueName(entryName);
 		List<Map<String,Object>> mappings = peptideMappingDao.findNaturalPeptideMappingAnnotationsByMasterId(masterId);
		Map<String,Integer> annotPepIsoRecCount = new TreeMap<>();
 		
 		for (Map<String,Object> mapping: mappings) {
 			Long annotationId = (Long) mapping.get(PeptideMappingDao.KEY_ANNOTATION_ID);
 			String qual = (String) mapping.get(PeptideMappingDao.KEY_QUALITY_QUALIFIER);
 			String iso = (String) mapping.get(PeptideMappingDao.KEY_ISO_UNIQUE_NAME);
 			String pep = (String) mapping.get(PeptideMappingDao.KEY_PEP_UNIQUE_NAME);
 			Integer p1 = (Integer) mapping.get(PeptideMappingDao.KEY_FIRST_POS);
 			Integer p2 = (Integer) mapping.get(PeptideMappingDao.KEY_LAST_POS);
 			Integer rank = (Integer) mapping.get(PeptideMappingDao.KEY_RANK);
 			
 			// basic check of the data returned by the DAO 
 			
 			assertTrue(annotationId>0);
 			assertTrue(qual.equals("GOLD") || qual.equals("SILVER"));
 			assertTrue(iso.startsWith(entryName));
 			assertTrue(pep.startsWith("NX_"));
 			assertTrue(p1>0);
 			assertTrue(p2>p1);
 			assertTrue(rank>0);
 			
 			// check that we NEVER have more than 1 mapping record for each peptide/isoform/annotation key
 			
			String key = "" + annotationId + "/" + pep + "/"+ iso;
			if (!annotPepIsoRecCount.containsKey(key)) annotPepIsoRecCount.put(key, new Integer(0));
			annotPepIsoRecCount.put(key, new Integer(annotPepIsoRecCount.get(key) + 1));
 		}
 		assertTrue(annotPepIsoRecCount.size()>0);
		for (Map.Entry<String, Integer> data: annotPepIsoRecCount.entrySet()) {
			//System.out.println(data.getKey() + " => " + data.getValue());
			assertTrue(data.getValue()==1);
		}
 	}	

}