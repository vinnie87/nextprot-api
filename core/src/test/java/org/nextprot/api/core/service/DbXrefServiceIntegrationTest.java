package org.nextprot.api.core.service;

import org.junit.Assert;
import org.junit.Test;
import org.nextprot.api.commons.constants.AnnotationCategory;
import org.nextprot.api.commons.service.MasterIdentifierService;
import org.nextprot.api.core.domain.DbXref;
import org.nextprot.api.core.domain.annotation.Annotation;
import org.nextprot.api.core.domain.annotation.AnnotationEvidence;
import org.nextprot.api.core.domain.annotation.AnnotationIsoformSpecificity;
import org.nextprot.api.core.test.base.CoreUnitBaseTest;
import org.springframework.beans.factory.annotation.Autowired;
import org.springframework.test.context.ActiveProfiles;

import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.PrintWriter;
import java.net.*;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

import static org.junit.Assert.assertTrue;

@ActiveProfiles({ "dev" })
public class DbXrefServiceIntegrationTest extends CoreUnitBaseTest {

	@Autowired private DbXrefService xrefService;
	@Autowired private MasterIdentifierService masterIdentifierService;
/*
 * This query finds entries having a single xref among 'Orphanet', 'KEGGPathway' , 'Reactome' and 'DrugBank'
 * It is convenient for tests: we know we get a single annotation from xrefs for a given entry
 * Example:
 * NX_A0AVF1 for Reactome
 * NX_A1L167 for Kegg
 * NX_A0PJY2 for Orphanet
 * NX_Q9Y2D1 for DrugBank

select a.unique_name, string_agg(a.acs, ',') as acs, string_agg(a.cv_name, ',') as dbs, count(*) as dbcount, sum(a.cnt) as xrcount from (
select si.unique_name, db.cv_name, count(*) as cnt, string_agg(x.accession, ',') as acs
from sequence_identifiers si
inner join identifier_resource_assoc ira on (si.identifier_id=ira.identifier_id)
inner join db_xrefs x on (ira.resource_id=x.resource_id)
inner join cv_databases db on (x.cv_database_id=db.cv_id)
where si.cv_type_id=1 and si.cv_status_id=1 
and db.cv_name in ('Orphanet', 'DrugBank','KEGGPathway','Reactome')
group by si.unique_name, db.cv_name
) a
group by a.unique_name
having sum(a.cnt)=1
;


 */
	
	@Test
	public void shouldReturn_1_ReactomeXrefAsAnnotation() {
		List<Annotation> annotations = this.xrefService.findDbXrefsAsAnnotationsByEntry("NX_A0AVF1");
		assertTrue(annotations.size() == 1);
		Annotation annot = annotations.get(0);
		assertTrue(annot.getCategory().equals(AnnotationCategory.PATHWAY.getDbAnnotationTypeName()));
		assertTrue(annot.getAPICategory()== AnnotationCategory.PATHWAY);
		assertTrue(annot.getQualityQualifier().equals("GOLD"));
		Assert.assertEquals("Intraflagellar transport", annot.getDescription());
		for (AnnotationIsoformSpecificity spec: annot.getTargetingIsoformsMap().values()) {
			assertTrue(spec.getSpecificity().equals("UNKNOWN"));
		}
		assertTrue(annot.getEvidences().size()==1);
		AnnotationEvidence evi = annot.getEvidences().get(0);
		assertTrue(evi.getAssignedBy().equals("Uniprot"));
		assertTrue(evi.getEvidenceCodeAC().equals("ECO:0000305"));
		assertTrue(evi.getResourceAccession().equals("R-HSA-5620924"));
		assertTrue(evi.getResourceDb().equals("Reactome"));

		Assert.assertTrue(annotations.get(0).getProperties().isEmpty());
	}
	
	@Test
	public void shouldReturn_1_KEGGPathwayXrefAsAnnotation() {
		List<Annotation> annotations = this.xrefService.findDbXrefsAsAnnotationsByEntry("NX_A1L167");
		assertTrue(annotations.size() == 1);
		Annotation annot = annotations.get(0);
		assertTrue(annot.getCategory().equals(AnnotationCategory.PATHWAY.getDbAnnotationTypeName()));
		assertTrue(annot.getAPICategory()== AnnotationCategory.PATHWAY);
		assertTrue(annot.getQualityQualifier().equals("GOLD"));
		Assert.assertEquals("Ubiquitin mediated proteolysis", annot.getDescription());
		for (AnnotationIsoformSpecificity spec: annot.getTargetingIsoformsMap().values()) {
			assertTrue(spec.getSpecificity().equals("UNKNOWN"));
		}
		assertTrue(annot.getEvidences().size()==1);
		AnnotationEvidence evi = annot.getEvidences().get(0);
		assertTrue(evi.getAssignedBy().equals("NextProt"));
		assertTrue(evi.getEvidenceCodeAC().equals("ECO:0000305"));
		assertTrue(evi.getResourceAccession().equals("hsa04120+134111"));
		assertTrue(evi.getResourceDb().equals("KEGGPathway"));

		Assert.assertTrue(annotations.get(0).getProperties().isEmpty());
	}
	
	@Test
	public void shouldReturn_1_OrphanetXrefAsAnnotation() {
		List<Annotation> annotations = this.xrefService.findDbXrefsAsAnnotationsByEntry("NX_A0PJY2");
		assertTrue(annotations.size() == 1);
		Annotation annot = annotations.get(0);
		assertTrue(annot.getCategory().equals(AnnotationCategory.DISEASE.getDbAnnotationTypeName()));
		assertTrue(annot.getAPICategory()== AnnotationCategory.DISEASE);
		assertTrue(annot.getQualityQualifier().equals("GOLD"));
		Assert.assertEquals("Kallmann syndrome", annot.getDescription());
		for (AnnotationIsoformSpecificity spec: annot.getTargetingIsoformsMap().values()) {
			assertTrue(spec.getSpecificity().equals("UNKNOWN"));
		}
		assertTrue(annot.getEvidences().size()==1);
		AnnotationEvidence evi = annot.getEvidences().get(0);
		assertTrue(evi.getAssignedBy().equals("Uniprot"));
		assertTrue(evi.getEvidenceCodeAC().equals("ECO:0000305"));
		assertTrue(evi.getResourceAccession().equals("478"));
		assertTrue(evi.getResourceDb().equals("Orphanet"));

		Assert.assertTrue(annotations.get(0).getProperties().isEmpty());
	}
	
	@Test
	public void shouldReturn_1_DrugBankXrefAsAnnotation() {
		List<Annotation> annotations = this.xrefService.findDbXrefsAsAnnotationsByEntry("NX_Q9Y2D1");
		assertTrue(annotations.size() == 1);
		Annotation annot = annotations.get(0);
		assertTrue(annot.getCategory().equals(AnnotationCategory.SMALL_MOLECULE_INTERACTION.getDbAnnotationTypeName()));
		assertTrue(annot.getAPICategory()== AnnotationCategory.SMALL_MOLECULE_INTERACTION);
		assertTrue(annot.getQualityQualifier().equals("GOLD"));
		Assert.assertEquals("Pseudoephedrine", annot.getDescription());
		for (AnnotationIsoformSpecificity spec: annot.getTargetingIsoformsMap().values()) {
			assertTrue(spec.getSpecificity().equals("UNKNOWN"));
		}
		assertTrue(annot.getEvidences().size()==1);
		AnnotationEvidence evi = annot.getEvidences().get(0);
		assertTrue(evi.getAssignedBy().equals("Uniprot"));
		assertTrue(evi.getEvidenceCodeAC().equals("ECO:0000305"));
		assertTrue(evi.getResourceAccession().equals("DB00852"));
		assertTrue(evi.getResourceDb().equals("DrugBank"));

		Assert.assertTrue(annotations.get(0).getProperties().isEmpty());
	}

	@Test
	public void reactomeXrefShouldHaveEmptyProperties() {

		assertEmptyProperties("NX_A0AVF1", 42610527);
	}

	@Test
	public void KEGGPathwayXrefShouldHaveEmptyProperties() {

		assertEmptyProperties("NX_A1L167", 14559832);
	}

	@Test
	public void orphanetXrefShouldHaveEmptyProperties() {

		assertEmptyProperties("NX_A0PJY2", 1077769);
	}

	@Test
	public void drugBankXrefShouldHaveEmptyProperties() {

		assertEmptyProperties("NX_Q9Y2D1", 983678);
	}

    @Test
    public void testPercentSignSTypeLinkHasUrlCorrectlyResolved() {

        List<DbXref> xrefs = this.xrefService.findDbXrefsByMaster("NX_P01308");

        Assert.assertEquals(1120, xrefs.size());

        for (DbXref xref : xrefs) {

            if (xref.getDbXrefId() == 1272250) {

                Assert.assertEquals("http://www.ncbi.nlm.nih.gov/protein/%s", xref.getLinkUrl());
                Assert.assertEquals("http://www.ncbi.nlm.nih.gov/protein/NP_000198.1", xref.getResolvedUrl());

                break;
            }
        }
    }

	@Test
	public void testPercentSignUTypeLinkHasUrlCorrectlyResolved() {

		List<DbXref> xrefs = this.xrefService.findDbXrefsByMaster("NX_P01308");

		Assert.assertEquals(1120, xrefs.size());

		for (DbXref xref : xrefs) {

			if (xref.getDbXrefId() == 16387756) {

				Assert.assertEquals("http://pbil.univ-lyon1.fr/cgi-bin/acnuc-ac2tree?query=%u&db=HOGENOM", xref.getLinkUrl());
				Assert.assertEquals("http://pbil.univ-lyon1.fr/cgi-bin/acnuc-ac2tree?query=P01308&db=HOGENOM", xref.getResolvedUrl());

                break;
			}
		}
	}

	@Test
	public void testBrendaTypeLinkHasUrlCorrectlyResolved() {

		List<DbXref> xrefs = this.xrefService.findDbXrefsByMaster("NX_Q9BXA6");

		Assert.assertEquals(520, xrefs.size());

		for (DbXref xref : xrefs) {

			if (xref.getDbXrefId() == 964246) {

				Assert.assertEquals("http://www.brenda-enzymes.org/enzyme.php?ecno=%s", xref.getLinkUrl());
				Assert.assertEquals("http://www.brenda-enzymes.org/enzyme.php?ecno=2.7.11.1", xref.getResolvedUrl());

                break;
			}
		}
	}

	//@Test
	public void testAllEntriesDbXrefs() {

		Set<String> allEntryNames = masterIdentifierService.findUniqueNames();

		for (String entryName : allEntryNames) {

			List<DbXref> xrefs = this.xrefService.findDbXrefsByMaster(entryName);

			for (DbXref xref : xrefs) {

				Assert.assertTrue(!xref.getAccession().isEmpty());
				Assert.assertTrue(!xref.getUrl().isEmpty());
				Assert.assertTrue(!xref.getLinkUrl().isEmpty());
				Assert.assertTrue(!xref.getResolvedUrl().isEmpty());
			}
		}
	}

	//@Test
	public void logAllEntriesXrefUrlStatus() throws FileNotFoundException {

		Set<String> visitedLinkedURLs = new HashSet<>();

		PrintWriter pw = new PrintWriter("/tmp/allentries-xrefs-url.tsv");

		Set<String> allEntryAcs = masterIdentifierService.findUniqueNames();

		pw.write("entry ac\tdb\txref ac\turl\thttp status\tresolved url\thttp status\n");

		int count=0;
		for (String entryAc : allEntryAcs) {

			List<DbXref> xrefs = this.xrefService.findDbXrefsByMaster(entryAc);

			for (DbXref xref : xrefs) {

				if ((count % 50) == 0) {
					pw.flush();
				}

				String resolvedUrl = xref.getResolvedUrl();

				String linkedURL = xref.getLinkUrl();

				if (!visitedLinkedURLs.contains(linkedURL)) {

					Response response = requestUrls(xref);

					int j=0;
					int tries = 3;
					while (response.getResolvedUrlHttpStatus().equals("TIMEOUT") && j<tries) {

						response = requestUrls(xref);
						j++;
					}

					String db = xref.getDatabaseName();
					String xrefAc = xref.getAccession();
					String url = xref.getUrl();

					pw.write(entryAc);
					pw.write("\t");
					pw.write(db);
					pw.write("\t");
					pw.write(xrefAc);
					pw.write("\t");
					pw.write(url);
					pw.write("\t");
					pw.write(response.getUrlHttpStatus());
					pw.write("\t");
					pw.write(resolvedUrl);
					pw.write("\t");
					pw.write(response.getResolvedUrlHttpStatus());
					pw.write("\n");

					visitedLinkedURLs.add(linkedURL);

					count++;
				}
			}
		}

		pw.flush();
		pw.close();
	}

	private Response requestUrls(DbXref xref) {

		String url = xref.getUrl();
		String urlHttpStatus = getResponseCode(url);
		String resolvedUrlHttpStatus = getResponseCode(xref.getResolvedUrl());

		return new Response(urlHttpStatus, resolvedUrlHttpStatus);
	}

	private static class Response {

		String urlHttpStatus;
		String resolvedUrlHttpStatus;

		public Response(String urlHttpStatus, String resolvedUrlHttpStatus) {
			this.urlHttpStatus = urlHttpStatus;
			this.resolvedUrlHttpStatus = resolvedUrlHttpStatus;
		}

		public String getUrlHttpStatus() {
			return urlHttpStatus;
		}

		public String getResolvedUrlHttpStatus() {
			return resolvedUrlHttpStatus;
		}

	}

	public String getResponseCode(String url)  {

		String response;
		HttpURLConnection con = null;

		try {
			URL obj = new URL(url);
			con = (HttpURLConnection) obj.openConnection();

			con.setRequestMethod("HEAD");
			con.setRequestProperty("User-Agent", "Mozilla/5.0");
			con.setConnectTimeout(5000);
			con.connect();

			System.out.println("Http HEAD request "+url);
			response = String.valueOf(con.getResponseCode());

		} catch (SocketTimeoutException e) {
			System.err.println(e.getMessage());
			response = "TIMEOUT";
		} catch (ProtocolException e) {
			System.err.println(e.getMessage());
			response = "PROTOCOL";
		} catch (MalformedURLException e) {
			System.err.println(e.getMessage());
			response = "MALFORMEDURL";
		} catch (IOException e) {
			System.err.println(e.getMessage());
			response = "IO";
		}

		if (con != null)
			con.disconnect();

		return response;
	}

	private void assertEmptyProperties(String entryName, long propertyId) {

		List<DbXref> dbxrefs = this.xrefService.findDbXrefsByMaster(entryName);

		for (DbXref xref : dbxrefs)
			if (xref.getDbXrefId() == propertyId)
				Assert.assertTrue(xref.getProperties().isEmpty());
	}
	
}
