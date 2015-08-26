package org.nextprot.api.web.service.impl;

import org.junit.Assert;
import org.junit.Test;
import org.nextprot.api.commons.service.MasterIdentifierService;
import org.nextprot.api.core.domain.DbXref;
import org.nextprot.api.core.domain.Entry;
import org.nextprot.api.core.service.EntryBuilderService;
import org.nextprot.api.core.service.fluent.EntryConfig;
import org.nextprot.api.web.dbunit.base.mvc.WebIntegrationBaseTest;
import org.springframework.beans.factory.annotation.Autowired;

public class BuildEntryTest extends WebIntegrationBaseTest {

	@Autowired
	private MasterIdentifierService masterIdentifierService;

	@Autowired
	private EntryBuilderService entryBuilderService;

	@Test
	public void testWithEnsemblGeneShouldBePresent() throws Exception {

		Entry entry = entryBuilderService.build(EntryConfig.newConfig("NX_P01308").withGenomicMappings().withChromosomalLocations());

		Assert.assertEquals(1, entry.getChromosomalLocations().size());
		Assert.assertEquals("ENSG00000254647", entry.getChromosomalLocations().get(0).getAccession());
		Assert.assertTrue(!entry.getGenomicMappings().isEmpty());
		Assert.assertEquals("ENSG00000254647", entry.getGenomicMappings().get(0).getAccession());
		Assert.assertEquals("Ensembl", entry.getGenomicMappings().get(0).getDatabase());
	}

	@Test
	public void testVirtualGeneShouldBeAbsent() throws Exception {

		Entry entry = entryBuilderService.build(EntryConfig.newConfig("NX_Q6ZTC4").withGenomicMappings().withChromosomalLocations().withXrefs());

		Assert.assertEquals(1, entry.getChromosomalLocations().size());
		Assert.assertTrue(entry.getChromosomalLocations().get(0).getAccession().isEmpty());
		Assert.assertTrue(entry.getGenomicMappings().isEmpty());
		for (DbXref xref : entry.getXrefs()) {
			Assert.assertTrue(!xref.getAccession().matches("NX_VG.+"));
		}
	}

	@Test
	public void testVirtualGene2ShouldBeAbsent() throws Exception {

		Entry entry = entryBuilderService.build(EntryConfig.newConfig("NX_O00370").withGenomicMappings().withChromosomalLocations().withXrefs());

		Assert.assertEquals(1, entry.getChromosomalLocations().size());
		Assert.assertTrue(entry.getChromosomalLocations().get(0).getAccession().isEmpty());
		Assert.assertTrue(entry.getGenomicMappings().isEmpty());
		for (DbXref xref : entry.getXrefs()) {
			Assert.assertTrue(!xref.getAccession().matches("VG.+"));
		}
	}
}