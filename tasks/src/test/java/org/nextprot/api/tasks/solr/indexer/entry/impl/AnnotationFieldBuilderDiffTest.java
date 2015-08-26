package org.nextprot.api.tasks.solr.indexer.entry.impl;

import static org.junit.Assert.assertEquals;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Iterator;
import java.util.List;
import java.util.Set;

import org.junit.Test;
import org.nextprot.api.commons.service.MasterIdentifierService;
import org.nextprot.api.core.domain.Entry;
import org.nextprot.api.core.service.EntryBuilderService;
import org.nextprot.api.core.service.fluent.EntryConfig;
import org.nextprot.api.solr.index.EntryIndex.Fields;
import org.nextprot.api.tasks.solr.indexer.entry.SolrDiffTest;
import org.springframework.beans.factory.annotation.Autowired;

public class AnnotationFieldBuilderDiffTest extends SolrDiffTest {

	@Autowired
	private EntryBuilderService entryBuilderService = null;
	@Autowired
	private MasterIdentifierService masterIdentifierService = null;

	@Test
	public void testAnnotationsAndFunctionalDescriptions() {

		Set<String> entries = masterIdentifierService.findUniqueNames();
		Iterator<String> entriesIt = entries.iterator();

		int i = 0;
		while (entriesIt.hasNext() && i < 10) {

			String entryName = entriesIt.next();
			System.out.print(i++);
			System.out.println(entryName);

			Entry entry = entryBuilderService.build(EntryConfig.newConfig(entryName).withAnnotations());

			testAnnotations(entryName, entry);
			testFunctionalDesc(entryName, entry);
		}

	}

	@SuppressWarnings("unchecked")
	public void testFunctionalDesc(String entryName, Entry entry) {

		Fields field = Fields.FUNCTION_DESC;

		AnnotationFieldBuilder afb = new AnnotationFieldBuilder(entry);
		List<String> functionalDescriptions = afb.build(entry, field, List.class);

		List<String> expectedValues = (List<String>) getValueForFieldInCurrentSolrImplementation(entryName, field);

		if (!((expectedValues == null) && (functionalDescriptions == null))) {
			assertEquals(functionalDescriptions.size(), expectedValues.size());
			if (!functionalDescriptions.isEmpty()) {
				assertEquals(functionalDescriptions.get(0), expectedValues.get(0));
			}
		}

	}

	@SuppressWarnings("unchecked")
	public void testAnnotations(String entryName, Entry entry) {

		Fields field = Fields.ANNOTATIONS;

		AnnotationFieldBuilder afb = new AnnotationFieldBuilder(entry);
		List<String> annotations = afb.build(entry, field, List.class);
		List<String> expectedRawValues = (List<String>) getValueForFieldInCurrentSolrImplementation(entryName, field);
		List<String> expectedValues = new ArrayList<String>();

		for (String s : expectedRawValues) {
			String aux = getDescriptionFromRawData(s);
			expectedValues.add(aux);
		}

		Collections.sort(annotations);
		Collections.sort(expectedValues);

		assertEquals(annotations.size(), expectedValues.size());

		// TODO remove "reference proteome"
		for (int i = 0; i < annotations.size(); i++) {
			assertEquals(annotations.get(i), expectedValues.get(i));
		}

	}

	@Test
	public void testCleanRawData() {
		String result = getDescriptionFromRawData("<p><b>anno_name : </b>caution</p><p><b>anno_qualname : </b>GOLD</p><p><b>description : </b>Product of a dubious CDS prediction.</p>");
		assertEquals("Product of a dubious CDS prediction.", result);
	}

	/**
	 * This removes all other elements, example:
	 * 
	 * 
	 * <b>anno_name : </b>mature protein <b>anno_qualname : </b>GOLD
	 * <b>description : </b>Putative uncharacterized protein HSD52
	 * <b>an_synonyms : </b>PRO_0000342680 <b>anno_name : </b>uniprot keyword
	 * <b>anno_qualname : </b>GOLD <b>description : </b>Complete proteome
	 * <b>cv_ancestors_acs : </b>KW-0181 <b>cv_ac : </b>KW-0181 <b>cv_name :
	 * </b>Complete proteome
	 */
	static String getDescriptionFromRawData(String html) {
		String aux = "";
		aux = html.replaceAll("<p>|</p>", "");
		String[] btags = aux.split("<b>");
		for (String bt : btags) {
			if (bt.startsWith("description")) {
				return bt.substring(bt.indexOf("</b>") + 4);
			}
		}
		return null;
	}

}
