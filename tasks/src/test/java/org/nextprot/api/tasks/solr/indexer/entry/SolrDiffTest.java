package org.nextprot.api.tasks.solr.indexer.entry;

import java.util.ArrayList;
import java.util.List;

import org.apache.solr.client.solrj.SolrQuery;
import org.apache.solr.client.solrj.SolrServerException;
import org.apache.solr.client.solrj.impl.HttpSolrServer;
import org.apache.solr.client.solrj.response.QueryResponse;
import org.apache.solr.common.SolrDocumentList;
import org.nextprot.api.commons.service.MasterIdentifierService;
import org.nextprot.api.core.domain.Entry;
import org.nextprot.api.core.service.EntryBuilderService;
import org.nextprot.api.core.service.fluent.EntryConfig;
import org.nextprot.api.solr.index.EntryIndex.Fields;
import org.springframework.beans.factory.InitializingBean;
import org.springframework.beans.factory.annotation.Autowired;


public class SolrDiffTest extends SolrBuildIntegrationTest implements InitializingBean {
	
	@Autowired	private EntryBuilderService entryBuilderService = null;
	@Autowired	private MasterIdentifierService masterIdentifierService = null;

	HttpSolrServer solr = new HttpSolrServer("http://kant.isb-sib.ch:8983/solr/npentries1");

	protected Object getValueForFieldInCurrentSolrImplementation(String entryName, Fields field) {

		List<Object> result = new ArrayList<Object>();

		try {
			SolrQuery query = new SolrQuery();
			query.setQuery("id:" + entryName);
			query.setFields(field.getName());
			QueryResponse response;
			response = solr.query(query);
			SolrDocumentList solrResults = response.getResults();

			for (int i = 0; i < solrResults.size(); ++i) {
				result.add(solrResults.get(i).get(field.getName()));
			}
		} catch (SolrServerException e) {
			e.printStackTrace();
			throw new RuntimeException(e);
		}
		
		return result.get(0);
	}
	

	private List<String> entries = null;
	
	protected Entry getEntry(String entryName){
		return entryBuilderService.build(EntryConfig.newConfig(entryName).withEverything());
	}
	
	protected Entry getEntry(EntryConfig config){
		return entryBuilderService.build(config);
	}
	
	protected Entry getEntry(int i){
		return entryBuilderService.build(EntryConfig.newConfig(entries.get(i)).withEverything());
	}

	@Override
	public void afterPropertiesSet() throws Exception {
		entries = new ArrayList<String>(masterIdentifierService.findUniqueNames());
	}

	
	public int getEntriesCount() {
		return entries.size();
	}

	/**
	 * Utility to parse html-like lines generated by current Solr implementation
	 * It removes all other elements, example:
	 * 
	 * 
	 * <b>anno_name : </b>mature protein <b>anno_qualname : </b>GOLD
	 * <b>description : </b>Putative uncharacterized protein HSD52
	 * <b>an_synonyms : </b>PRO_0000342680 <b>anno_name : </b>uniprot keyword
	 * <b>anno_qualname : </b>GOLD <b>description : </b>Complete proteome
	 * <b>cv_ancestors_acs : </b>KW-0181 <b>cv_ac : </b>KW-0181 <b>cv_name :
	 * </b>Complete proteome
	 */

	public static String getValueFromRawData(String html, String subfield) {
		String aux = "";
		aux = html.replaceAll("<p>|</p>", "");
		String[] btags = aux.split("<b>");
		for (String bt : btags) {
			if (bt.startsWith(subfield)) {
				return bt.substring(bt.indexOf("</b>") + 4);
			}
		}
		return null;
	}

}
