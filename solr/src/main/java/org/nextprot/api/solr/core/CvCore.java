package org.nextprot.api.solr.core;

import org.apache.solr.client.solrj.SolrQuery.ORDER;
import org.nextprot.api.commons.utils.Pair;
import org.nextprot.api.solr.query.config.AutocompleteConfiguration;
import org.nextprot.api.solr.query.config.FieldConfigSet;
import org.nextprot.api.solr.query.config.IndexConfiguration;
import org.nextprot.api.solr.query.config.IndexParameter;
import org.nextprot.api.solr.query.config.SortConfig;
import org.springframework.beans.factory.annotation.Value;
import org.springframework.stereotype.Component;

@Component
public class CvCore extends CoreTemplate {

	private static final String NAME = "npcvs1";
	private static final String SIMPLE = "simple";

	@Value("${solr.url}")
	private String solrServerUrl;

	public CvCore() {
		super(NAME, Entity.Term);
	}

	@Override
	protected IndexConfiguration newDefaultConfiguration() {

		IndexConfiguration defaultConfig = IndexConfiguration.SIMPLE();

		defaultConfig.addConfigSet(new FieldConfigSet(IndexParameter.FL)
				.add(CvSolrField.AC)
				.add(CvSolrField.NAME)
				.add(CvSolrField.SYNONYMS)
				.add(CvSolrField.DESCRIPTION)
				.add(CvSolrField.PROPERTIES)
				.add(CvSolrField.FILTERS));
		//.add(CvField.TEXT));

		defaultConfig.addConfigSet(new FieldConfigSet(IndexParameter.QF)
				.addWithBoostFactor(CvSolrField.AC, 64)
				.addWithBoostFactor(CvSolrField.NAME, 32)
				.addWithBoostFactor(CvSolrField.SYNONYMS, 32)
				.addWithBoostFactor(CvSolrField.DESCRIPTION, 16)
				.addWithBoostFactor(CvSolrField.PROPERTIES, 8)
				.addWithBoostFactor(CvSolrField.OTHER_XREFS, 8));

		defaultConfig.addConfigSet(new FieldConfigSet(IndexParameter.PF)
				.addWithBoostFactor(CvSolrField.AC, 640)
				.addWithBoostFactor(CvSolrField.NAME, 320)
				.addWithBoostFactor(CvSolrField.SYNONYMS, 320)
				.addWithBoostFactor(CvSolrField.DESCRIPTION, 160)
				.addWithBoostFactor(CvSolrField.PROPERTIES, 80)
				.addWithBoostFactor(CvSolrField.OTHER_XREFS, 80));

		defaultConfig.addOtherParameter("defType", "edismax")
				.addOtherParameter("facet", "true")
				.addOtherParameter("facet.field", "filters")
				.addOtherParameter("facet.limit", "10")
				.addOtherParameter("facet.method", "enum")
				.addOtherParameter("facet.mincount", "1")
				.addOtherParameter("facet.sort", "count");

		defaultConfig.addSortConfig(
				SortConfig.create("default", new Pair[] {
						Pair.create(CvSolrField.SCORE, ORDER.desc),
						Pair.create(CvSolrField.FILTERS, ORDER.asc)
				}),
				SortConfig.create("name", new Pair[] {
						Pair.create(CvSolrField.NAME_S, ORDER.asc),
						Pair.create(CvSolrField.FILTERS, ORDER.asc)
				})
		);
		defaultConfig.setDefaultSortName("default");

		return defaultConfig;
	}

	@Override
	protected AutocompleteConfiguration newAutoCompleteConfiguration(IndexConfiguration configuration) {
		AutocompleteConfiguration autocompleteConfig = new AutocompleteConfiguration(configuration);

		autocompleteConfig.addOtherParameter("defType", "edismax")
				.addOtherParameter("facet", "true")
				.addOtherParameter("facet.field", "text")
				.addOtherParameter("facet.limit", "10")
				.addOtherParameter("facet.method", "enum")
				.addOtherParameter("facet.mincount", "1")
				.addOtherParameter("facet.sort", "count")
				.addOtherParameter("stopwords", "true");

		return autocompleteConfig;
	}

	@Override
	protected SortConfig[] newSortConfigurations() {
		return new SortConfig[0];
	}

	@Override
	protected void setupConfigurations() {

		addConfiguration(defaultConfiguration);
		setConfigAsDefault(SIMPLE);
		addConfiguration(autocompleteConfiguration);
	}

	@Override
	protected String getServerUrl() {
		return solrServerUrl;
	}

	@Override
	public SolrField[] getSchema() {
		return CvSolrField.values();
	}
}
