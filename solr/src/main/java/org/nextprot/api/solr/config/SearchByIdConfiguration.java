package org.nextprot.api.solr.config;


import org.nextprot.api.solr.Query;

public class SearchByIdConfiguration extends IndexConfiguration {

	public static final String ID_SEARCH = "id";
	public static final String PL_SEARCH = "pl_search";

	public SearchByIdConfiguration(String name) {
		super(name);
	}
	
	@Override
	public String buildQuery(Query query) {
		return query.getQueryString(true);
	}
}
