package org.nextprot.api.tasks.solr.indexer.entry.impl;

import java.util.Arrays;
import java.util.Collection;

import org.nextprot.api.core.domain.Entry;
import org.nextprot.api.core.domain.EntryProperties;
import org.nextprot.api.solr.index.EntryIndex.Fields;
import org.nextprot.api.tasks.solr.indexer.entry.EntryFieldBuilder;
import org.nextprot.api.tasks.solr.indexer.entry.FieldBuilder;


@EntryFieldBuilder
public class FilterAndPropertiesFieldsBuilder extends FieldBuilder{
	
	@Override
	protected void init(Entry entry){
		
		// Filters and entry properties
		EntryProperties props = entry.getProperties();
		addField(Fields.ISOFORM_NUM, props.getIsoformCount());
		int cnt;
		cnt = props.getPtmCount();
		if(cnt > 0) {
			addField(Fields.PTM_NUM, cnt);
		}
		cnt = props.getVarCount();
		if(cnt > 0) {
			addField(Fields.VAR_NUM, cnt);
		}
		String filters = "";
		if(props.getFilterstructure()) filters += "filterstructure ";
		if(props.getFilterdisease()) filters += "filterdisease ";
		if(props.getFilterexpressionprofile()) filters += "filterexpressionprofile ";
		if(props.getFiltermutagenesis()) filters += "filtermutagenesis ";
		if(props.getFilterproteomics()) filters += "filterproteomics ";
		if(filters.length() > 0) {
			addField(Fields.FILTERS, filters.trim());
		}
		addField(Fields.AA_LENGTH, props.getMaxSeqLen()); // max length among all isoforms
				
	}

	@Override
	public Collection<Fields> getSupportedFields() {
		return Arrays.asList(Fields.ISOFORM_NUM, Fields.PTM_NUM, Fields.VAR_NUM, Fields.FILTERS, Fields.AA_LENGTH);
	}


}