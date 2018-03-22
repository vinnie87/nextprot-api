package org.nextprot.api.core.domain;

import java.io.Serializable;
import java.util.ArrayList;
import java.util.Collection;
import java.util.List;

public class GenomicMapping implements Serializable {

	private static final long serialVersionUID = 1L;

	private long geneSeqId;
	private String database;
	private String accession;
	private List<IsoformGeneMapping> isoformGeneMappings;

	public GenomicMapping() {

		isoformGeneMappings = new ArrayList<>();
	}

	public String getDatabase() {
		return database;
	}

	public void setDatabase(String database) {
		this.database = database;
	}

	public String getAccession() {
		return accession;
	}

	public void setAccession(String accession) {
		this.accession = accession;
	}

	public long getGeneSeqId() {
		return geneSeqId;
	}

	public void setGeneSeqId(long geneSeqId) {
		this.geneSeqId = geneSeqId;
	}

	public List<IsoformGeneMapping> getIsoformGeneMappings() {
		return isoformGeneMappings;
	}

	public boolean addAllIsoformGeneMappings(Collection<IsoformGeneMapping> mappings) {

		return isoformGeneMappings.addAll(mappings);
	}
}
