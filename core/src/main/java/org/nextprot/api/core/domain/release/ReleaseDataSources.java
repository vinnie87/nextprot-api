package org.nextprot.api.core.domain.release;

import java.util.HashSet;
import java.util.Set;

/**
 * Enum representing the datasources to display in the stats page
 * 
 * This enum was generated from the Grails code using http://jsfiddle.net/ddtxra/a2qmfx22/6/
 * 
 * @author Daniel Teixeira http://github.com/ddtxra
 *
 */
public enum ReleaseDataSources {
        
        UniProtKBSwissProt("UniProtKB/Swiss-Prot", "UniProtKB", "UniProt knowledgebase Swiss-Prot section", "http://www.uniprot.org/"),
        BGee("BGee", "BGee", "Database for gene expression evolution", "http://bgee.unil.ch"),
        COSMIC("COSMIC", "Cosmic", "Catalogue of somatic mutations in cancer", "http://www.sanger.ac.uk/genetics/CGP/cosmic/"),
        Ensembl("Ensembl", "Ensembl genome", "Ensembl human genome browser", "http://www.ensembl.org/"),
        ENZYME("ENZYME", "UniProtKB enzyme classification", "Enzyme nomenclature database", "http://www.expasy.org/enzyme/"),
        GO("GO", "GO", "Gene Ontology", "http://www.geneontology.org/"),
        HPA("HPA", "HPA", "Human Protein Atlas", "http://www.proteinatlas.org/"),
        IntAct("IntAct", "IntAct interactions", "Molecular interaction database", "http://www.ebi.ac.uk/intact/"),
        InterPro("InterPro", "InterPro", "Integrated resource of protein families, domains and functional sites", "http://www.ebi.ac.uk/interpro/"),
        MeSH("MeSH", "MeSH", "Medical Subject Headings", "http://www.nlm.nih.gov/mesh/MBrowser.html"),
        PeptideAtlas("PeptideAtlas", "Human PeptideAtlas", "Peptides identified by mass spectrometry", "http://www.peptideatlas.org/"),
        PROSITE("PROSITE", "PROSITEDOC", "Protein domain and family database", "http://www.expasy.org/prosite/"),
        PubMed("PubMed", "PubMed publications", "Citations for biomedical literature from MEDLINE, life science journals, and online books.", "http://www.ncbi.nlm.nih.gov/sites/entrez?db=pubmed"),
        SRMAtlas("SRMAtlas", "Human SRMAtlas", "Targeted proteomics assays", "http://www.srmatlas.org/"),
        UniProtGOA("UniProt-GOA", "GO annotations", "Gene Ontology annotations", "http://www.ebi.ac.uk/GOA");
         
        private String displayName, cvName, description, url; 

		ReleaseDataSources(String displayName, String cvName, String description, String url){
        	this.displayName = displayName;
        	this.cvName = cvName;
        	this.description = description;
        	this.url = url;
        }
		
        public String getDisplayName() {
			return displayName;
		}

		public void setDisplayName(String displayName) {
			this.displayName = displayName;
		}

		public String getCvName() {
			return cvName;
		}

		public void setCvName(String cvName) {
			this.cvName = cvName;
		}

		public String getDescription() {
			return description;
		}

		public void setDescription(String description) {
			this.description = description;
		}

		public String getUrl() {
			return url;
		}

		public void setUrl(String url) {
			this.url = url;
		}

		public static Set<String> getDistinctCvNames() {
			Set<String> cvNames = new HashSet<String>();
			for(ReleaseDataSources ds : ReleaseDataSources.values()){
				cvNames.add(ds.getCvName());
			}
			return cvNames;
		}

		public static ReleaseDataSources cvValueOf(String cvName) {
			for(ReleaseDataSources ds : ReleaseDataSources.values()){
				if(cvName.equals(ds.getCvName()))
					return ds;
			} return null;
		}

		
}

    