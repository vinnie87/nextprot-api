package org.nextprot.api.core.domain;

import java.io.Serializable;
import java.util.List;
import java.util.Set;

import org.nextprot.api.core.domain.annotation.Annotation;
import org.nextprot.api.core.service.fluent.EntryConfig;
import org.nextprot.api.core.utils.AnnotationUtils;
import org.nextprot.api.core.utils.ExperimentalContextUtil;
import org.nextprot.api.core.utils.PublicationUtils;
import org.nextprot.api.core.utils.XrefUtils;

public class EntryUtils implements Serializable{
	
	private static final long serialVersionUID = 3009334685615648172L;

	public static Entry filterEntryBySubPart(Entry entry, EntryConfig config) {
		
		
		List<Annotation> annotations;
		List<DbXref> xrefs;
		List<Publication> publications;
		List<ExperimentalContext> experimentalContexts;
		
		// Filter if necessary
		if (config.getSubpart() != null) {

			annotations = AnnotationUtils.filterAnnotationsByCategory(entry.getAnnotations(), config.getSubpart());
			entry.setAnnotations(annotations);
			
			if(!config.hasNoAdditionalReferences()){ //In case we don't care about xrefs, publications and experimental contexts (will be faster)
				Set<Long> xrefIds = AnnotationUtils.getXrefIdsForAnnotations(annotations);
				xrefIds.addAll(AnnotationUtils.getXrefIdsForInteractionsInteractants(annotations));
				xrefIds.addAll(AnnotationUtils.getXrefIdsForAnnotationsProperties(annotations));
				xrefs = XrefUtils.filterXrefsByIds(entry.getXrefs(), xrefIds);
				publications = PublicationUtils.filterPublicationsByIds(entry.getPublications(), AnnotationUtils.getPublicationIdsForAnnotations(annotations));

				experimentalContexts = ExperimentalContextUtil.filterExperimentalContextsByIds(entry.getExperimentalContexts(), AnnotationUtils.getExperimentalContextIdsForAnnotations(annotations));
				entry.setXrefs(xrefs);
				entry.setPublications(publications);
				entry.setExperimentalContexts(experimentalContexts);
			}

		}
		
		return entry;
	}

	
}
