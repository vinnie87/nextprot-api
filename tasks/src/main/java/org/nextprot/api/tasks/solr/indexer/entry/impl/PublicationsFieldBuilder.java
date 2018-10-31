package org.nextprot.api.tasks.solr.indexer.entry.impl;

import org.apache.log4j.Logger;
import org.nextprot.api.core.domain.Entry;
import org.nextprot.api.core.domain.Publication;
import org.nextprot.api.core.domain.PublicationAuthor;
import org.nextprot.api.core.domain.publication.GlobalPublicationStatistics;
import org.nextprot.api.core.domain.publication.JournalResourceLocator;
import org.nextprot.api.core.service.PublicationService;
import org.nextprot.api.solr.index.EntryField;
import org.nextprot.api.tasks.solr.indexer.entry.EntryFieldBuilder;
import org.springframework.beans.factory.annotation.Autowired;
import org.springframework.stereotype.Service;

import java.util.Arrays;
import java.util.Collection;
import java.util.List;
import java.util.SortedSet;

@Service
public class PublicationsFieldBuilder extends EntryFieldBuilder {
	
	protected Logger logger = Logger.getLogger(PublicationsFieldBuilder.class);

	@Autowired
	private PublicationService publicationService;

	@Override
	public void collect(Entry entry, boolean gold) {

		// Publications
		// Shouldn't Xrefs to PubMed and DOIs be also indexed here ?
		List<Publication> publications = entry.getPublications();
		int publi_computed_count = 0;
		int publi_curated_count = 0;
		int publi_large_scale_count = 0;
		String Jinfo = "";

		for (Publication currpubli : publications) {

			long pubId = currpubli.getPublicationId();
			logger.debug("looking for stats about pair " + entry.getUniqueName() + " - pubId:" + pubId);
            GlobalPublicationStatistics.PublicationStatistics publiStats = publicationService.getPublicationStatistics(pubId);
            if(publiStats.isComputed()) publi_computed_count++;
			if(publiStats.isCurated()) publi_curated_count++;
			if(publiStats.isLargeScale()) publi_large_scale_count++;

			if(currpubli.isLocatedInScientificJournal()) {

				JournalResourceLocator journalLocator = currpubli.getJournalResourceLocator();

				if (journalLocator.hasJournalId())
					addEntryFieldValue(EntryField.PUBLICATIONS, journalLocator.getNLMid());

				Jinfo = currpubli.getJournalResourceLocator().getName();
				if (journalLocator.hasJournalId())
					Jinfo += " - " + currpubli.getJournalResourceLocator().getMedAbbrev(); // Index name and abbrev in the same token

				addEntryFieldValue(EntryField.PUBLICATIONS,Jinfo);
			}
			String title = currpubli.getTitle();
			if(title.length() > 0) addEntryFieldValue(EntryField.PUBLICATIONS,title);
			SortedSet<PublicationAuthor> authors = currpubli.getAuthors();
			for (PublicationAuthor currauthor : authors) {
				String forename = currauthor.getForeName();
				if(forename.contains(".")) // Submission author
					addEntryFieldValue(EntryField.PUBLICATIONS, currauthor.getLastName() + "  " + currauthor.getInitials());
				else if(!forename.isEmpty() ) // trim not to add spaces when forename/initials are empty
					addEntryFieldValue(EntryField.PUBLICATIONS, (currauthor.getLastName() + " " + forename + " " + currauthor.getInitials()).trim());
				else
					addEntryFieldValue(EntryField.PUBLICATIONS, (currauthor.getLastName() + " " + currauthor.getInitials()).trim());
			}
		}
		
		addEntryFieldValue(EntryField.PUBLI_COMPUTED_COUNT, publi_computed_count);
		addEntryFieldValue(EntryField.PUBLI_CURATED_COUNT, publi_curated_count);
		addEntryFieldValue(EntryField.PUBLI_LARGE_SCALE_COUNT, publi_large_scale_count);

		// Based on the publications and the protein existence level we can compute informational score
		int pe_level = entry.getOverview().getProteinExistences().getProteinExistence().getLevel();
		
		float info_score = 0;
		if(pe_level == 1) info_score=12;
		else if(pe_level == 2) info_score=10;
		else if(pe_level == 3 || pe_level == 4) info_score=8;
		else if(pe_level == 5) info_score=5;
		float coeff = 100*publi_curated_count + 25*publi_computed_count + 10*publi_large_scale_count;
		info_score = coeff * info_score / 10;
		addEntryFieldValue(EntryField.INFORMATIONAL_SCORE, info_score);

	}


	@Override
	public Collection<EntryField> getSupportedFields() {
		return Arrays.asList(EntryField.PUBLICATIONS, EntryField.PUBLI_COMPUTED_COUNT, EntryField.PUBLI_CURATED_COUNT, EntryField.PUBLI_LARGE_SCALE_COUNT, EntryField.INFORMATIONAL_SCORE);
	}

}
