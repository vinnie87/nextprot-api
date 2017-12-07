package org.nextprot.api.web.controller;

import org.jsondoc.core.annotation.Api;
import org.jsondoc.core.annotation.ApiMethod;
import org.jsondoc.core.annotation.ApiPathParam;
import org.jsondoc.core.pojo.ApiVerb;
import org.nextprot.api.commons.exception.NextProtException;
import org.nextprot.api.core.domain.Publication;
import org.nextprot.api.core.domain.publication.*;
import org.nextprot.api.core.service.EntryPublicationService;
import org.nextprot.api.core.service.PublicationService;
import org.nextprot.api.core.service.PublicationStatisticsService;
import org.springframework.beans.factory.annotation.Autowired;
import org.springframework.http.MediaType;
import org.springframework.stereotype.Controller;
import org.springframework.web.bind.annotation.PathVariable;
import org.springframework.web.bind.annotation.RequestMapping;
import org.springframework.web.bind.annotation.RequestMethod;
import org.springframework.web.bind.annotation.ResponseBody;

import java.io.Serializable;
import java.util.*;
import java.util.function.Function;
import java.util.function.Predicate;
import java.util.stream.Collectors;

@Controller
@Api(name = "Entry Publications", description = "Method to retrieve a publications linked to a neXtProt entry")
public class EntryPublicationController {

	@Autowired
    private EntryPublicationService entryPublicationService;
    @Autowired
    private PublicationService publicationService;
    @Autowired
    private PublicationStatisticsService publicationStatisticsService;

	@ApiMethod(path = "/entry-publications/entry/{entry}/category/{category}", verb = ApiVerb.GET, description = "Exports publications associated with a neXtProt entry and a publication category",
			produces = { MediaType.APPLICATION_JSON_VALUE })
	@RequestMapping(value = "/entry-publications/entry/{entry}/category/{category}", method = { RequestMethod.GET })
	@ResponseBody
	public List<EntryPublicationView> getEntryPublicationList(
			@ApiPathParam(name = "entry", description = "The name of the neXtProt entry. For example, the insulin: NX_P01308",  allowedvalues = { "NX_P01308"})
			@PathVariable("entry") String entryName,
			@ApiPathParam(name = "category", description = "publication category (CURATED, SUBMISSION, ADDITIONAL, WEB_RESOURCE, PATENT or ALL)", allowedvalues = { "CURATED" })
			@PathVariable(value = "category") String publicationCategory) {

        String pubCategoryName = publicationCategory.toUpperCase();

		if (PublicationCategory.hasName(pubCategoryName)) {

            return buildView(entryPublicationService.findEntryPublications(entryName), PublicationCategory.valueOfName(pubCategoryName));
		}

		throw new NextProtException(publicationCategory + ": Unknown publication view");
	}

	@ApiMethod(path = "/entry-publications/entry/{entry}/pubid/{pubid}", verb = ApiVerb.GET, description = "Exports identified publication associated with a neXtProt entry",
			produces = { MediaType.APPLICATION_JSON_VALUE })
	@RequestMapping(value = "/entry-publications/entry/{entry}/pubid/{pubid}", method = { RequestMethod.GET })
	@ResponseBody
	public EntryPublication getEntryPublication(
			@ApiPathParam(name = "entry", description = "The name of the neXtProt entry. For example, the insulin: NX_P01308",  allowedvalues = { "NX_P01308"})
			@PathVariable("entry") String entryName,
			@ApiPathParam(name = "pubid", description = "A publication id", allowedvalues = { "630194" })
			@PathVariable("pubid") long publicationId) {

		return entryPublicationService.findEntryPublications(entryName).getEntryPublication(publicationId);
	}

    @ApiMethod(path = "/entry-publications/entry/{entry}/count", verb = ApiVerb.GET, description = "Count entry publications associated with a neXtProt entry by publication category",
            produces = { MediaType.APPLICATION_JSON_VALUE })
    @RequestMapping(value = "/entry-publications/entry/{entry}/count", method = { RequestMethod.GET })
    @ResponseBody
    public Map<PublicationCategory, Integer> countEntryPublication(
            @ApiPathParam(name = "entry", description = "The name of the neXtProt entry. For example, the insulin: NX_P01308",  allowedvalues = { "NX_P01308"})
            @PathVariable("entry") String entryAccession) {

	    Map<PublicationCategory, Integer> count = new HashMap<>();

	    for (PublicationCategory view : PublicationCategory.values()) {

            count.put(view, entryPublicationService.findEntryPublications(entryAccession).getEntryPublicationList(view).size());
        }

        return count;
    }

    @ApiMethod(path = "/publications/stats", verb = ApiVerb.GET, description = "Get overall statistics over publications",
            produces = { MediaType.APPLICATION_JSON_VALUE })
    @RequestMapping(value = "/publications/stats", method = { RequestMethod.GET })
    @ResponseBody
    public GlobalPublicationStatistics calcGlobalPublicationStats() {

        return publicationStatisticsService.getGlobalPublicationStatistics();
    }

    @ApiMethod(path = "/publications/pubid/{pubid}/stats", verb = ApiVerb.GET, description = "Get statistics over all publications linked with a neXtProt entry",
            produces = { MediaType.APPLICATION_JSON_VALUE })
    @RequestMapping(value = "/publications/stats/pubid/{pubid}", method = { RequestMethod.GET })
    @ResponseBody
    public GlobalPublicationStatistics.PublicationStatistics calcPublicationStats(
            @ApiPathParam(name = "pubid", description = "A publication id", allowedvalues = { "630194" })
            @PathVariable("pubid") long publicationId) {

        return publicationStatisticsService.getGlobalPublicationStatistics().getPublicationStatistics(publicationId);
    }

    @ApiMethod(path = "/publications/pubids/{statstype}", verb = ApiVerb.GET, description = "Get all publication ids by statistics type",
            produces = { MediaType.APPLICATION_JSON_VALUE })
    @RequestMapping(value = "/publications/pubids/{statstype}", method = { RequestMethod.GET })
    @ResponseBody
    public Set<Long> getPublicationIds(
            @ApiPathParam(name = "statstype", description = "a publication statistics type (ALL, CITED, COMPUTED, CURATED, LARGE_SCALE)", allowedvalues = { "ALL" })
            @PathVariable(value = "statstype") String statisticsType) {

        StatisticsType type = StatisticsType.valueOf(statisticsType.toUpperCase());
        Map<Long, GlobalPublicationStatistics.PublicationStatistics> map = publicationStatisticsService.getGlobalPublicationStatistics().getPublicationStatisticsById();

        if (type == StatisticsType.ALL) {
            return map.keySet();
        }

        return map.values().stream()
                .filter(type.getPredicate())
                .map(ps -> ps.getPublicationId())
                .collect(Collectors.toSet());
    }

    @ApiMethod(path = "/entry-publications/pubid/{pubid}", verb = ApiVerb.GET, description = "Exports identified publication associated with neXtProt entries",
            produces = { MediaType.APPLICATION_JSON_VALUE })
    @RequestMapping(value = "/entry-publications/pubid/{pubid}", method = { RequestMethod.GET })
    @ResponseBody
    public PublicationView getEntryPublicationsByPubId(
            @ApiPathParam(name = "pubid", description = "A publication id", allowedvalues = { "630194" })
            @PathVariable("pubid") long publicationId) {

        PublicationView view = new PublicationView();
        view.setPublication(publicationService.findPublicationById(publicationId));
        view.setEntryPublicationList(publicationService.getEntryPublications(publicationId));

        return view;
    }

    private List<EntryPublicationView> buildView(EntryPublications entryPublications, PublicationCategory publicationCategory) {

        List<EntryPublicationView> list = new ArrayList<>();

        Map<Long, EntryPublication> entryPublicationMap = entryPublications
                .getEntryPublicationList(publicationCategory).stream()
                .collect(Collectors.toMap(
                        EntryPublication::getPubId,
                        Function.identity(),
                        (entryPublication, entryPublication2) -> entryPublication
                ));

        List<Publication> publications = publicationService.findPublicationsByEntryName(entryPublications.getEntryAccession());

        for (Publication publication : publications) {

            if (entryPublicationMap.containsKey(publication.getPublicationId())) {
                EntryPublicationView view = new EntryPublicationView();

                EntryPublication entryPublication = entryPublicationMap.get(publication.getPublicationId());

                view.setCitedInViews(entryPublication.getCitedInViews());
                view.setDirectLinks(entryPublication.getDirectLinks());
                view.setPublication(publication);

                list.add(view);
            }
        }

        return list;
    }

    public static class EntryPublicationView implements Serializable {

        private static final long serialVersionUID = 1L;

        private Publication publication;
        private Map<String,String> citedInViews;
        private List<PublicationDirectLink> directLinks;

        public Publication getPublication() {
            return publication;
        }

        public void setPublication(Publication publication) {
            this.publication = publication;
        }

        public Map<String, String> getCitedInViews() {
            return citedInViews;
        }

        public void setCitedInViews(Map<String, String> citedInViews) {
            this.citedInViews = citedInViews;
        }

        public List<PublicationDirectLink> getDirectLinks() {
            return directLinks;
        }

        public void setDirectLinks(List<PublicationDirectLink> directLinks) {
            this.directLinks = directLinks;
        }
    }

    public static class PublicationView implements Serializable {

        private static final long serialVersionUID = 1L;

        private Publication publication;
        private List<EntryPublication> entryPublicationList;

        public Publication getPublication() {
            return publication;
        }

        public void setPublication(Publication publication) {
            this.publication = publication;
        }

        public List<EntryPublication> getEntryPublicationList() {
            return entryPublicationList;
        }

        public void setEntryPublicationList(List<EntryPublication> entryPublicationList) {
            this.entryPublicationList = entryPublicationList;
        }
    }

    private enum StatisticsType {

        ALL(publicationStatistics -> true),
        CITED(publicationStatistics -> publicationStatistics.isCited()),
        COMPUTED(publicationStatistics -> publicationStatistics.isComputed()),
        CURATED(publicationStatistics -> publicationStatistics.isCurated()),
        LARGE_SCALE(publicationStatistics -> publicationStatistics.isLargeScale())
        ;

        private final Predicate<GlobalPublicationStatistics.PublicationStatistics> predicate;

        StatisticsType(Predicate<GlobalPublicationStatistics.PublicationStatistics> predicate) {
            this.predicate = predicate;
        }

        public Predicate<GlobalPublicationStatistics.PublicationStatistics> getPredicate() {
            return predicate;
        }
    }
}
