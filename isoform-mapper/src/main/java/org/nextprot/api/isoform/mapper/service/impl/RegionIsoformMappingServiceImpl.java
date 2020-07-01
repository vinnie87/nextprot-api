package org.nextprot.api.isoform.mapper.service.impl;

import org.apache.commons.logging.Log;
import org.apache.commons.logging.LogFactory;
import org.apache.commons.text.similarity.LevenshteinDistance;
import org.nextprot.api.commons.app.ApplicationContextProvider;
import org.nextprot.api.core.domain.Isoform;
import org.nextprot.api.core.service.IsoformService;
import org.nextprot.api.core.service.impl.AnnotationServiceImpl;
import org.nextprot.api.core.utils.seqmap.IsoformSequencePositionMapper;
import org.nextprot.api.isoform.mapper.domain.query.RegionalFeatureQuery;
import org.nextprot.api.isoform.mapper.domain.query.result.impl.BaseFeatureQueryResult;
import org.nextprot.api.isoform.mapper.domain.query.result.impl.RegionFeatureQuerySuccessImpl;
import org.nextprot.api.isoform.mapper.service.RegionIsoformMappingService;
import org.springframework.beans.factory.annotation.Autowired;
import org.springframework.stereotype.Service;

import java.util.logging.Logger;


/**
 * Regional isoform mapping service, which maps regional features on to other isoforms
 */
@Service
public class RegionIsoformMappingServiceImpl implements RegionIsoformMappingService {

    private static final Log LOGGER = LogFactory.getLog(AnnotationServiceImpl.class);

    @Autowired
    private IsoformService isoformService;

    /**
     * Propagates feature on to other isoforms given the feature query
     * @param query
     * @return result
     */
    public BaseFeatureQueryResult propagateFeature(RegionalFeatureQuery query) {

        // Should build the isoform from the query
        Isoform isoform = ApplicationContextProvider.getApplicationContext().getBean(IsoformService.class)
                .getIsoformByNameOrCanonical(query.getAccession());

        // Validates
        if(validate(query, isoform)) {

            // Propagate
            return propagate(query, isoform);
        } else {
            RegionFeatureQuerySuccessImpl result = new RegionFeatureQuerySuccessImpl(query, isoform);
            return result;
        }
    }

    /**
     * Validates if a given subsequence exists in a given isoform
     * @param query
     * @param isoform
     * @return validity
     */
    private boolean validate(RegionalFeatureQuery query, Isoform isoform) {
        int regionStart = query.getRegionStart();
        int regionEnd = query.getRegionEnd();
        String regionFromQuery = query.getRegionSequence();
        String regionFromIsoform = isoform.getSequence().substring(regionStart - 1, regionEnd);

        if(regionFromQuery != null && regionFromIsoform != null) {
            // Allow tolerance
            return matchWithTolerance(regionFromQuery,regionFromIsoform, 0.05);
        } else {
            return false;
        }
    }

    /**
     * Propagates the given region on the protein sequence into other isoforms
     * @param query
     * @param isoform
     * @return result
     */
    private BaseFeatureQueryResult propagate(RegionalFeatureQuery query, Isoform isoform) {

        int regionStart = query.getRegionStart();
        int regionEnd = query.getRegionEnd();
        String regionFromQuery = query.getRegionSequence();
        int regionLength = regionFromQuery.length();

        RegionFeatureQuerySuccessImpl result = new RegionFeatureQuerySuccessImpl(query, isoform);
        for (Isoform targetIsoform : isoformService.getOtherIsoforms(isoform.getIsoformAccession())) {

            // Propagate the first position
            Integer targetIsoformRegionStart = IsoformSequencePositionMapper.getProjectedPosition(isoform, regionStart, targetIsoform);
            Integer targetIsoformRegionEnd = IsoformSequencePositionMapper.getProjectedPosition(isoform, regionEnd, targetIsoform);

            if(targetIsoformRegionStart == null || targetIsoformRegionEnd == null) {
                LOGGER.info("Project start/end position does not exist on " + targetIsoform.getIsoformAccession());
                continue;
            } else {

                String targetIsoformRegion = targetIsoform.getSequence().substring(targetIsoformRegionStart - 1, targetIsoformRegionEnd);
                LevenshteinDistance editDistanceCalculator = new LevenshteinDistance();
                int editDistance = editDistanceCalculator.apply(regionFromQuery, targetIsoformRegion);
                LOGGER.info("Sequence1: " + regionFromQuery + ",Sequence2:" + targetIsoformRegion + ",Levenshtein Distance:" + editDistance);
                float matchingScore = (float)(regionLength - editDistance)/regionLength;
                if(matchingScore >= 0.96) {
                    LOGGER.info("Match:true,Isoform:" + isoform.getIsoformAccession());
                    result.addMappedFeature(targetIsoform, targetIsoformRegionStart, targetIsoformRegionEnd);
                } else {
                    LOGGER.info("Match:true,Isoform:" + isoform.getIsoformAccession());
                    continue;
                }
            }
        }

        return result;
    }


    private boolean matchWithTolerance(String s, String r, double tolerance) {
        if(s.length() != r.length()) {
            LOGGER.warn("Sequences are not in the same length: " + s + " " + r);
            return false;
        }
        int unmatchCount = s.length();
        boolean matched = true;
        for(int i = 0; i < s.length(); i++) {
            if(s.charAt(i) != r.charAt(i)) {
                unmatchCount--;
                if((float)(unmatchCount/s.length()) > tolerance ) {
                    matched = false;
                    break;
                }
            }
        }

        return matched;
    }
}

