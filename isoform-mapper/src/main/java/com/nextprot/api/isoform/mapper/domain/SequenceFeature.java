package com.nextprot.api.isoform.mapper.domain;

import com.nextprot.api.isoform.mapper.domain.impl.exception.UnknownIsoformRuntimeException;
import org.nextprot.api.commons.bio.variation.SequenceVariation;
import org.nextprot.api.core.domain.Entry;
import org.nextprot.api.core.domain.Isoform;

/**
 * A sequence feature on an isoform sequence on a specific gene
 */
public interface SequenceFeature {

    /** @return the gene name */
    String getGeneName();

    /** @return isoform name or null if canonical */
    String getIsoformName();

    /** @return the variation as a string */
    String getFormattedVariation();

    /** format a feature specifically to isoform */
    String formatIsoSpecificFeature(Isoform isoform, int firstAAPos, int lastAAPos);

    /** @return the protein sequence variation */
    SequenceVariation getProteinVariation();

    /** @return true if gene name is referenced in entry */
    boolean isValidGeneName(Entry entry);

    /** @return the entry isoform where lie the feature or ab */
    Isoform getIsoform(Entry entry) throws UnknownIsoformRuntimeException;

    boolean isValidIsoform(Entry entry);
}