package org.nextprot.api.core.service;


import org.nextprot.api.core.domain.IsoformSequenceInfoPeff;
import org.nextprot.api.core.service.impl.peff.ModResPsiFormatter;

import java.util.List;

public interface PeffService {

    IsoformSequenceInfoPeff formatSequenceInfo(String isoformAccession);

    String formatIsoformAccession(String isoformAccession);

    String formatProteinName(String isoformAccession);

    String formatGeneName(String isoformAccession);

    String formatNcbiTaxonomyIdentifier(String isoformAccession);

    String formatTaxonomyName(String isoformAccession);

    String formatSequenceLength(String isoformAccession);

    String formatSequenceVersion(String isoformAccession);

    String formatEntryVersion(String isoformAccession);

    String formatProteinEvidence(String isoformAccession);

    String formatVariantSimple(String isoformAccession);

    String formatVariantComplex(String isoformAccession);

    String formatModResPsi(String isoformAccession, List<ModResPsiFormatter.ModResInfos> unmappedUniprotMods);

    String formatModRes(String isoformAccession);

    String formatProcessedMolecule(String isoformAccession);
}
