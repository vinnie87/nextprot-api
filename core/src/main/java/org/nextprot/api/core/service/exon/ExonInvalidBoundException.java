package org.nextprot.api.core.service.exon;

import org.nextprot.api.core.domain.exon.UncategorizedExon;

public class ExonInvalidBoundException extends InvalidExonException {

    public ExonInvalidBoundException(UncategorizedExon exon) {

        super(exon, "invalid bounds ["+exon.getFirstPositionOnGene()+"-"+exon.getLastPositionOnGene()+"]");
    }
}