package org.nextprot.api.core.utils.exon;

import com.google.common.base.Preconditions;
import org.nextprot.api.core.domain.Exon;

/**
 * Categorise exons according to gene coordinate of a protein isoform
 *
 * Created by fnikitin on 21/07/15.
 */
public class ExonCategorizer {

    private final int startPositionIsoform;
    private final int endPositionIsoform;

    public ExonCategorizer(int startPositionIsoform, int endPositionIsoform) {

        Preconditions.checkArgument(startPositionIsoform>0);
        Preconditions.checkArgument(endPositionIsoform>=startPositionIsoform);

        this.startPositionIsoform = startPositionIsoform;
        this.endPositionIsoform = endPositionIsoform;
    }

    public ExonCategory categorize(Exon exon) {

        int startPositionExon = exon.getFirstPositionOnGene();
        int endPositionExon = exon.getLastPositionOnGene();

        Preconditions.checkArgument(startPositionIsoform <= endPositionIsoform, "The start position of the isoform on the gene " + startPositionIsoform + " can not be bigger than the end " + endPositionIsoform);
        Preconditions.checkArgument(startPositionExon <= endPositionExon, "The start position of the exon on the gene " + startPositionIsoform + " can not be bigger than the end " + endPositionIsoform);

        ExonCategory codingStatus;

        // not coding exons in the beginning of the transcript
        if (endPositionExon < startPositionIsoform) {
            codingStatus = ExonCategory.NOT_CODING_PRE;
            // ************ SPI ******************* EPI *******************
            // **<SPE>***EPE***********************************************
        }

        // end codon or stop only exon
        else if (startPositionExon > endPositionIsoform) {
            // Some kind of hack has probably been done in the db here !!
            // We consider exon to be of kind STOP_ONLY if it is closed to the last coding exon !!
            if (startPositionExon - endPositionIsoform < 3) codingStatus = ExonCategory.STOP_ONLY;
            else codingStatus = ExonCategory.NOT_CODING_POST;

            // ************ SPI ******************* EPI *******************
            // ********************************************SPE*<EPE>*******
        }

        // start codon
        else if (startPositionExon <= startPositionIsoform && endPositionExon < endPositionIsoform) {
            codingStatus = ExonCategory.START;
            // ************ SPI ******************* EPI *******************
            // *******SPE**********<EPE>***********************************
        }

        // end codon
        else if (endPositionExon >= endPositionIsoform && startPositionExon > startPositionIsoform && startPositionExon <= endPositionIsoform) {
            codingStatus = ExonCategory.STOP;
            // ************ SPI ******************* EPI *******************
            // *********************<SPE>******************EPE*************
        }

        // Case where only one exon can translate the whole isoform
        else if (startPositionExon <= startPositionIsoform && endPositionExon >= endPositionIsoform) {
            codingStatus = ExonCategory.MONO;
            // ************ SPI ******************* EPI *******************
            // *************SPE**********************************EPE*******
        } else {

            // In the last case it must be a coding exon
            codingStatus = ExonCategory.CODING;
        }

        return codingStatus;
    }
}