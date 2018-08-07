package org.nextprot.api.commons.bio.variation.prot.impl.format;

import com.google.common.collect.Lists;
import org.nextprot.api.commons.bio.AminoAcidCode;
import org.nextprot.api.commons.bio.variation.prot.SequenceVariation;
import org.nextprot.api.commons.bio.variation.prot.SequenceVariationFormat;
import org.nextprot.api.commons.bio.variation.prot.impl.seqchange.AminoAcidModification;
import org.nextprot.api.commons.bio.variation.prot.impl.seqchange.Glycosylation;
import org.nextprot.api.commons.bio.variation.prot.impl.seqchange.format.SingleGenericModificationBEDFormat;
import org.nextprot.api.commons.bio.variation.prot.impl.varseq.format.AminoAcidModificationBEDFormatter;
import org.nextprot.api.commons.bio.variation.prot.seqchange.SequenceChange;
import org.nextprot.api.commons.bio.variation.prot.seqchange.SequenceChangeFormat;

import java.util.Collection;

public class SequenceModificationBedFormat extends SequenceVariationFormat {

    private final AminoAcidModificationBEDFormatter aminoAcidModificationFormatter;
    private final SingleGenericModificationBEDFormat changeFormat;

    public SequenceModificationBedFormat() {

        aminoAcidModificationFormatter = new AminoAcidModificationBEDFormatter();
        changeFormat = new SingleGenericModificationBEDFormat();
    }

    @Override
    protected AminoAcidModificationBEDFormatter getChangingSequenceFormatter() {

        return aminoAcidModificationFormatter;
    }

    @Override
    protected SequenceChangeFormat getSequenceChangeFormat(SequenceChange.Type changeType) {

        return changeFormat;
    }

    @Override
    protected Collection<SequenceChange.Type> getAvailableChangeTypes() {

        return Lists.newArrayList(SequenceChange.Type.PTM);
    }

    @Override
    public String format(SequenceVariation variation, AminoAcidCode.CodeType type) {

        StringBuilder sb = new StringBuilder();

        if (variation.getSequenceChange() instanceof Glycosylation) {
            changeFormat.format(sb, (AminoAcidModification) variation.getSequenceChange(), type);
            aminoAcidModificationFormatter.format(variation, type, sb);

            return sb.toString();
        }

        throw new IllegalArgumentException("Not a PTM: cannot format variation "+variation.getSequenceChange());
    }
}
