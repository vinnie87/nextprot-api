package org.nextprot.api.commons.bio.variation.prot.impl.seqchange.format;

import org.nextprot.api.commons.bio.AminoAcidCode;
import org.nextprot.api.commons.bio.variation.prot.SequenceVariation;
import org.nextprot.api.commons.bio.variation.prot.SequenceVariationBuilder;
import org.nextprot.api.commons.bio.variation.prot.impl.seqchange.Glycosylation;
import org.nextprot.api.commons.bio.variation.prot.seqchange.SequenceChangeFormat;

import java.text.ParseException;
import java.util.regex.Matcher;
import java.util.regex.Pattern;


/**
 * Parse single glycosylation with the format:
 *
 * PTM-id_pos (example: PTM-0253_Asn21 represents a glycosylation of an asparagine at position 21)
 */
public class SingleGlycosylationBEDFormat implements SequenceChangeFormat<Glycosylation> {

    private static final Pattern PATTERN = Pattern.compile("^(PTM-\\d{4})_([A-Z])([a-z]{2})?(\\d+)$");

    @Override
    public SequenceVariation parse(String source, SequenceVariationBuilder.FluentBuilding builder) throws ParseException {

        Matcher m = PATTERN.matcher(source);

        if (m.matches()) {

            Glycosylation aaChange = new Glycosylation(m.group(1));
            AminoAcidCode affectedAA = AminoAcidCode.parseAminoAcidCode(m.group(2) + ((m.group(3) != null) ? m.group(3) : ""));
            int affectedAAPos = Integer.parseInt(m.group(4));

            return builder.selectAminoAcid(affectedAA, affectedAAPos).thenAddModification(aaChange).build();
        }

        return null;
    }

    @Override
    public boolean matches(String source) {
        return source.matches(PATTERN.pattern());
    }

    @Override
    public void format(StringBuilder sb, Glycosylation change, AminoAcidCode.CodeType type) {

        sb
                .append(change.getPTMId())
                .append("_");
    }
}