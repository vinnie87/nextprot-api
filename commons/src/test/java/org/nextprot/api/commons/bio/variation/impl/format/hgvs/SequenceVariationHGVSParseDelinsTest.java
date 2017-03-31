package org.nextprot.api.commons.bio.variation.impl.format.hgvs;

import org.junit.Assert;
import org.junit.Test;
import org.nextprot.api.commons.bio.AminoAcidCode;
import org.nextprot.api.commons.bio.variation.SequenceChange;
import org.nextprot.api.commons.bio.variation.SequenceVariation;

public class SequenceVariationHGVSParseDelinsTest {

    SequenceVariantHGVSFormat format = new SequenceVariantHGVSFormat();

    ///// DELETION FOLLOWED BY INSERTIONS

    @Test
    public void testParseDeletion1AaAndInsertion1() throws Exception {

        SequenceVariation pm = format.parse("p.T399delinsL");

        Assert.assertEquals(AminoAcidCode.THREONINE, pm.getChangingSequence().getFirstAminoAcid());
        Assert.assertEquals(399, pm.getChangingSequence().getFirstAminoAcidPos());
        Assert.assertEquals(SequenceChange.Type.DELETION_INSERTION, pm.getSequenceChange().getType());
        Assert.assertArrayEquals(AminoAcidCode.valueOfAminoAcidCodeSequence("L"), (AminoAcidCode[]) pm.getSequenceChange().getValue());
    }

    @Test
    public void testParseDeletionMultiAasAndInsertion1() throws Exception {

        SequenceVariation pm = format.parse("p.L330_A331delinsF");

        Assert.assertEquals(AminoAcidCode.LEUCINE, pm.getChangingSequence().getFirstAminoAcid());
        Assert.assertEquals(AminoAcidCode.ALANINE, pm.getChangingSequence().getLastAminoAcid());
        Assert.assertEquals(330, pm.getChangingSequence().getFirstAminoAcidPos());
        Assert.assertEquals(331, pm.getChangingSequence().getLastAminoAcidPos());
        Assert.assertEquals(SequenceChange.Type.DELETION_INSERTION, pm.getSequenceChange().getType());
        Assert.assertArrayEquals(AminoAcidCode.valueOfAminoAcidCodeSequence("F"), (AminoAcidCode[]) pm.getSequenceChange().getValue());
    }

    @Test
    public void testParseDeletionMultiAndInsertionMulti() throws Exception {

        SequenceVariation pm = format.parse("p.D419_R420delinsSSDG");

        Assert.assertEquals(AminoAcidCode.ASPARTIC_ACID, pm.getChangingSequence().getFirstAminoAcid());
        Assert.assertEquals(AminoAcidCode.ARGININE, pm.getChangingSequence().getLastAminoAcid());
        Assert.assertEquals(419, pm.getChangingSequence().getFirstAminoAcidPos());
        Assert.assertEquals(420, pm.getChangingSequence().getLastAminoAcidPos());
        Assert.assertEquals(SequenceChange.Type.DELETION_INSERTION, pm.getSequenceChange().getType());
        Assert.assertArrayEquals(AminoAcidCode.valueOfAminoAcidCodeSequence("SSDG"), (AminoAcidCode[]) pm.getSequenceChange().getValue());
    }

    @Test
    public void testParseDeletion1AaAndInsertion1Code3() throws Exception {

        SequenceVariation pm = format.parse("p.Thr399delinsLeu");

        Assert.assertEquals(AminoAcidCode.THREONINE, pm.getChangingSequence().getFirstAminoAcid());
        Assert.assertEquals(399, pm.getChangingSequence().getFirstAminoAcidPos());
        Assert.assertEquals(SequenceChange.Type.DELETION_INSERTION, pm.getSequenceChange().getType());
        Assert.assertArrayEquals(AminoAcidCode.valueOfAminoAcidCodeSequence("L"), (AminoAcidCode[]) pm.getSequenceChange().getValue());
    }

    @Test
    public void testParseDeletionMultiAasAndInsertion1Code3() throws Exception {

        SequenceVariation pm = format.parse("p.Leu330_Ala331delinsPhe");

        Assert.assertEquals(AminoAcidCode.LEUCINE, pm.getChangingSequence().getFirstAminoAcid());
        Assert.assertEquals(AminoAcidCode.ALANINE, pm.getChangingSequence().getLastAminoAcid());
        Assert.assertEquals(330, pm.getChangingSequence().getFirstAminoAcidPos());
        Assert.assertEquals(331, pm.getChangingSequence().getLastAminoAcidPos());
        Assert.assertEquals(SequenceChange.Type.DELETION_INSERTION, pm.getSequenceChange().getType());
        Assert.assertArrayEquals(AminoAcidCode.valueOfAminoAcidCodeSequence("F"), (AminoAcidCode[]) pm.getSequenceChange().getValue());
    }

    @Test
    public void testParseDeletionMultiAndInsertionMultiCode3() throws Exception {

        SequenceVariation pm = format.parse("p.Asp419_Arg420delinsSerSerAspGly");

        Assert.assertEquals(AminoAcidCode.ASPARTIC_ACID, pm.getChangingSequence().getFirstAminoAcid());
        Assert.assertEquals(AminoAcidCode.ARGININE, pm.getChangingSequence().getLastAminoAcid());
        Assert.assertEquals(419, pm.getChangingSequence().getFirstAminoAcidPos());
        Assert.assertEquals(420, pm.getChangingSequence().getLastAminoAcidPos());
        Assert.assertEquals(SequenceChange.Type.DELETION_INSERTION, pm.getSequenceChange().getType());
        Assert.assertArrayEquals(AminoAcidCode.valueOfAminoAcidCodeSequence("SSDG"), (AminoAcidCode[]) pm.getSequenceChange().getValue());
    }

    @Test
    public void testParseAANonStandardDelins1() throws Exception {

        format = new SequenceVariantHGVSFormat(SequenceVariantHGVSFormat.ParsingMode.PERMISSIVE);
        SequenceVariation pm = format.parse("p.T399>L");

        Assert.assertEquals(AminoAcidCode.THREONINE, pm.getChangingSequence().getFirstAminoAcid());
        Assert.assertEquals(399, pm.getChangingSequence().getFirstAminoAcidPos());
        Assert.assertEquals(SequenceChange.Type.DELETION_INSERTION, pm.getSequenceChange().getType());
        Assert.assertArrayEquals(AminoAcidCode.valueOfAminoAcidCodeSequence("L"), (AminoAcidCode[]) pm.getSequenceChange().getValue());
    }

    @Test
    public void testParseAANonStandardDelins2() throws Exception {

        format = new SequenceVariantHGVSFormat(SequenceVariantHGVSFormat.ParsingMode.PERMISSIVE);
        SequenceVariation pm = format.parse("p.L330_A331>F");

        Assert.assertEquals(AminoAcidCode.LEUCINE, pm.getChangingSequence().getFirstAminoAcid());
        Assert.assertEquals(AminoAcidCode.ALANINE, pm.getChangingSequence().getLastAminoAcid());
        Assert.assertEquals(330, pm.getChangingSequence().getFirstAminoAcidPos());
        Assert.assertEquals(331, pm.getChangingSequence().getLastAminoAcidPos());
        Assert.assertEquals(SequenceChange.Type.DELETION_INSERTION, pm.getSequenceChange().getType());
        Assert.assertArrayEquals(AminoAcidCode.valueOfAminoAcidCodeSequence("F"), (AminoAcidCode[]) pm.getSequenceChange().getValue());
    }

    @Test
    public void testParseAANonStandardDelins3() throws Exception {

        format = new SequenceVariantHGVSFormat(SequenceVariantHGVSFormat.ParsingMode.PERMISSIVE);
        SequenceVariation pm = format.parse("p.W39_E40>*");

        Assert.assertEquals(AminoAcidCode.TRYPTOPHAN, pm.getChangingSequence().getFirstAminoAcid());
        Assert.assertEquals(AminoAcidCode.GLUTAMIC_ACID, pm.getChangingSequence().getLastAminoAcid());
        Assert.assertEquals(39, pm.getChangingSequence().getFirstAminoAcidPos());
        Assert.assertEquals(40, pm.getChangingSequence().getLastAminoAcidPos());
        Assert.assertEquals(SequenceChange.Type.DELETION_INSERTION, pm.getSequenceChange().getType());
        Assert.assertArrayEquals(AminoAcidCode.valueOfAminoAcidCodeSequence("*"), (AminoAcidCode[]) pm.getSequenceChange().getValue());
    }

    @Test
    public void testParseAANonStandardDelins4() throws Exception {

        format = new SequenceVariantHGVSFormat(SequenceVariantHGVSFormat.ParsingMode.PERMISSIVE);
        SequenceVariation pm = format.parse("p.D419_R420>SSDG");

        Assert.assertEquals(AminoAcidCode.ASPARTIC_ACID, pm.getChangingSequence().getFirstAminoAcid());
        Assert.assertEquals(AminoAcidCode.ARGININE, pm.getChangingSequence().getLastAminoAcid());
        Assert.assertEquals(419, pm.getChangingSequence().getFirstAminoAcidPos());
        Assert.assertEquals(420, pm.getChangingSequence().getLastAminoAcidPos());
        Assert.assertEquals(SequenceChange.Type.DELETION_INSERTION, pm.getSequenceChange().getType());
        Assert.assertArrayEquals(AminoAcidCode.valueOfAminoAcidCodeSequence("SSDG"), (AminoAcidCode[]) pm.getSequenceChange().getValue());
    }

    @Test
    public void testParseAmbiguousAminoAcidsFixCALIPHOMISC568() throws Exception {

        format = new SequenceVariantHGVSFormat(SequenceVariantHGVSFormat.ParsingMode.PERMISSIVE);
        SequenceVariation pm = format.parse("p.Phe709delinsXaaTrpXaa");

        Assert.assertEquals(AminoAcidCode.PHENYLALANINE, pm.getChangingSequence().getFirstAminoAcid());
        Assert.assertEquals(AminoAcidCode.PHENYLALANINE, pm.getChangingSequence().getLastAminoAcid());
        Assert.assertEquals(709, pm.getChangingSequence().getFirstAminoAcidPos());
        Assert.assertEquals(709, pm.getChangingSequence().getLastAminoAcidPos());
        Assert.assertEquals(SequenceChange.Type.DELETION_INSERTION, pm.getSequenceChange().getType());
        Assert.assertArrayEquals(AminoAcidCode.valueOfAminoAcidCodeSequence("XWX"), (AminoAcidCode[]) pm.getSequenceChange().getValue());
    }
}