package org.nextprot.api.isoform.mapper.domain.impl;

import com.google.common.collect.Sets;
import org.junit.Assert;
import org.junit.Test;
import org.mockito.Mockito;
import org.nextprot.api.commons.constants.AnnotationCategory;
import org.nextprot.api.core.domain.Entry;
import org.nextprot.api.core.domain.Isoform;
import org.nextprot.api.core.service.BeanService;
import org.nextprot.api.core.service.EntryBuilderService;
import org.nextprot.api.core.service.IsoformService;
import org.nextprot.api.core.service.MasterIdentifierService;
import org.nextprot.api.core.service.fluent.EntryConfig;
import org.nextprot.api.isoform.mapper.domain.SingleFeatureQuery;

import java.text.ParseException;

import static org.mockito.Matchers.any;

public class FeatureQuerySuccessTest {

    @Test
    public void testOnSuccess() throws ParseException {

        SingleFeatureQuery query = new SingleFeatureQuery("NX_Q9UI33", "SCN11A-p.Leu1158Pro", AnnotationCategory.VARIANT.getApiTypeName());

        Entry entry = mockNX_Q9UI33();

        SequenceVariant sequenceVariant = SequenceVariant.variant("SCN11A-p.Leu1158Pro", mockBeanService("SCN11A", entry));

        SingleFeatureQuerySuccessImpl result = new SingleFeatureQuerySuccessImpl(entry, query, sequenceVariant);
        result.addMappedFeature(entry.getIsoforms().get(0), 1158, 1158);
        result.addMappedFeature(entry.getIsoforms().get(1), 1158, 1158);
        result.addMappedFeature(entry.getIsoforms().get(2), 1120, 1120);

        result.addUnmappedFeature(SequenceVariantTest.mockIsoform("NX_Q9UI33-23", "Iso 23", false, ""));

        Assert.assertEquals("SCN11A-iso1-p.Leu1158Pro", result.getData().get("NX_Q9UI33-1").getIsoSpecificFeature());
        Assert.assertEquals("SCN11A-iso2-p.Leu1158Pro", result.getData().get("NX_Q9UI33-2").getIsoSpecificFeature());
        Assert.assertEquals("SCN11A-iso3-p.Leu1120Pro", result.getData().get("NX_Q9UI33-3").getIsoSpecificFeature());
        Assert.assertNull(result.getData().get("NX_Q9UI33-23").getIsoSpecificFeature());
    }

    private static Entry mockNX_Q9UI33() {

        String seqIso1 = "MDDRCYPVIFPDERNFRPFTSDSLAAIEKRIAIQKEKKKSKDQTGEVPQPRPQLDLKASRKLPKLYGDIPRELIGKPLEDLDPFYRNHKTFMVLNRKRTIYRFSAKHALFIFGPFNSIRSLAIRVSVHSLFSMFIIGTVIINCVFMATGPAKNSNSNNTDIAECVFTGIYIFEALIKILARGFILDEFSFLRDPWNWLDSIVIGIAIVSYIPGITIKLLPLRTFRVFRALKAISVVSRLKVIVGALLRSVKKLVNVIILTFFCLSIFALVGQQLFMGSLNLKCISRDCKNISNPEAYDHCFEKKENSPEFKMCGIWMGNSACSIQYECKHTKINPDYNYTNFDNFGWSFLAMFRLMTQDSWEKLYQQTLRTTGLYSVFFFIVVIFLGSFYLINLTLAVVTMAYEEQNKNVAAEIEAKEKMFQEAQQLLKEEKEALVAMGIDRSSLTSLETSYFTPKKRKLFGNKKRKSFFLRESGKDQPPGSDSDEDCQKKPQLLEQTKRLSQNLSLDHFDEHGDPLQRQRALSAVSILTITMKEQEKSQEPCLPCGENLASKYLVWNCCPQWLCVKKVLRTVMTDPFTELAITICIIINTVFLAMEHHKMEASFEKMLNIGNLVFTSIFIAEMCLKIIALDPYHYFRRGWNIFDSIVALLSFADVMNCVLQKRSWPFLRSFRVLRVFKLAKSWPTLNTLIKIIGNSVGALGSLTVVLVIVIFIFSVVGMQLFGRSFNSQKSPKLCNPTGPTVSCLRHWHMGDFWHSFLVVFRILCGEWIENMWECMQEANASSSLCVIVFILITVIGKLVVLNLFIALLLNSFSNEERNGNLEGEARKTKVQLALDRFRRAFCFVRHTLEHFCHKWCRKQNLPQQKEVAGGCAAQSKDIIPLVMEMKRGSETQEELGILTSVPKTLGVRHDWTWLAPLAEEEDDVEFSGEDNAQRITQPEPEQQAYELHQENKKPTSQRVQSVEIDMFSEDEPHLTIQDPRKKSDVTSILSECSTIDLQDGFGWLPEMVPKKQPERCLPKGFGCCFPCCSVDKRKPPWVIWWNLRKTCYQIVKHSWFESFIIFVILLSSGALIFEDVHLENQPKIQELLNCTDIIFTHIFILEMVLKWVAFGFGKYFTSAWCCLDFIIVIVSVTTLINLMELKSFRTLRALRPLRALSQFEGMKVVVNALIGAIPAILNVLLVCLIFWLVFCILGVYFFSGKFGKCINGTDSVINYTIITNKSQCESGNFSWINQKVNFDNVGNAYLALLQVATFKGWMDIIYAAVDSTEKEQQPEFESNSLGYIYFVVFIIFGSFFTLNLFIGVIIDNFNQQQKKLGGQDIFMTEEQKKYYNAMKKLGSKKPQKPIPRPLNKCQGLVFDIVTSQIFDIIIISLIILNMISMMAESYNQPKAMKSILDHLNWVFVVIFTLECLIKIFALRQYYFTNGWNLFDCVVVLLSIVSTMISTLENQEHIPFPPTLFRIVRLARIGRILRLVRAARGIRTLLFALMMSLPSLFNIGLLLFLIMFIYAILGMNWFSKVNPESGIDDIFNFKTFASSMLCLFQISTSAGWDSLLSPMLRSKESCNSSSENCHLPGIATSYFVSYIIISFLIVVNMYIAVILENFNTATEESEDPLGEDDFDIFYEVWEKFDPEATQFIKYSALSDFADALPEPLRVAKPNKYQFLVMDLPMVSEDRLHCMDILFAFTARVLGGSDGLDSMKAMMEEKFMEANPLKKLYEPIVTTTKRKEEERGAAIIQKAFRKYMMKVTKGDQGDQNDLENGPHSPLQTLCNGDLSSFGVAKGKVHCD";
        String seqIso2 = "MDDRCYPVIFPDERNFRPFTSDSLAAIEKRIAIQKEKKKSKDQTGEVPQPRPQLDLKASRKLPKLYGDIPRELIGKPLEDLDPFYRNHKTFMVLNRKRTIYRFSAKHALFIFGPFNSIRSLAIRVSVHSLFSMFIIGTVIINCVFMATGPAKNSNSNNTDIAECVFTGIYIFEALIKILARGFILDEFSFLRDPWNWLDSIVIGIAIVSYIPGITIKLLPLRTFRVFRALKAISVVSRLKVIVGALLRSVKKLVNVIILTFFCLSIFALVGQQLFMGSLNLKCISRDCKNISNPEAYDHCFEKKENSPEFKMCGIWMGNSACSIQYECKHTKINPDYNYTNFDNFGWSFLAMFRLMTQDSWEKLYQQTLRTTGLYSVFFFIVVIFLGSFYLINLTLAVVTMAYEEQNKNVAAEIEAKEKMFQEAQQLLKEEKEALVAMGIDRSSLTSLETSYFTPKKRKLFGNKKRKSFFLRESGKDQPPGSDSDEDCQKKPQLLEQTKRLSQNLSLDHFDEHGDPLQRQRALSAVSILTITMKEQEKSQEPCLPCGENLASKYLVWNCCPQWLCVKKVLRTVMTDPFTELAITICIIINTVFLAMEHHKMEASFEKMLNIGNLVFTSIFIAEMCLKIIALDPYHYFRRGWNIFDSIVALLSFADVMNCVLQKRSWPFLRSFRVLRVFKLAKSWPTLNTLIKIIGNSVGALGSLTVVLVIVIFIFSVVGMQLFGRSFNSQKSPKLCNPTGPTVSCLRHWHMGDFWHSFLVVFRILCGEWIENMWECMQEANASSSLCVIVFILITVIGKLVVLNLFIALLLNSFSNEERNGNLEGEARKTKVQLALDRFRRAFCFVRHTLEHFCHKWCRKQNLPQQKEVAGGCAAQSKDIIPLVMEMKRGSETQEELGILTSVPKTLGVRHDWTWLAPLAEEEDDVEFSGEDNAQRITQPEPEQQAYELHQENKKPTSQRVQSVEIDMFSEDEPHLTIQDPRKKSDVTSILSECSTIDLQDGFGWLPEMVPKKQPERCLPKGFGCCFPCCSVDKRKPPWVIWWNLRKTCYQIVKHSWFESFIIFVILLSSGALIFEDVHLENQPKIQELLNCTDIIFTHIFILEMVLKWVAFGFGKYFTSAWCCLDFIIVIVSVTTLINLMELKSFRTLRALRPLRALSQFEGMKVVVNALIGAIPAILNVLLVCLIFWLVFCILGVYFFSGKFGKCINGTDSVINYTIITNKSQCESGNFSWINQKVNFDNVGNAYLALLQVATFKGWMDIIYAAVDSTEKEQQPEFESNSLGYIYFVVFIIFGSFFTLNLFIGVIIDNFNQQQKKLGGQDIFMTEEQKKYYNAMKKLGSKKPQKPIPRPLNKCQGLVFDIVTSQIFDIIIISLIILNMISMMAESYNQPKAMKSILDHLNWVFVVIFTLECLIKIFALRQYYFTNGWNLFDCVVVLLSIVSK";
        String seqIso3 = "MDDRCYPVIFPDERNFRPFTSDSLAAIEKRIAIQKEKKKSKDQTGEVPQPRPQLDLKASRKLPKLYGDIPRELIGKPLEDLDPFYRNHKTFMVLNRKRTIYRFSAKHALFIFGPFNSIRSLAIRVSVHSLFSMFIIGTVIINCVFMATGPAKNSNSNNTDIAECVFTGIYIFEALIKILARGFILDEFSFLRDPWNWLDSIVIGIAIVSYIPGITIKLLPLRTFRVFRALKAISVVSRLKVIVGALLRSVKKLVNVIILTFFCLSIFALVGQQLFMGSLNLKCISRDCKNISNPEAYDHCFEKKENSPEFKMCGIWMGNSACSIQYECKHTKINPDYNYTNFDNFGWSFLAMFRLMTQDSWEKLYQQTLRTTGLYSVFFFIVVIFLGSFYLINLTLAVVTMAYEEQNKNVAAEIEAKEKMFQEAQQLLKEEKEALVAMGIDRSSLTSLETSYFTPKKRKLFGNKKRKSFFLRESGKDQPPGSDSDEDCQKKPQLLEQTKRLSQNLSLDHFDEHGDPLQRQRALSAVSILTITMKEQEKSQEPCLPCGENLASKYLVWNCCPQWLCVKKVLRTVMTDPFTELAITICIIINTVFLAMEHHKMEASFEKMLNIGNLVFTSIFIAEMCLKIIALDPYHYFRRGWNIFDSIVALLSFADVMNCVLQKRSWPFLRSFRVLRVFKLAKSWPTLNTLIKIIGNSVGALGSLTVVLVIVIFIFSVVGMQLFGRSFNSQKSPKLCNPTGPTVSCLRHWHMGDFWHSFLVVFRILCGEWIENMWECMQEANASSSLCVIVFILITVIGKLVVLNLFIALLLNSFSNEERNGNLEGEARKTKVQLALDRFRRAFCFVRHTLEHFCHKWCRKQNLPQQKEVAGGCAAQSKDIIPLVMEMKRGSETQEELGILTSVPKTLGVRHDWTWLAPLAEEEDDVEFSGEDNAQRITQPEPEQQKSDVTSILSECSTIDLQDGFGWLPEMVPKKQPERCLPKGFGCCFPCCSVDKRKPPWVIWWNLRKTCYQIVKHSWFESFIIFVILLSSGALIFEDVHLENQPKIQELLNCTDIIFTHIFILEMVLKWVAFGFGKYFTSAWCCLDFIIVIVSVTTLINLMELKSFRTLRALRPLRALSQFEGMKVVVNALIGAIPAILNVLLVCLIFWLVFCILGVYFFSGKFGKCINGTDSVINYTIITNKSQCESGNFSWINQKVNFDNVGNAYLALLQVATFKGWMDIIYAAVDSTEKEQQPEFESNSLGYIYFVVFIIFGSFFTLNLFIGVIIDNFNQQQKKLGGQDIFMTEEQKKYYNAMKKLGSKKPQKPIPRPLNKCQGLVFDIVTSQIFDIIIISLIILNMISMMAESYNQPKAMKSILDHLNWVFVVIFTLECLIKIFALRQYYFTNGWNLFDCVVVLLSIVSTMISTLENQEHIPFPPTLFRIVRLARIGRILRLVRAARGIRTLLFALMMSLPSLFNIGLLLFLIMFIYAILGMNWFSKVNPESGIDDIFNFKTFASSMLCLFQISTSAGWDSLLSPMLRSKESCNSSSENCHLPGIATSYFVSYIIISFLIVVNMYIAVILENFNTATEESEDPLGEDDFDIFYEVWEKFDPEATQFIKYSALSDFADALPEPLRVAKPNKYQFLVMDLPMVSEDRLHCMDILFAFTARVLGGSDGLDSMKAMMEEKFMEANPLKKLYEPIVTTTKRKEEERGAAIIQKAFRKYMMKVTKGDQGDQNDLENGPHSPLQTLCNGDLSSFGVAKGKVHCD";

        Isoform iso1 = SequenceVariantTest.mockIsoform("NX_Q9UI33-1", "Iso 1", true, seqIso1);
        Isoform iso2 = SequenceVariantTest.mockIsoform("NX_Q9UI33-2", "Iso 2", false, seqIso2);
        Isoform iso3 = SequenceVariantTest.mockIsoform("NX_Q9UI33-3", "Iso 3", false, seqIso3);

        return SequenceVariantTest.mockEntry("NX_Q9UI33", iso1, iso2, iso3);
    }

    private static BeanService mockBeanService(String geneName, Entry entry) {

        BeanService service = Mockito.mock(BeanService.class);

        MasterIdentifierService masterIdentifierService = mockMasterIdentifierService(geneName, entry.getUniqueName());
        EntryBuilderService entryBuilderService = mockEntryBuilderService(entry);
        IsoformService isoformService = mockIsoformService();

        Mockito.when(service.getBean(MasterIdentifierService.class)).thenReturn(masterIdentifierService);
        Mockito.when(service.getBean(EntryBuilderService.class)).thenReturn(entryBuilderService);
        Mockito.when(service.getBean(IsoformService.class)).thenReturn(isoformService);

        return service;
    }

    private static EntryBuilderService mockEntryBuilderService(Entry entry) {

        EntryBuilderService service = Mockito.mock(EntryBuilderService.class);

        Mockito.when(service.build(any(EntryConfig.class))).thenReturn(entry);

        return service;
    }

    private static MasterIdentifierService mockMasterIdentifierService(String geneName, String entryAccession) {

        MasterIdentifierService service = Mockito.mock(MasterIdentifierService.class);

        Mockito.when(service.findEntryAccessionByGeneName(geneName, false)).thenReturn(Sets.newHashSet(entryAccession));

        return service;
    }

    private static IsoformService mockIsoformService() {

        IsoformService service = Mockito.mock(IsoformService.class);


        return service;
    }
}