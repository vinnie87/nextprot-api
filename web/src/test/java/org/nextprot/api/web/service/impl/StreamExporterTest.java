package org.nextprot.api.web.service.impl;

import org.junit.Ignore;
import org.junit.Test;
import org.nextprot.api.web.dbunit.base.mvc.WebIntegrationBaseTest;
import org.nextprot.api.web.service.ExportService;
import org.nextprot.api.web.utils.XMLUnitUtils;
import org.skyscreamer.jsonassert.JSONAssert;
import org.w3c.dom.NodeList;

import java.io.ByteArrayOutputStream;
import java.io.PrintWriter;
import java.io.Writer;
import java.util.Arrays;
import java.util.HashMap;
import java.util.Map;

import static org.junit.Assert.assertEquals;

public class StreamExporterTest extends WebIntegrationBaseTest {

    @Test
    public void testXMLExportStream() throws Exception {

    	ByteArrayOutputStream out = new ByteArrayOutputStream();
        Writer writer = new PrintWriter(out);
        NPStreamExporter exporter = new XMLStreamExporter();
        exporter.export(Arrays.asList("NX_P06213", "NX_P01308"), writer, "overview", null);

        NodeList recommendedNodes = XMLUnitUtils.getMatchingNodes(out.toString(), "nextprot-export/entry-list/entry/overview/gene-list/gene/gene-name[@type='primary']");
        assertEquals(recommendedNodes.item(0).getTextContent(), "INSR");
        assertEquals(recommendedNodes.item(1).getTextContent(), "INS");
    }

    @Test
    public void testJSONExportStream() throws Exception {

        ByteArrayOutputStream out = new ByteArrayOutputStream();
        Writer writer = new PrintWriter(out);

        NPStreamExporter exporter = new JSONStreamExporter();

        exporter.export(Arrays.asList("NX_P06213", "NX_P01308"), writer, "overview", null);

        JSONAssert.assertEquals("{\"properties\":{\"proteinExistence\":\"Evidence at protein level\",\"ptmCount\":49,\"varCount\":242,\"isoformCount\":2,\"mutagenesisCount\":1,\"interactionCount\":46,\"maxSeqLen\":1382,\"filterstructure\":true,\"filterdisease\":true,\"filtermutagenesis\":true,\"filterproteomics\":true,\"filterexpressionprofile\":true,\"referencesCount\":513,\"referencesSubmissionsCount\":2,\"referencesAdditionalPublicationsCount\":341,\"referencesCuratedPublicationsCount\":190,\"referencesWebResourcesCount\":1,\"referencesPatentsCount\":0},\"uniqueName\":\"NX_P06213\",\"overview\":{\"history\":{\"proteinExistence\":\"Evidence_at_protein_level\",\"nextprotIntegrationDate\":1269233503237,\"nextprotUpdateDate\":1430705017403,\"uniprotIntegrationDate\":568029600000,\"uniprotUpdateDate\":1430301600000,\"uniprotVersion\":\"213\",\"lastSequenceUpdate\":\"2010-10-05\",\"sequenceVersion\":\"4\",\"formattedUniprotIntegrationDate\":\"1988-01-01\",\"formattedUniprotUpdateDate\":\"2015-04-29\",\"proteinExistenceRaw\":\"protein level\",\"formattedNextprotIntegrationDate\":\"2010-03-22\",\"formattedNextprotUpdateDate\":\"2015-05-04\",\"proteinExistenceLevel\":1},\"families\":[{\"familyId\":67489,\"accession\":\"FA-03128\",\"name\":\"Insulin receptor\",\"level\":\"Subfamily\",\"description\":\"Belongs to the protein kinase superfamily. Tyr protein kinase family. Insulin receptor subfamily.\",\"region\":null,\"parent\":{\"familyId\":67478,\"accession\":\"FA-03117\",\"name\":\"Tyr protein kinase\",\"level\":\"Family\",\"description\":null,\"region\":null,\"parent\":{\"familyId\":67418,\"accession\":\"FA-03057\",\"name\":\"Protein kinase\",\"level\":\"Superfamily\",\"description\":null,\"region\":null,\"parent\":null}}}],\"proteinNames\":[{\"clazz\":\"PROTEIN_NAMES\",\"type\":\"name\",\"qualifier\":\"full\",\"id\":\"PR_888984\",\"category\":\"protein\",\"name\":\"Insulin receptor\",\"parentId\":null,\"synonyms\":[{\"clazz\":\"PROTEIN_NAMES\",\"type\":\"name\",\"qualifier\":\"short\",\"id\":\"PR_888983\",\"category\":\"protein\",\"name\":\"IR\",\"parentId\":\"PR_888984\",\"synonyms\":null,\"main\":false,\"synonymId\":\"PR_888983\",\"composedName\":\"short name\",\"synonymName\":\"IR\"},{\"clazz\":\"PROTEIN_NAMES\",\"type\":\"enzyme name\",\"qualifier\":\"EC\",\"id\":\"PR_2472057\",\"category\":\"EC\",\"name\":\"2.7.10.1\",\"parentId\":\"PR_888984\",\"synonyms\":null,\"main\":true,\"synonymId\":\"PR_2472057\",\"composedName\":\"EC enzyme name\",\"synonymName\":\"2.7.10.1\"}],\"main\":true,\"synonymId\":\"PR_888984\",\"composedName\":\"full name\",\"synonymName\":\"Insulin receptor\"}],\"geneNames\":[{\"clazz\":\"GENE_NAMES\",\"type\":\"gene name\",\"qualifier\":null,\"id\":\"PR_1171491\",\"category\":\"gene name\",\"name\":\"INSR\",\"parentId\":null,\"synonyms\":null,\"main\":true,\"synonymId\":\"PR_1171491\",\"composedName\":\"gene name\",\"synonymName\":\"INSR\"}],\"functionalRegionNames\":null,\"cleavedRegionNames\":[{\"clazz\":\"CLEAVED_REGION_NAMES\",\"type\":\"name\",\"qualifier\":\"full\",\"id\":\"MP_10131410\",\"category\":\"protein\",\"name\":\"Insulin receptor subunit alpha\",\"parentId\":null,\"synonyms\":null,\"main\":true,\"synonymId\":\"MP_10131410\",\"composedName\":\"full name\",\"synonymName\":\"Insulin receptor subunit alpha\"},{\"clazz\":\"CLEAVED_REGION_NAMES\",\"type\":\"name\",\"qualifier\":\"full\",\"id\":\"MP_10131411\",\"category\":\"protein\",\"name\":\"Insulin receptor subunit beta\",\"parentId\":null,\"synonyms\":null,\"main\":true,\"synonymId\":\"MP_10131411\",\"composedName\":\"full name\",\"synonymName\":\"Insulin receptor subunit beta\"}],\"additionalNames\":[{\"clazz\":\"ADDITIONAL_NAMES\",\"type\":\"CD antigen\",\"qualifier\":\"CD antigen\",\"id\":\"PR_888982\",\"category\":\"CD antigen\",\"name\":\"CD220\",\"parentId\":null,\"synonyms\":null,\"main\":false,\"synonymId\":\"PR_888982\",\"composedName\":\"CD antigen CD antigen\",\"synonymName\":\"CD220\"}],\"recommendedProteinName\":{\"clazz\":\"PROTEIN_NAMES\",\"type\":\"name\",\"qualifier\":\"full\",\"id\":\"PR_888984\",\"category\":\"protein\",\"name\":\"Insulin receptor\",\"parentId\":null,\"synonyms\":[{\"clazz\":\"PROTEIN_NAMES\",\"type\":\"name\",\"qualifier\":\"short\",\"id\":\"PR_888983\",\"category\":\"protein\",\"name\":\"IR\",\"parentId\":\"PR_888984\",\"synonyms\":null,\"main\":false,\"synonymId\":\"PR_888983\",\"composedName\":\"short name\",\"synonymName\":\"IR\"},{\"clazz\":\"PROTEIN_NAMES\",\"type\":\"enzyme name\",\"qualifier\":\"EC\",\"id\":\"PR_2472057\",\"category\":\"EC\",\"name\":\"2.7.10.1\",\"parentId\":\"PR_888984\",\"synonyms\":null,\"main\":true,\"synonymId\":\"PR_2472057\",\"composedName\":\"EC enzyme name\",\"synonymName\":\"2.7.10.1\"}],\"main\":true,\"synonymId\":\"PR_888984\",\"composedName\":\"full name\",\"synonymName\":\"Insulin receptor\"},\"alternativeProteinNames\":[{\"clazz\":\"ADDITIONAL_NAMES\",\"type\":\"CD antigen\",\"qualifier\":\"CD antigen\",\"id\":\"PR_888982\",\"category\":\"CD antigen\",\"name\":\"CD220\",\"parentId\":null,\"synonyms\":null,\"main\":false,\"synonymId\":\"PR_888982\",\"composedName\":\"CD antigen CD antigen\",\"synonymName\":\"CD220\"}],\"mainProteinName\":\"Insulin receptor\",\"mainGeneName\":\"INSR\",\"proteinExistence\":\"Evidence_at_protein_level\",\"proteinExistenceLevel\":1},\"uniprotName\":\"P06213\",\"proteinExistence\":\"Evidence_at_protein_level\",\"proteinExistenceLevel\":1}{\"properties\":{\"proteinExistence\":\"Evidence at protein level\",\"ptmCount\":3,\"varCount\":41,\"isoformCount\":1,\"mutagenesisCount\":0,\"interactionCount\":4,\"maxSeqLen\":110,\"filterstructure\":true,\"filterdisease\":true,\"filtermutagenesis\":false,\"filterproteomics\":true,\"filterexpressionprofile\":true,\"referencesCount\":592,\"referencesSubmissionsCount\":3,\"referencesAdditionalPublicationsCount\":510,\"referencesCuratedPublicationsCount\":77,\"referencesWebResourcesCount\":3,\"referencesPatentsCount\":0},\"uniqueName\":\"NX_P01308\",\"overview\":{\"history\":{\"proteinExistence\":\"Evidence_at_protein_level\",\"nextprotIntegrationDate\":1267464499678,\"nextprotUpdateDate\":1430703122346,\"uniprotIntegrationDate\":522324000000,\"uniprotUpdateDate\":1430301600000,\"uniprotVersion\":\"203\",\"lastSequenceUpdate\":\"1986-07-21\",\"sequenceVersion\":\"1\",\"formattedUniprotIntegrationDate\":\"1986-07-21\",\"formattedUniprotUpdateDate\":\"2015-04-29\",\"proteinExistenceRaw\":\"protein level\",\"formattedNextprotIntegrationDate\":\"2010-03-01\",\"formattedNextprotUpdateDate\":\"2015-05-04\",\"proteinExistenceLevel\":1},\"families\":[{\"familyId\":66231,\"accession\":\"FA-01869\",\"name\":\"Insulin\",\"level\":\"Family\",\"description\":\"Belongs to the insulin family.\",\"region\":null,\"parent\":null}],\"proteinNames\":[{\"clazz\":\"PROTEIN_NAMES\",\"type\":\"name\",\"qualifier\":\"full\",\"id\":\"PR_720135\",\"category\":\"protein\",\"name\":\"Insulin\",\"parentId\":null,\"synonyms\":null,\"main\":true,\"synonymId\":\"PR_720135\",\"composedName\":\"full name\",\"synonymName\":\"Insulin\"}],\"geneNames\":[{\"clazz\":\"GENE_NAMES\",\"type\":\"gene name\",\"qualifier\":null,\"id\":\"PR_1171036\",\"category\":\"gene name\",\"name\":\"INS\",\"parentId\":null,\"synonyms\":null,\"main\":true,\"synonymId\":\"PR_1171036\",\"composedName\":\"gene name\",\"synonymName\":\"INS\"}],\"functionalRegionNames\":null,\"cleavedRegionNames\":[{\"clazz\":\"CLEAVED_REGION_NAMES\",\"type\":\"name\",\"qualifier\":\"full\",\"id\":\"MP_10132975\",\"category\":\"protein\",\"name\":\"Insulin B chain\",\"parentId\":null,\"synonyms\":null,\"main\":true,\"synonymId\":\"MP_10132975\",\"composedName\":\"full name\",\"synonymName\":\"Insulin B chain\"},{\"clazz\":\"CLEAVED_REGION_NAMES\",\"type\":\"name\",\"qualifier\":\"full\",\"id\":\"MP_10132976\",\"category\":\"protein\",\"name\":\"Insulin A chain\",\"parentId\":null,\"synonyms\":null,\"main\":true,\"synonymId\":\"MP_10132976\",\"composedName\":\"full name\",\"synonymName\":\"Insulin A chain\"}],\"additionalNames\":null,\"recommendedProteinName\":{\"clazz\":\"PROTEIN_NAMES\",\"type\":\"name\",\"qualifier\":\"full\",\"id\":\"PR_720135\",\"category\":\"protein\",\"name\":\"Insulin\",\"parentId\":null,\"synonyms\":null,\"main\":true,\"synonymId\":\"PR_720135\",\"composedName\":\"full name\",\"synonymName\":\"Insulin\"},\"alternativeProteinNames\":[],\"mainProteinName\":\"Insulin\",\"mainGeneName\":\"INS\",\"proteinExistence\":\"Evidence_at_protein_level\",\"proteinExistenceLevel\":1},\"uniprotName\":\"P01308\",\"proteinExistence\":\"Evidence_at_protein_level\",\"proteinExistenceLevel\":1}", out.toString(), true);
    }

    @Test
    public void testFastaExportStream() throws Exception {

        ByteArrayOutputStream out = new ByteArrayOutputStream();
        Writer writer = new PrintWriter(out);

        NPStreamExporter exporter = new FastaStreamExporter();

        exporter.export(Arrays.asList("NX_P06213", "NX_P01308"), writer, "overview", null);

        assertEquals(">nxp|NX_P06213-2|INSR|Insulin receptor|Short\n" +
                "MATGGRRGAAAAPLLVAVAALLLGAAGHLYPGEVCPGMDIRNNLTRLHELENCSVIEGHL\n" +
                "QILLMFKTRPEDFRDLSFPKLIMITDYLLLFRVYGLESLKDLFPNLTVIRGSRLFFNYAL\n" +
                "VIFEMVHLKELGLYNLMNITRGSVRIEKNNELCYLATIDWSRILDSVEDNYIVLNKDDNE\n" +
                "ECGDICPGTAKGKTNCPATVINGQFVERCWTHSHCQKVCPTICKSHGCTAEGLCCHSECL\n" +
                "GNCSQPDDPTKCVACRNFYLDGRCVETCPPPYYHFQDWRCVNFSFCQDLHHKCKNSRRQG\n" +
                "CHQYVIHNNKCIPECPSGYTMNSSNLLCTPCLGPCPKVCHLLEGEKTIDSVTSAQELRGC\n" +
                "TVINGSLIINIRGGNNLAAELEANLGLIEEISGYLKIRRSYALVSLSFFRKLRLIRGETL\n" +
                "EIGNYSFYALDNQNLRQLWDWSKHNLTITQGKLFFHYNPKLCLSEIHKMEEVSGTKGRQE\n" +
                "RNDIALKTNGDQASCENELLKFSYIRTSFDKILLRWEPYWPPDFRDLLGFMLFYKEAPYQ\n" +
                "NVTEFDGQDACGSNSWTVVDIDPPLRSNDPKSQNHPGWLMRGLKPWTQYAIFVKTLVTFS\n" +
                "DERRTYGAKSDIIYVQTDATNPSVPLDPISVSNSSSQIILKWKPPSDPNGNITHYLVFWE\n" +
                "RQAEDSELFELDYCLKGLKLPSRTWSPPFESEDSQKHNQSEYEDSAGECCSCPKTDSQIL\n" +
                "KELEESSFRKTFEDYLHNVVFVPRPSRKRRSLGDVGNVTVAVPTVAAFPNTSSTSVPTSP\n" +
                "EEHRPFEKVVNKESLVISGLRHFTGYRIELQACNQDTPEERCSVAAYVSARTMPEAKADD\n" +
                "IVGPVTHEIFENNVVHLMWQEPKEPNGLIVLYEVSYRRYGDEELHLCVSRKHFALERGCR\n" +
                "LRGLSPGNYSVRIRATSLAGNGSWTEPTYFYVTDYLDVPSNIAKIIIGPLIFVFLFSVVI\n" +
                "GSIYLFLRKRQPDGPLGPLYASSNPEYLSASDVFPCSVYVPDEWEVSREKITLLRELGQG\n" +
                "SFGMVYEGNARDIIKGEAETRVAVKTVNESASLRERIEFLNEASVMKGFTCHHVVRLLGV\n" +
                "VSKGQPTLVVMELMAHGDLKSYLRSLRPEAENNPGRPPPTLQEMIQMAAEIADGMAYLNA\n" +
                "KKFVHRDLAARNCMVAHDFTVKIGDFGMTRDIYETDYYRKGGKGLLPVRWMAPESLKDGV\n" +
                "FTTSSDMWSFGVVLWEITSLAEQPYQGLSNEQVLKFVMDGGYLDQPDNCPERVTDLMRMC\n" +
                "WQFNPKMRPTFLEIVNLLKDDLHPSFPEVSFFHSEENKAPESEELEMEFEDMENVPLDRS\n" +
                "SHCQREEAGGRDGGSSLGFKRSYEEHIPYTHMNGGKKNGRILTLPRSNPS\n" +
                ">nxp|NX_P06213-1|INSR|Insulin receptor|Long\n" +
                "MATGGRRGAAAAPLLVAVAALLLGAAGHLYPGEVCPGMDIRNNLTRLHELENCSVIEGHL\n" +
                "QILLMFKTRPEDFRDLSFPKLIMITDYLLLFRVYGLESLKDLFPNLTVIRGSRLFFNYAL\n" +
                "VIFEMVHLKELGLYNLMNITRGSVRIEKNNELCYLATIDWSRILDSVEDNYIVLNKDDNE\n" +
                "ECGDICPGTAKGKTNCPATVINGQFVERCWTHSHCQKVCPTICKSHGCTAEGLCCHSECL\n" +
                "GNCSQPDDPTKCVACRNFYLDGRCVETCPPPYYHFQDWRCVNFSFCQDLHHKCKNSRRQG\n" +
                "CHQYVIHNNKCIPECPSGYTMNSSNLLCTPCLGPCPKVCHLLEGEKTIDSVTSAQELRGC\n" +
                "TVINGSLIINIRGGNNLAAELEANLGLIEEISGYLKIRRSYALVSLSFFRKLRLIRGETL\n" +
                "EIGNYSFYALDNQNLRQLWDWSKHNLTITQGKLFFHYNPKLCLSEIHKMEEVSGTKGRQE\n" +
                "RNDIALKTNGDQASCENELLKFSYIRTSFDKILLRWEPYWPPDFRDLLGFMLFYKEAPYQ\n" +
                "NVTEFDGQDACGSNSWTVVDIDPPLRSNDPKSQNHPGWLMRGLKPWTQYAIFVKTLVTFS\n" +
                "DERRTYGAKSDIIYVQTDATNPSVPLDPISVSNSSSQIILKWKPPSDPNGNITHYLVFWE\n" +
                "RQAEDSELFELDYCLKGLKLPSRTWSPPFESEDSQKHNQSEYEDSAGECCSCPKTDSQIL\n" +
                "KELEESSFRKTFEDYLHNVVFVPRKTSSGTGAEDPRPSRKRRSLGDVGNVTVAVPTVAAF\n" +
                "PNTSSTSVPTSPEEHRPFEKVVNKESLVISGLRHFTGYRIELQACNQDTPEERCSVAAYV\n" +
                "SARTMPEAKADDIVGPVTHEIFENNVVHLMWQEPKEPNGLIVLYEVSYRRYGDEELHLCV\n" +
                "SRKHFALERGCRLRGLSPGNYSVRIRATSLAGNGSWTEPTYFYVTDYLDVPSNIAKIIIG\n" +
                "PLIFVFLFSVVIGSIYLFLRKRQPDGPLGPLYASSNPEYLSASDVFPCSVYVPDEWEVSR\n" +
                "EKITLLRELGQGSFGMVYEGNARDIIKGEAETRVAVKTVNESASLRERIEFLNEASVMKG\n" +
                "FTCHHVVRLLGVVSKGQPTLVVMELMAHGDLKSYLRSLRPEAENNPGRPPPTLQEMIQMA\n" +
                "AEIADGMAYLNAKKFVHRDLAARNCMVAHDFTVKIGDFGMTRDIYETDYYRKGGKGLLPV\n" +
                "RWMAPESLKDGVFTTSSDMWSFGVVLWEITSLAEQPYQGLSNEQVLKFVMDGGYLDQPDN\n" +
                "CPERVTDLMRMCWQFNPKMRPTFLEIVNLLKDDLHPSFPEVSFFHSEENKAPESEELEME\n" +
                "FEDMENVPLDRSSHCQREEAGGRDGGSSLGFKRSYEEHIPYTHMNGGKKNGRILTLPRSN\n" +
                "PS\n" +
                ">nxp|NX_P01308-1|INS|Insulin|Iso 1\n" +
                "MALWMRLLPLLALLALWGPDPAAAFVNQHLCGSHLVEALYLVCGERGFFYTPKTRREAED\n" +
                "LQVGQVELGGGPGAGSLQPLALEGSLQKRGIVEQCCTSICSLYQLENYCN\n", out.toString());
    }

    @Ignore
    @Test
    public void testPeffExportStream() throws Exception {

        ByteArrayOutputStream out = new ByteArrayOutputStream();

        Writer writer = new PrintWriter(out);

        NPStreamExporter exporter = new PeffStreamExporter();

        exporter.export(Arrays.asList("NX_P06213", "NX_P01308"), writer, "overview", null);

        assertEquals("#nb entries=2\nNX_P06213\nNX_P01308\n", out.toString());
    }

    @Test
    public void testTXTExportStream() throws Exception {

        ByteArrayOutputStream out = new ByteArrayOutputStream();

        Writer writer = new PrintWriter(out);

        NPStreamExporter exporter = new TXTStreamExporter();

        Map<String, Object> params = new HashMap<>();
        params.put(ExportService.ENTRIES_COUNT_PARAM, 2);

        exporter.export(Arrays.asList("NX_P06213", "NX_P01308"), writer, "overview", params);

        assertEquals("#nb entries=2\nNX_P06213\nNX_P01308\n", out.toString());
    }
}