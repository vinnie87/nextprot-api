package org.nextprot.api.core.domain;

import org.junit.Assert;
import org.junit.Test;

import java.util.Collections;

import static org.nextprot.api.core.utils.dbxref.DbXrefURLResolverTest.createDbXref;
import static org.nextprot.api.core.utils.dbxref.DbXrefURLResolverTest.createDbXrefProperty;

@Deprecated
public class DbXrefTest {

    @Test
    public void testResolveWEBINFO() throws Exception {

        DbXref xref = createDbXref("babebibobu", "WEBINFO", "");

        Assert.assertEquals("babebibobu", xref.resolveLinkTarget());
    }

    @Test
    public void testMissingHttpProtocolTemplate() throws Exception {

        DbXref xref = createDbXref("babebibobu", "Ensembl", "www.ensembl.org/id/%s");

        Assert.assertEquals("http://www.ensembl.org/id/babebibobu", xref.resolveLinkTarget());
    }

    @Test
    public void testUnknownDbNameAndEmptyURL() throws Exception {

        DbXref xref = createDbXref("babebibobu", "unknownDb", "");

        Assert.assertEquals("", xref.resolveLinkTarget());
    }

    @Test
    public void testResolveEmbl() throws Exception {

        DbXref xref = createDbXref("AF009225", "EMBL", "http://www.ebi.ac.uk/ena/data/view/%s");

        Assert.assertEquals("http://www.ebi.ac.uk/ena/data/view/AF009225", xref.resolveLinkTarget());
    }

    @Test
    public void testResolveEmblWithDotAccession() throws Exception {

        DbXref xref = createDbXref("CAH72401.1", "EMBL", "whatever");

        Assert.assertEquals("http://www.ebi.ac.uk/cgi-bin/dbfetch?db=emblcds&id=CAH72401", xref.resolveLinkTarget());
    }

    // entry/NX_Q9BXA6/xref.json
    @Test
    public void testResolveEnsemblENSG() throws Exception {

        DbXref xref = createDbXref("ENSG00000178093", "Ensembl", "whatever");

        Assert.assertEquals("http://www.ensembl.org/Homo_sapiens/Gene/Summary?db=core;g=ENSG00000178093", xref.resolveLinkTarget());
    }

    // entry/NX_Q9BXA6/xref.json
    @Test
    public void testResolveEnsemblENSP() throws Exception {

        DbXref xref = createDbXref("ENSP00000466056", "Ensembl", "whatever");

        Assert.assertEquals("http://www.ensembl.org/Homo_sapiens/Transcript/ProteinSummary?db=core;p=ENSP00000466056", xref.resolveLinkTarget());
    }

    // entry/NX_Q9BXA6/xref.json
    @Test
    public void testResolveEnsemblENST() throws Exception {

        DbXref xref = createDbXref("ENST00000587522", "Ensembl", "whatever");

        Assert.assertEquals("http://www.ensembl.org/Homo_sapiens/Transcript/Summary?db=core;t=ENST00000587522", xref.resolveLinkTarget());
    }

    @Test
    public void testResolveEnsemblBadPrimaryId() throws Exception {

        DbXref xref = createDbXref("EST00000587522", "Ensembl", "whatever");

        Assert.assertNull(xref.resolveLinkTarget());
    }

    // entry/NX_Q9BXA6/xref.json
    @Test
    public void testResolveCosmicCOSM() throws Exception {

        DbXref xref = createDbXref("COSM1172604", "Cosmic", "whatever");

        Assert.assertEquals("http://cancer.sanger.ac.uk/cosmic/mutation/overview?id=1172604", xref.resolveLinkTarget());
    }

    // entry/NX_?????/xref.json
    @Test
    public void testResolveCosmicCOSS() throws Exception {

        DbXref xref = createDbXref("COSS1172604", "Cosmic", "whatever");

        Assert.assertEquals("http://cancer.sanger.ac.uk/cosmic/sample/overview?id=1172604", xref.resolveLinkTarget());
    }

    @Test
    public void testResolveCosmicBadPrimaryId() throws Exception {

        DbXref xref = createDbXref("HCFC1", "Cosmic", "whatever");

        Assert.assertEquals("http://cancer.sanger.ac.uk/cosmic/gene/overview?ln=HCFC1", xref.resolveLinkTarget());
    }

    @Test
    public void testResolveClinvar() throws Exception {

        DbXref xref = createDbXref("HCFC1", "Clinvar", "url.whatever");

        Assert.assertEquals("http://www.ncbi.nlm.nih.gov/clinvar/?term=HCFC1", xref.resolveLinkTarget());
    }

    // entry/NX_P51610/xref.json
    @Test
    public void testResolvePIR() throws Exception {

        DbXref xref = createDbXref("babebibobu", "PIR", "http://pir.georgetown.edu/cgi-bin/nbrfget?uid=%s");

        xref.setProperties(Collections.singletonList(createDbXrefProperty("entry name", "A40718")));

        Assert.assertEquals("http://pir.georgetown.edu/cgi-bin/nbrfget?uid=A40718", xref.resolveLinkTarget());
    }

    // entry/NX_P51610/xref.json
    @Test
    public void testResolvePIRWithoutURLProtocol() throws Exception {

        DbXref xref = createDbXref("babebibobu", "PIR", "pir.georgetown.edu/cgi-bin/nbrfget?uid=%s");

        xref.setProperties(Collections.singletonList(createDbXrefProperty("entry name", "A40718")));

        Assert.assertEquals("http://pir.georgetown.edu/cgi-bin/nbrfget?uid=A40718", xref.resolveLinkTarget());
    }

    // unexpectedly throw NPE !!
    //@Test
    public void testResolvePIRShouldNotThrowsNPE() throws Exception {

        DbXref xref = createDbXref("babebibobu", "PIR", "pir.georgetown.edu/cgi-bin/nbrfget?uid=%s");

        Assert.assertEquals("http://pir.georgetown.edu/cgi-bin/nbrfget?uid=A40718", xref.resolveLinkTarget());
    }

    // entry/NX_Q9BXA6/xref.json
    @Test
    public void testResolveGermOnline() throws Exception {

        DbXref xref = createDbXref("ENSG00000178093", "GermOnline", "whatever");

        Assert.assertEquals("http://www.germonline.org/Homo_sapiens/geneview?gene=ENSG00000178093", xref.resolveLinkTarget());
    }

    @Test
    public void testResolveGenevestigator() throws Exception {

        DbXref xref = createDbXref("P01308", "Genevestigator", "whatever");

        Assert.assertEquals("http://genevisible.com/tissues/HS/UniProt/P01308", xref.resolveLinkTarget());
    }

    // entry/NX_Q9BXA6/xref.json
    @Test
    public void testResolveHPAGene() throws Exception {

        DbXref xref = createDbXref("ENSG00000254647", "HPA", "whatever");

        Assert.assertEquals("http://www.proteinatlas.org/ENSG00000254647", xref.resolveLinkTarget());
    }

    // entry/NX_P51610/xref.json
    @Test
    public void testResolveHPASubcellular() throws Exception {

        DbXref xref = createDbXref("ENSG00000254647/subcellular", "HPA", "whatever");

        Assert.assertEquals("http://www.proteinatlas.org/ENSG00000254647/subcellular", xref.resolveLinkTarget());
    }

    // entry/NX_P51610/xref.json
    @Test
    public void testResolveHPAAntibody() throws Exception {

        DbXref xref = createDbXref("HPA018312", "HPA", "whatever");

        Assert.assertEquals("http://www.proteinatlas.org/search/HPA018312", xref.resolveLinkTarget());
    }

    // entry/NX_P51610/xref.json
    @Test
    public void testResolveGenevisible() throws Exception {

        DbXref xref = createDbXref("P51610", "Genevisible", "http://genevisible.com/tissues/%s2/UniProt/%s1");

        Assert.assertEquals("http://genevisible.com/tissues/HS/UniProt/P51610", xref.resolveLinkTarget());
    }

    @Test
    public void testResolveUniGene() throws Exception {

        DbXref xref = createDbXref("Hs.83634", "UniGene", "http://www.ncbi.nlm.nih.gov/UniGene/clust.cgi?ORG=%s1&CID=%s2");

        Assert.assertEquals("http://www.ncbi.nlm.nih.gov/UniGene/clust.cgi?ORG=Hs&CID=83634", xref.resolveLinkTarget());
    }

    @Test
    public void testResolveUCSC() throws Exception {

        DbXref xref = createDbXref("uc004fjp.3", "UCSC", "http://genome.ucsc.edu/cgi-bin/hgGene?hgg_gene=%s1&org=%s2");

        Assert.assertEquals("http://genome.ucsc.edu/cgi-bin/hgGene?hgg_gene=uc004fjp.3&org=human", xref.resolveLinkTarget());
    }

    @Test
    public void testResolveIntAct() throws Exception {

        DbXref xref = createDbXref("EBI-1644164,EBI-396176", "IntAct", "whatever");

        Assert.assertEquals("http://www.ebi.ac.uk/intact/search/do/search?binary=EBI-1644164,EBI-396176", xref.resolveLinkTarget());
    }

    @Test
    public void testResolveIntActNoEBIId() throws Exception {

        DbXref xref = createDbXref("P51610", "IntAct", "http://www.ebi.ac.uk/intact/pages/interactions/interactions.xhtml?query=%s");

        Assert.assertEquals("http://www.ebi.ac.uk/intact/pages/interactions/interactions.xhtml?query=P51610", xref.resolveLinkTarget());
    }

    @Test
    public void testResolveProsite() throws Exception {

        DbXref xref = createDbXref("PS50853", "PROSITE", "whatever");

        Assert.assertEquals("http://prosite.expasy.org/cgi-bin/prosite/prosite-search-ac?PS50853", xref.resolveLinkTarget());
    }

    @Test
    public void testResolveHSSP() throws Exception {

        DbXref xref = createDbXref("PS50853", "HSSP", "%s");

        Assert.assertEquals("http://ps50853", xref.resolveLinkTarget());
    }

    @Test
    public void testResolveHSSPWithPDB() throws Exception {

        DbXref xref = createDbXref("PS50853", "HSSP", "%s");
        xref.setProperties(Collections.singletonList(createDbXrefProperty("PDB accession", "1A7F")));

        Assert.assertEquals("http://1a7f", xref.resolveLinkTarget());
    }

    // entry/P51610/xref.json
    @Test
    public void testResolveBgeeENSG() throws Exception {

        DbXref xref = createDbXref("ENSG00000172534", "Bgee", "http://bgee.unil.ch/bgee/bgee?uniprot_id=%s");

        Assert.assertEquals("http://bgee.unil.ch/bgee/bgee?page=expression&action=data&ENSG00000172534", xref.resolveLinkTarget());
    }

    // entry/P51610/xref.json
    @Test
    public void testResolveBgeeNoENSG() throws Exception {

        DbXref xref = createDbXref("P51610", "Bgee", "http://bgee.unil.ch/bgee/bgee?uniprot_id=%s");

        Assert.assertEquals("http://bgee.unil.ch/bgee/bgee?uniprot_id=P51610", xref.resolveLinkTarget());
    }

    // entry/NX_P01308/xref.json
    @Test
    public void testResolvePeptideAtlasPap() throws Exception {

        DbXref xref = createDbXref("PAp00001490", "PeptideAtlas", "whatever");

        Assert.assertEquals("https://db.systemsbiology.net/sbeams/cgi/PeptideAtlas/GetPeptide?searchWithinThis=Peptide+Name&searchForThis=PAp00001490;organism_name=Human", xref.resolveLinkTarget());
    }

    // entry/NX_P01308/xref.json
    @Test
    public void testResolvePeptideAtlasNoPap() throws Exception {

        DbXref xref = createDbXref("P01308", "PeptideAtlas", "whatever");

        Assert.assertEquals("https://db.systemsbiology.net/sbeams/cgi/PeptideAtlas/GetProtein?protein_name=P01308;organism_name=Human;action=GO", xref.resolveLinkTarget());
    }

    // entry/NX_P01308/xref.json
    @Test
    public void testResolveSRMAtlas() throws Exception {

        DbXref xref = createDbXref("PAp00968082", "SRMAtlas", "https://db.systemsbiology.net/sbeams/cgi/PeptideAtlas/GetTransitions?organism_name=Human;default_search=1;peptide_sequence_constraint=%s;apply_action=QUERY");
        xref.setProperties(Collections.singletonList(createDbXrefProperty("sequence", "GFFYTPK")));

        Assert.assertEquals("https://db.systemsbiology.net/sbeams/cgi/PeptideAtlas/GetTransitions?organism_name=Human;default_search=1;peptide_sequence_constraint=GFFYTPK;apply_action=QUERY", xref.resolveLinkTarget());
    }

    // entry/NX_P01308/xref.json
    @Test
    public void testResolvePDB() throws Exception {

        DbXref xref = createDbXref("1A7F", "PDB", "whatever");

        Assert.assertEquals("http://www.pdb.org/pdb/explore/explore.do?pdbId=1A7F", xref.resolveLinkTarget());
    }

    @Test
    public void testResolveTKG() throws Exception {

        DbXref xref = createDbXref("0377", "TKG", "http://www2.idac.tohoku.ac.jp/dep/ccr/TKGdate/TKGvo10%n/%s.html");

        Assert.assertEquals("http://www2.idac.tohoku.ac.jp/dep/ccr/TKGdate/TKGvo103/0377.html", xref.resolveLinkTarget());
    }

    @Test
    public void testResolveNIH_ARP() throws Exception {

        DbXref xref = createDbXref("11411-223", "NIH-ARP", "https://www.aidsreagent.org/reagentdetail.cfm?t=cell_lines&id=%s");

        Assert.assertEquals("https://www.aidsreagent.org/reagentdetail.cfm?t=cell_lines&id=223", xref.resolveLinkTarget());
    }

    @Test
    public void testResolveCGH_DB() throws Exception {

        DbXref xref = createDbXref("9029-4", "CGH-DB", "http://www.cghtmd.jp/CGHDatabase/mapViewer?hid=%s&aid=%t&lang=en");

        Assert.assertEquals("http://www.cghtmd.jp/CGHDatabase/mapViewer?hid=9029&aid=4&lang=en", xref.resolveLinkTarget());
    }

    @Test
    public void testResolveIFOMock() throws Exception {

        DbXref xref = createDbXref("IFO1234", "IFO", "whatever/%s");
        Assert.assertEquals("http://whatever/ifo1234", xref.resolveLinkTarget());
    }

    @Test
    public void testResolveJCRBMock() throws Exception {

        DbXref xref = createDbXref("JCRB1234", "JCRB", "whatever/%s");
        Assert.assertEquals("http://whatever/jcrb1234", xref.resolveLinkTarget());
    }

    @Test
    public void testResolveIFO() throws Exception {

        DbXref xref = createDbXref("IFO1234", "IFO", "http://cellbank.nibio.go.jp/~cellbank/cgi-bin/search_res_det.cgi?RNO=%s");
        Assert.assertEquals("http://cellbank.nibio.go.jp/~cellbank/cgi-bin/search_res_det.cgi?RNO=ifo1234", xref.resolveLinkTarget());
    }

    @Test
    public void testResolveJCRB() throws Exception {

        DbXref xref = createDbXref("JCRB1234", "JCRB", "http://cellbank.nibio.go.jp/~cellbank/en/search_res_list.cgi?KEYWOD=%s");
        Assert.assertEquals("http://cellbank.nibio.go.jp/~cellbank/en/search_res_list.cgi?KEYWOD=jcrb1234", xref.resolveLinkTarget());
    }

    // unexpectedly returns null !!
    // entry/NX_Q9BXA6/xref.json
    //@Test
    public void testResolvePRO() throws Exception {

        DbXref xref = createDbXref("PR:000028527", "PRO", "http://purl.obolibrary.org/obo/PR_%u");

        // cannot resolve %u here !!
        Assert.assertEquals("http://purl.obolibrary.org/obo/PR_000028527", xref.resolveLinkTarget());
    }

    @Test
    public void testResolveCLO() throws Exception {

        DbXref xref = createDbXref("CLO:0000031", "CLO", "purl.obolibrary.org/obo/%s");
        Assert.assertEquals("http://purl.obolibrary.org/obo/CLO_0000031", xref.resolveLinkTarget());
    }

    @Test
    public void testResolveFMA() throws Exception {

        DbXref xref = createDbXref("FMA:62955", "FMA", "http://purl.obolibrary.org/obo/%s");
        Assert.assertEquals("http://purl.obolibrary.org/obo/FMA_62955", xref.resolveLinkTarget());
    }

    @Test
    public void testResolveBrenda() throws Exception {

        DbXref xref = createDbXref("2.7.11.21", "BRENDA", "http://www.brenda-enzymes.org/enzyme.php?ecno=%s");
        Assert.assertEquals("http://www.brenda-enzymes.org/enzyme.php?ecno=2.7.11.21", xref.resolveLinkTarget());
    }

    @Test
    public void testResolveHogenom() throws Exception {

        DbXref xref = createDbXref("HOG000007899", "HOGENOM", "http://pbil.univ-lyon1.fr/cgi-bin/acnuc-ac2tree?query=%u&db=HOGENOM");
        Assert.assertEquals("http://pbil.univ-lyon1.fr/cgi-bin/acnuc-ac2tree?query=Q8NBS9&db=HOGENOM", DbXref.resolvePercentULinkTarget("NX_Q8NBS9", xref));
    }
}