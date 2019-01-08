package org.nextprot.api.core.controller;

import org.expasy.mzjava.proteomics.mol.digest.Protease;
import org.jsondoc.core.annotation.Api;
import org.jsondoc.core.annotation.ApiMethod;
import org.jsondoc.core.annotation.ApiQueryParam;
import org.jsondoc.core.pojo.ApiVerb;
import org.nextprot.api.commons.exception.NextProtException;
import org.nextprot.api.core.service.DigestionService;
import org.springframework.beans.factory.annotation.Autowired;
import org.springframework.http.MediaType;
import org.springframework.stereotype.Controller;
import org.springframework.web.bind.annotation.RequestMapping;
import org.springframework.web.bind.annotation.RequestMethod;
import org.springframework.web.bind.annotation.RequestParam;
import org.springframework.web.bind.annotation.ResponseBody;

import java.util.Set;

// See also sources of mzjava-proteomics are available at https://bitbucket.org/sib-pig/mzjava-proteomics
@Controller
@Api(name = "Protein digestion", description = "Digest proteins with proteases", group = "Tools")
public class ProteinDigestionController {

	@Autowired
	private DigestionService digestionService;

	@ResponseBody
	@RequestMapping(value = "/digestion/available-protease-list", method = { RequestMethod.GET }, produces = { MediaType.APPLICATION_JSON_VALUE })
	@ApiMethod(path = "/digestion/available-protease-list", verb = ApiVerb.GET, description = "list all available proteases")
	public Protease[] listAllProteases() {

		return digestionService.getProteases();
	}

	@ResponseBody
	@RequestMapping(value = "/digestion/digest-all-proteins", method = { RequestMethod.GET }, produces = { MediaType.APPLICATION_JSON_VALUE })
	@ApiMethod(path = "/digestion/digest-all-proteins", verb = ApiVerb.GET, description = "digest all neXtProt mature proteins with TRYPSIN (with number of missed cleavage max of 2)")
	public Set<String> digestAllProteins() {

		return digestionService.digestAllWithTrypsin();
	}

	@ResponseBody
	@RequestMapping(value = "/digestion/{entryAccession}", method = { RequestMethod.GET }, produces = { MediaType.APPLICATION_JSON_VALUE })
	@ApiMethod(path = "/digestion/{entryAccession}", verb = ApiVerb.GET, description = "digest mature parts of a protein with a specific protease")
	public Set<String> digestProtein(
			@ApiQueryParam(name = "entryAccession", description = "A neXtProt entry accession.", allowedvalues = { "NX_P01308" })
			@RequestParam(value = "entryAccession") String entryAccession,
			@ApiQueryParam(name = "protease", description = "chose a protease to digest a protein", allowedvalues = { "TRYPSIN" })
			@RequestParam(value = "protease") String protease,
			@ApiQueryParam(name = "minpeplen", description = "minimum peptide length", allowedvalues = { "7" })
			@RequestParam(value = "minpeplen", required = false) Integer minPepLen,
			@ApiQueryParam(name = "maxpeplen", description = "maximum peptide length", allowedvalues = { "77" })
			@RequestParam(value = "maxpeplen", required = false) Integer maxPepLen,
			@ApiQueryParam(name = "maxmissedcleavages", description = "maximum number of missed cleavages (cannot be greater than 2)", allowedvalues = { "2" })
			@RequestParam(value = "maxmissedcleavages", required = false) Integer maxMissedCleavages) {

		if (maxMissedCleavages > 2) {

			throw new NextProtException("Cannot configure digestion with such a high number of missed cleavages ("+maxMissedCleavages+")");
		}

		return digestionService.digest(entryAccession, Protease.valueOf(protease.toUpperCase()), minPepLen, maxPepLen, maxMissedCleavages);
	}
}
