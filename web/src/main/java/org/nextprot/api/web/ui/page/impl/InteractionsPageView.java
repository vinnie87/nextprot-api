package org.nextprot.api.web.ui.page.impl;

import org.nextprot.api.commons.constants.AnnotationCategory;
import org.nextprot.api.web.ui.page.EntryPage;

import javax.annotation.Nonnull;
import java.util.Arrays;
import java.util.List;

/**
 * Please keep this class in sync with specs in https://swissprot.isb-sib.ch/wiki/display/cal/neXtProt+Interactions+view+specs
 * @author pmichel
 *
 */
public class InteractionsPageView extends PageViewBase {

	InteractionsPageView() {
		super(EntryPage.INTERACTIONS);
	}

	@Nonnull
	@Override
	protected List<AnnotationCategory> getAnnotationCategoryWhiteList() {
		return Arrays.asList(
				
				AnnotationCategory.INTERACTION_INFO, 
				AnnotationCategory.ENZYME_REGULATION,
				AnnotationCategory.GO_MOLECULAR_FUNCTION, // further refined
				AnnotationCategory.BINARY_INTERACTION,
				AnnotationCategory.COFACTOR,
				AnnotationCategory.COFACTOR_INFO,
				AnnotationCategory.SMALL_MOLECULE_INTERACTION,
				AnnotationCategory.MISCELLANEOUS
		);
	}

	@Nonnull
	@Override
	protected List<AnnotationCategory> getFeatureCategoryWhiteList() {
		return Arrays.asList(
				AnnotationCategory.MISCELLANEOUS_REGION,
				AnnotationCategory.CALCIUM_BINDING_REGION,
				AnnotationCategory.DNA_BINDING_REGION,
				AnnotationCategory.NUCLEOTIDE_PHOSPHATE_BINDING_REGION,
				AnnotationCategory.INTERACTING_REGION,
				AnnotationCategory.BINDING_SITE,
				AnnotationCategory.METAL_BINDING_SITE
			);
	}

	@Nonnull
	@Override
	protected List<String> getXrefDbNameWhiteList() {
		return Arrays.asList(
				"BindingDB","DIP","IntAct","MINT",
				"STRING", "SignaLink", "BioGrid","SIGNOR"
			);
	}
}
