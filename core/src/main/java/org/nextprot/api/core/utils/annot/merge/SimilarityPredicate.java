package org.nextprot.api.core.utils.annot.merge;

import com.google.common.base.Preconditions;
import org.nextprot.api.commons.constants.AnnotationCategory;
import org.nextprot.api.core.domain.annotation.Annotation;
import org.nextprot.api.core.utils.annot.merge.impl.ObjectSimilarityPredicate;
import org.nextprot.api.core.utils.annot.merge.impl.SimilarityPredicateChain;
import org.nextprot.api.core.utils.annot.merge.impl.VariantPositionMatcher;

import java.util.Arrays;

/**
 * Defines the contract to evaluate similarity of 2 annotations
 *
 * Created by fnikitin on 02/08/16.
 */
public interface SimilarityPredicate {

    /**
     * @return true if annotations are similar else false
     */
    boolean isSimilar(Annotation annotation1, Annotation annotation2);

    /**
     * Static factory method that return a default predicate specific of the given category
     * @param category the annotation category to estimate similarity
     */
    static SimilarityPredicate newSimilarityPredicate(AnnotationCategory category) {

        Preconditions.checkNotNull(category);

        switch (category) {
            case GO_BIOLOGICAL_PROCESS:
            case GO_CELLULAR_COMPONENT:
            case GO_MOLECULAR_FUNCTION:
                // TODO: there is a 'is_negative' field to consider in the future
                return new ObjectSimilarityPredicate<>(Annotation::getCvTermAccessionCode);
            case VARIANT:
            case MUTAGENESIS:
                return new SimilarityPredicateChain(Arrays.asList(
                        new ObjectSimilarityPredicate<>(Annotation::getVariant,
                                (v1, v2) -> v1.getOriginal().equals(v2.getOriginal()) && v1.getVariant().equals(v2.getVariant())),
                        new ObjectSimilarityPredicate<>(Annotation::getTargetingIsoformsMap,
                                new VariantPositionMatcher()
                        )
                ));
            case BINARY_INTERACTION:
            case SMALL_MOLECULE_INTERACTION:
                return new ObjectSimilarityPredicate<>(Annotation::getBioObject,
                        (bo1, bo2) -> bo1.getAccession().equals(bo2.getAccession()) && bo1.getDatabase().equalsIgnoreCase(bo2.getDatabase())
                );
            default:
                return null;
        }
    }
}