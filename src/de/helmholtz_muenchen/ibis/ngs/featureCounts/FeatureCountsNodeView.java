package de.helmholtz_muenchen.ibis.ngs.featureCounts;

import de.helmholtz_muenchen.ibis.utils.abstractNodes.ExecutorNode.ExecutorNodeView;

/**
 * <code>NodeView</code> for the "FeatureCounts" Node.
 * featureCounts: an efficient general purpose program for assigning sequence reads to genomic features
 * For details, please see paper: 
 * Liao et al, Bioinformatics 2013; doi: 10.1093/bioinformatics/btt656
 *
 * @author Michael Kluge
 */
public class FeatureCountsNodeView extends ExecutorNodeView<FeatureCountsNodeModel> {

    /**
     * Creates a new view.
     * 
     * @param nodeModel The model (class: {@link FeatureCountsNodeModel})
     */
    protected FeatureCountsNodeView(final FeatureCountsNodeModel nodeModel) {
        super(nodeModel);
    }
}

