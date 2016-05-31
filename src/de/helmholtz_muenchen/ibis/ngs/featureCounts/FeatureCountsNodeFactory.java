package de.helmholtz_muenchen.ibis.ngs.featureCounts;

import org.knime.core.node.NodeDialogPane;
import org.knime.core.node.NodeFactory;
import org.knime.core.node.NodeView;

import de.helmholtz_muenchen.ibis.utils.abstractNodes.HTExecutorNode.HTENodeView;

/**
 * <code>NodeFactory</code> for the "FeatureCounts" Node.
 * featureCounts: an efficient general purpose program for assigning sequence reads to genomic features
 * For details, please see paper: 
 * Liao et al, Bioinformatics 2013; doi: 10.1093/bioinformatics/btt656
 *
 * @author Michael Kluge
 */
public class FeatureCountsNodeFactory 
        extends NodeFactory<FeatureCountsNodeModel> {

    /**
     * {@inheritDoc}
     */
    @Override
    public FeatureCountsNodeModel createNodeModel() {
        return new FeatureCountsNodeModel();
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public int getNrNodeViews() {
        return 1;
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public NodeView<FeatureCountsNodeModel> createNodeView(final int viewIndex, final FeatureCountsNodeModel nodeModel) {
        return new HTENodeView<FeatureCountsNodeModel>(nodeModel);
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public boolean hasDialog() {
        return true;
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public NodeDialogPane createNodeDialogPane() {
        return new FeatureCountsNodeDialog();
    }

}

