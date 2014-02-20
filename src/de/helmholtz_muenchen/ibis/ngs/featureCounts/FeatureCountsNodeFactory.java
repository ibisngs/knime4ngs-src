package de.helmholtz_muenchen.ibis.ngs.featureCounts;

import org.knime.core.node.NodeDialogPane;
import org.knime.core.node.NodeFactory;
import org.knime.core.node.NodeView;

/**
 * <code>NodeFactory</code> for the "Star" Node.
 * STAR aligns RNA-seq reads to a reference genome using uncompressed suffix arrays. 
 * For details, please see paper: 
 * Dobin et al, Bioinformatics 2012; doi: 10.1093/bioinformatics/bts635
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
        return new FeatureCountsNodeView(nodeModel);
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

