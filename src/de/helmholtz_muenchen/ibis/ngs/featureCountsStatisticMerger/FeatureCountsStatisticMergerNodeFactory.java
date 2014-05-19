package de.helmholtz_muenchen.ibis.ngs.featureCountsStatisticMerger;

import org.knime.core.node.NodeDialogPane;
import org.knime.core.node.NodeFactory;
import org.knime.core.node.NodeView;

/**
 * <code>NodeFactory</code> for the "FeatureCountsStatisticMergerNode" Node.
 * 
 *
 * @author Michael Kluge
 */
public class FeatureCountsStatisticMergerNodeFactory extends NodeFactory<FeatureCountsStatisticMergerNodeModel> {

    /**
     * {@inheritDoc}
     */
    @Override
    public FeatureCountsStatisticMergerNodeModel createNodeModel() {
        return new FeatureCountsStatisticMergerNodeModel();
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
    public NodeView<FeatureCountsStatisticMergerNodeModel> createNodeView(final int viewIndex, final FeatureCountsStatisticMergerNodeModel nodeModel) {
        return new FeatureCountsStatisticMergerNodeView(nodeModel);
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
        return new FeatureCountsStatisticMergerNodeDialog();
    }

}

