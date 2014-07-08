package de.helmholtz_muenchen.ibis.ngs.featureCountsMerger;

import org.knime.core.node.NodeDialogPane;
import org.knime.core.node.NodeFactory;
import org.knime.core.node.NodeView;

/**
 * <code>NodeFactory</code> for the "FeatureCountsMerger" Node.
 * 
 *
 * @author Michael Kluge
 */
public class FeatureCountsMergerNodeFactory 
        extends NodeFactory<FeatureCountsMergerNodeModel> {

    /**
     * {@inheritDoc}
     */
    @Override
    public FeatureCountsMergerNodeModel createNodeModel() {
        return new FeatureCountsMergerNodeModel();
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
    public NodeView<FeatureCountsMergerNodeModel> createNodeView(final int viewIndex,
            final FeatureCountsMergerNodeModel nodeModel) {
        return new FeatureCountsMergerNodeView(nodeModel);
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
        return new FeatureCountsMergerNodeDialog();
    }

}

