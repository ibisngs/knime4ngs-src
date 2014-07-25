package de.helmholtz_muenchen.ibis.ngs.depthofcoverage;

import org.knime.core.node.NodeDialogPane;
import org.knime.core.node.NodeFactory;
import org.knime.core.node.NodeView;

/**
 * <code>NodeFactory</code> for the "DepthOfCoverage" Node.
 * 
 *
 * @author 
 */
public class DepthOfCoverageNodeFactory 
        extends NodeFactory<DepthOfCoverageNodeModel> {

    /**
     * {@inheritDoc}
     */
    @Override
    public DepthOfCoverageNodeModel createNodeModel() {
        return new DepthOfCoverageNodeModel();
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
    public NodeView<DepthOfCoverageNodeModel> createNodeView(final int viewIndex,
            final DepthOfCoverageNodeModel nodeModel) {
        return new DepthOfCoverageNodeView(nodeModel);
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
        return new DepthOfCoverageNodeDialog();
    }

}

