package de.helmholtz_muenchen.ibis.ngs.plotdepthofcoverage;

import org.knime.core.node.NodeDialogPane;
import org.knime.core.node.NodeFactory;
import org.knime.core.node.NodeView;

import de.helmholtz_muenchen.ibis.utils.abstractNodes.RNode.RNodeView;

/**
 * <code>NodeFactory</code> for the "PlotDepthOfCoverage" Node.
 * 
 *
 * @author 
 */
public class PlotDepthOfCoverageNodeFactory 
        extends NodeFactory<PlotDepthOfCoverageNodeModel> {

    /**
     * {@inheritDoc}
     */
    @Override
    public PlotDepthOfCoverageNodeModel createNodeModel() {
        return new PlotDepthOfCoverageNodeModel();
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
    public NodeView<PlotDepthOfCoverageNodeModel> createNodeView(final int viewIndex,
            final PlotDepthOfCoverageNodeModel nodeModel) {
        return new RNodeView(nodeModel);
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
        return new PlotDepthOfCoverageNodeDialog();
    }

}

