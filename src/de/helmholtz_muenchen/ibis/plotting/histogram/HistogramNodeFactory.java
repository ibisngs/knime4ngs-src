package de.helmholtz_muenchen.ibis.plotting.histogram;

import org.knime.core.node.NodeDialogPane;
import org.knime.core.node.NodeFactory;
import org.knime.core.node.NodeView;

import de.helmholtz_muenchen.ibis.utils.abstractNodes.RNode.RNodeView;

/**
 * <code>NodeFactory</code> for the "Boxplot" Node.
 * 
 *
 * @author Jonas Zierer
 */
public class HistogramNodeFactory 
        extends NodeFactory<HistogramNodeModel> {

    /**
     * {@inheritDoc}
     */
    @Override
    public HistogramNodeModel createNodeModel() {
        return new HistogramNodeModel();
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
    public NodeView<HistogramNodeModel> createNodeView(final int viewIndex, final HistogramNodeModel nodeModel) {
        return new RNodeView<HistogramNodeModel>(nodeModel);
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
        return new HistogramNodeDialog();
    }

}

