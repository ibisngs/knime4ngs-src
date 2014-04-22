package de.helmholtz_muenchen.ibis.misc.adjustPvalues;

import org.knime.core.node.NodeDialogPane;
import org.knime.core.node.NodeFactory;
import org.knime.core.node.NodeView;

import de.helmholtz_muenchen.ibis.utils.abstractNodes.RNode.RNodeView;

/**
 * <code>NodeFactory</code> for the "RNodeTemplate" Node.
 * This is a template for a node, which executes a R-Script.
 *
 * @author Jonas Zierer
 */
public class PValAdjustNodeFactory extends NodeFactory<PValAdjustNodeModel> {

    /**
     * {@inheritDoc}
     */
    @Override
    public PValAdjustNodeModel createNodeModel() {
        return new PValAdjustNodeModel();
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
    public NodeView<PValAdjustNodeModel> createNodeView(final int viewIndex, final PValAdjustNodeModel nodeModel) {
        return new RNodeView<PValAdjustNodeModel>(nodeModel);
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
        return new PValAdjustNodeDialog();
    }

}

