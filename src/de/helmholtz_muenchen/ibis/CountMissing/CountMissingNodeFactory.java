package de.helmholtz_muenchen.ibis.CountMissing;

import org.knime.core.node.NodeDialogPane;
import org.knime.core.node.NodeFactory;
import org.knime.core.node.NodeView;

/**
 * <code>NodeFactory</code> for the "RNodeTemplate" Node.
 * This is a template for a node, which executes a R-Script.
 *
 * @author Jonas Zierer
 */
public class CountMissingNodeFactory 
        extends NodeFactory<CountMissingNodeModel> {

    /**
     * {@inheritDoc}
     */
    @Override
    public CountMissingNodeModel createNodeModel() {
        return new CountMissingNodeModel();
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public int getNrNodeViews() {
        return 0;
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public NodeView<CountMissingNodeModel> createNodeView(final int viewIndex, final CountMissingNodeModel nodeModel) {
        return null;
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
        return new CountMissingNodeDialog();
    }

}

