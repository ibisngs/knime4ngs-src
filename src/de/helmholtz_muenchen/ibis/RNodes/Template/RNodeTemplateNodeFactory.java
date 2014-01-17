package de.helmholtz_muenchen.ibis.RNodes.Template;

import org.knime.core.node.NodeDialogPane;
import org.knime.core.node.NodeFactory;
import org.knime.core.node.NodeView;

/**
 * <code>NodeFactory</code> for the "RNodeTemplate" Node.
 * This is a template for a node, which executes a R-Script.
 *
 * @author Jonas Zierer
 */
public class RNodeTemplateNodeFactory 
        extends NodeFactory<RNodeTemplateNodeModel> {

    /**
     * {@inheritDoc}
     */
    @Override
    public RNodeTemplateNodeModel createNodeModel() {
        return new RNodeTemplateNodeModel();
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
    public NodeView<RNodeTemplateNodeModel> createNodeView(final int viewIndex,
            final RNodeTemplateNodeModel nodeModel) {
        return new RNodeTemplateNodeView(nodeModel);
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
        return new RNodeTemplateNodeDialog();
    }

}

