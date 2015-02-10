package de.helmholtz_muenchen.ibis.htetrigger;

import org.knime.core.node.NodeDialogPane;
import org.knime.core.node.NodeFactory;
import org.knime.core.node.NodeView;

/**
 * <code>NodeFactory</code> for the "HTETrigger" Node.
 * 
 *
 * @author Tim Jeske
 */
public class HTETriggerNodeFactory 
        extends NodeFactory<HTETriggerNodeModel> {

    /**
     * {@inheritDoc}
     */
    @Override
    public HTETriggerNodeModel createNodeModel() {
        return new HTETriggerNodeModel();
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
    public NodeView<HTETriggerNodeModel> createNodeView(final int viewIndex,
            final HTETriggerNodeModel nodeModel) {
        return new HTETriggerNodeView(nodeModel);
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
        return new HTETriggerNodeDialog();
    }

}

