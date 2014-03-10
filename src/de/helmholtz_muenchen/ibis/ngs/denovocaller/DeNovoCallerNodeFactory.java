package de.helmholtz_muenchen.ibis.ngs.denovocaller;

import org.knime.core.node.NodeDialogPane;
import org.knime.core.node.NodeFactory;
import org.knime.core.node.NodeView;

/**
 * <code>NodeFactory</code> for the "DeNovoCaller" Node.
 * 
 *
 * @author Maximilian Hastreiter
 */
public class DeNovoCallerNodeFactory 
        extends NodeFactory<DeNovoCallerNodeModel> {

    /**
     * {@inheritDoc}
     */
    @Override
    public DeNovoCallerNodeModel createNodeModel() {
        return new DeNovoCallerNodeModel();
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
    public NodeView<DeNovoCallerNodeModel> createNodeView(final int viewIndex,
            final DeNovoCallerNodeModel nodeModel) {
        return new DeNovoCallerNodeView(nodeModel);
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
        return new DeNovoCallerNodeDialog();
    }

}

