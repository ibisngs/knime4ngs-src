package de.helmholtz_muenchen.ibis.utils.abstractNodes.HTExecutorNode;

import org.knime.core.node.NodeView;

/**
 * <code>NodeView</code> for the "HTExecutorNode" Node.
 * 
 *
 * @author Tim Jeske
 */
public abstract class HTExecutorNodeView extends NodeView<HTExecutorNodeModel> {

    /**
     * Creates a new view.
     * 
     * @param nodeModel The model (class: {@link HTExecutorNodeModel})
     */
    protected HTExecutorNodeView(final HTExecutorNodeModel nodeModel) {
        super(nodeModel);
        // TODO: generated method stub
    }

    /**
     * {@inheritDoc}
     */
    @Override
    protected void modelChanged() {
        // TODO: generated method stub
    }

    /**
     * {@inheritDoc}
     */
    @Override
    protected void onClose() {
        // TODO: generated method stub
    }

    /**
     * {@inheritDoc}
     */
    @Override
    protected void onOpen() {
        // TODO: generated method stub
    }

}

