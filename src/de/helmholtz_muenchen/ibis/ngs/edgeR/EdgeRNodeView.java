package de.helmholtz_muenchen.ibis.ngs.edgeR;

import org.knime.core.node.NodeView;

/**
 * <code>NodeView</code> for the "EdgeR" Node.
 * 
 *
 * @author Michael Kluge
 */
public class EdgeRNodeView extends NodeView<EdgeRNodeModel> {

    /**
     * Creates a new view.
     * 
     * @param nodeModel The model (class: {@link EdgeRNodeModel})
     */
    protected EdgeRNodeView(final EdgeRNodeModel nodeModel) {
        super(nodeModel);
    }

    /**
     * {@inheritDoc}
     */
    @Override
    protected void modelChanged() {
    }

    /**
     * {@inheritDoc}
     */
    @Override
    protected void onClose() {
    }

    /**
     * {@inheritDoc}
     */
    @Override
    protected void onOpen() {
    }

}

